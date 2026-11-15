"""批次样本配对的核心逻辑（仅依赖标准库）。

目标：怀疑两批测序结果的样本标签互换时，用两份多样本 VCF 中共有的
双等位单碱基位点（SNP）的基因型核对样本身份。

规则：
1. 两份 VCF 各 1–6 个样本，各不超过 256 KiB；
2. 只使用两侧都出现、REF 与参考序列相符、ALT 相同且为双等位单碱基的位点；
3. GT 仅接受未分相的 ``0/0``、``0/1``、``1/1`` 与缺失调用 ``./.``；
   ``./.`` 不计入可比位点；其他 GT 写法一律判失败；
4. 两个样本至少有 8 个可比位点且不一致位点不超过 1 个才可配对；
5. 全局先最大化配对数，再最小化所有配对的累计不一致位点数；
   页面需要列出每个样本在所有最优方案中的可能对应对象（含未配对）。

任何重复位点、非法 GT、REF 与参考不符都会让整次分析失败，不输出部分映射。
"""

from dataclasses import dataclass, field
from typing import Optional

from .bio import CheckError, Reference

# 每份 VCF 的样本数上限
MAX_SAMPLES_PER_VCF = 6

# 配对阈值
MIN_COMPARABLE_SITES = 8
MAX_INCONSISTENT_SITES = 1

_BASES = frozenset("ACGT")
_ALLOWED_GTS = frozenset(("0/0", "0/1", "1/1", "./."))

# GT 在内部用剂量表示：0 = 0/0，1 = 0/1，2 = 1/1，None = ./.
GT_MISSING: Optional[int] = None


@dataclass(frozen=True)
class SnpRecord:
    """一条可参与配对的双等位 SNV 记录（每个样本一份 GT 剂量）。"""

    ref: str
    alt: str
    gts: tuple[Optional[int], ...]


@dataclass
class BatchVcf:
    """一份多样本 VCF 的解析结果。"""

    filename: str
    samples: list[str]
    # 仅保留双等位单碱基记录；键为 1-based POS。
    sites: dict[int, SnpRecord] = field(default_factory=dict)


@dataclass(frozen=True)
class SharedSite:
    pos: int
    ref: str
    alt: str

    @property
    def label(self) -> str:
        return f"{self.pos}:{self.ref}>{self.alt}"


@dataclass
class PairStats:
    """批次 A 样本与批次 B 样本之间逐位点的统计。"""

    a: str
    b: str
    comparable: int
    inconsistent: int
    eligible: bool
    in_optimal: bool = False


@dataclass
class Counterpart:
    sample: str
    comparable: int
    inconsistent: int


@dataclass
class SamplePossibilities:
    """某个原批次样本在全部最优方案中的可能对应对象。"""

    sample: str
    batch: str  # "A" | "B"
    counterparts: list[Counterpart] = field(default_factory=list)
    unpaired_possible: bool = False


@dataclass
class OptimalScheme:
    """一个最优配对方案。pairs 为 (批次 A 样本, 批次 B 样本)。"""

    pairs: list[tuple[str, str]]
    inconsistent: int


@dataclass
class PairingResult:
    reference_name: str
    reference_length: int
    samples_a: list[str]
    samples_b: list[str]
    shared_sites: list[SharedSite]
    # 所有 A×B 组合的逐位点统计（含不可配对组合），便于页面完整展示。
    pair_stats: list[PairStats]
    max_pairs: int
    min_inconsistent: int
    schemes: list[OptimalScheme]
    possibilities_a: list[SamplePossibilities]
    possibilities_b: list[SamplePossibilities]

    @property
    def scheme_count(self) -> int:
        return len(self.schemes)

    @property
    def max_paired_samples(self) -> int:
        return self.max_pairs * 2


# --------------------------------------------------------------------------- #
# 多样本 VCF 解析
# --------------------------------------------------------------------------- #
def _parse_gt(
    gt: str,
    filename: str,
    line_no: int,
    sample: str,
) -> Optional[int]:
    """校验 GT 并转换成剂量；``./.`` 返回 None。

    只接受未分相的 ``0/0``、``0/1``、``1/1`` 与缺失调用 ``./.``。
    其他任何写法（``1/2``、``0|1``、``./1``、非二倍体、空值等）都属于非法 GT，
    整次分析失败。
    """
    if gt not in _ALLOWED_GTS:
        raise CheckError(
            f"{filename} 第 {line_no} 行：样本 {sample!r} 的 GT={gt!r} 非法；"
            "仅接受未分相的 0/0、0/1、1/1 或缺失调用 ./."
        )
    if gt == "./.":
        return GT_MISSING
    return {"0/0": 0, "0/1": 1, "1/1": 2}[gt]


def parse_batch_vcf(data: bytes, filename: str, refobj: Reference) -> BatchVcf:
    """解析一份 1–6 样本的 VCF。

    所有记录都会校验 CHROM、POS、REF 与每个样本的 GT（仅允许 ``0/0``、
    ``0/1``、``1/1``、``./.``）；只有满足双等位单碱基条件的记录才会进入
    ``sites``。任何硬性错误抛出 CheckError，调用方不得展示部分结果。
    """
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise CheckError(
            f"{filename} 不是 UTF-8 文本（字节 {exc.start} 处解码失败）。"
        )

    lines = text.splitlines()
    header_idx = -1
    for i, line in enumerate(lines):
        if line.startswith("#CHROM"):
            header_idx = i
            break
    if header_idx < 0:
        raise CheckError(f"{filename}：缺少 '#CHROM' 表头行，不是合法 VCF。")

    header_fields = lines[header_idx].rstrip("\r").split("\t")
    if len(header_fields) < 10:
        raise CheckError(
            f"{filename}：'#CHROM' 表头只有 {len(header_fields)} 列，"
            "缺少样本列；配对分析要求 1–6 个样本。"
        )
    n_samples = len(header_fields) - 9
    if n_samples > MAX_SAMPLES_PER_VCF:
        raise CheckError(
            f"{filename}：表头有 {n_samples} 个样本列，超过上限 "
            f"{MAX_SAMPLES_PER_VCF} 个。"
        )

    sample_names = [name.strip() for name in header_fields[9:]]
    if any(not name for name in sample_names):
        raise CheckError(f"{filename}：样本列中存在空样本名。")
    duplicate = next(
        (name for name in sample_names if sample_names.count(name) > 1), None
    )
    if duplicate is not None:
        raise CheckError(f"{filename}：样本名 {duplicate!r} 重复，无法唯一标识样本。")

    batch = BatchVcf(filename=filename, samples=sample_names)
    seen_positions: dict[int, int] = {}

    for offset, raw in enumerate(lines[header_idx + 1:], start=header_idx + 2):
        line = raw.rstrip("\r")
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 9 + n_samples:
            raise CheckError(
                f"{filename} 第 {offset} 行：列数为 {len(fields)}，"
                f"少于表头要求的 {9 + n_samples} 列。"
            )
        if len(fields) > 9 + n_samples:
            raise CheckError(
                f"{filename} 第 {offset} 行：列数为 {len(fields)}，"
                f"多于表头要求的 {9 + n_samples} 列。"
            )

        chrom, pos_s, _vid, ref, alt_s, _qual, _filt, _info = fields[:8]
        if chrom != refobj.name:
            raise CheckError(
                f"{filename} 第 {offset} 行：CHROM={chrom!r} 与参考 FASTA 名称 "
                f"{refobj.name!r} 不一致。"
            )

        try:
            pos = int(pos_s)
        except ValueError:
            raise CheckError(f"{filename} 第 {offset} 行：POS={pos_s!r} 不是整数。")
        if pos < 1 or pos > len(refobj.sequence):
            raise CheckError(
                f"{filename} 第 {offset} 行：POS={pos} 超出参考范围"
                f"（参考长度 {len(refobj.sequence)}）。"
            )

        if not ref or any(ch not in _BASES for ch in ref):
            raise CheckError(
                f"{filename} 第 {offset} 行：REF 等位 {ref!r} 不是纯 A/C/G/T 序列。"
            )
        end = pos + len(ref) - 1
        if end > len(refobj.sequence):
            raise CheckError(
                f"{filename} 第 {offset} 行：REF 延伸至 {end}，"
                f"超出参考长度 {len(refobj.sequence)}。"
            )
        ref_sub = refobj.sequence[pos - 1:end]
        if ref_sub != ref:
            raise CheckError(
                f"{filename} 第 {offset} 行：REF={ref} 与参考序列在 pos={pos} 处的 "
                f"{ref_sub} 不一致，整次配对分析失败。"
            )

        if pos in seen_positions:
            raise CheckError(
                f"{filename} 第 {seen_positions[pos]} 行与第 {offset} 行重复出现"
                f"同一位点 POS={pos}，整次配对分析失败。"
            )
        seen_positions[pos] = offset

        alts = alt_s.split(",")
        is_biallelic_snp = (
            len(ref) == 1
            and len(alts) == 1
            and len(alts[0]) == 1
            and alts[0] in _BASES
            and alts[0] != ref
        )

        fmt = fields[8]
        fmt_keys = fmt.split(":")
        if "GT" not in fmt_keys:
            raise CheckError(
                f"{filename} 第 {offset} 行：FORMAT={fmt!r} 未包含 GT。"
            )
        gt_idx = fmt_keys.index("GT")

        gts: list[Optional[int]] = []
        for sample_idx, sample_name in enumerate(sample_names):
            sample_field = fields[9 + sample_idx]
            values = sample_field.split(":")
            if len(values) > len(fmt_keys):
                raise CheckError(
                    f"{filename} 第 {offset} 行：样本 {sample_name!r} 的子字段多于 "
                    f"FORMAT 定义（{fmt!r}）。"
                )
            if gt_idx >= len(values):
                raise CheckError(
                    f"{filename} 第 {offset} 行：样本 {sample_name!r} 缺少 GT 子字段"
                    f"（FORMAT={fmt!r}）。"
                )
            gts.append(
                _parse_gt(values[gt_idx], filename, offset, sample_name)
            )

        # indel / MNV / 多 ALT / 符号化 ALT 等记录不参与配对，
        # 但仍需通过上面的 GT 严格校验。
        if is_biallelic_snp:
            batch.sites[pos] = SnpRecord(ref=ref, alt=alts[0], gts=tuple(gts))

    return batch


# --------------------------------------------------------------------------- #
# 全局配对
# --------------------------------------------------------------------------- #
def _collect_shared_sites(
    batch_a: BatchVcf, batch_b: BatchVcf
) -> list[SharedSite]:
    shared: list[SharedSite] = []
    for pos in sorted(set(batch_a.sites) & set(batch_b.sites)):
        rec_a = batch_a.sites[pos]
        rec_b = batch_b.sites[pos]
        if rec_a.ref == rec_b.ref and rec_a.alt == rec_b.alt:
            shared.append(SharedSite(pos=pos, ref=rec_a.ref, alt=rec_a.alt))
    return shared


def _pair_stats(
    batch_a: BatchVcf,
    batch_b: BatchVcf,
    shared: list[SharedSite],
) -> list[PairStats]:
    stats: list[PairStats] = []
    for ai, sample_a in enumerate(batch_a.samples):
        for bi, sample_b in enumerate(batch_b.samples):
            comparable = 0
            inconsistent = 0
            for site in shared:
                gt_a = batch_a.sites[site.pos].gts[ai]
                gt_b = batch_b.sites[site.pos].gts[bi]
                if gt_a is None or gt_b is None:
                    continue
                comparable += 1
                if gt_a != gt_b:
                    inconsistent += 1
            eligible = (
                comparable >= MIN_COMPARABLE_SITES
                and inconsistent <= MAX_INCONSISTENT_SITES
            )
            stats.append(
                PairStats(
                    a=sample_a,
                    b=sample_b,
                    comparable=comparable,
                    inconsistent=inconsistent,
                    eligible=eligible,
                )
            )
    return stats


def _enumerate_optimal(
    n_a: int,
    n_b: int,
    edges: dict[tuple[int, int], PairStats],
) -> tuple[int, int, list[tuple[tuple[int, int], ...]]]:
    """枚举全部最优匹配。

    最优含义：先最大化配对数，再最小化所有配对的累计不一致位点数。
    批次规模最多 6×6，全部匹配数上界 13327，直接穷举即可。
    """
    best_pairs = -1
    best_cost = -1
    schemes: list[tuple[tuple[int, int], ...]] = []

    def visit(ai: int, used_b: set[int], pairs: list[tuple[int, int]]) -> None:
        nonlocal best_pairs, best_cost, schemes
        if ai == n_a:
            cost = sum(edges[(i, j)].inconsistent for i, j in pairs)
            if (
                len(pairs) > best_pairs
                or (len(pairs) == best_pairs and cost < best_cost)
            ):
                best_pairs = len(pairs)
                best_cost = cost
                schemes = [tuple(pairs)]
            elif len(pairs) == best_pairs and cost == best_cost:
                schemes.append(tuple(pairs))
            return

        # 该样本不配对
        visit(ai + 1, used_b, pairs)

        # 与任一尚未使用的批次 B 样本配对（仅限可配对边）
        for bi in range(n_b):
            if bi in used_b or (ai, bi) not in edges:
                continue
            used_b.add(bi)
            pairs.append((ai, bi))
            visit(ai + 1, used_b, pairs)
            pairs.pop()
            used_b.remove(bi)

    visit(0, set(), [])
    # 按样本名排序，保证同一输入的结果稳定、便于人工比对
    schemes = [
        tuple(sorted(scheme, key=lambda p: (p[0], p[1]))) for scheme in schemes
    ]
    schemes.sort()
    return best_pairs, best_cost, schemes


def pair_batches(
    refobj: Reference, batch_a: BatchVcf, batch_b: BatchVcf
) -> PairingResult:
    """在已解析的两批样本之间求全部最优配对方案及每个样本的可能对应对象。"""
    shared = _collect_shared_sites(batch_a, batch_b)
    stats = _pair_stats(batch_a, batch_b, shared)
    edges = {
        (ai, bi): stats[ai * len(batch_b.samples) + bi]
        for ai in range(len(batch_a.samples))
        for bi in range(len(batch_b.samples))
        if stats[ai * len(batch_b.samples) + bi].eligible
    }

    max_pairs, min_inconsistent, schemes_idx = _enumerate_optimal(
        len(batch_a.samples), len(batch_b.samples), edges
    )

    # 标记哪些边真正出现在至少一个最优方案中，以及哪些样本可能未配对。
    edges_in_optimal: set[tuple[int, int]] = set()
    unpaired_a: set[int] = set()
    unpaired_b: set[int] = set()
    for scheme in schemes_idx:
        matched_a = {i for i, _ in scheme}
        matched_b = {j for _, j in scheme}
        edges_in_optimal.update(scheme)
        unpaired_a.update(set(range(len(batch_a.samples))) - matched_a)
        unpaired_b.update(set(range(len(batch_b.samples))) - matched_b)

    for ai, bi in edges_in_optimal:
        stats[ai * len(batch_b.samples) + bi].in_optimal = True

    schemes = [
        OptimalScheme(
            pairs=[
                (batch_a.samples[i], batch_b.samples[j]) for i, j in scheme
            ],
            inconsistent=sum(edges[(i, j)].inconsistent for i, j in scheme),
        )
        for scheme in schemes_idx
    ]

    def possibilities(
        samples: list[str], batch_letter: str, other: list[str], is_a: bool
    ) -> list[SamplePossibilities]:
        result: list[SamplePossibilities] = []
        for idx, sample in enumerate(samples):
            candidates: list[Counterpart] = []
            for other_idx, other_name in enumerate(other):
                key = (idx, other_idx) if is_a else (other_idx, idx)
                if key not in edges_in_optimal:
                    continue
                pair = edges[key]
                candidates.append(
                    Counterpart(
                        sample=other_name,
                        comparable=pair.comparable,
                        inconsistent=pair.inconsistent,
                    )
                )
            result.append(
                SamplePossibilities(
                    sample=sample,
                    batch=batch_letter,
                    counterparts=candidates,
                    unpaired_possible=(
                        idx in (unpaired_a if is_a else unpaired_b)
                    ),
                )
            )
        return result

    return PairingResult(
        reference_name=refobj.name,
        reference_length=len(refobj.sequence),
        samples_a=list(batch_a.samples),
        samples_b=list(batch_b.samples),
        shared_sites=shared,
        pair_stats=stats,
        max_pairs=max_pairs,
        min_inconsistent=min_inconsistent,
        schemes=schemes,
        possibilities_a=possibilities(
            batch_a.samples, "A", batch_b.samples, is_a=True
        ),
        possibilities_b=possibilities(
            batch_b.samples, "B", batch_a.samples, is_a=False
        ),
    )
