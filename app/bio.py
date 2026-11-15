"""VCF 核对的核心逻辑（仅依赖标准库）。

处理流程：
1. 解析单条染色体参考 FASTA（仅允许 A/C/G/T）。
2. 解析两份单样本 VCF，仅接受：
   - A/C/G/T 构成的单碱基替换；
   - 长度 <= 20 碱基的纯插入 / 纯缺失；
   - 非缺失、未分相的二倍体 GT。
3. 对照参考做 REF 校验，再把 indel 左对齐到唯一的归一位置
   （重复区域中写法不同的同一 indel 归为同一等位变异）。
4. 多 ALT 记录按 GT 分别计算每个 ALT 的 0/1/2 剂量，忽略剂量 0。
5. 同一文件内同一归一变异出现互相矛盾的剂量（单条记录多值或合并后 >2）即整体失败。
"""

from dataclasses import dataclass, field

# 上传体积上限（字节）
MAX_FASTA_BYTES = 1 * 1024 * 1024      # 1 MiB
MAX_VCF_BYTES = 256 * 1024             # 256 KiB

# indel 事件长度上限（碱基）
MAX_INDEL_LEN = 20

_BASES = frozenset("ACGT")


class CheckError(Exception):
    """可直接展示给用户的核对错误；任何一次抛出都代表整次核对失败。"""


@dataclass(frozen=True)
class Reference:
    name: str
    sequence: str


@dataclass(frozen=True)
class Variant:
    """归一后的等位变异。

    pos 为 1-based 位置：
    - SNV：被替换碱基的位置；
    - 插入：插入点左侧（锚定）碱基的位置，ref 是该锚定碱基，alt 为锚定碱基 + 插入序列；
    - 缺失：被删除序列前一个（锚定）碱基的位置，ref 为锚定碱基 + 被删序列，alt 为锚定碱基。
    """

    pos: int
    ref: str
    alt: str
    kind: str  # "SNV" | "INS" | "DEL"

    @property
    def label(self) -> str:
        arrow = ">" if self.kind == "SNV" else "ins" if self.kind == "INS" else "del"
        return f"{self.pos}:{self.ref}{arrow}{self.alt}"


@dataclass
class Call:
    """某文件中一个归一变异的核对结果。"""

    dosage: int
    line_numbers: list[int] = field(default_factory=list)
    raw_descriptions: list[str] = field(default_factory=list)


# --------------------------------------------------------------------------- #
# FASTA
# --------------------------------------------------------------------------- #
def parse_fasta(data: bytes) -> Reference:
    """解析单条染色体 FASTA。"""
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise CheckError(f"参考 FASTA 不是 UTF-8 文本（字节 {exc.start} 处解码失败）。")

    lines = text.splitlines()
    header_idx = next((i for i, ln in enumerate(lines) if ln.startswith(">")), -1)
    if header_idx < 0:
        raise CheckError("参考 FASTA 缺少以 '>' 开头的描述行。")

    name = lines[header_idx][1:].strip()
    if not name:
        raise CheckError("参考 FASTA 描述行 '>' 后缺少染色体名称。")
    name = name.split()[0]

    seq_parts: list[str] = []
    for i, line in enumerate(lines):
        if i == header_idx:
            continue
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            raise CheckError("参考 FASTA 含有多条记录，本工具仅接受单条染色体参考。")
        seq_parts.append(stripped)

    sequence = "".join(seq_parts).upper()
    if not sequence:
        raise CheckError("参考 FASTA 不包含任何序列。")
    bad = next((b for b in sequence if b not in _BASES), None)
    if bad is not None:
        raise CheckError(
            f"参考序列仅允许 A/C/G/T，但发现字符 {bad!r}；请先确认参考为纯 ACGT 序列。"
        )
    return Reference(name=name, sequence=sequence)


# --------------------------------------------------------------------------- #
# 变异归一化
# --------------------------------------------------------------------------- #
def normalize(pos: int, ref: str, alt: str, sequence: str) -> tuple[int, str, str]:
    """把 REF/ALT 等位变异按参考序列左对齐，返回归一后的 (1-based pos, ref, alt)。

    与 bcftools norm / vt normalize 等价的左对齐：先去共同前缀，
    再在锚定碱基可保留时不断把等位整体左移，最后去共同后缀。
    调用方需保证 REF 与参考在 pos 处一致，且 pos/长度不越界。
    """
    # 0-based 起点
    p = pos - 1
    r = ref
    a = alt

    # 去除共同前缀，各自至少保留 1 个锚定碱基
    while len(r) > 1 and len(a) > 1 and r[0] == a[0]:
        r = r[1:]
        a = a[1:]
        p += 1

    # 左移：等位末端碱基相同且左侧仍有参考碱基时，整体向左滑动
    while p > 0 and r[-1] == a[-1]:
        p -= 1
        r = sequence[p] + r[:-1]
        a = sequence[p] + a[:-1]

    # 去除共同后缀，各自至少保留 1 个锚定碱基
    while len(r) > 1 and len(a) > 1 and r[-1] == a[-1]:
        r = r[:-1]
        a = a[:-1]

    return p + 1, r, a


def classify(ref: str, alt: str) -> str:
    """按归一后的形态分类；不支持的形态抛 CheckError。"""
    if ref == alt:
        raise CheckError("存在 REF 与 ALT 相同的无效记录。")
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) == 1 and len(alt) > 1:
        event_len = len(alt) - 1
        if event_len > MAX_INDEL_LEN:
            raise CheckError(f"插入长度 {event_len} bp 超过上限 {MAX_INDEL_LEN} bp。")
        return "INS"
    if len(alt) == 1 and len(ref) > 1:
        event_len = len(ref) - 1
        if event_len > MAX_INDEL_LEN:
            raise CheckError(f"缺失长度 {event_len} bp 超过上限 {MAX_INDEL_LEN} bp。")
        return "DEL"
    raise CheckError(
        "仅接受单碱基替换或纯插入/纯缺失，发现多核苷酸替换等复杂变异"
        f"（REF={ref}, ALT={alt}）。"
    )


# --------------------------------------------------------------------------- #
# GT
# --------------------------------------------------------------------------- #
def parse_diploid_gt(gt: str, line_no: int, n_alt: int) -> dict[int, int]:
    """解析未分相、非缺失的二倍体 GT，返回 {alt 1-based 索引: 剂量}。

    剂量为 0 的 ALT 不写入结果。多 ALT 时同一 ALT 出现两次则剂量为 2。
    """
    if not gt:
        raise CheckError(f"第 {line_no} 行：GT 为空。")
    if "|" in gt:
        raise CheckError(f"第 {line_no} 行：仅接受未分相 GT（'/'），发现分相标记 '|'：{gt}")

    parts = gt.split("/")
    if len(parts) != 2:
        raise CheckError(f"第 {line_no} 行：仅接受二倍体 GT（两个等位），收到 {gt!r}。")

    alleles: list[int] = []
    for part in parts:
        if part in ("", ".", ".*"):
            raise CheckError(f"第 {line_no} 行：GT 含缺失等位，本工具不接受缺失 GT：{gt}")
        try:
            idx = int(part)
        except ValueError:
            raise CheckError(f"第 {line_no} 行：GT 等位 {part!r} 不是整数索引。")
        if idx < 0:
            raise CheckError(f"第 {line_no} 行：GT 等位索引不能为负：{gt}")
        if idx > n_alt:
            raise CheckError(
                f"第 {line_no} 行：GT 等位索引 {idx} 超出 ALT 数量（{n_alt} 个 ALT）。"
            )
        alleles.append(idx)

    dosages: dict[int, int] = {}
    for idx in alleles:
        if idx == 0:  # REF 等位，不计入任何 ALT 剂量
            continue
        dosages[idx] = dosages.get(idx, 0) + 1
    return dosages


# --------------------------------------------------------------------------- #
# VCF
# --------------------------------------------------------------------------- #
def _parse_record(
    fields: list[str], line_no: int, refobj: Reference
) -> dict[Variant, int]:
    """解析一条 VCF 数据行，返回该行内 {归一变异: 剂量}（已剔除剂量 0 并在记录内聚合）。"""
    if len(fields) < 8:
        raise CheckError(f"第 {line_no} 行：列数不足 8 列，不是合法 VCF 记录。")

    chrom, pos_s, _vid, ref, alt_s, _qual, _filt, info = fields[:8]
    if chrom != refobj.name:
        raise CheckError(
            f"第 {line_no} 行：CHROM={chrom!r} 与参考 FASTA 名称 {refobj.name!r} 不一致。"
        )

    try:
        pos = int(pos_s)
    except ValueError:
        raise CheckError(f"第 {line_no} 行：POS={pos_s!r} 不是整数。")
    if pos < 1 or pos > len(refobj.sequence):
        raise CheckError(
            f"第 {line_no} 行：POS={pos} 超出参考范围（参考长度 {len(refobj.sequence)}）。"
        )

    if not ref or any(ch not in _BASES for ch in ref):
        raise CheckError(
            f"第 {line_no} 行：REF 等位 {ref!r} 不是纯 A/C/G/T 序列。"
        )
    alts = alt_s.split(",")
    if any(not a or any(ch not in _BASES for ch in a) for a in alts):
        raise CheckError(
            f"第 {line_no} 行：ALT {alt_s!r} 中存在空值或非 A/C/G/T 字符"
            "（本工具不接受 *、<DEL>、符号化等位及缺失占位）。"
        )

    end = pos + len(ref) - 1
    if end > len(refobj.sequence):
        raise CheckError(
            f"第 {line_no} 行：REF 延伸至 {end}，超出参考长度 {len(refobj.sequence)}。"
        )
    ref_sub = refobj.sequence[pos - 1:end]
    if ref_sub != ref:
        raise CheckError(
            f"第 {line_no} 行：REF={ref} 与参考序列在 pos={pos} 处的 {ref_sub} 不一致，"
            "整次核对失败。"
        )

    # FORMAT/GT
    if len(fields) == 8:
        raise CheckError(f"第 {line_no} 行：缺少 FORMAT 与样本列，本工具要求单样本二倍体 GT。")
    fmt = fields[8]
    fmt_keys = fmt.split(":")
    if not fmt_keys or fmt_keys[0] != "GT":
        raise CheckError(f"第 {line_no} 行：FORMAT 必须以 GT 开头，实际为 {fmt!r}。")
    sample = fields[9]
    sample_values = sample.split(":")
    if len(sample_values) > len(fmt_keys):
        raise CheckError(f"第 {line_no} 行：样本列的子字段多于 FORMAT 定义。")
    gt = sample_values[0] if sample_values else ""
    dosages = parse_diploid_gt(gt, line_no, len(alts))

    # 同一记录中不同 ALT 归一后可能指向同一变异（如重复区里的等价写法），
    # 先在记录内聚合剂量；二倍体单条记录的剂量不可能超过 2。
    per_record: dict[Variant, int] = {}
    for alt_idx, dose in dosages.items():
        alt = alts[alt_idx - 1]
        npos, nref, nalt = normalize(pos, ref, alt, refobj.sequence)
        kind = classify(nref, nalt)  # 同时承担长度与形态校验
        variant = Variant(npos, nref, nalt, kind)
        per_record[variant] = per_record.get(variant, 0) + dose
    return per_record


def parse_vcf(data: bytes, filename: str, refobj: Reference) -> dict[Variant, Call]:
    """解析单样本 VCF，返回 {归一变异: Call}。任何硬性错误都抛出 CheckError。"""
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise CheckError(f"{filename} 不是 UTF-8 文本（字节 {exc.start} 处解码失败）。")

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
            "本工具要求恰好一个样本列（共 10 列）。"
        )
    if len(header_fields) > 10:
        raise CheckError(
            f"{filename}：'#CHROM' 表头有 {len(header_fields) - 9} 个样本列，"
            "本工具仅接受单样本 VCF。"
        )
    sample_name = header_fields[9]

    calls: dict[Variant, Call] = {}
    for offset, raw in enumerate(lines[header_idx + 1:], start=header_idx + 2):
        line = raw.rstrip("\r")
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        record_calls = _parse_record(fields, offset, refobj)
        raw_desc = f"{fields[1]}:{fields[3]}>{fields[4]}"
        for variant, dose in record_calls.items():
            existing = calls.get(variant)
            if existing is None:
                calls[variant] = Call(
                    dosage=dose,
                    line_numbers=[offset],
                    raw_descriptions=[raw_desc],
                )
            else:
                # 同一归一变异：单条记录内已可能是剂量 2；
                # 若多条记录的剂量无法在二倍体中并存（合计 > 2）即互相矛盾。
                if existing.dosage + dose > 2:
                    raise CheckError(
                        f"{filename}：归一变异 {variant.label} 在第 "
                        f"{', '.join(map(str, existing.line_numbers))} 行与第 {offset} 行"
                        f"出现相互矛盾的剂量（{existing.dosage} + {dose} 超过二倍体上限 2），"
                        "整次核对失败。"
                    )
                existing.dosage += dose
                existing.line_numbers.append(offset)
                existing.raw_descriptions.append(raw_desc)

    return calls, sample_name


# --------------------------------------------------------------------------- #
# 双侧比对
# --------------------------------------------------------------------------- #
def compare(
    refobj: Reference,
    left: dict[Variant, Call],
    right: dict[Variant, Call],
) -> dict[str, list[dict]]:
    """按归一变异归并两份结果，分为一致 / 仅单侧 / 剂量不符三组。"""
    consistent: list[dict] = []
    left_only: list[dict] = []
    right_only: list[dict] = []
    dosage_mismatch: list[dict] = []

    for variant in sorted(set(left) | set(right), key=lambda v: (v.pos, v.ref, v.alt)):
        lc = left.get(variant)
        rc = right.get(variant)
        entry = {
            "variant": variant,
            "left": lc,
            "right": rc,
        }
        if lc is not None and rc is not None:
            if lc.dosage == rc.dosage:
                consistent.append(entry)
            else:
                dosage_mismatch.append(entry)
        elif lc is not None:
            left_only.append(entry)
        else:
            right_only.append(entry)

    return {
        "consistent": consistent,
        "left_only": left_only,
        "right_only": right_only,
        "dosage_mismatch": dosage_mismatch,
    }
