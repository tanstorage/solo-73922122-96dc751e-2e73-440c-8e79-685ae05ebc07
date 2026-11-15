"""从 samples/demo.fa 生成两份多样本 VCF，演示批次样本配对。

运行：python3 scripts/make_pair_demo.py

场景覆盖：
  - 4 个个体在 20 个共有双等位 SNP 位点上的基因型，两两差异明显；
  - 两批 VCF 样本标签同为 S1–S4，但批次 B 的样本列数据被错位
    （模拟样本标签互换），配对结果会指出真实的对应关系；
  - 其中一个个体在两侧有 1 处基因型不一致（仍在阈值内）；
  - 每个样本都有 2 个 ./. 缺失调用（不计入可比位点）；
  - 一条 indel 记录（非双等位单碱基，不参与配对）。

预期配对（唯一最优方案，4 对，累计不一致 1）：
  A.S1 ↔ B.S2、A.S2 ↔ B.S4、A.S3 ↔ B.S1、A.S4 ↔ B.S3
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from app.bio import parse_fasta  # noqa: E402

HEADER = (
    "##fileformat=VCFv4.2\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
)

# 20 个共享 SNP 位点（1-based），彼此间隔 6 bp，落在 132 bp 参考内。
POSITIONS = [2 + 6 * i for i in range(20)]

# 每个个体的基因型剂量由 (斜率, 截距) 决定：dosage = (slope * site + intercept) % 3
# 四组参数两两不同，保证任意两个个体在多数位点上基因型不同。
PROFILES = [(1, 0), (1, 1), (2, 0), (2, 1)]
GT_BY_DOSAGE = {0: "0/0", 1: "0/1", 2: "1/1"}

# 每个个体固定缺失的位点序号（两侧一致，均不计入可比位点）。
MISSING = {0: {5, 6}, 1: {11, 12}, 2: {15, 16}, 3: {18, 19}}

# 批次 B 中被改写的基因型：(个体序号, 位点序号) -> 剂量。
# 只给个体 2 改 1 处，演示“不一致 ≤ 1 仍可配对”。
PERTURBED = {(2, 0): 2}


def dosage(ind: int, site: int) -> int:
    slope, intercept = PROFILES[ind]
    return (slope * site + intercept) % 3


def sample_gt(ind: int, site: int, altered: bool) -> str:
    if site in MISSING[ind]:
        return "./."
    dose = dosage(ind, site)
    if altered and (ind, site) in PERTURBED:
        dose = PERTURBED[(ind, site)]
    return GT_BY_DOSAGE[dose]


def build(sequence: str, order: list[int], altered: bool) -> str:
    # 表头始终使用 S1–S4；order 决定每一列实际装的是哪个个体的数据，
    # 以此模拟“标签与数据错位”。
    names = [f"S{i + 1}" for i in range(len(order))]
    rows = [
        HEADER
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(names)
    ]
    for site, pos in enumerate(POSITIONS):
        ref = sequence[pos - 1]
        alt = next(b for b in "ACGT" if b != ref)
        gts = [sample_gt(ind, site, altered) for ind in order]
        rows.append(
            f"chrDemo\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(gts)
        )
    # 一条 indel 记录：不参与配对（不会被当成 SNP）
    pos = 60
    ref = sequence[pos - 1]
    rows.append(
        f"chrDemo\t{pos}\t.\t{ref}\t{ref}A\t.\tPASS\t.\tGT\t"
        + "\t".join("0/1" for _ in order)
    )
    return "\n".join(rows) + "\n"


def main() -> None:
    ref = parse_fasta((ROOT / "samples" / "demo.fa").read_bytes())
    # 批次 A：第 i 列装个体 i 的数据（标签与数据一致）。
    # 批次 B：标签仍为 S1–S4，但第 i 列装的是 order_b[i] 号个体的数据，
    # 即标签与数据的对应关系被整体错位（模拟样本标签互换）。
    order_a = [0, 1, 2, 3]
    order_b = [2, 0, 3, 1]

    (ROOT / "samples" / "pair_batch_a.vcf").write_text(
        build(ref.sequence, order_a, altered=False), encoding="utf-8"
    )
    (ROOT / "samples" / "pair_batch_b.vcf").write_text(
        build(ref.sequence, order_b, altered=True), encoding="utf-8"
    )
    print("wrote 21 records -> samples/pair_batch_a.vcf（样本 S1–S4）")
    print("wrote 21 records -> samples/pair_batch_b.vcf（样本 S1–S4，数据列被错位）")
    print("预期：A.S1↔B.S2、A.S2↔B.S4、A.S3↔B.S1、A.S4↔B.S3，累计不一致 1")


if __name__ == "__main__":
    main()
