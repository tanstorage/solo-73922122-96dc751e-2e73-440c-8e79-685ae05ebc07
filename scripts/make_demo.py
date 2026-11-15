"""从 samples/demo.fa 生成两份演示 VCF（坐标由参考序列自动计算，保证 REF 合法）。

运行：python3 scripts/make_demo.py
场景覆盖：
  SNV 一致 / 同一插入的不同写法 / 同一缺失的不同写法 / 多 ALT /
  剂量不符 / 仅单侧出现 / 剂量 0 被忽略
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.bio import parse_fasta  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
)


def find(seq: str, motif: str, start: int = 1) -> int:
    """返回 motif 在参考中的 1-based 起点。"""
    idx = seq.find(motif, start - 1)
    if idx < 0:
        raise SystemExit(f"motif {motif!r} not found in reference")
    return idx + 1


def row(pos, ref, alt, gt, extra_fmt=None):
    fmt = "GT" if extra_fmt is None else "GT:DP"
    sample = gt if extra_fmt is None else f"{gt}:{extra_fmt}"
    return f"chrDemo\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\t{fmt}\t{sample}"


def main() -> None:
    ref = parse_fasta((ROOT / "samples" / "demo.fa").read_bytes())
    seq = ref.sequence

    # 找一段纯 A 重复区，演示插入/缺失在重复区内的等价写法
    a_run = find(seq, "AAAA")
    # 找一个 AT 二碱基重复区，演示不同锚点写法的缺失
    at_run = find(seq, "ATATATAT")

    rows_a: list[str] = []
    rows_b: list[str] = []

    # 1) SNV 两侧一致（0/1）
    p = find(seq, "ACGTAC", 1)
    snv_pos = p + 1  # 参考碱基 C -> G
    rows_a.append(row(snv_pos, "C", "G", "0/1"))
    rows_b.append(row(snv_pos, "C", "G", "0/1"))

    # 2) 同一插入的两种重复区写法：A 在 A 串前插入 A（左对齐后应相同）
    anchor_pos = a_run - 1          # A 串前一个碱基的 1-based 位置
    anchor = seq[anchor_pos - 1]    # 该锚定碱基本身
    rows_a.append(row(anchor_pos, anchor, anchor + "A", "0/1"))
    rows_b.append(row(a_run, "A", "AA", "0/1"))

    # 3) 同一缺失的两种写法（AT 重复区，删一个 AT）
    rows_a.append(row(at_run - 1, seq[at_run - 2] + "AT", seq[at_run - 2], "0/1"))
    rows_b.append(row(at_run + 1, "TAT", "T", "0/1"))

    # 4) 多 ALT 记录：A 侧 1/2 两个等位各剂量 1；B 侧只报第一个 ALT
    p = find(seq, "GGATCC")
    multi_pos = p  # 参考 G
    rows_a.append(row(multi_pos, "G", "A,T", "1/2"))
    rows_b.append(row(multi_pos, "G", "A", "0/1"))

    # 5) 剂量不符：A 纯合 1/1，B 杂合 0/1
    p = find(seq, "TTTT")
    mm_pos = p + 1  # 参考 T
    rows_a.append(row(mm_pos, "T", "C", "1/1", "42"))
    rows_b.append(row(mm_pos, "T", "C", "0/1", "40"))

    # 6) 仅 A 侧出现
    p = find(seq, "GCGCGCGC")
    a_only_pos = p  # 参考 G
    rows_a.append(row(a_only_pos, "G", "T", "0/1"))

    # 7) 仅 B 侧出现
    b_only_pos = p + 2  # 同区域内另一个 G
    rows_b.append(row(b_only_pos, "G", "A", "0/1"))

    # 8) B 侧一条 0/0 记录：剂量 0，应被完全忽略
    rows_b.append(row(snv_pos, "C", "G", "0/0"))

    (ROOT / "samples" / "caller_a.vcf").write_text(
        HEADER + "\n".join(rows_a) + "\n", encoding="utf-8"
    )
    (ROOT / "samples" / "caller_b.vcf").write_text(
        HEADER + "\n".join(rows_b) + "\n", encoding="utf-8"
    )
    print(f"wrote {len(rows_a)} records -> samples/caller_a.vcf")
    print(f"wrote {len(rows_b)} records -> samples/caller_b.vcf")


if __name__ == "__main__":
    main()
