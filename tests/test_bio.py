"""核心逻辑测试，仅依赖标准库：python -I -m unittest -v"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.bio import (  # noqa: E402
    CheckError,
    compare,
    normalize,
    parse_fasta,
    parse_vcf,
)

# 1-based: 1A 2C 3G 4T | 5A 6A 7A 8A | 9C 10G 11T 12T
SEQ = "ACGTAAAACGTT"
FASTA = f">chrT\n{SEQ}\n".encode()


def vcf(records, sample="S1", extra_header_cols=""):
    lines = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample
        + extra_header_cols,
    ]
    lines.extend(records)
    return "\n".join(lines) + "\n"


def rec(pos, ref, alt, gt, chrom="chrT", fmt="GT"):
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\t{fmt}\t{gt}"


def parse(data):
    refobj = parse_fasta(FASTA)
    calls, _ = parse_vcf(data.encode(), "test.vcf", refobj)
    return refobj, calls


class NormalizeTests(unittest.TestCase):
    def test_insertion_in_homopolymer_different_writing(self):
        # 在 5-8 的 A 串中插入一个 A：pos4 T>TA 与 pos5 A>AA 应归一为同一变异
        a = normalize(4, "T", "TA", SEQ)
        b = normalize(5, "A", "AA", SEQ)
        self.assertEqual(a, b)
        self.assertEqual(a, (4, "T", "TA"))

    def test_deletion_in_homopolymer_different_writing(self):
        a = normalize(4, "TA", "T", SEQ)   # 删 pos5 的 A
        b = normalize(5, "AA", "A", SEQ)   # 删 pos6 的 A，同属 A 串
        self.assertEqual(a, b)
        self.assertEqual(a, (4, "TA", "T"))

    def test_snv_not_moved(self):
        self.assertEqual(normalize(2, "C", "G", SEQ), (2, "C", "G"))

    def test_different_indel_lengths_stay_distinct(self):
        one = normalize(4, "TA", "T", SEQ)
        two = normalize(4, "TAA", "T", SEQ)
        self.assertNotEqual(one, two)

    def test_repeat_run_tt_at_end_deletion(self):
        # 末端 TT 串（pos11-12）中删一个 T：两种写法归一到最左锚点
        a = normalize(10, "GT", "G", SEQ)    # 删 pos11 的 T
        b = normalize(11, "TT", "T", SEQ)    # 删 pos12 的 T
        self.assertEqual(a, b)
        self.assertEqual(a[0], 10)

    def test_insertion_shift_stops_when_it_would_become_mnv(self):
        # pos12 后插入 G：再左移就不再是纯插入，规范结果保持在原位
        self.assertEqual(normalize(12, "T", "TG", SEQ), (12, "T", "TG"))


class ParseAndCompareTests(unittest.TestCase):
    def test_snv_consistent_with_line_numbers(self):
        _, ca = parse(vcf([rec(2, "C", "G", "0/1")]))
        _, cb = parse(vcf([rec(2, "C", "G", "0/1")]))
        refobj = parse_fasta(FASTA)
        groups = compare(refobj, ca, cb)
        self.assertEqual(len(groups["consistent"]), 1)
        entry = groups["consistent"][0]
        self.assertEqual(entry["left"].dosage, 1)
        self.assertEqual(entry["left"].line_numbers, [3])
        self.assertEqual(entry["right"].line_numbers, [3])

    def test_indel_written_differently_treated_as_same_variant(self):
        _, ca = parse(vcf([rec(4, "T", "TA", "0/1")]))
        _, cb = parse(vcf([rec(5, "A", "AA", "0/1")]))
        refobj = parse_fasta(FASTA)
        groups = compare(refobj, ca, cb)
        self.assertEqual(len(groups["consistent"]), 1)
        self.assertEqual(groups["consistent"][0]["variant"].label, "4:TinsTA")
        # 原始行号与原始写法都要保留
        self.assertEqual(groups["consistent"][0]["right"].raw_descriptions, ["5:A>AA"])

    def test_multi_alt_dosages_per_allele(self):
        _, calls = parse(vcf([rec(2, "C", "G,T", "1/2")]))
        keys = sorted(calls, key=lambda v: v.alt)
        self.assertEqual([k.alt for k in keys], ["G", "T"])
        self.assertTrue(all(c.dosage == 1 for c in calls.values()))

    def test_multi_alt_homozygous_alt_and_ignore_zero_dose(self):
        _, calls = parse(vcf([rec(2, "C", "G,T", "1/1")]))
        self.assertEqual(len(calls), 1)  # T 剂量 0 被忽略
        (only,) = calls.values()
        self.assertEqual(only.dosage, 2)

    def test_dosage_mismatch_group(self):
        _, ca = parse(vcf([rec(2, "C", "G", "1/1")]))
        _, cb = parse(vcf([rec(2, "C", "G", "0/1")]))
        groups = compare(parse_fasta(FASTA), ca, cb)
        self.assertEqual(len(groups["dosage_mismatch"]), 1)
        self.assertEqual(groups["consistent"], [])

    def test_single_side_groups(self):
        _, ca = parse(vcf([rec(2, "C", "G", "0/1")]))
        _, cb = parse(vcf([rec(3, "G", "T", "0/1")]))
        groups = compare(parse_fasta(FASTA), ca, cb)
        self.assertEqual(len(groups["left_only"]), 1)
        self.assertEqual(len(groups["right_only"]), 1)

    def test_two_insertion_alts_of_different_length_are_distinct(self):
        # 同一记录的两个 ALT：插 1 个 A 与插 2 个 A，归一后是两个不同变异，
        # 各剂量 1，行号各自只记一次
        data = vcf(["chrT\t5\t.\tA\tAA,AAA\t.\t.\t.\tGT\t1/2"])
        _, calls = parse(data)
        labels = {v.label: c for v, c in calls.items()}
        self.assertEqual(labels["4:TinsTA"].dosage, 1)
        self.assertEqual(labels["4:TinsTA"].line_numbers, [3])
        ins2 = [v for v in calls if v.kind == "INS" and len(v.alt) - 1 == 2]
        self.assertEqual(len(ins2), 1)
        self.assertEqual(calls[ins2[0]].dosage, 1)

    def test_same_variant_two_heterozygous_records_merge_to_dosage_two(self):
        data = vcf([rec(4, "T", "TA", "0/1"), rec(5, "A", "AA", "0/1")])
        _, calls = parse(data)
        self.assertEqual(len(calls), 1)
        (only,) = calls.values()
        self.assertEqual(only.dosage, 2)
        self.assertEqual(only.line_numbers, [3, 4])

    def test_ref_mismatch_fails_whole_check(self):
        with self.assertRaises(CheckError) as ctx:
            parse(vcf([rec(2, "G", "T", "0/1")]))  # pos2 参考实为 C
        self.assertIn("REF", str(ctx.exception))

    def test_contradictory_dosage_in_same_file_fails(self):
        data = vcf([rec(2, "C", "G", "0/1"), rec(2, "C", "G", "1/1")])
        with self.assertRaises(CheckError) as ctx:
            parse(data)
        self.assertIn("矛盾", str(ctx.exception))

    def test_homozygous_plus_extra_record_contradiction(self):
        data = vcf([rec(2, "C", "G", "1/1"), rec(2, "C", "G", "0/1")])
        with self.assertRaises(CheckError):
            parse(data)

    def test_phased_gt_rejected(self):
        with self.assertRaises(CheckError):
            parse(vcf([rec(2, "C", "G", "0|1")]))

    def test_missing_gt_rejected(self):
        for bad in ("./1", "0/.", "./."):
            with self.subTest(gt=bad):
                with self.assertRaises(CheckError):
                    parse(vcf([rec(2, "C", "G", bad)]))

    def test_wrong_ploidy_rejected(self):
        for bad in ("0", "0/1/1"):
            with self.subTest(gt=bad):
                with self.assertRaises(CheckError):
                    parse(vcf([rec(2, "C", "G", bad)]))

    def test_gt_index_out_of_range(self):
        with self.assertRaises(CheckError):
            parse(vcf([rec(2, "C", "G,T", "0/3")]))

    def test_insertion_too_long_rejected(self):
        with self.assertRaises(CheckError):
            parse(vcf([rec(1, "A", "A" + "C" * 21, "0/1")]))  # 插入 21 bp

    def test_deletion_too_long_rejected(self):
        # 构造 21 bp 缺失：需要足够参考长度，直接用 classify 路径验证
        long_ref = "A" + "C" * 21
        seq2 = long_ref + "GG"
        fasta2 = b">chrL\n" + seq2.encode()
        refobj = parse_fasta(fasta2)
        data = vcf([rec(1, long_ref, "A", "0/1")])
        with self.assertRaises(CheckError):
            parse_vcf(data.encode(), "x.vcf", refobj)

    def test_mnv_rejected(self):
        with self.assertRaises(CheckError):
            parse(vcf([rec(1, "AC", "GT", "0/1")]))

    def test_symbolic_and_non_acgt_alleles_rejected(self):
        for alt in (".", "G*", "<DEL>", "g"):
            with self.subTest(alt=alt):
                with self.assertRaises(CheckError):
                    parse(vcf([rec(2, "C", alt, "0/1")]))

    def test_chrom_mismatch_rejected(self):
        with self.assertRaises(CheckError):
            parse(vcf([rec(2, "C", "G", "0/1", chrom="chrOther")]))

    def test_pos_out_of_range_rejected(self):
        with self.assertRaises(CheckError):
            parse(vcf([rec(99, "C", "G", "0/1")]))

    def test_multi_sample_vcf_rejected(self):
        data = vcf([rec(2, "C", "G", "0/1")], extra_header_cols="\tS2")
        with self.assertRaises(CheckError):
            parse(data)

    def test_missing_header_rejected(self):
        with self.assertRaises(CheckError):
            parse(rec(2, "C", "G", "0/1") + "\n")

    def test_format_not_leading_gt_rejected(self):
        data = vcf([rec(2, "C", "G", "30", fmt="DP:GT")])
        with self.assertRaises(CheckError):
            parse(data)

    def test_format_with_gt_plus_other_fields_ok(self):
        data = vcf(["chrT\t2\t.\tC\tG\t.\t.\t.\tGT:DP\t0/1:30"])
        _, calls = parse(data)
        self.assertEqual(len(calls), 1)

    def test_ref_only_gt_produces_no_variant(self):
        _, calls = parse(vcf([rec(2, "C", "G", "0/0")]))
        self.assertEqual(calls, {})

    def test_comment_and_blank_lines_ignored(self):
        data = vcf(["##source=unit-test", "", rec(2, "C", "G", "0/1")])
        _, calls = parse(data)
        self.assertEqual(len(calls), 1)
        (only,) = calls.values()
        self.assertEqual(only.line_numbers, [5])  # 含两条 ## 头与空行的物理行号

    def test_crlf_line_endings_keep_physical_line_numbers(self):
        data = ("##fileformat=VCFv4.2\r\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\r\n"
                + rec(2, "C", "G", "0/1") + "\r\n")
        _, calls = parse(data)
        (only,) = calls.values()
        self.assertEqual(only.line_numbers, [3])

    def test_missing_gt_in_multi_alt_rejected(self):
        with self.assertRaises(CheckError):
            parse(vcf([rec(2, "C", "G,T", "./2")]))

    def test_phased_multi_alt_rejected(self):
        with self.assertRaises(CheckError):
            parse(vcf([rec(2, "C", "G,T", "1|2")]))


class FastaTests(unittest.TestCase):
    def test_bad_base_fails(self):
        with self.assertRaises(CheckError):
            parse_fasta(b">c\nACGN\n")

    def test_missing_header_fails(self):
        with self.assertRaises(CheckError):
            parse_fasta(b"ACGT\n")

    def test_two_records_fails(self):
        with self.assertRaises(CheckError):
            parse_fasta(b">a\nAC\n>b\nGT\n")

    def test_wrapped_sequence_ok(self):
        refobj = parse_fasta(b">chrT\nACGT\nAAAA\nCGTT\n")
        self.assertEqual(refobj.sequence, SEQ)


if __name__ == "__main__":
    unittest.main()
