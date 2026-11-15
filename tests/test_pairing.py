"""批次样本配对核心逻辑测试，仅依赖标准库。

python3 -I -m unittest discover -s tests -t . -p 'test_pairing.py' -v
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.bio import CheckError, parse_fasta  # noqa: E402
from app.pairing import (  # noqa: E402
    MAX_SAMPLES_PER_VCF,
    PairStats,
    parse_batch_vcf,
    pair_batches,
)

# 参考："ACGT" * 10（40 bp），所有 C 位点为 4k+2。
SEQ = "ACGT" * 10
FASTA = f">chrP\n{SEQ}\n".encode()
C_POSITIONS = list(range(2, 41, 4))  # 2, 6, 10, ..., 38
T_POSITIONS = list(range(4, 41, 4))  # 4, 8, ..., 40（ALT 用于不同位点集合）


def refobj():
    return parse_fasta(FASTA)


def vcf(samples, rows):
    header = (
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        "\tFORMAT\t" + "\t".join(samples)
    )
    lines = [header]
    for pos, ref, alt, gts in rows:
        assert len(gts) == len(samples)
        lines.append(
            f"chrP\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(gts)
        )
    return ("\n".join(lines) + "\n").encode()


def parse(data, name="batch"):
    return parse_batch_vcf(data, name, refobj())


def snp_rows(gts_per_sample, n=10, alt="G"):
    """在 n 个 C 位点上生成同一 ALT、各样本指定基因型的记录。"""
    return [
        (C_POSITIONS[i], "C", alt, [g[i] for g in gts_per_sample])
        for i in range(n)
    ]


class ParsingTests(unittest.TestCase):
    def test_multi_sample_parsing_keeps_samples(self):
        batch = parse(vcf(["S1", "S2"], [(2, "C", "G", ["0/0", "1/1"])]))
        self.assertEqual(batch.samples, ["S1", "S2"])
        self.assertEqual(batch.sites[2].gts, (0, 2))

    def test_all_allowed_gt_forms(self):
        batch = parse(
            vcf(["S1", "S2", "S3", "S4"],
                [(2, "C", "G", ["0/0", "0/1", "1/1", "./."])])
        )
        self.assertEqual(batch.sites[2].gts, (0, 1, 2, None))

    def test_illegal_gt_fails_whole_analysis(self):
        for bad in ("0|1", "./1", "0/.", "1|1", "1/2", "2/2", "0", "0/0/1",
                    ".", ""):
            with self.subTest(gt=bad):
                with self.assertRaises(CheckError):
                    parse(vcf(["S1"], [(2, "C", "G", [bad])]))

    def test_multiallelic_record_ignored_but_gt_must_be_allowed(self):
        # 多 ALT 记录不参与配对；但 GT 仍须在允许范围内，1/2 会整次失败。
        with self.assertRaises(CheckError):
            parse(vcf(["S1"], [(2, "C", "G,T", ["1/2"])]))
        batch = parse(vcf(["S1"], [(2, "C", "G,T", ["0/1"])]))
        self.assertEqual(batch.sites, {})

    def test_non_snp_records_ignored_but_gt_validated(self):
        header = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        )
        rows = [
            "chrP\t2\t.\tC\tG\t.\tPASS\t.\tGT\t0/1",      # 双等位 SNP，保留
            "chrP\t6\t.\tC\tCT\t.\tPASS\t.\tGT\t0/1",     # 插入，忽略
            "chrP\t10\t.\tCG\tC\t.\tPASS\t.\tGT\t1/1",    # 缺失，忽略
        ]
        batch = parse((header + "\n".join(rows) + "\n").encode())
        self.assertEqual(sorted(batch.sites), [2])

        # 非法 GT 即使在不参与配对的记录上也会整次失败
        for gt in ("0|1", "1/2", "./1", "0", ""):
            with self.subTest(gt=gt):
                bad = header + f"chrP\t6\t.\tC\tCT\t.\tPASS\t.\tGT\t{gt}\n"
                with self.assertRaises(CheckError):
                    parse(bad.encode())

    def test_duplicate_position_fails(self):
        with self.assertRaises(CheckError) as ctx:
            parse(vcf(["S1"], [(2, "C", "G", ["0/1"]),
                              (2, "C", "T", ["0/1"])]))
        self.assertIn("重复", str(ctx.exception))

    def test_duplicate_position_fails_even_for_ignored_record_shapes(self):
        # 重复判定在“是否参与配对”的筛选之前完成
        with self.assertRaises(CheckError):
            parse(vcf(["S1"], [(2, "C", "G,T", ["0/1"]),
                              (2, "C", "G,T", ["0/1"])]))

    def test_ref_mismatch_fails(self):
        with self.assertRaises(CheckError) as ctx:
            parse(vcf(["S1"], [(3, "A", "G", ["0/1"])]))  # pos3 实为 G
        self.assertIn("REF", str(ctx.exception))

    def test_chrom_mismatch_fails(self):
        data = (
            b"##fileformat=VCFv4.2\n"
            b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            b"chrX\t2\t.\tC\tG\t.\tPASS\t.\tGT\t0/1\n"
        )
        with self.assertRaises(CheckError):
            parse(data)

    def test_sample_count_limits(self):
        ok = vcf([f"S{i}" for i in range(MAX_SAMPLES_PER_VCF)],
                 [(2, "C", "G", ["0/1"] * MAX_SAMPLES_PER_VCF)])
        self.assertEqual(len(parse(ok).samples), MAX_SAMPLES_PER_VCF)

        too_many = vcf([f"S{i}" for i in range(MAX_SAMPLES_PER_VCF + 1)],
                       [(2, "C", "G",
                         ["0/1"] * (MAX_SAMPLES_PER_VCF + 1))])
        with self.assertRaises(CheckError):
            parse(too_many)

    def test_zero_samples_fails(self):
        data = (
            b"##fileformat=VCFv4.2\n"
            b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )
        with self.assertRaises(CheckError):
            parse(data)

    def test_duplicate_sample_name_fails(self):
        with self.assertRaises(CheckError):
            parse(vcf(["S1", "S1"], [(2, "C", "G", ["0/1", "1/1"])]))

    def test_non_snp_records_ignored_but_still_validated(self):
        header = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
        )
        rows = [
            "chrP\t2\t.\tC\tG\t.\tPASS\t.\tGT\t0/1\t0/1",      # 保留
            "chrP\t6\t.\tC\tG,T\t.\tPASS\t.\tGT\t0/1\t0/1",    # 多 ALT，忽略
            "chrP\t10\t.\tCG\tC\t.\tPASS\t.\tGT\t0/1\t0/1",    # 缺失，忽略
            "chrP\t14\t.\tC\tG,C\t.\tPASS\t.\tGT\t0/1\t0/1",   # 重复 ALT，忽略
        ]
        batch = parse((header + "\n".join(rows) + "\n").encode())
        self.assertEqual(sorted(batch.sites), [2])

    def test_non_snp_record_with_bad_ref_still_fails(self):
        header = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        )
        bad = header + "chrP\t10\t.\tTT\tT\t.\tPASS\t.\tGT\t0/1\n"
        with self.assertRaises(CheckError):
            parse(bad.encode())

    def test_row_column_count_must_match_header(self):
        header = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        )
        with self.assertRaises(CheckError):
            parse((header + "chrP\t2\t.\tC\tG\t.\tPASS\t.\tGT\n").encode())

    def test_gt_not_first_in_format_is_accepted(self):
        header = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        )
        data = header + "chrP\t2\t.\tC\tG\t.\tPASS\t.\tDP:GT\t30:0/1\n"
        self.assertEqual(parse(data.encode()).sites[2].gts, (1,))

    def test_format_without_gt_fails(self):
        header = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        )
        with self.assertRaises(CheckError):
            parse((header + "chrP\t2\t.\tC\tG\t.\tPASS\t.\tDP\t30\n").encode())


class PairingTests(unittest.TestCase):
    def test_perfect_match(self):
        a = parse(vcf(["A1"], snp_rows([["0/1"] * 10])))
        b = parse(vcf(["B1"], snp_rows([["0/1"] * 10])))
        result = pair_batches(refobj(), a, b)
        self.assertEqual(len(result.shared_sites), 10)
        self.assertEqual(result.max_pairs, 1)
        self.assertEqual(result.min_inconsistent, 0)
        (stat,) = result.pair_stats
        self.assertEqual((stat.comparable, stat.inconsistent), (10, 0))
        self.assertTrue(stat.eligible)
        self.assertTrue(stat.in_optimal)
        counterpart = result.possibilities_a[0].counterparts[0]
        self.assertEqual(counterpart.sample, "B1")
        self.assertEqual(counterpart.comparable, 10)
        self.assertEqual(counterpart.inconsistent, 0)
        self.assertFalse(result.possibilities_a[0].unpaired_possible)

    def test_eight_comparable_is_threshold(self):
        a = parse(vcf(["A1"], snp_rows([["0/1"] * 8], n=8)))
        b = parse(vcf(["B1"], snp_rows([["0/1"] * 8], n=8)))
        self.assertEqual(pair_batches(refobj(), a, b).max_pairs, 1)

        a7 = parse(vcf(["A1"], snp_rows([["0/1"] * 7], n=7)))
        b7 = parse(vcf(["B1"], snp_rows([["0/1"] * 7], n=7)))
        result7 = pair_batches(refobj(), a7, b7)
        self.assertEqual(result7.max_pairs, 0)
        self.assertFalse(result7.pair_stats[0].eligible)

    def test_one_inconsistency_allowed_two_not(self):
        a = parse(vcf(["A1"], snp_rows([["0/1"] * 10])))
        b1 = parse(vcf(["B1"], snp_rows([["0/1"] * 9 + ["0/0"]])))
        result = pair_batches(refobj(), a, b1)
        self.assertEqual(result.min_inconsistent, 1)
        self.assertEqual(result.max_pairs, 1)

        b2 = parse(vcf(["B2"], snp_rows([["0/1"] * 8 + ["0/0"] * 2])))
        self.assertEqual(pair_batches(refobj(), a, b2).max_pairs, 0)

    def test_missing_calls_not_comparable(self):
        # 两侧共有 10 个位点，其中 2 个 ./. 不计，剩 8 个可比且一致
        a = parse(vcf(["A1"], snp_rows([["0/1"] * 10])))
        b = parse(vcf(["B1"], snp_rows([["0/1"] * 8 + ["./."] * 2])))
        stat = pair_batches(refobj(), a, b).pair_stats[0]
        self.assertEqual((stat.comparable, stat.inconsistent), (8, 0))
        self.assertTrue(stat.eligible)

        # 1 个缺失 + 1 处不一致：可比 9，不一致 1，仍可配对
        b2 = parse(vcf(["B2"], snp_rows([["0/1"] * 8 + ["0/0", "./."]])))
        stat2 = pair_batches(refobj(), a, b2).pair_stats[0]
        self.assertEqual((stat2.comparable, stat2.inconsistent), (9, 1))
        self.assertTrue(stat2.eligible)

    def test_zero_dose_calls_are_comparable(self):
        a = parse(vcf(["A1"], snp_rows([["0/0"] * 8], n=8)))
        b = parse(vcf(["B1"], snp_rows([["0/0"] * 8], n=8)))
        stat = pair_batches(refobj(), a, b).pair_stats[0]
        self.assertEqual((stat.comparable, stat.inconsistent), (8, 0))
        self.assertTrue(stat.eligible)

    def test_sites_must_be_in_both_files_with_same_alt(self):
        # A 报告 10 个位点、B 只报告 6 个 → 共有 6 个，不可配对
        a = parse(vcf(["A1"], snp_rows([["0/1"] * 10])))
        b = parse(vcf(["B1"], snp_rows([["0/1"] * 6], n=6)))
        result = pair_batches(refobj(), a, b)
        self.assertEqual(len(result.shared_sites), 6)
        self.assertEqual(result.max_pairs, 0)

        # 同一位置但 ALT 不同 → 不是共有位点
        b_alt = parse(vcf(["B1"], snp_rows([["0/1"] * 10], alt="T")))
        result_alt = pair_batches(refobj(), a, b_alt)
        self.assertEqual(result_alt.shared_sites, [])
        self.assertEqual(result_alt.max_pairs, 0)

    def test_max_pairs_priority_over_total_inconsistency(self):
        # a1、a2 基因型相同；b1 与它们完全一致，b2 各有 1 处不一致。
        # 任意两对方案必须用到 b2，累计不一致 = 1；
        # 单对方案 a1-b1 累计不一致 = 0，但配对数更少。
        # 全局目标要求先最大化配对数，所以最优是 2 对、累计 1。
        a = parse(vcf(["a1", "a2"],
                      snp_rows([["0/1"] * 10, ["0/1"] * 10])))
        b = parse(vcf(["b1", "b2"],
                      snp_rows([["0/1"] * 10,
                                ["0/1"] * 9 + ["0/0"]])))
        result = pair_batches(refobj(), a, b)
        self.assertEqual(result.max_pairs, 2)
        self.assertEqual(result.min_inconsistent, 1)
        self.assertEqual(len(result.schemes), 2)
        for scheme in result.schemes:
            self.assertEqual(len(scheme.pairs), 2)
            self.assertEqual(scheme.inconsistent, 1)

        stats = {(s.a, s.b): s for s in result.pair_stats}
        self.assertEqual(stats[("a1", "b1")].inconsistent, 0)
        self.assertEqual(stats[("a1", "b2")].inconsistent, 1)
        self.assertEqual(stats[("a2", "b2")].inconsistent, 1)

    def test_ties_are_all_reported(self):
        # 两个 B 样本都能和 A1 完美配对 → 两个同分方案
        a = parse(vcf(["A1"], snp_rows([["0/1"] * 10])))
        b = parse(vcf(["B1", "B2"],
                      snp_rows([["0/1"] * 10, ["0/1"] * 10])))
        result = pair_batches(refobj(), a, b)
        self.assertEqual(result.max_pairs, 1)
        self.assertEqual(len(result.schemes), 2)
        names = {c.sample for c in result.possibilities_a[0].counterparts}
        self.assertEqual(names, {"B1", "B2"})
        self.assertFalse(result.possibilities_a[0].unpaired_possible)

    def test_sample_may_be_unpaired_in_some_optimal_scheme(self):
        # B 只有 1 个样本，两个 A 样本都能与它完美配对 → 两个方案各配一个，
        # 因此两个 A 样本都可能未配对；B 样本在任一方案中都配对。
        a = parse(vcf(["A1", "A2"],
                      snp_rows([["0/1"] * 10, ["0/1"] * 10])))
        b = parse(vcf(["B1"], snp_rows([["0/1"] * 10])))
        result = pair_batches(refobj(), a, b)
        self.assertEqual(result.max_pairs, 1)
        self.assertEqual(len(result.schemes), 2)
        for poss in result.possibilities_a:
            self.assertTrue(poss.unpaired_possible)
            self.assertEqual([c.sample for c in poss.counterparts], ["B1"])
        self.assertFalse(result.possibilities_b[0].unpaired_possible)

    def test_no_shared_site_means_nothing_paired(self):
        a = parse(vcf(["A1"], snp_rows([["0/1"] * 10])))
        b = parse(vcf(["B1"], []))
        result = pair_batches(refobj(), a, b)
        self.assertEqual(result.shared_sites, [])
        self.assertEqual(result.max_pairs, 0)
        self.assertEqual(result.schemes[0].pairs, [])
        self.assertTrue(result.possibilities_a[0].unpaired_possible)
        self.assertTrue(result.possibilities_b[0].unpaired_possible)
        self.assertEqual(result.possibilities_a[0].counterparts, [])

    def test_ineligible_pairs_still_listed_with_stats(self):
        a = parse(vcf(["A1"], snp_rows([["0/1"] * 7], n=7)))
        b = parse(vcf(["B1"], snp_rows([["0/1"] * 7], n=7)))
        result = pair_batches(refobj(), a, b)
        stat = result.pair_stats[0]
        self.assertEqual((stat.comparable, stat.inconsistent), (7, 0))
        self.assertFalse(stat.eligible)
        self.assertFalse(stat.in_optimal)
        self.assertTrue(result.possibilities_a[0].unpaired_possible)

    def test_matching_is_independent_of_sample_order(self):
        a = parse(vcf(["A1", "A2"],
                      snp_rows([["0/1"] * 10, ["0/0"] * 10])))
        b = parse(vcf(["B2", "B1"],
                      snp_rows([["0/1"] * 10, ["0/0"] * 10])))
        result = pair_batches(refobj(), a, b)
        self.assertEqual(result.max_pairs, 2)
        self.assertEqual(
            sorted(result.schemes[0].pairs), [("A1", "B2"), ("A2", "B1")]
        )

    def test_counterpart_stats_are_per_pair(self):
        a = parse(vcf(["A1"], snp_rows([["0/1"] * 10])))
        b = parse(vcf(["B1", "B2"],
                      snp_rows([["0/1"] * 10, ["0/1"] * 10])))
        result = pair_batches(refobj(), a, b)
        self.assertEqual(len(result.possibilities_a[0].counterparts), 2)
        for c in result.possibilities_a[0].counterparts:
            self.assertEqual((c.comparable, c.inconsistent), (10, 0))

    def test_eligible_pair_not_in_any_optimal_scheme_is_not_listed(self):
        # 20 个位点，两位点 s0/s1 是“区分位点”。
        # a1 与 a2 仅在 s0/s1 不同；b1 在 s0/s1 缺失，其余与 a1 相同；
        # b2 在 s0/s1 与 a1 相同，但在另外 12 个位点缺失，只与 a1/a2 共有 8 个可比位点。
        #   a1-b1：18 可比、0 不一致 → 可配对；
        #   a2-b1：18 可比、0 不一致 → 可配对；
        #   a1-b2：8 可比、0 不一致 → 可配对；
        #   a2-b2：8 可比、2 不一致 → 不可配对。
        # 因此最大配对数 2 只能由 a1-b2 + a2-b1 达成；a1-b1 虽可配对，
        # 却不出现在任何最优方案中，不应被列为“可能对应对象”。
        seq = "ACGT" * 20
        ref = parse_fasta(f">chrP\n{seq}\n".encode())
        positions = list(range(2, len(seq) + 1, 4))  # 20 个 C 位点
        n = len(positions)
        base = [i % 3 for i in range(n)]
        s0, s1 = 0, 1
        missing_b1 = {s0, s1}
        missing_b2 = set(range(2, 14))

        def gt(dose):
            return {0: "0/0", 1: "0/1", 2: "1/1"}[dose]

        a_rows = []
        for i in range(n):
            d1 = base[i]
            d2 = (base[i] + 1) % 3 if i in (s0, s1) else base[i]
            a_rows.append((positions[i], "C", "G", [gt(d1), gt(d2)]))
        b_rows = []
        for i in range(n):
            v1 = "./." if i in missing_b1 else gt(base[i])
            v2 = "./." if i in missing_b2 else gt(base[i])
            b_rows.append((positions[i], "C", "G", [v1, v2]))

        batch_a = parse_batch_vcf(
            vcf(["A1", "A2"], a_rows), "A", ref
        )
        batch_b = parse_batch_vcf(vcf(["B1", "B2"], b_rows), "B", ref)
        result = pair_batches(ref, batch_a, batch_b)
        stats = {(s.a, s.b): s for s in result.pair_stats}

        self.assertTrue(stats[("A1", "B1")].eligible)
        self.assertTrue(stats[("A2", "B1")].eligible)
        self.assertTrue(stats[("A1", "B2")].eligible)
        self.assertFalse(stats[("A2", "B2")].eligible)
        self.assertEqual(result.max_pairs, 2)
        self.assertEqual(
            sorted(result.schemes[0].pairs), [("A1", "B2"), ("A2", "B1")]
        )
        self.assertFalse(stats[("A1", "B1")].in_optimal)
        # A1 的可能对象只有 B2，B1 不出现
        names = {c.sample for c in result.possibilities_a[0].counterparts}
        self.assertEqual(names, {"B2"})

    def test_six_by_six_is_fast_and_complete(self):
        # 6×6 全部可配对且完全一致 → 方案数 = 6! = 720
        n = 6
        a = parse(vcf([f"A{i}" for i in range(n)],
                      snp_rows([["0/1"] * 10] * n)))
        b = parse(vcf([f"B{i}" for i in range(n)],
                      snp_rows([["0/1"] * 10] * n)))
        result = pair_batches(refobj(), a, b)
        self.assertEqual(result.max_pairs, 6)
        self.assertEqual(result.min_inconsistent, 0)
        self.assertEqual(len(result.schemes), 720)
        self.assertEqual(result.max_paired_samples, 12)
        # 每个样本都能与任一对方样本配对
        for poss in result.possibilities_a:
            self.assertEqual(len(poss.counterparts), 6)
            self.assertFalse(poss.unpaired_possible)

    def test_fully_connected_six_by_six_enumeration_is_bounded(self):
        # 完全图（任意一对都可配对）时穷举所有匹配，上界为 13327；
        # 这里确认在可接受时间内完成且候选对象覆盖全部对方样本。
        n = 6
        a = parse(vcf([f"A{i}" for i in range(n)],
                      snp_rows([["0/1"] * 10] * n)))
        b = parse(vcf([f"B{i}" for i in range(n)],
                      snp_rows([["0/1"] * 10] * n)))
        result = pair_batches(refobj(), a, b)
        for poss in result.possibilities_b:
            self.assertEqual(len(poss.counterparts), 6)
        self.assertEqual(len(result.pair_stats), 36)
        self.assertEqual(result.max_pairs, min(n, n))

    def test_enumeration_matches_independent_bruteforce(self):
        """用独立的暴力穷举交叉验证全局最优配对与全部同分方案。"""
        import itertools
        import random

        from app.pairing import _enumerate_optimal

        def brute(n_a, n_b, edges):
            best_pairs, best_cost, schemes = -1, -1, []
            for k in range(0, min(n_a, n_b) + 1):
                for asub in itertools.combinations(range(n_a), k):
                    for bsub in itertools.permutations(range(n_b), k):
                        pairs = tuple(zip(asub, bsub))
                        if any(p not in edges for p in pairs):
                            continue
                        cost = sum(edges[p].inconsistent for p in pairs)
                        if k > best_pairs or (
                            k == best_pairs and cost < best_cost
                        ):
                            best_pairs, best_cost, schemes = k, cost, [pairs]
                        elif k == best_pairs and cost == best_cost:
                            schemes.append(pairs)
            return best_pairs, best_cost, schemes

        rng = random.Random(20261115)
        for _ in range(200):
            n_a = rng.randint(1, 3)
            n_b = rng.randint(1, 3)
            edges = {
                (i, j): PairStats(
                    str(i), str(j), comparable=8,
                    inconsistent=rng.randint(0, 1), eligible=True,
                )
                for i in range(n_a)
                for j in range(n_b)
                if rng.random() < 0.7
            }
            got = _enumerate_optimal(n_a, n_b, edges)
            want = brute(n_a, n_b, edges)
            self.assertEqual(got[:2], want[:2])
            self.assertEqual(
                sorted(map(sorted, got[2])), sorted(map(sorted, want[2]))
            )


if __name__ == "__main__":
    unittest.main()
