"""批次配对 Web 层端到端测试，需要 fastapi/httpx 依赖。

PYTHONPATH=/tmp/pylibs python3 -I -m unittest -v tests.test_pairing_web
"""

import sys
import unittest
from pathlib import Path

for _cand in (Path("/tmp/pylibs"), Path("/workspace/.pylibs")):
    if _cand.is_dir() and str(_cand) not in sys.path:
        sys.path.insert(0, str(_cand))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from fastapi.testclient import TestClient  # noqa: E402

from app.main import app  # noqa: E402
from app.bio import MAX_FASTA_BYTES, MAX_VCF_BYTES  # noqa: E402

client = TestClient(app)

FASTA = b">chrP\n" + ("ACGT" * 10).encode() + b"\n"
C_POSITIONS = list(range(2, 41, 4))


def batch_vcf(samples, rows):
    header = (
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        "\tFORMAT\t" + "\t".join(samples)
    )
    lines = [header]
    for pos, alt, gts in rows:
        lines.append(
            f"chrP\t{pos}\t.\tC\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(gts)
        )
    return ("\n".join(lines) + "\n").encode()


def rows(gts_per_sample, n=10, alt="G"):
    return [
        (C_POSITIONS[i], alt, [g[i] for g in gts_per_sample]) for i in range(n)
    ]


def upload(fasta=FASTA, a=b"", b=b"", endpoint="/pair"):
    return client.post(
        endpoint,
        files={
            "fasta": ("ref.fa", fasta, "text/plain"),
            "vcf_a": ("a.vcf", a, "text/plain"),
            "vcf_b": ("b.vcf", b, "text/plain"),
        },
    )


GOOD_A = batch_vcf(["A1", "A2"], rows([["0/1"] * 10, ["0/0"] * 10]))
GOOD_B = batch_vcf(["B1", "B2"], rows([["0/0"] * 10, ["0/1"] * 10]))


class PairWebTests(unittest.TestCase):
    def test_pair_page_renders(self):
        resp = client.get("/pair")
        self.assertEqual(resp.status_code, 200)
        body = resp.text
        self.assertIn("批次样本配对", body)
        self.assertIn("1–6 个样本", body)
        # 原有核对入口仍存在
        self.assertIn("双 VCF 变异核对", body)

    def test_index_links_both_tabs(self):
        resp = client.get("/")
        self.assertIn('href="/pair"', resp.text)
        self.assertIn('href="/"', resp.text)

    def test_pairing_result_shows_possibilities_and_stats(self):
        resp = upload(a=GOOD_A, b=GOOD_B)
        self.assertEqual(resp.status_code, 200)
        body = resp.text
        self.assertIn("配对概览", body)
        self.assertIn("全部最优方案", body)
        self.assertIn("每个样本的可能对应对象", body)
        self.assertIn("可比位点 10", body)
        self.assertIn("不一致位点 0", body)
        # 两个 A 样本各自只能对应一个 B 样本
        self.assertIn("A1", body)
        self.assertIn("B2", body)
        # 候选组合表包含不可配对组合的统计
        self.assertIn("全部候选组合", body)

    def test_unpaired_possible_marker(self):
        # 2 个 A 样本、1 个 B 样本，且都能完美配对 → 两个 A 都可能未配对
        a = batch_vcf(["A1", "A2"], rows([["0/1"] * 10, ["0/1"] * 10]))
        b = batch_vcf(["B1"], rows([["0/1"] * 10]))
        resp = upload(a=a, b=b)
        self.assertEqual(resp.status_code, 200)
        self.assertIn("未配对", resp.text)

    def test_no_eligible_pair_shows_warning(self):
        a = batch_vcf(["A1"], rows([["0/1"] * 7], n=7))
        b = batch_vcf(["B1"], rows([["0/1"] * 7], n=7))
        resp = upload(a=a, b=b)
        self.assertEqual(resp.status_code, 200)
        self.assertIn("没有任何一对样本", resp.text)

    def test_illegal_gt_fails_whole_analysis_without_partial_map(self):
        bad = batch_vcf(["A1"], rows([["0/1"] * 9 + ["1|1"]]))
        resp = upload(a=bad, b=GOOD_B)
        self.assertEqual(resp.status_code, 422)
        self.assertIn("配对分析失败", resp.text)
        self.assertNotIn("配对概览", resp.text)
        self.assertNotIn("每个样本的可能对应对象", resp.text)

    def test_ref_mismatch_fails_whole_analysis(self):
        bad = batch_vcf(["A1"], [(3, "G", ["0/1"])])  # pos3 实为 G，REF=C 不符
        resp = upload(a=bad, b=GOOD_B)
        self.assertEqual(resp.status_code, 422)
        self.assertIn("REF", resp.text)
        self.assertNotIn("配对概览", resp.text)

    def test_duplicate_position_fails_whole_analysis(self):
        a = batch_vcf(["A1"], rows([["0/1"] * 10]) + rows([["0/1"] * 10])[:1])
        # 上面构造出重复 POS=2 的记录
        resp = upload(a=a, b=GOOD_B)
        self.assertEqual(resp.status_code, 422)
        self.assertIn("重复", resp.text)
        self.assertNotIn("配对概览", resp.text)

    def test_size_limit_enforced(self):
        big = b">chrP\n" + b"A" * (MAX_FASTA_BYTES + 10)
        resp = upload(fasta=big, a=GOOD_A, b=GOOD_B)
        self.assertEqual(resp.status_code, 422)
        self.assertIn("超过大小上限", resp.text)

        filler = b"#\n" * (MAX_VCF_BYTES // 2)
        resp = upload(a=GOOD_A + filler, b=GOOD_B)
        self.assertEqual(resp.status_code, 422)
        self.assertIn("超过大小上限", resp.text)

    def test_sample_count_over_limit_fails(self):
        samples = [f"S{i}" for i in range(7)]
        too_many = batch_vcf(samples, [(2, "G", ["0/1"] * 7)])
        resp = upload(a=too_many, b=GOOD_B)
        self.assertEqual(resp.status_code, 422)
        self.assertIn("超过上限", resp.text)

    def test_api_pair_json(self):
        resp = upload(a=GOOD_A, b=GOOD_B, endpoint="/api/pair")
        self.assertEqual(resp.status_code, 200)
        data = resp.json()
        self.assertTrue(data["ok"])
        result = data["result"]
        self.assertEqual(result["samples"], {"a": ["A1", "A2"], "b": ["B1", "B2"]})
        self.assertEqual(result["shared_sites"]["count"], 10)
        self.assertEqual(result["optimal"]["pair_count"], 2)
        self.assertEqual(result["optimal"]["total_inconsistent"], 0)
        self.assertEqual(result["optimal"]["scheme_count"], 1)
        pairs = result["optimal"]["schemes"][0]["pairs"]
        self.assertEqual(
            sorted((p["a"], p["b"]) for p in pairs),
            [("A1", "B2"), ("A2", "B1")],
        )
        poss = {p["sample"]: p for p in result["possibilities"]["a"]}
        self.assertEqual(poss["A1"]["counterparts"], [
            {"sample": "B2", "comparable": 10, "inconsistent": 0}
        ])
        self.assertFalse(poss["A1"]["unpaired_possible"])
        stats = {(s["a"], s["b"]): s for s in result["pair_stats"]}
        self.assertFalse(stats[("A1", "B1")]["eligible"])

    def test_api_pair_failure_json(self):
        bad = batch_vcf(["A1"], [(2, "G", ["1/2"])])
        resp = upload(a=bad, b=GOOD_B, endpoint="/api/pair")
        self.assertEqual(resp.status_code, 422)
        self.assertFalse(resp.json()["ok"])
        self.assertIn("GT", resp.json()["error"])

    def test_existing_compare_endpoint_still_works(self):
        # 原有单样本核对入口不受影响
        header = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        )
        body = (header + "chrP\t2\t.\tC\tG\t.\tPASS\t.\tGT\t0/1\n").encode()
        resp = upload(a=body, b=body, endpoint="/compare")
        self.assertEqual(resp.status_code, 200)
        self.assertIn("核对概览", resp.text)


if __name__ == "__main__":
    unittest.main()
