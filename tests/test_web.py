"""Web 层端到端测试，需要 fastapi/httpx 依赖：
PYTHONPATH=/tmp/pylibs python -I -m unittest -v tests.test_web
"""

import io
import sys
import unittest
from pathlib import Path

# 本地无 pip 环境时，依赖安装在隔离目录；容器内该目录不存在则跳过
for _cand in (Path("/tmp/pylibs"), Path("/workspace/.pylibs")):
    if _cand.is_dir() and str(_cand) not in sys.path:
        sys.path.insert(0, str(_cand))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from fastapi.testclient import TestClient  # noqa: E402

from app.main import app  # noqa: E402
from app.bio import MAX_FASTA_BYTES, MAX_VCF_BYTES  # noqa: E402

client = TestClient(app)

FASTA = b">chrT\nACGTAAAACGTT\n"
HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"


def vcf(rows):
    return (HEADER + "".join(r + "\n" for r in rows)).encode()


def upload(fasta=FASTA, a=None, b=None):
    return client.post(
        "/compare",
        files={
            "fasta": ("ref.fa", fasta, "text/plain"),
            "vcf_a": ("a.vcf", a if a is not None else b"", "text/plain"),
            "vcf_b": ("b.vcf", b if b is not None else b"", "text/plain"),
        },
    )


class WebTests(unittest.TestCase):
    def test_index_renders(self):
        resp = client.get("/")
        self.assertEqual(resp.status_code, 200)
        self.assertIn("双 VCF 变异核对".encode(), resp.content)

    def test_happy_path_groups_and_line_numbers(self):
        a = vcf([
            "chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1",       # SNV 一致
            "chrT\t4\t.\tT\tTA\t.\t.\t.\tGT\t0/1",      # 插入（写法一）
            "chrT\t10\t.\tG\tT\t.\t.\t.\tGT\t1/1",      # 仅 A 有
        ])
        b = vcf([
            "chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1",
            "chrT\t5\t.\tA\tAA\t.\t.\t.\tGT\t0/1",      # 同一插入的另一种写法
            "chrT\t3\t.\tG\tC\t.\t.\t.\tGT\t0/1",       # 仅 B 有
        ])
        resp = upload(a=a, b=b)
        self.assertEqual(resp.status_code, 200)
        body = resp.text
        self.assertIn("两侧一致", body)
        self.assertIn("4:TinsTA", body)          # 归一变异
        self.assertIn("原始行号：4", body)        # 插入记录的物理行号
        self.assertIn("5:A&gt;AA", body)         # B 侧保留原始写法
        self.assertIn("仅单侧出现", body)
        self.assertIn("剂量不符", body)

    def test_ref_mismatch_shows_failure_page_without_partial_results(self):
        a = vcf(["chrT\t2\t.\tG\tT\t.\t.\t.\tGT\t0/1"])  # pos2 参考是 C
        b = vcf(["chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1"])
        resp = upload(a=a, b=b)
        self.assertEqual(resp.status_code, 422)
        self.assertIn("核对失败", resp.text)
        self.assertNotIn("核对概览", resp.text)

    def test_size_limit_enforced(self):
        big = b">chrT\n" + b"A" * (MAX_FASTA_BYTES + 10)
        resp = upload(fasta=big, a=vcf([]), b=vcf([]))
        self.assertEqual(resp.status_code, 422)
        self.assertIn("超过大小上限", resp.text)

        big_vcf = HEADER.encode() + b"#\n" * (MAX_VCF_BYTES // 2)
        resp = upload(a=big_vcf, b=vcf([]))
        self.assertEqual(resp.status_code, 422)

    def test_empty_file_rejected(self):
        resp = upload(a=b"", b=vcf([]))
        self.assertEqual(resp.status_code, 422)
        self.assertIn("为空", resp.text)

    def test_api_json_happy_path(self):
        a = vcf([
            "chrT\t2\t.\tC\tG,T\t.\t.\t.\tGT\t1/2",
            "chrT\t4\t.\tT\tTA\t.\t.\t.\tGT\t1/1",
        ])
        b = vcf([
            "chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1",  # G 剂量不符（1 vs 2）
            "chrT\t4\t.\tT\tTA\t.\t.\t.\tGT\t1/1",
        ])
        resp = client.post(
            "/api/compare",
            files={
                "fasta": ("ref.fa", FASTA, "text/plain"),
                "vcf_a": ("a.vcf", a, "text/plain"),
                "vcf_b": ("b.vcf", b, "text/plain"),
            },
        )
        self.assertEqual(resp.status_code, 200)
        data = resp.json()
        self.assertTrue(data["ok"])
        labels = {e["variant"]["label"] for e in data["result"]["consistent"]}
        self.assertIn("4:TinsTA", labels)
        # 多 ALT 记录 1/2：G、T 各剂量 1；B 侧只有 G 剂量 1 → T 仅 A 有，G 一致
        mismatch = data["result"]["dosage_mismatch"]
        self.assertEqual(mismatch, [])
        left_only = {e["variant"]["label"]: e for e in data["result"]["left_only"]}
        self.assertIn("2:C>T", left_only)
        self.assertEqual(left_only["2:C>T"]["left"]["dosage"], 1)

    def test_api_json_dosage_mismatch(self):
        a = vcf(["chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t1/1"])
        b = vcf(["chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1"])
        resp = client.post(
            "/api/compare",
            files={
                "fasta": ("ref.fa", FASTA, "text/plain"),
                "vcf_a": ("a.vcf", a, "text/plain"),
                "vcf_b": ("b.vcf", b, "text/plain"),
            },
        )
        self.assertEqual(resp.status_code, 200)
        mismatch = resp.json()["result"]["dosage_mismatch"]
        self.assertEqual(len(mismatch), 1)
        self.assertEqual(mismatch[0]["left"]["dosage"], 2)
        self.assertEqual(mismatch[0]["right"]["dosage"], 1)

    def test_api_json_failure(self):
        resp = client.post(
            "/api/compare",
            files={
                "fasta": ("ref.fa", FASTA, "text/plain"),
                "vcf_a": ("a.vcf", vcf(["chrT\t2\t.\tG\tT\t.\t.\t.\tGT\t0/1"]), "text/plain"),
                "vcf_b": ("b.vcf", vcf([]), "text/plain"),
            },
        )
        self.assertEqual(resp.status_code, 422)
        self.assertFalse(resp.json()["ok"])


if __name__ == "__main__":
    unittest.main()
