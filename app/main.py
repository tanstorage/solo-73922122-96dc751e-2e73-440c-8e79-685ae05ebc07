"""FastAPI 入口：

1. ``/`` 上传参考 FASTA + 两份单样本 VCF，核对归一化变异差异（原功能）；
2. ``/pair`` 上传同染色体参考 FASTA + 两份多样本 VCF（各 1–6 样本），
   用共有双等位 SNP 的基因型核对样本身份、给出全部最优配对方案。
"""

from pathlib import Path

from fastapi import FastAPI, File, Request, UploadFile
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.templating import Jinja2Templates

from .bio import (
    MAX_FASTA_BYTES,
    MAX_VCF_BYTES,
    CheckError,
    compare,
    parse_fasta,
    parse_vcf,
)
from .pairing import (
    MAX_INCONSISTENT_SITES,
    MAX_SAMPLES_PER_VCF,
    MIN_COMPARABLE_SITES,
    PairingResult,
    parse_batch_vcf,
    pair_batches,
)

app = FastAPI(
    title="VCF 变异核对与批次样本配对",
    description=(
        "上传单条染色体参考 FASTA 与 VCF，既可以核对两份单样本 VCF 的变异差异，"
        "也可以对两份多样本 VCF 做基于共有 SNP 基因型的样本配对。"
    ),
)

templates = Jinja2Templates(directory=str(Path(__file__).parent / "templates"))

_PAIRING_CONTEXT = {
    "fasta_limit_kib": MAX_FASTA_BYTES // 1024,
    "vcf_limit_kib": MAX_VCF_BYTES // 1024,
    "max_samples": MAX_SAMPLES_PER_VCF,
    "min_comparable": MIN_COMPARABLE_SITES,
    "max_inconsistent": MAX_INCONSISTENT_SITES,
}


def _read_limited(file: UploadFile, limit: int, label: str) -> bytes:
    data = file.file.read(limit + 1)
    if len(data) > limit:
        raise CheckError(
            f"{label} 超过大小上限 {limit} 字节"
            f"（FASTA 上限 {MAX_FASTA_BYTES // 1024} KiB，VCF 上限 {MAX_VCF_BYTES // 1024} KiB）。"
        )
    if not data:
        raise CheckError(f"{label} 为空文件。")
    return data


def _run_compare(fasta_bytes: bytes, vcf_a_bytes: bytes, vcf_b_bytes: bytes):
    refobj = parse_fasta(fasta_bytes)
    calls_a, sample_a = parse_vcf(vcf_a_bytes, "VCF A", refobj)
    calls_b, sample_b = parse_vcf(vcf_b_bytes, "VCF B", refobj)
    groups = compare(refobj, calls_a, calls_b)
    return refobj, sample_a, sample_b, groups


@app.get("/", response_class=HTMLResponse)
def index(request: Request):
    return templates.TemplateResponse(
        request,
        "index.html",
        {
            "fasta_limit_kib": MAX_FASTA_BYTES // 1024,
            "vcf_limit_kib": MAX_VCF_BYTES // 1024,
            "max_samples": MAX_SAMPLES_PER_VCF,
        },
    )


@app.post("/compare", response_class=HTMLResponse)
def compare_page(
    request: Request,
    fasta: UploadFile = File(...),
    vcf_a: UploadFile = File(...),
    vcf_b: UploadFile = File(...),
):
    try:
        fasta_bytes = _read_limited(fasta, MAX_FASTA_BYTES, "参考 FASTA")
        vcf_a_bytes = _read_limited(vcf_a, MAX_VCF_BYTES, "VCF A")
        vcf_b_bytes = _read_limited(vcf_b, MAX_VCF_BYTES, "VCF B")
        refobj, sample_a, sample_b, groups = _run_compare(
            fasta_bytes, vcf_a_bytes, vcf_b_bytes
        )
    except CheckError as exc:
        # 整次核对失败：只展示错误，不展示任何局部结论
        return templates.TemplateResponse(
            request,
            "index.html",
            {
                "fasta_limit_kib": MAX_FASTA_BYTES // 1024,
                "vcf_limit_kib": MAX_VCF_BYTES // 1024,
                "error": str(exc),
            },
            status_code=422,
        )

    return templates.TemplateResponse(
        request,
        "result.html",
        {
            "ref_name": refobj.name,
            "ref_length": len(refobj.sequence),
            "sample_a": sample_a,
            "sample_b": sample_b,
            "groups": groups,
        },
    )


@app.post("/api/compare")
def compare_api(
    fasta: UploadFile = File(...),
    vcf_a: UploadFile = File(...),
    vcf_b: UploadFile = File(...),
):
    """JSON 接口，便于脚本化调用；失败返回 422 与错误信息。"""
    try:
        fasta_bytes = _read_limited(fasta, MAX_FASTA_BYTES, "参考 FASTA")
        vcf_a_bytes = _read_limited(vcf_a, MAX_VCF_BYTES, "VCF A")
        vcf_b_bytes = _read_limited(vcf_b, MAX_VCF_BYTES, "VCF B")
        refobj, sample_a, sample_b, groups = _run_compare(
            fasta_bytes, vcf_a_bytes, vcf_b_bytes
        )
    except CheckError as exc:
        return JSONResponse(status_code=422, content={"ok": False, "error": str(exc)})

    def dump_entry(entry):
        def dump_call(call):
            if call is None:
                return None
            return {
                "dosage": call.dosage,
                "line_numbers": call.line_numbers,
                "raw": call.raw_descriptions,
            }

        v = entry["variant"]
        return {
            "variant": {
                "pos": v.pos,
                "ref": v.ref,
                "alt": v.alt,
                "kind": v.kind,
                "label": v.label,
            },
            "left": dump_call(entry["left"]),
            "right": dump_call(entry["right"]),
        }

    return {
        "ok": True,
        "reference": {"name": refobj.name, "length": len(refobj.sequence)},
        "samples": {"left": sample_a, "right": sample_b},
        "result": {key: [dump_entry(e) for e in rows] for key, rows in groups.items()},
    }


# --------------------------------------------------------------------------- #
# 批次样本配对
# --------------------------------------------------------------------------- #
def _run_pairing(fasta_bytes: bytes, vcf_a_bytes: bytes, vcf_b_bytes: bytes) -> PairingResult:
    refobj = parse_fasta(fasta_bytes)
    batch_a = parse_batch_vcf(vcf_a_bytes, "VCF A", refobj)
    batch_b = parse_batch_vcf(vcf_b_bytes, "VCF B", refobj)
    return pair_batches(refobj, batch_a, batch_b)


def _pairing_page(request: Request, error: str, status_code: int = 422):
    return templates.TemplateResponse(
        request,
        "pair.html",
        {**_PAIRING_CONTEXT, "error": error},
        status_code=status_code,
    )


@app.get("/pair", response_class=HTMLResponse)
def pair_index(request: Request):
    return templates.TemplateResponse(
        request,
        "pair.html",
        _PAIRING_CONTEXT,
    )


@app.post("/pair", response_class=HTMLResponse)
def pair_page(
    request: Request,
    fasta: UploadFile = File(...),
    vcf_a: UploadFile = File(...),
    vcf_b: UploadFile = File(...),
):
    try:
        fasta_bytes = _read_limited(fasta, MAX_FASTA_BYTES, "参考 FASTA")
        vcf_a_bytes = _read_limited(vcf_a, MAX_VCF_BYTES, "VCF A")
        vcf_b_bytes = _read_limited(vcf_b, MAX_VCF_BYTES, "VCF B")
        result = _run_pairing(fasta_bytes, vcf_a_bytes, vcf_b_bytes)
    except CheckError as exc:
        # 整次分析失败：只展示错误，不展示任何局部映射
        return _pairing_page(request, str(exc))

    return templates.TemplateResponse(
        request,
        "pair_result.html",
        {
            **_PAIRING_CONTEXT,
            "result": result,
            "scheme_limit": 200,
            "schemes": result.schemes[:200],
            "schemes_truncated": len(result.schemes) > 200,
        },
    )


def _dump_pairing(result: PairingResult) -> dict:
    def dump_counterpart(c):
        return {
            "sample": c.sample,
            "comparable": c.comparable,
            "inconsistent": c.inconsistent,
        }

    def dump_possibility(p):
        return {
            "sample": p.sample,
            "batch": p.batch,
            "counterparts": [dump_counterpart(c) for c in p.counterparts],
            "unpaired_possible": p.unpaired_possible,
        }

    return {
        "reference": {
            "name": result.reference_name,
            "length": result.reference_length,
        },
        "samples": {
            "a": result.samples_a,
            "b": result.samples_b,
        },
        "shared_sites": {
            "count": len(result.shared_sites),
            "positions": [s.pos for s in result.shared_sites],
        },
        "pair_stats": [
            {
                "a": s.a,
                "b": s.b,
                "comparable": s.comparable,
                "inconsistent": s.inconsistent,
                "eligible": s.eligible,
                "in_optimal": s.in_optimal,
            }
            for s in result.pair_stats
        ],
        "optimal": {
            "pair_count": result.max_pairs,
            "total_inconsistent": result.min_inconsistent,
            "scheme_count": len(result.schemes),
            "schemes": [
                {
                    "pairs": [
                        {"a": a, "b": b} for a, b in scheme.pairs
                    ],
                    "inconsistent": scheme.inconsistent,
                }
                for scheme in result.schemes
            ],
        },
        "possibilities": {
            "a": [dump_possibility(p) for p in result.possibilities_a],
            "b": [dump_possibility(p) for p in result.possibilities_b],
        },
    }


@app.post("/api/pair")
def pair_api(
    fasta: UploadFile = File(...),
    vcf_a: UploadFile = File(...),
    vcf_b: UploadFile = File(...),
):
    """JSON 接口，便于脚本化调用；失败返回 422 与错误信息。"""
    try:
        fasta_bytes = _read_limited(fasta, MAX_FASTA_BYTES, "参考 FASTA")
        vcf_a_bytes = _read_limited(vcf_a, MAX_VCF_BYTES, "VCF A")
        vcf_b_bytes = _read_limited(vcf_b, MAX_VCF_BYTES, "VCF B")
        result = _run_pairing(fasta_bytes, vcf_a_bytes, vcf_b_bytes)
    except CheckError as exc:
        return JSONResponse(status_code=422, content={"ok": False, "error": str(exc)})

    return {"ok": True, "result": _dump_pairing(result)}
