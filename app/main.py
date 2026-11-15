"""FastAPI 入口：上传参考 FASTA + 两份单样本 VCF，展示归一化核对结果。"""

from pathlib import Path

from fastapi import FastAPI, File, Form, Request, UploadFile
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

app = FastAPI(
    title="双 VCF 变异核对",
    description="上传单条染色体参考 FASTA 与两份单样本 VCF，按归一化等位变异核对差异。",
)

templates = Jinja2Templates(directory=str(Path(__file__).parent / "templates"))


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
