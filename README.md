# VCF 变异核对与批次样本配对（FastAPI）

同一参考下提供两个入口：

1. **双 VCF 变异核对**（原功能）：核对两套变异识别程序输出的单样本 VCF，按**归一化后的
   等位变异**列出两侧一致、仅单侧出现与剂量不符的结果，并保留每条结论的 VCF 原始行号；
2. **批次样本配对**：怀疑两批测序结果的样本标签互换时，用两份多样本 VCF 中共有的
   双等位 SNP 的基因型核对样本身份，列出全部最优配对方案。

## 一、双 VCF 变异核对

- **输入限制**
  - 单条染色体参考 FASTA，纯 A/C/G/T 序列，不超过 **1 MiB**；
  - 两份单样本 VCF（恰好 10 列、1 个样本），各不超过 **256 KiB**。
- **接受的变异**
  - A/C/G/T 构成的**单碱基替换**（SNV）；
  - 事件长度不超过 **20 bp** 的**纯插入 / 纯缺失**；
  - GT 必须是**非缺失、未分相的二倍体**（`0/1`、`1/1`，多 ALT 可为 `1/2`）。
  - 其他形态（MNV、符号化 ALT、`*`、分相 `|`、缺失 `.`、非二倍体等）一律判失败。
- **归一化**：对照参考序列对 indel 做左对齐（与 `bcftools norm` / `vt normalize` 等价），
  重复区域中位置写法不同的同一插入 / 缺失归为同一等位变异；SNV 位置不变。
- **剂量计算**：多 ALT 记录按 GT 分别计算每个 ALT 的 0/1/2 剂量，剂量为 0 的变异忽略；
  同一文件内同一归一变异的多条记录剂量累加，超过二倍体上限 2 即判定为**互相矛盾**。
- **整次失败**（返回 422，只显示错误、不展示任何局部结论）：
  - 任一记录 REF 与参考序列不符（含坐标越界、CHROM 不一致）；
  - 同一文件中同一归一变异存在互相矛盾的剂量；
  - 以及任何不满足上述输入规则的情形。
- **结果分组**（均按归一变异排序，并标注两侧原始行号与原始 `POS:REF>ALT` 写法）：
  两侧一致 / 仅 VCF A 或 B 出现 / 剂量不符。

### 试用双 VCF 核对

打开 <http://localhost:8000/>，分别选择：

- 参考：`samples/demo.fa`
- VCF A：`samples/caller_a.vcf`
- VCF B：`samples/caller_b.vcf`

演示数据覆盖：SNV 一致、重复区内插入 / 缺失的两种等价写法归一、多 ALT 按 GT 拆分、
`0/0` 剂量 0 被忽略、剂量不符（1/1 vs 0/1）、仅单侧出现。

也可以走 JSON 接口：

```bash
curl -F 'fasta=@samples/demo.fa' \
     -F 'vcf_a=@samples/caller_a.vcf' \
     -F 'vcf_b=@samples/caller_b.vcf' \
     http://localhost:8000/api/compare
```

## 二、批次样本配对

怀疑两批测序结果的样本标签互换时，用两份多样本 VCF 中共有 SNP 的基因型核对身份。

- **输入限制**
  - 同一条染色体的参考 FASTA，纯 A/C/G/T，不超过 **1 MiB**；
  - 两份**多样本 VCF**，每份 **1–6 个样本**，各不超过 **256 KiB**；
    样本名在文件内必须唯一。
- **参与配对的位点**
  - 只使用**两侧都出现**、REF 与参考相符、ALT 相同且为**双等位单碱基**的位点；
  - indel、MNV、多 ALT、符号化 ALT 等记录不参与配对，但其 GT 仍须为
    `0/0`、`0/1`、`1/1` 或 `./.`（多 ALT 记录中的 `1/2` 等写法会判失败）；
  - **重复位点**（同一文件内同一 POS 出现多次）使整次分析失败。
- **GT 接受范围**：所有记录（含不参与配对的记录）的 GT 都必须是未分相的
  `0/0`、`0/1`、`1/1` 或缺失调用 `./.`；`./.` 视为缺失调用，**不计入可比位点**；
  其他写法（分相 `|`、部分缺失 `./1`、多等位 `1/2`、非二倍体、空值等）使整次分析失败。
- **可配对条件**：两个样本**至少 8 个可比位点**且不一致位点**不超过 1 处**。
- **全局目标**：先**最大化配对数**（配成最多对），配对对数相同时再**最小化所有配对的
  累计不一致位点数**；求出的所有同分方案都视为最优。
- **结果展示**
  - 每个原批次样本在**所有最优方案**中的可能对应对象（及“未配对”可能性），
    每对列出**可比位点数与不一致位点数**；
  - 全部最优方案的逐对明细，以及全部 A×B 候选组合（含不可配对者）的统计；
  - 任一文件出现重复位点、非法 GT 或 REF 不符时**整次分析失败**（422），
    只显示错误，不展示任何部分映射。

### 试用配对功能

打开 <http://localhost:8000/pair>，选择：

- 参考：`samples/demo.fa`
- 批次 A：`samples/pair_batch_a.vcf`
- 批次 B：`samples/pair_batch_b.vcf`

演示数据含 4 个样本（S1–S4）、20 个共有 SNP 位点与 2 个 `./.` 缺失调用；
两批的标签与数据被错位，预期唯一最优方案为
`S1↔S2、S2↔S4、S3↔S1、S4↔S3`，累计不一致 1。

JSON 接口：

```bash
curl -F 'fasta=@samples/demo.fa' \
     -F 'vcf_a=@samples/pair_batch_a.vcf' \
     -F 'vcf_b=@samples/pair_batch_b.vcf' \
     http://localhost:8000/api/pair
```

## 目录结构

```
app/            FastAPI 应用与核心逻辑（bio.py / pairing.py 仅用标准库）
  main.py       路由与上传校验（/compare 与 /pair、/api/compare 与 /api/pair）
  bio.py        双 VCF 核对：FASTA/VCF 解析、左对齐归一化、剂量与分组
  pairing.py    批次配对：多样本 VCF 解析、共有 SNP 统计、全局最优配对
  templates/    上传页与结果页（Jinja2）
tests/          单元 / 端到端测试
samples/        可直接试用的演示数据（核对与配对各一套）
scripts/        演示 VCF 生成器（坐标由参考自动计算）
Dockerfile, docker-compose.yml, requirements.txt
```

## 用 Docker Compose 启动（推荐）

```bash
docker compose up --build
```

浏览器入口：

- 双 VCF 变异核对：<http://localhost:8000/>
- 批次样本配对：<http://localhost:8000/pair>
- 交互式 API 文档（Swagger UI）：<http://localhost:8000/docs>

停止：`docker compose down`。容器以只读根文件系统、`no-new-privileges` 运行，
内置 `/` 健康检查。

## 本地运行（无 Docker 时）

```bash
python3 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

## 测试

```bash
# 核心逻辑（仅标准库）
python3 -I -m unittest discover -s tests -t . -p 'test_bio.py' -v
python3 -I -m unittest discover -s tests -t . -p 'test_pairing.py' -v

# 全部测试（Web 层需要先安装 requirements.txt）
python3 -I -m unittest discover -s tests -t . -v
```

当前共 95 个单元 / 端到端测试（含原有 48 个与批次配对的 47 个）。
