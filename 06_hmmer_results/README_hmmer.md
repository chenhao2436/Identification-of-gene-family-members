# CaTRM 家族鉴定流程（HMMER 部分）

Zunla 辣椒 TRM 基因家族鉴定的第二半程：**MAFFT → hmmbuild → hmmsearch → 与 BLAST 合并 → Pfam 结构域核查 → 最终 CaTRM 列表**。

前半程（34 条 AtTRM 提取、BLASTP 初筛）在 `../02_script_output/` 与 `../05_blast_results/`。

---

## 一、快速开始

```bash
cd 06_hmmer_results/scripts

bash run.sh                    # 列出 7 个阶段与运行状态
bash run.sh 01                 # 只跑阶段 01（环境与输入文件检查）← 建议先跑这个
bash run.sh 02                 # MAFFT 比对
bash run.sh 03                 # hmmbuild 建 TRM.hmm
bash run.sh 04                 # hmmsearch 扫描 Zunla 全蛋白
bash run.sh 05                 # 与 BLAST 候选合并
bash run.sh 06                 # Pfam 结构域核查（需要 Pfam-A.hmm）
bash run.sh 07                 # 汇总最终 CaTRM 列表
bash run.sh all                # 依次全跑
```

**分阶段可单独执行**，每步独立可重跑，产物落在对应阶段目录，互不覆盖。

透传参数：`--gff3 <file>` `--genome <file>` `--pfam <Pfam-A.hmm>` `--threads N`

```bash
bash run.sh 06 --pfam /data/db/Pfam-A.hmm --threads 16
```

---

## 二、阶段与产物

| 阶段 | 做什么 | 输出目录 | 依赖 |
|---|---|---|---|
| 01 | 环境与输入检查（**只读**） | `01_env_check/` | 无 |
| 02 | `mafft --auto` 比对 34 条 AtTRM | `02_msa/` | mafft |
| 03 | `hmmbuild` 建 `TRM.hmm` | `03_hmmbuild/` | 阶段 02 |
| 04 | `hmmsearch` 扫描 Zunla 全蛋白 | `04_hmmsearch/` | 阶段 03 + Zunla 蛋白 |
| 05 | 与 BLAST 候选合并去重 | `05_merge/` | 阶段 04 + `05_blast_results/` |
| 06 | `hmmscan` Pfam 结构域核查 | `06_domain/` | 阶段 05 + **Pfam-A.hmm** |
| 07 | 汇总最终 CaTRM 列表 | `07_final/` | 阶段 05 + 06 |

每个阶段产出三件套：`<阶段>.log`（完整日志）、`<阶段>.subset.tsv`（机读小结）、`<阶段>_report.txt`（中文结论）。

---

## 三、方法、参数与判定口径（可直接写进论文 Methods）

### 3.1 多序列比对（阶段 02）

```bash
mafft --auto --anysymbol --thread <N> AtTRM_34.fasta > AtTRM_34.aln.fasta
```

- `--auto`：MAFFT 按序列数与相似度自动选择算法（FFT-NS-2 / L-INS-i），官方推荐的通用默认，**未针对本家族调优**。
- **未做 gap 裁剪**（trimAl/Gblocks）：HMMER 建库对 gap 自带加权，裁剪会丢信息。如需裁剪属方法变更。

### 3.2 家族 HMM 构建（阶段 03）

```bash
hmmbuild --amino --cpu <N> TRM.hmm AtTRM_34.aln.fasta
```

自建家族 HMM，未引入外部 seed。hmmbuild 默认使用 Henikoff 位置特异性权重 + 有效序列数校正（EFFN 会记录在 `TRM.hmm.header.txt`）。

### 3.3 全蛋白组扫描（阶段 04）

```bash
hmmsearch --domtblout TRM_vs_Zunla.domtblout --tblout TRM_vs_Zunla.tblout \
          -E 1e-5 --cpu <N> TRM.hmm Zunla_Canz.pep.fa > TRM_vs_Zunla.hmmsearch.txt
```

**判定口径（两档，均按固定阈值统计，不凑数）**：

| 档位 | 条件 | 含义 |
|---|---|---|
| `strict` | 全序列 E-value ≤ 1e-5 | 主口径 |
| `inclusive_only` | 全序列 E > 1e-5 但最佳 domain 的 i-Evalue ≤ 1e-2 | Rfam/Pfam 惯例的包含阈值，最宽松档 |
| `below_inclusive` | 连 i-Evalue ≤ 1e-2 都不到 | 保留在 `parsed.tsv`，不计入候选 |

> `-E` 只决定**输出哪些命中**，不改变打分。所有判定依据都取自原始 `domtblout` 的 E-value / i-Evalue。**没有为了让候选数落在某个区间而调整任何阈值。**

### 3.4 BLASTP 初筛（已有结果）

沿用 `05_blast_results/`，**未重跑**。实际使用的阈值会从日志中提取并记录（注意：实测为 identity ≥25%、coverage ≥30%，与早期文档写的 30%/50% 不同，报告里以实际值为准）。

### 3.5 结构域核查（阶段 06）

```bash
hmmscan --domtblout CaTRM_vs_Pfam.domtblout --cpu <N> Pfam-A.hmm CaTRM_candidates_proteins.fasta
```

显著性阈值 `i-Evalue ≤ 1e-5`。

**重点关注的两个域（已用 InterPro 官方注释实测核实）**：

| Pfam | 名称 | 位置 | 本流程的态度 |
|---|---|---|---|
| **PF14309** | LONGIFOLIA 1/2-like C-terminal domain (DUF4378) | **C 端**（约 180 aa） | **主判定依据**。3/3 条实测的 AtTRM 都命中，信号很强（1e-30~1e-35） |
| **PF14383** | DUF761-associated sequence motif (VARLMGL) | **N 端**（仅 22–31 aa） | **辅助证据**。3 条里仅 1 条命中，且 score 仅 0.01（边缘） |

**实测证据（InterPro API，逐条核对）**：

| 蛋白 | 长度 | PANTHER 家族模型 | PF14309 (C 端) | PF14383 (N 端) |
|---|---|---|---|---|
| LNG2 / TRM1 (Q9S823) | 905 | PTHR31680:SF15，**覆盖 1–902，score 0** | 722–901，**4.1e-30** | **未命中** |
| LNG1 / TRM2 (Q9LF24) | 927 | PTHR31680:SF22，**覆盖 1–925，score 0** | 721–904，**7.8e-35** | **未命中** |
| TRM4 (Q0WNQ5) | 1025 | PTHR31680:SF4，**覆盖 1–1022，score 0** | 820–1002，**1.8e-33** | 274–295，22 aa，**score 0.01（边缘）** |

### ⚠️ 关键结论：Pfam 无法证明家族归属

- **PANTHER 的 `PTHR31680`（IPR033334，"Protein LONGIFOLIA 1/2"，GO:0051513）覆盖蛋白质全长**，
  是唯一能代表 TRM 家族全长的模型。
- **Pfam 在 TRM 家族中只有 C 端域 PF14309**，**没有全长家族模型**；
  N 端在 Pfam 里只有一条 22–31 aa 的边缘短基序 PF14383。
- 因此：**阶段 06 的 hmmscan(Pfam) 只能验证"C 端域是否存在"，
  不能证明"这条蛋白是 TRM 家族成员"**。这正是必须自建全家族 HMM
  （阶段 02–04）的根本原因，两件事不可互相替代。

**N 端的真实组成**（此前我判断不完整，现已查清）：

```
N 端:  PF14383 (VARLMGL 基序, 22-31 aa)  +  PF12552 / DUF3741 (约 50 aa, IPR022212)
       └─ PF14383 的官方定位就是「位于 DUF3741 的 N 端」
C 端:  PF14309 (LONGIFOLIA 1/2-like C-terminal domain, 约 180 aa)
```

**判读建议**：
- PF14309 位于 C 端、域长约 180 aa、E 值 1e-30 ~ 1e-35 —— 用作 `conservative` 层的必要条件是合理的。
- PF14383 是短基序、信号从 `0.01` 到 `3.1e-10` 波动、多数成员未命中 —— **不宜单独作为判定依据**。
- 注意 Pfam 对 PF14383 有已知的**重复命中抑制规则**，因此"未报 PF14383"**不等于**没有该段序列。
- 若要把"N 端模块"也纳入证据，应查 **PF12552/DUF3741**，而不是只看 PF14383。

脚本在阶段 06 会用**你自己的 34 条 AtTRM 实测**这两个域的命中情况（`focus_hits.tsv`），报告里给出结论，不预设。

### 可选补充证据：PANTHER 全长家族模型

由于 Pfam 无全长模型，**建议在阶段 07 之后加一步 PANTHER 核查**作为独立佐证：

- 服务：InterProScan（在线或本地），勾选 **PANTHER**（或同时勾 CDD/SMART）
- 目标条目：`PTHR31680`（IPR033334，Protein LONGIFOLIA 1/2）
- 预期：真正的 CaTRM 成员应被 PTHR31680 覆盖大部分序列长度、score 很低（0 量级）
- 提交对象：`07_final/CaTRM_final.proteins.fasta`

> **为什么本流程没把它自动化**：本机沙箱的 shell 只允许 HTTP（HTTPS 因 TLS 凭证
> 无法建立而失败），而 InterProScan 的提交接口需要 HTTPS POST，`web_fetch` 只支持 GET，
> 因此**无法在此环境自动提交序列**。故 PANTHER 核查列为可选人工步骤。
> 若你的服务器能访问外网，可自行批量提交；或用本地 InterProScan 跑（见第五节）。

### 3.6 候选合并（阶段 05）

取 BLAST 与 HMM 候选的**并集**，按 **gene**（剥离 isoform 版本号）去重，标注四类证据：

| 标签 | 含义 |
|---|---|
| `BLAST+HMM` | 两种方法都显著支持 |
| `BLAST-only` | 仅 BLAST 支持（HMM 未显著命中） |
| `HMM-only` | 仅 HMM 严格命中（BLAST 漏检——profile HMM 对远缘同源更敏感） |
| `HMM-inclusive-only` | 仅达 HMM 包含阈值 |

**不删除任何候选。**

### 3.7 分层（阶段 07）

| 层 | 条件 |
|---|---|
| `conservative` | BLAST 与 HMM 双重支持 **且** 检出 PF14309 |
| `standard` | BLAST 与 HMM 双重支持；或单方法支持但 HMM 达严格阈值 |
| `review` | 仅达 HMM 包含阈值 / 仅单方法且无结构域支持 / 覆盖度偏低 |

`review` 层单列于 `review.tsv`，**仅标记不删除**。若论文只用 conservative + standard，需说明 review 层的处置理由。

---

## 四、ID 命名体系（三者关系）

Zunla 的 GFF3 与蛋白 FASTA 使用同一套编号：

```
gene        ZLC01G0003110        ← GFF3 的 gene 特征；最终统计单位
  └ mRNA    ZLC01G0003110.1      ← GFF3 的 mRNA
      └ 蛋白 ZLC01G0003110.1      ← 蛋白 FASTA 的序列 ID
```

剥离版本号（`.1`/`.2`）即得 gene ID。**所有基因级统计都以 gene 为单位。**

---

## 五、需要准备的数据库

| 资源 | 用途 | 获取方式 |
|---|---|---|
| **Pfam-A.hmm** | 阶段 06 结构域核查 | `wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz` → `gunzip` → `hmmpress Pfam-A.hmm`（约 1.5–2 GB） |
| InterProScan（可选） | 含 CDD/SMART/PROSITE，覆盖最全 | 约 10–20 GB；适合需要正式 InterPro 注释时 |
| CDD / SMART（可选） | 交叉验证 | NCBI CDD 或 SMART 在线批量提交 |

阶段 06 若找不到 Pfam 库，会**明确报告为阻塞项**并给出三种方案，不会静默跳过；候选蛋白 FASTA 仍会产出，可直接用于在线提交。

---

## 六、运行前检查（阶段 01）

```bash
cd 06_hmmer_results/scripts
bash run.sh 01
cat ../01_env_check/01_env_check_report.txt
```

检查项：软件（mafft/hmmbuild/hmmsearch/hmmscan/hmmpress/interproscan/blast/python）、Pfam 库位置与 hmmpress 索引、输入文件（存在性/格式/条数）、**ID 交叉核对**（候选基因对 GFF3、候选蛋白对蛋白库）、BLAST 结果自洽性与实际阈值、磁盘余量。

输出 `01_env_check_summary.tsv` 为机读表（检查项·期望·实测·结论），阻塞项与警示项计数在报告头部。

### 6.1 已有 BLAST 候选的覆盖审计（已离线完成）

```bash
python scripts/audit_blast_coverage.py
# 产出 05_merge/blast_query_coverage.tsv（At 基因 -> 候选 Zunla 基因）与 .txt 报告
```

实测结果（基于现有 24 基因 / 39 蛋白候选）：

| 项 | 值 |
|---|---|
| AtTRM query 总数 | 34 |
| 至少有一条 E≤1e-5 命中的 | **32** |
| 进入主候选的 | 32（**没有 query 被 I/C 阈值刷掉**） |
| 主候选覆盖的 At 基因 | 32 / 34 |
| 每个 At 基因的候选数 | 1–4 个，**中位数 2**，无异常离群 |

**2 条完全零命中的 query（重点观察对象）**：

| query | At 基因 | TRM 符号 | 说明 |
|---|---|---|---|
| `AT5G01370.1` | At5g01370 | **TRM29, ACI1** | 蛋白较短（427 aa），且兼有 ALC-interacting protein 1 的身份 |
| `AT5G58630.1` | At5g58630 | TRM31 | 全家族最短（372 aa） |

这 2 条**正是 HMMER 该补检的对象**：BLAST 单条比对对短/远缘序列不够敏感，而 34 条建成的 profile HMM 覆盖了全家族保守模式，很可能把它们捞回来。若阶段 04 后仍有 query 无对应候选，才可判断为"辣椒中丢失/拟南芥特有"，并写入结论。

> 注 `AT5G01370` 同时是 ACI1：它是否属经典 TRM 家族值得复核。若 HMMER 也不支持，建议在最终列表里单列说明，而非直接计入家族。

---

## 七、验证情况与已知限制（重要）

### 已在本机验证的部分

```bash
python scripts/test_hmmer_parsers.py                                    # 53 条断言
python scripts/check_id_mapping.py                                      # 10 条断言（真实数据）
python ../02_script_output/06_verify_bash_structure.py scripts/*.sh     # 9 个脚本结构检查
```

| 验证项 | 结果 |
|---|---|
| 四个解析/合并/汇总脚本（合成 HMMER 输出） | 53 条断言全通过 |
| **Zunla ID 映射（真实 GFF3 + 真实候选）** | 10 条断言全通过：候选基因 24/24 命中 GFF3 `gene` 特征；候选蛋白 39/39 命中蛋白库；39 个候选蛋白同时是 GFF3 的 mRNA（gene/transcript/protein 三级命名一致） |
| GFF3 规模 | 12 条染色体、39,068 个 gene、52,385 个 mRNA（平均 1.34 转录本/基因） |
| bash 结构（括号/关键字/heredoc/CRLF/`grep -c` 陷阱） | 9 个脚本全通过 |

### ⚠️ 已知限制：7 个阶段脚本尚未真正执行过

本机沙箱**禁止 bash 创建 signal pipe**（`bash: fatal error - couldn't create signal pipe`），
因此**无法在本机跑任何 shell 脚本**，也做不了端到端冒烟测试。

我尝试过用假二进制（mock mafft/hmmbuild/hmmsearch/hmmscan）搭建离线冒烟测试，
但同样被该限制挡住，故已删除该测试脚本（不留跑不起来的测试）。

**这意味着**：阶段之间的产物衔接、外部程序调用语法，**只经过静态检查与人工审查，
未经实机执行**。请在服务器上**先跑阶段 01**（纯只读，风险最低），确认无误后再逐步推进。

若某阶段报错，优先看该阶段的 `.log` 与 `_report.txt`；阶段脚本都设计成失败即明确报错，
不会静默产出空结果。

### 6.2 顺带修复：BLAST 流水线的候选汇总表丢失表头

审计时实测发现 `05_blast_results/01_candidate_summary.tsv` **没有表头**（10 列纯数据）。
定位到 `02_script_output/03_trm_blast_pipeline.sh` 阶段 3 的缺陷：

```bash
{
  printf 'query\ttrm_gene\t...\ttier\n'      # 表头写到 stdout（即文件）
  if [[ -n "$MAPFILE" && -f "$MAPFILE" ]]; then
    awk ... > "$SUMMARY"                     # ← awk 再次打开同一文件并截断，表头被冲掉
  ...
```

**最小修复**：改为子 shell 单次重定向，表头与 awk 输出共用同一个 fd：

```bash
(
  printf 'query\ttrm_gene\t...\ttier\n'
  if [[ -n "$MAPFILE" && -f "$MAPFILE" ]]; then
    awk ...                                  # 不再自带 > "$SUMMARY"
  ...
) > "$SUMMARY"
```

同时补上了 else 分支 awk 缺失的 `-F'\t'`（原先字段切分会错，导致输出退化为单列）。
已用 Python 复现并验证修复语义（表头保留 + 数据完整 + 列数一致）。

> 你现有的 `01_candidate_summary.tsv` 是**修复前**产物（无表头）。它只影响人工阅读，
> 下游解析（`merge_blast_hmm.py` 等）按位置列读取、对两种格式都兼容，因此**无需重跑**。
> 下次重跑 BLAST 时会自动带上表头。

---

## 八、设计上的几个取舍（供审阅）

1. **`-E 1e-5` 作初筛但另记 inclusive 档**：既保留标准严格口径，又不在数据里丢掉边界候选，两档数量都进报告。
2. **不做 gap 裁剪**：见 3.1。
3. **PF14383 用实测说话**：见 3.5，已用 InterPro 官方注释核实其共现性与弱信号特征。
4. **零删除原则**：所有候选保留，用 `tier` 与 `manual_review` 列表达置信度。
5. **每步留痕**：软件版本、完整参数、候选数量都写入各阶段报告与 `07_final/CaTRM_Methods.txt`。
