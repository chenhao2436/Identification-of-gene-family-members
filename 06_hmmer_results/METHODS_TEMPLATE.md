# CaTRM 家族鉴定 —— Methods 素材（模板 + 已填内容）

本文件把所有**不依赖服务器**的方法学信息预先填好；标注 `【待填】` 的项需在服务器跑完对应阶段后补齐。
服务器端跑完阶段 07 会自动生成带实际数字的 `07_final/CaTRM_Methods.txt`。

---

## 1. 数据来源

| 项 | 内容 |
|---|---|
| 辣椒参考蛋白组 | Zunla（`Zunla_Canz.pep.fa`），**52,385** 条序列 |
| 辣椒基因组注释 | `Canz.genome.gff3`（证据行 `Zunla3.0`），12 条染色体、**39,068** 个 gene、**52,385** 个 mRNA |
| 平均可变剪接体 | **1.34** 个转录本/gene（52,385 ÷ 39,068） |
| 拟南芥来源 | Araport11 蛋白组（`Araport11_pep_20250411.gz`） |
| AtTRM query | 34 个 AtTRM 基因，各取**最长转录本**作代表蛋白 |

**ID 体系**（三级一致，已用真实数据验证 39/39、24/24 命中）：

```
gene        ZLC01G0003110      ← GFF3 gene 特征；最终统计单位
  └ mRNA    ZLC01G0003110.1    ← GFF3 mRNA
      └ 蛋白 ZLC01G0003110.1    ← 蛋白 FASTA 序列 ID
```

---

## 2. query 构建与完整性验证

- 提取：从 Araport11 按基因 ID 提取全部转录本 → 逐基因比较氨基酸长度 → 取最长者
  （并列最长且序列相同时取转录本编号最小者）
- **完整性验证（已通过）**：

| 检查项 | 结果 |
|---|---|
| query 条数 | 34 |
| 与 Araport11 往返逐残基比对 | **34/34 一致，0 不一致** |
| 长度与注释表一致 | 34/34 |
| 字母表 | 仅 20 种标准氨基酸，无 `*`/`X`/空序列 |
| 长度范围 | 372 – 1025 aa，平均 693.9 aa |
| UniProt 抽查 | AT3G02170.1 与 Q9S823 **逐残基完全一致（905/905）** |

> 脚本：`06_hmmer_results/scripts/check_query_integrity.py`

---

## 3. 软件与版本

| 软件 | 用途 | 版本 |
|---|---|---|
| BLAST+ | 第一轮同源初筛 | 【待填】≥2.9（`qcovhsp` 要求） |
| MAFFT | 多序列比对 | 【待填】 |
| HMMER（hmmbuild/hmmsearch/hmmscan/hmmpress） | HMM 建库、扫描、结构域 | 【待填】 |
| Pfam-A | 结构域数据库 | 【待填】 |

---

## 4. 方法与参数

### 4.1 第一轮：BLASTP 初筛

```
blastp -query <34_AtTRM.fasta> -db <Zunla_Canz.pep> \
       -evalue 1e-5 \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend \
                sstart send evalue bitscore qlen slen qcovhsp" \
       -max_target_seqs 50 -num_threads <N> -seg no
```

- **实际使用的过滤阈值：identity ≥ 25%，query coverage ≥ 30%，E ≤ 1e-5**
  （由 `05_blast_results/01_summary_stats.txt` 实测确认）
- `-seg no`：TRM 富含重复序列，开启低复杂度屏蔽会屏蔽掉最需比对的区域
- 结果：**24 个候选基因 / 39 条候选蛋白**

### 4.2 第二轮：家族 HMM 建库与扫描

```bash
# 2a 多序列比对
mafft --auto --anysymbol --thread <N> AtTRM_34.fasta > AtTRM_34.aln.fasta

# 2b 建立家族 HMM
hmmbuild --amino --cpu <N> TRM.hmm AtTRM_34.aln.fasta

# 2c 扫描全蛋白组
hmmsearch --domtblout TRM_vs_Zunla.domtblout --tblout TRM_vs_Zunla.tblout \
          -E 1e-5 --cpu <N> TRM.hmm Zunla_Canz.pep.fa > TRM_vs_Zunla.hmmsearch.txt
```

- `--auto`：MAFFT 按序列数与相似度自动选算法（FFT-NS-2 / L-INS-i），**未针对本家族调优**
- **未做 gap 裁剪**：HMMER 建库对 gap 自带加权，裁剪会丢失信息
- 判定口径（两档，均按固定阈值统计，**未为凑数量调整**）：
  - `strict`：全序列 E-value ≤ 1e-5
  - `inclusive_only`：全序列 E > 1e-5 但最佳 domain 的 i-Evalue ≤ 1e-2
- `-E` 只决定**报告范围**，不改变打分；判定依据全部取自原始 `domtblout`
- 结果：strict 【待填】 基因 / 【待填】 蛋白；inclusive 【待填】 基因

### 4.3 候选合并

- 取 BLAST 与 HMM 候选的**并集**，按 **gene**（剥离 isoform 版本号）去重
- 标注四类证据：`BLAST+HMM` / `BLAST-only` / `HMM-only` / `HMM-inclusive-only`
- **不删除任何候选**；可疑者以 `tier` 与 `manual_review` 标记
- 结果：合并并集 【待填】 基因

### 4.4 结构域核查

```bash
hmmscan --domtblout CaTRM_vs_Pfam.domtblout --cpu <N> Pfam-A.hmm CaTRM_candidates_proteins.fasta
```

- 显著性阈值：i-Evalue ≤ 1e-5

**⚠️ 结构域解释的关键前提（务必写进 Methods 的局限讨论）**

经 InterPro 官方注释逐条核实：

| 条目 | accession | 长度 | 位置 | 家族判定价值 |
|---|---|---|---|---|
| LONGIFOLIA 1/2-like C-terminal domain | Pfam **PF14309** / InterPro **IPR025486** | ~180 aa | C 端 | ✅ 可用（3/3 成员命中，1e-30 量级） |
| DUF761-associated sequence motif (VARLMGL) | Pfam **PF14383** / InterPro **IPR032795** | 22–31 aa | N 端 | ⚠️ 不宜单独使用（3 条仅 1 条命中，score 0.01） |
| DUF3741 | Pfam **PF12552** / InterPro **IPR022212** | ~50 aa | N 端 | ⚠️ 未知功能域 |
| Protein LONGIFOLIA 1/2 | PANTHER **PTHR31680** / InterPro **IPR033334** | 近全长 | 1→C 端 | ✅ 最佳家族标记（GO:0051513） |

**要点：Pfam 在 TRM 家族中没有全长家族模型**，只有 C 端域 PF14309。
因此 **hmmscan(Pfam) 只能验证"C 端域存在"，不能证明"该蛋白属 TRM 家族"**；
全长家族归属依赖于 4.2 节自建的全家族 HMM（或 PANTHER PTHR31680）。

**N 端结构**：`PF14383 (VARLMGL 基序) + PF12552/DUF3741`——PF14383 的官方定位即"位于 DUF3741 的 N 端"。

结果：检出 PF14309 的候选 【待填】 / 【待填】

### 4.5 分层标准

| 层 | 条件 |
|---|---|
| `conservative` | BLAST 与 HMM 双重支持 **且** 检出 PF14309 |
| `standard` | BLAST 与 HMM 双重支持；或单方法支持但 HMM 达严格阈值 |
| `review` | 仅达 HMM 包含阈值 / 仅单方法且无结构域支持 / 覆盖度偏低 |

---

## 5. 各步候选数量

| 步骤 | 数量 |
|---|---|
| AtTRM 输入代表蛋白 | **34** |
| BLAST 主流（E≤1e-5）命中 | 32 条 query 有命中 |
| BLAST 候选基因（I≥25% / C≥30%） | **24** |
| BLAST 候选蛋白 | **39** |
| HMM strict 候选基因 | 【待填】 |
| HMM inclusive 候选基因 | 【待填】 |
| 合并去重后候选基因 | 【待填】 |
| ├ conservative | 【待填】 |
| ├ standard | 【待填】 |
| └ review（待人工确认） | 【待填】 |
| 最终代表蛋白序列 | 【待填】 |

**BLAST 侧覆盖情况（已离线审计）**：
34 条 query 中 32 条至少有一条 E≤1e-5 命中；**2 条完全零命中**：

| query | At 基因 | TRM | 说明 |
|---|---|---|---|
| `AT5G01370.1` | At5g01370 | TRM29, ACI1 | 较短（427 aa），兼有 ALC-interacting protein 1 身份，家族归属需复核 |
| `AT5G58630.1` | At5g58630 | TRM31 | 全家族最短（372 aa） |

每条 At 基因贡献 1–4 个候选，中位数 2，无异常离群。
> 脚本：`06_hmmer_results/scripts/audit_blast_coverage.py`

---

## 6. 已记录的方法学限制（建议写入讨论）

1. **Pfam 无全长家族模型**（见 4.4）：结构域阳性不能单独作为家族成员证据。
2. **BLAST 阈值与早期文档不一致**：实际为 I≥25%/C≥30%（非 30%/50%），
   所有结果均以实测值为准；若需两法阈值统一，需重跑 BLAST 并声明。
3. **候选汇总表曾有缺表头缺陷**：`01_candidate_summary.tsv` 因流水线缺陷丢失表头，
   已定位并修复；现有文件为修复前产物，仅影响人工阅读，不影响解析。
4. **E 值阈值的潜伏缺陷（已修复，本数据集无影响）**：
   流水线曾以 `awk -v E="$EVALUE"` 配合 `(E+0)` 比较，awk 的 strtod 会把字符串
   `"1e-5"` 按十进制解析为 `1`，使条件退化为 `E ≤ 1`，**E 值过滤实际失效**。
   已改为传数值字面量并去掉 `+0`。
   **对本数据集的影响已量化：三档（主候选/高可信/片段）共 245 行，最大 E 值 7.14e-06，
   无任何行因该缺陷被错误放行**，故 24 基因 / 39 蛋白的结论不受影响。
5. **PANTHER 全长核查未自动完成**：因本机网络限制（HTTPS 不可用）无法提交
   InterProScan，列为可选人工步骤。

---

## 7. 可引用文献

| PMID | 期刊/年份 | 要点 |
|---|---|---|
| 22286137 | Plant Cell, 2012 | TRM1–TON1 互作，植物皮层微管阵列与真核中心体的共同招募网络 |
| 17038516 | Development, 2006 | LONGIFOLIA1/2 调控拟南芥纵向细胞伸长 |
