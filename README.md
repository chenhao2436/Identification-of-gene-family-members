# 拟南芥 TRM 基因家族成员鉴定（Zunla 辣椒）

用拟南芥 TRM 家族 34 条代表蛋白作 query，鉴定 Zunla 辣椒中的 TRM 家族成员。
本目录是该工作的**独立项目包**，自带输入、脚本与全部产出。

**两阶段流程**：

```
第一阶段：BLASTP 初筛                     第二阶段：HMMER 精筛与确认
  34 条 AtTRM                              34 条 AtTRM
      │                                        │  MAFFT 比对
      ▼  blastp vs Zunla 全蛋白                ▼  hmmbuild → TRM.hmm
  BLAST 候选（24 基因 / 39 蛋白）           hmmsearch 扫描 Zunla 全蛋白
      │                                        │
      └──────────────┬─────────────────────────┘
                     ▼  并集合并去重（按 gene）
              标注 BLAST+HMM / BLAST-only / HMM-only / HMM-inclusive-only
                     ▼  Pfam 结构域核查（PF14309 等）
                     ▼  分层：conservative / standard / review
                  最终 CaTRM 基因/蛋白列表
```

- **第一阶段** → `02_script_output/` + `05_blast_results/`
- **第二阶段** → `06_hmmer_results/`（详见 [README_hmmer.md](06_hmmer_results/README_hmmer.md)）

---

## 一、目录结构

```
Identification_of_gene_family_members/
│
├── 01_At_TRM_query/                  输入数据
│   ├── Zunla_Canz.pep.fa               ★ Zunla 辣椒蛋白组（52,385 条，21 MB）
│   ├── Araport11_pep_20250411.gz       Araport11 全蛋白组（8.5 MB，来源数据）
│   └── TRM基因家族注释表.xlsx           你提供的 TRM 家族注释表
│                                      （BLAST query 在 02_script_output/①）
│
├── 02_script_output/                 ★ 第一阶段：六个核心文件 + 运行脚本
│   ├── 01_At_TRM_query_34.fasta        ★① BLAST query：34 条 TRM 代表蛋白
│   ├── 02_At_TRM_annotation.tsv        ★② query 注释表（基因/转录本/TRM符号/长度）
│   ├── 03_trm_blast_pipeline.sh        ★③ BLAST 鉴定流水线（一键运行）
│   ├── 04_summarize_candidates.py      ★④ 候选统计与出图
│   ├── 05_verify_pipeline_logic.py     ★⑤ 流水线逻辑验证（22 条断言）
│   └── 06_verify_bash_structure.py     ★⑥ bash 结构校验器（含默认路径检查）
│
├── 05_blast_results/                 第一阶段产出：BLAST 候选
│   ├── TRM_vs_Zunla.blast.tsv          原始结果（一对多）
│   ├── TRM_candidates_gene_level.txt   ★ 候选基因 ID（24）
│   ├── TRM_candidates_protein_level.txt ★ 候选蛋白 ID（39）
│   ├── Zunla_TRM_candidate_proteins.fasta
│   └── 01_*.tsv|txt|png                统计表/报告/分布图
│
├── 06_hmmer_results/                 ★ 第二阶段：HMMER 精筛与家族确认
│   ├── README_hmmer.md                 流程说明（方法/参数/判读，可写进论文）
│   ├── scripts/                        7 个阶段脚本 + 主入口（分阶段可单独执行）
│   │   ├── run.sh                      主入口：bash run.sh <01-07|all>
│   │   ├── common.sh                   共享配置与工具函数
│   │   ├── stage_01_env_check.sh       环境与输入文件检查（只读）
│   │   ├── stage_02_msa.sh             MAFFT 比对
│   │   ├── stage_03_hmmbuild.sh        建立 TRM.hmm
│   │   ├── stage_04_hmmsearch.sh       扫描 Zunla 全蛋白
│   │   ├── stage_05_merge.sh           与 BLAST 候选合并
│   │   ├── stage_06_domain.sh          Pfam 结构域核查
│   │   ├── stage_07_final.sh           汇总最终 CaTRM 列表
│   │   ├── parse_hmmsearch_domtbl.py   解析 hmmsearch（纯标准库）
│   │   ├── parse_hmmscan_domtbl.py     解析 hmmscan（纯标准库）
│   │   ├── merge_blast_hmm.py          合并去重（四类证据标签）
│   │   ├── final_summary.py            最终列表 + 证据与理由
│   │   └── test_hmmer_parsers.py       合成数据测试（53 条断言）
│   └── 01_env_check/ … 07_final/       各阶段产出（运行后生成）
│
├── 03_auxiliary/                     辅助材料（前序步骤与备查）
│   ├── 00_peek_fasta_format.py         FASTA 头格式查看
│   ├── 01_extract_trm_proteins.py      按基因 ID 提取 At TRM 蛋白
│   ├── 02_select_representative.py     每基因挑最长转录本作代表蛋白
│   ├── At_TRM_all_isoforms_87.fasta    全部 87 条转录本（备用）
│   ├── At_TRM_all_isoforms_87.tsv      全转录本汇总表
│   ├── select_representative_report.txt 代表蛋白挑选过程报告
│   ├── not_found_genes.txt             未匹配基因（当前为空 = 34 个全部匹配）
│   └── TRM基因家族注释表_copy.xlsx       注释表副本
│
├── 04_temp/                          ★ 临时文件（可随时清理，不入 git）
│   ├── README.md                       本目录说明与 SnapGene 文件注意事项
│   └── SnapGene_prot_files/            34 个 SnapGene 二进制 .prot（非 FASTA，不可用于 BLAST）
```

> 版本控制忽略规则在本仓库根的 `.gitignore`（忽略两个大参考数据文件、BLAST/HMMER 运行产物与 `04_temp/`）。

---

## 二、第一阶段：六个核心文件说明

| # | 文件 | 作用 | 在哪运行 |
|---|---|---|---|
| ① | `01_At_TRM_query_34.fasta` | **BLAST 的 query**：34 个 At TRM 基因各取一条最长转录本 | 输入 |
| ② | `02_At_TRM_annotation.tsv` | query 注释表，含 `gene_id` 列用于把结果映射回 At 基因号 | 输入 |
| ③ | `03_trm_blast_pipeline.sh` | **主流水线**：校验→建库→比对→过滤分档→去重→抽序列→统计出图 | 零参数一键运行 |
| ④ | `04_summarize_candidates.py` | 被 ③ 自动调用；I/C 分布统计、断崖检测、三联图 | 无需手动运行 |
| ⑤ | `05_verify_pipeline_logic.py` | 22 条断言验证阈值边界/三档/去重逻辑 | 本地或服务器 |
| ⑥ | `06_verify_bash_structure.py` | 静态检查 bash 结构 + 核对默认路径 | 本地或服务器 |

> 注：③④⑤⑥ 无第三方依赖（④ 的绘图在无 matplotlib 时自动降级为 CSV+文本）。
> ③ 会自动在同目录下查找 ④，因此两者**必须放在同一目录**。

---

## 三、鉴定流程

```
01_At_TRM_query/01_At_TRM_query_34.fasta   (34 条)
              │  03_trm_blast_pipeline.sh
              ▼  blastp  vs  Zunla_Canz.pep.fa (52,385 条)
         ① 原始 BLAST 结果（一对多）
              │  过滤 + 分档
              ▼
         ② 候选 Zunla 蛋白 ID（去重）  →  ③ 候选蛋白 FASTA
```

### 阈值设定与理由

| 指标 | 取值 | 理由 |
|---|---|---|
| `-evalue` | **1e-5** | **特异度主力**，家族鉴定的常规起点 |
| `pident` | **≥25%** | 长比对 pident 天然偏低；卡 30% 会挡掉 E 值最好的一批全长命中 |
| `qcovhsp` | **≥30%** | 只作粗筛，剔除"仅局部重复区匹配"的假阳性 |

**关键理念：特异度由 E-value 承担，coverage 只做粗筛。** TRM 成员互为同源时只在
C 端 TRM 结构域那一段对齐，N 端 MORN 重复高度发散，真直系同源的全长 query 覆盖率
只有 26–46%——所以 coverage 卡在 50% 会稳定漏掉真成员。

这三个取值是从本数据集实测反推的：coverage 50%→30%、identity 30%→25%，候选基因数
从 18 → 22 → 24 后即**饱和**（再把 E 放宽到 1e-3、coverage 到 20%、identity 到 20%
均不再新增），说明 24 是搜索的饱和点而非阈值人为截断。注意饱和只说明"搜不到了"，
不等于"都是真成员"——交付前仍需用 TRM 结构域 + MORN 重复做复核。

收紧回旧口径：`--ident 30 --cov 50`（得 18 个）。

### 三档输出

| 档位 | 条件 | 含义 | 用途 |
|---|---|---|---|
| **候选（主）** | C≥30% 且 I≥25% 且 E≤1e-5 | ★ 主交付集 | 建树、家族规模统计 |
| 高可信 | C≥70% 且 E≤1e-10 | 全长直系同源 | 核心成员 |
| 片段待查 | E≤1e-5 但 C<30% | 仅局部匹配，疑似假阳性 | 人工复查 |

去重给两份：**蛋白级**（保留 `.1/.2`，用于取序列建树）+ **基因级**（剥离 isoform
后缀，用于统计真实"候选基因"数）。

---

## 四、运行步骤

脚本**自动识别本项目的文件路径**，不需要先搬文件、也不需要写一堆 `--query/--db-fasta`。
把项目整个拷到服务器（或用 git clone），加载 BLAST+ 后直接跑即可。

### 1. 加载 BLAST+（三选一）

```bash
module load blast+/2.14.0        # HPC
conda activate blast             # conda
# 或下载 NCBI 静态包后 export PATH=<...>/ncbi-blast-*/bin:$PATH
blastp -version                  # 必须 >= 2.9（否则无 qcovhsp 字段）
```

### 2. 一键运行（零参数）

```bash
cd 02_script_output
bash 03_trm_blast_pipeline.sh
```

就这样。默认值全部相对脚本自身定位，**在任意目录下调用都能正确解析**：

| 参数 | 默认值（自动识别） |
|---|---|
| `--query` | `02_script_output/01_At_TRM_query_34.fasta`（34 条代表蛋白） |
| `--db-fasta` | `01_At_TRM_query/Zunla_Canz.pep.fa`（52,385 条） |
| `--map` | `02_script_output/02_At_TRM_annotation.tsv`（Query 注释表） |
| `--outdir` | `05_blast_results/`（结果输出，脚本自动创建） |
| `--threads` | `4` |

**只调线程数**（常见需求）：

```bash
bash 03_trm_blast_pipeline.sh --threads 16
```

**文件放在别处时**再显式指定（其余仍走默认）：

```bash
bash 03_trm_blast_pipeline.sh \
  --query    /path/to/TRM_34.fasta \
  --db-fasta /path/to/Canz.pep.fa \
  --outdir   /path/to/out
```

`--help` 看全部参数。流水线每阶段打印进度，全程落盘到
`05_blast_results/01_diagnostics.txt`，结束给出核心数字与判读建议。

> 前序脚本支持 `TRM_INPUT_DIR` / `TRM_OUTPUT_DIR` 环境变量覆盖输入输出目录，
> 便于隔离测试；流水线本身用命令行参数即可。

### 3. 产出文件（均在 `05_blast_results/` 下）

| # | 文件 | 作用 |
|---|---|---|
| ① | `TRM_vs_Zunla.blast.tsv` | 原始 BLAST 结果（15 列，一对多） |
| ② | `TRM_candidates_protein_level.txt` | 去重候选**蛋白** ID（含 isoform） |
| ② | `TRM_candidates_gene_level.txt` | 去重候选**基因** ID ← **家族规模看这个** |
| ③ | `Zunla_TRM_candidate_proteins.fasta` | 去重候选蛋白 FASTA（建树直接可用） |
| A1 | `01_candidate_summary.tsv` | 候选表：query / At基因 / TRM符号 / Zunla ID / I / E / bitscore / C / 档位 |
| A2 | `01_candidate_counts.tsv` | 每 query 候选数 |
| A3 | `01_best_hit_per_query.tsv` | 每 query 最佳命中 |
| A4 | `01_zero_hit_queries.txt` | 零命中 query |
| A5 | `01_distribution.tsv` | I/C 分箱频数表 |
| A6 | `01_plots.png` | 三联图（覆盖率/一致性直方图 + 散点图，带阈值线） |
| A7 | `01_summary_stats.txt` | 统计报告 + 断崖检测 + 中文判读建议 |
| A8 | `Zunla_TRM_candidate_proteins_genelevel.fasta` | 基因级候选全部 isoform |
| A9 | `01_diagnostics.txt` | 全流程日志 |

---

## 五、本地验证

两个验证脚本都在 `02_script_output/`，可独立运行（把脚本路径作为参数即可）：

```bash
cd 02_script_output
python 05_verify_pipeline_logic.py                            # 22 条断言
python 06_verify_bash_structure.py 03_trm_blast_pipeline.sh   # bash 结构静态检查
python 06_verify_bash_structure.py --check-defaults 03_trm_blast_pipeline.sh
                                                              # 核对零参数默认路径是否可用
```

`--check-defaults` 会复刻流水线的路径探测逻辑，逐项确认 `--query`（34 条）、
`--db-fasta`（52,385 条）、`--map` 的默认路径在磁盘上真实存在、条数正确，
**确保"一键运行"不会因找不到文件而失败**。

已验证通过：阈值边界（`I=30/C=50/E=1e-5` 恰好通过，`29.99/49/1e-4` 被排除）、
三档互斥、isoform 剥离（含两位数 `.12`）、零命中识别、最佳命中按 bitscore 取值、
默认路径解析。

**端到端复现测试**（隔离临时目录重跑前序流水线）已通过：`01_extract` + `02_select`
的四个产出与项目内 `01_At_TRM_query_34.fasta`、`02_At_TRM_annotation.tsv`、
`At_TRM_all_isoforms_87.fasta/tsv`、`select_representative_report.txt`
**SHA256 完全一致**，说明归位后的路径与命名可靠、结果可复现。

> 前序脚本（`03_auxiliary/`）支持用环境变量指定输入/输出目录，便于隔离测试：
> `TRM_INPUT_DIR`（默认 `../01_At_TRM_query`）、`TRM_OUTPUT_DIR`（默认 `../02_script_output`）。

---

## 六、结果判读与调整

1. **优先看 coverage 分布**（`01_distribution.tsv` / `01_plots.png`）：重复序列家族里
   coverage 的断崖通常比 identity 干净，断崖位置比固定数字更可靠。
2. 候选基因数 **<20** → 偏少，用 `--evalue 1e-3` 重跑。
3. 候选基因数 **>100** → 偏多，用 `--cov 70 --evalue 1e-10` 重跑。
4. 某 query **零命中** → 该亚家族可能拟南芥特有或辣椒中丢失，本身是有意义的结论。
5. 高分但 **coverage <30%** 的命中 → 典型低复杂度假阳性。

---

## 七、注意事项

- **BLAST+ 必须 ≥2.9**：`qcovhsp` 字段在更老版本不存在，流水线阶段 0 会提前拦下。
- **`-seg no` 不能省**：TRM 富含重复序列，开启低复杂度屏蔽会屏蔽掉最该比对的区域。
- **query 用 34 条代表蛋白**，不要用 87 条全转录本——多转录本会互相抢命中，
  导致 top hit 混乱、候选数虚高。
- **04_temp/SnapGene_prot_files/** 里 34 个 `.prot` 是 **SnapGene 二进制格式**，
  不是文本 FASTA，**不能用于 BLAST**（已实测确认文件头为 `SnapGene` 魔数）。
  如需使用请在 SnapGene 里 `File → Export → Protein FASTA` 重新导出。

---

## 八、数据来源

| 数据 | 来源 | 说明 |
|---|---|---|
| At TRM 序列 | `input/Araport11_pep_20250411.gz`（Araport11） | 34 基因 → 87 转录本 → 34 代表蛋白 |
| Zunla 辣椒蛋白组 | `E:\华为\fsdownload\REF\Canz.pep.fa` | 52,385 条，20,285,222 残基，ID 形如 `ZLC01G0000010.1` |

> 本目录内的 `Zunla_Canz.pep.fa` 与参考目录下的原文件是**硬链接**（共享同一份数据，
> 不额外占空间）；修改其一即同时影响两者，请勿在本目录内改动该文件。
