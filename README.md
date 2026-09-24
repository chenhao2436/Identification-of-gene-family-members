# 拟南芥 TRM 基因家族成员鉴定（Zunla 辣椒）

用拟南芥 TRM 家族 34 条代表蛋白作 query，比对 Zunla 辣椒蛋白组，鉴定辣椒中的
TRM 候选家族成员。本目录是该工作的**独立项目包**，自带输入、脚本与全部产出。

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
├── 02_script_output/                 ★ 六个核心文件 + 运行脚本
│   ├── 01_At_TRM_query_34.fasta        ★① BLAST query：34 条 TRM 代表蛋白
│   ├── 02_At_TRM_annotation.tsv        ★② query 注释表（基因/转录本/TRM符号/长度）
│   ├── 03_trm_blast_pipeline.sh        ★③ BLAST 鉴定流水线（一键运行）
│   ├── 04_summarize_candidates.py      ★④ 候选统计与出图
│   ├── 05_verify_pipeline_logic.py     ★⑤ 流水线逻辑验证（22 条断言）
│   └── 06_verify_bash_structure.py     ★⑥ bash 结构校验器（含默认路径检查）
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
│
└── 05_blast_results/                 ★ 运行产物（首次运行时自动创建）
    └── ...                             原始结果 / 候选 ID / 候选 FASTA / 统计图
```

> 版本控制忽略规则在本仓库根的 `.gitignore`（忽略两个大参考数据文件、BLAST 产物与 `04_temp/`）。

---

## 二、六个核心文件说明

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
| `-evalue` | **1e-5** | 家族鉴定的常规起点 |
| `pident` | **≥30%** | TRM 富含重复序列，identity 天然偏低；卡 50% 会漏真成员 |
| `qcovhsp` | **≥50%** | **主力过滤器**，排除"只匹配到局部重复区"的假阳性 |

**关键理念：不要两个都卡死。** E 值 1e-5 + identity 30% + coverage 50% 三者同时满足，
是"宁可稍宽、但不放水"的平衡点。重复序列家族里 **coverage 比 identity 更可靠**。

### 三档输出

| 档位 | 条件 | 含义 | 用途 |
|---|---|---|---|
| **候选（主）** | C≥50% 且 I≥30% 且 E≤1e-5 | ★ 主交付集 | 建树、家族规模统计 |
| 高可信 | C≥70% 且 E≤1e-10 | 全长直系同源 | 核心成员 |
| 片段待查 | E≤1e-5 但 C<50% | 仅局部匹配，疑似假阳性 | 人工复查 |

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
