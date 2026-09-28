# TRM 家族结构域注释 — 权威参考（供论文 Methods 与结果解释使用）

数据来源：InterPro API（EBI），逐条核对。查询日期见文件末尾。
用途：确定哪些结构域可以作为 TRM 家族的判定依据，避免误用。

---

## 1. 结论速览

| | 结构域 | 数据库条目 | 长度 | 位置 | 能否作为家族判定依据 |
|---|---|---|---|---|---|
| C 端 | LONGIFOLIA 1/2-like C-terminal domain | **Pfam PF14309** = InterPro **IPR025486** | ~180 aa | C 端 | ✅ **可以**（信号 1e-30 量级，全部成员命中） |
| N 端 | DUF761-associated sequence motif (VARLMGL) | **Pfam PF14383** = InterPro **IPR032795** | 22–31 aa | N 端 | ⚠️ **不建议单独使用**（短基序、信号弱、多数成员未命中） |
| N 端 | Protein of unknown function DUF3741 | **Pfam PF12552** = InterPro **IPR022212** | ~50 aa | N 端 | ⚠️ 未知功能域，仅作位置参考 |
| 全长 | Protein LONGIFOLIA 1/2 | **PANTHER PTHR31680** = InterPro **IPR033334** | 近全长 | 1 → C 端 | ✅ **最佳家族标记**（GO:0051513） |

**关键点：Pfam 在 TRM 家族中没有全长家族模型。**
Pfam 侧只有 C 端域 PF14309；因此 **hmmscan(Pfam) 只能验证"C 端域存在"，
不能证明"该蛋白是 TRM 家族成员"**。全长家族归属需靠
① 自建 TRM.hmm（本流程阶段 02–04）或 ② PANTHER PTHR31680。

---

## 2. 逐条实测证据

| 蛋白 | UniProt | 长度 | PANTHER | PF14309 (C 端) | PF14383 (N 端) |
|---|---|---|---|---|---|
| LNG2 / TRM1 (At3g02170) | Q9S823 | 905 | PTHR31680**SF15**（PROTEIN LONGIFOLIA 2），覆盖 **1–902**，score 0 | **722–901，4.1e-30** | 未命中 |
| LNG1 / TRM2 (At5g15580) | Q9LF24 | 927 | PTHR31680**SF22**（PROTEIN LONGIFOLIA 1），覆盖 **1–925**，score 0 | **721–904，7.8e-35** | 未命中 |
| TRM4 (At1g74160) | Q0WNQ5 | 1025 | PTHR31680**SF4**，覆盖 **1–1022**，score 0 | **820–1002，1.8e-33** | 274–295（22 aa），**score 0.01** |
| 辣椒同源蛋白 (Capsicum) | A0A022S2L0 | 805 | — | **654–797，1.7e-30** | 81–111（31 aa），3.1e-10 |

> 三条已复核的拟南芥成员中，**PF14309 全部命中（3/3）**，
> 而 **PF14383 仅 1/3 命中且 score 仅 0.01**。

---

## 3. 结构域排布（以 TRM4 为例）

```
N 端                                                     C 端
|--------------------------------------------------------------|
   PF14383        PF12552/DUF3741                    PF14309
   (22 aa)        (~50 aa)                           (~180 aa)
   274-295        (PF14383 定位于 DUF3741 的 N 端)    820-1002
   └──────── N 端模块 ────────┘                      └─ 家族特征域 ─┘

   └──────────────── PTHR31680（PANTHER 全长家族模型, 1-1022）─────────┘
```

---

## 4. 参考条目官方描述（原文要点）

- **IPR025486 / PF14309**："This entry represents the C-terminal domain in
  Centrosome-associated protein 350 (Cep350), **Protein TRM32** and
  **Protein LONGIFOLIA 1/2** found in eukaryotes. This domain contains a helical fold."
  并引用了 TRM1–TON1 那篇 *Plant Cell* 2012（PMID 22286137）与 LNG1/LNG2
  的 *Development* 2006（PMID 17038516）。

- **IPR032795 / PF14383**："This **short domain** is found frequently at the
  **N terminus** of domain DUF3741 (IPR022212). It contains the sequence motif VARLMGL."

- **IPR022212 / PF12552 (DUF3741)**："This domain is found in plant proteins,
  and is approximately 50 amino acids in length."

- **IPR033334 / PTHR31680**："Protein LONGIFOLIA 1/2"，GO:0051513
  （regulation of monopolar cell growth）。

---

## 5. 对判读的具体建议

1. **`conservative` 层用 PF14309 作必要条件** —— 依据充分（3/3 命中、1e-30 量级）。
2. **不要把 PF14383 缺失当作否证** —— 它多数成员不命中，且 Pfam 对该条目有
   已知的重复命中抑制规则。
3. **若候选只支持 C 端域、无 HMM 全长支持，应标为 `review`**，不要直接计入家族。
4. **论文中若要写"家族成员确认"**，应同时报告：
   - 自建 TRM.hmm 的全长命中（本流程阶段 04）
   - PF14309 结构域支持（本流程阶段 06）
   - 可选：PANTHER PTHR31680 支持
5. **N 端模块若需描述**，用 **PF12552/DUF3741**，不要只写 PF14383。

---

## 6. 如何复核本文件的数字

```bash
# 取某 AtTRM 的 UniProt 编号
curl -s "https://rest.uniprot.org/uniprotkb/search?query=gene:AT3G02170+AND+organism_id:3702&format=json&fields=accession" | head

# 取该蛋白在全部数据库中的注释（含 Panthera/Pfam 位置与分数）
curl -s "https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/Q9S823/?page_size=30"
```

> 核对日期：2026-09（本轮）。Pfam / InterPro 会定期更新，正式发表前建议重查一次。

## 7. 相关文献（可写入 Methods 引用）

| PMID | 期刊/年份 | 标题要点 |
|---|---|---|
| 22286137 | Plant Cell, 2012 | TRM1–TON1 互作，揭示植物皮层微管阵列与真核中心体的共同招募网络 |
| 17038516 | Development, 2006 | LONGIFOLIA1/2 调控拟南芥纵向细胞伸长 |
