# -*- coding: utf-8 -*-
"""
用真实数据核对 Zunla 的 ID 映射关系（无需 bash/HMMER）。

为什么先做这个
    阶段 06/07 要用 GFF3 给候选基因补坐标，并与蛋白库对照。
    如果 ID 体系对不上，服务器上会白跑一轮才在阶段 01 报警。
    这里用真实 GFF3 + 真实候选列表先验证。

验证项
    1. 候选基因 ID -> GFF3 gene 特征 的命中率
    2. 候选蛋白 ID -> Zunla 蛋白库 的命中率
    3. GFF3 mRNA ID 与蛋白库 ID 的交集（确认三级命名一致）
    4. 候选基因的染色体分布与坐标可提取性
    5. Zunla 蛋白库的 ID 唯一性与空序列
"""
import os
import re
import sys
from collections import Counter

PROJECT = (r"E:\华为\deepseek harness\gene_famaily_analysis"
           r"\Identification_of_gene_family_members")
BLAST = os.path.join(PROJECT, "05_blast_results")
PEP = os.path.join(PROJECT, "01_At_TRM_query", "Zunla_Canz.pep.fa")
GFF3 = r"E:\华为\fsdownload\Chr11-12_region_BAM\PP64_ref\Canz.genome.gff3"


def read_lines(p):
    with open(p, "r", encoding="utf-8", errors="replace") as fh:
        return [l.strip() for l in fh if l.strip()]


def main():
    fails = []
    n = [0]

    def chk(cond, label, detail=""):
        n[0] += 1
        mark = "OK  " if cond else "FAIL"
        print("  [%s] %s%s" % (mark, label, ("  " + detail) if detail else ""))
        if not cond:
            fails.append(label)

    print("=" * 74)
    print("Zunla ID 映射真实数据核对")
    print("=" * 74)

    for p in (PEP, GFF3,
              os.path.join(BLAST, "TRM_candidates_gene_level.txt"),
              os.path.join(BLAST, "TRM_candidates_protein_level.txt")):
        if not os.path.exists(p):
            print("缺少文件: %s" % p)
            return 2

    cand_genes = read_lines(os.path.join(BLAST, "TRM_candidates_gene_level.txt"))
    cand_prots = read_lines(os.path.join(BLAST, "TRM_candidates_protein_level.txt"))
    print("\n候选: %d 基因 / %d 蛋白" % (len(cand_genes), len(cand_prots)))

    # ---------------- 蛋白库 ----------------
    print("\n[1] Zunla 蛋白库")
    pep_ids = []
    empty = 0
    cur_len = None
    with open(PEP, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_len == 0:
                    empty += 1
                pep_ids.append(line[1:].strip().split()[0])
                cur_len = 0
            elif cur_len is not None:
                cur_len += len(line.strip())
    if cur_len == 0:
        empty += 1
    n_gt = len(pep_ids)
    n_uniq = len(set(pep_ids))
    chk(n_gt > 0, "蛋白库有序列", "%d 条" % n_gt)
    chk(n_gt == n_uniq, "蛋白库 ID 唯一",
        "重复 %d 个" % (n_gt - n_uniq) if n_gt != n_uniq else "唯一")
    chk(empty == 0, "无空序列", "空序列 %d 条" % empty)
    print("       序列 ID 样例: %s" % ", ".join(pep_ids[:3]))

    # ---------------- 候选蛋白 -> 蛋白库 ----------------
    print("\n[2] 候选蛋白 -> 蛋白库")
    pep_set = set(pep_ids)
    hit = [p for p in cand_prots if p in pep_set]
    chk(len(hit) == len(cand_prots), "全部候选蛋白都能在库中找到",
        "%d/%d" % (len(hit), len(cand_prots)))
    missing = [p for p in cand_prots if p not in pep_set]
    if missing:
        print("       未命中: %s" % ", ".join(missing[:8]))

    # ---------------- GFF3 ----------------
    print("\n[3] GFF3 解析（66 MB，单遍扫描）")
    gene_c = {}
    mrna_ids = set()
    n_gene = n_mrna = 0
    with open(GFF3, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            t = f[2]
            if t == "gene":
                m = re.search(r"(?:^|;)ID=([^;]+)", f[8])
                if m:
                    n_gene += 1
                    gene_c[m.group(1)] = (f[0], f[3], f[4], f[6])
            elif t == "mRNA":
                m = re.search(r"(?:^|;)ID=([^;]+)", f[8])
                if m:
                    n_mrna += 1
                    mrna_ids.add(m.group(1))
    chk(n_gene > 0, "GFF3 有 gene 特征", "%d 个" % n_gene)
    chk(n_mrna > 0, "GFF3 有 mRNA 特征", "%d 个" % n_mrna)

    print("\n[4] 候选基因 -> GFF3 gene 特征")
    g_hit = [g for g in cand_genes if g in gene_c]
    chk(len(g_hit) == len(cand_genes), "全部候选基因都能映射到 gene 特征",
        "%d/%d" % (len(g_hit), len(cand_genes)))
    gmiss = [g for g in cand_genes if g not in gene_c]
    if gmiss:
        print("       未命中: %s" % ", ".join(gmiss[:8]))

    print("\n[5] 候选基因的坐标与染色体分布")
    chrom = Counter(gene_c[g][0] for g in g_hit)
    print("       染色体分布: %s" % dict(sorted(chrom.items())))
    chk(len(chrom) > 0, "能提取到坐标（阶段06/07 需要）",
        "覆盖 %d 条染色体" % len(chrom))
    if g_hit:
        g0 = sorted(g_hit)[0]
        print("       样例: %s -> %s" % (g0, gene_c[g0]))

    # ---------------- 三级命名一致性 ----------------
    print("\n[6] GFF3 mRNA 与蛋白库 ID 的交集（三级命名一致性）")
    common = mrna_ids & pep_set
    chk(len(common) > 0, "mRNA ID 与蛋白 ID 有交集",
        "%d 个" % len(common))
    if common:
        c0 = sorted(common)[0]
        print("       样例: mRNA=%s  蛋白=%s" % (c0, c0))

    # 候选蛋白是否都在 mRNA 集合里（即 protein == mRNA ID）
    p_in_mrna = [p for p in cand_prots if p in mrna_ids]
    print("       候选蛋白中同时是 GFF3 mRNA 的: %d/%d"
          % (len(p_in_mrna), len(cand_prots)))
    chk(len(p_in_mrna) == len(cand_prots),
        "候选蛋白 ID 与 GFF3 mRNA ID 一一对应（gene/transcript/protein 三级一致）")

    print("\n" + "=" * 74)
    print("断言数: %d   失败: %d" % (n[0], len(fails)))
    if fails:
        print("失败项:")
        for f in fails:
            print("  - %s" % f)
        return 1
    print("全部通过 ✅  阶段 06/07 的 ID 映射与坐标提取有保障")
    return 0


if __name__ == "__main__":
    sys.exit(main())
