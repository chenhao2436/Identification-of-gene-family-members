#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
merge_blast_hmm.py —— 阶段 05：合并 BLAST 与 HMM 候选，按 gene 去重并标注证据来源

只用标准库。

证据标签（你要求的四类）
  BLAST+HMM        两种方法都支持
  BLAST-only       仅 BLAST 支持（HMM 未显著命中）
  HMM-only         仅 HMM 显著命中（BLAST 未显著命中）
  HMM-inclusive-only  仅达 HMM 包含阈值(i-Evalue<=1e-2)，未达全序列 strict 阈值

主统计单位: gene（剥离 isoform 版本号）

输入
  --blast-gene   05_blast_results/TRM_candidates_gene_level.txt
  --blast-prots  05_blast_results/TRM_candidates_protein_level.txt
  --blast-main   05_blast_results/TRM_candidates_main.tsv（15 列 BLAST outfmt）
  --hmm-gene     04_hmmsearch/TRM.hmmsearch.gene_level.tsv
  --hmm-gene-strict    ...gene_ids.txt
  --hmm-gene-incl      ...gene_ids_inclusive.txt
  --anno         02_script_output/02_At_TRM_annotation.tsv（把 At 转录本映射回 At 基因）

产出（前缀由 --prefix 指定）
  <prefix>.merged.tsv        合并后的完整候选表（每行一个基因）
  <prefix>.evidence_counts.tsv  四类证据标签计数
  <prefix>.by_evidence.tsv   按证据标签分组列表
  <prefix>.stats.txt         统计与判读要点
"""
import argparse
import os
import re
import sys
from collections import defaultdict


def parse_args():
    ap = argparse.ArgumentParser(description="合并 BLAST 与 HMM 候选（纯标准库）")
    ap.add_argument("--blast-gene", required=True)
    ap.add_argument("--blast-prots", default="")
    ap.add_argument("--blast-main", default="")
    ap.add_argument("--hmm-gene", default="")
    ap.add_argument("--hmm-gene-strict", default="")
    ap.add_argument("--hmm-gene-incl", default="")
    ap.add_argument("--anno", default="")
    ap.add_argument("--prefix", required=True)
    return ap.parse_args()


def gene_of(x):
    return re.sub(r"\.\d+$", "", x.strip())


def read_lines(path):
    if not path or not os.path.exists(path):
        return []
    out = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                out.append(line)
    return out


def read_anno(path):
    """At 转录本 -> (At基因, TRM符号)"""
    m = {}
    if not path or not os.path.exists(path):
        return m
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i == 0:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4:
                tid = f[1].strip()
                at_gene = f[0].strip()
                sym = f[3].strip()
                m[tid] = (at_gene, sym)
    return m


def read_blast_main(path):
    """15 列 BLAST outfmt -> gene -> 最佳指标 + 支持的 At query 集合"""
    per_gene = defaultdict(lambda: {"pident": 0.0, "qcov": 0.0, "evalue": None,
                                    "bitscore": 0.0, "queries": set(), "prots": set()})
    if not path or not os.path.exists(path):
        return per_gene
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 15:
                continue
            try:
                q, s = f[0], f[1]
                pid = float(f[2]); ev = float(f[10]); bs = float(f[11]); cov = float(f[14])
            except ValueError:
                continue
            g = gene_of(s)
            d = per_gene[g]
            d["queries"].add(q)
            d["prots"].add(s)
            if pid > d["pident"]:
                d["pident"] = pid
            if cov > d["qcov"]:
                d["qcov"] = cov
            if bs > d["bitscore"]:
                d["bitscore"] = bs
            if d["evalue"] is None or ev < d["evalue"]:
                d["evalue"] = ev
    return per_gene


def read_hmm_gene_level(path):
    """TRM.hmmsearch.gene_level.tsv -> gene -> 指标"""
    per_gene = {}
    if not path or not os.path.exists(path):
        return per_gene
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i == 0:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            g = f[0].strip()
            try:
                per_gene[g] = {
                    "n_proteins": int(f[1]),
                    "best_protein": f[2].strip(),
                    "n_models": int(f[3]),
                    "full_e": float(f[4]),
                    "dom_score": float(f[5]),
                    "hmm_cov": float(f[6]),
                    "ali_cov": float(f[7]),
                    "significance": f[8].strip(),
                }
            except ValueError:
                continue
    return per_gene


def main():
    args = parse_args()
    odir = os.path.dirname(args.prefix) or "."
    if odir and not os.path.isdir(odir):
        os.makedirs(odir)

    blast_genes = {gene_of(g) for g in read_lines(args.blast_gene)}
    blast_prots = {p.strip() for p in read_lines(args.blast_prots)}
    hmm_strict = {gene_of(g) for g in read_lines(args.hmm_gene_strict)}
    hmm_incl = {gene_of(g) for g in read_lines(args.hmm_gene_incl)}
    blast_metrics = read_blast_main(args.blast_main)
    hmm_metrics = read_hmm_gene_level(args.hmm_gene)
    anno = read_anno(args.anno)

    # 若 strict 列表缺失，退化为从 gene_level 表里按 significance 判断
    if not hmm_strict and hmm_metrics:
        hmm_strict = {g for g, d in hmm_metrics.items() if d["significance"] == "strict"}
    if not hmm_incl and hmm_metrics:
        hmm_incl = {g for g, d in hmm_metrics.items()
                    if d["significance"] in ("strict", "inclusive_only")}

    only_blast = blast_genes - hmm_incl
    only_hmm_strict = hmm_strict - blast_genes
    only_hmm_incl = (hmm_incl - hmm_strict) - blast_genes
    both = blast_genes & hmm_strict
    union = blast_genes | hmm_incl

    def evidence_of(g):
        if g in both:
            return "BLAST+HMM"
        if g in blast_genes:
            return "BLAST-only"
        if g in hmm_strict:
            return "HMM-only"
        if g in hmm_incl:
            return "HMM-inclusive-only"
        return "?"

    merged = os.path.join(odir, os.path.basename(args.prefix) + ".merged.tsv")
    rows = []
    with open(merged, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("gene\tevidence\tn_blast_proteins\tblast_best_pident\tblast_best_qcov\t"
                 "blast_min_evalue\tblast_best_bitscore\tn_hmm_proteins\thmm_best_full_evalue\t"
                 "hmm_best_dom_score\thmm_best_hmm_cov\thmm_significance\t"
                 "n_at_queries\tat_genes\tat_symbols\tblast_proteins\thmm_best_protein\tmanual_review\n")
        for g in sorted(union):
            b = blast_metrics.get(g, {"pident": 0.0, "qcov": 0.0, "evalue": None,
                                      "bitscore": 0.0, "queries": set(), "prots": set()})
            h = hmm_metrics.get(g, {})
            ev = "%.3g" % b["evalue"] if b["evalue"] is not None else "NA"
            # At query -> At 基因 / TRM 符号
            at_genes, at_syms = set(), set()
            for q in sorted(b["queries"]):
                if q in anno:
                    at_genes.add(anno[q][0])
                    if anno[q][1]:
                        at_syms.add(anno[q][1])
            # 待人工确认标记：证据单薄或仅达包含阈值
            flags = []
            if g in only_hmm_incl:
                flags.append("仅达HMM包含阈值")
            if g in only_blast:
                flags.append("仅BLAST支持")
            if g in only_hmm_strict:
                flags.append("仅HMM支持")
            if not b["queries"] and not h:
                flags.append("无有效证据")
            if h and h.get("hmm_cov", 0) < 50:
                flags.append("HMM覆盖不足50%")
            if b["qcov"] and b["qcov"] < 50:
                flags.append("BLAST覆盖不足50%")
            mr = ";".join(flags) if flags else "-"
            rows.append((g, evidence_of(g), mr))
            fh.write("%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n" % (
                g, evidence_of(g),
                len(b["prots"]),
                ("%.1f" % b["pident"]) if b["pident"] else "NA",
                ("%.0f" % b["qcov"]) if b["qcov"] else "NA",
                ev,
                ("%.1f" % b["bitscore"]) if b["bitscore"] else "NA",
                h.get("n_proteins", "NA"),
                ("%.3g" % h["full_e"]) if h else "NA",
                ("%.1f" % h["dom_score"]) if h else "NA",
                ("%.1f" % h["hmm_cov"]) if h else "NA",
                h.get("significance", "NA"),
                len(b["queries"]),
                ",".join(sorted(at_genes)) if at_genes else "NA",
                ",".join(sorted(at_syms)) if at_syms else "NA",
                ",".join(sorted(b["prots"])) if b["prots"] else "NA",
                h.get("best_protein", "NA"),
                mr))

    # ---------- 证据计数 ----------
    counts = os.path.join(odir, os.path.basename(args.prefix) + ".evidence_counts.tsv")
    tally = defaultdict(int)
    for g in union:
        tally[evidence_of(g)] += 1
    with open(counts, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("evidence\tn_genes\n")
        for k in ("BLAST+HMM", "BLAST-only", "HMM-only", "HMM-inclusive-only"):
            fh.write("%s\t%d\n" % (k, tally.get(k, 0)))
        fh.write("TOTAL\t%d\n" % len(union))

    # ---------- 按证据分组 ----------
    byev = os.path.join(odir, os.path.basename(args.prefix) + ".by_evidence.tsv")
    groups = defaultdict(list)
    for g, ev, mr in rows:
        groups[ev].append(g)
    with open(byev, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("evidence\tgene\n")
        for ev in ("BLAST+HMM", "BLAST-only", "HMM-only", "HMM-inclusive-only"):
            for g in groups.get(ev, []):
                fh.write("%s\t%s\n" % (ev, g))

    # ---------- 统计报告 ----------
    stats = os.path.join(odir, os.path.basename(args.prefix) + ".stats.txt")
    with open(stats, "w", encoding="utf-8", newline="\n") as fh:
        A = lambda s: fh.write(s + "\n")
        A("阶段 05：BLAST + HMM 候选合并统计")
        A("=" * 70)
        A("主统计单位: gene（isoform 版本号已剥离）")
        A("")
        A("【各自规模】")
        A("  BLAST 候选基因            : %d" % len(blast_genes))
        A("  HMM strict 候选基因        : %d" % len(hmm_strict))
        A("  HMM inclusive 候选基因     : %d" % len(hmm_incl))
        A("")
        A("【合并结果（并集）】")
        A("  ★ 合并候选基因总数         : %d" % len(union))
        A("    BLAST+HMM 双重支持       : %d" % tally.get("BLAST+HMM", 0))
        A("    BLAST-only               : %d" % tally.get("BLAST-only", 0))
        A("    HMM-only (strict)        : %d" % tally.get("HMM-only", 0))
        A("    HMM-inclusive-only       : %d" % tally.get("HMM-inclusive-only", 0))
        A("")
        A("【值得注意】")
        if tally.get("HMM-only", 0) > 0:
            A("  * 有 %d 个基因仅 HMM 命中、BLAST 未命中 ——" % tally.get("HMM-only", 0))
            A("    profile HMM 对远缘同源比 BLAST 敏感，这类是 BLAST 漏检的候选，")
            A("    正是引入 HMMER 的价值所在，不应直接丢弃。")
        if tally.get("BLAST-only", 0) > 0:
            A("  * 有 %d 个基因仅 BLAST 命中 —— 需在阶段 06 用结构域核查判断" % tally.get("BLAST-only", 0))
            A("    是否含 TRM 特征结构域；无结构域支持者标为待人工确认。")
        if tally.get("HMM-inclusive-only", 0) > 0:
            A("  * 有 %d 个基因仅达 HMM 包含阈值(i-Evalue<=1e-2)，属最宽松档，"
              % tally.get("HMM-inclusive-only", 0))
            A("    默认保留并标记，不删除。")
        A("")
        A("【产出】")
        A("  merged.tsv           每基因一行，含双方法指标与证据标签")
        A("  evidence_counts.tsv  四类证据计数")
        A("  by_evidence.tsv      按证据分组的基因清单")
        A("")
        A("说明: 本步骤不做任何删除。可疑候选在 merged.tsv 的 manual_review 列标记，")
        A("      留给阶段 06/07 结合结构域证据与人工判断。")

    print("  BLAST 候选基因 : %d" % len(blast_genes))
    print("  HMM strict     : %d" % len(hmm_strict))
    print("  HMM inclusive  : %d" % len(hmm_incl))
    print("  ★ 合并并集     : %d" % len(union))
    for k in ("BLAST+HMM", "BLAST-only", "HMM-only", "HMM-inclusive-only"):
        print("      %-20s %d" % (k, tally.get(k, 0)))
    print("    %s" % merged)
    return 0


if __name__ == "__main__":
    sys.exit(main())
