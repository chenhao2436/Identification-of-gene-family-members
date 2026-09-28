# -*- coding: utf-8 -*-
"""
audit_blast_coverage.py —— BLAST 候选的 query 覆盖审计（离线，不依赖服务器）

为什么做这个
    1. 最终 CaTRM 列表要给出"每个成员的证据来源与筛选理由"，
       而 BLAST 侧的线索是"它对应哪个 AtTRM query"。
       先离线把这张对应表建出来，阶段 05/07 一到手就能直接用。
    2. 34 条 query 里若有的没有被任何候选覆盖，说明 BLAST 漏了，
       这正是 HMMER 该补回来的部分——提前标出来，便于判断 HMMER 是否生效。
    3. 某条 query 若贡献了异常多的候选，需警惕假阳性。

输入（均为 BLAST 阶段已有产物）
    05_blast_results/01_candidate_summary.tsv      query/At基因/TRM符号/Zunla基因/I/E/bitscore/C/档位
    05_blast_results/TRM_candidates_main.tsv       15 列原始 outfmt
    05_blast_results/01_zero_hit_queries.txt       零命中 query
    02_script_output/02_At_TRM_annotation.tsv      34 条 query 全集

产出
    06_hmmer_results/05_merge/blast_query_coverage.tsv   每 At 基因 -> 候选 Zunla 基因
    06_hmmer_results/05_merge/blast_query_coverage.txt   可读报告与判读
"""
import os
import re
import sys
from collections import defaultdict

PROJECT = (r"E:\华为\deepseek harness\gene_famaily_analysis"
           r"\Identification_of_gene_family_members")
BLAST = os.path.join(PROJECT, "05_blast_results")
ANNO = os.path.join(PROJECT, "02_script_output", "02_At_TRM_annotation.tsv")
OUTDIR = os.path.join(PROJECT, "06_hmmer_results", "05_merge")


def read_tsv(path):
    """读带表头的 TSV -> [dict]"""
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        header = None
        for i, line in enumerate(fh):
            f = line.rstrip("\n").split("\t")
            if i == 0:
                header = f
                continue
            if header and len(f) == len(header):
                rows.append(dict(zip(header, f)))
    return rows


def read_tsv_positional(path):
    """
    读 TSV -> [[字段,...]]，**自动跳过表头**（若存在）。

    背景: 01_candidate_summary.tsv 早先因流水线缺陷丢失了表头（已修复）。
    为兼容修复前后两种产物，这里按首个字段是否为已知列名来判断表头。
    """
    header_tokens = {"query", "qseqid", "gene", "sseqid", "protein", "at_gene"}
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            line = line.rstrip("\n")
            if not line:
                continue
            f = line.split("\t")
            if i == 0 and f[0].strip().lower() in header_tokens:
                continue          # 跳过表头
            rows.append(f)
    return rows


# 01_candidate_summary.tsv 的列序（无表头）
# query trm_gene trm_symbol zunla_protein zunla_gene pident evalue bitscore qcovhsp tier
SUM_COLS = ["query", "trm_gene", "trm_symbol", "zunla_protein", "zunla_gene",
            "pident", "evalue", "bitscore", "qcovhsp", "tier"]


def read_lines(p):
    if not os.path.exists(p):
        return []
    with open(p, "r", encoding="utf-8", errors="replace") as fh:
        return [l.strip() for l in fh if l.strip()]


def main():
    for p in (ANNO,
              os.path.join(BLAST, "01_candidate_summary.tsv"),
              os.path.join(BLAST, "TRM_candidates_main.tsv")):
        if not os.path.exists(p):
            print("缺少输入: %s" % p)
            return 2

    # 34 条 query 全集（含未命中的）
    all_q = {}
    for r in read_tsv(ANNO):
        tid = r.get("representative_transcript", "").strip()
        if tid:
            all_q[tid] = (r.get("gene", ""), r.get("symbols", ""), r.get("length_aa", ""))

    zero_q = set(read_lines(os.path.join(BLAST, "01_zero_hit_queries.txt")))

    # 候选汇总表：主候选(query -> Zunla 候选)。无表头，按位置列解析。
    summ_rows = read_tsv_positional(os.path.join(BLAST, "01_candidate_summary.tsv"))
    summ = []
    bad_rows = 0
    for f in summ_rows:
        if len(f) != len(SUM_COLS):
            bad_rows += 1
            continue
        summ.append(dict(zip(SUM_COLS, f)))
    if not summ:
        print("[错误] 01_candidate_summary.tsv 解析为空（列数不匹配？）")
        return 1

    # 15 列原始 outfmt，同样无表头
    main_rows = read_tsv_positional(os.path.join(BLAST, "TRM_candidates_main.tsv"))

    hit_any = defaultdict(set)   # query -> 所有命中过的 Zunla 蛋白（不论是否过阈值）
    for f in main_rows:
        if len(f) >= 2:
            hit_any[f[0]].add(f[1])

    cand_by_q = defaultdict(set)  # query -> 主候选 Zunla 基因
    for r in summ:
        q = r["query"].strip()
        g = r["zunla_gene"].strip()
        if q and g and g != "NA":
            cand_by_q[q].add(g)

    at_gene_of = {tid: v[0] for tid, v in all_q.items()}
    sym_of = {tid: v[1] for tid, v in all_q.items()}

    # 汇总到 At 基因级
    cand_by_at = defaultdict(set)
    q_of_at = defaultdict(set)
    for tid, genes in cand_by_q.items():
        ag = at_gene_of.get(tid, tid)
        cand_by_at[ag] |= genes
        q_of_at[ag].add(tid)

    # ---------------- 报告 ----------------
    lines = []
    A = lines.append
    A("BLAST 候选的 query 覆盖审计")
    A("=" * 74)
    A("输入 : %s" % os.path.join(BLAST, "01_candidate_summary.tsv"))
    A("")
    A("【规模】")
    A("  AtTRM query 总数        : %d" % len(all_q))
    A("  至少 R1 命中过的 query  : %d" % len(hit_any))
    A("  零 R1 命中 (E<=1e-5) 的 : %d  -> %s"
      % (len(zero_q), ", ".join(sorted(zero_q)) if zero_q else "无"))
    A("  进入主候选的 query      : %d" % len(cand_by_q))
    A("  主候选覆盖的 At 基因    : %d / %d" % (len(cand_by_at), len(all_q)))
    A("")

    # 未被主候选覆盖的 query
    uncovered = [t for t in sorted(all_q) if t not in cand_by_q]
    A("【未进入主候选的 query（%d 条）】" % len(uncovered))
    A("  说明: 这些 query 的候选要么 E 值不达标，要么被 I/C 阈值刷掉。")
    A("        它们是最需要 HMMER 补检的部分。")
    for t in uncovered:
        ag, sym, _ = all_q[t]
        n_any = len(hit_any.get(t, ()))
        reason = "R1 完全无命中" if n_any == 0 else "命中 %d 条但未过阈值" % n_any
        A("  %-14s %-10s %-16s %s" % (t, ag, sym, reason))
    A("")

    # 零命中（E<=1e-5 都没有）
    A("【E<=1e-5 内完全无命中的 query】")
    if zero_q:
        for t in sorted(zero_q):
            ag, sym, _ = all_q.get(t, ("?", "?", "?"))
            A("  %-14s %-10s %-16s -> 该亚家族可能拟南芥特有或辣椒中丢失" % (t, ag, sym))
    else:
        A("  无")
    A("")

    # 候选数分布：找异常
    A("【每个 At 基因贡献的候选数（降序）】")
    ranked = sorted(cand_by_at.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    for ag, genes in ranked:
        sym = ""
        for t in q_of_at[ag]:
            sym = sym_of.get(t, "")
            break
        A("  %-12s %-18s %3d 个候选" % (ag, sym, len(genes)))
    A("")

    if ranked:
        top_ag, top_genes = ranked[0]
        n = len(top_genes)
        sizes = sorted(len(genes) for genes in cand_by_at.values())
        median = sizes[len(sizes) // 2]
        A("【异常检查】")
        A("  候选数最多的 At 基因: %s (%d 个候选)，中位数 %d"
          % (top_ag, n, median))
        if n >= 3 * max(median, 1):
            A("  [注意] 该基因贡献的候选数远超中位数，建议在阶段 06 重点核查其")
            A("         候选的 PF14309 结构域是否成立——重复序列家族常见单点")
            A("         多命中的假阳性。")
        else:
            A("  候选数分布未见明显异常离群。")
        A("")

    # 写出 At 基因 -> 候选 Zunla 基因 对应表
    os.makedirs(OUTDIR, exist_ok=True)
    out_tsv = os.path.join(OUTDIR, "blast_query_coverage.tsv")
    with open(out_tsv, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("at_gene\tat_symbol\tat_queries\tn_candidates\tzunla_genes\n")
        for ag, genes in sorted(cand_by_at.items()):
            syms = sorted({sym_of.get(t, "") for t in q_of_at[ag]} - {""})
            fh.write("%s\t%s\t%s\t%d\t%s\n" % (
                ag, ",".join(syms), ",".join(sorted(q_of_at[ag])),
                len(genes), ",".join(sorted(genes))))

    rpt = os.path.join(OUTDIR, "blast_query_coverage.txt")
    with open(rpt, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")

    print("\n".join(lines))
    print("=" * 74)
    print("产出:")
    print("  %s" % out_tsv)
    print("  %s" % rpt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
