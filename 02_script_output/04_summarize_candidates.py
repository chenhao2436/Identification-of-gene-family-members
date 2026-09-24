#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
01_summarize_candidates.py
    BLAST 结果的分布统计与出图（trm_blast_pipeline.sh 阶段 7 调用）。

设计原则
    * 零第三方依赖也能跑完：只用标准库完成全部统计；
    * matplotlib 可选：未安装则跳过 PNG，CSV/文本报告照常产出，不阻塞主流程。

输入
    --blast   out/TRM_vs_Zunla.blast.tsv   15 列 TSV（阶段 2 产出）
    --query   query/TRM_34.fasta           只用于统计 query 总数
    --outdir  out

产出
    01_distribution.tsv      覆盖率/一致性分箱频数表
    01_summary_stats.txt     关键统计与中文判读要点
    01_plots.png             三联图（若 matplotlib 可用）
"""

import argparse
import os
import re
import sys
from collections import Counter, OrderedDict

# 15 列 outfmt 的列序（1-based 对应 awk 的 $1..$15）
COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        "qlen", "slen", "qcovhsp"]

TIER_HIGH_COV = 70.0
TIER_HIGH_E = 1e-10


def parse_args():
    ap = argparse.ArgumentParser(
        description="BLAST 结果分布统计与出图（零依赖降级）")
    ap.add_argument("--blast", required=True, help="15 列 BLAST TSV")
    ap.add_argument("--query", default="", help="query FASTA（统计总条数）")
    ap.add_argument("--outdir", default="out", help="输出目录")
    ap.add_argument("--ident", type=float, default=30.0, help="identity 阈值")
    ap.add_argument("--cov", type=float, default=50.0, help="coverage 阈值")
    ap.add_argument("--evalue", default="1e-5", help="E-value 阈值")
    return ap.parse_args()


def load_hits(path):
    """读取 BLAST TSV -> [dict, ...]。容忍末尾空行与列数不足的行。"""
    hits = []
    skipped = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            f = line.split("\t")
            if len(f) < len(COLS):
                # 有些 BLAST 版本在 outfmt 失败时会输出注释行
                skipped += 1
                continue
            try:
                hits.append({
                    "qseqid": f[0],
                    "sseqid": f[1],
                    "pident": float(f[2]),
                    "length": int(f[3]),
                    "evalue": float(f[10]),
                    "bitscore": float(f[11]),
                    "qlen": int(f[12]),
                    "slen": int(f[13]),
                    "qcovhsp": float(f[14]),
                })
            except ValueError:
                skipped += 1
    return hits, skipped


def gene_of(sseqid):
    """ZLC01G0000010.1 -> ZLC01G0000010（剥离 isoform 后缀）"""
    return re.sub(r"\.\d+$", "", sseqid)


def binned(values, bins, label):
    """分箱频数表 -> [(区间下界, 上界, 计数), ...]"""
    rows = []
    for lo in bins:
        hi = lo + (bins[1] - bins[0] if len(bins) > 1 else 1)
        n = sum(1 for v in values if lo <= v < hi)
        rows.append((lo, hi, n))
    return rows


def cliff(values, gap_min=10.0):
    """
    在已排序的取值里找“最大间隙分割点”。
    返回 (阈值, 间隙大小, 低簇范围, 高簇范围) 或 None。
    用于判读分布是否存在天然断崖。
    """
    if len(values) < 4:
        return None
    u = sorted(set(values))
    if len(u) < 2:
        return None
    best = (0.0, None)
    for a, b in zip(u, u[1:]):
        d = b - a
        if d > best[0]:
            best = (d, (a, b))
    d, pair = best
    if pair is None or d < gap_min:
        return None
    lo_vals = [v for v in values if v <= pair[0]]
    hi_vals = [v for v in values if v >= pair[1]]
    return (pair[1], d, (min(lo_vals), max(lo_vals)), (min(hi_vals), max(hi_vals)))


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    if not os.path.exists(args.blast):
        print("[错误] 找不到 BLAST 结果: %s" % args.blast, file=sys.stderr)
        return 1

    hits, skipped = load_hits(args.blast)
    if not hits:
        print("[错误] BLAST 结果里没有可解析的数据行: %s" % args.blast,
              file=sys.stderr)
        return 1

    n_query_total = None
    if args.query and os.path.exists(args.query):
        with open(args.query, "r", encoding="utf-8", errors="replace") as fh:
            n_query_total = sum(1 for l in fh if l.startswith(">"))

    # ---------------- 基础统计 ----------------
    qcov = [h["qcovhsp"] for h in hits]
    pident = [h["pident"] for h in hits]
    evalues = [h["evalue"] for h in hits]

    queries = sorted({h["qseqid"] for h in hits})
    prot_ids = {h["sseqid"] for h in hits}
    gene_ids = {gene_of(h["sseqid"]) for h in hits}

    # 主候选集（与 pipeline 阈值一致）
    e_thr = float(args.evalue)
    main_hits = [h for h in hits
                 if h["evalue"] <= e_thr
                 and h["pident"] >= args.ident
                 and h["qcovhsp"] >= args.cov]
    main_prot = {h["sseqid"] for h in main_hits}
    main_gene = {gene_of(h["sseqid"]) for h in main_hits}

    high_hits = [h for h in hits
                 if h["qcovhsp"] >= TIER_HIGH_COV and h["evalue"] <= TIER_HIGH_E]
    frag_hits = [h for h in hits
                 if h["evalue"] <= e_thr and h["qcovhsp"] < args.cov]

    # 每 query 候选数
    per_q = Counter(h["qseqid"] for h in main_hits)

    # 零命中 query
    zero_q = []
    if n_query_total:
        all_q = set()
        with open(args.query, "r", encoding="utf-8", errors="replace") as fh:
            for l in fh:
                if l.startswith(">"):
                    all_q.add(l[1:].strip().split()[0])
        zero_q = sorted(all_q - set(queries))

    # ---------------- 分箱 ----------------
    cov_bins = list(range(0, 100, 10))
    id_bins = list(range(0, 100, 5))
    cov_table = binned(qcov, cov_bins, "coverage")
    id_table = binned(pident, id_bins, "identity")

    dist_path = os.path.join(args.outdir, "01_distribution.tsv")
    with open(dist_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("metric\tbin_start\tbin_end\tcount\tpercent\n")
        for lo, hi, n in cov_table:
            fh.write("qcovhsp\t%d\t%d\t%d\t%.2f\n" % (lo, hi, n, 100.0 * n / len(hits)))
        for lo, hi, n in id_table:
            fh.write("pident\t%d\t%d\t%d\t%.2f\n" % (lo, hi, n, 100.0 * n / len(hits)))

    # ---------------- 断崖检测 ----------------
    cov_cliff = cliff(qcov, gap_min=10.0)
    id_cliff = cliff(pident, gap_min=10.0)

    # ---------------- 文本报告 ----------------
    lines = []
    A = lines.append
    A("BLAST 结果分布统计与判读要点")
    A("=" * 70)
    A("输入文件            : %s" % args.blast)
    A("阈值设定            : E<=%s, identity>=%g%%, coverage>=%g%%"
      % (args.evalue, args.ident, args.cov))
    A("")
    A("【总量】")
    A("  HSP 行数           : %d" % len(hits))
    if skipped:
        A("  跳过的不可解析行   : %d" % skipped)
    if n_query_total:
        A("  query 总数         : %d（其中 %d 条至少有一条命中，%d 条零命中）"
          % (n_query_total, len(queries), len(zero_q)))
    else:
        A("  有命中的 query 数  : %d" % len(queries))
    A("  命中到的蛋白 ID 数 : %d（去 isoform 后 %d 个基因）"
      % (len(prot_ids), len(gene_ids)))
    A("")
    A("【主候选集  E<=%s 且 I>=%g%% 且 C>=%g%%】" % (args.evalue, args.ident, args.cov))
    A("  候选 HSP 行数      : %d" % len(main_hits))
    A("  ★ 候选蛋白 ID 数   : %d" % len(main_prot))
    A("  ★ 候选基因数       : %d   <-- 家族规模" % len(main_gene))
    A("")
    A("【分档】")
    A("  高可信 (C>=%g%% 且 E<=%g) : %d HSP" % (TIER_HIGH_COV, TIER_HIGH_E, len(high_hits)))
    A("  主候选 (见上)             : %d HSP" % len(main_hits))
    A("  片段待查 (E 显著但 C<%g%%) : %d HSP  <-- 疑似假阳性，人工复查"
      % (args.cov, len(frag_hits)))
    A("")
    A("【分布特征】")
    A("  coverage  最小/中位/最大 : %.1f / %.1f / %.1f"
      % (min(qcov), sorted(qcov)[len(qcov) // 2], max(qcov)))
    A("  identity  最小/中位/最大 : %.1f / %.1f / %.1f"
      % (min(pident), sorted(pident)[len(pident) // 2], max(pident)))
    A("  覆盖率 <30%% 的 HSP 占比  : %.1f%%  （越高说明低复杂度假阳性越多）"
      % (100.0 * sum(1 for v in qcov if v < 30) / len(qcov)))
    A("")
    if cov_cliff:
        A("  [coverage 断崖] 在 %.1f%% 处最大间隙 %.1f 个百分点"
          % (cov_cliff[0], cov_cliff[1]))
        A("      低簇 %.1f-%.1f%%，高簇 %.1f-%.1f%%"
          % (cov_cliff[2][0], cov_cliff[2][1], cov_cliff[3][0], cov_cliff[3][1]))
        A("      -> 若断崖明显，可考虑把第二轮 coverage 阈值提到 %.0f%%"
          % cov_cliff[0])
    else:
        A("  [coverage 断崖] 未发现 >=10 个百分点的明显间隙")
        A("      -> 分布连续，沿用经验阈值 coverage>=%g%% 即可" % args.cov)
    if id_cliff:
        A("  [identity 断崖] 在 %.1f%% 处最大间隙 %.1f 个百分点"
          % (id_cliff[0], id_cliff[1]))
    else:
        A("  [identity 断崖] 未发现明显间隙")
        A("      -> 重复序列家族常见现象，不必强找拐点，用 coverage 把关")
    A("")
    A("【判读建议】")
    if len(main_gene) < 20:
        A("  * 候选基因数 %d < 20，偏少：可放宽 --evalue 1e-3 重跑" % len(main_gene))
    elif len(main_gene) > 100:
        A("  * 候选基因数 %d > 100，偏多：可收紧到 --cov 70 --evalue 1e-10 重跑"
          % len(main_gene))
    else:
        A("  * 候选基因数 %d 落在预期区间(20-100)，阈值设置合理" % len(main_gene))
    if zero_q:
        A("  * 零命中 query（%d 条）：%s" % (len(zero_q), ", ".join(zero_q)))
        A("      -> 该亚家族可能为拟南芥特有或辣椒中丢失，本身是有意义的结论")
    if frag_hits:
        A("  * 片段待查 %d 条 HSP：建议人工检查其 sseqid 是否只匹配到重复区"
          % len(frag_hits))
    A("  * 优先看 coverage 分布而非 identity：重复序列家族里 coverage 的断崖更干净")

    rpt_path = os.path.join(args.outdir, "01_summary_stats.txt")
    with open(rpt_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")

    # 每 query 候选数（pipeline 已产出 01_candidate_counts.tsv，这里补全零命中者）
    cnt_path = os.path.join(args.outdir, "01_candidate_counts.tsv")
    with open(cnt_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("query\tn_candidates\n")
        for q in queries:
            fh.write("%s\t%d\n" % (q, per_q.get(q, 0)))
        for q in zero_q:
            fh.write("%s\t0\n" % q)

    # ---------------- 出图（可选） ----------------
    png_path = os.path.join(args.outdir, "01_plots.png")
    made_png = False
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

        axes[0].hist(qcov, bins=cov_bins + [100], edgecolor="black")
        axes[0].axvline(args.cov, color="red", linestyle="--",
                        label="threshold %g%%" % args.cov)
        axes[0].set_xlabel("query coverage (%)")
        axes[0].set_ylabel("HSP count")
        axes[0].set_title("Coverage distribution")
        axes[0].legend()

        axes[1].hist(pident, bins=id_bins + [100], edgecolor="black", color="tab:green")
        axes[1].axvline(args.ident, color="red", linestyle="--",
                        label="threshold %g%%" % args.ident)
        axes[1].set_xlabel("identity (%)")
        axes[1].set_ylabel("HSP count")
        axes[1].set_title("Identity distribution")
        axes[1].legend()

        axes[2].scatter(qcov, pident, s=12, alpha=0.4)
        axes[2].axvline(args.cov, color="red", linestyle="--")
        axes[2].axhline(args.ident, color="red", linestyle="--")
        axes[2].set_xlabel("query coverage (%)")
        axes[2].set_ylabel("identity (%)")
        axes[2].set_title("Coverage vs Identity")

        fig.tight_layout()
        fig.savefig(png_path, dpi=150)
        plt.close(fig)
        made_png = True
    except ImportError:
        made_png = False
    except Exception as exc:                                  # noqa: BLE001
        print("[警告] 出图失败（不影响统计结果）: %s" % exc, file=sys.stderr)
        made_png = False

    # ---------------- 屏幕摘要 ----------------
    print("  分布统计完成:")
    print("    %s" % dist_path)
    print("    %s" % rpt_path)
    print("    %s" % cnt_path)
    if made_png:
        print("    %s" % png_path)
    else:
        print("    (未生成 PNG：matplotlib 不可用；可用 01_distribution.tsv 自行绘图)")
    print("    ★ 候选基因数 %d，候选蛋白数 %d"
          % (len(main_gene), len(main_prot)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
