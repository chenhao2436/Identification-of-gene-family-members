#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parse_hmmscan_domtbl.py —— 阶段 06：解析 hmmscan(Pfam-A) 的 --domtblout

只用标准库。

重要：hmmscan 与 hmmsearch 的 domtblout 列序相同（23 列），但语义相反：
    hmmsearch:  target = 蛋白,  query = HMM 模型
    hmmscan  :  target = Pfam 模型, query = 蛋白
因此本脚本按下标取列，与 parse_hmmsearch_domtbl.py 相同。

Pfam 关注域（可配置）
    PF14309  LONGIFOLIA 1/2-like C-terminal domain (DUF4378)
             —— 官方注释明确包含 Protein TRM32 / LONGIFOLIA 1/2，是 TRM 家族核心域
    PF14383  DUF761-associated sequence motif (VARLMGL)
             —— 官方注释与 TRM 无已知关联，本脚本用实测数据检验它是否真的命中

产出（前缀 --prefix）
  <prefix>.all_domains.tsv        全部显著 Pfam 命中（protein × model）
  <prefix>.per_protein.tsv        每蛋白一行：结构域清单 + 关注域 1/0 + TRM域判定
  <prefix>.focus_matrix.tsv       候选 × 关注域 的 1/0 矩阵
  <prefix>.focus_hits.tsv         关注域命中的逐条明细
  <prefix>.stats.txt              统计与判读（含对 PF14383 的实测结论）
"""
import argparse
import os
import re
import sys
from collections import defaultdict

FOCUS = {
    "PF14309": "LONGIFOLIA1/2-like C-terminal domain (DUF4378) — TRM 家族核心域",
    "PF14383": "DUF761-associated sequence motif (VARLMGL) — 与 TRM 关联未知，待检验",
}


def parse_args():
    ap = argparse.ArgumentParser(description="解析 hmmscan Pfam 结果（纯标准库）")
    ap.add_argument("--domtbl", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--dom-eval", type=float, default=1e-5,
                    help="domain 显著性阈值（用 i-Evalue 判定）")
    ap.add_argument("--focus", default="PF14309,PF14383",
                    help="关注的 Pfam 登录号，逗号分隔")
    return ap.parse_args()


def gene_of(x):
    return re.sub(r"\.\d+$", "", x.strip())


def read_domtbl(path):
    rows = []
    bad = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 22:
                bad += 1
                continue
            # hmmscan: f[0]=Pfam model, f[3]=蛋白
            model_full = f[0]
            model = model_full.split(".")[0]        # PF14309.14 -> PF14309
            protein = f[3]
            try:
                rows.append({
                    "model_full": model_full,
                    "model": model,
                    "protein": protein,
                    "gene": gene_of(protein),
                    "tlen_model": int(f[2]),
                    "plen": int(f[5]),
                    "full_e": float(f[6]),
                    "full_score": float(f[7]),
                    "dom_ie": float(f[11]),
                    "dom_score": float(f[12]),
                    "hmm_from": int(f[15]),
                    "hmm_to": int(f[16]),
                    "ali_from": int(f[17]),
                    "ali_to": int(f[18]),
                    "desc": " ".join(f[22:]) if len(f) > 22 else "",
                })
            except (ValueError, IndexError):
                bad += 1
    return rows, bad


def main():
    args = parse_args()
    if not os.path.exists(args.domtbl):
        print("[错误] 找不到 domtblout: %s" % args.domtbl, file=sys.stderr)
        return 1

    odir = os.path.dirname(args.prefix) or "."
    if odir and not os.path.isdir(odir):
        os.makedirs(odir)

    focus = [x.strip() for x in args.focus.split(",") if x.strip()]
    rows, bad = read_domtbl(args.domtbl)
    sig = [r for r in rows if r["dom_ie"] <= args.dom_eval]

    # ---------- 全部显著命中 ----------
    allf = os.path.join(odir, os.path.basename(args.prefix) + ".all_domains.tsv")
    with open(allf, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("protein\tgene\tpfam\tpfam_full\tplen\tmodel_len\tfull_evalue\t"
                 "dom_ievalue\tdom_score\tali_from\tali_to\tali_cov\tpfam_desc\n")
        for r in sorted(sig, key=lambda x: (x["protein"], x["dom_ie"])):
            ali_cov = 100.0 * (r["ali_to"] - r["ali_from"] + 1) / r["plen"] if r["plen"] else 0.0
            fh.write("%s\t%s\t%s\t%s\t%d\t%d\t%.3g\t%.3g\t%.1f\t%d\t%d\t%.1f\t%s\n" % (
                r["protein"], r["gene"], r["model"], r["model_full"], r["plen"],
                r["tlen_model"], r["full_e"], r["dom_ie"], r["dom_score"],
                r["ali_from"], r["ali_to"], ali_cov, r["desc"]))

    # ---------- 每蛋白汇总 ----------
    per_prot = defaultdict(lambda: {"models": set(), "focus": set(), "desc": {}})
    for r in sig:
        d = per_prot[r["protein"]]
        d["models"].add(r["model"])
        d["desc"][r["model"]] = r["desc"]
        if r["model"] in focus:
            d["focus"].add(r["model"])

    ppf = os.path.join(odir, os.path.basename(args.prefix) + ".per_protein.tsv")
    focus_cols = "".join("\t%s" % f for f in focus)
    with open(ppf, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("protein\tgene\tn_pfam_domains\tpfam_models\thas_TRM_core_domain" + focus_cols + "\n")
        for p in sorted(per_prot):
            d = per_prot[p]
            models = sorted(d["models"])
            has_core = "yes" if "PF14309" in d["focus"] else "no"
            cells = "".join("\t%d" % (1 if f in d["focus"] else 0) for f in focus)
            fh.write("%s\t%s\t%d\t%s\t%s%s\n" % (
                p, gene_of(p), len(models), ",".join(models), has_core, cells))

    # ---------- 关注域矩阵（按 gene） ----------
    gene_focus = defaultdict(set)
    for p, d in per_prot.items():
        for f in d["focus"]:
            gene_focus[gene_of(p)].add(f)

    fm = os.path.join(odir, os.path.basename(args.prefix) + ".focus_matrix.tsv")
    with open(fm, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("gene" + "".join("\t%s" % f for f in focus) + "\n")
        for g in sorted(gene_focus):
            fh.write(g + "".join("\t%d" % (1 if f in gene_focus[g] else 0) for f in focus) + "\n")

    # ---------- 关注域明细 ----------
    fhf = os.path.join(odir, os.path.basename(args.prefix) + ".focus_hits.tsv")
    with open(fhf, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("pfam\tprotein\tgene\tfull_evalue\tdom_ievalue\tdom_score\tali_from\tali_to\n")
        for r in sorted([r for r in sig if r["model"] in focus],
                        key=lambda x: (x["model"], x["dom_ie"])):
            fh.write("%s\t%s\t%s\t%.3g\t%.3g\t%.1f\t%d\t%d\n" % (
                r["model"], r["protein"], r["gene"], r["full_e"], r["dom_ie"],
                r["dom_score"], r["ali_from"], r["ali_to"]))

    # ---------- 统计 ----------
    n_prot = len(per_prot)
    stats = os.path.join(odir, os.path.basename(args.prefix) + ".stats.txt")
    with open(stats, "w", encoding="utf-8", newline="\n") as fh:
        A = lambda s: fh.write(s + "\n")
        A("阶段 06：Pfam 结构域核查统计")
        A("=" * 70)
        A("domtblout      : %s" % args.domtbl)
        A("显著性阈值     : i-Evalue <= %g" % args.dom_eval)
        A("")
        A("【规模】")
        A("  输入候选蛋白数           : (见 hmmscan 输入 FASTA)")
        A("  至少命中 1 个 Pfam 域的  : %d" % n_prot)
        A("  显著 domain 命中行数     : %d" % len(sig))
        if bad:
            A("  跳过的不可解析行         : %d" % bad)
        A("")
        A("【关注域实测结果（重点）】")
        for f in focus:
            hit_prots = sorted({r["protein"] for r in sig if r["model"] == f})
            hit_genes = sorted({r["gene"] for r in sig if r["model"] == f})
            A("  %s : %d 蛋白 / %d 基因" % (f, len(hit_prots), len(hit_genes)))
            A("      注释: %s" % FOCUS.get(f, "(未收录说明)"))
            if hit_genes:
                A("      基因: %s" % ", ".join(hit_genes[:20]) +
                  (" ..." if len(hit_genes) > 20 else ""))
        A("")
        # 对 PF14383 的实测结论
        n_8309 = sum(1 for r in sig if r["model"] == "PF14309")
        n_8383 = sum(1 for r in sig if r["model"] == "PF14383")
        A("【结论】")
        if n_8309 > 0:
            A("  * PF14309 在本候选集中有 %d 条显著命中 -> 确认可用作 TRM 家族结构域证据。" % n_8309)
        else:
            A("  * PF14309 在本候选集中 0 命中 -> 需检查 Pfam 版本是否含该模型，")
            A("    或候选是否确实不含该域（这本身是重要信息）。")
        if n_8383 > 0:
            A("  * PF14383 有 %d 条显著命中 -> 实测显示它与本家族确有重叠，" % n_8383)
            A("    可作为辅助标记，但其 Pfam 注释未提及 TRM，解释时需谨慎。")
        else:
            A("  * PF14383 在本候选集中 0 命中 -> 与 InterPro 官方注释一致")
            A("    （该模型描述为 DUF761 相关序列基序，未提及 TRM）。")
            A("    建议不要把 PF14383 作为 TRM 家族判定依据。")
        A("")
        A("【产出】")
        A("  all_domains.tsv   全部显著 Pfam 命中明细")
        A("  per_protein.tsv   每蛋白的结构域清单与 has_TRM_core_domain 标志")
        A("  focus_matrix.tsv  候选 × 关注域 矩阵")
        A("  focus_hits.tsv    关注域逐条命中")

    print("  显著 domain 命中: %d 行, 覆盖 %d 个蛋白" % (len(sig), n_prot))
    for f in focus:
        n = sum(1 for r in sig if r["model"] == f)
        print("    %s : %d 条" % (f, n))
    print("    %s" % ppf)
    print("    %s" % stats)
    return 0


if __name__ == "__main__":
    sys.exit(main())
