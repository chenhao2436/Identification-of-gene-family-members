#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parse_hmmsearch_domtbl.py —— 解析 hmmsearch 的 --domtblout 输出

只用标准库，服务器无需装 Biopython/pandas。

输入（HMMER 3.x domtblout 列序，共 23 列）:
   1 target name            2 target accession     3 tlen
   4 query name             5 query accession      6 qlen
   7 full-seq E-value       8 full-seq score       9 full-seq bias
  10 # of domains
  11 c-Evalue              12 i-Evalue            13 domain score
  14 domain bias           15 hmm from             16 hmm to
  17 ali from              18 ali to               19 env from
  20 env to                21 acc                  22 description
  23 (空)

产出:
  <prefix>.parsed.tsv            全部蛋白×模型命中（每行一个 domain）
  <prefix>.best_per_protein.tsv  每个蛋白的最佳命中
  <prefix>.gene_level.tsv        基因级汇总（主统计口径）
  <prefix>.protein_ids.txt       命中蛋白 ID（含 isoform）
  <prefix>.gene_ids.txt          命中基因 ID（剥 isoform）
  <prefix>.hit_proteins.fasta    命中蛋白序列（从 Zunla 蛋白库抽取）
  <prefix>.stats.txt             统计与判读要点

判定口径（重要，写论文时照抄）
  * 初筛保留:  全序列 E-value <= HMM_EVAL（默认 1e-5）→ 标记 strict
  * 另列“仅达包含阈值”: 全序列 E-value > HMM_EVAL 但最佳 domain 的 i-Evalue <= 1e-2
                        → 标记 inclusive_only（Rfam/Pfam 惯例的包含阈值，属最宽松档）
  * 不做任何“为凑数量”的阈值调整；两档数量都会写入 stats。
"""
import argparse
import os
import re
import sys
from collections import defaultdict

DOM_COLS = 23


def parse_args():
    ap = argparse.ArgumentParser(description="解析 hmmsearch domtblout（纯标准库）")
    ap.add_argument("--domtbl", required=True, help="hmmsearch --domtblout 输出")
    ap.add_argument("--pep", required=True, help="Zunla 蛋白 FASTA（抽序列用）")
    ap.add_argument("--prefix", required=True, help="输出文件前缀（含目录）")
    ap.add_argument("--hmm-eval", type=float, default=1e-5, help="初筛全序列 E 值")
    ap.add_argument("--incl-eval", type=float, default=1e-2, help="包含阈值(i-Evalue)")
    return ap.parse_args()


def gene_of(pid):
    """ZLC01G0003110.1 -> ZLC01G0003110"""
    return re.sub(r"\.\d+$", "", pid)


def best_record_of(prots):
    """prots 为 [(pid, record), ...]；取该基因/集合的最佳记录。
    排序优先级: 全序列 E 值最小 -> domain score 最大。"""
    ordered = sorted(prots, key=lambda x: (x[1]["full_e"], -x[1]["dom_score"]))
    return ordered[0][1]


def read_domtbl(path):
    """返回 [(target, tlen, query, qlen, full_e, full_score, dom_ce, dom_ie, dom_score,
              hmm_from, hmm_to, ali_from, ali_to, acc, desc), ...]"""
    rows = []
    bad = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 22:
                bad += 1
                continue
            try:
                rows.append({
                    "target": f[0],
                    "tlen": int(f[2]),
                    "query": f[3],
                    "qlen": int(f[5]),
                    "full_e": float(f[6]),
                    "full_score": float(f[7]),
                    "n_dom": int(f[9]),
                    "dom_ce": float(f[10]),
                    "dom_ie": float(f[11]),
                    "dom_score": float(f[12]),
                    "hmm_from": int(f[15]),
                    "hmm_to": int(f[16]),
                    "ali_from": int(f[17]),
                    "ali_to": int(f[18]),
                    "acc": float(f[20]) if f[20] not in ("-",) else 0.0,
                    "desc": " ".join(f[22:]) if len(f) > 22 else "",
                })
            except (ValueError, IndexError):
                bad += 1
    return rows, bad


def read_fasta(path):
    """返回 OrderedDict id -> seq（保持文件顺序）"""
    seqs = {}
    order = []
    cur = None
    buf = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf)
                cur = line[1:].strip().split()[0]
                order.append(cur)
                buf = []
            else:
                buf.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(buf)
    return seqs, order


def main():
    args = parse_args()
    for p in (args.domtbl, args.pep):
        if not os.path.exists(p):
            print("[错误] 找不到文件: %s" % p, file=sys.stderr)
            return 1

    rows, bad = read_domtbl(args.domtbl)
    if not rows:
        print("[错误] domtblout 里没有可解析的数据行: %s" % args.domtbl, file=sys.stderr)
        return 1

    # ---------- 全部命中明细 ----------
    parsed = os.path.join(os.path.dirname(args.prefix) or ".", os.path.basename(args.prefix) + ".parsed.tsv")
    with open(parsed, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("target\tgene\ttlen\tquery_model\tqlen\tfull_evalue\tfull_score\t"
                 "dom_cevalue\tdom_ievalue\tdom_score\tacc\thmm_from\thmm_to\t"
                 "ali_from\tali_to\thmm_cov\tali_cov\n")
        for r in rows:
            hmm_cov = 100.0 * (r["hmm_to"] - r["hmm_from"] + 1) / r["qlen"] if r["qlen"] else 0.0
            ali_cov = 100.0 * (r["ali_to"] - r["ali_from"] + 1) / r["tlen"] if r["tlen"] else 0.0
            fh.write("%s\t%s\t%d\t%s\t%d\t%.3g\t%.1f\t%.3g\t%.3g\t%.1f\t%.1f\t"
                     "%d\t%d\t%d\t%d\t%.1f\t%.1f\n" % (
                         r["target"], gene_of(r["target"]), r["tlen"], r["query"], r["qlen"],
                         r["full_e"], r["full_score"], r["dom_ce"], r["dom_ie"],
                         r["dom_score"], r["acc"], r["hmm_from"], r["hmm_to"],
                         r["ali_from"], r["ali_to"], hmm_cov, ali_cov))

    # ---------- 每蛋白最佳命中（按 domain score 最大） ----------
    best = {}
    for r in rows:
        k = r["target"]
        if k not in best or r["dom_score"] > best[k]["dom_score"]:
            best[k] = r

    bp = os.path.join(os.path.dirname(args.prefix) or ".", os.path.basename(args.prefix) + ".best_per_protein.tsv")
    with open(bp, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("protein\tgene\ttlen\tquery_model\ttlen_query\tfull_evalue\tfull_score\t"
                 "best_dom_ievalue\tbest_dom_score\thmm_cov\tali_cov\tsignificance\n")
        for pid in sorted(best):
            r = best[pid]
            hmm_cov = 100.0 * (r["hmm_to"] - r["hmm_from"] + 1) / r["qlen"] if r["qlen"] else 0.0
            ali_cov = 100.0 * (r["ali_to"] - r["ali_from"] + 1) / r["tlen"] if r["tlen"] else 0.0
            if r["full_e"] <= args.hmm_eval:
                sig = "strict"
            elif r["dom_ie"] <= args.incl_eval:
                sig = "inclusive_only"
            else:
                sig = "below_inclusive"
            fh.write("%s\t%s\t%d\t%s\t%d\t%.3g\t%.1f\t%.3g\t%.1f\t%.1f\t%.1f\t%s\n" % (
                pid, gene_of(pid), r["tlen"], r["query"], r["qlen"], r["full_e"],
                r["full_score"], r["dom_ie"], r["dom_score"], hmm_cov, ali_cov, sig))

    # ---------- 基因级汇总（主统计口径） ----------
    gene_prots = defaultdict(list)
    for pid, r in best.items():
        gene_prots[gene_of(pid)].append((pid, r))
    n_models = len({r["query"] for r in rows})

    gl = os.path.join(os.path.dirname(args.prefix) or ".", os.path.basename(args.prefix) + ".gene_level.tsv")
    with open(gl, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("gene\tn_proteins\tbest_protein\tn_models_hit\tbest_full_evalue\t"
                 "best_dom_score\tbest_hmm_cov\tbest_ali_cov\tsignificance\n")
        for g in sorted(gene_prots):
            prots = gene_prots[g]
            # 该基因的最佳蛋白 = 全序列 E 最小、其次 domain score 最大
            prots.sort(key=lambda x: (x[1]["full_e"], -x[1]["dom_score"]))
            bpid, br = prots[0]
            models = {p[1]["query"] for p in prots}
            hmm_cov = 100.0 * (br["hmm_to"] - br["hmm_from"] + 1) / br["qlen"] if br["qlen"] else 0.0
            ali_cov = 100.0 * (br["ali_to"] - br["ali_from"] + 1) / br["tlen"] if br["tlen"] else 0.0
            if br["full_e"] <= args.hmm_eval:
                sig = "strict"
            elif br["dom_ie"] <= args.incl_eval:
                sig = "inclusive_only"
            else:
                sig = "below_inclusive"
            fh.write("%s\t%d\t%s\t%d\t%.3g\t%.1f\t%.1f\t%.1f\t%s\n" % (
                g, len(prots), bpid, len(models), br["full_e"], br["dom_score"],
                hmm_cov, ali_cov, sig))

    # ---------- ID 列表 ----------
    strict_prots = sorted(p for p, r in best.items() if r["full_e"] <= args.hmm_eval)
    incl_prots = sorted(p for p, r in best.items() if r["dom_ie"] <= args.incl_eval)
    strict_genes = sorted({gene_of(p) for p in strict_prots})
    incl_genes = sorted({gene_of(p) for p in incl_prots})

    pid_f = os.path.join(os.path.dirname(args.prefix) or ".", os.path.basename(args.prefix) + ".protein_ids.txt")
    gid_f = os.path.join(os.path.dirname(args.prefix) or ".", os.path.basename(args.prefix) + ".gene_ids.txt")
    with open(pid_f, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(strict_prots) + ("\n" if strict_prots else ""))
    with open(gid_f, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(strict_genes) + ("\n" if strict_genes else ""))

    pid_i = os.path.join(os.path.dirname(args.prefix) or ".", os.path.basename(args.prefix) + ".protein_ids_inclusive.txt")
    gid_i = os.path.join(os.path.dirname(args.prefix) or ".", os.path.basename(args.prefix) + ".gene_ids_inclusive.txt")
    with open(pid_i, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(incl_prots) + ("\n" if incl_prots else ""))
    with open(gid_i, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(incl_genes) + ("\n" if incl_genes else ""))

    # ---------- 抽取命中蛋白序列 ----------
    seqs, order = read_fasta(args.pep)
    fa = os.path.join(os.path.dirname(args.prefix) or ".", os.path.basename(args.prefix) + ".hit_proteins.fasta")
    n_out = 0
    with open(fa, "w", encoding="utf-8", newline="\n") as fh:
        for pid in order:
            if pid in best:
                fh.write(">%s\n" % pid)
                s = seqs[pid]
                for i in range(0, len(s), 60):
                    fh.write(s[i:i + 60] + "\n")
                n_out += 1

    # ---------- 统计 ----------
    n_below = sum(1 for r in best.values() if r["dom_ie"] > args.incl_eval)
    stats = os.path.join(os.path.dirname(args.prefix) or ".", os.path.basename(args.prefix) + ".stats.txt")
    with open(stats, "w", encoding="utf-8", newline="\n") as fh:
        A = lambda s: fh.write(s + "\n")
        A("hmmsearch 结果解析统计")
        A("=" * 70)
        A("domtblout        : %s" % args.domtbl)
        A("Zunla 蛋白库     : %s" % args.pep)
        A("判定口径         : 初筛 全序列 E<=%g ; 包含阈值 i-Evalue<=%g"
          % (args.hmm_eval, args.incl_eval))
        A("")
        A("【原始规模】")
        A("  domain 命中行数        : %d" % len(rows))
        if bad:
            A("  跳过的不可解析行       : %d" % bad)
        A("  命中的蛋白数           : %d" % len(best))
        A("  命中的基因数           : %d" % len(gene_prots))
        A("  被命中的 HMM 模型数    : %d" % n_models)
        A("")
        A("【按判定口径】")
        A("  strict  (E<=%g)                  : %d 蛋白 / %d 基因 <- 主口径"
          % (args.hmm_eval, len(strict_prots), len(strict_genes)))
        A("  inclusive (i-Evalue<=%g)         : %d 蛋白 / %d 基因"
          % (args.incl_eval, len(incl_prots), len(incl_genes)))
        A("  仅达 inclusive、未达 strict       : %d 蛋白"
          % (len(set(incl_prots) - set(strict_prots))))
        A("  连 inclusive 都不到               : %d 蛋白（已单独保留在 parsed.tsv）" % n_below)
        A("")
        A("【分布特征（基因级，最佳命中）】")
        if gene_prots:
            best_recs = [best_record_of(prots) for prots in gene_prots.values()]
            fes = sorted(r["full_e"] for r in best_recs)
            scs = sorted(r["dom_score"] for r in best_recs)
            A("  全序列 E 值 最小/中位/最大 : %.3g / %.3g / %.3g"
              % (fes[0], fes[len(fes) // 2], fes[-1]))
            A("  domain score 最小/中位/最大: %.1f / %.1f / %.1f"
              % (scs[0], scs[len(scs) // 2], scs[-1]))
        A("")
        A("【说明】")
        A("  * 未为了凑某个预设数量调整阈值：以上两档数量均按固定口径统计。")
        A("  * hmmsearch 的 -E 参数只影响「输出哪些命中」，不改变打分；")
        A("    真正的判定依据是这里的 E-value / i-Evalue，均取自原始 domtblout。")
        A("  * 命中蛋白序列已抽取到 %s（%d 条）"
          % (os.path.basename(fa), n_out))

    print("  解析完成:")
    print("    strict    : %d 蛋白 / %d 基因" % (len(strict_prots), len(strict_genes)))
    print("    inclusive : %d 蛋白 / %d 基因" % (len(incl_prots), len(incl_genes)))
    print("    %s" % parsed)
    print("    %s" % gl)
    print("    %s" % stats)
    return 0


if __name__ == "__main__":
    sys.exit(main())
