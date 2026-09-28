#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
final_summary.py —— 阶段 07：汇总最终 CaTRM 列表，附每个成员的证据来源与筛选理由

只用标准库。

输入
  --merged       05_merge/CaTRM_candidates.merged.tsv
  --per-protein  06_domain/CaTRM_vs_Pfam.per_protein.tsv（可缺，缺则结构域列填 NA）
  --focus-matrix 06_domain/CaTRM_vs_Pfam.focus_matrix.tsv（可缺）
  --gene-level   04_hmmsearch/TRM.hmmsearch.gene_level.tsv（可缺）
  --gff3         Zunla GFF3（可缺；有则补 gene 坐标与染色体）
  --pep          Zunla 蛋白 FASTA
  --prefix       输出前缀

分层（conservative / standard / review，可同时出现在 standard 与 review）
  conservative  双方法支持 + 含 PF14309 核心域
  standard      双方法支持，或缺核心域但双方法支持，或单方法+HMM严格命中
  review        仅达 HMM 包含阈值 / 仅单方法且无结构域支持 / 覆盖度偏低等

产出
  <prefix>.gene_list.tsv     最终基因列表（主交付）
  <prefix>.gene_ids.txt      仅基因 ID（34 行格式）
  <prefix>.protein_list.tsv  蛋白级列表
  <prefix>.proteins.fasta    最终代表蛋白序列（每基因取最佳蛋白）
  <prefix>.evidence.tsv      每成员证据来源与筛选理由
  <prefix>.review.tsv        待人工确认清单（不删除，只列出）
  <prefix>.stats.txt         分层计数与说明
"""
import argparse
import os
import re
import sys
from collections import defaultdict


def parse_args():
    ap = argparse.ArgumentParser(description="汇总最终 CaTRM 列表（纯标准库）")
    ap.add_argument("--merged", required=True)
    ap.add_argument("--per-protein", default="")
    ap.add_argument("--focus-matrix", default="")
    ap.add_argument("--gene-level", default="")
    ap.add_argument("--gff3", default="")
    ap.add_argument("--pep", default="")
    ap.add_argument("--prefix", required=True)
    return ap.parse_args()


def gene_of(x):
    return re.sub(r"\.\d+$", "", x.strip())


def read_tsv(path):
    """返回 (表头列表, [dict,...])"""
    if not path or not os.path.exists(path):
        return [], []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        rows = []
        header = None
        for i, line in enumerate(fh):
            f = line.rstrip("\n").split("\t")
            if i == 0:
                header = f
                continue
            if len(f) == len(header):
                rows.append(dict(zip(header, f)))
        return header or [], rows


def read_text_lines(path):
    if not path or not os.path.exists(path):
        return []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        return [l.strip() for l in fh if l.strip() and not l.startswith("#")]


def read_fasta(path):
    seqs, order = {}, []
    if not path or not os.path.exists(path):
        return seqs, order
    cur, buf = None, []
    opener = open
    if path.endswith(".gz"):
        import gzip
        opener = gzip.open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
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


def read_gff3_coords(path):
    """gene ID -> (chr, start, end, strand, n_mrna)"""
    out = {}
    if not path or not os.path.exists(path):
        return out
    import gzip
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            m = re.search(r"(?:^|;)ID=([^;]+)", f[8])
            if not m:
                continue
            gid = m.group(1)
            out[gid] = (f[0], int(f[3]), int(f[4]), f[6], 0)
    # 补 mRNA 计数
    with op(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "mRNA":
                continue
            m = re.search(r"(?:^|;)Parent=([^;]+)", f[8])
            if not m:
                continue
            p = m.group(1)
            if p in out:
                c = out[p]
                out[p] = (c[0], c[1], c[2], c[3], c[4] + 1)
    return out


def main():
    args = parse_args()
    odir = os.path.dirname(args.prefix) or "."
    if odir and not os.path.isdir(odir):
        os.makedirs(odir)

    _, merged = read_tsv(args.merged)
    if not merged:
        print("[错误] 合并表为空或不存在: %s" % args.merged, file=sys.stderr)
        return 1

    _, per_prot = read_tsv(args.per_protein)
    prot_dom = {r["protein"]: r for r in per_prot if "protein" in r}

    _, gene_level = read_tsv(args.gene_level)
    hmm_gene = {r["gene"]: r for r in gene_level if "gene" in r}

    coords = read_gff3_coords(args.gff3)
    seqs, order = read_fasta(args.pep)

    # ---------- 分层判定 ----------
    have_domain_data = bool(prot_dom)   # 没有结构域数据时不能声称"含核心域"

    def tier_of(r):
        ev = r.get("evidence", "")
        hbp = r.get("hmm_best_protein", "NA")
        has_core = "no"
        if hbp != "NA" and hbp in prot_dom:
            has_core = prot_dom[hbp].get("has_TRM_core_domain", "no")
        if not have_domain_data:
            has_core = "NA"

        reasons = []
        if ev == "BLAST+HMM" and has_core == "yes":
            tier = "conservative"
            reasons.append("BLAST 与 HMM 双重支持")
            reasons.append("含 PF14309 核心域")
        elif ev == "BLAST+HMM":
            tier = "standard"
            reasons.append("BLAST 与 HMM 双重支持")
            if has_core == "NA":
                reasons.append("未做结构域核查（缺少阶段 06 结果），核心域未经验证")
            else:
                reasons.append("未检出 PF14309（可能因该域短/数据库版本差异）")
        elif ev in ("HMM-only",):
            if has_core == "yes":
                tier = "standard"
                reasons.append("HMM 严格命中且含 PF14309 核心域（BLAST 漏检，profile HMM 更敏感）")
            else:
                tier = "standard"
                reasons.append("HMM 严格命中（BLAST 漏检）")
                if has_core == "NA":
                    reasons.append("未做结构域核查，核心域未经验证")
                else:
                    reasons.append("未检出 PF14309 -> 建议人工确认")
        elif ev == "BLAST-only":
            tier = "review"
            reasons.append("仅 BLAST 支持，HMM 未显著命中")
            reasons.append("需人工确认是否为远缘/片段同源")
        elif ev == "HMM-inclusive-only":
            tier = "review"
            reasons.append("仅达 HMM 包含阈值(i-Evalue<=1e-2)，未达严格阈值")
        else:
            tier = "review"
            reasons.append("证据不足")

        mr = r.get("manual_review", "-")
        if mr and mr != "-":
            reasons.append("自动标记: %s" % mr)
        return tier, has_core, reasons

    # ---------- 写基因列表 ----------
    gl = os.path.join(odir, os.path.basename(args.prefix) + ".gene_list.tsv")
    evf = os.path.join(odir, os.path.basename(args.prefix) + ".evidence.tsv")
    rf = os.path.join(odir, os.path.basename(args.prefix) + ".review.tsv")
    tier_count = defaultdict(int)
    rows_for_fasta = []

    with open(gl, "w", encoding="utf-8", newline="\n") as fh, \
         open(evf, "w", encoding="utf-8", newline="\n") as ef, \
         open(rf, "w", encoding="utf-8", newline="\n") as rfh:
        fh.write("gene\tchrom\tstart\tend\tstrand\tn_transcripts\ttier\tevidence\t"
                 "has_PF14309\thmm_significance\thmm_best_protein\thmm_best_full_evalue\t"
                 "blast_best_pident\tblast_best_qcov\tn_at_queries\tat_genes\tat_symbols\n")
        ef.write("gene\ttier\tevidence\tselection_reasons\n")
        rfh.write("gene\ttier\tevidence\treview_reasons\n")

        for r in sorted(merged, key=lambda x: x.get("gene", "")):
            g = r.get("gene", "")
            tier, has_core, reasons = tier_of(r)
            tier_count[tier] += 1

            c = coords.get(g)
            if c is None:
                chrom = start = end = strand = nmrna = "NA"
            else:
                chrom, start, end, strand, nmrna = c
            has_core_out = has_core if has_core != "NA" else "NA"
            fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                g, chrom, start, end, strand, nmrna,
                tier, r.get("evidence", "NA"), has_core_out,
                r.get("hmm_significance", "NA"), r.get("hmm_best_protein", "NA"),
                r.get("hmm_best_full_evalue", "NA"),
                r.get("blast_best_pident", "NA"), r.get("blast_best_qcov", "NA"),
                r.get("n_at_queries", "NA"), r.get("at_genes", "NA"),
                r.get("at_symbols", "NA")))
            ef.write("%s\t%s\t%s\t%s\n" % (g, tier, r.get("evidence", "NA"),
                                           "; ".join(reasons)))
            if tier == "review":
                rfh.write("%s\t%s\t%s\t%s\n" % (g, tier, r.get("evidence", "NA"),
                                                "; ".join(reasons)))
            if r.get("hmm_best_protein", "NA") != "NA":
                rows_for_fasta.append(r["hmm_best_protein"])
            for p in (r.get("blast_proteins", "") or "").split(","):
                if p and p != "NA":
                    rows_for_fasta.append(p)

    # ---------- 基因 ID 清单 ----------
    gid_f = os.path.join(odir, os.path.basename(args.prefix) + ".gene_ids.txt")
    with open(gid_f, "w", encoding="utf-8", newline="\n") as fh:
        for r in sorted(merged, key=lambda x: x.get("gene", "")):
            fh.write(r.get("gene", "") + "\n")

    # ---------- 蛋白列表 ----------
    pl = os.path.join(odir, os.path.basename(args.prefix) + ".protein_list.tsv")
    prot_rows = []
    for r in merged:
        g = r.get("gene", "")
        bp = r.get("hmm_best_protein", "NA")
        if bp != "NA":
            prot_rows.append((g, bp, "HMM_best", r.get("evidence", "")))
        for p in (r.get("blast_proteins", "") or "").split(","):
            if p and p != "NA":
                prot_rows.append((g, p, "BLAST_hit", r.get("evidence", "")))
    seen = set()
    with open(pl, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("gene\tprotein\tsource\tevidence\n")
        for g, p, src, ev in prot_rows:
            if (g, p) in seen:
                continue
            seen.add((g, p))
            fh.write("%s\t%s\t%s\t%s\n" % (g, p, src, ev))

    # ---------- 代表蛋白 FASTA（每基因取 HMM 最佳蛋白） ----------
    fa = os.path.join(odir, os.path.basename(args.prefix) + ".proteins.fasta")
    n_fa = 0
    with open(fa, "w", encoding="utf-8", newline="\n") as fh:
        for r in sorted(merged, key=lambda x: x.get("gene", "")):
            bp = r.get("hmm_best_protein", "NA")
            if bp == "NA" or bp not in seqs:
                # 退而取该基因第一个 BLAST 蛋白
                cands = [p for p in (r.get("blast_proteins", "") or "").split(",")
                         if p and p != "NA"]
                bp = cands[0] if cands else None
            if not bp or bp not in seqs:
                continue
            s = seqs[bp]
            fh.write(">%s | gene=%s | evidence=%s\n" % (bp, r.get("gene", ""),
                                                          r.get("evidence", "")))
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
            n_fa += 1

    # ---------- 统计 ----------
    stats = os.path.join(odir, os.path.basename(args.prefix) + ".stats.txt")
    with open(stats, "w", encoding="utf-8", newline="\n") as fh:
        A = lambda s: fh.write(s + "\n")
        A("阶段 07：最终 CaTRM 候选汇总")
        A("=" * 70)
        A("主统计单位: gene")
        A("")
        A("【分层结果】")
        A("  conservative (双方法 + PF14309)  : %d" % tier_count.get("conservative", 0))
        A("  standard     (双方法 或 单方法+HMM严格) : %d" % tier_count.get("standard", 0))
        A("  review       (待人工确认)        : %d" % tier_count.get("review", 0))
        A("  ------------------------------------------")
        A("  ★ 候选基因总数                   : %d" % len(merged))
        A("")
        A("【证据来源分布】")
        ev_c = defaultdict(int)
        for r in merged:
            ev_c[r.get("evidence", "NA")] += 1
        for k in ("BLAST+HMM", "BLAST-only", "HMM-only", "HMM-inclusive-only"):
            A("  %-22s %d" % (k, ev_c.get(k, 0)))
        A("")
        A("【结构域】")
        n_core = sum(1 for r in merged
                     if r.get("hmm_best_protein", "NA") in prot_dom
                     and prot_dom[r["hmm_best_protein"]].get("has_TRM_core_domain") == "yes")
        A("  检出 PF14309(TRM 核心域) 的基因: %d" % n_core)
        A("  未检出 PF14309 的基因        : %d" % (len(merged) - n_core))
        A("  说明: 未检出不等于不是家族成员——PF14309 是较短的 C 端域，")
        A("        在远缘同源或部分序列上可能低于 Pfam 的显著阈值。")
        A("")
        A("【未删除任何候选】")
        A("  review 层 %d 个候选已单列到 review.tsv，仅标记不删除。" % tier_count.get("review", 0))
        A("  论文中如只采用 conservative + standard，请明确说明 review 层的处置理由。")
        A("")
        A("【产出】")
        A("  gene_list.tsv     最终基因列表（主交付，含坐标/证据/分层）")
        A("  gene_ids.txt      仅基因 ID")
        A("  protein_list.tsv  蛋白级列表")
        A("  proteins.fasta    每基因代表蛋白序列（%d 条）" % n_fa)
        A("  evidence.tsv      每成员证据来源与筛选理由")
        A("  review.tsv        待人工确认清单")

    print("  ★ 最终候选基因 : %d" % len(merged))
    print("      conservative : %d" % tier_count.get("conservative", 0))
    print("      standard     : %d" % tier_count.get("standard", 0))
    print("      review       : %d" % tier_count.get("review", 0))
    print("    代表蛋白 FASTA : %d 条" % n_fa)
    print("    %s" % gl)
    return 0


if __name__ == "__main__":
    sys.exit(main())
