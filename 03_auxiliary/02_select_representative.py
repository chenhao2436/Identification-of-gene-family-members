# -*- coding: utf-8 -*-
"""
【脚本 2】每个基因挑选一条代表蛋白（最长转录本）

策略（不写死 .1）：
    1 个 Gene -> N 个 protein isoform -> 逐一比较氨基酸长度 -> 取最长者作为代表蛋白
    并列最长时：优先“序列完全相同”者取转录本编号最小；
                若序列不同但长度相同，取编号最小，并在报告中列出全部并列项。

输入：02_script_output/TRM_proteins.fasta   （脚本 1 的全转录本结果）
输出：02_script_output/TRM_representative_proteins.fasta   代表蛋白序列（34 条）
      02_script_output/TRM_representative_proteins.tsv     代表蛋白汇总表
      02_script_output/TRM_representative_report.txt       挑选过程报告

注意：输出的中间名文件请按项目命名规范改名归位为
      01_At_TRM_query_34.fasta / 02_At_TRM_annotation.tsv /
      select_representative_report.txt。
"""
import os
import re
from collections import OrderedDict

WORK = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.environ.get(
    "TRM_OUTPUT_DIR", os.path.join(WORK, "..", "02_script_output"))
if not os.path.isdir(OUTDIR):
    os.makedirs(OUTDIR)
IN_FA = os.path.join(OUTDIR, "TRM_proteins.fasta")
OUT_FA = os.path.join(OUTDIR, "TRM_representative_proteins.fasta")
OUT_TSV = os.path.join(OUTDIR, "TRM_representative_proteins.tsv")
OUT_RPT = os.path.join(OUTDIR, "TRM_representative_report.txt")

HEAD_RE = re.compile(r"^>(\S+)\s*\|\s*gene=(\S+)\s*\|\s*Symbols:([^|]*)\|\s*(.*)$")


def read_fasta(path):
    """读入脚本1产出的 FASTA -> [(tid, gene, symbols, desc, seq), ...]"""
    recs = []
    cur = None
    buf = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            if line.startswith(">"):
                if cur:
                    recs.append((cur[0], cur[1], cur[2], cur[3], "".join(buf)))
                m = HEAD_RE.match(line)
                if not m:
                    raise ValueError("无法解析的 FASTA 头: %s" % line)
                cur = (m.group(1), m.group(2), m.group(3).strip(), m.group(4).strip())
                buf = []
            else:
                buf.append(line.strip())
    if cur:
        recs.append((cur[0], cur[1], cur[2], cur[3], "".join(buf)))
    return recs


def tid_key(tid):
    """AT3G02170.2 -> (3, 2170, 2)，用于稳定排序"""
    base, _, ver = tid.partition(".")
    m = re.match(r"AT(\d)G(\d+)", base.upper())
    chrom = int(m.group(1)) if m else 0
    num = int(m.group(2)) if m else 0
    return (chrom, num, int(ver) if ver.isdigit() else 0)


def main():
    recs = read_fasta(IN_FA)

    # 按基因分组（保持首次出现顺序）
    genes = OrderedDict()
    for tid, gene, sym, desc, seq in recs:
        genes.setdefault(gene, []).append(
            {"tid": tid, "symbols": sym, "desc": desc, "seq": seq, "len": len(seq)}
        )

    rep_rows = []
    report = []
    report.append("每个基因代表蛋白（最长转录本）挑选报告")
    report.append("=" * 78)
    report.append("规则：1 Gene -> N isoform -> 比较氨基酸长度 -> 取最长者")
    report.append("并列最长：序列相同者取转录本编号最小；序列不同则在下方标注 [并列-序列不同]")
    report.append("")

    n_tie_same = 0
    n_tie_diff = 0

    for gene, iso in genes.items():
        iso.sort(key=lambda r: tid_key(r["tid"]))
        max_len = max(r["len"] for r in iso)
        cands = [r for r in iso if r["len"] == max_len]

        # 并列时优先序列相同（去重后只剩 1 条独特序列即视为完全一致）
        uniq_seqs = {r["seq"] for r in cands}
        chosen = cands[0]
        tie_note = ""

        if len(cands) > 1:
            if len(uniq_seqs) == 1:
                n_tie_same += 1
                tie_note = "并列最长 %d 条且序列完全相同，取编号最小 %s" % (
                    len(cands), chosen["tid"])
            else:
                n_tie_diff += 1
                tie_note = "[并列-序列不同] 并列最长 %d 条: %s，取编号最小 %s" % (
                    len(cands), ", ".join(r["tid"] for r in cands), chosen["tid"])

        rep_rows.append({
            "gene": gene,
            "tid": chosen["tid"],
            "symbols": chosen["symbols"],
            "desc": chosen["desc"],
            "len": chosen["len"],
            "n_iso": len(iso),
            "max_len": max_len,
            "tie_note": tie_note,
        })

        report.append("Gene %-10s isoform=%d  长度: %s" % (
            gene, len(iso), ", ".join("%s(%d)" % (r["tid"], r["len"]) for r in iso)))
        report.append("    代表蛋白 -> %s  (%d aa)%s" % (
            chosen["tid"], chosen["len"], ("   " + tie_note) if tie_note else ""))

    # ---- 写出 FASTA ----
    with open(OUT_FA, "w", encoding="utf-8", newline="\n") as fa:
        for r in rep_rows:
            fa.write(">%s | gene=%s | Symbols:%s | %s\n" % (
                r["tid"], r["gene"], r["symbols"], r["desc"]))
            seq = next(x["seq"] for x in genes[r["gene"]] if x["tid"] == r["tid"])
            for i in range(0, len(seq), 60):
                fa.write(seq[i:i + 60] + "\n")

    # ---- 写出 TSV ----
    with open(OUT_TSV, "w", encoding="utf-8", newline="\n") as tsv:
        # 第 3 列 gene_id 为 transcript 去掉版本号后的形式，专供下游
        # BLAST 结果映射（如 trm_blast_pipeline.sh 把 query 归回 At 基因号）
        tsv.write("gene\trepresentative_transcript\tgene_id\tsymbols\tdescription\t"
                  "length_aa\tisoform_count\tselected_is_longest\tnote\n")
        for r in rep_rows:
            tsv.write("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\n" % (
                r["gene"], r["tid"], r["tid"].split(".")[0], r["symbols"],
                r["desc"], r["len"], r["n_iso"], "yes", r["tie_note"]))

    # ---- 报告尾部统计 ----
    report.append("")
    report.append("=" * 78)
    report.append("基因总数: %d" % len(rep_rows))
    report.append("代表蛋白总数: %d" % len(rep_rows))
    n_changed = sum(1 for r in rep_rows if not r["tid"].endswith(".1"))
    report.append("代表蛋白不是 .1 的基因数: %d" % n_changed)
    if n_changed:
        for r in rep_rows:
            if not r["tid"].endswith(".1"):
                report.append("    %s -> %s (%d aa)" % (r["gene"], r["tid"], r["len"]))
    report.append("并列最长且序列完全相同: %d 个基因" % n_tie_same)
    report.append("并列最长但序列不同: %d 个基因" % n_tie_diff)
    report.append("")
    report.append("氨基酸长度范围: %d - %d aa；平均 %.1f aa" % (
        min(r["len"] for r in rep_rows), max(r["len"] for r in rep_rows),
        sum(r["len"] for r in rep_rows) / len(rep_rows)))

    with open(OUT_RPT, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(report) + "\n")

    print("基因数: %d -> 代表蛋白数: %d" % (len(genes), len(rep_rows)))
    print("代表蛋白非 .1 的基因: %d 个" % n_changed)
    print("并列(序列相同): %d 个；并列(序列不同): %d 个" % (n_tie_same, n_tie_diff))
    print("输出: %s" % OUT_FA)


if __name__ == "__main__":
    main()
