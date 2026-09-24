# -*- coding: utf-8 -*-
"""
从 Araport11 蛋白组 FASTA 中，按基因 ID 提取对应蛋白序列。
- 基因 ID 不区分大小写（At3g02170 / AT3G02170 均可）
- 一个基因有多个转录本时，默认全部保留（可用 --isoform N 只取第 N 个）
输出（写入 02_script_output/，再经 02_select_representative.py 命名归位）：
  TRM_proteins.fasta      —— 提取到的序列（FASTA）
  TRM_proteins.tsv        —— 汇总表（基因ID/转录本/符号/描述/长度）
  TRM_not_found.txt       —— 未找到的基因ID

注意：本脚本会把中间结果写到 02_script_output/。若该目录已存在正式的
      At_TRM_all_isoforms_87.fasta 等成品，重跑前请先备份，避免覆盖。
"""
import gzip
import os
import re
import sys

WORK = os.path.dirname(os.path.abspath(__file__))
# 输入目录默认指向项目内的 01_At_TRM_query/，可用环境变量覆盖：
#   Windows:  set TRM_INPUT_DIR=D:\path\to\dir
#   Linux  :  export TRM_INPUT_DIR=/path/to/dir
INPUT_DIR = os.environ.get(
    "TRM_INPUT_DIR", os.path.join(WORK, "..", "01_At_TRM_query"))
FASTA = os.path.join(INPUT_DIR, "Araport11_pep_20250411.gz")
# 输出目录默认 02_script_output/，同样可用 TRM_OUTPUT_DIR 覆盖
OUTDIR = os.environ.get(
    "TRM_OUTPUT_DIR", os.path.join(WORK, "..", "02_script_output"))
if not os.path.isdir(OUTDIR):
    os.makedirs(OUTDIR)

GENE_IDS = [
    "At3g02170", "At5g15580", "At1g18620", "At1g74160", "At3g63430",
    "At3g05750", "At3g58650", "At5g26910", "At4g00770", "At4g11780",
    "At4g23020", "At2g36420", "At5g03670", "At2g45900", "At3g61380",
    "At4g00440", "At1g01695", "At2g20240", "At3g53540", "At4g28760",
    "At5g43880", "At2g17550", "At4g25430", "At5g51850", "At5g62170",
    "At5g01370", "At5g58630", "At5g02390", "At1g07620", "At1g63670",
    "At2g39435", "At1g67040", "At5g42710", "At3g24630",
]

header_re = re.compile(
    r"^>(\S+)"                      # 1: transcript id, e.g. AT1G01010.1
    r"(?:\s*\|\s*Symbols:([^|]*))?"  # 2: symbols
    r"(?:\s*\|\s*([^|]*))?"          # 3: description
    r"(?:\s*\|\s*(.*))?$"            # 4: location / LENGTH
)


def main():
    if not os.path.exists(FASTA):
        sys.exit("找不到输入文件: %s" % FASTA)

    wanted = {g.lower(): g for g in GENE_IDS}
    hits = {}          # gene_lower -> list of record dicts
    order = []         # 记录首次出现的顺序
    cur = None

    with gzip.open(FASTA, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            if line[0] == ">":
                if cur is not None:
                    hits.setdefault(cur["gene_lower"], []).append(cur)
                    if cur["gene_lower"] not in order:
                        order.append(cur["gene_lower"])
                    cur = None
                m = header_re.match(line)
                if not m:
                    continue
                tid = m.group(1)
                gene_lower = tid.split(".")[0].lower()
                if gene_lower not in wanted:
                    continue
                length = None
                loc = m.group(4) or ""
                lm = re.search(r"LENGTH=(\d+)", loc)
                if lm:
                    length = int(lm.group(1))
                cur = {
                    "gene_lower": gene_lower,
                    "gene": wanted[gene_lower],
                    "tid": tid,
                    "symbols": (m.group(2) or "").strip(),
                    "desc": (m.group(3) or "").strip(),
                    "loc": loc.strip(),
                    "declared_len": length,
                    "seq": [],
                }
            elif cur is not None:
                cur["seq"].append(line.strip())

    if cur is not None:
        hits.setdefault(cur["gene_lower"], []).append(cur)
        if cur["gene_lower"] not in order:
            order.append(cur["gene_lower"])

    # ---- 输出 FASTA（按用户给定基因ID顺序，多转录本按编号排序）----
    out_fa = os.path.join(OUTDIR, "TRM_proteins.fasta")
    out_tsv = os.path.join(OUTDIR, "TRM_proteins.tsv")
    out_nf = os.path.join(OUTDIR, "TRM_not_found.txt")

    found_genes = []
    n_rec = 0
    with open(out_fa, "w", encoding="utf-8", newline="\n") as fa, \
         open(out_tsv, "w", encoding="utf-8", newline="\n") as tsv:
        tsv.write("query_gene\ttranscript\tsymbols\tdescription\tlength_aa\tdeclared_length\tlocation\n")
        for g in GENE_IDS:
            recs = hits.get(g.lower())
            if not recs:
                continue
            found_genes.append(g)
            recs.sort(key=lambda r: r["tid"])
            for r in recs:
                seq = "".join(r["seq"])
                fa.write(">%s | gene=%s | Symbols:%s | %s\n" % (
                    r["tid"], r["gene"], r["symbols"], r["desc"]))
                for i in range(0, len(seq), 60):
                    fa.write(seq[i:i + 60] + "\n")
                tsv.write("%s\t%s\t%s\t%s\t%d\t%s\t%s\n" % (
                    r["gene"], r["tid"], r["symbols"], r["desc"], len(seq),
                    r["declared_len"] if r["declared_len"] else "", r["loc"]))
                n_rec += 1

    missing = [g for g in GENE_IDS if g.lower() not in hits]
    with open(out_nf, "w", encoding="utf-8", newline="\n") as f:
        for g in missing:
            f.write(g + "\n")

    # ---- 控制台报告 ----
    print("输入基因数: %d" % len(GENE_IDS))
    print("匹配到基因数: %d" % len(found_genes))
    print("提取到转录本数: %d" % n_rec)
    multi = {g: len(hits[g.lower()]) for g in found_genes}
    multi = {k: v for k, v in multi.items() if v > 1}
    if multi:
        print("\n多转录本基因 (共 %d 个):" % len(multi))
        for k, v in multi.items():
            print("  %s -> %d 个转录本" % (k, v))
    if missing:
        print("\n未找到 (%d 个):" % len(missing))
        for g in missing:
            print("  " + g)
    else:
        print("\n全部 %d 个基因 ID 均已匹配。" % len(GENE_IDS))


if __name__ == "__main__":
    main()
