#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_hmmer_parsers.py —— 用合成 HMMER 输出测试三个解析脚本

为什么需要: 服务器上无法先跑一遍验证，本地又没有 HMMER。
因此用严格按 HMMER 3.x --domtblout 格式构造的合成数据，
把三个解析脚本与合并/汇总脚本跑一遍并断言结果。

覆盖的关键场景
  * 同一蛋白的多个 domain 命中 -> 取最佳
  * 同一基因的多个 isoform -> 基因级去重
  * strict vs inclusive_only vs below_inclusive 三档判定
  * 关注域 PF14309 / PF14383 的矩阵与 has_TRM_core_domain 标志
  * BLAST 与 HMM 的四种证据标签（BLAST+HMM / BLAST-only / HMM-only / HMM-inclusive-only）
  * 无 BLAST 文件时的降级（BLAST 列填 NA）

用法: python test_hmmer_parsers.py
退出码 0 = 全部通过
"""
import os
import shutil
import subprocess
import sys

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
TMP = os.path.join(os.environ.get("TEMP", "."), "trm_hmmer_test")

# ---------------------------------------------------------------------------
# 合成 Zunla 蛋白库
# ---------------------------------------------------------------------------
PROTEINS = {
    "ZLC01G0003110.1": "MSAKLLYNLSDENPNLNKQFGCMNGIFQVFYRQHCPATPVTVSGGAEKSLPPGERRGSVG",
    "ZLC01G0003110.2": "MSAKLLYNLSDENPNLNKQFGCMNGIFQVFYRQHCPATPVTVSGGAEKSLPPGERRG",
    "ZLC02G0005000.1": "MEDQVGFGFRPNDEELVGHYLRNKIEGNTSRDVEVAISEVNICSYDPWNLRFQSKYKS",
    "ZLC03G0007000.1": "MTRMTRMTRMTRMTRMTRMTRMTRMTRMTRMTRMTRMTRMTRMTRMTRMTRMTRMTR",
    "ZLC04G0009000.1": "MAAAAGGGGKLPPSSTTDDEEKKRRQQNNVVYYWWFFAALLLIIIVVVCCCGGGSSSTT",
}


def write_pep(path):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for pid, seq in PROTEINS.items():
            fh.write(">%s\n" % pid)
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


# ---------------------------------------------------------------------------
# 合成 hmmsearch --domtblout（23 列）
#   target=蛋白  query=HMM模型
# ---------------------------------------------------------------------------
def dom_line(target, tlen, query, qlen, full_e, full_score, ndom,
             ce, ie, dscore, hf, ht, af, at, acc, desc):
    return ("%-24s %-10s %6d  %-18s %-10s %6d  %9.2e %8.1f %8.1f %6d  "
            "%9.2e %9.2e %8.1f %8.1f %5d %5d %5d %5d %5d %5d %6.2f %s" % (
                target, "-", tlen, query, "-", qlen, full_e, full_score, 0.0, ndom,
                ce, ie, dscore, 0.0, hf, ht, af, at, af, at, acc, desc))


def write_hmmsearch_domtbl(path):
    lines = [
        "#                                                               --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord",
        "# target        name  tlen  query  acc  qlen  E-value  score  bias  #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target",
    ]
    # 蛋白1: 同一蛋白两个 domain 命中，取高分那个
    lines.append(dom_line("ZLC01G0003110.1", 60, "TRM.hmm", 200, 1e-40, 150.0, 2,
                          1e-40, 1e-42, 148.0, 10, 60, 5, 55, 0.95, "TRM domain"))
    lines.append(dom_line("ZLC01G0003110.1", 60, "TRM.hmm", 200, 1e-40, 150.0, 2,
                          1e-30, 1e-32, 120.0, 80, 150, 10, 58, 0.90, "TRM domain"))
    # 蛋白2: strict 命中
    lines.append(dom_line("ZLC01G0003110.2", 55, "TRM.hmm", 200, 1e-20, 90.0, 1,
                          1e-20, 1e-22, 88.0, 10, 50, 3, 45, 0.88, "TRM domain"))
    # 蛋白3: 仅达 inclusive (full E 不够，i-Evalue 够)
    lines.append(dom_line("ZLC02G0005000.1", 60, "TRM.hmm", 200, 3.5e-4, 25.0, 1,
                          1e-4, 5e-3, 24.0, 20, 60, 8, 48, 0.70, "TRM domain"))
    # 蛋白4: 连 inclusive 都不到
    lines.append(dom_line("ZLC03G0007000.1", 58, "TRM.hmm", 200, 2.0e-1, 8.0, 1,
                          5e-2, 3e-2, 7.0, 30, 55, 12, 40, 0.40, "weak"))
    # 蛋白5: strict 命中（独立基因）
    lines.append(dom_line("ZLC04G0009000.1", 60, "TRM.hmm", 200, 1e-15, 70.0, 1,
                          1e-15, 1e-16, 68.0, 15, 58, 4, 50, 0.85, "TRM domain"))
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# 合成 hmmscan --domtblout（列序相同，但 target=Pfam 模型、query=蛋白）
# ---------------------------------------------------------------------------
def scan_line(model, mlen, protein, plen, full_e, full_score, ndom,
              ce, ie, dscore, hf, ht, af, at, acc, desc):
    return ("%-14s %-10s %6d  %-18s %-10s %6d  %9.2e %8.1f %8.1f %6d  "
            "%9.2e %9.2e %8.1f %8.1f %5d %5d %5d %5d %5d %5d %6.2f %s" % (
                model, "-", mlen, protein, "-", plen, full_e, full_score, 0.0, ndom,
                ce, ie, dscore, 0.0, hf, ht, af, at, af, at, acc, desc))


def write_hmmscan_domtbl(path):
    lines = [
        "# target  acc  tlen  query  acc  qlen  E-value score bias # of c-Evalue i-Evalue score bias hmm from to ali from to env from to acc desc",
    ]
    # 蛋白1、蛋白2、蛋白5 含 PF14309（TRM 核心域）
    lines.append(scan_line("PF14309.20", 220, "ZLC01G0003110.1", 60, 1e-30, 100.0, 1,
                           1e-30, 1e-32, 99.0, 30, 200, 18, 58, 0.92,
                           "LONGIFOLIA 1/2-like C-terminal domain"))
    lines.append(scan_line("PF14309.20", 220, "ZLC01G0003110.2", 55, 1e-18, 65.0, 1,
                           1e-18, 1e-20, 64.0, 40, 190, 15, 50, 0.88,
                           "LONGIFOLIA 1/2-like C-terminal domain"))
    lines.append(scan_line("PF14309.20", 220, "ZLC04G0009000.1", 60, 1e-12, 50.0, 1,
                           1e-12, 1e-14, 49.0, 50, 195, 20, 59, 0.80,
                           "LONGIFOLIA 1/2-like C-terminal domain"))
    # 蛋白3 含 PF14383（检验用）
    lines.append(scan_line("PF14383.9", 90, "ZLC02G0005000.1", 60, 2e-8, 30.0, 1,
                           2e-8, 3e-9, 29.0, 5, 70, 3, 40, 0.75,
                           "DUF761-associated sequence motif"))
    # 蛋白4 一个无关域
    lines.append(scan_line("PF00001.30", 300, "ZLC03G0007000.1", 58, 1e-6, 20.0, 1,
                           1e-6, 1e-7, 19.0, 10, 100, 5, 35, 0.65,
                           "Some unrelated domain"))
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# 合成 BLAST 15 列结果
# ---------------------------------------------------------------------------
def write_blast(path, rows):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def blast_row(q, s, pident, evalue, bits, qcov, qlen=900, slen=900):
    return [q, s, pident, 300, 50, 2, 1, 300, 5, 304, evalue, bits, qlen, slen, qcov]


# ---------------------------------------------------------------------------
# 断言框架
# ---------------------------------------------------------------------------
class Checker:
    def __init__(self):
        self.n = 0
        self.fails = []

    def eq(self, got, want, label):
        self.n += 1
        if got != want:
            self.fails.append("%s\n      期望: %r\n      实际: %r" % (label, want, got))

    def true(self, cond, label):
        self.n += 1
        if not cond:
            self.fails.append(label)


def lines_of(path):
    if not os.path.exists(path):
        return []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        return [l.rstrip("\n") for l in fh if l.strip()]


def read_tsv(path):
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        header = None
        for i, l in enumerate(fh):
            f = l.rstrip("\n").split("\t")
            if i == 0:
                header = f
                continue
            if header and len(f) == len(header):
                rows.append(dict(zip(header, f)))
    return rows


def main():
    c = Checker()
    if os.path.isdir(TMP):
        shutil.rmtree(TMP)
    os.makedirs(TMP)

    pep = os.path.join(TMP, "Zunla.pep.fa")
    hmm_dom = os.path.join(TMP, "TRM.domtblout")
    scan_dom = os.path.join(TMP, "Pfam.domtblout")
    write_pep(pep)
    write_hmmsearch_domtbl(hmm_dom)
    write_hmmscan_domtbl(scan_dom)

    print("=" * 72)
    print("HMMER 解析脚本测试（合成数据）")
    print("=" * 72)

    # ---------------- 1. hmmsearch 解析 ----------------
    print("\n[1] parse_hmmsearch_domtbl.py")
    pre = os.path.join(TMP, "TRM.hmmsearch")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "parse_hmmsearch_domtbl.py"),
                        "--domtbl", hmm_dom, "--pep", pep, "--prefix", pre,
                        "--hmm-eval", "1e-5", "--incl-eval", "1e-2"],
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
    c.eq(r.returncode, 0, "解析器退出码")
    if r.returncode != 0:
        print(r.stdout, r.stderr)
        return report(c)

    best = read_tsv(pre + ".best_per_protein.tsv")
    by_prot = {x["protein"]: x for x in best}
    c.eq(len(best), 5, "每蛋白最佳表行数（应保留全部 5 个蛋白，含未达标者）")

    # strict: 全序列 E<=1e-5 -> ZLC01G0003110.1 / .2 / ZLC04G0009000.1
    strict = set(lines_of(pre + ".protein_ids.txt"))
    c.eq(strict, {"ZLC01G0003110.1", "ZLC01G0003110.2", "ZLC04G0009000.1"},
         "strict 蛋白集（E<=1e-5）")
    c.eq(set(lines_of(pre + ".gene_ids.txt")),
         {"ZLC01G0003110", "ZLC04G0009000"},
         "strict 基因集（isoform 已合并）")

    # inclusive: 加上 ZLC02G0005000.1（i-Evalue 5e-3 <= 1e-2）
    incl = set(lines_of(pre + ".protein_ids_inclusive.txt"))
    c.eq(incl, {"ZLC01G0003110.1", "ZLC01G0003110.2",
                "ZLC02G0005000.1", "ZLC04G0009000.1"},
         "inclusive 蛋白集（i-Evalue<=1e-2）")
    c.true("ZLC03G0007000.1" not in incl, "低于包含阈值的蛋白不应出现在 inclusive")

    # 最佳 domain 选择：ZLC01G0003110.1 应取 score 148 而非 120
    c.eq(by_prot["ZLC01G0003110.1"]["best_dom_score"], "148.0",
         "同一蛋白多 domain 应取最高分")
    c.eq(by_prot["ZLC01G0003110.1"]["significance"], "strict", "蛋白1 判定 strict")
    c.eq(by_prot["ZLC02G0005000.1"]["significance"], "inclusive_only",
         "蛋白3 判定 inclusive_only")
    c.eq(by_prot["ZLC03G0007000.1"]["significance"], "below_inclusive",
         "蛋白4 判定 below_inclusive")

    # 基因级：ZLC01G0003110 应有 2 个蛋白
    gl = read_tsv(pre + ".gene_level.tsv")
    by_gene = {x["gene"]: x for x in gl}
    c.eq(by_gene["ZLC01G0003110"]["n_proteins"], "2", "基因级 isoform 计数")
    c.eq(len(gl), 4, "基因级表行数")

    # 抽序列
    n_fa = sum(1 for l in lines_of(pre + ".hit_proteins.fasta") if l.startswith(">"))
    c.eq(n_fa, 5, "命中蛋白 FASTA 条数")

    # ---------------- 2. hmmscan 解析 ----------------
    print("\n[2] parse_hmmscan_domtbl.py")
    spre = os.path.join(TMP, "CaTRM_vs_Pfam")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "parse_hmmscan_domtbl.py"),
                        "--domtbl", scan_dom, "--prefix", spre,
                        "--dom-eval", "1e-5", "--focus", "PF14309,PF14383"],
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
    c.eq(r.returncode, 0, "hmmscan 解析器退出码")
    if r.returncode != 0:
        print(r.stdout, r.stderr)
        return report(c)

    pp = read_tsv(spre + ".per_protein.tsv")
    byp = {x["protein"]: x for x in pp}
    # 5 个蛋白都至少命中 1 个显著域（蛋白4 命中 PF00001 无关域），故应为 5 行
    c.eq(len(pp), 5, "每蛋白表行数（5 个蛋白都有显著域命中）")
    c.eq(byp["ZLC01G0003110.1"]["has_TRM_core_domain"], "yes", "蛋白1 含 PF14309")
    c.eq(byp["ZLC01G0003110.1"]["PF14309"], "1", "蛋白1 PF14309 标志=1")
    c.eq(byp["ZLC01G0003110.1"]["PF14383"], "0", "蛋白1 PF14383 标志=0")
    c.eq(byp["ZLC02G0005000.1"]["has_TRM_core_domain"], "no", "蛋白3 不含 PF14309")
    c.eq(byp["ZLC02G0005000.1"]["PF14383"], "1", "蛋白3 PF14383 标志=1")

    fm = read_tsv(spre + ".focus_matrix.tsv")
    # 矩阵只含「有任一关注域命中」的基因: ZLC01G0003110 / ZLC02G0005000 / ZLC04G0009000
    c.eq(len(fm), 3, "关注域矩阵行数（仅含命中关注域的基因）")
    fh = lines_of(spre + ".focus_hits.tsv")
    c.eq(len(fh) - 1, 4, "关注域命中明细条数（3×PF14309 + 1×PF14383）")

    # ---------------- 3. 合并 ----------------
    print("\n[3] merge_blast_hmm.py")
    blast_main = os.path.join(TMP, "TRM_candidates_main.tsv")
    # BLAST 支持 ZLC01G0003110、ZLC04G0009000（与 HMM 重叠）+ ZLC05G0001234（仅 BLAST）
    write_blast(blast_main, [
        blast_row("AT3G02170.1", "ZLC01G0003110.1", 45.0, 1e-30, 200.0, 70),
        blast_row("AT3G02170.1", "ZLC01G0003110.2", 44.0, 1e-28, 190.0, 68),
        blast_row("AT5G15580.1", "ZLC04G0009000.1", 38.0, 1e-20, 150.0, 60),
        blast_row("AT1G18620.2", "ZLC05G0001234.1", 35.0, 1e-10, 100.0, 55),
    ])
    bg = os.path.join(TMP, "blast_gene.txt")
    bp = os.path.join(TMP, "blast_prot.txt")
    with open(bg, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("ZLC01G0003110\nZLC04G0009000\nZLC05G0001234\n")
    with open(bp, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("ZLC01G0003110.1\nZLC01G0003110.2\nZLC04G0009000.1\nZLC05G0001234.1\n")

    anno = os.path.join(TMP, "anno.tsv")
    with open(anno, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("gene\trepresentative_transcript\tgene_id\tsymbols\tdescription\t"
                 "length_aa\tisoform_count\tselected_is_longest\tnote\n")
        fh.write("At3g02170\tAT3G02170.1\tAT3G02170\tTRM1,LNG2\tdesc\t905\t3\tyes\t\n")
        fh.write("At5g15580\tAT5G15580.1\tAT5G15580\tTRM2,LNG1\tdesc\t927\t1\tyes\t\n")
        fh.write("At1g18620\tAT1G18620.2\tAT1G18620\tTRM3\tdesc\t1014\t5\tyes\t\n")

    mpre = os.path.join(TMP, "CaTRM_candidates")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "merge_blast_hmm.py"),
                        "--blast-gene", bg, "--blast-prots", bp, "--blast-main", blast_main,
                        "--hmm-gene", pre + ".gene_level.tsv",
                        "--hmm-gene-strict", pre + ".gene_ids.txt",
                        "--hmm-gene-incl", pre + ".gene_ids_inclusive.txt",
                        "--anno", anno, "--prefix", mpre],
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
    c.eq(r.returncode, 0, "合并脚本退出码")
    if r.returncode != 0:
        print(r.stdout, r.stderr)
        return report(c)

    mrows = read_tsv(mpre + ".merged.tsv")
    byg = {x["gene"]: x for x in mrows}
    c.eq(len(mrows), 4, "合并后基因总数（并集）")
    c.eq(byg["ZLC01G0003110"]["evidence"], "BLAST+HMM", "ZLC01G0003110 -> BLAST+HMM")
    c.eq(byg["ZLC04G0009000"]["evidence"], "BLAST+HMM", "ZLC04G0009000 -> BLAST+HMM")
    c.eq(byg["ZLC05G0001234"]["evidence"], "BLAST-only", "ZLC05G0001234 -> BLAST-only")
    c.eq(byg["ZLC02G0005000"]["evidence"], "HMM-inclusive-only",
         "ZLC02G0005000 -> HMM-inclusive-only")
    c.eq(byg["ZLC01G0003110"]["at_genes"], "At3g02170",
         "At 基因映射（沿用注释表里的写法）")
    c.eq(byg["ZLC01G0003110"]["at_symbols"], "TRM1,LNG2", "TRM 符号映射")
    c.eq(byg["ZLC01G0003110"]["blast_best_pident"], "45.0", "最佳 identity")

    ev = {x["evidence"]: x["n_genes"] for x in read_tsv(mpre + ".evidence_counts.tsv")}
    c.eq(ev.get("BLAST+HMM"), "2", "证据计数 BLAST+HMM")
    c.eq(ev.get("BLAST-only"), "1", "证据计数 BLAST-only")
    c.eq(ev.get("HMM-inclusive-only"), "1", "证据计数 HMM-inclusive-only")

    # ---------------- 4. 最终汇总 ----------------
    print("\n[4] final_summary.py")
    fpre = os.path.join(TMP, "CaTRM_final")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "final_summary.py"),
                        "--merged", mpre + ".merged.tsv",
                        "--per-protein", spre + ".per_protein.tsv",
                        "--focus-matrix", spre + ".focus_matrix.tsv",
                        "--gene-level", pre + ".gene_level.tsv",
                        "--pep", pep,
                        "--prefix", fpre],
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
    c.eq(r.returncode, 0, "汇总脚本退出码")
    if r.returncode != 0:
        print(r.stdout, r.stderr)
        return report(c)

    final = read_tsv(fpre + ".gene_list.tsv")
    byf = {x["gene"]: x for x in final}
    c.eq(len(final), 4, "最终基因列表行数")
    c.eq(byf["ZLC01G0003110"]["tier"], "conservative",
         "双方法+PF14309 -> conservative")
    c.eq(byf["ZLC04G0009000"]["tier"], "conservative",
         "双方法+PF14309 -> conservative")
    c.eq(byf["ZLC05G0001234"]["tier"], "review", "仅BLAST -> review")
    c.eq(byf["ZLC02G0005000"]["tier"], "review", "仅inclusive -> review")
    c.eq(byf["ZLC01G0003110"]["has_PF14309"], "yes", "最终表 PF14309 标志")

    ev_rows = read_tsv(fpre + ".evidence.tsv")
    c.eq(len(ev_rows), 4, "证据表行数（每成员一条理由）")
    c.true(all(x["selection_reasons"].strip() for x in ev_rows),
           "每个成员都必须有筛选理由")
    rev = read_tsv(fpre + ".review.tsv")
    c.eq(len(rev), 2, "review 层应为 2 个（仅BLAST + 仅inclusive）")

    # ZLC05G0001234 只有 BLAST 命中，其蛋白不在本测试的合成蛋白库里，
    # 因此最终 FASTA 只含 3 条（该基因仍出现在 gene_list 中并标 review）。
    n_final_fa = sum(1 for l in lines_of(fpre + ".proteins.fasta") if l.startswith(">"))
    c.eq(n_final_fa, 3, "最终代表蛋白 FASTA 条数（缺序列的基因跳过）")
    c.true(all(x["gene"] == "ZLC05G0001234" or x["gene"] in
               {"ZLC01G0003110", "ZLC02G0005000", "ZLC04G0009000"} for x in final),
           "所有候选都应出现在最终基因列表")

    # ---------------- 5. 降级路径 ----------------
    print("\n[5] 降级：缺 BLAST / 缺结构域")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "merge_blast_hmm.py"),
                        "--blast-gene", "", "--blast-prots", "", "--blast-main", "",
                        "--hmm-gene", pre + ".gene_level.tsv",
                        "--hmm-gene-strict", pre + ".gene_ids.txt",
                        "--hmm-gene-incl", pre + ".gene_ids_inclusive.txt",
                        "--anno", "", "--prefix", os.path.join(TMP, "noblast")],
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
    c.eq(r.returncode, 0, "缺 BLAST 时合并脚本应正常退出")
    nb = read_tsv(os.path.join(TMP, "noblast.merged.tsv"))
    c.eq(len(nb), 3, "缺 BLAST 时应只有 HMM 候选（3 基因）")
    c.true(all(x["evidence"] in ("HMM-only", "HMM-inclusive-only") for x in nb),
           "缺 BLAST 时证据标签应为 HMM-only / HMM-inclusive-only")

    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "final_summary.py"),
                        "--merged", mpre + ".merged.tsv",
                        "--per-protein", "", "--gene-level", "",
                        "--pep", pep, "--prefix", os.path.join(TMP, "nodom")],
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
    c.eq(r.returncode, 0, "缺结构域结果时汇总脚本应正常退出")
    nd = read_tsv(os.path.join(TMP, "nodom.gene_list.tsv"))
    c.eq(len(nd), 4, "缺结构域时仍应输出全部候选")
    c.true(all(x["has_PF14309"] == "NA" for x in nd), "缺结构域时 has_PF14309 应为 NA")
    c.eq(sum(1 for x in nd if x["tier"] == "conservative"), 0,
         "缺结构域时不应有 conservative（无法确认核心域）")

    return report(c)


def report(c):
    print("\n" + "=" * 72)
    print("断言数: %d   失败: %d" % (c.n, len(c.fails)))
    if c.fails:
        print("\n失败明细:")
        for f in c.fails:
            print("  [FAIL] %s" % f)
        return 1
    print("全部通过 ✅")
    shutil.rmtree(TMP, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
