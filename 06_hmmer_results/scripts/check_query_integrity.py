# -*- coding: utf-8 -*-
"""
check_query_integrity.py —— 34 条 AtTRM query 的完整性校验（离线，不需服务器）

为什么做
    TRM.hmm 由这 34 条序列建成。query 若有缺残基/错序列/非标准字符，
    HMM 的质量会整体受损，而且是后续所有步骤都无法发现的问题。
    因此必须在下游浪费算力之前先验证。

验证项
    1. 条数是否为 34
    2. 与 Araport11 源蛋白组**往返比对**：每条 query 的序列必须与
       Araport11_pep_20250411.gz 中同 ID 的序列逐残基一致
    3. 长度与注释表（02_At_TRM_annotation.tsv）核对
    4. 字母表检查：仅标准氨基酸，无 *、无 X/B/Z/U
    5. ID 与注释表一一对应（无多、无缺）

产出：屏幕报告；退出码 0 = 全部通过
"""
import gzip
import os
import re
import sys

PROJECT = (r"E:\华为\deepseek harness\gene_famaily_analysis"
           r"\Identification_of_gene_family_members")
QUERY = os.path.join(PROJECT, "02_script_output", "01_At_TRM_query_34.fasta")
ANNO = os.path.join(PROJECT, "02_script_output", "02_At_TRM_annotation.tsv")
ARAPORT = os.path.join(PROJECT, "01_At_TRM_query", "Araport11_pep_20250411.gz")

STD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_fasta(path, gz=False):
    """返回 (ids_in_order, {id: seq})"""
    opener = gzip.open if gz else open
    ids, seqs = [], {}
    cur, buf = None, []
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf)
                cur = line[1:].strip().split()[0]
                ids.append(cur)
                buf = []
            else:
                buf.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(buf)
    return ids, seqs


class Checker:
    def __init__(self):
        self.n = 0
        self.fails = []

    def chk(self, cond, label, detail=""):
        self.n += 1
        print("  [%s] %s%s" % ("OK  " if cond else "FAIL", label,
                               ("  " + detail) if detail else ""))
        if not cond:
            self.fails.append(label)


def main():
    for p in (QUERY, ANNO, ARAPORT):
        if not os.path.exists(p):
            print("缺少输入: %s" % p)
            return 2

    c = Checker()
    print("=" * 74)
    print("34 条 AtTRM query 完整性校验")
    print("=" * 74)

    q_ids, q_seqs = read_fasta(QUERY)
    a_ids, a_seqs = read_fasta(ARAPORT, gz=True)

    print("\n[1] 条数")
    c.chk(len(q_ids) == 34, "query 条数 = 34", "实际 %d" % len(q_ids))
    c.chk(len(q_ids) == len(set(q_ids)), "query ID 无重复")

    print("\n[2] 与 Araport11 源蛋白组往返比对（逐残基）")
    n_match = n_miss = n_diff = 0
    diffs = []
    for tid in q_ids:
        ref = a_seqs.get(tid)
        if ref is None:
            n_miss += 1
            diffs.append("%s: Araport11 中不存在此 ID" % tid)
            continue
        if ref == q_seqs[tid]:
            n_match += 1
        else:
            n_diff += 1
            diffs.append("%s: 长度或残基不一致 (query %d vs Araport11 %d)"
                         % (tid, len(q_seqs[tid]), len(ref)))
    c.chk(n_match == len(q_ids),
          "全部 query 与 Araport11 逐残基一致",
          "%d/%d 一致, %d 不一致, %d 缺失" % (n_match, len(q_ids), n_diff, n_miss))
    for d in diffs[:10]:
        print("        %s" % d)

    print("\n[3] 长度与注释表核对")
    anno = {}
    with open(ANNO, "r", encoding="utf-8", errors="replace") as fh:
        header = None
        for i, line in enumerate(fh):
            f = line.rstrip("\n").split("\t")
            if i == 0:
                header = f
                continue
            if header and len(f) == len(header):
                anno[f[1].strip()] = int(f[5])
    n_len_ok = sum(1 for t in q_ids if anno.get(t) == len(q_seqs[t]))
    c.chk(n_len_ok == len(q_ids),
          "query 长度与注释表一致", "%d/%d" % (n_len_ok, len(q_ids)))

    print("\n[4] ID 一致性")
    c.chk(set(q_ids) == set(anno.keys()),
          "query ID 集合与注释表完全相同",
          "query %d 个, 注释表 %d 个" % (len(q_ids), len(anno)))

    print("\n[5] 字母表检查")
    bad_alpha = {}
    for tid in q_ids:
        s = q_seqs[tid]
        bad = sorted(set(s) - STD_AA)
        if bad:
            bad_alpha[tid] = bad
    c.chk(not bad_alpha, "全部序列仅含 20 种标准氨基酸",
          "异常: %s" % bad_alpha if bad_alpha else "")
    lens = [len(q_seqs[t]) for t in q_ids]
    c.chk(all(l > 0 for l in lens), "无空序列")
    print("        长度范围: %d - %d aa，平均 %.1f aa"
          % (min(lens), max(lens), sum(lens) / len(lens)))

    print("\n" + "=" * 74)
    print("断言数: %d   失败: %d" % (c.n, len(c.fails)))
    if c.fails:
        print("失败项:")
        for f in c.fails:
            print("  - %s" % f)
        return 1
    print("全部通过 ✅  query 可作为 hmmbuild 的可靠输入")
    return 0


if __name__ == "__main__":
    sys.exit(main())
