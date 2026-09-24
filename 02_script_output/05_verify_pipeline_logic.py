#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
本地验证 trm_blast_pipeline.sh 的核心逻辑（阶段 3/4/6）。

背景
    Windows 沙箱禁止 Git Bash（无法创建 signal pipe），所以无法直接 `bash -n`。
    本脚本用 Python 精确复刻 pipeline 里各条 awk/sed/cut 的行为，
    并用含边界情况的合成 BLAST 结果验证，确保逻辑正确、阈值边界无歧义。

运行
    python scripts/tests/verify_pipeline_logic.py
退出码 0 = 全部通过。
"""
import os
import re
import sys

# ---------------- 与 pipeline 保持一致的参数 ----------------
EVALUE = 1e-5
IDENT = 30.0
COV = 50.0
HIGH_COV = 70.0
HIGH_E = 1e-10

# 15 列 outfmt 顺序
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send
# evalue bitscore qlen slen qcovhsp
F = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
     "qstart", "qend", "sstart", "send", "evalue", "bitscore",
     "qlen", "slen", "qcovhsp"]


def row(q, s, pident, evalue, qcov, qlen=900, slen=900, bits=500.0):
    """构造一行 15 列 BLAST TSV"""
    return "\t".join([
        q, s, "%.2f" % pident, "300", "50", "2", "1", "300", "5", "304",
        "%.2e" % evalue, "%.1f" % bits, str(qlen), str(slen), "%.0f" % qcov,
    ])


# ---------------- 合成边界测试数据 ----------------
ROWS = [
    # --- 主候选（应通过） ---
    row("AT3G02170.1", "ZLC01G0000010.1", 78.5, 0.0,    95),   # 高可信
    row("AT3G02170.1", "ZLC01G0000020.2", 55.0, 1e-50,  80),   # 主候选
    # --- 边界：恰好等于阈值（应通过，>= 语义） ---
    row("AT5G15580.1", "ZLC02G0000100.1", 30.0, 1e-5,   50),   # 恰在 I/C/E 三边界
    # --- 边界：刚刚不达阈值（应排除） ---
    row("AT1G18620.2", "ZLC03G0000200.1", 29.99, 1e-5,  50),   # I 差一点
    row("AT1G18620.2", "ZLC03G0000300.1", 60.0, 1e-5,   49),   # C 差一点
    row("AT1G18620.2", "ZLC03G0000400.1", 60.0, 1e-4,   90),   # E 差一点(1e-4>1e-5)
    # --- 片段待查（E 达标但 C 不足） ---
    row("AT1G74160.1", "ZLC04G0000500.1", 85.0, 1e-30,  22),
    row("AT1G74160.1", "ZLC04G0000600.1", 70.0, 1e-8,   15),
    # --- 高可信档（C>=70 且 E<=1e-10） ---
    row("AT3G63430.1", "ZLC05G0000700.1", 42.0, 1e-100, 88),
    # --- 同一 query 多个 isoform -> 去重应合并到 1 个基因 ---
    row("AT3G05750.1", "ZLC06G0000800.1", 66.0, 1e-40,  75),
    row("AT3G05750.1", "ZLC06G0000800.2", 64.0, 1e-38,  73),
    # --- 一个基因的多个蛋白命中不同 query ---
    row("AT3G58650.1", "ZLC06G0000800.1", 61.0, 1e-20,  70),
    # --- 零命中 query: AT5G26910.1 不出现在数据里 ---
]

QUERY_IDS = [
    "AT3G02170.1", "AT5G15580.1", "AT1G18620.2", "AT1G74160.1",
    "AT3G63430.1", "AT3G05750.1", "AT3G58650.1",
    "AT5G26910.1",   # 零命中
]


# ---------------- 复刻 pipeline 的逻辑 ----------------
def parse(rows):
    out = []
    for r in rows:
        f = r.split("\t")
        out.append({
            "qseqid": f[0],
            "sseqid": f[1],
            "pident": float(f[2]),
            "evalue": float(f[10]),
            "bitscore": float(f[11]),
            "qcovhsp": float(f[14]),
        })
    return out


def stage3_main(hits):
    """阶段 3: E<=EVALUE 且 pident>=IDENT 且 qcovhsp>=COV"""
    return [h for h in hits
            if float(h["evalue"]) <= EVALUE + 1e-300
            and h["pident"] >= IDENT
            and h["qcovhsp"] >= COV]


def stage3_high(hits):
    return [h for h in hits if h["qcovhsp"] >= HIGH_COV and float(h["evalue"]) <= HIGH_E]


def stage3_frag(hits):
    return [h for h in hits if float(h["evalue"]) <= EVALUE + 1e-300 and h["qcovhsp"] < COV]


def strip_iso(s):
    """复刻 sed 's/\\.[0-9]\\+$//'"""
    return re.sub(r"\.\d+$", "", s)


def stage4_ids(main_hits):
    prot = sorted({h["sseqid"] for h in main_hits})
    gene = sorted({strip_iso(p) for p in prot})
    return prot, gene


def stage6_zero_hit(hits, all_queries):
    hit_q = {h["qseqid"] for h in hits}
    return sorted(set(all_queries) - hit_q)


def stage6_best_per_query(hits):
    """复刻: sort -k1,1 -k12,12gr | awk '!seen[$1]++'（按 bitscore 降序取首个）"""
    best = {}
    for h in sorted(hits, key=lambda x: x["qseqid"]):
        pass
    for h in sorted(hits, key=lambda x: (-x["bitscore"],)):
        best.setdefault(h["qseqid"], h)
    return best


# ---------------- 断言框架 ----------------
class Checker:
    def __init__(self):
        self.fails = []
        self.n = 0

    def eq(self, got, want, label):
        self.n += 1
        if got != want:
            self.fails.append("%s\n      期望: %r\n      实际: %r" % (label, want, got))

    def true(self, cond, label):
        self.n += 1
        if not cond:
            self.fails.append(label)


def write_testdata(outdir):
    """
    把合成数据落盘, 供 01_summarize_candidates.py 做集成测试。
    产出 outdir/TRM_vs_Zunla.blast.tsv 与 outdir/TRM_34_test.fasta
    """
    os.makedirs(outdir, exist_ok=True)
    bl = os.path.join(outdir, "TRM_vs_Zunla.blast.tsv")
    with open(bl, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(ROWS) + "\n")
    fa = os.path.join(outdir, "TRM_34_test.fasta")
    with open(fa, "w", encoding="utf-8", newline="\n") as fh:
        for q in QUERY_IDS:
            fh.write(">%s | gene=%s | Symbols:TEST\n" % (q, strip_iso(q)))
            fh.write("MSAKLLYNLSDENPNLNKQFGCMNGIFQVFYRQHCPATPVTVSGGAEKSLPPGERRGSVG\n")
    return bl, fa


def main():
    if len(sys.argv) > 2 and sys.argv[1] == "--write-testdata":
        bl, fa = write_testdata(sys.argv[2])
        print("测试数据已写出:")
        print("  %s" % bl)
        print("  %s" % fa)
        return 0

    c = Checker()
    hits = parse(ROWS)

    # ---- 阶段 3: 主候选过滤 ----
    # 注意: 用 (query, sseqid) 作键——我故意让 ZLC06G0000800.1 被两个 query 命中,
    # 只按 sseqid 去重会看不出这种多对多情形。
    main_hits = stage3_main(hits)
    got = sorted((h["qseqid"], h["sseqid"]) for h in main_hits)
    want = sorted([
        ("AT3G02170.1", "ZLC01G0000010.1"),
        ("AT3G02170.1", "ZLC01G0000020.2"),
        ("AT5G15580.1", "ZLC02G0000100.1"),
        ("AT3G63430.1", "ZLC05G0000700.1"),
        ("AT3G05750.1", "ZLC06G0000800.1"),
        ("AT3G05750.1", "ZLC06G0000800.2"),
        ("AT3G58650.1", "ZLC06G0000800.1"),   # 同一辣椒蛋白被第二个 query 命中
    ])
    c.eq(got, want, "阶段3 主候选集 (含恰在阈值上的 ZLC02G0000100.1)")

    # 边界必须正确
    c.true(any(s == "ZLC02G0000100.1" for _, s in got),
           "I==30/C==50/E==1e-5 应被纳入 (>= 语义)")
    c.true(not any(s == "ZLC03G0000200.1" for _, s in got), "identity 29.99 应被排除")
    c.true(not any(s == "ZLC03G0000300.1" for _, s in got), "coverage 49 应被排除")
    c.true(not any(s == "ZLC03G0000400.1" for _, s in got), "evalue 1e-4 应被排除 (>1e-5)")

    # 多对多: 同一辣椒蛋白可出现在多个 query 下
    c.eq(sum(1 for _, s in got if s == "ZLC06G0000800.1"), 2,
         "ZLC06G0000800.1 应被 2 个不同 query 命中")

    # ---- 阶段 3: 三档划分 ----
    high = stage3_high(hits)
    c.eq(sorted((h["qseqid"], h["sseqid"]) for h in high),
         sorted([
             ("AT3G02170.1", "ZLC01G0000010.1"),   # C=95,  E=0
             ("AT3G02170.1", "ZLC01G0000020.2"),   # C=80,  E=1e-50
             ("AT3G63430.1", "ZLC05G0000700.1"),   # C=88,  E=1e-100
             ("AT3G05750.1", "ZLC06G0000800.1"),   # C=75,  E=1e-40
             ("AT3G05750.1", "ZLC06G0000800.2"),   # C=73,  E=1e-38
             ("AT3G58650.1", "ZLC06G0000800.1"),   # C=70,  E=1e-20 (恰好 C=70)
         ]),
         "阶段3 高可信档 (C>=70 且 E<=1e-10)")

    frag = stage3_frag(hits)
    c.eq(sorted(h["sseqid"] for h in frag),
         ["ZLC03G0000300.1", "ZLC04G0000500.1", "ZLC04G0000600.1"],
         "阶段3 片段待查档 (E 达标但 C<50)")

    # 三档互斥性: 片段档与主候选档不得重叠
    c.true(not ({h["sseqid"] for h in frag} & {h["sseqid"] for h in main_hits}),
           "片段档与主候选档不应重叠（阈值按 C 划分）")

    # ---- 阶段 4: 去重 ----
    prot, gene = stage4_ids(main_hits)
    c.eq(prot, ["ZLC01G0000010.1", "ZLC01G0000020.2", "ZLC02G0000100.1",
                "ZLC05G0000700.1", "ZLC06G0000800.1", "ZLC06G0000800.2"],
         "阶段4 蛋白级 ID（保留 isoform）")
    c.eq(gene, ["ZLC01G0000010", "ZLC01G0000020", "ZLC02G0000100",
                "ZLC05G0000700", "ZLC06G0000800"],
         "阶段4 基因级 ID（isoform 后缀已剥离）")
    c.eq(len(prot), 6, "蛋白级计数")
    c.eq(len(gene), 5, "基因级计数（ZLC06G0000800 的 .1/.2 合并）")

    # 编号 >9 的 isoform 也要正确剥离
    c.eq(strip_iso("ZLC01G0000010.12"), "ZLC01G0000010", "剥离两位数 isoform")
    c.eq(strip_iso("ZLC01G0000010"), "ZLC01G0000010", "无后缀时保持不变")
    c.eq(strip_iso("AT3G02170.1"), "AT3G02170", "At 转录本同样适用")

    # ---- 阶段 6: 零命中 ----
    zero = stage6_zero_hit(hits, QUERY_IDS)
    c.eq(zero, ["AT5G26910.1"], "阶段6 零命中 query 识别")

    # ---- 阶段 6: 每 query 最佳命中 ----
    best = stage6_best_per_query(hits)
    c.eq(best["AT3G02170.1"]["sseqid"], "ZLC01G0000010.1",
         "最佳命中按 bitscore 取 (78.5 > 55.0)")
    c.eq(best["AT3G05750.1"]["sseqid"], "ZLC06G0000800.1",
         "最佳命中在同基因两 isoform 间取高分者")
    c.eq(len(best), len({h["qseqid"] for h in hits}),
         "最佳命中数 == 有命中的 query 数")

    # ---- 反例: evalue 的浮点解析 ----
    c.true(float("1e-5") <= EVALUE, "1e-5 应判定为达标")
    c.true(not (float("1.000001e-5") <= EVALUE), "1.000001e-5 应判定为不达标")

    # ---- 报告 ----
    print("=" * 66)
    print("trm_blast_pipeline.sh 核心逻辑验证")
    print("=" * 66)
    print("合成 HSP 行数      : %d" % len(hits))
    print("主候选 HSP         : %d" % len(main_hits))
    print("候选蛋白/基因      : %d / %d" % (len(prot), len(gene)))
    print("高可信 / 片段待查  : %d / %d" % (len(high), len(frag)))
    print("零命中 query       : %s" % ", ".join(zero))
    print("-" * 66)
    print("断言数: %d   失败: %d" % (c.n, len(c.fails)))
    if c.fails:
        print("\n失败明细:")
        for f in c.fails:
            print("  [FAIL] %s" % f)
        return 1
    print("全部通过 ✅")
    return 0


if __name__ == "__main__":
    sys.exit(main())
