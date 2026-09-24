#!/usr/bin/env bash
# =============================================================================
# trm_blast_pipeline.sh
#   拟南芥 TRM 家族(34 条代表蛋白)比对 Zunla 辣椒蛋白组 -> 鉴定候选家族
#
# 用法:
#   bash trm_blast_pipeline.sh \
#       --query   query/TRM_34.fasta \
#       --db-fasta db/Canz.pep.fa \
#       --outdir  out \
#       --threads 8 \
#       --evalue  1e-5 \
#       --ident   30 \
#       --cov     50
#
# 阈值(采纳方案定稿):
#   E-value <= 1e-5 且 identity >= 30% 且 query coverage >= 50%
#   理由: TRM 富含重复序列, identity 天然偏低(卡 50% 会漏真成员);
#         coverage 才是主力过滤器(排除"仅局部重复区匹配"的假阳性)。
#
# 产出:
#   ① out/TRM_vs_Zunla.blast.tsv                     原始 BLAST 结果
#   ② out/TRM_candidates_protein_level.txt            去重候选蛋白 ID(含 isoform)
#      out/TRM_candidates_gene_level.txt               去重候选基因 ID(剥离后缀)
#   ③ out/Zunla_TRM_candidate_proteins.fasta          去重候选蛋白 FASTA
#   辅助 out/01_*.tsv|txt|png, out/Zunla_TRM_candidate_proteins_genelevel.fasta
# =============================================================================
set -euo pipefail

VERSION="1.0"

# ------------------------------ 默认参数 -------------------------------------
QUERY=""
DB_FASTA=""
OUTDIR="out"
THREADS=4
EVALUE="1e-5"
IDENT=30
COV=50
MAX_TARGET=50
DB_PREFIX=""
FORCE=0
MAPFILE=""          # 可选: At 基因映射表(推荐传 TRM_representative_proteins.tsv)

# ------------------------------ 工具函数 -------------------------------------
log()  { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }
step() { printf '\n========== 阶段 %s: %s ==========\n' "$1" "$2"; }
die()  { printf '\n[错误] %s\n' "$*" >&2; exit 1; }
warn() { printf '[警告] %s\n' "$*" >&2; }

usage() {
  # 打印脚本头部注释块（第 2-27 行）作为帮助；不依赖行号以外的东西
  sed -n '2,27p' "$0" | sed 's/^# \{0,1\}//'
  exit "${1:-0}"
}

# ------------------------------ 参数解析 -------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --query)       QUERY="$2";       shift 2 ;;
    --db-fasta)    DB_FASTA="$2";    shift 2 ;;
    --outdir)      OUTDIR="$2";      shift 2 ;;
    --threads)     THREADS="$2";     shift 2 ;;
    --evalue)      EVALUE="$2";      shift 2 ;;
    --ident)       IDENT="$2";       shift 2 ;;
    --cov)         COV="$2";         shift 2 ;;
    --max-target)  MAX_TARGET="$2";  shift 2 ;;
    --db-prefix)   DB_PREFIX="$2";   shift 2 ;;
    --map)         MAPFILE="$2";     shift 2 ;;
    --force)       FORCE=1;          shift ;;
    -h|--help)     usage 0 ;;
    *) usage 1 >&2; echo "未知参数: $1 (用 --help 查看用法)" >&2; exit 2 ;;
  esac
done

[[ -n "$QUERY"    ]] || die "缺少 --query"
[[ -n "$DB_FASTA" ]] || die "缺少 --db-fasta"

mkdir -p "$OUTDIR"
DB_PREFIX="${DB_PREFIX:-${OUTDIR}/Canz.pep}"

BL_TSV="${OUTDIR}/TRM_vs_Zunla.blast.tsv"
DIAG="${OUTDIR}/01_diagnostics.txt"

# 诊断日志同时写屏幕和文件。
# 主流程放在 run_pipeline() 里, 最后用 `run_pipeline | tee` 输出。
# 这样 stdout/stderr 全程落盘又上屏, 且不像 exec >(tee ...) 那样有末尾截断风险。
run_pipeline() {

log "trm_blast_pipeline.sh v${VERSION}"
log "query=$QUERY  db=$DB_FASTA  outdir=$OUTDIR  threads=$THREADS"
log "阈值: E<=$EVALUE  identity>=$IDENT%  coverage>=$COV%"

# =============================================================================
step 0 "输入校验与环境检查"
# =============================================================================
[[ -f "$QUERY"    ]] || die "找不到 query 文件: $QUERY"
[[ -f "$DB_FASTA" ]] || die "找不到蛋白库文件: $DB_FASTA"

command -v blastp >/dev/null 2>&1     || die "blastp 未找到。请先加载环境:
    module load blast+/2.14.0     或
    conda activate blast          或
    下载 NCBI 静态包后 export PATH=<...>/ncbi-blast-*/bin:\$PATH"
command -v makeblastdb >/dev/null 2>&1 || die "makeblastdb 未找到(同上处理)"

# BLAST+ 版本必须 >= 2.9, 否则 outfmt 不支持 qcovhsp
BL_VERSION_RAW="$(blastp -version 2>&1 | head -1)"
BL_VER="$(printf '%s' "$BL_VERSION_RAW" | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"
[[ -n "$BL_VER" ]] || die "无法解析 BLAST 版本: $BL_VERSION_RAW"
BL_MAJ="${BL_VER%%.*}"
BL_MIN="$(printf '%s' "$BL_VER" | cut -d. -f2)"
if (( BL_MAJ < 2 )) || { (( BL_MAJ == 2 )) && (( BL_MIN < 9 )); }; then
  die "BLAST+ 版本 $BL_VER 过旧, outfmt 不支持 qcovhsp(需 >= 2.9)。请升级后重跑。"
fi
log "BLAST+ 版本: $BL_VER  (满足 >= 2.9)"

N_QUERY="$(grep -c '^>' "$QUERY" || true)"
N_DB="$(grep -c '^>' "$DB_FASTA" || true)"
log "query 序列数: $N_QUERY"
log "db    序列数: $N_DB"
(( N_QUERY > 0 )) || die "query 文件里没有序列(不以 '>' 开头?)"
(( N_DB    > 0 )) || die "db 文件里没有序列"

if (( N_QUERY != 34 )); then
  warn "query 不是预期的 34 条(实际 $N_QUERY)。若你确实用的是 34 条代表蛋白请检查文件。"
fi
if (( N_DB != 52385 )); then
  warn "db 不是预期的 52385 条(实际 $N_DB)。若参考蛋白组不同属正常。"
fi

# 二进制文件检查(防止把 SnapGene 导出的 .prot 当 FASTA 用)
if command -v file >/dev/null 2>&1; then
  QT="$(file -b "$QUERY" || true)"
  case "$QT" in
    *text*|*ASCII*|*Unicode*|*UTF-8*) : ;;
    *) die "query 看起来不是纯文本 FASTA(file: $QT)。
    如果你是 SnapGene 导出的文件, 请重新用 File -> Export -> Protein FASTA 导出。" ;;
  esac
  log "query 文件类型: $QT"
fi

# =============================================================================
step 1 "建立 BLAST 蛋白库"
# =============================================================================
if [[ -f "${DB_PREFIX}.pdb" || -f "${DB_PREFIX}.phr" || -f "${DB_PREFIX}.psq" ]]; then
  if (( FORCE == 1 )); then
    log "已存在库, --force 指定, 重新建库"
    rm -f "${DB_PREFIX}".{pdb,phr,pin,psq,pot,ptf,pto,pjs}
    makeblastdb -in "$DB_FASTA" -dbtype prot -out "$DB_PREFIX" \
      -title "Zunla_capsicum_pep" -parse_seqids
  else
    log "已存在库 ${DB_PREFIX}.* , 跳过建库(需要重建请加 --force)"
  fi
else
  log "建库中..."
  if ! makeblastdb -in "$DB_FASTA" -dbtype prot -out "$DB_PREFIX" \
        -title "Zunla_capsicum_pep" -parse_seqids; then
    warn "-parse_seqids 建库失败(可能 SeqId 重复), 回退为不解析 SeqId"
    rm -f "${DB_PREFIX}".{pdb,phr,pin,psq,pot,ptf,pto,pjs}
    makeblastdb -in "$DB_FASTA" -dbtype prot -out "$DB_PREFIX" \
      -title "Zunla_capsicum_pep"
  fi
fi

DB_INFO="$(blastdbcmd -db "$DB_PREFIX" -info 2>&1 | head -3 || true)"
log "库信息: $(printf '%s' "$DB_INFO" | tr '\n' ' ')"

# =============================================================================
step 2 "第一轮 blastp -> ① 原始 BLAST 结果"
# =============================================================================
log "运行 blastp (线程 $THREADS) ..."
blastp \
  -query "$QUERY" \
  -db    "$DB_PREFIX" \
  -out   "$BL_TSV" \
  -evalue "$EVALUE" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp" \
  -max_target_seqs "$MAX_TARGET" \
  -num_threads "$THREADS" \
  -seg no

N_HSP="$(wc -l < "$BL_TSV" | tr -d ' ')"
log "① 原始结果: $BL_TSV  ($N_HSP 行 HSP)"
(( N_HSP > 0 )) || die "BLAST 结果为空。请检查 query 是否为蛋白序列、db 是否为蛋白库。"

# =============================================================================
step 3 "过滤 + 分档 -> A1 候选汇总表"
# =============================================================================
MAIN_TSV="${OUTDIR}/TRM_candidates_main.tsv"
HIGH_TSV="${OUTDIR}/TRM_candidates_highconf.tsv"
FRAG_TSV="${OUTDIR}/TRM_candidates_fragment.tsv"
SUMMARY="${OUTDIR}/01_candidate_summary.tsv"

# 主候选: E<=EVALUE 且 pident>=IDENT 且 qcovhsp>=COV
awk -F'\t' -v OFS='\t' -v E="$EVALUE" -v I="$IDENT" -v C="$COV" \
  '($11+0) <= (E+0) && ($3+0) >= I && ($15+0) >= C' "$BL_TSV" > "$MAIN_TSV"

# 高可信档: cov>=70 且 E<=1e-10
awk -F'\t' -v OFS='\t' '($15+0) >= 70 && ($11+0) <= 1e-10' "$BL_TSV" > "$HIGH_TSV"

# 片段待查档: 满足 E 值但覆盖度不足
awk -F'\t' -v OFS='\t' -v E="$EVALUE" -v C="$COV" \
  '($11+0) <= (E+0) && ($15+0) < C' "$BL_TSV" > "$FRAG_TSV"

N_MAIN="$(wc -l < "$MAIN_TSV" | tr -d ' ')"
N_HIGH="$(wc -l < "$HIGH_TSV" | tr -d ' ')"
N_FRAG="$(wc -l < "$FRAG_TSV" | tr -d ' ')"
log "主候选 HSP: $N_MAIN   高可信 HSP: $N_HIGH   片段待查 HSP: $N_FRAG"
(( N_MAIN > 0 )) || die "主候选集为空(阈值可能过严)。请查看 $BL_TSV 的 pident/qcovhsp 分布后调 --ident/--cov/--evalue。"

# 汇总表: 拼上 At 基因号与 TRM 符号(可选映射表)
# 映射表格式(scripts/02 产出): gene  representative_transcript  gene_id  symbols ...
{
  printf 'query\ttrm_gene\ttrm_symbol\tzunla_protein\tzunla_gene\tpident\tevalue\tbitscore\tqcovhsp\ttier\n'
  if [[ -n "$MAPFILE" && -f "$MAPFILE" ]]; then
    awk -F'\t' -v OFS='\t' -v E="$EVALUE" -v I="$IDENT" -v C="$COV" '
      NR==FNR { if (FNR>1) { g[$2]=$1; s[$2]=$4 } next }
      {
        q=$1
        tier = (($15+0)>=70 && ($11+0)<=1e-10) ? "highconf" : "candidate"
        gid = $2; sub(/\.[0-9]+$/, "", gid)
        print q, (q in g ? g[q] : "NA"), (q in s ? s[q] : "NA"), $2, gid, $3, $11, $12, $15, tier
      }' "$MAPFILE" "$MAIN_TSV" > "$SUMMARY"
  else
    warn "未提供 --map 映射表, At 基因号列填 NA(建议传 TRM_representative_proteins.tsv)"
    awk -F'\t' -v OFS='\t' '
      {
        tier = (($15+0)>=70 && ($11+0)<=1e-10) ? "highconf" : "candidate"
        gid = $2; sub(/\.[0-9]+$/, "", gid)
        print $1, "NA", "NA", $2, gid, $3, $11, $12, $15, tier
      }' "$MAIN_TSV" > "$SUMMARY"
  fi
} 
log "A1 候选汇总表: $SUMMARY"

# =============================================================================
step 4 "去重 -> ② 候选 ID(蛋白级 + 基因级)"
# =============================================================================
PROT_TXT="${OUTDIR}/TRM_candidates_protein_level.txt"
GENE_TXT="${OUTDIR}/TRM_candidates_gene_level.txt"

cut -f4 "$SUMMARY" | grep -v '^query$' | sort -u > "$PROT_TXT"
# 基因级: 剥离 isoform 后缀 (.1/.2 ...)
sed 's/\.[0-9]\+$//' "$PROT_TXT" | sort -u > "$GENE_TXT"

N_PROT="$(wc -l < "$PROT_TXT" | tr -d ' ')"
N_GENE="$(wc -l < "$GENE_TXT" | tr -d ' ')"
log "② 候选蛋白(含 isoform): $N_PROT 条 -> $PROT_TXT"
log "② 候选基因(去 isoform): $N_GENE 个 -> $GENE_TXT"

# =============================================================================
step 5 "抽取候选蛋白序列 -> ③ 去重 FASTA"
# =============================================================================
CAND_FA="${OUTDIR}/Zunla_TRM_candidate_proteins.fasta"
CAND_FA_GENE="${OUTDIR}/Zunla_TRM_candidate_proteins_genelevel.fasta"

# 5a: 基因级(列出全部 isoform, 供建树时挑选)
if ! blastdbcmd -db "$DB_PREFIX" -entry_batch "$GENE_TXT" \
      -out "$CAND_FA_GENE" 2>/dev/null || ! grep -q '^>' "$CAND_FA_GENE" 2>/dev/null; then
  warn "blastdbcmd 基因级取序列失败(库可能未加 -parse_seqids), 回退纯 awk 提取"
  awk 'NR==FNR{want[$1]=1; next}
       /^>/ { keep = ($1 in want); if (keep) print; next }
       keep' FS='[ \t]+' "$GENE_TXT" "$DB_FASTA" > "$CAND_FA_GENE"
fi

# 5b: 蛋白级(精确到 isoform)
ok=0
if blastdbcmd -db "$DB_PREFIX" -entry_batch "$PROT_TXT" \
     -out "$CAND_FA" 2>/dev/null && grep -q '^>' "$CAND_FA" 2>/dev/null; then
  ok=1
  n_got="$(grep -c '^>' "$CAND_FA" || true)"
  n_got="${n_got//[^0-9]/}"
  if (( ${n_got:-0} != N_PROT )); then
    warn "蛋白级取到 $n_got 条, 但候选 ID 有 $N_PROT 条, 改用按蛋白 ID 精确提取"
    ok=0
  fi
fi
if (( ok == 0 )); then
  log "回退: 按蛋白 ID(含 isoform 后缀)从原始 FASTA 精确提取"
  awk '
    BEGIN{FS="[ \t]+"}
    NR==FNR{want[$1]=1; next}
    /^>/ { id=$0; sub(/^>/,"",id); sub(/[ \t].*$/,"",id); keep=(id in want); if(keep) print; next }
    keep
  ' "$PROT_TXT" "$DB_FASTA" > "$CAND_FA"
  n_got="$(grep -c '^>' "$CAND_FA" || true)"
  n_got="${n_got//[^0-9]/}"
  if (( ${n_got:-0} == 0 )); then
    warn "蛋白级精确提取仍为 0 条, ③ 退化为基因级结果(含 isoform)"
    cp "$CAND_FA_GENE" "$CAND_FA"
  fi
fi

N_CAND="$(grep -c '^>' "$CAND_FA" || true)";  N_CAND="${N_CAND//[^0-9]/}"
N_CANDG="$(grep -c '^>' "$CAND_FA_GENE" || true)"; N_CANDG="${N_CANDG//[^0-9]/}"
log "③ 候选蛋白 FASTA: $CAND_FA (${N_CAND:-0} 条)"
log "   基因级 FASTA  : $CAND_FA_GENE (${N_CANDG:-0} 条)"

# =============================================================================
step 6 "辅助统计"
# =============================================================================
# A2: 每 query 候选数
awk -F'\t' 'NR>1{c[$2]++} END{for(q in c) printf "%s\t%d\n", q, c[q]}' "$SUMMARY" \
  | sort > "${OUTDIR}/01_candidate_counts.tsv"

# A3: 每 query 最佳命中(按 bitscore)
sort -k1,1 -k12,12gr "$BL_TSV" | awk -F'\t' '!seen[$1]++' \
  > "${OUTDIR}/01_best_hit_per_query.tsv"

# A4: 零命中 query
cut -f1 "$BL_TSV" | sort -u > "${OUTDIR}/.hit_queries.tmp"
grep '^>' "$QUERY" | sed 's/^>//' | awk '{print $1}' | sort -u > "${OUTDIR}/.all_queries.tmp"
comm -23 "${OUTDIR}/.all_queries.tmp" "${OUTDIR}/.hit_queries.tmp" \
  > "${OUTDIR}/01_zero_hit_queries.txt"
rm -f "${OUTDIR}/.hit_queries.tmp" "${OUTDIR}/.all_queries.tmp"

N_ZERO="$(wc -l < "${OUTDIR}/01_zero_hit_queries.txt" | tr -d ' ')"
log "A2 每query候选数 / A3 最佳命中 / A4 零命中($N_ZERO 条)"

# =============================================================================
step 7 "分布统计与出图"
# =============================================================================
PY=""
for c in python3 python; do
  if command -v "$c" >/dev/null 2>&1; then PY="$c"; break; fi
done

if [[ -n "$PY" ]]; then
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  # 统计脚本可能叫这几个名字(随项目命名规范变化), 依次探测取第一个存在的
  SUMMARIZER=""
  for cand in \
      "${SCRIPT_DIR}/04_summarize_candidates.py" \
      "${SCRIPT_DIR}/01_summarize_candidates.py" \
      "${SCRIPT_DIR}/summarize_candidates.py"; do
    if [[ -f "$cand" ]]; then SUMMARIZER="$cand"; break; fi
  done

  if [[ -n "$SUMMARIZER" ]]; then
    log "调用统计脚本: $SUMMARIZER"
    "$PY" "$SUMMARIZER" \
      --blast "$BL_TSV" \
      --query "$QUERY" \
      --outdir "$OUTDIR" \
      --ident "$IDENT" --cov "$COV" --evalue "$EVALUE" || \
      warn "统计脚本返回非零, 但主结果已产出, 请检查上面报错"
  else
    warn "未找到统计脚本(04_summarize_candidates.py / 01_summarize_candidates.py)"
    warn "请确认它与本脚本放在同一目录; 跳过分布统计与出图"
  fi
else
  warn "未找到 python3/python, 跳过分布统计与出图(主结果不受影响)"
fi

# =============================================================================
step 8 "完成"
# =============================================================================
cat <<EOF

=================== 流水线完成 ===================
① 原始 BLAST 结果 : $BL_TSV
② 候选蛋白 ID     : $PROT_TXT  ($N_PROT 条)
   候选基因 ID     : $GENE_TXT  ($N_GENE 个)
③ 候选蛋白 FASTA  : $CAND_FA
辅助文件           : ${OUTDIR}/01_*.tsv|txt|png

核心数字
  候选基因数(去 isoform) : $N_GENE      <-- 你的候选家族规模
  候选蛋白数(含 isoform) : $N_PROT
  高可信 HSP             : $N_HIGH
  片段待查 HSP           : $N_FRAG
  零命中 query           : $N_ZERO

判读建议
  * 候选基因数 <20 -> 偏少, 用 --evalue 1e-3 重跑
  * 候选基因数 >100 -> 偏多, 用 --ident 30 --cov 70 --evalue 1e-10 重跑
  * 优先看 01_distribution.tsv 里 coverage 的分布断崖, 它比 identity 更干净
  * 零命中的 query 记录为该亚家族可能拟南芥特有/辣椒丢失, 本身是结论
==================================================
EOF

}   # run_pipeline() 结束

# 主流程执行: stdout+stderr 同时上屏并落盘到诊断日志。
# set -o pipefail 保证 run_pipeline 里 die 的非零退出能传递给脚本。
run_pipeline 2>&1 | tee "$DIAG"
exit "${PIPESTATUS[0]}"
