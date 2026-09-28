#!/usr/bin/env bash
# =============================================================================
# stage_05_merge.sh —— 阶段 05：BLAST 与 HMM 候选合并去重（按 gene）
#
# 方法
#   并集合并，不删除任何候选；按 gene（剥 isoform）去重，标注四类证据：
#     BLAST+HMM / BLAST-only / HMM-only / HMM-inclusive-only
#
# 用法:
#   bash stage_05_merge.sh
#
# 产出（06_hmmer_results/05_merge/）
#   05_merge.log
#   CaTRM_candidates.merged.tsv            每基因一行（主结果）
#   CaTRM_candidates.evidence_counts.tsv   四类证据计数
#   CaTRM_candidates.by_evidence.tsv       按证据分组
#   CaTRM_candidates.stats.txt             统计与判读要点
#   05_merge.subset.tsv                    机读小结
#   05_merge_report.txt                    中文结论
# =============================================================================
set -euo pipefail

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPTS_DIR}/common.sh"

parse_ref_args "$@" || { sed -n '2,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0; }

OUTDIR="$D05"
mkdir -p "$OUTDIR"
LOG="${OUTDIR}/05_merge.log"
PREFIX="${OUTDIR}/CaTRM_candidates"
SUMMARY="${OUTDIR}/05_merge.subset.tsv"
REPORT="${OUTDIR}/05_merge_report.txt"
MERGER="${SCRIPTS_DIR}/merge_blast_hmm.py"

# HMM 侧输入（阶段 04 产出）
HMM_GENE_LEVEL="$D04/TRM.hmmsearch.gene_level.tsv"
HMM_GENE_STRICT="$D04/TRM.hmmsearch.gene_ids.txt"
HMM_GENE_INCL="$D04/TRM.hmmsearch.gene_ids_inclusive.txt"

subset_init "$SUMMARY"

run_stage() {

step "阶段 05：合并 BLAST 与 HMM 候选"
log "BLAST 候选 : $BLAST_GENE_LIST"
log "HMM 候选   : $HMM_GENE_STRICT"
log "输出前缀   : $PREFIX"

# ------------------------------ 前置检查 ------------------------------------
PY="$(find_python || true)"
[[ -n "$PY" ]] || die "需要 python3，但未找到"
[[ -f "$MERGER" ]] || die "找不到合并脚本: $MERGER"

if [[ ! -f "$BLAST_GENE_LIST" ]]; then
  warn "找不到 BLAST 候选列表: $BLAST_GENE_LIST"
  warn "  将只输出 HMM 候选（BLAST 侧记 NA）。若 BLAST 结果在别处，设 BLAST_DIR=<目录> 重跑。"
  printf '输入\tBLAST候选列表\t存在\t未找到\tWARN:仅输出HMM候选\t\n' >> "$SUMMARY"
else
  n_b="$(count_lines "$BLAST_GENE_LIST")"
  log "BLAST 候选基因数: $n_b"
  printf '输入\tBLAST候选基因\t-\t%d\tOK\t\n' "$n_b" >> "$SUMMARY"
fi

if [[ ! -f "$HMM_GENE_STRICT" && ! -f "$HMM_GENE_LEVEL" ]]; then
  die "找不到阶段 04 的 HMM 结果。
    请先运行: bash run.sh 04"
fi
if [[ -f "$HMM_GENE_STRICT" ]]; then
  n_h="$(count_lines "$HMM_GENE_STRICT")"
  log "HMM strict 候选基因数: $n_h"
  printf '输入\tHMM strict候选基因\t-\t%d\tOK\t\n' "$n_h" >> "$SUMMARY"
fi

# ------------------------------ 运行合并 ------------------------------------
step "运行合并脚本"
"$PY" "$MERGER" \
  --blast-gene  "$BLAST_GENE_LIST" \
  --blast-prots "$BLAST_PROT_LIST" \
  --blast-main  "$BLAST_MAIN" \
  --hmm-gene    "$HMM_GENE_LEVEL" \
  --hmm-gene-strict "$HMM_GENE_STRICT" \
  --hmm-gene-incl   "$HMM_GENE_INCL" \
  --anno        "$AT_ANNO_TSV" \
  --prefix      "$PREFIX"

MERGED="${PREFIX}.merged.tsv"
[[ -f "$MERGED" ]] || die "合并未产出: $MERGED"

N_UNION="$(awk 'NR>1' "$MERGED" | wc -l | tr -d ' ')"
awk -F'\t' 'NR>1{c[$2]++} END{for(k in c) printf "%s\t%d\n", k, c[k]}' "$MERGED" \
  | sort > "${OUTDIR}/.ev.tmp" || true
while IFS=$'\t' read -r ev n; do
  printf '合并\t%s\t-\t%d\tINFO\t\n' "$ev" "$n" >> "$SUMMARY"
done < "${OUTDIR}/.ev.tmp"
rm -f "${OUTDIR}/.ev.tmp"
printf '合并\t候选基因总数(并集)\t-\t%d\tINFO\t\n' "$N_UNION" >> "$SUMMARY"

# ------------------------------ 报告 ----------------------------------------
{
  echo "阶段 05：BLAST + HMM 合并候选 —— 报告"
  echo "======================================================================"
  echo "时间        : $(date '+%F %T')"
  echo "主统计单位  : gene（isoform 版本号已剥离）"
  echo "合并策略    : 并集，不删除任何候选"
  echo ""
  cat "${PREFIX}.stats.txt"
  echo ""
  echo "【列说明（merged.tsv）】"
  echo "  gene                   Zunla 基因 ID（主键）"
  echo "  evidence               证据标签（四类之一）"
  echo "  n_blast_proteins       BLAST 命中的该基因蛋白数"
  echo "  blast_best_pident      最佳 identity(%)"
  echo "  blast_best_qcov        最佳 query coverage(%)"
  echo "  blast_min_evalue       最小 E 值"
  echo "  blast_best_bitscore    最大 bitscore"
  echo "  n_hmm_proteins         HMM 命中的该基因蛋白数"
  echo "  hmm_best_full_evalue   HMM 最佳全序列 E 值"
  echo "  hmm_best_dom_score     HMM 最佳 domain score"
  echo "  hmm_best_hmm_cov       命中的 HMM 模型覆盖度(%)"
  echo "  hmm_significance       strict / inclusive_only"
  echo "  n_at_queries           支持的 AtTRM query 数"
  echo "  at_genes / at_symbols  对应的 At 基因与 TRM 符号"
  echo "  blast_proteins         BLAST 命中的蛋白 ID 列表"
  echo "  hmm_best_protein       HMM 最佳蛋白"
  echo "  manual_review          待人工确认标记（不删除，只标记）"
  echo ""
  echo "下一步: bash run.sh 06   （Pfam 结构域核查）"
} > "$REPORT"

log "产出: $MERGED"
log "      ${PREFIX}.evidence_counts.tsv"
log "      ${PREFIX}.stats.txt"
log "      $REPORT"

}   # run_stage 结束

run_stage 2>&1 | tee "$LOG"
exit "${PIPESTATUS[0]}"
