#!/usr/bin/env bash
# =============================================================================
# stage_07_final.sh —— 阶段 07：汇总最终 CaTRM 基因/蛋白列表
#
# 汇总 BLAST + HMM + 结构域 + Zunla 注释，输出最终列表，
# 并为每个成员给出证据来源与筛选理由；可疑者单列"待人工确认"，不删除。
#
# 分层口径
#   conservative  双方法支持 + 含 PF14309 核心域
#   standard      双方法支持，或 单方法但 HMM 严格命中
#   review        仅达包含阈值 / 仅单方法无结构域支持 / 覆盖度偏低 等
#
# 用法:
#   bash stage_07_final.sh [--gff3 /path/Canz.genome.gff3]
#
# 产出（06_hmmer_results/07_final/）
#   07_final.log
#   CaTRM_final.gene_list.tsv     最终基因列表（主交付）
#   CaTRM_final.gene_ids.txt      仅基因 ID
#   CaTRM_final.protein_list.tsv  蛋白级列表
#   CaTRM_final.proteins.fasta    每基因代表蛋白
#   CaTRM_final.evidence.tsv      每成员证据与理由
#   CaTRM_final.review.tsv        待人工确认清单
#   CaTRM_final.stats.txt         分层统计
#   CaTRM_Methods.txt             论文 Methods 素材（软件/版本/参数/各步候选数）
#   07_final.subset.tsv / 07_final_report.txt
# =============================================================================
set -euo pipefail

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPTS_DIR}/common.sh"

parse_ref_args "$@" || { sed -n '2,24p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0; }

OUTDIR="$D07"
mkdir -p "$OUTDIR"
LOG="${OUTDIR}/07_final.log"
PREFIX="${OUTDIR}/CaTRM_final"
SUMMARY="${OUTDIR}/07_final.subset.tsv"
REPORT="${OUTDIR}/07_final_report.txt"
METHODS="${OUTDIR}/CaTRM_Methods.txt"
SUMMARIZER="${SCRIPTS_DIR}/final_summary.py"

MERGED="$D05/CaTRM_candidates.merged.tsv"
PER_PROT="$D06/CaTRM_vs_Pfam.per_protein.tsv"
FOCUS_MX="$D06/CaTRM_vs_Pfam.focus_matrix.tsv"
HMM_GENE_LEVEL="$D04/TRM.hmmsearch.gene_level.tsv"

subset_init "$SUMMARY"

run_stage() {

step "阶段 07：汇总最终 CaTRM 列表"
log "合并表   : $MERGED"
log "结构域表 : $PER_PROT"

# ------------------------------ 前置检查 ------------------------------------
PY="$(find_python || true)"
[[ -n "$PY" ]] || die "需要 python3，但未找到"
[[ -f "$SUMMARIZER" ]] || die "找不到汇总脚本: $SUMMARIZER"
[[ -f "$MERGED" ]] || die "找不到阶段 05 的合并结果: $MERGED
    请先运行: bash run.sh 05"

if [[ ! -f "$PER_PROT" ]]; then
  warn "未找到阶段 06 的结构域结果: $PER_PROT"
  warn "  has_PF14309 列将填 NA，分层只能依据 BLAST+HMM 证据。"
  warn "  如需完整分层，请先补齐阶段 06（需要 Pfam-A.hmm）。"
  printf '输入\t结构域结果\t存在\t未找到\tWARN:分层缺少结构域证据\t\n' >> "$SUMMARY"
else
  n_dom="$(awk 'NR>1' "$PER_PROT" | wc -l | tr -d ' ')"
  printf '输入\t结构域结果\t-\t%d蛋白\tOK\t\n' "$n_dom" >> "$SUMMARY"
fi

GFF3_USE="$(probe_gff3 || true)"
if [[ -n "$GFF3_USE" ]]; then
  log "GFF3     : $GFF3_USE（用于补 gene 坐标）"
  printf '输入\tGFF3坐标\t存在\t%s\tOK\t\n' "$(basename "$GFF3_USE")" >> "$SUMMARY"
else
  warn "未找到 GFF3，坐标列将填 NA（用 --gff3 指定）"
  printf '输入\tGFF3坐标\t存在\t未找到\tWARN:坐标为NA\t\n' >> "$SUMMARY"
fi

# ------------------------------ 运行汇总 ------------------------------------
step "运行汇总脚本"
"$PY" "$SUMMARIZER" \
  --merged       "$MERGED" \
  --per-protein  "$PER_PROT" \
  --focus-matrix "$FOCUS_MX" \
  --gene-level   "$HMM_GENE_LEVEL" \
  --gff3         "$GFF3_USE" \
  --pep          "$ZUNLA_PEP" \
  --prefix       "$PREFIX"

GENE_LIST="${PREFIX}.gene_list.tsv"
[[ -f "$GENE_LIST" ]] || die "汇总未产出基因列表: $GENE_LIST"

N_FINAL="$(awk 'NR>1' "$GENE_LIST" | wc -l | tr -d ' ')"
N_CONS="$(awk -F'\t' 'NR>1 && $7=="conservative"' "$GENE_LIST" | wc -l | tr -d ' ')"
N_STD="$(awk -F'\t' 'NR>1 && $7=="standard"' "$GENE_LIST" | wc -l | tr -d ' ')"
N_REV="$(awk -F'\t' 'NR>1 && $7=="review"' "$GENE_LIST" | wc -l | tr -d ' ')"
printf '最终\t候选基因总数\t-\t%d\tINFO\t\n' "$N_FINAL" >> "$SUMMARY"
printf '最终\tconservative\t-\t%d\tINFO\t\n' "$N_CONS" >> "$SUMMARY"
printf '最终\tstandard\t-\t%d\tINFO\t\n' "$N_STD" >> "$SUMMARY"
printf '最终\treview(待人工确认)\t-\t%d\tINFO\t\n' "$N_REV" >> "$SUMMARY"

# -------------------------- Methods 素材（可写进论文） -----------------------
step "生成 Methods 素材"
{
  echo "CaTRM 基因家族鉴定 —— Methods 素材"
  echo "======================================================================"
  echo "生成时间: $(date '+%F %T')"
  echo ""
  echo "1. 数据来源"
  echo "   1.1 拟南芥 TRM 家族: 以 Araport11 蛋白组为来源，取 34 个 AtTRM 基因"
  echo "       各一条最长转录本作代表蛋白（$AT_QUERY_FA）。"
  echo "   1.2 辣椒 Zunla 全蛋白: $ZUNLA_PEP"
  echo "       序列数: $(count_fasta "$ZUNLA_PEP")"
  echo "   1.3 基因组与注释: ${GFF3_USE:-未指定}"
  echo ""
  echo "2. 软件与版本"
  echo "   MAFFT     : $( [[ -f "$D02/02_msa_report.txt" ]] && grep -m1 '^软件' "$D02/02_msa_report.txt" | sed 's/软件 *: *//' || echo '见 02_msa_report.txt' )"
  echo "   HMMER     : $( [[ -f "$D03/03_hmmbuild_report.txt" ]] && grep -m1 '^软件' "$D03/03_hmmbuild_report.txt" | sed 's/软件 *: *//' || echo '见 03_hmmbuild_report.txt' )"
  echo "   Pfam      : $( [[ -f "$D06/06_domain_report.txt" ]] && grep -m1 '^Pfam 数据库' "$D06/06_domain_report.txt" | sed 's/Pfam 数据库 *: *//' || echo '见 06_domain_report.txt' )"
  echo ""
  echo "3. 方法与参数"
  echo "   3.1 多序列比对"
  echo "       mafft --auto --anysymbol --thread <N> AtTRM_34.fasta > AtTRM_34.aln.fasta"
  echo "   3.2 家族 HMM 构建"
  echo "       hmmbuild --amino --cpu <N> TRM.hmm AtTRM_34.aln.fasta"
  echo "   3.3 全蛋白组扫描"
  echo "       hmmsearch --domtblout TRM_vs_Zunla.domtblout -E ${HMM_EVAL} --cpu <N> \\"
  echo "                 TRM.hmm Zunla_Canz.pep.fa"
  echo "       判定口径: 全序列 E-value <= ${HMM_EVAL} 记为 strict;"
  echo "                 另按包含阈值 i-Evalue <= 1e-2 记录 inclusive 档（最宽松）"
  echo "   3.4 BLASTP 初筛（本项目已有）"
  echo "       $( [[ -f "$BLAST_SUMMARY" ]] && grep -m1 '阈值设定' "$BLAST_SUMMARY" | sed 's/^ *//' || echo '见 05_blast_results/01_summary_stats.txt' )"
  echo "   3.5 结构域核查"
  echo "       hmmscan --domtblout CaTRM_vs_Pfam.domtblout --cpu <N> Pfam-A.hmm \\"
  echo "               CaTRM_candidates_proteins.fasta"
  echo "       显著性阈值 i-Evalue <= ${DOM_EVAL}"
  echo "       重点域: PF14309 (LONGIFOLIA1/2-like C-terminal domain, DUF4378)"
  echo "               PF14383 (DUF761-associated sequence motif; 与 TRM 关联未经证实)"
  echo "   3.6 候选合并"
  echo "       取 BLAST 与 HMM 候选的并集，按 gene（剥离 isoform 版本号）去重，"
  echo "       标注 BLAST+HMM / BLAST-only / HMM-only / HMM-inclusive-only 四类证据。"
  echo ""
  echo "4. 各步候选数量"
  printf "   %-38s %s\n" "AtTRM 输入代表蛋白" "$(count_fasta "$AT_QUERY_FA")"
  printf "   %-38s %s\n" "BLAST 候选基因" "$( [[ -f "$BLAST_GENE_LIST" ]] && count_lines "$BLAST_GENE_LIST" || echo NA )"
  printf "   %-38s %s\n" "HMM strict 候选基因" "$( [[ -f "$D04/TRM.hmmsearch.gene_ids.txt" ]] && count_lines "$D04/TRM.hmmsearch.gene_ids.txt" || echo NA )"
  printf "   %-38s %s\n" "HMM inclusive 候选基因" "$( [[ -f "$D04/TRM.hmmsearch.gene_ids_inclusive.txt" ]] && count_lines "$D04/TRM.hmmsearch.gene_ids_inclusive.txt" || echo NA )"
  printf "   %-38s %s\n" "合并去重后候选基因" "$N_FINAL"
  printf "   %-38s %s\n" "  其中 conservative" "$N_CONS"
  printf "   %-38s %s\n" "  其中 standard" "$N_STD"
  printf "   %-38s %s\n" "  其中 review(待人工确认)" "$N_REV"
  printf "   %-38s %s\n" "最终代表蛋白序列" "$(count_fasta "${PREFIX}.proteins.fasta" 2>/dev/null || echo NA)"
  echo ""
  echo "5. 分层标准"
  echo "   conservative: BLAST 与 HMM 双重支持 且 检出 PF14309 核心域"
  echo "   standard    : BLAST 与 HMM 双重支持; 或 单方法支持但 HMM 达严格阈值"
  echo "   review      : 仅达 HMM 包含阈值 / 仅单方法且无结构域支持 / 覆盖度偏低"
  echo "   说明: 所有候选均保留，review 层单列于 review.tsv，未做任何删除。"
} > "$METHODS"

printf '产出\tMethods 素材\t-\t%s\tOK\t\n' "$(basename "$METHODS")" >> "$SUMMARY"

# ------------------------------ 报告 ----------------------------------------
{
  echo "阶段 07：最终 CaTRM 列表 —— 报告"
  echo "======================================================================"
  echo "时间: $(date '+%F %T')"
  echo ""
  cat "${PREFIX}.stats.txt"
  echo ""
  echo "【最终基因列表（前 40 行）】"
  head -41 "$GENE_LIST" | column -t 2>/dev/null | cut -c1-200 | sed 's/^/  /' \
    || head -41 "$GENE_LIST" | sed 's/^/  /'
  echo ""
  echo "【待人工确认（review 层）】"
  if [[ -s "${PREFIX}.review.tsv" ]]; then
    cat "${PREFIX}.review.tsv" | sed 's/^/  /'
  else
    echo "  (无)"
  fi
  echo ""
  echo "【论文 Methods 素材】"
  echo "  $METHODS"
  echo ""
  echo "流程全部完成。主要交付:"
  echo "  $(basename "$GENE_LIST")       最终 CaTRM 基因列表"
  echo "  $(basename "${PREFIX}.proteins.fasta")         每基因代表蛋白序列"
  echo "  $(basename "${PREFIX}.evidence.tsv")           每成员证据与筛选理由"
} > "$REPORT"

cat "$METHODS"
echo ""
log "产出: $GENE_LIST"
log "      ${PREFIX}.proteins.fasta"
log "      $METHODS"
log "      $REPORT"

}   # run_stage 结束

run_stage 2>&1 | tee "$LOG"
exit "${PIPESTATUS[0]}"
