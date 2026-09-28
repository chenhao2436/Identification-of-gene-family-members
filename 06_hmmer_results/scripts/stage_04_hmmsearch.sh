#!/usr/bin/env bash
# =============================================================================
# stage_04_hmmsearch.sh —— 阶段 04：hmmsearch 扫描 Zunla 全蛋白
#
# 方法（标准参数，保留完整原始结果）
#   hmmsearch --domtblout <domtbl> -E 1e-5 --cpu N TRM.hmm Zunla_Canz.pep.fa > <stdout>
#     -E 1e-5 : 只影响「报告哪些命中」，不改变打分；
#               真正的判定依据是 domtblout 里的 E-value / i-Evalue，两档都会统计。
#     --domtblout : 保留逐 domain 的完整原始结果（含 i-Evalue、坐标）
#
# 用法:
#   bash stage_04_hmmsearch.sh [--threads N] [--hmm-eval 1e-5]
#
# 产出（06_hmmer_results/04_hmmsearch/）
#   04_hmmsearch.log                完整日志
#   TRM_vs_Zunla.hmmsearch.txt      hmmsearch 原始 stdout（完整保留）
#   TRM_vs_Zunla.domtblout          逐 domain 原始结果（完整保留）
#   TRM_vs_Zunla.tblout             逐序列汇总原始结果
#   TRM.hmmsearch.*                 解析后的多张表（见 parse_hmmsearch_domtbl.py）
#   04_hmmsearch.subset.tsv         机读小结
#   04_hmmsearch_report.txt         中文结论
# =============================================================================
set -euo pipefail

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPTS_DIR}/common.sh"

# 支持 --hmm-eval 覆盖初筛阈值
HMM_EVAL_ARG=""
PASS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --hmm-eval) HMM_EVAL="$2"; shift 2 ;;
    *) PASS+=("$1"); shift ;;
  esac
done
parse_ref_args "${PASS[@]+"${PASS[@]}"}" || { sed -n '2,24p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0; }

OUTDIR="$D04"
mkdir -p "$OUTDIR"
LOG="${OUTDIR}/04_hmmsearch.log"
HMM="$D03/TRM.hmm"
STDOUT_F="${OUTDIR}/TRM_vs_Zunla.hmmsearch.txt"
DOMTBL="${OUTDIR}/TRM_vs_Zunla.domtblout"
TBLOUT="${OUTDIR}/TRM_vs_Zunla.tblout"
PREFIX="${OUTDIR}/TRM.hmmsearch"
SUMMARY="${OUTDIR}/04_hmmsearch.subset.tsv"
REPORT="${OUTDIR}/04_hmmsearch_report.txt"
PARSER="${SCRIPTS_DIR}/parse_hmmsearch_domtbl.py"

subset_init "$SUMMARY"

run_stage() {

step "阶段 04：hmmsearch 扫描 Zunla 全蛋白"
log "HMM 输入   : $HMM"
log "蛋白库     : $ZUNLA_PEP"
log "初筛 E 值  : $HMM_EVAL"
log "线程数     : $THREADS"

# ------------------------------ 前置检查 ------------------------------------
have hmmsearch || die "未找到 hmmsearch。请先安装/加载 HMMER（阶段 01 会报告）。"
[[ -f "$HMM" ]] || die "找不到阶段 03 的 HMM: $HMM
    请先运行: bash run.sh 03"
[[ -f "$ZUNLA_PEP" ]] || die "找不到 Zunla 蛋白库: $ZUNLA_PEP"

N_DB="$(count_fasta "$ZUNLA_PEP")"
log "蛋白库序列数: $N_DB"
if (( N_DB < 10000 )); then
  warn "蛋白库只有 $N_DB 条，明显少于完整蛋白组。"
  warn "候选数会被系统性低估——请先确认这份 FASTA 是否完整。"
fi
printf '输入\tZunla蛋白库\t完整蛋白组\t%d条\t%s\t\n' "$N_DB" \
  "$( (( N_DB < 10000 )) && echo 'WARN:条数偏少' || echo OK )" >> "$SUMMARY"
printf '输入\tTRM.hmm\t存在\t%s\tOK\t\n' "$(basename "$HMM")" >> "$SUMMARY"

# ------------------------------ 运行 hmmsearch ------------------------------
step "运行 hmmsearch"
HMMER_VER="$(hmmsearch -h 2>&1 | grep -oE '# HMMER [0-9.]+' | head -1 || echo 'unknown')"
log "HMMER 版本: $HMMER_VER"

CMD=(hmmsearch --domtblout "$DOMTBL" --tblout "$TBLOUT" -E "$HMM_EVAL" --cpu "$THREADS" "$HMM" "$ZUNLA_PEP")
log "命令: ${CMD[*]}"
"${CMD[@]}" > "$STDOUT_F" 2>&1

[[ -s "$DOMTBL" ]] || die "hmmsearch 未产出 domtblout: $DOMTBL（查看 $STDOUT_F）"
log "hmmsearch 完成"
log "  stdout       : $STDOUT_F ($(wc -l < "$STDOUT_F" | tr -d ' ') 行)"
log "  domtblout    : $DOMTBL ($(grep -vc '^#' "$DOMTBL" || true) 条 domain 命中)"
log "  tblout       : $TBLOUT ($(grep -vc '^#' "$TBLOUT" || true) 条序列命中)"
printf 'hmmsearch\t原始结果\t保留\tdomtblout+stdout+tblout\tOK\t\n' >> "$SUMMARY"

# ------------------------------ 解析 ----------------------------------------
step "解析原始结果（纯标准库）"
PY="$(find_python || true)"
[[ -n "$PY" ]] || die "需要 python3 来解析 domtblout，但未找到"
[[ -f "$PARSER" ]] || die "找不到解析脚本: $PARSER"

"$PY" "$PARSER" \
  --domtbl "$DOMTBL" \
  --pep "$ZUNLA_PEP" \
  --prefix "$PREFIX" \
  --hmm-eval "$HMM_EVAL" \
  --incl-eval 1e-2

GL="${PREFIX}.gene_level.tsv"
[[ -f "$GL" ]] || die "解析未产出基因级汇总表: $GL"

N_GENE_STRICT="$( [[ -f "${PREFIX}.gene_ids.txt" ]] && count_lines "${PREFIX}.gene_ids.txt" || echo 0 )"
N_PROT_STRICT="$( [[ -f "${PREFIX}.protein_ids.txt" ]] && count_lines "${PREFIX}.protein_ids.txt" || echo 0 )"
N_GENE_INCL="$( [[ -f "${PREFIX}.gene_ids_inclusive.txt" ]] && count_lines "${PREFIX}.gene_ids_inclusive.txt" || echo 0 )"
N_PROT_INCL="$( [[ -f "${PREFIX}.protein_ids_inclusive.txt" ]] && count_lines "${PREFIX}.protein_ids_inclusive.txt" || echo 0 )"

printf 'HMM结果\tstrict(E<=%s) 蛋白\t-\t%d\tINFO\t\n' "$HMM_EVAL" "$N_PROT_STRICT" >> "$SUMMARY"
printf 'HMM结果\tstrict(E<=%s) 基因\t-\t%d\tINFO\t\n' "$HMM_EVAL" "$N_GENE_STRICT" >> "$SUMMARY"
printf 'HMM结果\tinclusive(i-E<=1e-2) 蛋白\t-\t%d\tINFO\t\n' "$N_PROT_INCL" >> "$SUMMARY"
printf 'HMM结果\tinclusive(i-E<=1e-2) 基因\t-\t%d\tINFO\t\n' "$N_GENE_INCL" >> "$SUMMARY"

# ------------------------------ 报告 ----------------------------------------
{
  echo "阶段 04：hmmsearch 扫描 Zunla 全蛋白 —— 报告"
  echo "======================================================================"
  echo "时间          : $(date '+%F %T')"
  echo "软件          : $HMMER_VER"
  echo "参数          : --domtblout --tblout -E $HMM_EVAL --cpu $THREADS"
  echo "HMM           : $HMM"
  echo "蛋白库        : $ZUNLA_PEP ($N_DB 条)"
  echo ""
  echo "【原始结果（完整保留）】"
  echo "  $(basename "$STDOUT_F")     hmmsearch stdout"
  echo "  $(basename "$DOMTBL")       逐 domain（含 i-Evalue / 坐标）"
  echo "  $(basename "$TBLOUT")       逐序列汇总"
  echo ""
  echo "【候选数量（两档口径，均按固定阈值统计，未凑数）】"
  echo "  strict    全序列 E <= $HMM_EVAL : $N_PROT_STRICT 蛋白 / $N_GENE_STRICT 基因"
  echo "  inclusive i-Evalue <= 1e-2      : $N_PROT_INCL 蛋白 / $N_GENE_INCL 基因"
  echo ""
  echo "【解析产物】"
  echo "  $(basename "$PREFIX").parsed.tsv            全部命中明细"
  echo "  $(basename "$PREFIX").best_per_protein.tsv  每蛋白最佳"
  echo "  $(basename "$GL")   基因级汇总(主统计口径)"
  echo "  $(basename "$PREFIX").hit_proteins.fasta    命中蛋白序列"
  echo ""
  echo "说明: -E 只决定报告范围，不改变打分；判定依据全部取自原始 domtblout 的"
  echo "      E-value 与 i-Evalue。未为了获得某个预设候选数量而调整阈值。"
  echo ""
  echo "下一步: bash run.sh 05   （与 BLAST 候选合并去重）"
} > "$REPORT"

cat "${PREFIX}.stats.txt"

log "产出目录: ${OUTDIR#"$PROJECT_ROOT"/}"

}   # run_stage 结束

run_stage 2>&1 | tee "$LOG"
exit "${PIPESTATUS[0]}"
