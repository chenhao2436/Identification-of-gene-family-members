#!/usr/bin/env bash
# =============================================================================
# stage_02_msa.sh —— 阶段 02：MAFFT 多序列比对（34 条 AtTRM）
#
# 方法（标准参数，不调优、不裁剪）
#   mafft --auto --anysymbol --thread N AtTRM.fasta > AtTRM.aln.fasta
#   --auto      : 按序列数与相似度自动选 FFT-NS-2 / L-INS-i（MAFFT 官方推荐默认）
#   --anysymbol : 允许非标准残基（B/Z/X/U）而不报错
#
# 用法:
#   bash stage_02_msa.sh [--threads N]
#
# 产出（06_hmmer_results/02_msa/）
#   02_msa.log                 完整日志
#   AtTRM_34.aln.fasta         比对结果（构建 HMM 的输入）
#   AtTRM_34.aln.stats.txt     esl-alistat 统计（列数/序列数/平均一致性等）
#   02_msa.subset.tsv          机读小结
#   02_msa_report.txt          中文结论
# =============================================================================
set -euo pipefail

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPTS_DIR}/common.sh"

parse_ref_args "$@" || { sed -n '2,24p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0; }

OUTDIR="$D02"
mkdir -p "$OUTDIR"
LOG="${OUTDIR}/02_msa.log"
ALN="${OUTDIR}/AtTRM_34.aln.fasta"
STATS="${OUTDIR}/AtTRM_34.aln.stats.txt"
SUMMARY="${OUTDIR}/02_msa.subset.tsv"
REPORT="${OUTDIR}/02_msa_report.txt"

subset_init "$SUMMARY"

# 日志同时上屏与落盘：用「函数 + 管道 tee」而非 exec >(tee ...)，
# 后者是后台进程，脚本退出时可能截断末尾输出。
run_stage() {

step "阶段 02：MAFFT 多序列比对"
log "输入 query : $AT_QUERY_FA"
log "输出比对   : $ALN"
log "线程数     : $THREADS"

# ------------------------------ 前置检查 ------------------------------------
have mafft || die "未找到 mafft。请先安装/加载（阶段 01 会报告）。"
[[ -f "$AT_QUERY_FA" ]] || die "找不到 AtTRM query FASTA: $AT_QUERY_FA"
is_text_fasta "$AT_QUERY_FA" || die "AtTRM query 不是纯文本 FASTA: $AT_QUERY_FA"

N_IN="$(count_fasta "$AT_QUERY_FA")"
log "query 序列数: $N_IN (期望 $EXPECT_AT_N)"
if (( N_IN != EXPECT_AT_N )); then
  warn "query 不是 $EXPECT_AT_N 条（实际 $N_IN）。继续执行，但请确认用的是 34 条代表蛋白。"
fi
printf '输入\tAtTRM query\t%d\t%d\t%s\t\n' "$EXPECT_AT_N" "$N_IN" \
  "$( (( N_IN == EXPECT_AT_N )) && echo OK || echo 'WARN:条数不符' )" >> "$SUMMARY"

# ------------------------------ 运行 MAFFT ----------------------------------
step "运行 mafft --auto"
MAFFT_VER="$(mafft --version 2>&1 | head -1 || echo unknown)"
log "mafft 版本: $MAFFT_VER"

MAFFT_CMD="mafft --auto --anysymbol --thread $THREADS"
log "命令: $MAFFT_CMD '$AT_QUERY_FA' > '$ALN'"

# shellcheck disable=SC2086
mafft --auto --anysymbol --thread "$THREADS" "$AT_QUERY_FA" > "$ALN"

[[ -s "$ALN" ]] || die "MAFFT 输出为空: $ALN"
N_OUT="$(count_fasta "$ALN")"
log "比对完成: $N_OUT 条序列"
(( N_OUT == N_IN )) || die "比对前后序列数不一致 ($N_IN -> $N_OUT)"
printf '比对\t输出序列数\t%d\t%d\tOK\t\n' "$N_IN" "$N_OUT" >> "$SUMMARY"

# ------------------------------ 比对统计 ------------------------------------
step "比对统计"
ALN_COLS="$(awk '/^>/{next}{n+=length($0)}END{print n}' "$ALN")"
log "未 gapped 残基总数: $ALN_COLS"
# 比对列数（取第一条序列的 gapped 长度）
ALN_LEN="$(awk '/^>/{if(seq!=""){print length(seq); exit} next}{seq=seq $0}' "$ALN")"
log "比对列数(含 gap): ${ALN_LEN:-未知}"

if have esl-alistat; then
  esl-alistat "$ALN" > "$STATS" 2>&1 || true
  log "esl-alistat 统计已写入: $STATS"
  # 抽取关键行便于阅读
  grep -iE 'Number of sequences|Alignment length|Average identity|Most related|Least related|Total residues' \
    "$STATS" | sed 's/^/    /' || true
  printf '比对\tesl-alistat\t运行成功\t成功\tOK\t\n' >> "$SUMMARY"
else
  warn "未找到 esl-alistat（HMMER 套件自带），跳过详细统计"
  printf '比对\tesl-alistat\t可用\t缺失\tWARN:跳过详细统计\t\n' >> "$SUMMARY"
fi
printf '比对\t比对列数\t-\t%s\tINFO\t\n' "${ALN_LEN:-NA}" >> "$SUMMARY"

# ------------------------------ 报告 ----------------------------------------
{
  echo "阶段 02：MAFFT 多序列比对报告"
  echo "======================================================================"
  echo "时间        : $(date '+%F %T')"
  echo "软件        : mafft $MAFFT_VER"
  echo "参数        : --auto --anysymbol --thread $THREADS"
  echo "输入        : $AT_QUERY_FA"
  echo "输出        : $ALN"
  echo "输入序列数  : $N_IN"
  echo "输出序列数  : $N_OUT"
  echo "比对列数    : ${ALN_LEN:-NA}"
  echo ""
  echo "说明: --auto 让 MAFFT 按序列数与相似度自动选择算法(FFT-NS-2 或 L-INS-i)，"
  echo "      这是 MAFFT 官方推荐的通用默认，未针对本家族调优。"
  echo "      未做 gap 裁剪(trimAl/Gblocks)：HMMER 建库本身对 gap 有加权，"
  echo "      裁剪会丢弃信息；若你希望裁剪，请告知，属方法变更。"
  echo ""
  echo "下一步: bash run.sh 03   （hmmbuild 建立 TRM.hmm）"
} > "$REPORT"

echo ""
log "产出: $ALN"
log "      $STATS"
log "      $REPORT"

}   # run_stage 结束

run_stage 2>&1 | tee "$LOG"
exit "${PIPESTATUS[0]}"
