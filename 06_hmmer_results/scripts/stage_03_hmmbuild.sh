#!/usr/bin/env bash
# =============================================================================
# stage_03_hmmbuild.sh —— 阶段 03：hmmbuild 建立 TRM.hmm
#
# 方法（标准参数）
#   hmmbuild --amino --cpu N TRM.hmm AtTRM_34.aln.fasta
#   --amino : 明确按蛋白处理（默认会按序列自动判定，显式指定更可复现）
#   --cpu   : 使用指定线程（hmmbuild 单序列建库，线程主要影响内部并行）
#
# 用法:
#   bash stage_03_hmmbuild.sh [--threads N]
#
# 产出（06_hmmer_results/03_hmmbuild/）
#   03_hmmbuild.log       完整日志
#   TRM.hmm               用于 hmmsearch 的 profile HMM
#   TRM.hmm.header.txt    HMM 头部统计（序列数/比对长度/有效序列数）
#   03_hmmbuild.subset.tsv 机读小结
#   03_hmmbuild_report.txt 中文结论
# =============================================================================
set -euo pipefail

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPTS_DIR}/common.sh"

parse_ref_args "$@" || { sed -n '2,22p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0; }

OUTDIR="$D03"
mkdir -p "$OUTDIR"
LOG="${OUTDIR}/03_hmmbuild.log"
ALN="$D02/AtTRM_34.aln.fasta"
HMM="${OUTDIR}/TRM.hmm"
HDR="${OUTDIR}/TRM.hmm.header.txt"
SUMMARY="${OUTDIR}/03_hmmbuild.subset.tsv"
REPORT="${OUTDIR}/03_hmmbuild_report.txt"

subset_init "$SUMMARY"

run_stage() {

step "阶段 03：hmmbuild 建立 TRM.hmm"
log "输入比对 : $ALN"
log "输出 HMM : $HMM"

# ------------------------------ 前置检查 ------------------------------------
have hmmbuild || die "未找到 hmmbuild。请先安装/加载 HMMER（阶段 01 会报告）。"
[[ -f "$ALN" ]] || die "找不到阶段 02 的比对结果: $ALN
    请先运行: bash run.sh 02"
is_text_fasta "$ALN" || die "比对结果不是纯文本 FASTA: $ALN"

N_SEQ="$(count_fasta "$ALN")"
log "比对序列数: $N_SEQ"
(( N_SEQ >= 4 )) || die "比对序列太少($N_SEQ)，无法建立可靠的 profile HMM"
printf '输入\t比对序列数\t%d\t%d\tOK\t\n' "$EXPECT_AT_N" "$N_SEQ" >> "$SUMMARY"

# ------------------------------ 运行 hmmbuild -------------------------------
step "运行 hmmbuild"
HMMBUILD_VER="$(hmmbuild -h 2>&1 | grep -oE '# HMMER [0-9.]+' | head -1 || echo 'unknown')"
log "HMMER 版本: $HMMBUILD_VER"

log "命令: hmmbuild --amino --cpu $THREADS '$HMM' '$ALN'"
hmmbuild --amino --cpu "$THREADS" "$HMM" "$ALN"

[[ -s "$HMM" ]] || die "hmmbuild 输出为空: $HMM"
log "hmmbuild 完成，HMM 大小: $(du -h "$HMM" | cut -f1)"
printf 'HMM\tTRM.hmm\t生成成功\t%s\tOK\t\n' "$(du -h "$HMM" | cut -f1)" >> "$SUMMARY"

# ------------------------------ HMM 头部统计 --------------------------------
step "HMM 统计"
grep -E '^(NAME|LENG|ALPH|NSEQ|EFFN|CKSUM)' "$HMM" > "$HDR" 2>/dev/null || true
if [[ -s "$HDR" ]]; then
  sed 's/^/    /' "$HDR"
  HMM_LENG="$(awk '/^LENG/{print $2}' "$HDR")"
  HMM_NSEQ="$(awk '/^NSEQ/{print $2}' "$HDR")"
  HMM_EFFN="$(awk '/^EFFN/{print $2}' "$HDR")"
  printf 'HMM\t比对长度(LENG)\t-\t%s\tINFO\t\n' "${HMM_LENG:-NA}" >> "$SUMMARY"
  printf 'HMM\t建库序列数(NSEQ)\t%d\t%s\tOK\t\n' "$N_SEQ" "${HMM_NSEQ:-NA}" >> "$SUMMARY"
  # EFFN 远低于 NSEQ 说明序列高度冗余，建库权重被大幅下调
  if [[ -n "${HMM_EFFN:-}" ]]; then
    log "有效序列数(EFFN) = $HMM_EFFN"
    printf 'HMM\t有效序列数(EFFN)\t-\t%s\tINFO\t\n' "$HMM_EFFN" >> "$SUMMARY"
  fi
else
  warn "未能解析 HMM 头部统计"
  printf 'HMM\t头部统计\t可解析\t失败\tWARN\t\n' >> "$SUMMARY"
fi

# ------------------------------ 报告 ----------------------------------------
{
  echo "阶段 03：hmmbuild 建立 TRM.hmm 报告"
  echo "======================================================================"
  echo "时间        : $(date '+%F %T')"
  echo "软件        : $HMMBUILD_VER"
  echo "参数        : --amino --cpu $THREADS"
  echo "输入比对    : $ALN  ($N_SEQ 条序列)"
  echo "输出 HMM    : $HMM"
  echo ""
  echo "【HMM 头部统计】"
  sed 's/^/  /' "$HDR" 2>/dev/null || echo "  (未解析到)"
  echo ""
  echo "说明: profile HMM 由 34 条 AtTRM 的 MAFFT 比对直接建立，"
  echo "      未引入外部 seed、未做序列加权以外的调参（hmmbuild 默认使用"
  echo "      Henikoff 位置特异性权重 + 有效序列数校正）。"
  echo ""
  echo "注意: 本 HMM 是「自建家族 HMM」，不等同于 Pfam 官方模型。"
  echo "      Pfam 中与 TRM 相关的模型(如 PF14309)在阶段 06 用 hmmscan 单独核查。"
  echo ""
  echo "下一步: bash run.sh 04   （hmmsearch 扫描 Zunla 全蛋白）"
} > "$REPORT"

log "产出: $HMM"
log "      $HDR"
log "      $REPORT"

}   # run_stage 结束

run_stage 2>&1 | tee "$LOG"
exit "${PIPESTATUS[0]}"
