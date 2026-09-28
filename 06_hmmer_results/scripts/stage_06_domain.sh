#!/usr/bin/env bash
# =============================================================================
# stage_06_domain.sh —— 阶段 06：Pfam 结构域核查
#
# 方法
#   hmmscan --domtblout <domtbl> --cpu N <Pfam-A.hmm> <候选蛋白.fasta>
#   重点核查 PF14309（TRM 家族核心域，官方注释含 TRM32/LONGIFOLIA1/2）
#   与 PF14383（与 TRM 关联未知，用实测数据检验）
#
# 降级路径（缺 Pfam 数据库时）
#   脚本会明确报告为阻塞项，并给出三种可选方案（见报告），不静默跳过。
#
# 用法:
#   bash stage_06_domain.sh [--pfam /path/Pfam-A.hmm] [--threads N]
#
# 产出（06_hmmer_results/06_domain/）
#   06_domain.log
#   CaTRM_candidates_proteins.fasta      送检的候选蛋白
#   CaTRM_vs_Pfam.domtblout              hmmscan 原始结果
#   CaTRM_vs_Pfam.domtblout.parsed.tsv   全部显著命中
#   CaTRM_vs_Pfam.per_protein.tsv        每蛋白结构域清单 + has_TRM_core_domain
#   CaTRM_vs_Pfam.focus_matrix.tsv       候选 × 关注域 矩阵
#   CaTRM_vs_Pfam.focus_hits.tsv         关注域明细
#   CaTRM_vs_Pfam.stats.txt              统计（含 PF14383 实测结论）
#   06_domain.subset.tsv / 06_domain_report.txt
# =============================================================================
set -euo pipefail

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPTS_DIR}/common.sh"

parse_ref_args "$@" || { sed -n '2,22p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0; }

OUTDIR="$D06"
mkdir -p "$OUTDIR"
LOG="${OUTDIR}/06_domain.log"
CAND_FA="${OUTDIR}/CaTRM_candidates_proteins.fasta"
DOMTBL="${OUTDIR}/CaTRM_vs_Pfam.domtblout"
PREFIX="${OUTDIR}/CaTRM_vs_Pfam"
SUMMARY="${OUTDIR}/06_domain.subset.tsv"
REPORT="${OUTDIR}/06_domain_report.txt"
PARSER="${SCRIPTS_DIR}/parse_hmmscan_domtbl.py"
MERGED="$D05/CaTRM_candidates.merged.tsv"
HMM_HIT_FA="$D04/TRM.hmmsearch.hit_proteins.fasta"

subset_init "$SUMMARY"

run_stage() {

step "阶段 06：Pfam 结构域核查"
log "候选来源 : $MERGED"
log "关注域   : $PFAM_FOCUS"

# ------------------------------ 组装候选蛋白 --------------------------------
step "组装候选蛋白 FASTA"
[[ -f "$MERGED" ]] || die "找不到阶段 05 的合并结果: $MERGED
    请先运行: bash run.sh 05"

# 从 merged.tsv 收集蛋白 ID（blast_proteins 逗号分隔列 + hmm_best_protein 列）
awk -F'\t' '
  NR==1 { for(i=1;i<=NF;i++){ if($i=="blast_proteins") b=i; if($i=="hmm_best_protein") h=i } next }
  {
    if(b && $b!="NA"){ n=split($b,arr,","); for(k=1;k<=n;k++) print arr[k] }
    if(h && $h!="NA") print $h
  }
' "$MERGED" | sed 's/[[:space:]]//g' | grep -v '^$' | sort -u > "${OUTDIR}/.want_prot.tmp" || true

N_WANT="$(count_lines "${OUTDIR}/.want_prot.tmp")"
log "merged.tsv 中引用的蛋白 ID 数: $N_WANT"

if [[ ! -f "$ZUNLA_PEP" ]]; then
  die "找不到 Zunla 蛋白库，无法抽取候选序列: $ZUNLA_PEP"
fi
awk '
  BEGIN{FS="[ \t]+"}
  NR==FNR{ want[$1]=1; next }
  /^>/ { id=$0; sub(/^>/,"",id); sub(/[ \t].*$/,"",id); keep=(id in want); if(keep) print; next }
  keep
' "${OUTDIR}/.want_prot.tmp" "$ZUNLA_PEP" > "$CAND_FA"

N_CAND="$(count_fasta "$CAND_FA")"
log "抽取到候选蛋白序列: $N_CAND 条"
printf '候选\t送检蛋白数\t%d\t%d\t%s\t\n' "$N_WANT" "$N_CAND" \
  "$( (( N_CAND == N_WANT )) && echo OK || echo 'WARN:有蛋白未在库中找到' )" >> "$SUMMARY"
if (( N_CAND < N_WANT )); then
  comm -23 "${OUTDIR}/.want_prot.tmp" \
    <(grep '^>' "$ZUNLA_PEP" | sed 's/^>//' | awk '{print $1}' | sort -u) \
    | sed 's/^/[未在Zunla库中找到] /' > "${OUTDIR}/06_missing_proteins.txt" || true
  warn "有 $((N_WANT - N_CAND)) 个蛋白在 Zunla 库中未找到，已列出到 06_missing_proteins.txt"
fi
rm -f "${OUTDIR}/.want_prot.tmp"
(( N_CAND > 0 )) || die "候选蛋白 FASTA 为空，无法进行结构域核查"

# ------------------------------ Pfam 数据库 ---------------------------------
step "定位 Pfam 数据库"
have hmmscan || die "未找到 hmmscan。结构域核查需要 HMMER 的 hmmscan。"

PFAM="$(probe_pfam || true)"
if [[ -z "$PFAM" ]]; then
  warn "未找到 Pfam-A.hmm —— 结构域核查无法进行（这是阻塞项，不是可以跳过的步骤）"
  printf 'Pfam\tPfam-A.hmm\t存在\t未找到\tBLOCK:无法核查结构域\t\n' >> "$SUMMARY"
  cat > "$REPORT" <<EOF
阶段 06：Pfam 结构域核查 —— 未能执行（缺 Pfam 数据库）
======================================================================
候选蛋白已就绪: $CAND_FA ($N_CAND 条)

阻塞原因: 未找到 Pfam-A.hmm，hmmscan 无法运行。

三种可选方案（请选一种后告知，或自行补齐）:

  A. 本地 Pfam-A.hmm（推荐，可批量、可复现）
     wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
     gunzip Pfam-A.hmm.gz && hmmpress Pfam-A.hmm
     然后: bash run.sh 06 --pfam /path/to/Pfam-A.hmm
     注意: 解压约 1.5-2 GB，hmmpress 后更大；请确认磁盘余量。

  B. InterProScan（含 Pfam/CDD/SMART/ProSite，覆盖最全，但体积与耗时大）
     需下载约 10-20 GB 数据；适合需要正式 InterPro 注释的场景。

  C. 在线提交（无需本地安装，但不适合大批量、且结果不可脚本复现）
     InterPro:  https://www.ebi.ac.uk/interpro/search/sequence/
     NCBI CDD:  https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi
     SMART:     http://smart.embl-heidelberg.de/

已产出的中间文件（不浪费）:
  $CAND_FA   <- 可直接用于上述任一方案

关键提醒: PF14309 是 TRM 家族核心域（官方注释含 TRM32/LONGIFOLIA1/2），
          PF14383 与 TRM 关联未知，需用实测结果判断，勿预设。
EOF
  log "报告已写: $REPORT"
  return 0
fi

PFAM_VER_LINE="$(head -20 "$PFAM" | grep -iE '^#=GF (DE|PI)' | head -1 || true)"
log "Pfam 数据库: $PFAM"
log "  $PFAM_VER_LINE"

IDX_MISS=""
for ext in h3f h3i h3m h3p; do
  [[ -f "${PFAM}.${ext}" ]] || IDX_MISS="$IDX_MISS .${ext}"
done
if [[ -n "$IDX_MISS" ]]; then
  warn "Pfam 库缺少 hmmpress 索引:$IDX_MISS"
  if have hmmpress; then
    log "尝试自动运行 hmmpress（需要写权限与磁盘空间）..."
    if hmmpress -f "$PFAM"; then
      log "hmmpress 完成"
    else
      die "hmmpress 失败。请手动运行: hmmpress -f $PFAM"
    fi
  else
    die "缺 hmmpress 索引且系统无 hmmpress 命令。请先 hmmpress -f $PFAM"
  fi
fi
printf 'Pfam\t数据库\t存在\t%s\tOK\t\n' "$(basename "$PFAM")" >> "$SUMMARY"

# ------------------------------ 运行 hmmscan --------------------------------
step "运行 hmmscan"
HMMER_VER="$(hmmscan -h 2>&1 | grep -oE '# HMMER [0-9.]+' | head -1 || echo 'unknown')"
log "HMMER 版本: $HMMER_VER"
log "命令: hmmscan --domtblout $DOMTBL --cpu $THREADS $PFAM $CAND_FA"

hmmscan --domtblout "$DOMTBL" --cpu "$THREADS" "$PFAM" "$CAND_FA" \
  > "${OUTDIR}/CaTRM_vs_Pfam.hmmscan.txt" 2>&1

[[ -s "$DOMTBL" ]] || die "hmmscan 未产出 domtblout: $DOMTBL"
log "hmmscan 完成，domtblout 命中行数: $(grep -vc '^#' "$DOMTBL" || true)"
printf 'hmmscan\t原始结果\t保留\tdomtblout\tOK\t\n' >> "$SUMMARY"

# ------------------------------ 解析 ----------------------------------------
step "解析结构域结果"
PY="$(find_python || true)"
[[ -n "$PY" ]] || die "需要 python3，但未找到"
[[ -f "$PARSER" ]] || die "找不到解析脚本: $PARSER"

"$PY" "$PARSER" \
  --domtbl "$DOMTBL" \
  --prefix "$PREFIX" \
  --dom-eval "$DOM_EVAL" \
  --focus "$(printf '%s' "$PFAM_FOCUS" | tr ' ' ',')"

# 记录关注域命中数
FOCUS_HITS="${PREFIX}.focus_hits.tsv"
if [[ -f "$FOCUS_HITS" ]]; then
  for f in $PFAM_FOCUS; do
    n="$(awk -F'\t' -v F="$f" 'NR>1 && $1==F' "$FOCUS_HITS" | wc -l | tr -d ' ')"
    printf '结构域\t%s 显著命中\t-\t%d\tINFO\t\n' "$f" "$n" >> "$SUMMARY"
  done
fi
N_CORE="$(awk -F'\t' 'NR>1 && $5=="yes"' "${PREFIX}.per_protein.tsv" 2>/dev/null | wc -l | tr -d ' ')"
printf '结构域\t含 PF14309(TRM核心域) 的蛋白\t-\t%d\tINFO\t\n' "$N_CORE" >> "$SUMMARY"

# ------------------------------ 报告 ----------------------------------------
{
  echo "阶段 06：Pfam 结构域核查 —— 报告"
  echo "======================================================================"
  echo "时间        : $(date '+%F %T')"
  echo "软件        : $HMMER_VER"
  echo "参数        : hmmscan --domtblout --cpu $THREADS"
  echo "Pfam 数据库 : $PFAM"
  echo "  $PFAM_VER_LINE"
  echo "候选蛋白    : $CAND_FA ($N_CAND 条)"
  echo ""
  cat "${PREFIX}.stats.txt"
  echo ""
  echo "【关注域矩阵（候选 × 关注域）】"
  if [[ -f "${PREFIX}.focus_matrix.tsv" ]]; then
    column -t "${PREFIX}.focus_matrix.tsv" 2>/dev/null | sed 's/^/  /' \
      || sed 's/^/  /' "${PREFIX}.focus_matrix.tsv"
  else
    echo "  (无命中)"
  fi
  echo ""
  echo "【每蛋白结构域（前 30 行）】"
  head -31 "${PREFIX}.per_protein.tsv" 2>/dev/null | column -t 2>/dev/null | sed 's/^/  /' \
    || head -31 "${PREFIX}.per_protein.tsv" | sed 's/^/  /'
  echo ""
  echo "【CDD / SMART 补充核查（如需更全面）】"
  echo "  本步骤用 Pfam-A 完成主核查。若需 CDD/SMART 交叉验证，可选:"
  echo "    - NCBI CDD 在线批量: https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"
  echo "    - InterProScan 本地版（含 CDD/SMART/PROSITE/PRINTS）"
  echo "  候选蛋白 FASTA 已备好: $CAND_FA"
  echo ""
  echo "下一步: bash run.sh 07   （汇总最终 CaTRM 列表）"
} > "$REPORT"

log "产出: $CAND_FA"
log "      $DOMTBL"
log "      ${PREFIX}.per_protein.tsv"
log "      $REPORT"

}   # run_stage 结束

run_stage 2>&1 | tee "$LOG"
exit "${PIPESTATUS[0]}"
