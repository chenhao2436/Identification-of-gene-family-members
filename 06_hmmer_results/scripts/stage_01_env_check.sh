#!/usr/bin/env bash
# =============================================================================
# stage_01_env_check.sh —— 阶段 01：环境与输入文件检查（只读，不做分析）
#
# 检查内容
#   A 软件环境      mafft / hmmbuild / hmmsearch / hmmscan / hmmpress / interproscan / blast
#   B Pfam 库       位置、版本、hmmpress 索引是否齐全
#   C 输入文件      存在性、格式、序列条数（失败样本会被抽样展示）
#   D ID 交叉核对   24 个候选基因对 GFF3、39 个候选蛋白对 Zunla 蛋白库
#   E BLAST 核查    从日志提取实际使用的阈值、候选数自洽性
#   F 磁盘余量
#
# 用法:
#   bash stage_01_env_check.sh
#   bash stage_01_env_check.sh --gff3 /path/Canz.genome.gff3 --genome /path/Canz.genome.fa
#
# 产出（06_hmmer_results/01_env_check/）
#   01_env_check.log            完整日志
#   01_env_check_report.txt     中文结论（齐备/缺失/阻塞项）
#   01_env_check_summary.tsv    机读小结
# =============================================================================
set -euo pipefail

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPTS_DIR}/common.sh"

# 允许 --gff3/--genome/--pfam/--threads 覆盖
parse_ref_args "$@" || { sed -n '2,30p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0; }

OUTDIR="$D01"
mkdir -p "$OUTDIR"
LOG="${OUTDIR}/01_env_check.log"
REPORT="${OUTDIR}/01_env_check_report.txt"
SUMMARY="${OUTDIR}/01_env_check_summary.tsv"
SAMPLES="${OUTDIR}/01_failed_id_samples.txt"

subset_init "$SUMMARY"
: > "$SAMPLES"

BLOCKERS=0      # 阻塞项计数
WARNINGS=0      # 警示项计数

# 把检查结论同时写入小结与屏幕
rec() {   # rec 类别 条目 期望 实测 结论 [是否阻塞]
  local cat="$1" item="$2" exp="$3" obs="$4" concl="$5" block="${6:-}"
  local last="${concl##*:}"
  printf '%s\t%s\t%s\t%s\t%s\n' "$cat" "$item" "$exp" "$obs" "$concl" >> "$SUMMARY"
  if [[ "$concl" == *"OK"* ]]; then
    ok "$item = $obs"
  elif [[ "$concl" == BLOCK* ]]; then
    bad "$item = $obs   <-- 阻塞"
    (( BLOCKERS++ )) || true
  else
    warn "$item = $obs"
    (( WARNINGS++ )) || true
  fi
}

# =============================================================================
step "阶段 01：环境与输入文件检查"
echo "项目根目录 : $PROJECT_ROOT"
echo "流程目录   : ${HMMER_DIR#"$PROJECT_ROOT"/}"
echo "输出目录   : ${OUTDIR#"$PROJECT_ROOT"/}"
echo "时间       : $(date '+%F %T')"
echo "主机       : $(uname -n 2>/dev/null || echo '?')"

# =============================================================================
step "A. 软件环境"
# =============================================================================
SW_RAW="${OUTDIR}/01_software_versions.txt"
: > "$SW_RAW"

check_sw() {   # check_sw <命令> <必需?> <版本参数...>
  local cmd="$1" required="$2"; shift 2
  if have "$cmd"; then
    local p v
    p="$(command -v "$cmd")"
    v="$(show_version "$cmd" "$@")"
    printf '%-14s %s\t%s\n' "$cmd" "$p" "$v" >> "$SW_RAW"
    echo "  [OK]   $cmd = $v"
    printf '软件\t%s\t-\t%s\tOK\t\n' "$cmd" "$p" >> "$SUMMARY"
    return 0
  else
    printf '%-14s MISSING\n' "$cmd" >> "$SW_RAW"
    echo "  [缺失] $cmd"
    if [[ "$required" == "required" ]]; then
      printf '软件\t%s\t需要安装\t缺失\tBLOCK:必需软件缺失\t\n' "$cmd" >> "$SUMMARY"
      (( BLOCKERS++ )) || true
    else
      printf '软件\t%s\t可选\t缺失\tWARN:可选软件缺失\t\n' "$cmd" >> "$SUMMARY"
      (( WARNINGS++ )) || true
    fi
    return 1
  fi
}

check_sw mafft      required --version
check_sw hmmbuild   required -h
check_sw hmmsearch  required -h
check_sw hmmscan    optional -h
check_sw hmmpress   optional -h
check_sw esl-alistat optional -h
check_sw blastp     optional -version
check_sw makeblastdb optional -version
check_sw interproscan.sh optional --version
check_sw emapper.py optional --version

PY="$(find_python || true)"
if [[ -n "$PY" ]]; then
  echo "  [OK]   $PY = $("$PY" --version 2>&1)"
  printf '软件\t%s\t-\t%s\tOK\t\n' "$PY" "$("$PY" --version 2>&1)" >> "$SUMMARY"
else
  echo "  [缺失] python3（解析脚本需要）"
  printf '软件\tpython3\t需要安装\t缺失\tBLOCK:解析脚本依赖\t\n' >> "$SUMMARY"
  (( BLOCKERS++ )) || true
fi

# =============================================================================
step "B. Pfam 库"
# =============================================================================
PFAM_FOUND="$(probe_pfam || true)"
if [[ -n "$PFAM_FOUND" ]]; then
  PFAM_SIZE="$(du -h "$PFAM_FOUND" 2>/dev/null | cut -f1)"
  PFAM_VER="$(head -20 "$PFAM_FOUND" 2>/dev/null | grep -iE '^#=GF (DE|PI)' | head -1 | sed 's/^#=GF //')"
  echo "  位置: $PFAM_FOUND  ($PFAM_SIZE)"
  echo "  $PFAM_VER"
  printf 'Pfam\tPfam-A.hmm\t存在\t%s\tOK:%s\t\n' "$PFAM_FOUND" "$PFAM_SIZE" >> "$SUMMARY"

  # hmmpress 索引检查
  idx_missing=""
  for ext in h3f h3i h3m h3p; do
    [[ -f "${PFAM_FOUND}.${ext}" ]] || idx_missing="${idx_missing} .${ext}"
  done
  if [[ -z "$idx_missing" ]]; then
    echo "  [OK]   hmmpress 索引齐全 (.h3f/.h3i/.h3m/.h3p)"
    printf 'Pfam\thmmpress索引\t齐全\t齐全\tOK\t\n' >> "$SUMMARY"
  else
    echo "  [缺失] hmmpress 索引:$idx_missing  → 需先运行 hmmpress"
    printf 'Pfam\thmmpress索引\t齐全\t缺少%s\tBLOCK:需先hmmpress\t\n' "$idx_missing" >> "$SUMMARY"
    (( BLOCKERS++ )) || true
  fi
else
  echo "  [缺失] 未找到 Pfam-A.hmm"
  echo "         阶段 06 结构域核查将无法运行 hmmscan（可改用 InterProScan）"
  printf 'Pfam\tPfam-A.hmm\t存在\t未找到\tBLOCK:阶段06依赖\t\n' >> "$SUMMARY"
  (( BLOCKERS++ )) || true
fi

# =============================================================================
step "C. 输入文件"
# =============================================================================
check_fasta() {   # check_fasta <路径> <期望条数|-> <说明> <必需?>
  local f="$1" exp="$2" label="$3" req="${4:-required}"
  if [[ ! -f "$f" ]]; then
    echo "  [缺失] $label -> $f"
    printf '输入\t%s\t存在\t未找到\tBLOCK:文件缺失\t\n' "$label" >> "$SUMMARY"
    (( BLOCKERS++ )) || true
    return 1
  fi
  if ! is_text_fasta "$f"; then
    echo "  [异常] $label 不是纯文本 FASTA -> $f"
    printf '输入\t%s\t纯文本FASTA\t格式异常\tBLOCK:格式错误\t\n' "$label" >> "$SUMMARY"
    (( BLOCKERS++ )) || true
    return 1
  fi
  local n sz
  n="$(count_fasta "$f")"
  sz="$(du -h "$f" 2>/dev/null | cut -f1)"
  if [[ "$exp" == "-" ]]; then
    echo "  [OK]   $label = $n 条 ($sz)"
    printf '输入\t%s\t有序列\t%d条/%s\tOK\t\n' "$label" "$n" "$sz" >> "$SUMMARY"
  elif (( n == exp )); then
    echo "  [OK]   $label = $n 条 ($sz, 期望 $exp)"
    printf '输入\t%s\t%d\t%d\tOK\t\n' "$label" "$exp" "$n" >> "$SUMMARY"
  else
    echo "  [警告] $label = $n 条 (期望 $exp, $sz)"
    printf '输入\t%s\t%d\t%d\tWARN:条数不符\t\n' "$label" "$exp" "$n" >> "$SUMMARY"
    (( WARNINGS++ )) || true
  fi
  return 0
}

check_fasta "$AT_QUERY_FA" "$EXPECT_AT_N" "AtTRM query FASTA" required || true
check_fasta "$ZUNLA_PEP"  "-"            "Zunla 全蛋白 FASTA" required || true

# 注释表
if [[ -f "$AT_ANNO_TSV" ]]; then
  cols="$(head -1 "$AT_ANNO_TSV" | awk -F'\t' '{print NF}')"
  rows="$(wc -l < "$AT_ANNO_TSV" | tr -d ' ')"
  echo "  [OK]   AtTRM 注释表 = $((rows-1)) 行数据 × ${cols} 列"
  printf '输入\tAtTRM注释表\t9列×34行\t%d列×%d行\tOK\t\n' "$cols" "$((rows-1))" >> "$SUMMARY"
else
  echo "  [缺失] AtTRM 注释表 -> $AT_ANNO_TSV"
  printf '输入\tAtTRM注释表\t存在\t未找到\tWARN:At基因映射会填NA\t\n' >> "$SUMMARY"
  (( WARNINGS++ )) || true
fi

# GFF3 / genome
GFF3_FOUND="$(probe_gff3 || true)"
GENOME_FOUND="$(probe_genome || true)"
if [[ -n "$GFF3_FOUND" ]]; then
  GZ="否"; CAT="cat"
  case "$GFF3_FOUND" in *.gz) GZ="是"; CAT="zcat" ;; esac
  echo "  [OK]   GFF3 = $GFF3_FOUND (压缩: $GZ, $(du -h "$GFF3_FOUND" 2>/dev/null | cut -f1))"
  printf '输入\tZunla GFF3\t存在\t%s\tOK\t\n' "$GFF3_FOUND" >> "$SUMMARY"
else
  echo "  [缺失] 未找到 Zunla GFF3（用 --gff3 指定）"
  printf '输入\tZunla GFF3\t存在\t未找到\tWARN:gene/transcript映射受限\t\n' >> "$SUMMARY"
  (( WARNINGS++ )) || true
fi
if [[ -n "$GENOME_FOUND" ]]; then
  echo "  [OK]   genome.fa = $GENOME_FOUND ($(du -h "$GENOME_FOUND" 2>/dev/null | cut -f1))"
  [[ -f "${GENOME_FOUND}.fai" ]] && echo "  [OK]   .fai 索引存在" \
    || echo "  [警告] 缺 .fai 索引"
  printf '输入\tZunla genome\t存在\t%s\tOK\t\n' "$GENOME_FOUND" >> "$SUMMARY"
else
  echo "  [警告] 未找到 Zunla genome.fa（用 --genome 指定；本流程不直接用，仅登记）"
  printf '输入\tZunla genome\t存在\t未找到\tWARN:本流程不直接使用\t\n' >> "$SUMMARY"
  (( WARNINGS++ )) || true
fi

# Zunla 蛋白库内部一致性(上一轮出现过 grep 计数与 makeblastdb 计数不一致)
if [[ -f "$ZUNLA_PEP" ]]; then
  n_gt="$(count_fasta "$ZUNLA_PEP")"
  n_uniq="$(grep '^>' "$ZUNLA_PEP" 2>/dev/null | sed 's/^>//' | awk '{print $1}' | sort -u | wc -l | tr -d ' ')"
  echo "  Zunla 蛋白: '>' 计数=$n_gt, 唯一 ID 数=$n_uniq"
  if (( n_gt == n_uniq )); then
    printf '输入\tZunla蛋白ID唯一性\t%d\t%d\tOK\t\n' "$n_gt" "$n_uniq" >> "$SUMMARY"
  else
    echo "  [警告] 存在重复 ID($((n_gt-n_uniq)) 个), makeblastdb 会合并它们"
    printf '输入\tZunla蛋白ID唯一性\t唯一\t%d个重复\tWARN:建库会合并\t\n' "$((n_gt-n_uniq))" >> "$SUMMARY"
    (( WARNINGS++ )) || true
  fi
  # 空序列检查
  n_empty="$(awk '/^>/{if(n!=""&&len==0)c++; n=$0; len=0; next}{len+=length($0)}END{if(n!=""&&len==0)c++; print c+0}' "$ZUNLA_PEP")"
  echo "  Zunla 蛋白: 空序列数=$n_empty"
  if (( n_empty > 0 )); then
    printf '输入\tZunla蛋白空序列\t0\t%d\tWARN:会导致建库条数不一致\t\n' "$n_empty" >> "$SUMMARY"
    (( WARNINGS++ )) || true
  else
    printf '输入\tZunla蛋白空序列\t0\t0\tOK\t\n' >> "$SUMMARY"
  fi
fi

# =============================================================================
step "D. ID 交叉核对（protein / transcript / gene）"
# =============================================================================
# D1: 候选基因 ID -> GFF3 gene 特征
if [[ -n "$GFF3_FOUND" && -f "$BLAST_GENE_LIST" ]]; then
  n_gene="$(count_lines "$BLAST_GENE_LIST")"
  : > "${OUTDIR}/.gff3_gene_ids.tmp"
  case "$GFF3_FOUND" in
    *.gz) zcat "$GFF3_FOUND" ;;
    *)    cat "$GFF3_FOUND" ;;
  esac | awk -F'\t' '$3=="gene"{ for(i=1;i<=NF;i++) if($i ~ /^ID=/){ sub(/^ID=/,"",$i); sub(/;.*/,"",$i); print $i } }' \
       > "${OUTDIR}/.gff3_gene_ids.tmp" || true
  n_in_gff="$(count_lines "${OUTDIR}/.gff3_gene_ids.tmp")"
  n_match="$(comm -12 <(sort -u "$BLAST_GENE_LIST") <(sort -u "${OUTDIR}/.gff3_gene_ids.tmp") | wc -l | tr -d ' ')"
  echo "  GFF3 gene 特征数 = $n_in_gff"
  echo "  候选基因($n_gene) 中能在 GFF3 找到 gene 特征: $n_match"
  if (( n_match == n_gene )); then
    printf 'ID映射\t候选基因->GFF3 gene\t%d\t%d\tOK\t\n' "$n_gene" "$n_match" >> "$SUMMARY"
  else
    printf 'ID映射\t候选基因->GFF3 gene\t%d\t%d\tWARN:有未匹配\t\n' "$n_gene" "$n_match" >> "$SUMMARY"
    (( WARNINGS++ )) || true
    comm -23 <(sort -u "$BLAST_GENE_LIST") <(sort -u "${OUTDIR}/.gff3_gene_ids.tmp") \
      | sed 's/^/[未匹配GFF3] /' >> "$SAMPLES"
  fi
  rm -f "${OUTDIR}/.gff3_gene_ids.tmp"
else
  printf 'ID映射\t候选基因->GFF3 gene\t-\t跳过\tWARN:缺GFF3或候选列表\t\n' >> "$SUMMARY"
  (( WARNINGS++ )) || true
fi

# D2: 候选蛋白 ID -> Zunla 蛋白库
if [[ -f "$ZUNLA_PEP" && -f "$BLAST_PROT_LIST" ]]; then
  n_prot="$(count_lines "$BLAST_PROT_LIST")"
  grep '^>' "$ZUNLA_PEP" | sed 's/^>//' | awk '{print $1}' | sort -u > "${OUTDIR}/.zunla_ids.tmp"
  n_match="$(comm -12 <(sort -u "$BLAST_PROT_LIST") "${OUTDIR}/.zunla_ids.tmp" | wc -l | tr -d ' ')"
  echo "  候选蛋白($n_prot) 中能在 Zunla 蛋白库找到: $n_match"
  if (( n_match == n_prot )); then
    printf 'ID映射\t候选蛋白->Zunla库\t%d\t%d\tOK\t\n' "$n_prot" "$n_match" >> "$SUMMARY"
  else
    printf 'ID映射\t候选蛋白->Zunla库\t%d\t%d\tWARN:有未匹配\t\n' "$n_prot" "$n_match" >> "$SUMMARY"
    (( WARNINGS++ )) || true
    comm -23 <(sort -u "$BLAST_PROT_LIST") "${OUTDIR}/.zunla_ids.tmp" \
      | sed 's/^/[未匹配Zunla] /' >> "$SAMPLES"
  fi
  rm -f "${OUTDIR}/.zunla_ids.tmp"
else
  printf 'ID映射\t候选蛋白->Zunla库\t-\t跳过\tWARN:缺文件\t\n' >> "$SUMMARY"
fi

# D3: GFF3 mRNA 与蛋白 ID 是否一一对应
if [[ -n "$GFF3_FOUND" && -f "$ZUNLA_PEP" ]]; then
  case "$GFF3_FOUND" in
    *.gz) zcat "$GFF3_FOUND" ;;
    *)    cat "$GFF3_FOUND" ;;
  esac | awk -F'\t' '$3=="mRNA"{ for(i=1;i<=NF;i++) if($i ~ /^ID=/){ sub(/^ID=/,"",$i); sub(/;.*/,"",$i); print $i } }' \
       | sort -u > "${OUTDIR}/.gff3_mrna.tmp" || true
  n_mrna="$(count_lines "${OUTDIR}/.gff3_mrna.tmp")"
  n_pep="$(count_fasta "$ZUNLA_PEP")"
  n_common="$(comm -12 "${OUTDIR}/.gff3_mrna.tmp" <(grep '^>' "$ZUNLA_PEP" | sed 's/^>//' | awk '{print $1}' | sort -u) | wc -l | tr -d ' ')"
  echo "  GFF3 mRNA 数 = $n_mrna ; 蛋白库序列数 = $n_pep ; 两者 ID 交集 = $n_common"
  if (( n_common > 0 )); then
    printf 'ID映射\tmRNA<->蛋白ID\t>0\t%d\tOK:可一一对应\t\n' "$n_common" >> "$SUMMARY"
  else
    printf 'ID映射\tmRNA<->蛋白ID\t>0\t0\tWARN:命名体系可能不一致\t\n' >> "$SUMMARY"
    (( WARNINGS++ )) || true
  fi
  rm -f "${OUTDIR}/.gff3_mrna.tmp"
fi

# =============================================================================
step "E. 已有 BLAST 结果核查"
# =============================================================================
if [[ -f "$BLAST_MAIN" ]]; then
  n_hsp="$(wc -l < "$BLAST_MAIN" | tr -d ' ')"
  echo "  [OK]   主候选 HSP 表 = $n_hsp 行 -> $BLAST_MAIN"
  printf 'BLAST\t主候选HSP表\t存在\t%d行\tOK\t\n' "$n_hsp" >> "$SUMMARY"

  # 从 HSP 表反推实际使用的阈值（这是唯一权威来源）
  awk -F'\t' '
    NR==1{next}
    { if($3+0<imin||imin==""){imin=$3+0} if($3+0>imax)imax=$3+0
      if($15+0<cmin||cmin==""){cmin=$15+0} if($15+0>cmax)cmax=$15+0
      if($11+0>emax||emax==""){emax=$11+0} }
    END{ printf "  BLAST 实际过滤后: identity %.1f-%.1f%%, coverage %.0f-%.0f%%, E 最大 %.2g\n", imin,imax,cmin,cmax,emax }
  ' "$BLAST_MAIN"
else
  echo "  [缺失] 主候选 HSP 表 -> $BLAST_MAIN"
  printf 'BLAST\t主候选HSP表\t存在\t未找到\tBLOCK:阶段05依赖\t\n' >> "$SUMMARY"
  (( BLOCKERS++ )) || true
fi

# 候选数三者自洽性
n_gene_l="$( [[ -f "$BLAST_GENE_LIST" ]] && count_lines "$BLAST_GENE_LIST" || echo 0 )"
n_prot_l="$( [[ -f "$BLAST_PROT_LIST" ]] && count_lines "$BLAST_PROT_LIST" || echo 0 )"
n_prot_fa="$( [[ -f "${BLAST_DIR}/Zunla_TRM_candidate_proteins.fasta" ]] && count_fasta "${BLAST_DIR}/Zunla_TRM_candidate_proteins.fasta" || echo 0)"
echo "  BLAST 候选: 基因=$n_gene_l, 蛋白列表=$n_prot_l, 蛋白FASTA=$n_prot_fa"
if (( n_prot_l == n_prot_fa && n_prot_l > 0 )); then
  printf 'BLAST\t候选数自洽\t蛋白列表=FASTA\t%d=%d\tOK\t\n' "$n_prot_l" "$n_prot_fa" >> "$SUMMARY"
else
  printf 'BLAST\t候选数自洽\t蛋白列表=FASTA\t%d vs %d\tWARN:不一致\t\n' "$n_prot_l" "$n_prot_fa" >> "$SUMMARY"
  (( WARNINGS++ )) || true
fi
printf 'BLAST\t候选基因数\t-\t%d\tINFO\t\n' "$n_gene_l" >> "$SUMMARY"
printf 'BLAST\t候选蛋白数\t-\t%d\tINFO\t\n' "$n_prot_l" >> "$SUMMARY"

# 提取日志中记录的实际参数
if [[ -f "$BLAST_DIAG" ]]; then
  echo "  从流水线日志提取的实际参数:"
  grep -E '阈值|db    序列数|query 序列数|库信息|候选蛋白\(|候选基因\(' "$BLAST_DIAG" 2>/dev/null \
    | tail -8 | sed 's/^/    /' || true
else
  echo "  [警告] 未找到 $BLAST_DIAG（无法回溯实际参数）"
  (( WARNINGS++ )) || true
fi

# =============================================================================
step "F. 磁盘余量"
# =============================================================================
AVAIL="$(df -h "$HMMER_DIR" 2>/dev/null | awk 'NR==2{print $4}')"
AVAIL_KB="$(df -Pk "$HMMER_DIR" 2>/dev/null | awk 'NR==2{print $4}')"
echo "  可用空间: ${AVAIL:-未知}"
if [[ -n "${AVAIL_KB:-}" ]] && (( AVAIL_KB > 20*1024*1024 )); then
  printf '磁盘\t可用空间\t>20GB\t%s\tOK\t\n' "$AVAIL" >> "$SUMMARY"
else
  printf '磁盘\t可用空间\t>20GB\t%s\tWARN:可能不足以装Pfam/InterProScan\t\n' "${AVAIL:-未知}" >> "$SUMMARY"
  (( WARNINGS++ )) || true
fi

# =============================================================================
step "G. 结论"
# =============================================================================
{
  echo "CaTRM 家族鉴定 —— 阶段 01 环境与输入文件检查报告"
  echo "======================================================================"
  echo "时间       : $(date '+%F %T')"
  echo "项目根目录 : $PROJECT_ROOT"
  echo "阻塞项     : $BLOCKERS"
  echo "警示项     : $WARNINGS"
  echo ""
  echo "【软件】"
  cat "$SW_RAW" | sed 's/\t/  /' | sed 's/^/  /'
  echo ""
  echo "【阻塞项与警示项明细】"
  awk -F'\t' 'NR>1 && $5 !~ /^OK/ { printf "  [%s] %s: 期望 %s, 实测 %s (%s)\n", $1, $2, $3, $4, $5 }' "$SUMMARY"
  echo ""
  echo "【ID 映射】"
  awk -F'\t' 'NR>1 && $1=="ID映射" { printf "  %s: %s -> %s (%s)\n", $2, $3, $4, $5 }' "$SUMMARY"
  echo ""
  echo "【BLAST 现状】"
  awk -F'\t' 'NR>1 && $1=="BLAST" { printf "  %s: %s (%s)\n", $2, $4, $5 }' "$SUMMARY"
  echo ""
  if [[ -s "$SAMPLES" ]]; then
    echo "【未匹配 ID 抽样(最多 20 条)】"
    head -20 "$SAMPLES" | sed 's/^/  /'
    echo ""
  fi
  echo "【下一步】"
  if (( BLOCKERS == 0 )); then
    echo "  无阻塞项，可执行:  bash run.sh 02"
  else
    echo "  存在 $BLOCKERS 个阻塞项，请先解决后再进入阶段 02。"
    echo "  缺失软件获取建议见 06_hmmer_results/README_hmmer.md"
  fi
} > "$REPORT"

echo "阻塞项: $BLOCKERS   警示项: $WARNINGS"
echo ""
echo "产出:"
echo "  ${LOG#"$PROJECT_ROOT"/}"
echo "  ${REPORT#"$PROJECT_ROOT"/}"
echo "  ${SUMMARY#"$PROJECT_ROOT"/}"
[[ -s "$SAMPLES" ]] && echo "  ${SAMPLES#"$PROJECT_ROOT"/}"

# 只读检查：始终以 0 退出，让调用方按报告内容决定是否继续
exit 0
