#!/usr/bin/env bash
# =============================================================================
# common.sh —— CaTRM 家族鉴定流程的共享配置与工具函数
#
# 被 scripts/run.sh 与 scripts/stage_*.sh 通过 source 加载，不单独执行。
#
# 设计要点
#   1. 路径全部自动识别：以 PROJECT_ROOT 为基准，在任意目录调用都能工作
#   2. 每个阶段产出 *.subset.tsv（机读小结），便于逐轮追踪候选数变化
#   3. 所有阶段默认只读已有结果，不覆盖、不删改输入
# =============================================================================

# 本文件是库，设计上只被 source。若被直接执行则给出明确提示（也便于 shellcheck）。
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  echo "common.sh 是共享库，请勿直接执行。用法:" >&2
  echo "  bash run.sh <01-07|all>" >&2
  exit 2
fi

# ------------------------------ 路径 ----------------------------------------
SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HMMER_DIR="$(dirname "$SCRIPTS_DIR")"
PROJECT_ROOT="$(dirname "$HMMER_DIR")"

# 输入（可用环境变量覆盖，便于服务器上文件在别处时调整）
TRM_INPUT_DIR="${TRM_INPUT_DIR:-${PROJECT_ROOT}/01_At_TRM_query}"
SCRIPT_OUT_DIR="${SCRIPT_OUT_DIR:-${PROJECT_ROOT}/02_script_output}"
BLAST_DIR="${BLAST_DIR:-${PROJECT_ROOT}/05_blast_results}"

AT_QUERY_FA="${AT_QUERY_FA:-${SCRIPT_OUT_DIR}/01_At_TRM_query_34.fasta}"
AT_ANNO_TSV="${AT_ANNO_TSV:-${SCRIPT_OUT_DIR}/02_At_TRM_annotation.tsv}"
ZUNLA_PEP="${ZUNLA_PEP:-${TRM_INPUT_DIR}/Zunla_Canz.pep.fa}"

# 参考注释（服务器上若在别处，用 --gff3/--genome 传入或设 GFF3_FILE/GENOME_FA）
GFF3_FILE="${GFF3_FILE:-}"
GENOME_FA="${GENOME_FA:-}"

# 输出阶段目录
D01="${HMMER_DIR}/01_env_check"
D02="${HMMER_DIR}/02_msa"
D03="${HMMER_DIR}/03_hmmbuild"
D04="${HMMER_DIR}/04_hmmsearch"
D05="${HMMER_DIR}/05_merge"
D06="${HMMER_DIR}/06_domain"
D07="${HMMER_DIR}/07_final"

# ------------------------------ 可调参数 ------------------------------------
EXPECT_AT_N=34                 # AtTRM query 条数
# hmmsearch 初筛 E 值（标准做法；不为了凑某个预设数量而调整）
HMM_EVAL="${HMM_EVAL:-1e-5}"
# hmmscan 结构域初筛
DOM_EVAL="${DOM_EVAL:-1e-5}"
# 关注的 Pfam 结构域（PF14383 待实测验证，见 README）
PFAM_FOCUS="PF14309 PF14383"
PFAM_DB="${PFAM_DB:-}"
# BLAST 候选（现有结果，不重跑）
BLAST_MAIN="${BLAST_DIR}/TRM_candidates_main.tsv"
BLAST_GENE_LIST="${BLAST_DIR}/TRM_candidates_gene_level.txt"
BLAST_PROT_LIST="${BLAST_DIR}/TRM_candidates_protein_level.txt"
BLAST_SUMMARY="${BLAST_DIR}/01_summary_stats.txt"
BLAST_DIAG="${BLAST_DIR}/01_diagnostics.txt"

# ------------------------------ 日志 ----------------------------------------
log()  { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }
step() { printf '\n========== %s ==========\n' "$*"; }
warn() { printf '[警告] %s\n' "$*" >&2; }
die()  { printf '\n[错误] %s\n' "$*" >&2; exit 1; }
ok()   { printf '  [OK]   %s\n' "$*"; }
bad()  { printf '  [缺失] %s\n' "$*"; }

# ------------------------------ 机读小结 ------------------------------------
# 每个阶段写 <outdir>/<stage>.subset.tsv，列固定，便于逐轮追踪
subset_init() {
  local f="$1"
  printf 'check\titem\texpected\tobserved\tconclusion\n' > "$f"
}
# subset_add <小结文件> <检查类别> <条目> <期望> <实测> <结论>
subset_add() {
  printf '%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5" "${6:-}" >> "$7"
}

# ------------------------------ 通用工具 ------------------------------------
have() { command -v "$1" >/dev/null 2>&1; }

# 取版本：依次尝试多种 --version 形式
show_version() {
  local cmd="$1"; shift
  local out
  out="$("$cmd" "$@" 2>&1 | head -2 | tr '\n' ' ' | sed 's/  */ /g')" || true
  printf '%s' "$out"
}

# 数 FASTA 序列条数（空文件返回 0）
count_fasta() {
  local n
  n="$(grep -c '^>' "$1" 2>/dev/null || true)"
  n="${n//[^0-9]/}"
  printf '%s' "${n:-0}"
}

# 数非空行数（空文件/无匹配都返回 0）
# 必须内部吞掉 grep 的退出码 1：在 set -e 下，$(grep -c ...) 在无匹配时
# 返回 1 会让整个脚本提前中止（这是踩过的坑，勿改）。
count_lines() {
  local n
  n="$(grep -c . "$1" 2>/dev/null || true)"
  n="${n//[^0-9]/}"
  printf '%s' "${n:-0}"
}

# 检查文件是否为纯文本 FASTA（防 SnapGene 二进制）
is_text_fasta() {
  local f="$1"
  [[ -f "$f" ]] || return 1
  if have file; then
    file -b "$f" | grep -qiE 'text|ascii|utf-8|unicode' || return 1
  fi
  head -c 1 "$f" | grep -q '>' || return 1
  return 0
}

# 剥掉 isoform 版本号：ZLC01G0003110.1 -> ZLC01G0003110
strip_ver() { printf '%s' "$1" | sed 's/\.[0-9]\+$//'; }
export -f strip_ver 2>/dev/null || true

# 检查 Python 可用性，输出解释器路径
find_python() {
  local c
  for c in python3 python; do
    if have "$c"; then printf '%s' "$c"; return 0; fi
  done
  return 1
}

# ------------------------------ 参数解析helper ------------------------------
# 各阶段脚本共用的 --gff3/--genome 解析
parse_ref_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --gff3)   GFF3_FILE="$2";   shift 2 ;;
      --genome) GENOME_FA="$2";   shift 2 ;;
      --pfam)   PFAM_DB="$2";     shift 2 ;;
      --threads) THREADS="$2";    shift 2 ;;
      -h|--help) return 1 ;;
      *) shift ;;
    esac
  done
  return 0
}

THREADS="${THREADS:-4}"

# ------------------------------ 自动探测参考文件 ----------------------------
# 在常见位置找 Zunla GFF3（含 .gz）
probe_gff3() {
  [[ -n "$GFF3_FILE" && -f "$GFF3_FILE" ]] && { printf '%s' "$GFF3_FILE"; return 0; }
  local c
  for c in \
      "${TRM_INPUT_DIR}/Canz.genome.gff3" \
      "${TRM_INPUT_DIR}/Canz.genome.gff3.gz" \
      "${TRM_INPUT_DIR}/Canz.gff3" \
      "${PROJECT_ROOT}/Canz.genome.gff3" \
      "${PROJECT_ROOT}/01_At_TRM_query/Canz.genome.gff3"; do
    [[ -f "$c" ]] && { printf '%s' "$c"; return 0; }
  done
  return 1
}

probe_genome() {
  [[ -n "$GENOME_FA" && -f "$GENOME_FA" ]] && { printf '%s' "$GENOME_FA"; return 0; }
  local c
  for c in \
      "${TRM_INPUT_DIR}/Canz.genome.fa" \
      "${PROJECT_ROOT}/Canz.genome.fa" \
      "${PROJECT_ROOT}/01_At_TRM_query/Canz.genome.fa"; do
    [[ -f "$c" ]] && { printf '%s' "$c"; return 0; }
  done
  return 1
}

# Pfam-A.hmm 探测（含 hmmpress 索引检查）
probe_pfam() {
  if [[ -n "$PFAM_DB" && -f "$PFAM_DB" ]]; then printf '%s' "$PFAM_DB"; return 0; fi
  local c
  for c in \
      "${PFAM_DB:-}" \
      "$HOME/Pfam-A.hmm" \
      "/data/Pfam-A.hmm" \
      "/data/db/Pfam-A.hmm" \
      "/data2/db/Pfam-A.hmm" \
      "/opt/Pfam-A.hmm" \
      "/usr/local/share/Pfam-A.hmm"; do
    [[ -n "$c" && -f "$c" ]] && { printf '%s' "$c"; return 0; }
  done
  # 最后尝试在常见数据盘浅层搜索
  local hit
  hit="$(find /data /data2 /opt /home -maxdepth 4 -name 'Pfam-A.hmm' -type f 2>/dev/null | head -1 || true)"
  [[ -n "$hit" ]] && { printf '%s' "$hit"; return 0; }
  return 1
}
