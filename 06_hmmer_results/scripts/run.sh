#!/usr/bin/env bash
# =============================================================================
# run.sh —— CaTRM 家族鉴定流程主入口（分阶段执行）
#
# 用法:
#   bash run.sh              # 列出全部阶段与状态
#   bash run.sh 01           # 只跑阶段 01（环境与输入文件检查）
#   bash run.sh 02           # 只跑阶段 02（MAFFT 比对）
#   bash run.sh 01 02 03     # 按顺序跑多个阶段
#   bash run.sh all          # 依次跑 01→07
#
# 阶段:
#   01  环境与输入文件检查        -> 01_env_check/
#   02  MAFFT 比对 34 条 AtTRM    -> 02_msa/
#   03  hmmbuild 建 TRM.hmm       -> 03_hmmbuild/
#   04  hmmsearch 扫描 Zunla 全蛋白 -> 04_hmmsearch/
#   05  与 BLAST 候选合并去重     -> 05_merge/
#   06  Pfam 结构域核查           -> 06_domain/
#   07  最终 CaTRM 列表           -> 07_final/
#
# 额外参数会透传给阶段脚本（如 --gff3 / --genome / --pfam / --threads）
# =============================================================================
set -euo pipefail

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPTS_DIR}/common.sh"

STAGES=(
  "01:stage_01_env_check.sh:环境与输入文件检查"
  "02:stage_02_msa.sh:MAFFT 比对 34 条 AtTRM"
  "03:stage_03_hmmbuild.sh:hmmbuild 建立 TRM.hmm"
  "04:stage_04_hmmsearch.sh:hmmsearch 扫描 Zunla 全蛋白"
  "05:stage_05_merge.sh:与 BLAST 候选合并去重"
  "06:stage_06_domain.sh:Pfam 结构域核查"
  "07:stage_07_final.sh:最终 CaTRM 基因/蛋白列表"
)

stage_desc() {
  local want="$1" entry
  for entry in "${STAGES[@]}"; do
    if [[ "${entry%%:*}" == "$want" ]]; then
      printf '%s' "${entry#*:}"
      return 0
    fi
  done
  return 1
}
stage_script() {
  local want="$1" entry
  for entry in "${STAGES[@]}"; do
    if [[ "${entry%%:*}" == "$want" ]]; then
      local rest="${entry#*:}"
      printf '%s' "${rest%%:*}"
      return 0
    fi
  done
  return 1
}

list_stages() {
  echo "CaTRM 家族鉴定流程 —— 阶段一览"
  echo "======================================================================"
  local entry id script desc outdir n
  for entry in "${STAGES[@]}"; do
    id="${entry%%:*}"
    local rest="${entry#*:}"
    script="${rest%%:*}"
    desc="${rest#*:}"
    case "$id" in
      01) outdir="$D01" ;; 02) outdir="$D02" ;; 03) outdir="$D03" ;;
      04) outdir="$D04" ;; 05) outdir="$D05" ;; 06) outdir="$D06" ;;
      07) outdir="$D07" ;;
    esac
    if [[ -d "$outdir" ]]; then
      n="$(find "$outdir" -maxdepth 1 -type f 2>/dev/null | wc -l | tr -d ' ')"
    else
      n=0
    fi
    if (( n > 0 )); then
      printf '  %s  %-28s %s\n      -> %s  (%s 个文件)\n' \
        "$id" "$desc" "[已运行]" "${outdir#"$HMMER_DIR"/}" "$n"
    else
      printf '  %s  %-28s %s\n      -> %s\n' \
        "$id" "$desc" "[未运行]" "${outdir#"$HMMER_DIR"/}"
    fi
  done
  echo "======================================================================"
  echo "运行方式:  bash run.sh <阶段号> [阶段号...]     例: bash run.sh 01"
  echo "透传参数:  --gff3 <file> --genome <file> --pfam <Pfam-A.hmm> --threads N"
}

run_stage() {
  local id="$1"; shift
  local script
  script="$(stage_script "$id")" || die "未知阶段号: $id (可选 01-07 或 all)"
  local path="${SCRIPTS_DIR}/${script}"
  [[ -f "$path" ]] || die "阶段脚本不存在: $path"

  printf '\n'
  printf '######################################################################\n'
  printf '#  阶段 %s: %s\n' "$id" "$(stage_desc "$id")"
  printf '#  脚本: %s\n' "${path#"$SCRIPTS_DIR"/}"
  printf '######################################################################\n'
  bash "$path" "$@"
}

# ------------------------------ 主流程 --------------------------------------
if [[ $# -eq 0 ]]; then
  list_stages
  exit 0
fi

# 分离阶段号与透传参数
STAGE_ARGS=()
PASS_ARGS=()
for a in "$@"; do
  if [[ "$a" == "all" ]]; then
    STAGE_ARGS+=(01 02 03 04 05 06 07)
  elif [[ "$a" =~ ^0[1-7]$ ]]; then
    STAGE_ARGS+=("$a")
  else
    PASS_ARGS+=("$a")
  fi
done

[[ ${#STAGE_ARGS[@]} -gt 0 ]] || die "没有指定有效阶段号（01-07 或 all）"

for id in "${STAGE_ARGS[@]}"; do
  run_stage "$id" "${PASS_ARGS[@]}"
done

printf '\n全部指定阶段执行完毕。\n'
