#!/usr/bin/env bash

# If invoked via sh, restart with bash so bash-only syntax works.
if [ -z "${BASH_VERSION:-}" ]; then
  exec bash "$0" "$@"
fi

set -euo pipefail

# -------------------------------
# Simple configuration
# -------------------------------
SIM_FILE="sim.c"
SIM_BIN="./sim"
SIM_ARGS="-s0"
PYTHON_CMD="python3"
OPT_SCRIPT="optimize_source_params.py"
HISTORY_FILE="source_optimization_history.csv"

M0_TARGET="3.172138e+20"
M2_TARGET="7.682495e+32"
TOL_REL="0.000001"
MAX_ITERS=10

MIN_AMP="1e-30"
MIN_TEMP_EV="1e-6"

MC2P_FILE="gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl"
M0_FILE="gk_lorentzian_mirror-ion_source_M0_0.gkyl"
M2_FILE="gk_lorentzian_mirror-ion_source_M2_0.gkyl"
Z_RANGE="-0.98:0.98"

# Source line update controls.
# Defaults below target:
#   double gamma0 = ...;
#   double E_beam = ... * GKYL_ELEMENTARY_CHARGE;
# - Set *_LINE_NUM to a positive integer to target an exact line in SIM_FILE.
# - Leave *_LINE_NUM empty to target the first line matching *_LINE_MATCH regex.
# - *_LINE_FORMAT receives two printf placeholders: coeff then target.
AMP_LINE_NUM=""
TEMP_LINE_NUM=""
AMP_LINE_MATCH='^[[:space:]]*double[[:space:]]+gamma0[[:space:]]*='
TEMP_LINE_MATCH='^[[:space:]]*double[[:space:]]+E_beam[[:space:]]*='
AMP_LINE_FORMAT='double gamma0 = %s; // Beam intM0 = %s'
TEMP_LINE_FORMAT='double E_beam = %s * GKYL_ELEMENTARY_CHARGE; // Beam intM2 = %s'

usage() {
  cat <<'EOF'
Usage: ./optimize_source_params.sh [options]

Core options:
  --m0-target VALUE
  --m2-target VALUE
  --tol-rel VALUE
  --max-iters N

File/command options:
  --sim-file PATH
  --sim-args "ARGS"
  --python-cmd CMD
  --opt-script PATH
  --history-file PATH
  --mc2p FILE
  --z-range A:B
  -h, --help

Flow per iteration:
  1) run sim
  2) read M0, M2 with pgkyl
  3) append (amp,temp,M0,M2) to history csv
  4) ask Python optimizer for next (amp,temp)
  5) rewrite the two source lines in sim.c
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

is_numeric() {
  [[ "$1" =~ ^[+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?$ ]]
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --m0-target) M0_TARGET="$2"; shift 2 ;;
    --m2-target) M2_TARGET="$2"; shift 2 ;;
    --tol-rel) TOL_REL="$2"; shift 2 ;;
    --max-iters) MAX_ITERS="$2"; shift 2 ;;
    --python-cmd) PYTHON_CMD="$2"; shift 2 ;;
    --opt-script) OPT_SCRIPT="$2"; shift 2 ;;
    --history-file) HISTORY_FILE="$2"; shift 2 ;;
    --sim-file) SIM_FILE="$2"; shift 2 ;;
    --sim-args) SIM_ARGS="$2"; shift 2 ;;
    --z-range) Z_RANGE="$2"; shift 2 ;;
    --mc2p) MC2P_FILE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1" ;;
  esac
done

[[ "$MAX_ITERS" =~ ^[1-9][0-9]*$ ]] || die "max-iters must be a positive integer"
[[ -f "$SIM_FILE" ]] || die "sim file not found: $SIM_FILE"

need_cmd awk
need_cmd sed
need_cmd make
need_cmd pgkyl
need_cmd "$PYTHON_CMD"

[[ -f "$OPT_SCRIPT" ]] || die "optimizer script not found: $OPT_SCRIPT"

if [[ -n "$AMP_LINE_NUM" ]]; then
  [[ "$AMP_LINE_NUM" =~ ^[1-9][0-9]*$ ]] || die "AMP_LINE_NUM must be a positive integer or empty"
fi
if [[ -n "$TEMP_LINE_NUM" ]]; then
  [[ "$TEMP_LINE_NUM" =~ ^[1-9][0-9]*$ ]] || die "TEMP_LINE_NUM must be a positive integer or empty"
fi

get_source_line() {
  local kind="$1"
  local line_num match label line

  case "$kind" in
    amp)
      line_num="$AMP_LINE_NUM"
      match="$AMP_LINE_MATCH"
      label="gamma0"
      ;;
    temp)
      line_num="$TEMP_LINE_NUM"
      match="$TEMP_LINE_MATCH"
      label="E_beam"
      ;;
    *)
      die "unknown source kind: $kind"
      ;;
  esac

  if [[ -n "$line_num" ]]; then
    line=$(sed -n "${line_num}p" "$SIM_FILE")
    [[ -n "$line" ]] || die "line $line_num for $label does not exist in $SIM_FILE"
  else
    line=$(grep -E "$match" "$SIM_FILE" | head -n 1 || true)
    [[ -n "$line" ]] || die "could not find $label line in $SIM_FILE with regex: $match"
  fi

  echo "$line"
}

amp_line=$(get_source_line amp)
temp_line=$(get_source_line temp)

amp_coeff=$(echo "$amp_line" | sed -E 's/.*=[[:space:]]*([0-9.eE+-]+)[[:space:]]*;.*/\1/')
temp_coeff_ev=$(echo "$temp_line" | sed -E 's/.*=[[:space:]]*([0-9.eE+-]+)[[:space:]]*\*[[:space:]]*GKYL_ELEMENTARY_CHARGE.*/\1/')

[[ -n "$amp_coeff" ]] || die "failed to parse amplitude coefficient from sim.c"
[[ -n "$temp_coeff_ev" ]] || die "failed to parse temp coefficient from sim.c"

[[ "$M0_TARGET" =~ ^[0-9.eE+-]+$ ]] || die "m0-target must be numeric"
[[ "$M2_TARGET" =~ ^[0-9.eE+-]+$ ]] || die "m2-target must be numeric"
[[ "$TOL_REL" =~ ^[0-9.eE+-]+$ ]] || die "tol-rel must be numeric"

rel_err() {
  awk -v v="$1" -v t="$2" 'BEGIN { if (t==0) { print 0; exit } d=v-t; if (d<0) d=-d; print d/t }'
}

append_history() {
  local iter="$1" amp="$2" temp="$3" m0="$4" m2="$5"

  if [[ ! -f "$HISTORY_FILE" ]]; then
    echo "iter,amp,temp_eV,m0,m2" > "$HISTORY_FILE"
  fi

  echo "$iter,$amp,$temp,$m0,$m2" >> "$HISTORY_FILE"
}

suggest_next_inputs() {
  local current_amp="$1" current_temp="$2"
  local out s_amp s_temp s_method

  out=$("$PYTHON_CMD" "$OPT_SCRIPT" \
    --history "$HISTORY_FILE" \
    --m0-target "$M0_TARGET" \
    --m2-target "$M2_TARGET" \
    --current-amp "$current_amp" \
    --current-temp "$current_temp" \
    --min-amp "$MIN_AMP" \
    --min-temp "$MIN_TEMP_EV")

  read -r s_amp s_temp s_method <<< "$out"
  is_numeric "$s_amp" || die "optimizer returned non-numeric amp: '$s_amp'"
  is_numeric "$s_temp" || die "optimizer returned non-numeric temp: '$s_temp'"
  [[ -n "$s_method" ]] || s_method="unknown"
  echo "$s_amp $s_temp $s_method"
}

replace_source_line() {
  local kind="$1" coeff="$2" target="$3"
  local tmp_file line_num line_match fmt mode replacement
  tmp_file=$(mktemp)

  case "$kind" in
    amp)
      line_num="$AMP_LINE_NUM"
      line_match="$AMP_LINE_MATCH"
      fmt="$AMP_LINE_FORMAT"
      ;;
    temp)
      line_num="$TEMP_LINE_NUM"
      line_match="$TEMP_LINE_MATCH"
      fmt="$TEMP_LINE_FORMAT"
      ;;
    *)
      die "unknown source kind: $kind"
      ;;
  esac

  mode="regex"
  if [[ -n "$line_num" ]]; then
    mode="line"
  fi

  printf -v replacement "$fmt" "$coeff" "$target"

  awk -v mode="$mode" -v line_num="$line_num" -v line_match="$line_match" -v replacement="$replacement" '
    BEGIN { updated=0 }
    {
      if (!updated && ((mode=="line" && NR==line_num) || (mode=="regex" && $0 ~ line_match))) {
        print replacement;
        updated=1;
        next;
      }
      print $0
    }
    END {
      if (!updated) exit 7;
    }
  ' "$SIM_FILE" > "$tmp_file" || {
    rm -f "$tmp_file"
    if [[ "$mode" == "line" ]]; then
      die "failed to replace $kind line at line number $line_num in $SIM_FILE"
    fi
    die "failed to replace $kind line matching regex '$line_match' in $SIM_FILE"
  }

  mv "$tmp_file" "$SIM_FILE"
}

get_moment_max() {
  local moment_file="$1"
  local output value

  output=$(pgkyl --c2p "$MC2P_FILE" "$moment_file" interp sel --z0 "$Z_RANGE" integ 0 info)
  value=$(echo "$output" | awk '
    /Maximum:/ {
      for (i=1; i<=NF; i++) {
        if ($i ~ /^[+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?$/) { print $i; exit }
      }
    }
  ')

  [[ -n "$value" ]] || die "failed to parse Maximum from $moment_file"
  is_numeric "$value" || die "parsed non-numeric Maximum from $moment_file: '$value'"
  echo "$value"
}

run_sim() {
  echo "[run] make && $SIM_BIN $SIM_ARGS"
  make > /dev/null 2>&1
  # shellcheck disable=SC2086
  $SIM_BIN $SIM_ARGS > /dev/null 2>&1
}

echo "Initial coefficients: amplitude = $amp_coeff, temp_eV = $temp_coeff_ev"
echo "Targets: M0 = $M0_TARGET, M2 = $M2_TARGET"
echo "Stop when both relative errors are <= $TOL_REL, max_iters = $MAX_ITERS"
echo "Optimizer: $PYTHON_CMD $OPT_SCRIPT"
echo "History file: $HISTORY_FILE"

iter=0
err_m0=1
err_m2=1
while [ "$iter" -lt "$MAX_ITERS" ]; do
  iter=$((iter + 1))
  echo ""
  echo "=== Iteration $iter/$MAX_ITERS ==="
  echo "[state] amp=$amp_coeff temp_eV=$temp_coeff_ev"

  run_sim

  m0=$(get_moment_max "$M0_FILE")
  m2=$(get_moment_max "$M2_FILE")
  err_m0=$(rel_err "$m0" "$M0_TARGET")
  err_m2=$(rel_err "$m2" "$M2_TARGET")

  echo "[diag] M0=$m0 (rel_err=$err_m0)"
  echo "[diag] M2=$m2 (rel_err=$err_m2)"

  append_history "$iter" "$amp_coeff" "$temp_coeff_ev" "$m0" "$m2"

  if awk -v e1="$err_m0" -v e2="$err_m2" -v tol="$TOL_REL" 'BEGIN {exit !((e1<=tol) && (e2<=tol))}'; then
    echo "Converged."
    break
  fi

  read -r next_amp next_temp next_method <<< "$(suggest_next_inputs "$amp_coeff" "$temp_coeff_ev")"
  echo "[opt] $next_method -> amp=$next_amp temp_eV=$next_temp"

  amp_coeff="$next_amp"
  temp_coeff_ev="$next_temp"
  replace_source_line amp "$amp_coeff" "$M0_TARGET"
  replace_source_line temp "$temp_coeff_ev" "$M2_TARGET"
done

if [ "$iter" -ge "$MAX_ITERS" ]; then
  echo "Reached max iterations ($MAX_ITERS)."
fi

echo ""
echo "Done. Final source lines in $SIM_FILE:"
grep -E '^[[:space:]]*double[[:space:]]+ion_source_(amplitude|temp)[[:space:]]*=' "$SIM_FILE"
