#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  printf 'usage: %s ABSOLUTE_OUTPUT_ROOT [WORKERS]\n' "$0" >&2
  exit 64
fi

OUTPUT_ROOT=$(/usr/bin/realpath -m -- "$1")
WORKERS=${2:-4}
readonly A9B_SYSTEMD_UNIT_NAME=a9b-sensitivity-campaign
readonly A9B_RUNTIME_MAX_SEC=86400
export A9B_RUNTIME_MAX_SEC
SCRIPT_DIR=$(cd -- "$(/usr/bin/dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
SCRIPT_PATH="$SCRIPT_DIR/$(/usr/bin/basename -- "${BASH_SOURCE[0]}")"
PETRA_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
REPO_ROOT=$(cd -- "$PETRA_ROOT/.." && pwd)
CARGO_TARGET_DIR=/home/vsletten/.hermes/cache/cargo-target/vsletten-dissertation-da68ebaefa3a
export CARGO_TARGET_DIR

if [[ "$OUTPUT_ROOT" != /* ]]; then
  printf 'output root must be absolute\n' >&2
  exit 64
fi
case "$OUTPUT_ROOT/" in
  "$REPO_ROOT"/*)
    printf 'output root must be outside the git worktree: %s\n' "$OUTPUT_ROOT" >&2
    exit 64
    ;;
esac
if [[ ! "$WORKERS" =~ ^[1-4]$ ]]; then
  printf 'workers must be an integer in [1,4]\n' >&2
  exit 64
fi

# Only a live cgroup + D-Bus proof of exact unit membership and RuntimeMaxUSec
# can enter the workload. Caller-controlled environment markers are ignored.
SYSTEMD_PROOF=$(
  /usr/bin/python3 -I "$SCRIPT_DIR/a9b_bounded_launch_guard.py" \
    "$A9B_SYSTEMD_UNIT_NAME.service" "$A9B_RUNTIME_MAX_SEC" 2>/dev/null
) || {
  unset A9B_SYSTEMD_INNER A9B_SYSTEMD_UNIT INVOCATION_ID
  exec systemd-run --user --unit="$A9B_SYSTEMD_UNIT_NAME" --collect --wait --pipe \
    --property="RuntimeMaxSec=$A9B_RUNTIME_MAX_SEC" \
    /usr/bin/nice -n 10 "$SCRIPT_PATH" "$@"
}
read -r A9B_LIVE_RUNTIME_MAX_USEC A9B_LIVE_INVOCATION_ID A9B_LIVE_NICENESS <<<"$SYSTEMD_PROOF"
if [[ ! $A9B_LIVE_RUNTIME_MAX_USEC =~ ^[1-9][0-9]*$ || ! $A9B_LIVE_INVOCATION_ID =~ ^[0-9a-f]{32}$ || ! $A9B_LIVE_NICENESS =~ ^-?[0-9]+$ ]]; then
  printf 'refusing workload: malformed live systemd proof\n' >&2
  exit 77
fi

# Re-exec if needed, then require a fresh /proc attestation before workload.
if (( A9B_LIVE_NICENESS < 10 )); then
  exec /usr/bin/nice -n 10 "$SCRIPT_PATH" "$@"
fi
SYSTEMD_PROOF=$(
  /usr/bin/python3 -I "$SCRIPT_DIR/a9b_bounded_launch_guard.py" \
    "$A9B_SYSTEMD_UNIT_NAME.service" "$A9B_RUNTIME_MAX_SEC" 10
) || {
  printf 'refusing workload: live niceness attestation failed\n' >&2
  exit 77
}
read -r A9B_LIVE_RUNTIME_MAX_USEC A9B_LIVE_INVOCATION_ID A9B_LIVE_NICENESS <<<"$SYSTEMD_PROOF"
if [[ ! $A9B_LIVE_RUNTIME_MAX_USEC =~ ^[1-9][0-9]*$ || ! $A9B_LIVE_INVOCATION_ID =~ ^[0-9a-f]{32}$ || ! $A9B_LIVE_NICENESS =~ ^-?[0-9]+$ ]] || (( A9B_LIVE_NICENESS < 10 )); then
  printf 'refusing workload: malformed live niceness attestation\n' >&2
  exit 77
fi
readonly A9B_LIVE_RUNTIME_MAX_USEC A9B_LIVE_INVOCATION_ID A9B_LIVE_NICENESS
export A9B_LIVE_RUNTIME_MAX_USEC A9B_LIVE_INVOCATION_ID A9B_LIVE_NICENESS

THREADS=$((16 / WORKERS))
export OMP_NUM_THREADS=$THREADS
export OPENBLAS_NUM_THREADS=$THREADS
export MKL_NUM_THREADS=$THREADS
export RAYON_NUM_THREADS=$THREADS
export NUMEXPR_NUM_THREADS=$THREADS

OPERATOR_ROOT="${OUTPUT_ROOT}-operator"
mkdir -p -- "$OPERATOR_ROOT"
ATTEMPT=1
while [[ -e "$OPERATOR_ROOT/launch-attempt-${ATTEMPT}.log" || -e "$OPERATOR_ROOT/launch-attempt-${ATTEMPT}.json" ]]; do
  ATTEMPT=$((ATTEMPT + 1))
done
LOG="$OPERATOR_ROOT/launch-attempt-${ATTEMPT}.log"
RECEIPT="$OPERATOR_ROOT/launch-attempt-${ATTEMPT}.json"
STATUS=running
STARTED="$(date --iso-8601=seconds)"

finish() {
  local code=$?
  trap - EXIT
  local ended
  ended="$(date --iso-8601=seconds)"
  if [[ $code -eq 0 ]]; then STATUS=complete; else STATUS=failed; fi
  python3 -c 'import hashlib,json,os,pathlib,sys,tempfile
receipt=pathlib.Path(sys.argv[1]); log=pathlib.Path(sys.argv[2]); root=pathlib.Path(sys.argv[3]); status=sys.argv[4]; code=int(sys.argv[5]); attempt=int(sys.argv[6]); workers=int(sys.argv[7]); threads=int(sys.argv[8]); started=sys.argv[9]; ended=sys.argv[10]
def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
runtime_max_usec=int(os.environ["A9B_LIVE_RUNTIME_MAX_USEC"])
niceness=int(os.environ["A9B_LIVE_NICENESS"])
payload={"schema":"a9b-campaign-launch-receipt-v2","status":status,"exit_code":code,"attempt":attempt,"started_at":started,"ended_at":ended,"workers":workers,"threads_per_worker":threads,"niceness":niceness,"systemd_unit":"a9b-sensitivity-campaign.service","systemd_invocation_id":os.environ["A9B_LIVE_INVOCATION_ID"],"runtime_max_usec":runtime_max_usec,"runtime_max_sec":runtime_max_usec/1000000,"output_root":str(root),"log":str(log),"log_sha256":digest(log),"manifest_sha256":digest(root/"manifest.json"),"analysis_sha256":digest(root/"analysis.json"),"verification_sha256":digest(root/"verification.json")}
receipt.parent.mkdir(parents=True,exist_ok=True)
fd,tmp=tempfile.mkstemp(prefix="."+receipt.name+".",dir=receipt.parent)
try:
  with os.fdopen(fd,"w") as f: json.dump(payload,f,indent=2,sort_keys=True,allow_nan=False); f.write("\n"); f.flush(); os.fsync(f.fileno())
  os.link(tmp,receipt); os.unlink(tmp)
  directory=os.open(receipt.parent,os.O_RDONLY|os.O_DIRECTORY)
  try: os.fsync(directory)
  finally: os.close(directory)
except BaseException:
  try: os.unlink(tmp)
  except FileNotFoundError: pass
  raise' "$RECEIPT" "$LOG" "$OUTPUT_ROOT" "$STATUS" "$code" "$ATTEMPT" "$WORKERS" "$THREADS" "$STARTED" "$ended"
  exit "$code"
}
trap finish EXIT

# Fail closed on the v16 workstation stability gate. Sample swap-out growth
# over two seconds rather than treating a static nonzero swap total as growth.
read -r _ _ LOAD15 _ < /proc/loadavg
MEM_KB=$(python3 -c 'import pathlib; rows=dict(line.split(":",1) for line in pathlib.Path("/proc/meminfo").read_text().splitlines()); print(rows["MemAvailable"].split()[0])')
SWAP_BEFORE=$(python3 -c 'import pathlib; print(next(line.split()[1] for line in pathlib.Path("/proc/vmstat").read_text().splitlines() if line.startswith("pswpout ")))')
sleep 2
SWAP_AFTER=$(python3 -c 'import pathlib; print(next(line.split()[1] for line in pathlib.Path("/proc/vmstat").read_text().splitlines() if line.startswith("pswpout ")))')
python3 -c 'import sys; load=float(sys.argv[1]); mem=int(sys.argv[2]); before=int(sys.argv[3]); after=int(sys.argv[4]);
if load > 24: raise SystemExit(f"15-minute loadavg {load} exceeds 24")
if mem < 12*1024*1024: raise SystemExit(f"MemAvailable {mem} kB is below 12 GiB")
if after > before: raise SystemExit(f"swap-out grew from {before} to {after}")' "$LOAD15" "$MEM_KB" "$SWAP_BEFORE" "$SWAP_AFTER"

{
  printf 'A9b survey-tier campaign; never production\n'
  printf 'output_root=%s workers=%s threads_per_worker=%s niceness=%s\n' "$OUTPUT_ROOT" "$WORKERS" "$THREADS" "$A9B_LIVE_NICENESS"
  printf 'cargo_target=%s\n' "$CARGO_TARGET_DIR"
  cargo build --release -p petra-deck --example a9b_importance_sampling --manifest-path "$PETRA_ROOT/Cargo.toml"
  RUNNER="$CARGO_TARGET_DIR/release/examples/a9b_importance_sampling"
  python3 "$SCRIPT_DIR/a9b_sensitivity_ranking.py" validate "$OUTPUT_ROOT"
  python3 "$SCRIPT_DIR/a9b_sensitivity_ranking.py" run "$OUTPUT_ROOT" --runner "$RUNNER" --workers "$WORKERS"
  python3 "$SCRIPT_DIR/a9b_sensitivity_ranking.py" analyze "$OUTPUT_ROOT"
  python3 "$SCRIPT_DIR/verify_a9b_sensitivity_ranking.py" "$OUTPUT_ROOT" --write-verification
} 2>&1 | tee "$LOG"
sync -f "$LOG"
