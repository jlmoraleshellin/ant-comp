#!/usr/bin/env bash
# scripts/vib.sh — the ONE way this repo talks to the vib cluster.
#
#   bash scripts/vib.sh help
#
# Every call: reads the three personal values from .env (host, workspace, account),
# validates them, wraps ssh in `timeout`, starts in the workspace with the venv active,
# and refuses commands that would leave it. Nothing else should call ssh/rsync for vib.
#
# The guard is an ACCIDENT guard, not a sandbox: it catches `cd`, `..`, ~, absolute paths
# outside the workspace and an `rm` outside it. Code running inside an `srun` allocation is
# still your own user. The hard limit is your Unix permissions on group storage.
#
# Needs locally: bash, ssh, rsync, timeout, python3 (guard only). No uv locally; `bootstrap`
# builds the venv on vib with `module load uv`.

set -u -o pipefail

SSH_OPTS=(-n -x -o BatchMode=yes -o ConnectTimeout=20)
RSYNC_SSH="ssh -x -o BatchMode=yes -o ConnectTimeout=20"
MARKER=".sapia-workspace"
# What the cluster runs, and nothing else. `outputs/` and `.venv/` are not listed, so the
# --delete below can never reach them. `--delete` only acts inside the listed directories.
SYNC_PATHS=(tools activation utils .env)

die() { echo "vib: $*" >&2; exit 2; }

# POSIX single-quote a string so it survives one more shell layer unchanged.
q() { local s=${1//\'/\'\\\'\'}; printf "'%s'" "$s"; }

usage() {
  cat <<'EOF'
usage: bash scripts/vib.sh <command> [args]

  exec [--sync] [--timeout S] '<command>'
        Run <command> in the workspace, venv active, .env loaded. --sync pushes the
        working tree first (use it right before `sapia run`; not before polls).
        Default timeout 120 s. $SAPIA_VIB_WORKSPACE / $SAPIA_VIB_ACCOUNT expand ON vib,
        so write them in single quotes.
  exec --on-node [--cpus N] [--mem SIZE] [--walltime H:MM:SS] [...as above] '<command>'
        The same, but inside a small SLURM allocation instead of the login node: use it for
        `sapia collect`, pandas, anything beyond a scheduler query. Defaults: 2 cpus, 8G,
        0:30:00, ssh timeout 1800 s. --cpus is capped at 8; heavier work: write the srun yourself.
        Not needed for `squeue`/`sacct`/`sinfo` polls.
  sync [--dry-run]     push tools/ activation/ utils/ .env to the workspace
  check '<command>'    run the guard only; touches nothing
  init [--adopt]       create the workspace + its marker; --adopt for an existing one
  bootstrap [--commit SHA]
                       init + sync + build the venv with `uv` on a compute node (~5 min).
                       Default: `uv sync --frozen` from uv.lock, so vib runs the exact set Modal
                       does. --commit SHA: install that pushed prosapia commit instead (needed
                       while uv.lock points at the editable ../prosapia).
  status               marker, venv, prosapia commit vs the pin, and `module load uv`

Set VIB_PRINT=1 to print the ssh/rsync command instead of running it.
Set VIB_DEBUG=1 on `sync` for rsync -vv --stats and ssh -v (long: redirect stderr to a file).
Exit 3 = the guard refused. ssh rc 255 or 124 usually means the certificate lapsed.
EOF
}

# ── config ──────────────────────────────────────────────────────────────────
load_config() {
  cd "$(dirname "${BASH_SOURCE[0]}")/.." || die "cannot cd to the repo root"
  [ -f .env ] || die "no .env in $(pwd)"
  local b; for b in ssh rsync timeout python3; do
    command -v "$b" >/dev/null || die "'$b' is required locally"
  done
  set +u; set -a; . ./.env; set +a; set -u
  HOST=${SAPIA_VIB_HOST:-}
  WS=${SAPIA_VIB_WORKSPACE:-}
  ACCOUNT=${SAPIA_VIB_ACCOUNT:-}
  local re='^[A-Za-z0-9._-]+$'
  [[ $HOST =~ $re ]] || die "SAPIA_VIB_HOST must be an ssh alias ([A-Za-z0-9._-]), got '$HOST'"
  WS=${WS%/}
  re='^(/[A-Za-z0-9._@+-]+){4,}$'
  [[ $WS =~ $re ]] || die "SAPIA_VIB_WORKSPACE must be an absolute path of at least 4 components ([A-Za-z0-9._@+-] only), got '$WS'"
  case "$WS/" in *"/../"*) die "SAPIA_VIB_WORKSPACE must not contain '..'" ;; esac
}

need_account() {
  local re='^[A-Za-z0-9._-]+$'
  [[ $ACCOUNT =~ $re ]] || die "SAPIA_VIB_ACCOUNT is empty or malformed ('$ACCOUNT'); ask the user, never guess"
}

# ── ssh ─────────────────────────────────────────────────────────────────────
# run_ssh <timeout-s> <remote command string>
run_ssh() {
  local t=$1 remote=$2 err rc
  if [ -n "${VIB_PRINT:-}" ]; then
    printf 'ssh %s %s %s\n' "${SSH_OPTS[*]}" "$HOST" "$remote"
    return 0
  fi
  err=$(mktemp)
  timeout -k 5 "$t" ssh "${SSH_OPTS[@]}" "$HOST" "$remote" 2>"$err" && rc=0 || rc=$?
  grep -v -e Provisioner -e X11 -e xauth "$err" >&2 || true
  rm -f "$err"
  case $rc in
    255|124|137)
      echo "vib: ssh did not complete (rc=$rc). Most likely the certificate lapsed: ask the user to sign in (\`! ssh $HOST true\`). Do not retry in a loop." >&2 ;;
  esac
  return $rc
}

# remote <command> <timeout-s>: the command, run from the workspace with .env + venv.
remote() {
  local script
  script="cd $(q "$WS") || { echo 'vib: workspace directory is missing; run: bash scripts/vib.sh init' >&2; exit 97; }
set -a; [ -f .env ] && . ./.env; set +a
if [ -f .venv/bin/activate ]; then . .venv/bin/activate; else echo 'vib: no .venv in the workspace; run: bash scripts/vib.sh bootstrap' >&2; fi
$1"
  run_ssh "$2" "bash -c $(q "$script")"
}

# ── the guard ───────────────────────────────────────────────────────────────
IFS= read -r -d '' GUARD_PY <<'PY' || true
import re, shlex, sys

ws, cmd = sys.argv[1], sys.argv[2]
WSV = r"\$SAPIA_VIB_WORKSPACE|\$\{SAPIA_VIB_WORKSPACE\}"
bad = []

if not cmd.strip():
    bad.append("empty command")

# every command already starts in the workspace; there is nothing to cd to
if re.search(r"(?<![\w./$-])(?:cd|pushd|popd)(?![\w.-])", cmd) or "chdir" in cmd:
    bad.append("directory change (cd / pushd / popd / chdir): every command already starts in the workspace")

if re.search(r"""(?:^|[\s/'"=:(,;&|<>])\.\.(?:[/\s'"),;&|<>]|$)""", cmd):
    bad.append("'..' path component")

if re.search(r"""(?:^|[\s'"=:(])~""", cmd) or re.search(r"\$\{?(?:HOME|OLDPWD)\b", cmd):
    bad.append("home-relative path (~, $HOME, $OLDPWD)")

if re.search(r"(?<![\w./$-])(?:ssh|scp|sftp|rsync)(?![\w.-])", cmd):
    bad.append("ssh / scp / sftp / rsync from the cluster")

# the filesystem root as a bare argument (`find / ...`, `chmod -R 777 /`): PATH_RE below needs a
# name after the slash, so a lone `/` is caught here, only after a path-taking command
if re.search(r"(?<![\w./$-])(?:find|ls|du|chmod|chown|chgrp|tree|stat|tar|cat|cp|mv|ln|grep)(?:\s+[^\s;&|]+){0,3}?\s+/(?=[\s;&|)'\"]|$)", cmd):
    bad.append("the filesystem root '/' as an argument")

# absolute paths must be inside the workspace. A path reached THROUGH a variable
# ($SAPIA_VIB_WORKSPACE/..., $SAPIA_VIB_SHARED_ACTIVATION/...) is not seen here by design.
ALLOW = {"/dev/null", "/dev/stdin", "/dev/stdout", "/dev/stderr"}
PATH_RE = re.compile(r"(?<![\w.:/~$}\])*?\[{-])(/[\w.@%+=,*?-]+(?:/[\w.@%+=,*?-]*)*)")
for m in PATH_RE.finditer(cmd):
    p = m.group(1)
    if p in ALLOW or p == ws or p.startswith(ws + "/"):
        continue
    bad.append("absolute path outside the workspace: " + p)

# rm may delete anything INSIDE the workspace (outputs/ and .venv/ included); it may not
# reach outside it, nor remove the workspace root itself. `..` and absolute paths are already
# refused above, so a relative target is always inside.
WSRE = re.compile(r"^(?:" + WSV + "|" + re.escape(ws) + r")(?=/|$)")
for m in re.finditer(r"(?<![\w./$-])rm(?=\s)([^;&|\n)]*)", cmd):
    toks = None
    for seg in (m.group(1), m.group(1).rstrip("'\" ")):  # raw first; the 2nd drops a closing quote of an enclosing bash -lc '...'
        try:
            toks = shlex.split(seg)
            break
        except ValueError:
            continue
    if toks is None:
        bad.append("cannot parse the arguments of rm")
        continue
    targets = [t for t in toks if not t.startswith("-")]
    if not targets:
        bad.append("rm with no parsable target")
    for t in targets:
        wm = WSRE.match(t)
        if wm:
            rel = t[wm.end():]
        elif t.startswith("/"):
            bad.append("rm target outside the workspace: " + t)
            continue
        elif t.startswith(("$", "`")):
            bad.append("rm target behind a variable the guard cannot resolve: " + t)
            continue
        else:
            rel = t
        if rel.strip("/") in ("", "."):
            bad.append("rm of the workspace root itself: " + t)

if bad:
    for b in bad:
        sys.stderr.write("vib guard: " + b + "\n")
    sys.stderr.write("vib guard: refused (workspace " + ws + "). An accident guard, not a sandbox.\n")
    sys.exit(3)
PY

# guard <command text>: silent on pass, rc 3 and reasons on stderr on refusal
guard() { python3 -I -c "$GUARD_PY" "$WS" "$1"; }

# ── commands ────────────────────────────────────────────────────────────────
require_marker() {
  run_ssh 60 "test -f $(q "$WS/$MARKER")"
  case $? in
    0) return 0 ;;
    1) echo "vib: $WS has no $MARKER marker, so sync (which deletes) refuses to touch it. New workspace: bash scripts/vib.sh bootstrap. Existing and yours: bash scripts/vib.sh init --adopt" >&2
       return 1 ;;
    *) return 1 ;;
  esac
}

cmd_sync() {
  local extra=() paths=() p rc
  [ "${1:-}" = "--dry-run" ] && extra=(-n -i)
  for p in "${SYNC_PATHS[@]}"; do
    if [ -e "$p" ]; then paths+=("$p"); else echo "vib: skipping missing local path: $p" >&2; fi
  done
  require_marker || return 1
  # --timeout bounds I/O inactivity, so a stuck connection fails in ~90 s instead of eating the
  # whole 300 s. VIB_DEBUG=1 adds rsync -vv --stats and ssh -v (long; redirect stderr to a file).
  local ssh_cmd=$RSYNC_SSH dbg=()
  if [ -n "${VIB_DEBUG:-}" ]; then dbg=(-vv --stats); ssh_cmd="$RSYNC_SSH -v"; fi
  if [ -n "${VIB_PRINT:-}" ]; then
    echo "rsync -az --delete --relative --timeout=90 ${extra[*]:-} ${dbg[*]:-} -e '$ssh_cmd' ${paths[*]} $HOST:$WS/"
    return 0
  fi
  timeout -k 5 300 rsync -az --delete --relative --timeout=90 ${extra[@]+"${extra[@]}"} ${dbg[@]+"${dbg[@]}"} \
    --exclude='__pycache__' --exclude='*.pyc' --exclude='.ruff_cache' --exclude='results_shard_0' \
    -e "$ssh_cmd" "${paths[@]}" "$HOST:$WS/" && rc=0 || rc=$?
  case $rc in
    0) ;;
    255|124|137) echo "vib: rsync did not complete (rc=$rc). Certificate lapsed? Ask the user to sign in; do not retry in a loop." >&2 ;;
    *) echo "vib: rsync failed (rc=$rc)" >&2 ;;
  esac
  return $rc
}

cmd_exec() {
  local do_sync=0 t="" cmd on_node=0 cpus=2 mem=8G wall=0:30:00 re
  while [ $# -gt 0 ]; do
    case $1 in
      --sync) do_sync=1; shift ;;
      --timeout) t=${2:-}; shift 2 || die "--timeout needs seconds" ;;
      --on-node) on_node=1; shift ;;
      --cpus) cpus=${2:-}; shift 2 || die "--cpus needs a number" ;;
      --mem) mem=${2:-}; shift 2 || die "--mem needs a size, e.g. 8G" ;;
      --walltime) wall=${2:-}; shift 2 || die "--walltime needs H:MM:SS" ;;
      --) shift; break ;;
      *) break ;;
    esac
  done
  # a call on a compute node can wait in the queue first, so it gets the allocation's length
  if [ -z "$t" ]; then if [ "$on_node" = 1 ]; then t=1800; else t=120; fi; fi
  [[ $t =~ ^[0-9]+$ ]] || die "--timeout must be whole seconds, got '$t'"
  [ $# -ge 1 ] || die "exec needs a command (see: bash scripts/vib.sh help)"
  cmd="$*"
  guard "$cmd" || exit 3
  if [ "$on_node" = 1 ]; then
    need_account
    re='^[0-9]+$'; [[ $cpus =~ $re ]] && [ "$cpus" -ge 1 ] && [ "$cpus" -le 8 ] \
      || die "--cpus must be 1-8 (anything heavier: write the srun yourself), got '$cpus'"
    re='^[0-9]+[KMGT]?$'; [[ $mem =~ $re ]] || die "--mem must look like 8G, got '$mem'"
    re='^[0-9]+:[0-9]{2}:[0-9]{2}$'; [[ $wall =~ $re ]] || die "--walltime must be H:MM:SS, got '$wall'"
    # The step starts in the workspace (srun keeps the cwd). The venv is re-activated because
    # `bash -l` re-reads the login profile. The job's own SLURM_* variables are dropped so an
    # `sbatch` run from here (sapia run) cannot inherit this small allocation's cpus / memory.
    cmd=". .venv/bin/activate; unset \$(compgen -v SLURM_); $cmd"
    cmd="srun -A $(q "$ACCOUNT") -c $cpus --mem=$mem -t $wall --job-name=sapia_node bash -lc $(q "$cmd")"
  fi
  if [ "$do_sync" = 1 ]; then cmd_sync || exit $?; fi
  remote "$cmd" "$t"
}

cmd_init() {
  local adopt=0 script
  [ "${1:-}" = "--adopt" ] && adopt=1
  IFS= read -r -d '' script <<'EOF' || true
set -e
WS=@WS@
if [ -f "$WS/.sapia-workspace" ]; then echo "marker already present in $WS"; exit 0; fi
if [ -d "$WS" ] && [ -n "$(ls -A "$WS")" ] && [ @ADOPT@ != 1 ]; then
  echo "vib: $WS exists, is not empty and has no .sapia-workspace marker; refusing. If it IS your workspace: bash scripts/vib.sh init --adopt" >&2
  exit 2
fi
mkdir -p "$WS"
echo "created by scripts/vib.sh init" > "$WS/.sapia-workspace"
echo "marker created in $WS"
EOF
  script=${script//@WS@/$(q "$WS")}
  script=${script//@ADOPT@/$adopt}
  run_ssh 60 "bash -c $(q "$script")"
}

# The prosapia commit the venv must carry: --commit, else the git source in uv.lock.
resolve_commit() {
  local c=${1:-}
  if [ -z "$c" ] && [ -f uv.lock ]; then
    c=$(grep -oE 'prosapia\?branch=[A-Za-z0-9._-]+#[0-9a-f]{40}' uv.lock | head -1 | cut -d'#' -f2)
  fi
  [ -n "$c" ] || die "no prosapia commit to pin: uv.lock has no git source for prosapia (it is probably set to the editable ../prosapia checkout). Restore the git source, or pass --commit <sha of a PUSHED prosapia commit>."
  local re='^[0-9a-f]{7,40}$'
  [[ $c =~ $re ]] || die "not a commit sha: '$c'"
  echo "$c"
}

cmd_bootstrap() {
  local commit_arg="" commit script build
  while [ $# -gt 0 ]; do
    case $1 in
      --commit) commit_arg=${2:-}; shift 2 || die "--commit needs a sha" ;;
      *) die "unknown bootstrap option '$1'" ;;
    esac
  done
  need_account
  commit=$(resolve_commit "$commit_arg") || exit 2
  # `uv sync` reads pyproject.toml + uv.lock from the workspace. A per-run sync never needs
  # them, so only bootstrap sends them. (--delete cannot touch a file path.)
  local SYNC_PATHS=("${SYNC_PATHS[@]}")
  if [ -z "$commit_arg" ]; then
    [ -f pyproject.toml ] && [ -f uv.lock ] || die "pyproject.toml and uv.lock must both exist here: uv sync reads them"
    SYNC_PATHS+=(pyproject.toml uv.lock)
    echo "prosapia pinned to ${commit:0:8} (uv.lock)"
    build='uv sync --frozen'
  else
    echo "prosapia pinned to ${commit:0:8} (--commit)"
    build="uv venv --python '>=3.11' .venv
uv pip install --python .venv/bin/python \"prosapia[modal] @ git+https://github.com/jlmoraleshellin/prosapia@$commit\""
  fi
  echo "==> workspace + marker";   cmd_init   || exit $?
  echo "==> syncing the working tree"; cmd_sync || exit $?
  echo "==> building the venv on a compute node (several minutes)"
  IFS= read -r -d '' script <<'EOF' || true
set -e
cd @WS@
test -f .sapia-workspace
module load uv
# Take the interpreter from the EasyBuild Python that `module load uv` brings along. By default
# uv prefers a managed Python, downloaded into $HOME, and the venv would then link into one
# user's home directory.
export UV_PYTHON_PREFERENCE=only-system
# the uv cache lives in $HOME and the venv on group storage: copy, do not warn about hardlinks
export UV_LINK_MODE=copy
echo "building with $(uv --version), $(python3 -V) on $(hostname)"
rm -rf .venv
@BUILD@
EOF
  script=${script//@WS@/$(q "$WS")}
  script=${script//@BUILD@/$build}
  run_ssh 1800 "srun -A $(q "$ACCOUNT") -c 8 -t 2:0:0 --job-name=venv_bootstrap bash -lc $(q "$script")" || exit $?
  echo "==> verifying"
  remote "sapia run --help >/dev/null && echo 'sapia OK'" 120
}

cmd_status() {
  local pin out script
  pin=$(resolve_commit 2>/dev/null || true)
  IFS= read -r -d '' script <<'EOF' || true
cd @WS@ || { echo "workspace: MISSING"; exit 0; }
if [ -f .sapia-workspace ]; then echo "marker:    ok"; else echo "marker:    MISSING (sync will refuse)"; fi
if [ -x .venv/bin/python ]; then
  echo "venv:      $(.venv/bin/python -V 2>&1) ($(readlink -f .venv/bin/python))"
  # a uv venv has no pip: read the commit pip/uv recorded in the dist-info instead
  echo "prosapia:  $(.venv/bin/python -I -c 'import importlib.metadata as m,json;d=m.distribution("prosapia");u=json.loads(d.read_text("direct_url.json") or "{}");v=u.get("vcs_info");print(u.get("url","?")+"@"+v["commit_id"] if v else "not a git install: "+d.version+" "+str(u))' 2>/dev/null || echo 'not installed')"
else
  echo "venv:      MISSING"
fi
if uvv=$(bash -lc 'module load uv >/dev/null 2>&1 && uv --version' 2>/dev/null); then echo "uv:        $uvv (module load uv)"; else echo "uv:        NOT available through 'module load uv': bootstrap cannot run"; fi
EOF
  script=${script//@WS@/$(q "$WS")}
  out=$(run_ssh 60 "bash -c $(q "$script")") || exit $?
  echo "$out"
  if [ -n "$pin" ]; then
    if grep -qi "@$pin" <<<"$out"; then echo "pin:       venv matches ${pin:0:8}"
    else echo "pin:       MISMATCH with ${pin:0:8}: re-run bash scripts/vib.sh bootstrap"; fi
  else
    echo "pin:       cannot be read from uv.lock here (editable?); pass --commit to bootstrap"
  fi
}

main() {
  local sub=${1:-}
  shift || true
  case $sub in
    exec)       load_config; cmd_exec "$@" ;;
    sync)       load_config; cmd_sync "$@" ;;
    check)      load_config; guard "$*" && echo "guard: ok" ;;
    init)       load_config; cmd_init "$@" ;;
    bootstrap)  load_config; cmd_bootstrap "$@" ;;
    status)     load_config; cmd_status ;;
    ''|-h|--help|help) usage ;;
    *) die "unknown command '$sub' (see: bash scripts/vib.sh help)" ;;
  esac
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then main "$@"; fi
