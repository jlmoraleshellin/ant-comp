# vib: how this repo reaches the cluster

The Claude session stays on your machine. Compute runs on vib. `scripts/vib.sh` is the only path between them.

```
your machine                                     vib
------------                                     ---
.env, tools/, activation/, utils/  --rsync-->    $SAPIA_VIB_WORKSPACE/   (own .venv, outputs/)
vib-worker ---- bash scripts/vib.sh exec ---->   login node: squeue, sacct, sinfo, tail a log
           └─── exec --on-node ──────────────>   small srun (2 cpus, 8G, 30 min): sapia run / collect
                                                  └── sbatch → compute nodes → tool scripts
```

## The commands

| Command | Does |
| --- | --- |
| `exec [--sync] [--timeout S] '<cmd>'` | runs `<cmd>` in the workspace, venv active, `.env` loaded; `--sync` pushes first |
| `exec --on-node [--cpus N] [--mem S] [--walltime H:MM:SS] '<cmd>'` | the same, inside a small `srun` (defaults 2 cpus, 8G, 0:30:00; `--cpus` capped at 8; ssh timeout 1800 s) instead of the login node. Drops the job's own `SLURM_*` variables so a `sapia run` inside it does not inherit the small allocation's resources. Use it for `sapia run`, `sapia collect` and anything that reads tables; keep `squeue`/`sacct` polls on the login node |
| `sync [--dry-run]` | pushes `tools/`, `activation/`, `utils/` and `.env` |
| `check '<cmd>'` | runs the guard only |
| `init [--adopt]` | creates the workspace and its `.sapia-workspace` marker; `--adopt` for an existing one |
| `bootstrap [--commit SHA]` | init, sync (plus `pyproject.toml` and `uv.lock`), then builds the venv with `uv` on a compute node (`rm -rf .venv` first, so run it only when no job of yours is in flight). Default is `uv sync --frozen` from `uv.lock`; `--commit SHA` installs that prosapia commit instead, which must be pushed |
| `status` | marker, venv (and the interpreter it links to), prosapia commit against the pin, and whether `module load uv` works |

`$SAPIA_VIB_WORKSPACE` and `$SAPIA_VIB_ACCOUNT` expand **on vib**, so write them in single quotes.

## What it guarantees, and what it does not

- Every `exec` starts in the workspace. The guard refuses `cd`, `pushd`, `popd`, `chdir`, `..`, `~`, `$HOME`, ssh/scp/rsync from the cluster, absolute paths outside the workspace, and `rm` of anything outside it or of the workspace root itself (anything inside, `outputs/` and `.venv/` included, is allowed). A path reached through `$SAPIA_VIB_WORKSPACE` or `$SAPIA_VIB_SHARED_ACTIVATION` is allowed.
- `sync` copies only `tools/ activation/ utils/ .env`, and `--delete` acts only inside those directories. `outputs/`, `.venv/` and everything else in the workspace are unreachable by construction. It also refuses unless the remote directory holds the `.sapia-workspace` marker, so a mistyped path cannot be wiped.
- It is an accident guard, not a sandbox. Code inside `srun ... bash -lc '...'` is still your user. The hard limit is your Unix permissions on group storage.
- Every ssh runs under `timeout`. rc 255 or 124 usually means the certificate lapsed: sign in, do not retry in a loop.

## `.env`: three readers

| Section | Read by |
| --- | --- |
| 1. host, workspace, account | `scripts/vib.sh`, on your machine |
| 2-4. shared activation dir, `SAPIA_ACTIVATE_*` | prosapia on vib: at submit time via the workspace venv's `load_dotenv()`, at task time via the prelude's `source .env` in the submit directory |

Both prosapia reads must reach the same file, which is why the venv lives **inside** the workspace. Do not move it.

## The venv, `uv`, and the pin

- **`bootstrap` uses the cluster's `uv`**, loaded with `module load uv` inside the compute-node allocation (`bash -lc` is what makes `module` available). Default: `uv sync --frozen` from the `pyproject.toml` and `uv.lock` that `bootstrap` pushes, so the install matches Modal's. Nothing else uses `uv` on vib, and your machine still needs none for vib.
- **The interpreter is the EasyBuild Python that `module load uv` brings along** (3.13.5 when measured 2026-10-17). `bootstrap` sets `UV_PYTHON_PREFERENCE=only-system`, because otherwise `uv` prefers a managed Python downloaded into `$HOME` and the venv would link into one user's home directory. The venv's `python` is a symlink into the module tree and runs with an empty environment, so tasks do not need the module loaded. `UV_LINK_MODE=copy` stops `uv` warning that its cache (in `$HOME`) and the venv (group storage) are on different filesystems.
- **A `uv` venv has no `pip`.** `status` reads the installed prosapia commit from the package's `direct_url.json` instead of `pip freeze`.
- rsync never touches the venv. After a pin change, re-run `bootstrap`; `status` shows a mismatch.
- If `uv.lock` points at the editable `../prosapia` checkout, the default `bootstrap` stops. Pass `--commit <sha>` with a **pushed** commit: that path runs `uv venv` and `uv pip install` of the commit, with no lock, so transitive dependencies are resolved fresh and can differ from Modal's.

## The shared activation scripts

`SAPIA_VIB_SHARED_ACTIVATION` points at the lab's scripts for the built-in tools. The software they activate (rfd3, ProteinMPNN, Boltz, AF3...) is installed on the cluster and is not in this repo, so syncing the repo does not replace them. What the workspace venv replaced is the shared **python env** that `sapia` itself used to run from. (The scripts' contents live on the cluster and are not versioned here; copying them into `activation/vib/` would version them but not the software.)
