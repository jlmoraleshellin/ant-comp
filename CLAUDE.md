# ant-comp: a protein-design workspace

Computational design of protein binders (and other designed proteins) on top of [`prosapia`](https://github.com/jlmoraleshellin/prosapia) (CLI: `sapia`). In prosapia **a design is a row and a generation of designs is a table**: every tool reads columns and writes columns, and each step is one `run` plus one `collect` against a `run_dir`. Nothing is designed or predicted locally: every tool runs remotely, on Modal (a container per task) or on the VIB DataCore (SLURM array jobs over ssh), and the data stays where it was produced. A `run_dir` lives on one backend.

## How work gets done

Claude drives the campaign and the backend does the compute. The main session is the **`thinker`** (`claude --agent thinker`): it owns the science, decides what to run, reads the results and delegates everything else. It has no shell. Execution belongs to **`vib-worker`** or **`modal-worker`**, one per backend, and the specialists are **`tracker`** (step ledger), **`campaign-writer`** and **`campaign-condenser`** (decision log), **`table-analyst`** (batch statistics), **`target-scout`** (target recon), **`tool-creator`** and **`tool-reviewer`** (building and reviewing tools).

## Where things live

| path | what |
| --- | --- |
| `.claude/agents/` | the agent definitions, including the standing protocol of each worker |
| `.claude/skills/` | one skill per tool (premise, flags, columns, traps) plus the workflow skills: `prosapia`, `running-a-step`, `verifying-a-tool`, `binder-campaign`, `ph-switch`, `authoring-a-tool`, `editing-a-tool` |
| `tools/` | project-local tools and forks of bundled ones; they shadow the bundled tool of the same name |
| `docs/vib.md` | how the cluster is reached (`scripts/vib.sh`), partitions, workarounds |
| `.env` | backend configuration and one `SAPIA_ACTIVATE_<TOOL>` entry per tool that needs an activation script on the cluster |
| `campaigns/` | per campaign: the decision log `<run_dir>.md` and the step ledger `<run_dir>_progress.html`; analysis exports and reports under `campaigns/.analysis/` |

When a campaign is resumed, the thinker starts from the ledger and the decision log, then asks a worker for the registry of the `run_dir`.
