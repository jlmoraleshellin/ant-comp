# Overnight run — handover

**Started:** 2026-10-07 18:24 · **run_dir:** `outputs/20261007_155900_tnf_1000`
on volume `sapia-runs-luca` · **expected finish:** ~23:30

## What is running

`campaigns/tnf/run_campaign.sh`, detached with `setsid` (own session, survives
the terminal). 5 epitopes x 25 jobs x 8 designs = **1000 backbones**, one GPU at
a time.

Confirmed working at 18:28: first task returned **8 structures**.

Progress log: `/tmp/claude-1000/-home-lpignatti/9763510a-.../scratchpad/campaign2.log`
— one `tasks N/25  pdbs N` line per 5 minutes per epitope.

| label | target | epitope |
| --- | --- | --- |
| `table0_f1` | TNF_E1 | A32,33,145 + C87,91,92 — 6/6 exposed, composite, 0 mouse diffs, basic |
| `table0_f2` | TNF_E2 | A87,88,90,92,131,135 — de-buried rework, single protomer, basic |
| `table0_f3` | TNF_E3 | A113,115 + C72,75,77,97 — LPC's benchmark region, mostly buried |
| `table0_f4` | TNF_E4 | A24,65,67,138,140,141 — **control**: acidic + 3 exposed mouse diffs |
| `table0_f5` | TNF_E5 | none — unbiased |

Ignore `table0_e1`, `table0_e2`, `table0_e1b` — reserved by failed attempts, empty.

## First thing in the morning

```bash
cd ~/ant-comp
# how far it got
tail -40 /tmp/claude-1000/-home-lpignatti/9763510a-*/scratchpad/campaign2.log

# collect each epitope that finished
for L in f1 f2 f3 f4 f5; do
  .venv/bin/sapia modal-shell --cmd \
    "sapia collect laproteina outputs/20261007_155900_tnf_1000 -t table0_$L"
done
```

Then compare epitopes on `laproteina_binder_length`, and re-run
`campaigns/tnf/scripts/contacts.py` over a sample of each table to see **which
designs actually reached their hotspots** — on the smoke test only 2 of 6
hotspots were contacted, and whether that holds at n=200 is the first real
scientific question.

## Three open items

### 1. Histidines — unsolved, and it is the assignment

The goal is binding at pH 7.4 and release at pH 6.0, via His protonation. But
**all four atomium sequences yesterday contained zero histidines.** MPNN-family
designers strongly under-use His on idealised de novo bundles.

Filtering cannot fix this: a filter selects from what the designer produces, and
it produces none. The lever is `--bias-aa 'H:<w>'` **upstream**. And a global
bias scatters His through the fold where it does nothing — what the mechanism
needs is His *at the interface*, which argues for position-restricted bias over
a blanket one.

Also worth checking once designs exist: the epitope should present **Arg/Lys**
to the binder, so a protonated His+ is repelled. E1/E2/E3 are basic; **E4 is
acidic**, where a His+ would form a salt bridge and *strengthen* binding at low
pH. That is why E4 is the control.

### 2. `batch_size` override ignored in the container — sidestepped, not explained

`++generation.dataloader.batch_size=8` is accepted on the command line and
Hydra honours it locally with the identical override sequence, yet the task
resolved `batch_size: 16` and OOM-ed a 40 GiB A100 on a ~520-residue complex.

Current workaround: 25 jobs x 8 designs, so the **dataset** per job is 8 and the
batch cannot exceed it. Robust, but it costs 25 container starts per epitope
(~60 min instead of ~15). Worth finding the real cause before the next campaign.

### 3. The manifest race

The first E1 attempt died with `LPC_NAME: manifest field 1 (name) is empty` and
the prelude warning `ignored null byte in input` — a sparse read of a manifest
that is correct on disk. A second, identical submission a minute later read it
fine.

`laproteina.sh` now re-reads its line from `$MANIFEST` with backoff and strips
NULs. That makes the tool resilient; it does not fix the underlying visibility
race, which is prosapia/Modal territory and worth reporting upstream.

## Cost so far

Negligible. Every failure died at build or import; the one task that reached a
GPU before tonight lasted ~90 seconds. Tonight's 1000 designs are a few dollars.

## Not yet done

- `--num-jobs` > 1 had never run before tonight — this campaign is its first
  real exercise. Check the per-job `.exit` codes when collecting.
- Nothing has been validated with Boltz. No design is yet known to fold.
- Mouse is deferred to a final screen of the best designs.
