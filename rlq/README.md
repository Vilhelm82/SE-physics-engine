# rlq -- reflection-loop QEC harness

Phase 1 (2026-09-06). A driver is a page; the conventions are structural; every run gets a folder.

## Layout

```
rlq/                 the library (this package)
  receipt.py         Receipt: named checks, forks, held-out list, sealed predictions, corrections, provenance, artifacts
  instruments.py     word x hold policy, propagated as 14-dim density matrices; Model(capture=False) IS FD
  channels.py        Pauli twirl from a PTM, stacked-leak metrics (DB, polar isometry), Codex's native_channel
  precision.py       numpy/mpmath backends, Theorem-2 (Cella U-0784) paths, Richardson
  paper/figures.py         Codex's plot grammar: 3x2 constrained grid, palette, finish()
  words.py           Codex's reflection_loop_* re-exported at their current paths (phase 2 moves them here)
  decoder.py         Codex's qec_distance_stack re-exported
experiments/<date>/  drivers; one file per experiment, ~30-70 lines
runs/<date>/<exp>/   artifacts written by Receipt: receipt.json, run.log, *.png, predictions.md (if sealed)
runs/INDEX.md        one line per run, appended automatically
```

## A driver

```python
from rlq.receipt import Receipt
r = Receipt('rh4_capture', 'capture model at the absorbing interface', driver=__file__)
r.held_out('occupancy readout as an operation', 'hardware')      # required: close() refuses an empty list
r.predictions('docs/prereg/RH-4/DESIGNER-PREDICTIONS.md')          # optional: copied in, hash recorded
r.check('capture off == FD herald', value, 1e-9)                   # must pass, or the run exits nonzero
r.fork('ideal register', 'a', 'herald 0, unheralded = background') # an outcome; never pass/fail
r.record('grid', rows); r.figure(fig, 'register_requirement')      # data and pictures land in the run folder
sys.exit(r.close())
```

The four statements a run can make, and what each is:

| statement | meaning | pass/fail? |
|---|---|---|
| `check` | an identity, regression or sanity condition; if false the numbers cannot be trusted | yes -- the count in the commit message |
| `fork` | a physical outcome the spec pre-declared as (a)/(b)/...; a verdict is never written before a run | no -- recorded |
| `held_out` | what this run does not cover; where the next run comes from | required |
| `predictions` | the sealed designer guesses, so being wrong is on the record | hashed |

## Rules kept from before

- Old root runners are not moved in phase 1. `rh4_capture.py` is now a shim over `rlq.instruments` and reproduces its committed
  receipt with zero difference. Replay lines in docs/ still work.
- Receipts hash the driver, every `rlq/*.py`, and any module the driver names, plus the git head and a dirty flag.
- Phase 2 (coordinated with Codex): move `reflection_loop_*` and `qec_distance_stack` under `rlq/`, retire root drivers into
  `experiments/<date>/` with updated replay lines, and give the theorem suites (`thm_*`, `prim*`, ...) their own harness.

Acceptance test: `experiments/2026-09-06/rh4_capture.py` reproduces `docs/rh4_capture-checks.json` to 1.3e-16 on all 40 grid rows.
