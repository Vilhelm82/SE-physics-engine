# seated-root

A derivational-economy study: build relativistic and quantum structure from minimal geometric primitives — one root,
three lines, three planes, one seated observer — without importing algebras that already contain the target results.

## Read these three

| file | what it is |
|---|---|
| **`MODEL.md`** | The model. Part I: the labelled octahedron. Part II: the dynamics tier — what is derived, at what tier, with the runner that shows it. Self-contained. |
| **`HANDOFF.md`** | Only what `MODEL.md` cannot say: what to do next, and the state of work in progress. |
| **`AGENTS.md`** | Working instructions: how dependencies are recorded, how things are named, how results are labelled and committed. |

## The tree

| path | contents |
|---|---|
| `suites/` | Runners. Each is a self-contained check with its own header stating inputs, bans and kill conditions; each has a `.log` beside it from the run that was committed. Mixed eras: pre-pivot (Cl(3)) studies sit alongside the post-pivot rebuild, and a runner's header states which it is. |
| `docs/` | Dated records. `results/<date>/` is what was found on that date and is not edited afterwards; `prereg/` holds pre-registered protocols; `archive/` is superseded material. |
| `docs/EXTERNAL-MATHEMATICS-DEBT.md` | The dependency ledger. Every adopted external result with its source, exact use, assumptions, and whether it is retained standard mathematics or leaves a native construction owed. |
| `docs/CONJECTURE-COSMOLOGY.md` | The white-sheet conjecture, quarantined: CONJECTURE tier, with its kill conditions and their status. |
| `handoffs/` | Superseded handoffs, provenance only. The current one is `HANDOFF.md` in the root. |
| `paper/` | Draft and build scripts. |
| `rlq/`, `experiments/`, `runs/`, `tools/`, `external/` | Library code, sweeps, run records, utilities, and material brought in from outside. |
| `tests/` | Test suite. |

## Where the work stands

Residue: one declared line (KMS) and two quantities pinned to measurement (d = 3; K from LAGEOS). Everything else in
Part II is derived from the Newtonian potential, Laplace, Euler, induction and the primitives. Next open construction: whether a second propagating channel exists, and whether it is spin-2. See `MODEL.md` §13
and `HANDOFF.md`.
