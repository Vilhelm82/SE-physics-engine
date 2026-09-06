"""rlq -- reflection-loop QEC harness (phase 1, 2026-09-06).

Libraries that experiments share, so that a driver is a page and the conventions are structural:
  receipt      named checks, forks, held-out lists, sealed predictions, provenance, artifact folders
  figures      Codex's plot grammar (3x2 grid, palette) as helpers
  channels     Pauli twirl, herald, stacked-leak metrics
  instruments  words x hold policies (absorbing hold = FD; capture = register), density-matrix propagation
  precision    multi-precision backends, Theorem-2 (Cella U-0784) fit, Richardson
  words        Codex's reflection_loop_* re-exported at their current paths (phase 2 moves them here)
  decoder      Codex's qec_distance_stack re-exported

Nothing of Codex's is moved in phase 1.  Old root runners keep working; new drivers live in experiments/<date>/.
"""
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
