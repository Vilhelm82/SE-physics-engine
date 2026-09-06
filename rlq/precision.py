"""Multi-precision backends and Cella U-0784 machinery, re-exported from rh2_precision; plus Richardson."""
from rlq.rh2_precision import loop_np, loop_mp, stack_np, stack_mp, paths_np, paths_mp, EPS_GRID, U32, U64  # noqa


def richardson(values, order=3):
    """Repeated Richardson extrapolation of a sequence sampled at halving step (v[i] at h/2^i)."""
    seq = list(values)
    for j in range(1, order + 1):
        seq = [((2**j)*seq[i+1] - seq[i])/(2**j - 1) for i in range(len(seq) - 1)]
    return seq
