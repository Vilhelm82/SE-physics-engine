"""Channels: Pauli twirl from a Pauli transfer matrix, stacked-leak metrics, and Codex's native_channel for compatibility."""
import numpy as np
from scipy.linalg import polar
import rlq.qec_distance_stack as _Q

PAULIS = _Q.PAULIS; G = _Q.G
CHAR = np.array([[1, 1, 1, 1], [1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]])   # chi(P, a): commute +1, anticommute -1
native_channel = _Q.native_channel; operating_point = _Q.operating_point; unpack = _Q.unpack


def twirl_from_ptm(ptm):
    """Pauli weights (I, X, Y, Z) of the twirl of a map with PTM diagonal ptm (relative to the target)."""
    return CHAR@np.asarray(ptm)/4


def average_infidelity(weights, accepted=1.):
    """(2/3) (w_X + w_Y + w_Z)/accepted for a twirled qubit channel."""
    w = np.asarray(weights); return float(2*w[1:].sum()/3/accepted)


def stacked_leak(loops):
    """From chronological 4x4 loop propagators (P = 2,3; Q = 0,1): survivor K, leaks F_k, stacked F, its polar isometry, DB, c."""
    K = np.eye(2, dtype=complex); leaks = []
    for L in loops:
        leaks.append(L[:2, 2:]@K); K = L[2:, 2:]@K
    F = np.vstack(leaks); U, H = polar(F); s = np.linalg.svd(F, compute_uv=False); p = s*s
    return dict(K=K, leaks=leaks, F=F, U=U, H=H, singular_values=s, c=float(p.sum()/2), DB=float((p[0] - p[1])/p.sum()),
                identity_defect=float(np.linalg.norm(F.conj().T@F + K.conj().T@K - np.eye(2))))
