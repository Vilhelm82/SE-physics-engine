"""Instruments: a word times a hold policy, propagated as 14-dim density matrices with loss attribution.

FD's instrument (absorbing holds, monitored) is `Model(capture=False)`; the register version adds a coherent capture before each
absorbing hold.  Moved verbatim from rh4_capture.py (2026-09-06) so that the numbers are bit-identical; the old module re-exports
these names.  Conventions: P = indices 2,3 (d, r); Q = 0,1 (p, q); register slot k = 4+2k, 5+2k.
"""
import numpy as np
from scipy.linalg import expm, polar
import rlq.reflection_loop_finite_dump as FD
import rlq.qec_distance_stack as Q

N = 5; DIM = 4 + 2*N
CASE = (0.01, 1e-4, 0.0); GAMMA = 1e-6; KAPPA, ACTION = 10., 40.
PAULI2 = Q.PAULIS
CHAR = np.array([[1, 1, 1, 1], [1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]])   # chi(P, a): +1 commute, -1 anticommute
IP, IQ = [2, 3], [0, 1]
IR = lambda k: [4 + 2*k, 5 + 2*k]
ALLR = [m for k in range(N) for m in IR(k)]


def embed4(m):
    out = np.eye(DIM, dtype=complex); out[:4, :4] = m; return out


def embed4z(m):
    """Zero-padded embedding for EFFECT operators (an effect on P+Q must not act on the register)."""
    out = np.zeros((DIM, DIM), dtype=complex); out[:4, :4] = m; return out


class Model:
    """Word x hold policy; capture=False is FD exactly."""
    """Steps are no-jump matrices M (rho -> M rho M^dag); dephasing is a separate elementwise factor."""
    def __init__(self, e_c=0., t_c=1., gamma_r=0., gamma_phi=0., register_heralded=True, capture=True,
                 case=CASE, gamma=GAMMA, n=1, loops=None, effects=None):
        self.th = (np.pi/2)*(1 + e_c); self.t_c = t_c if capture else 0.; self.capture = capture
        self.gr, self.gphi, self.rh = gamma_r, gamma_phi, register_heralded
        self.case, self.gamma, self.n = tuple(case), gamma, n
        if loops is None:
            loops, effects = FD.loop_endpoints_with_effects([case], [gamma], n=n); loops, effects = loops[0], effects[0]
        self.loops, self.effects = loops, effects
        stages = FD.primitive_stages(n); T0 = sum(s.duration for s in stages)
        _, _, stretches = FD.controls('five'); self.loop_T = [T0*float(s) for s in stretches]
        self.td = ACTION/KAPPA; self.chi = np.exp(-ACTION/2); self.delta = np.array([case[1], case[2]])
        self.T_total = sum(self.loop_T) + N*(self.td + self.t_c)

    # ---- matrices
    def M_phases_P(self, t):
        d = np.ones(DIM, dtype=complex); d[IP] = np.exp(-1j*self.delta*t); return np.diag(d)
    def M_regdecay(self, t, filled):
        d = np.ones(DIM); d[[m for k in filled for m in IR(k)]] = np.exp(-self.gr*t/2); return np.diag(d)
    def M_loop(self, k): return embed4(self.loops[k])
    def M_capture(self, k):
        M = np.eye(DIM, dtype=complex)
        if not self.capture or self.t_c == 0: return M
        gen = np.array([[-self.gamma/2, -1j*self.th/self.t_c], [-1j*self.th/self.t_c, -self.gr/2]])
        blk = expm(gen*self.t_c)
        for m, r in zip(IQ, IR(k)): M[np.ix_([m, r], [m, r])] = blk
        return M@self.M_phases_P(self.t_c)
    def M_absorb(self):
        d = np.ones(DIM, dtype=complex); d[IQ] = self.chi; return np.diag(d)@self.M_phases_P(self.td)
    def dephase(self, rho, t, filled):
        if not self.gphi or not filled: return rho
        f = np.zeros(DIM); f[[m for k in filled for m in IR(k)]] = 1.
        decay = np.exp(-self.gphi*t*(f[:, None] + f[None, :])); np.fill_diagonal(decay, 1.)
        return rho*decay

    # ---- pure linear map (valid for gamma_phi = 0): total no-jump matrix, chronological
    def total_matrix(self):
        M = np.eye(DIM, dtype=complex); filled = []
        for k in range(N):
            M = self.M_regdecay(self.loop_T[k], filled)@self.M_loop(k)@M
            M = self.M_regdecay(self.t_c, filled)@self.M_capture(k)@M; filled.append(k)
            M = self.M_regdecay(self.td, filled)@self.M_absorb()@M
        return M

    # ---- density-matrix propagation with loss attribution
    def run(self, rho_P):
        rho = np.zeros((DIM, DIM), dtype=complex); rho[np.ix_(IP, IP)] = rho_P
        her = unh = 0.; filled = []
        def apply(M, rho, heralded):
            new = M@rho@M.conj().T; lost = float(np.trace(rho).real - np.trace(new).real)
            return new, (lost if heralded else 0.), (0. if heralded else lost)
        for k in range(N):
            unh += float(np.trace(embed4z(self.effects[k])@rho).real)                     # background during transport (exact effect)
            rho = self.M_loop(k)@rho@self.M_loop(k).conj().T
            rho, a, b = apply(self.M_regdecay(self.loop_T[k], filled), rho, self.rh); her += a; unh += b
            rho = self.dephase(rho, self.loop_T[k], filled)
            rho, a, b = apply(self.M_capture(k), rho, False); her += a; unh += b            # Q background + register loss during t_c: unheralded
            rho, a, b = apply(self.M_regdecay(self.t_c, filled), rho, self.rh); her += a; unh += b
            rho = self.dephase(rho, self.t_c, filled); filled = filled + [k]
            rho, a, b = apply(self.M_absorb(), rho, True); her += a; unh += b               # monitored absorption: heralded (FD)
            rho, a, b = apply(self.M_regdecay(self.td, filled), rho, self.rh); her += a; unh += b
            rho = self.dephase(rho, self.td, filled)
        return rho, her, unh


def channel(model, U_cal):
    """Herald and exact Pauli twirl of the accepted map, qec_distance_stack format, plus FD's decomposition."""
    Gm = Q.G; parts = {'survivor': np.zeros(4), 'recovered': np.zeros(4), 'replaced': np.zeros(4)}; heralds = []; unhs = []
    base = {}
    for a, sig in enumerate(PAULI2):
        rho_in = sig/2 if a == 0 else sig/2 + np.eye(2)/2
        rho, her, unh = model.run(rho_in)
        rP = rho[np.ix_(IP, IP)]; rR = rho[np.ix_(ALLR, ALLR)]; rQ = float(np.trace(rho[np.ix_(IQ, IQ)]).real)
        rec = U_cal.conj().T@rR@U_cal; leftover = float(np.trace(rR).real - np.trace(rec).real)
        outs = {'survivor': Gm.conj().T@rP@Gm, 'recovered': rec, 'replaced': (unh + rQ)*np.eye(2)/2}
        heralds.append(her + leftover); unhs.append(unh + rQ)
        for name, out in outs.items():
            if a == 0: base[name] = out; parts[name][0] = float(np.trace(out).real)
            else: parts[name][a] = float(np.trace(sig@(out - base[name])).real)
    h = float(heralds[0]); spread = float(np.ptp(heralds)); unh0 = float(unhs[0])
    w = {name: CHAR@ptm/4 for name, ptm in parts.items()}; wt = sum(w.values())
    assert abs(wt.sum() + h - 1) < 1e-9, (wt.sum(), h)
    cond = wt/(1 - h); cond[0] = 1 - cond[1:].sum()
    coh = {name: float(2*w[name][1:].sum()/3/(1 - h)) for name in ('survivor', 'recovered')}   # average infidelity contributions
    return dict(herald=h, herald_state_spread=spread, pauli=cond.tolist(), unheralded_pauli_rate=float(wt[1:].sum()), total_time=model.T_total,
                unobserved_absorption=unh0/(1 - h), surviving_error=coh['survivor'], recovered_error=coh['recovered'],
                unheralded_residual=unh0/(1 - h) + coh['survivor'] + coh['recovered'],          # FD's raw q_u: loss counts as a full error
                repaired_residual=float(2*wt[1:].sum()/3/(1 - h)), recovered_fraction=float(np.trace(base['recovered']).real))


def calibrate(t_c, **kw):
    ideal = Model(e_c=0., t_c=t_c, gamma_r=0., gamma_phi=0., **kw)
    F = ideal.total_matrix()[np.ix_(ALLR, IP)]; U, H = polar(F)
    s = np.linalg.svd(F, compute_uv=False)
    return U, dict(captured_fraction=float((s**2).sum()/2), DB_stack=float((s[0]**2 - s[1]**2)/(s**2).sum()))


