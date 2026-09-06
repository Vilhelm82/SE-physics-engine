"""RH-1 numerical instruments on the declared instantaneous coherent P/Q interface.

Trine: bare RF arcs. FD: existing finite ramps and stretched physical inverses.
No absorbing hold, extra bath, fitted recovery, or Pauli twirl is inserted.
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))  # repo root on sys.path (reorg 2026-09-06)
import argparse
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import polar

from rlq.reflection_loop_composite import primitive, return_time
from rlq.reflection_loop_dynamics import embed, frame, segment, hamiltonian
from rlq.reflection_loop_finite_dump import controls, loop_endpoints
from rlq.reflection_loop_reference_echo import rotation

I2 = np.eye(2, dtype=complex)
I4 = np.eye(4, dtype=complex)
G = np.diag([-1., 1.])
Z = np.diag([1., -1.])
TAU = return_time()
TRINE = np.pi*np.array([1/3, 2/3, 1/3])


def db(e):
    s = np.linalg.svd(e, compute_uv=False)
    p = s*s
    return dict(DB=float((p[0]-p[1])/sum(p)) if sum(p)>1e-28 else None,
                singular_values=s.tolist(), mean=float(sum(p)/2),
                extrema=p[::-1].tolist())


def foldback(s, e):
    # A unitary polar completion is chosen even at rank one. Its null-space
    # freedom has no effect on C E at calibration; it cannot restore rank.
    return polar(s)[0]@polar(e)[0].conj().T


def pattern_kraus(branches):
    result = [I2.copy()]
    for choices in branches:
        result = [b@k for k in result for b in choices]
    return result


def channel_metrics(kraus, target=G):
    completeness = sum((k.conj().T@k for k in kraus), np.zeros((2,2),complex))
    err = 0.
    for k in kraus:
        a = target.conj().T@k
        err += np.linalg.norm(a-np.trace(a)*I2/2)**2/3
    fe = sum(abs(np.trace(target.conj().T@k))**2 for k in kraus)/4
    return dict(infidelity=float(err), direct_infidelity=float((2-2*fe)/3),
                tp_defect=float(np.linalg.norm(completeness-I2)),
                mean_trace=float(np.trace(completeness).real/2))


@lru_cache(None)
def haar_quadrature(nz=24, nphi=48):
    z, w = np.polynomial.legendre.leggauss(nz)
    phi = np.arange(nphi)*2*np.pi/nphi
    states = np.stack(np.broadcast_arrays(np.sqrt((1+z[:,None])/2),
                     np.sqrt((1-z[:,None])/2)*np.exp(1j*phi)), axis=-1).reshape(-1,2)
    return states, np.repeat(w/(2*nphi), nphi)


def conditional_metrics(k, target=G):
    a = target.conj().T@k
    tf = a-np.trace(a)*I2/2
    states, weights = haar_quadrature()
    out = states@a.T
    perpendicular = np.column_stack((-states[:,1].conj(), states[:,0].conj()))
    orth = np.sum(perpendicular.conj()*(states@tf.T), axis=1)
    survive = np.sum(abs(out)**2, axis=1)
    numerator = float(np.linalg.norm(tf)**2/3)
    mean_success = float(np.linalg.norm(k)**2/2)
    return dict(conditional=float(weights@(abs(orth)**2/survive)),
                success_weighted=numerator/mean_success, success=mean_success,
                error_numerator=numerator,
                success_extrema=np.linalg.eigvalsh(k.conj().T@k).tolist())


def replacement_kraus(amplitude):
    # CP map rho -> tr(amplitude rho amplitude†) I/2.
    return [np.outer(I2[:,a], amplitude[b,:])/np.sqrt(2)
            for a in range(2) for b in range(amplitude.shape[0])]


def replay_branches(s, e, calibrated_u):
    back = calibrated_u.conj().T@e
    return [s, s@back] + replacement_kraus(e@back)


def compose_discard(loops):
    k = I2.copy()
    leaks = []
    for loop in loops:
        s, e = loop[2:,2:], loop[:2,2:]
        leaks.append(e@k)
        k = s@k
    r = sum((f.conj().T@f for f in leaks), np.zeros((2,2),complex))
    out = conditional_metrics(k)
    out.update(erasure=float(np.trace(r).real/2),
               erasure_spread=float(np.ptp(np.linalg.eigvalsh(r))),
               unconditional=out['error_numerator']+float(np.trace(r).real/4),
               completeness_defect=float(np.linalg.norm(r+k.conj().T@k-I2)))
    return k, leaks, out


@lru_cache(None)
def trine_loop(alpha, eps=0., dd=0., dr=0., dt=0., lock=(0.,0.,0.,0.),
               frame_sign=-1, start=0., max_step=.12):
    rot = rotation(alpha)
    if dd == dr == 0 and not any(lock):
        u = np.eye(3, dtype=complex)
        for j in range(3):
            u = segment(j, TAU, 1+eps, TAU+(dt if j==2 else 0))@u
        return rot@embed(frame(2, np.pi/2*(1+dt/TAU))@u)@rot.T
    noise = np.diag([0.,0.,dd,dr])
    lock = np.array(lock)
    u = I4.copy()
    for j in range(3):
        def rhs(t, flat):
            h = (1+eps)*rot@hamiltonian(j,t,TAU)@rot.T
            if np.any(lock):
                phases = np.exp(frame_sign*1j*(start+j*TAU+t)*lock)
                h = phases[:,None]*h*phases.conj()[None,:]
            return (-1j*(h+noise)@flat.reshape(4,4)).ravel()
        sol = solve_ivp(rhs, (0,TAU+(dt if j==2 else 0)), u.ravel(),
                        method='DOP853', rtol=2e-13, atol=2e-15, max_step=max_step)
        if not sol.success:
            raise RuntimeError(sol.message)
        u = sol.y[:,-1].reshape(4,4)
    return u


def word_loops(kind, cases):
    if kind == 'trine':
        return np.array([[trine_loop(float(a),*case) for a in TRINE] for case in cases])
    if kind != 'five':
        raise ValueError(kind)
    out = loop_endpoints(cases)
    # Closed form at pure gain avoids integration noise at 1e-15 infidelity.
    for i,(eps,dd,dr) in enumerate(cases):
        if dd == dr == 0:
            f = primitive(eps)
            angles, signs, _ = controls()
            out[i] = [rotation(a)@(f if sign>0 else f.conj().T)@rotation(a).T
                      for a,sign in zip(angles,signs)]
    return out


def norm_json(value):
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return dict(real=value.real.tolist(), imag=value.imag.tolist())
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def receipt(test, result):
    files = [f'experiments/2026-09-06/{test}.py','rlq/rh1_common.py','experiments/2026-09-05/split1_split_loop.py',
             'rlq/reflection_loop_finite_dump.py','rlq/reflection_loop_dynamics.py',
             'docs/results/2026-09-06/2026-09-06-RH-1-spec.md','docs/prereg/RH-1/DESIGNER-PREDICTIONS.md']
    result['provenance'] = dict(test=test, date='2026-09-06',
        head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        unblinded=True, reason='User supplied the predictions path; read during context refresh before runs.',
        interface='Instantaneous coherent P/Q projection; FD finite ramps retained, dump duration zero.',
        sha256={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in files})
    checks = result.get('checks',[])
    result['check_count'] = len(checks)
    result['checks_passed'] = sum(c['passed'] for c in checks)
    parser=argparse.ArgumentParser()
    parser.add_argument('--json',nargs='?',const=f'docs/{test}-checks.json',default=None)
    args=parser.parse_args()
    if args.json:
        Path(args.json).write_text(json.dumps(result,default=norm_json,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k in ('summary','check_count','checks_passed')},default=norm_json,indent=2))
    if any(not c['passed'] for c in checks):
        raise SystemExit(1)


def check(checks, name, value, tolerance):
    checks.append(dict(name=name,value=float(value),tolerance=tolerance,
                       passed=bool(np.isfinite(value) and abs(value)<=tolerance)))
