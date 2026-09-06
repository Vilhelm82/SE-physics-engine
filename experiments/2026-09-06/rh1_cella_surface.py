"""Cella three-channel analysis of a declared trine DSP angle chart.

Offline geometry/constitutive analysis, not an autonomous physical controller.
The native Cella exact referee evaluates nominal exact jets and rationalized
finite-difference jets. The latter retain numerical, not exact-jet, status.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.linalg import polar
from scipy.optimize import brentq

from rh1_common import G, TRINE, compose_discard, trine_loop

CELLA = Path('/home/williaml/Cella Framework')
REFEREE = CELLA / 'engine/src/cella/reference_lift.py'
Y = np.array([[0, -1j], [1j, 0]])


def load_referee():
    spec = importlib.util.spec_from_file_location('cella_reference', REFEREE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.three_channel_referee


def signed_constraint(x, eps, dd):
    k, _, metrics = compose_discard([trine_loop(float(a), eps, dd) for a in x])
    if np.linalg.svd(k, compute_uv=False)[-1] < 1e-8:
        raise ValueError('singular surviving map outside the local chart')
    w = polar(k)[0] @ G
    if abs(np.trace(w)) < 1e-8:
        raise ValueError('trace branch outside the local chart')
    value = .5 * np.arctan(np.real(1j * np.trace(Y @ w) / np.trace(w)))
    return float(value), metrics


def jet(x, eps, dd, h):
    basis = np.eye(3) * h
    f = lambda y: signed_constraint(y, eps, dd)[0]
    center = f(x)
    g, hessian = np.zeros(3), np.zeros((3, 3))
    for i in range(3):
        plus, minus = f(x + basis[i]), f(x - basis[i])
        g[i] = (plus - minus) / (2 * h)
        hessian[i, i] = (plus - 2 * center + minus) / h**2
        for j in range(i):
            hessian[i, j] = hessian[j, i] = (
                f(x + basis[i] + basis[j]) - f(x + basis[i] - basis[j])
                - f(x - basis[i] + basis[j]) + f(x - basis[i] - basis[j])) / (4 * h**2)
    return center, g, hessian


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json', type=Path)
    args = parser.parse_args()
    referee = load_referee()
    d, s, p = sp.symbols('D S P', real=True)
    f = d - s + p
    g = sp.Matrix([sp.diff(f, x) for x in (d, s, p)])
    hessian = sp.hessian(f, (d, s, p))
    exact = referee([str(x) for x in g], [[str(x) for x in row] for row in hessian.tolist()])
    assert all(v == 0 for v in exact['channels'].values())
    a = sp.ones(3, 1)
    assert (g.T * a)[0] == 1
    assert (g.T * a * a.T * g)[0] == 1
    keystone = referee([3, 1, 2], [[2, 1, 0], [1, 0, 0], [0, 0, 2]])
    assert [str(keystone['channels'][k]) for k in ('kappa_c', 'kappa_s', 'kappa_int', 'K_G')] == ['-1/49', '1/49', '-3/49', '-3/49']
    cases = []
    for eps, dd in ((.01, 1e-4), (.01, -1e-4), (.01, 1e-3)):
        shift = brentq(lambda u: signed_constraint(TRINE + u, eps, dd)[0], -.003, .003, xtol=1e-14)
        point = TRINE + shift
        levels = []
        for h in (.002, .001, .0005):
            value, grad, hes = jet(point, eps, dd, h)
            account = referee([str(v) for v in grad], [[str(v) for v in row] for row in hes])
            assert account['bordered_matches_partition']
            assert abs(value) < 1e-12
            assert abs(sum(grad)) > .99
            levels.append(dict(step=h, F=value, gradient=grad.tolist(), hessian=hes.tolist(),
                channels={k: float(v) for k, v in account['channels'].items()},
                common_angle_derivative=float(sum(grad)),
                exact_partition_for_supplied_decimal_jet=True))
        fine, medium = levels[-1], levels[-2]
        for key in fine['channels']:
            assert abs(fine['channels'][key] - medium['channels'][key]) < 1e-3 * max(abs(fine['channels'][key]), 1e-14)
        cases.append(dict(eps=eps, dd=dd, root_shift=shift, point=point.tolist(),
            metrics=signed_constraint(point, eps, dd)[1], levels=levels))
    paths = [Path(__file__), Path('rh1_common.py'), REFEREE,
        CELLA / 'Papers_Library/05_expository_companions_and_research_maps/dbp_role_channel_and_orbit_geometry/Three_Channel_KG_Strong_Spec.md',
        CELLA / 'Papers_Library/02_theorems_and_lemmas/dbp_role_channel_and_orbit_geometry/Canonical_Invariant_Reduction_Theorem.md']
    result = dict(date='2026-09-06',
        interpretation='Ordered loop angles (alpha_1,alpha_2,alpha_3) instantiate (D,S,P); other role realizations are not inferred.',
        defining_function='F=0.5*atan(real(i*trace(Y*polar(K)*G)/trace(polar(K)*G))) on the local trace-nonzero branch',
        exact_nominal=dict(F=str(f), gradient=[str(v) for v in g],
            hessian=[[str(v) for v in row] for row in hessian.tolist()],
            channels={k: str(v) for k, v in exact['channels'].items()}, common_angle_derivative='1'),
        scope='Signed Y-component closure only under noise; not zero erasure or complete logical correction.',
        arithmetic='Nominal jet exact. Noisy jets use convergent finite differences; exact rational referee arithmetic does not make those jets exact.',
        response='Declared E=k*F^2/2, M=m*a*a^T, x_dot=-M*grad(E), a=(1,1,1). A physical force/bath implementation is not supplied.',
        cases=cases, verification='Nominal exact channel zeros, keystone signs, actuator transversality, all determinant partitions and 0.1 percent fine/medium channel consistency passed.',
        sha256={str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths})
    if args.json:
        args.json.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(dict(exact_nominal=result['exact_nominal'], verification=result['verification'],
        cases=[dict(eps=c['eps'], dd=c['dd'], root_shift=c['root_shift'],
                    channels=c['levels'][-1]['channels'],
                    common_angle_derivative=c['levels'][-1]['common_angle_derivative']) for c in cases]), indent=2))


if __name__ == '__main__':
    main()
