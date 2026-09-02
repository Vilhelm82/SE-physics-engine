#!/usr/bin/env python3
# =============================================================================
# THM-TARGET J -- does the model's own geometry produce the field law?  (schedule item 0b, 2026-09-03)
# sympy + mpmath only.  Exit 0 iff every check passes.
#
# WHAT THE DATA DEMAND (thm_i_field.py, 10/10): the harmonic variable is sech^2 l = Delta/det G', the reciprocal
#   presented volume; every rival profile is dead at Mercury.
# WHAT A STATIC ENERGY CAN DO (D-1, proved below): any gradient-quadratic energy sum_ab M_ab(g) grad g_a . grad g_b of
#   click-invariant fields g_a(x), restricted to the spherically symmetric pinned family where every g_a is a function
#   of lambda(r), is F(lambda) lambda'^2, and its Euler-Lagrange equation makes the ARC LENGTH s = int sqrt(F) d lambda
#   harmonic.  So "produces KIN-2a" means "s is affine in sech^2 l": to fourth order, s4/s2 = -2/3.
# THE MODEL'S OWN CANDIDATE (D-3): the harmonic-map (Dirichlet) energy of x -> presented state in the FROZEN state
#   metric of docs/GRAM_SUBMERSION_CURVATURE: product round metric on (S^2)^3, mechanical connection, horizontal lift
#   by the 9x9 system  v_i . a_i = 0,  sum a_i x v_i = 0,  a_i . v_j + a_j . v_i = dG_ij,  quotient metric |v|^2.
#   No free choice.  Kill condition: s4/s2 != -2/3  ->  c_2 = -(s4/s2 + 2/3) != 0  ->  perihelion 3 + 2 c_2 != 3.
# CARRIED: H-5/H-7 presented Gram; KIN-2a as the target profile; thm_i_field B-2e (perihelion = 3 + 2 c_2).
# =============================================================================
import sys
import sympy as sp
import mpmath as mp
mp.mp.dps = 30
CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" -- {n}" if n else "")); sys.stdout.flush()
def z(e):
    e = sp.sympify(e)
    for rt in (sp.simplify, lambda q: sp.simplify(sp.expand(q))):
        try:
            if rt(e) == 0: return True
        except Exception: pass
    return False

print("=" * 100); print("D-1  THEOREM: a gradient-quadratic static energy makes the arc length along the family harmonic"); print("=" * 100)
r = sp.Symbol('r', positive=True); lam = sp.Function('lambda')(r); F = sp.Function('F'); s = sp.Function('s')
Lag = F(lam)*sp.diff(lam, r)**2*r**2                            # any sum_ab M_ab grad g_a . grad g_b reduces to this on the family
EL = sp.diff(Lag, lam) - sp.diff(sp.diff(Lag, sp.diff(lam, r)), r)
# with F = (ds/d lambda)^2 the EL equation is (r^2 ds/dr)' = 0
sub = {F(lam): sp.diff(s(lam), lam)**2}
EL_s = sp.simplify(EL.subs(sub).doit())
harmonic = sp.diff(r**2*sp.diff(s(lam), r), r)                # (r^2 s')' with s = s(lambda(r))
ratio = sp.simplify(EL_s/harmonic)
check("D-1a Euler-Lagrange of int F(l) l'^2 r^2 dr equals -2 s_l (r^2 ds/dr)' when F = s_l^2: stationary iff r^2 ds/dr = const,"
      " i.e. the arc length s(lambda) is HARMONIC, for every F", z(ratio + 2*sp.diff(s(lam), lam)))
check("D-1b so the data's 'sech^2 l harmonic' is the statement s = a + b sech^2 l along the family; to fourth order in lambda,"
      " s = s2 l^2 + s4 l^4 + ... must have s4/s2 = -2/3 (since tanh^2 l = l^2 - (2/3) l^4 + ...)",
      z(sp.series(sp.tanh(sp.Symbol('L'))**2, sp.Symbol('L'), 0, 7).removeO() - (sp.Symbol('L')**2 - sp.Rational(2, 3)*sp.Symbol('L')**4 + sp.Rational(17, 45)*sp.Symbol('L')**6)))

print(); print("=" * 100); print("D-2  OBSTRUCTION: no polynomial invariant of the presented cell is affine in sech^2 l"); print("=" * 100)
L = sp.Symbol('lambda', positive=True); k = sp.Symbol('k', positive=True); H = 1/sp.cosh(L)**2
cands = {'cosh^2 l (rho_K)': sp.cosh(L)**2, 'sinh^2 l': sp.sinh(L)**2, 'presented rod rho_i = 1 + k^2 sinh^2 l': 1 + k**2*sp.sinh(L)**2,
         'presented length sqrt(rho_i)': sp.sqrt(1 + k**2*sp.sinh(L)**2), "det G'/Delta = cosh^2 l": sp.cosh(L)**2}
def affine_in_H(g):                                            # d^2 g / dH^2 == 0 ?
    dgdH = sp.diff(g, L)/sp.diff(H, L); return z(sp.simplify(sp.diff(dgdH, L)/sp.diff(H, L)))
check("D-2a each polynomial-type invariant along the family is NOT affine in sech^2 l (second derivative w.r.t. sech^2 l nonzero): " +
      ", ".join(nm for nm in cands), not any(affine_in_H(g) for g in cands.values()))
check("D-2b ...while Delta/det G' = sech^2 l is (trivially), and it has a POLE where the presented volume vanishes (the centre,"
      " cosh l = i sinh mu -> 0) and is BOUNDED at the horizon (l -> oo): the opposite of every polynomial invariant",
      affine_in_H(H) and sp.limit(H, L, sp.oo) == 0 and sp.limit(sp.cosh(L)**2, L, sp.oo) == sp.oo)
print("  READING: the static energy that produces KIN-2a cannot be polynomial in the presented cell invariants; it must be built")
print("  from the RECIPROCAL of the volume invariant det G' (or another function bounded at the horizon with a pole at the centre).")

print(); print("=" * 100); print("D-3  THE MODEL'S OWN CANDIDATE: harmonic-map energy in the frozen quotient metric of the elliptope"); print("=" * 100)
Delta_of = lambda g: 1 - g[0]**2 - g[1]**2 - g[2]**2 + 2*g[0]*g[1]*g[2]
def frame_vectors(g):
    a, b, c = g; D = Delta_of(g)
    return [sp.Matrix([1, 0, 0]), sp.Matrix([a, sp.sqrt(1 - a**2), 0]), sp.Matrix([b, (c - a*b)/sp.sqrt(1 - a**2), sp.sqrt(D/(1 - a**2))])]
def presented_dirs(g, layer):
    """unit presented directions vhat_i(lambda) and the normalised presented Gram Gam(lambda), exact sympy in lambda"""
    A = frame_vectors(g); K = A[2] if layer == 'A' else (A[0] + A[1] + A[2])/sp.sqrt(((A[0] + A[1] + A[2]).T*(A[0] + A[1] + A[2]))[0, 0])
    vs = [ai + (sp.cosh(L) - 1)*(ai.T*K)[0, 0]*K for ai in A]
    vh = [v/sp.sqrt((v.T*v)[0, 0]) for v in vs]
    return vh
def lift_norm_fn(g, layer):
    """returns f(lambda) = ||d Gam/d lambda||_B (quotient metric) and the full (S^2)^3 speed of the actual presented lift"""
    vh = presented_dirs(g, layer)
    dvh = [sp.diff(v, L) for v in vh]
    Gam = lambda i, j: (vh[i].T*vh[j])[0, 0]
    dG = {(i, j): sp.diff(Gam(i, j), L) for (i, j) in ((0, 1), (0, 2), (1, 2))}
    vh_f = [sp.lambdify(L, v, 'mpmath') for v in vh]; dvh_f = [sp.lambdify(L, v, 'mpmath') for v in dvh]
    dG_f = {key: sp.lambdify(L, e, 'mpmath') for key, e in dG.items()}
    def cross(a, b): return mp.matrix([a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]])
    def quotient_speed(lv):
        a = [mp.matrix(vh_f[i](lv)) for i in range(3)]
        M = mp.matrix(9, 9); rhs = mp.matrix(9, 1); row = 0
        for i in range(3):                                   # v_i . a_i = 0
            for c in range(3): M[row, 3*i + c] = a[i][c]
            row += 1
        for c in range(3):                                   # sum_i a_i x v_i = 0, component c
            for i in range(3):
                for d in range(3):
                    e = mp.matrix(3, 1); e[d] = 1; M[row, 3*i + d] = cross(a[i], e)[c]
            row += 1
        for (i, j) in ((0, 1), (0, 2), (1, 2)):               # a_i . v_j + a_j . v_i = dG_ij
            for c in range(3): M[row, 3*j + c] += a[i][c]; M[row, 3*i + c] += a[j][c]
            rhs[row] = dG_f[(i, j)](lv); row += 1
        w = mp.lu_solve(M, rhs)
        return mp.sqrt(sum(w[q]**2 for q in range(9)))
    def total_speed(lv):
        return mp.sqrt(sum(sum(mp.mpf(dvh_f[i](lv)[c])**2 for c in range(3)) for i in range(3)))
    return quotient_speed, total_speed
def series_coeffs(speed, hs=('0.02', '0.01', '0.005')):
    """s(l) = s2 l^2 + s4 l^4 + ...  from  ds/dl = 2 s2 l + 4 s4 l^3 + ...:  three-point fit on small lambda"""
    vals = {}
    for h in (mp.mpf(hh) for hh in hs):
        vals[h] = speed(h)/(2*h)                                # = s2 + 2 s4 h^2 + O(h^4)
    h1, h2, h3 = sorted(vals)
    s2 = (4*vals[h1] - vals[h2])/3                              # h1 = h2/2: eliminate h^2 (crude); refine below
    # fit s2 + 2 s4 h^2 + 3 s6 h^4 through the three points (exact for the truncated series)
    A_ = mp.matrix([[1, 2*h**2, 3*h**4] for h in (h1, h2, h3)]); c = mp.lu_solve(A_, mp.matrix([vals[h] for h in (h1, h2, h3)]))
    return c[0], c[1]
GRAMS = {'N1': (sp.Rational(1, 4), -sp.Rational(1, 3), sp.Rational(1, 5)), 'N3': (sp.Rational(1, 20), -sp.Rational(1, 30), sp.Rational(1, 50))}
rows = {}
for nm, g in GRAMS.items():
    for layer in ('A', 'C'):
        qs, ts = lift_norm_fn(g, layer)
        s2, s4 = series_coeffs(qs); c2 = -(s4/s2 + mp.mpf(2)/3); per = 3 + 2*c2
        t2, t4 = series_coeffs(ts); c2t = -(t4/t2 + mp.mpf(2)/3)
        # exactness probe: is s(l)/tanh^2 l constant?  compare the quotient arc length at two lambdas
        sfun = lambda lv: mp.quad(qs, [0, lv])
        ratio1, ratio2 = sfun(mp.mpf('0.5'))/mp.tanh(mp.mpf('0.5'))**2, sfun(mp.mpf('1.5'))/mp.tanh(mp.mpf('1.5'))**2
        rows[(nm, layer)] = dict(s2=s2, s4=s4, c2=c2, per=per, c2t=c2t, r1=ratio1, r2=ratio2)
        print(f"  {nm} layer {layer}: quotient metric  s2 = {mp.nstr(s2, 8)}  s4/s2 = {mp.nstr(s4/s2, 8)}  (need -0.6667)  ->  c2 = {mp.nstr(c2, 6)}  perihelion coeff {mp.nstr(per, 6)}"
              f"   | full (S^2)^3 lift: c2 = {mp.nstr(c2t, 6)}   | s/tanh^2 at l = 0.5, 1.5: {mp.nstr(ratio1, 6)}, {mp.nstr(ratio2, 6)}")
# the fork, decided by the numbers (declared: DERIVED if |c_2| < 1e-4 everywhere; NOT DERIVED if |c_2| > 1e-2 everywhere)
derived = all(abs(v['c2']) < mp.mpf(10)**-4 for v in rows.values()); notder = all(abs(v['c2']) > mp.mpf('0.01') for v in rows.values())
check("D-3a the fit is converged: refitting s2 on a coarser triple of lambdas (0.04, 0.02, 0.01) changes it by less than 1e-8 (every frame, layer)",
      all(abs(series_coeffs(lift_norm_fn(g, ly)[0], hs=('0.04', '0.02', '0.01'))[0] - rows[(nm, ly)]['s2']) < mp.mpf(10)**-8 for nm, g in GRAMS.items() for ly in ('A', 'C')))
check("D-3b the fork is decided unambiguously across frames and layers (all DERIVED or all NOT DERIVED): branch = " +
      ("DERIVED (c_2 = 0): the frozen metric's harmonic-map energy reproduces KIN-2a" if derived else
       "NOT DERIVED: c_2 != 0, perihelion coefficients " + ", ".join(f"{nm}/{ly}: {mp.nstr(v['per'], 5)}" for (nm, ly), v in rows.items())),
      derived or notder)
check("D-3c the arc length of the frozen metric along the family is (or is not) proportional to tanh^2 l: ratios at l = 0.5 and 1.5 "
      + ("agree" if all(abs(v['r1'] - v['r2']) < mp.mpf(10)**-6 for v in rows.values()) else "differ") + " -- consistent with D-3b",
      (all(abs(v['r1'] - v['r2']) < mp.mpf(10)**-6 for v in rows.values())) == derived)
print("  ARGUMENT (not counted): D-1 + D-2 + D-3 together say the dynamics tier's static energy, if gradient-quadratic, must be")
print("  the Dirichlet energy of a function with a pole at the presented volume's zero -- and the only click-invariant with that")
print("  property among the objects on the table is the reciprocal presented volume Delta/det G'.  The frozen state metric is")
print("  a metric on DIRECTIONS; it cannot see the volume's reciprocal, and the numbers say it does not.")

print(); n_ = sum(CH); print(f"RESULT: {n_}/{len(CH)} checks passed.")
if n_ != len(CH): print("VERDICT: check failures above."); sys.exit(1)
if derived:
    print("VERDICT (THM-J, item 0b): DERIVED.  The harmonic-map energy of the frozen elliptope metric reproduces KIN-2a: the field")
    print("  law follows from the bare tier's own geometry.  Conditional on the harmonic-map form of the static energy.")
else:
    print("VERDICT (THM-J, item 0b): NOT DERIVED.  The model's frozen state-space geometry does not produce KIN-2a as a harmonic map.")
    print("  What IS established: any gradient-quadratic static energy makes its own arc length harmonic (D-1), so the energy that")
    print("  produces KIN-2a is exactly the gradient energy of Delta/det G' -- a reciprocal, not a polynomial, of the volume invariant")
    print("  (D-2).  The model has a Newton's law to declare, in one variable: E_static = int |grad(Delta/det G')|^2, sourced by the")
    print("  blind-mass measure (KIN-2a'').  Not an Einstein yet.")
print("COMPARISON STAGE: the Dirichlet energy of the metric N^2 = 1 - 2 Phi/c^2 is the Newtonian field energy; the elliptope")
print("  quotient metric is the shape-space metric of three unit vectors (mechanical connection, Littlejohn-Reinsch); harmonic maps")
print("  from R^3 into a target with the round quotient metric are sigma-model configurations, and sigma models do not give 1/r.")
sys.exit(0)
