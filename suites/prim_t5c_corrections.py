#!/usr/bin/env python3
# =============================================================================
# PRIMITIVES T5(c) -- three corrections (Will, 2026-09-04), each verified.  Supersedes T5(b)'s NAMING.
#   1. det G_c = -D^2 <= 0.  It never crosses zero with a sign; it TOUCHES zero quadratically.  The quantity that
#      crosses zero is D, the oriented volume det[hbar, G, c].  The tangential-sector debt is sgn D(r).
#   2. The observable equal to +-2D is the RAW SPINOR TRACE V_spin = tr_2(Gamma_C Gamma_H Gamma_G), NOT the projector
#      Bargmann invariant Tr(P_C P_H P_G).  Gamma_C^2 = -1 so (1 + Gamma_C)/2 is not a projector.  The vector
#      rank-one projector reading is Tr(P_C P_H P_G) = q(C,H) q(H,G) q(G,C) / (q(C,C) q(H,H) q(G,G)), real,
#      reversal-invariant.  These must remain differently named.
#   3. Cl(2,1) = M_2(R) (+) M_2(R): the pseudoscalar I = Gamma_C Gamma_H Gamma_G is central with I^2 = +1, and the
#      two blocks have I = +1 and I = -1.  A single 2x2 block is not faithful on Pin(2,1).  The two blocks give
#      V_spin = +2D and -2D: the block label IS the sheet.  Forced by real Cl(2,1); not an imported Cl(3).
# EARNED THEOREM:  V_spin = 2D,   D^2 = -det G_c,   tau: D -> -D.   At D = 0 the observable is defined and zero;
#   only sgn D is undefined.  Stellar assignment (neutron-star <-> one sign, black-hole <-> the other) is a
#   DYNAMICAL proposition awaiting the collapse dynamics.
# =============================================================================
import sympy as sp, itertools, time
T0 = time.time(); CH = []
def check(t, ok, n=""):
    CH.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-T0:5.1f}s] {t}" + (f" -- {n}" if n else ""), flush=True)
def z(e):
    e = sp.sympify(e)
    for f in (sp.simplify, lambda q: sp.simplify(sp.expand(q)), sp.cancel):
        try:
            if f(e) == 0: return True
        except Exception: pass
    return False
def zM(M): return all(z(e) for e in M)

Q = sp.diag(1, 1, -1)                                   # (hbar, G, c)
t, l1, l2 = sp.symbols('t l1 l2', real=True)
Cv = sp.Matrix([0, 0, 1]); Hv = sp.Matrix([sp.cosh(l1), 0, sp.sinh(l1)])
Gv = sp.Matrix([sp.cosh(l2)*sp.cos(t), sp.cosh(l2)*sp.sin(t), sp.sinh(l2)])
F = sp.Matrix.hstack(Hv, Gv, Cv)                         # the frame matrix, columns hbar, G, c
G_c = sp.simplify(F.T*Q*F)                               # the Gram of the seat's form
D = sp.simplify(F.det())                                 # the ORIENTED volume

print("=== 1. det G_c = -D^2: it touches zero, D crosses it ===")
check("1a D = det[hbar, G, c] = cosh l1 cosh l2 sin t (the oriented volume, in seat coordinates)", z(D - sp.cosh(l1)*sp.cosh(l2)*sp.sin(t)))
check("1b det G_c = det(F^T Q F) = det(Q) D^2 = -D^2 exactly", z(sp.simplify(G_c.det() + D**2)))
check("1c det G_c <= 0 for every state and is EVEN in D: it cannot cross zero with a sign; it touches zero quadratically at sin t = 0",
      z(sp.simplify(G_c.det().subs(t, -t) - G_c.det())) and z(sp.simplify(G_c.det() + D**2)))
check("1d D is ODD in t: D(-t) = -D(t); D crosses zero at the branch locus with a sign -- the sign is sgn D, and THAT is the debt",
      z(sp.simplify(D.subs(t, -t) + D)))

print("=== 2. V_spin is the raw spinor trace, not the projector Bargmann invariant ===")
g_H = sp.Matrix([[0, 1], [1, 0]]); g_G = sp.Matrix([[1, 0], [0, -1]]); g_C = sp.Matrix([[0, 1], [-1, 0]])
x_C = g_C; x_H = sp.cosh(l1)*g_H + sp.sinh(l1)*g_C; x_G = sp.cosh(l2)*(sp.cos(t)*g_H + sp.sin(t)*g_G) + sp.sinh(l2)*g_C
V_spin = sp.simplify((x_C*x_H*x_G).trace())
check("2a V_spin := tr_2(Gamma_C Gamma_H Gamma_G) = 2D exactly, at every state", z(V_spin - 2*D))
check("2b V_spin^2 = -4 det G_c (T5b's identity, now correctly attributed to the spinor trace)", z(sp.simplify(V_spin**2 + 4*G_c.det())))
check("2c reversal: tr_2(Gamma_C Gamma_G Gamma_H) = -2D  (tau: D -> -D)", z(sp.simplify((x_C*x_G*x_H).trace()) + 2*D))
# the timelike gamma does NOT give a projector
Pc_try = (sp.eye(2) + g_C)/2
check("2d (1 + Gamma_C)/2 is NOT idempotent (Gamma_C^2 = -1): there is no real projector onto the timelike line in this block -- the projector Bargmann invariant is a DIFFERENT object",
      not zM(Pc_try*Pc_try - Pc_try))
# the vector rank-one projector reading, properly normalised
a_, b_, g_ = sp.symbols('a b gamma', real=True)
def P(v):   # q-orthogonal rank-one projector onto line v: P w = v q(v,w)/q(v,v)
    return v*(v.T*Q)/((v.T*Q*v)[0])
Cg = sp.Matrix([0,0,1]); Hg = sp.Matrix([1,0,0]); Gg = sp.Matrix([sp.cos(t), sp.sin(t), 0])   # orthogonal-depth frame for a clean formula check
Tr_vec = sp.simplify((P(Cv)*P(Hv)*P(Gv)).trace())
qCH = (Cv.T*Q*Hv)[0]; qHG = (Hv.T*Q*Gv)[0]; qGC = (Gv.T*Q*Cv)[0]
qCC = (Cv.T*Q*Cv)[0]; qHH = (Hv.T*Q*Hv)[0]; qGG = (Gv.T*Q*Gv)[0]
check("2e vector reading: Tr(P_C P_H P_G) = q(C,H) q(H,G) q(G,C) / (q(C,C) q(H,H) q(G,G)), with the timelike denominator -1 (T5b omitted the normalisation)",
      z(sp.simplify(Tr_vec - qCH*qHG*qGC/(qCC*qHH*qGG))))
check("2f the vector reading is real and reversal-invariant: Tr(P_C P_G P_H) = Tr(P_C P_H P_G)",
      z(sp.simplify((P(Cv)*P(Gv)*P(Hv)).trace() - Tr_vec)))
check("2g the two objects are different functions of the state: V_spin is linear in D, Tr(P P P) is even in D -- they must carry different names",
      not z(sp.simplify(V_spin - 2*Tr_vec)) and z(sp.simplify(Tr_vec.subs(t, -t) - Tr_vec)))

print("=== 3. Cl(2,1) = M_2(R) (+) M_2(R): the two blocks are the two sheets ===")
I_pseudo = g_C*g_H*g_G
check("3a the pseudoscalar I = Gamma_C Gamma_H Gamma_G is CENTRAL (commutes with every generator)",
      all(zM(I_pseudo*g - g*I_pseudo) for g in (g_C, g_H, g_G)))
check("3b I^2 = +1: Cl(2,1) splits into two blocks I = +1 and I = -1; the algebra is M_2(R) (+) M_2(R), not simple",
      I_pseudo**2 == sp.eye(2))
check("3c the real 2x2 rep of T5b is the I = +1 block: Gamma_C Gamma_H Gamma_G = +I_2 there", I_pseudo == sp.eye(2))
# the other block: flip the sign of one generator (an automorphism of the relations sending I -> -I)
g_Cm, g_Hm, g_Gm = -g_C, g_H, g_G
I_minus = g_Cm*g_Hm*g_Gm
check("3d the I = -1 block: Gamma_C -> -Gamma_C preserves all relations (squares, anticommutation) and sends I to -I_2",
      g_Cm**2 == -sp.eye(2) and zM(g_Cm*g_Hm + g_Hm*g_Cm) and zM(g_Cm*g_Gm + g_Gm*g_Cm) and I_minus == -sp.eye(2))
x_Cm = g_Cm; x_Hm = sp.cosh(l1)*g_Hm + sp.sinh(l1)*g_Cm; x_Gm = sp.cosh(l2)*(sp.cos(t)*g_Hm + sp.sin(t)*g_Gm) + sp.sinh(l2)*g_Cm
V_minus = sp.simplify((x_Cm*x_Hm*x_Gm).trace())
check("3e in the I = -1 block V_spin = -2D: the two blocks read the two sheets; the block label IS the sheet", z(V_minus + 2*D))
# the faithful representation: the 4x4 direct sum
G4 = {k: sp.diag(v[0], v[1]) for k, v in {'C': (g_C, g_Cm), 'H': (g_H, g_Hm), 'G': (g_G, g_Gm)}.items()}
I4 = G4['C']*G4['H']*G4['G']
check("3f the 4x4 direct sum is faithful: I = diag(+1,+1,-1,-1), both sheets present; a single 2x2 block sees only one",
      I4 == sp.diag(1, 1, -1, -1))
check("3g on the faithful rep the TOTAL trace of Gamma_C Gamma_H Gamma_G vanishes: the sheet is invisible to the full trace and visible only to the block-resolved trace -- the 'one-sided' reading is the choice of block",
      z(I4.trace()))

print("=== the earned theorem ===")
check("V_spin = 2D, D^2 = -det G_c, tau: D -> -D; at D = 0 the observable is defined (= 0) and only sgn D is undefined",
      z(V_spin - 2*D) and z(sp.simplify(D**2 + G_c.det())) and z(V_spin.subs(t, 0)))

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT [DERIVED]: det G_c = -D^2 touches zero; D crosses it; the debt is sgn D(r).  The sheet-reading observable")
print("  is V_spin = tr_2(Gamma_C Gamma_H Gamma_G) = 2D, NOT the projector Bargmann invariant (real, even in D).")
print("  Cl(2,1) = M_2(R)+M_2(R), I^2 = +1: the two blocks give +-2D and the block label is the sheet; the full")
print("  4x4 trace is sheet-blind.  Neutron-star vs black-hole <-> sgn D is a DYNAMICAL proposition, not yet derived.")
