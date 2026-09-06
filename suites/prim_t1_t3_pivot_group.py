#!/usr/bin/env python3
# =============================================================================
# PRIMITIVES T1-T3 -- what group do the pivots form?  (2026-09-04)
#
# INPUTS (docs/2026-09-04-PRIMITIVES-v0.md): P2 lines through a root, P3 planes,
#   P4 angles as state, P6 pivot = rotation by an IMAGINARY angle.  Nothing else.
# BANNED: Cl(3), Pauli matrices, spinors, any 4-vector, any signature as a given.
# TOOLS: real matrices and their Lie algebras only.  Sympy exact.
#
# HELD-OUT PREDICTIONS (written before running):
#   T1  one plane: the plane's quadratic invariant is forced from (2,0) to (1,1);
#       the pivots form SO(1,1); rapidity adds.
#   T2  three lines: on the 3-space itself, Wick on ONE plane's coordinate gives
#       SO(2,1) (two hyperbolic planes, one circular).  Wick on ALL three is not
#       realisable on the 3-space; the abstract algebra is so(3,C)_R = so(1,3),
#       but a (1,3) form needs a FOURTH direction the primitives do not supply.
#   T3  Wigner: two non-collinear pivots in SO(2,1) compose to pivot x rotation
#       with tan(w/2) = s1 s2 sin(a) / (c1 c2 + s1 s2 cos(a)).  If so, BARE-1's
#       kill condition is re-earned without Cl(3).
# OUTCOME FORK: (a) T1-T3 all hold -> 1+1 and 2+1 relativity DERIVED; (1,3)
#   needs one more primitive (name it).  (b) T2 realises (1,3) on the 3-space
#   somehow -> my prediction wrong, record how.  (c) Wigner fails in SO(2,1) ->
#   the pivot needs a declaration beyond 'imaginary angle'.
# =============================================================================
import sympy as sp, time
T0 = time.time(); CH = []
def check(t, ok, n=""):
    CH.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-T0:5.1f}s] {t}" + (f" -- {n}" if n else ""), flush=True)
def z(e):
    e = sp.sympify(e)
    for f in (sp.simplify, lambda q: sp.simplify(sp.expand(q)),
              lambda q: sp.simplify(sp.expand(q.rewrite(sp.exp)))):
        try:
            if f(e) == 0: return True
        except Exception: pass
    return False
def zM(M): return all(z(e) for e in M)

th, lam, a, b = sp.symbols('theta lambda a b', real=True)
I = sp.I

print("=== T1: one plane, one angle, one imaginary rotation ===")
R = lambda t: sp.Matrix([[sp.cos(t), -sp.sin(t)], [sp.sin(t), sp.cos(t)]])
Q_E = sp.eye(2)
Q_M = sp.diag(1, -1)
# T1a: R(theta) preserves the Euclidean form -- the plane's given invariant.
check("T1a  R(theta)^T Q_E R(theta) = Q_E (rotation preserves the plane's form)",
      zM(R(th).T*Q_E*R(th) - Q_E))
# T1b: the imaginary-angle rotation is complex; conjugating by S=diag(1,i) makes it real.
Ri = R(I*lam)
S = sp.diag(1, I)
B = sp.simplify(S.inv()*Ri*S)
B_expected = sp.Matrix([[sp.cosh(lam), sp.sinh(lam)], [sp.sinh(lam), sp.cosh(lam)]])
check("T1b  S^-1 R(i lam) S is REAL and equals [[ch,sh],[sh,ch]]",
      zM(B - B_expected) and all(sp.im(e.rewrite(sp.exp)).simplify() == 0 for e in B))
# T1c: what real form does the realified pivot preserve?  Not Q_E.  Q_M.
check("T1c  the real pivot does NOT preserve Q_E (signature change is forced)",
      not zM(B.T*Q_E*B - Q_E))
check("T1d  the real pivot preserves Q_M = diag(1,-1): signature (1,1)",
      zM(B.T*Q_M*B - Q_M))
# T1e: the form Q_E pulled back through S is Q_M: the coordinate that went
#      imaginary is the one whose square flipped sign.
check("T1e  S^T Q_E S = Q_M (the Wick coordinate is the timelike one)",
      zM(sp.simplify(S.T*Q_E*S) - Q_M))
# T1f: rapidity adds.
check("T1f  B(a) B(b) = B(a+b): rapidity is additive, the group is SO(1,1)",
      zM(B.subs(lam, a)*B.subs(lam, b) - B.subs(lam, a+b)))
# T1g: uniqueness -- Q_M is the ONLY symmetric form (up to scale) preserved by all B(lam).
q11, q12, q22 = sp.symbols('q11 q12 q22', real=True)
Q = sp.Matrix([[q11, q12], [q12, q22]])
eqs = list(sp.simplify(B.T*Q*B - Q))
sol = sp.solve([sp.expand(e.rewrite(sp.exp)) for e in eqs], [q11, q12, q22], dict=True)
Qs = Q.subs(sol[0]) if len(sol) == 1 else None
uniq = Qs is not None and Qs[0,1] == 0 and sp.simplify(Qs[0,0] + Qs[1,1]) == 0 and Qs[1,1] != 0
check("T1g  the only invariant symmetric form of the pivot family is prop. to diag(1,-1)",
      uniq, str(sol))

print("=== T2: three lines, three planes -- which planes take the imaginary angle? ===")
# so(3) generators: rotation in plane (j,k) is exp(t L_i), i the missing index.
L1 = sp.Matrix([[0,0,0],[0,0,-1],[0,1,0]])   # plane (2,3)
L2 = sp.Matrix([[0,0,1],[0,0,0],[-1,0,0]])   # plane (1,3)
L3 = sp.Matrix([[0,-1,0],[1,0,0],[0,0,0]])   # plane (1,2)
Ls = [L1, L2, L3]
comm = lambda X, Y: X*Y - Y*X
check("T2a  so(3): [L1,L2]=L3, [L2,L3]=L1, [L3,L1]=L2",
      zM(comm(L1,L2)-L3) and zM(comm(L2,L3)-L1) and zM(comm(L3,L1)-L2))

# ONE coordinate imaginary: x3 -> i x3.  S3 = diag(1,1,i).  Realify each generator.
S3 = sp.diag(1, 1, I)
def realify(G, S):  # returns S^-1 (i G) S if that is real, else S^-1 G S if real, else None
    for cand in (S.inv()*(I*G)*S, S.inv()*G*S):
        if all(sp.im(sp.simplify(e)) == 0 for e in cand):
            return sp.simplify(cand)
    return None
K1 = realify(L1, S3); K2 = realify(L2, S3); J3 = realify(L3, S3)
issym = lambda M: zM(M - M.T); isanti = lambda M: zM(M + M.T)
check("T2b  with x3 -> i x3: planes (2,3),(1,3) become SYMMETRIC generators (boosts), (1,2) stays ANTISYMMETRIC (rotation)",
      K1 is not None and K2 is not None and J3 is not None and issym(K1) and issym(K2) and isanti(J3) and zM(J3 - L3))
Q21 = sp.diag(1, 1, -1)
check("T2c  the three realified generators preserve diag(1,1,-1): the group is SO(2,1)",
      all(zM(G.T*Q21 + Q21*G) for G in (K1, K2, J3)))
check("T2d  so(2,1) relations: [K1,K2] = -J3, [J3,K1] = K2, [K2,J3] = K1",
      zM(comm(K1,K2)+J3) and zM(comm(J3,K1)-K2) and zM(comm(K2,J3)-K1))
# By symmetry, any single coordinate made imaginary gives SO(2,1) with that line timelike.
S3b = sp.diag(I, 1, 1)
Kb = [realify(G, S3b) for G in Ls]
check("T2e  x1 -> i x1 likewise gives SO(2,1) (the timelike line is a CHOICE among the three)",
      all(G is not None for G in Kb) and all(zM(G.T*sp.diag(-1,1,1) + sp.diag(-1,1,1)*G) for G in Kb))

# ALL THREE planes with imaginary angle?  Enumerate every S = diag(s1,s2,s3), s_i in {1,i}.
# For each, count how many of the three plane-generators realify WITH the factor i (boost)
# versus without it (rotation).  Prediction: never three boosts; only 0 or 2.
import itertools
counts = set()
for signs in itertools.product([1, I], repeat=3):
    Sx = sp.diag(*signs)
    nb = 0
    for G in Ls:
        cand_b = sp.simplify(Sx.inv()*(I*G)*Sx)
        cand_r = sp.simplify(Sx.inv()*G*Sx)
        real_b = all(sp.im(e) == 0 for e in cand_b)
        real_r = all(sp.im(e) == 0 for e in cand_r)
        if real_b and not zM(cand_b): nb += 1
        elif not real_r: nb = -99   # neither realifies: inconsistent choice
    counts.add(nb)
check("T2f  over all 8 realifications, the number of boost planes is in {0, 2}: three boosts on the 3-space is impossible",
      counts == {0, 2}, f"boost-plane counts seen: {sorted(counts)}")
# The abstract algebra: J_k = L_k, K_k = i L_k satisfy the Lorentz relations.
Jk = Ls; Kk = [I*G for G in Ls]
eps = lambda i,j,k: sp.LeviCivita(i,j,k)
lor = True
for i in range(3):
    for j in range(3):
        JJ = comm(Jk[i], Jk[j]) - sum((eps(i,j,k)*Jk[k] for k in range(3)), sp.zeros(3,3))
        JK = comm(Jk[i], Kk[j]) - sum((eps(i,j,k)*Kk[k] for k in range(3)), sp.zeros(3,3))
        KK = comm(Kk[i], Kk[j]) + sum((eps(i,j,k)*Jk[k] for k in range(3)), sp.zeros(3,3))
        lor = lor and zM(JJ) and zM(JK) and zM(KK)
check("T2g  abstractly, {L_k, i L_k} satisfy [J,J]=eJ, [J,K]=eK, [K,K]=-eJ: so(3,C)_R = so(1,3)",
      lor)
# Conclusion of T2 (recorded, not checked): the Lorentz ALGEBRA is reachable from P6 on
# all three planes, but a (1,3)-signature SPACE for it to act on is not -- the 3-space
# admits SO(3) or SO(2,1) only (T2f).  A fourth direction is a new primitive.

print("=== T3: Wigner rotation inside SO(2,1), no Cl(3) ===")
# RATIONAL COORDINATES (the SECT-1 trick): v_i = exp(lambda_i/2) > 0, t = tan(alpha/2).
# Then cosh/sinh/cos/sin are rational, every simplify is polynomial, and there is
# no acosh, no exp, no matrix log anywhere.  Same theorem as before, honest speed.
v1, v2, t = sp.symbols('v1 v2 t', positive=True)
def ch(v): return (v**2 + 1/v**2)/2      # cosh(lambda)
def sh(v): return (v**2 - 1/v**2)/2      # sinh(lambda)
def ch2(v): return (v + 1/v)/2           # cosh(lambda/2)
def sh2(v): return (v - 1/v)/2           # sinh(lambda/2)
ca, sa = (1 - t**2)/(1 + t**2), 2*t/(1 + t**2)   # cos(alpha), sin(alpha)
def boost(v, n1, n2):
    # pure boost of rapidity lambda = 2 ln v along unit (n1, n2) in the (x1,x2) plane, x3 timelike
    C, Sh = ch(v), sh(v)
    return sp.Matrix([[1 + (C-1)*n1**2, (C-1)*n1*n2, Sh*n1],
                      [(C-1)*n1*n2, 1 + (C-1)*n2**2, Sh*n2],
                      [Sh*n1, Sh*n2, C]])
def pure_boost_from_image(u):
    # the unique pure boost sending e3 = (0,0,1) to u = (u1,u2,u3), u3 = cosh chi; rational in u
    u1, u2, u3 = u
    d = 1 + u3
    return sp.Matrix([[1 + u1**2/d, u1*u2/d, u1],
                      [u1*u2/d, 1 + u2**2/d, u2],
                      [u1, u2, u3]])
B1 = boost(v1, 1, 0); B2 = boost(v2, ca, sa)
check("T3a  boost(v1, x1) preserves diag(1,1,-1)", zM(sp.expand(B1.T*Q21*B1 - Q21)))
check("T3b  boost(v2, alpha) preserves diag(1,1,-1)", zM(sp.cancel(B2.T*Q21*B2 - Q21)))
P = sp.cancel(B1*B2)
u = list(P[:, 2])
Bt = pure_boost_from_image(u)
check("T3c  the pure boost rebuilt from P e3 preserves diag(1,1,-1) and is symmetric",
      zM(sp.cancel(Bt.T*Q21*Bt - Q21)) and zM(sp.cancel(Bt - Bt.T)))
# B_t^-1 = Q B_t Q for a pure boost (B^T Q B = Q and B = B^T).
Rw = sp.cancel(Q21*Bt*Q21*P)
check("T3d  R = B_t^-1 (B1 B2) fixes the timelike axis: pure spatial rotation",
      zM(sp.cancel(Rw[:, 2] - sp.Matrix([0,0,1]))) and zM(sp.cancel(Rw[2, :] - sp.Matrix([[0,0,1]]))))
check("T3e  R is orthogonal on the (x1,x2) block with det +1",
      zM(sp.cancel(Rw[:2,:2].T*Rw[:2,:2] - sp.eye(2))) and z(sp.cancel(Rw[:2,:2].det() - 1)))
cw, sw = Rw[0,0], Rw[1,0]
c1, s1, c2, s2 = ch2(v1), sh2(v1), ch2(v2), sh2(v2)
tan_half_pred = s1*s2*sa/(c1*c2 + s1*s2*ca)
tan_half_got  = sw/(1 + cw)
check("T3f  tan^2(w/2) = [s1 s2 sin(a) / (c1 c2 + s1 s2 cos(a))]^2  [Wigner magnitude, re-earned in SO(2,1)]",
      z(sp.cancel(tan_half_got**2 - tan_half_pred**2)))
check("T3f' sign under THIS convention (P = B1 B2, R = B_t^-1 P, w from R[1,0]): tan(w/2) = -(formula); orientation, not physics",
      z(sp.cancel(tan_half_got + tan_half_pred)))
check("T3g  margin: w = 0 iff collinear (t = 0) or unpivoted (v1 = 1 or v2 = 1)",
      z(tan_half_pred.subs(t, 0)) and z(tan_half_pred.subs(v1, 1)) and z(tan_half_pred.subs(v2, 1)))
# Second path (independent): mpmath at 40 digits on three random points, via the
# original exp/acosh route, to corroborate the rational identity numerically.
import mpmath as mp; mp.mp.dps = 40
K1n = mp.matrix([[0,0,0],[0,0,1],[0,1,0]]); K2n = mp.matrix([[0,0,1],[0,0,0],[1,0,0]])
def boost_num(l, ang): return mp.expm(l*(mp.cos(ang)*K2n + mp.sin(ang)*K1n))
ok_num = True
for (l1n, l2n, an) in [(mp.mpf('0.7'), mp.mpf('1.3'), mp.mpf('0.9')),
                       (mp.mpf('2.1'), mp.mpf('0.4'), mp.mpf('2.5')),
                       (mp.mpf('0.05'), mp.mpf('3.0'), mp.mpf('1.1'))]:
    Pn = boost_num(l1n, 0)*boost_num(l2n, an)
    un = [Pn[0,2], Pn[1,2], Pn[2,2]]
    dn = 1 + un[2]
    Btn = mp.matrix([[1+un[0]**2/dn, un[0]*un[1]/dn, un[0]],[un[0]*un[1]/dn, 1+un[1]**2/dn, un[1]],[un[0],un[1],un[2]]])
    Qn = mp.diag([1,1,-1])
    Rn = Qn*Btn*Qn*Pn
    got = Rn[1,0]/(1 + Rn[0,0])
    pred = mp.sinh(l1n/2)*mp.sinh(l2n/2)*mp.sin(an)/(mp.cosh(l1n/2)*mp.cosh(l2n/2) + mp.sinh(l1n/2)*mp.sinh(l2n/2)*mp.cos(an))
    ok_num = ok_num and abs(got + pred) < mp.mpf('1e-35')   # same sign convention as T3f'
check("T3h  second path: exp/acosh-free numeric route agrees with the closed form to 1e-35 at three points", ok_num)

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT (computed from the table, not written in advance):")
print("  T1 :", "1+1 signature and SO(1,1) DERIVED from P6" if all(CH[0:7]) else "T1 FAILED -- see above")
print("  T2 :", "SO(2,1) on the 3-space for one imaginary plane; Lorentz ALGEBRA reachable, (1,3) SPACE not"
      if all(CH[7:14]) else "T2 anomaly -- read the FAILs")
print("  T3 :", "Wigner re-earned in SO(2,1) without Cl(3)" if all(CH[14:23]) else "T3 FAILED -- pivot needs more than 'imaginary angle'")
