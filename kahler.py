#!/usr/bin/env python3
# =============================================================================
# THE KAHLER DIAL - g and omega as two faces of one tensor; J rotates the doors
# 2026-08-28.  Will's synthesis; verification + the e^{J theta} interpolation.
# Frame space (S^2)^3 = (CP^1)^3 is Kahler; J_a v = a x v (the plane P_a = I a).
# Door C (first order): da_i = a_i x grad_i V     (precession; omega-face)
# Door R (first order): da_i = -P_i grad_i V      (descent;    g-face)
# Dial:  da_i = sin(th) a_i x grad_i V - cos(th) P_i grad_i V  = -e^{J th} grad V
# =============================================================================
import sympy as sp
import math

CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" - {n}" if n else ""), flush=True)

print("=" * 78); print("PART 1 - THEOREMS (symbolic)"); print("=" * 78)
a = sp.Matrix(sp.symbols('a1:4')); x = sp.Matrix(sp.symbols('x1:4'))
y = sp.Matrix(sp.symbols('y1:4')); s = a.dot(a)
u = x - (a.dot(x)/s)*a                     # generic tangent vectors at a
v = y - (a.dot(y)/s)*a
check("K-1  BAC-CAB, exact polynomial: a x (a x w) = (a.w)a - (a.a)w",
      sp.simplify(sp.expand(a.cross(a.cross(x)) - (a.dot(x)*a - a.dot(a)*x))) == sp.zeros(3, 1),
      "with |a|=1 and w tangent: J^2 = -1 on the tangent space [corollary]")
check("K-2  Lagrange, exact polynomial: (a x w).(a x z) = (a.a)(w.z) - (a.w)(a.z)",
      sp.simplify(sp.expand((a.cross(x)).dot(a.cross(y)) - (a.dot(a)*x.dot(y) - a.dot(x)*a.dot(y)))) == 0,
      "with |a|=1, tangent w,z: g(Jw,Jz) = g(w,z) - J is g-orthogonal [corollary]")
check("K-3  omega(u,v) := g(Ju,v) = a.(u x v): antisymmetric - the area form",
      sp.simplify(sp.expand((a.cross(u)).dot(v) - a.dot(u.cross(v)))) == 0 and
      sp.simplify(sp.expand((a.cross(u)).dot(v) + (a.cross(v)).dot(u))) == 0,
      "one Hermitian tensor h = g + i*omega; the two doors are its two faces")

# flows with a generic click-invariant V: grad_i V = sum_j c_ij a_j (chain rule)
A1 = sp.Matrix(sp.symbols('p1:4')); A2 = sp.Matrix(sp.symbols('q1:4')); A3 = sp.Matrix(sp.symbols('r1:4'))
c12, c13, c23 = sp.symbols('c12 c13 c23')
G1 = c12*A2 + c13*A3; G2 = c12*A1 + c23*A3; G3 = c13*A1 + c23*A2
prec = [A1.cross(G1), A2.cross(G2), A3.cross(G3)]
Vdot_prec = sp.simplify(G1.dot(prec[0]) + G2.dot(prec[1]) + G3.dot(prec[2]))
check("K-4a Door C (precession) conserves V EXACTLY: dV/dt = Sum grad.(a x grad) = 0",
      Vdot_prec == 0, "rotates along the level sets - the omega-face signature")
check("K-4b Door C conserves the total axis vector: d/dt Sum a_i = Sum a_i x grad_i V = 0",
      sp.simplify(prec[0] + prec[1] + prec[2]) == sp.zeros(3, 1),
      "the spin-picture conserved momentum, same cancellation as S1")
th = sp.Symbol('theta', real=True)
def P(w, aa): return w - (aa.dot(w)/aa.dot(aa))*aa
dial = [sp.sin(th)*A1.cross(G1) - sp.cos(th)*P(G1, A1),
        sp.sin(th)*A2.cross(G2) - sp.cos(th)*P(G2, A2),
        sp.sin(th)*A3.cross(G3) - sp.cos(th)*P(G3, A3)]
Vdot = sp.simplify(G1.dot(dial[0]) + G2.dot(dial[1]) + G3.dot(dial[2]))
target = -sp.cos(th)*(P(G1, A1).dot(P(G1, A1)) + P(G2, A2).dot(P(G2, A2)) + P(G3, A3).dot(P(G3, A3)))
check("K-5  the DIAL: dV/dt = -cos(theta) * Sum |P grad_i V|^2, exactly",
      sp.simplify(sp.expand(Vdot - target)) == 0,
      "descent rate = cos(theta); pure rotation at pi/2; e^{J theta} interpolation")
print("      K-0 [identity, elementary, uncounted]: J's generator at axis a is the")
print("      bivector P_a = I a  (sec 12.2's plane); quarter-turn exp(-P_a pi/4) acts")
print("      as a x (.) - the operator that turns omega into g is the axis's own plane.")

print(); print("=" * 78)
print("PART 2 - THE DIAL, RUN  [numeric evidence, RK4]")
print("=" * 78)
def dot(u,v): return sum(x*y for x,y in zip(u,v))
def sub(u,v): return [x-y for x,y in zip(u,v)]
def add(u,v): return [x+y for x,y in zip(u,v)]
def scl(s,u): return [s*x for x in u]
def crs(u,v): return [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]]
def nrm(u):
    n = math.sqrt(dot(u,u)); return [x/n for x in u]
def gam(A): return (dot(A[0],A[1]), dot(A[0],A[2]), dot(A[1],A[2]))
def Dv(g): return 1 - g[0]**2 - g[1]**2 - g[2]**2 + 2*g[0]*g[1]*g[2]
def Vv(A): return -0.5*math.log(Dv(gam(A)))
def grads(A):
    g = gam(A); D = Dv(g)
    c = (-(-2*g[0]+2*g[1]*g[2])/(2*D), -(-2*g[1]+2*g[0]*g[2])/(2*D), -(-2*g[2]+2*g[0]*g[1])/(2*D))
    return [add(scl(c[0],A[1]), scl(c[1],A[2])),
            add(scl(c[0],A[0]), scl(c[2],A[2])),
            add(scl(c[1],A[0]), scl(c[2],A[1]))]
def rhs(A, thv):
    Gr = grads(A); out = []
    for i in range(3):
        Pg = sub(Gr[i], scl(dot(Gr[i],A[i]), A[i]))
        out.append(sub(scl(math.sin(thv), crs(A[i],Gr[i])), scl(math.cos(thv), Pg)))
    return out
def run(thv, N=3000, h=0.002):
    A = [[1.0,0,0],[0.55, math.sqrt(1-0.3025),0],[-0.30,0.35,math.sqrt(1-0.09-0.1225)]]
    V0 = Vv(A); flips = 0; last = None
    for k in range(N):
        k1 = rhs(A,thv); A2=[nrm(add(A[i],scl(h/2,k1[i]))) for i in range(3)]
        k2 = rhs(A2,thv); A3=[nrm(add(A[i],scl(h/2,k2[i]))) for i in range(3)]
        k3 = rhs(A3,thv); A4=[nrm(add(A[i],scl(h,k3[i]))) for i in range(3)]
        k4 = rhs(A4,thv)
        A = [nrm(add(A[i], scl(h/6, add(add(k1[i],scl(2,k2[i])), add(scl(2,k3[i]),k4[i]))))) for i in range(3)]
        g = gam(A)
        if last is not None and (g[0]-lastg)*lastd < 0: flips += 1
        if last is not None: lastd = g[0]-lastg
        else: lastd = 0.0
        lastg = g[0]; last = True
    return V0, Vv(A), flips
res = {t: run(t) for t in (0.0, math.pi/6, math.pi/3, math.pi/2)}
for t,(V0,V1,fl) in res.items():
    print(f"    theta = {t:6.3f}:  V {V0:.4f} -> {V1:.6f}   turning points of g12: {fl}")
check("K-6a theta = pi/2: V conserved (|dV| < 1e-8) - pure omega-face rotation",
      abs(res[math.pi/2][1]-res[math.pi/2][0]) < 1e-8)
check("K-6b V_end strictly increases with theta: descent rate = cos(theta)",
      res[0.0][1] < res[math.pi/6][1] < res[math.pi/3][1] < res[math.pi/2][1])
check("K-6c intermediate theta: many oscillations WHILE descending - EVANESCENCE",
      res[math.pi/3][2] > 10 and res[math.pi/3][1] < res[math.pi/2][1] - 1e-4,
      "damped rotation = complex frequency: the propagating/evanescent band realised")
print(f"\nRESULT: {sum(CH)}/{len(CH)} checks passed.")
print("""
ANSWER TO THE QUESTION: the operator turning omega into g is J, the complex
structure - concretely the axis's own plane P_a = I a (sec 12.2), applied to the
flow as e^{J theta}, i.e. a WICK ROTATION OF THE GENERATOR.  Door R is Door C in
imaginary time; LBL-1's thermal face IS the dissipative face; the band between
is a complex frequency - the model's own word for it is EVANESCENT.
WHO SETS theta: nothing internal can.  At theta = pi/2 the flow is norm- and
V-preserving; rotating the phase requires a symmetric positive part in the
generator, which is a coupling to unresolved degrees of freedom (a dissipator).
Closed structure => theta = pi/2 exactly, and every apparent descent is
coarse-graining.  The yield-point operator is therefore THE OPENING itself -
the same ruling as the two doors, now as a phase.  [status: the Kahler geometry
and the dial are theorems/verified; the physics assignments (quantum = omega
face, thermal/settling = g face) are conjecture-tier, labelled.]""")
