#!/usr/bin/env python3
# =============================================================================
# PRIMITIVES T8(a) -- the fourth direction is the operation that exchanges the two Cl(2,1) blocks.   (Will, 2026-09-04)
#
#   Cl(2,1) = M_2(R) (+) M_2(R) with central pseudoscalar I, I^2 = +1 (T5c).  The two blocks are the two sheets.
#   Claim: a fourth generator Gamma_4 anticommuting with Gamma_C, Gamma_H, Gamma_G satisfies Gamma_4 I = -I Gamma_4,
#   hence maps the I = +1 block to the I = -1 block.  Without Gamma_4 the sheet is a SUPERSELECTION label (no element
#   of Cl(2,1) mixes the blocks); with Gamma_4 the blocks can superpose coherently and an arriving configuration can
#   resolve the superposition.  This derives the algebraic JOB of the fourth object without assuming a fourth
#   spacetime coordinate.
# CHECKS:
#   a1  in the faithful 4x4 rep every element of Cl(2,1) is block-diagonal: commutes with I_4 = diag(1,1,-1,-1).
#       SUPERSELECTION: no matrix element connects the sheets.
#   a2  Gamma_4 anticommuting with all three generators exists, is OFF-diagonal, and Gamma_4 I_4 = -I_4 Gamma_4.
#   a3  the general solution: Gamma_4 = [[0, alpha Gamma_C], [beta Gamma_C, 0]]; Gamma_4^2 = -alpha beta I_4.
#       Gamma_4^2 = +1 (alpha beta = -1): signature (+,+,+,-) = Cl(3,1).  Gamma_4^2 = -1 (alpha beta = +1): (+,+,-,-) = Cl(2,2).
#   a4  with Gamma_4 adjoined the pseudoscalar I is no longer central; the centre of the generated algebra is the scalars:
#       the algebra is SIMPLE (M_4(R) in both signatures), the sectors are gone.
#   a5  superposition made concrete: for spinors psi_+ (block +) and psi_- (block -), every Cl(2,1) bilinear
#       <psi_+ | X | psi_-> vanishes -- no interference -- while <psi_+ | Gamma_4 | psi_-> does not.
#   a6  the two blocks are INEQUIVALENT irreps (no intertwiner); Gamma_4 conjugation is the grade involution x -> -x, made inner.
# SIGNATURE OF THE FOURTH DIRECTION: left OPEN here.  Cl(3,1) gives 3+1 with c timelike; Cl(2,2) gives 2+2.
#   The physical reading (the fourth direction is spatial) requires Gamma_4^2 = +1.  Kill: if the dynamics forces
#   Gamma_4^2 = -1, the fourth direction is timelike and the 3+1 reading is dead.
# =============================================================================
import sympy as sp, itertools, time
T0 = time.time(); CH = []
def check(t, ok, n=""):
    CH.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-T0:5.1f}s] {t}" + (f" -- {n}" if n else ""), flush=True)
def zM(M): return all(sp.simplify(e) == 0 for e in M)

g_H = sp.Matrix([[0, 1], [1, 0]]); g_G = sp.Matrix([[1, 0], [0, -1]]); g_C = sp.Matrix([[0, 1], [-1, 0]])
Z = sp.zeros(2)
blk = lambda A, B: sp.Matrix(sp.BlockMatrix([[A, Z], [Z, B]]))
off = lambda A, B: sp.Matrix(sp.BlockMatrix([[Z, A], [B, Z]]))
# faithful 4x4 rep: block + uses (g_C, g_H, g_G); block - uses (-g_C, g_H, g_G)  (T5c 3d)
G_C = blk(g_C, -g_C); G_H = blk(g_H, g_H); G_G = blk(g_G, g_G)
I4 = G_C*G_H*G_G
gens = {'C': G_C, 'H': G_H, 'G': G_G}

print("=== a1: superselection inside Cl(2,1) ===")
check("I_4 = Gamma_C Gamma_H Gamma_G = diag(1,1,-1,-1), central in Cl(2,1)", I4 == sp.diag(1, 1, -1, -1) and all(zM(I4*g - g*I4) for g in gens.values()))
# every element of Cl(2,1) is a real combination of the 8 basis monomials; all commute with I4
monomials = [sp.eye(4)] + list(gens.values()) + [G_C*G_H, G_C*G_G, G_H*G_G, I4]
check("all 8 basis monomials of Cl(2,1) commute with I_4: every element is BLOCK-DIAGONAL -- no matrix element connects the two sheets (superselection)",
      all(zM(m*I4 - I4*m) for m in monomials) and all(zM(m[0:2, 2:4]) and zM(m[2:4, 0:2]) for m in monomials))

print("=== a2-a3: the fourth generator ===")
al, be = sp.symbols('alpha beta', real=True)
G4 = off(al*g_C, be*g_C)
check("a2 Gamma_4 := [[0, alpha Gamma_C],[beta Gamma_C, 0]] ANTICOMMUTES with all three seat generators",
      all(zM(G4*g + g*G4) for g in gens.values()))
check("a2' Gamma_4 I_4 = -I_4 Gamma_4: it anticommutes with the pseudoscalar and is OFF-diagonal -- it maps block + to block - and back",
      zM(G4*I4 + I4*G4) and zM(G4[0:2, 0:2]) and zM(G4[2:4, 2:4]))
# uniqueness: any 4x4 anticommuting with all three must be of this form (solve generally)
X = sp.Matrix(4, 4, sp.symbols('x0:16', real=True))
eqs = []
for g in gens.values(): eqs += list(X*g + g*X)
sol = sp.solve(eqs, list(X), dict=True)
Xsol = X.subs(sol[0])
free = sorted(Xsol.free_symbols, key=str)
check("a3 the general solution of {Gamma_4, Gamma_X} = 0 for all three X is exactly the 2-parameter family above (off-diagonal, proportional to Gamma_C in each corner)",
      len(free) == 2 and zM(Xsol[0:2, 0:2]) and zM(Xsol[2:4, 2:4]) and zM(Xsol[0:2, 2:4] - (Xsol[0:2, 2:4][0, 1])*g_C) and zM(Xsol[2:4, 0:2] - (Xsol[2:4, 0:2][0, 1])*g_C),
      f"free parameters: {free}")
check("a3' Gamma_4^2 = -alpha beta I_4: sign of the fourth direction is alpha beta.  alpha beta = -1 -> Gamma_4^2 = +1, signature (+,+,+,-) = Cl(3,1);  alpha beta = +1 -> Gamma_4^2 = -1, (+,+,-,-) = Cl(2,2)",
      zM(G4*G4 + al*be*sp.eye(4)))
G4_space = G4.subs({al: 1, be: -1}); G4_time = G4.subs({al: 1, be: 1})
check("a3'' both realised: spacelike Gamma_4 (square +1) and timelike Gamma_4 (square -1) each anticommute with the three and exchange the blocks",
      G4_space**2 == sp.eye(4) and G4_time**2 == -sp.eye(4) and all(zM(G4_space*g + g*G4_space) and zM(G4_time*g + g*G4_time) for g in gens.values()))

print("=== a4: adjoining Gamma_4 lifts the superselection ===")
def centre_dim(generators):
    # dimension of the commutant of the generated algebra inside M_4(R), computed on a generic element
    Y = sp.Matrix(4, 4, sp.symbols('y0:16', real=True))
    eqs = []
    for g in generators: eqs += list(Y*g - g*Y)
    s = sp.solve(eqs, list(Y), dict=True)
    return len(Y.subs(s[0]).free_symbols)
check("a4 the commutant of Cl(2,1) in M_4(R) is 2-dimensional (spanned by 1 and I_4): two sectors",
      centre_dim(list(gens.values())) == 2)
check("a4' with the spacelike Gamma_4 adjoined the commutant is 1-dimensional (scalars only): the algebra is SIMPLE, the sectors are gone",
      centre_dim(list(gens.values()) + [G4_space]) == 1)
check("a4'' with the timelike Gamma_4 adjoined, likewise simple", centre_dim(list(gens.values()) + [G4_time]) == 1)

print("=== a5-a6: superposition made concrete; what Gamma_4 is (not an intertwiner) ===")
u1, u2, v1, v2 = sp.symbols('u1 u2 v1 v2', real=True)
psi_p = sp.Matrix([u1, u2, 0, 0]); psi_m = sp.Matrix([0, 0, v1, v2])
check("a5 every Cl(2,1) bilinear between the sheets vanishes: <psi_+| X |psi_-> = 0 for all 8 basis monomials -- no interference without Gamma_4",
      all(sp.simplify((psi_p.T*m*psi_m)[0]) == 0 for m in monomials))
check("a5' <psi_+| Gamma_4 |psi_-> = alpha (u1 v2 - u2 v1) != 0 generically: Gamma_4 creates the cross term -- the sheets can interfere",
      sp.simplify((psi_p.T*G4*psi_m)[0] - al*(u1*v2 - u2*v1)) == 0)
rho_p = {'C': g_C, 'H': g_H, 'G': g_G}; rho_m = {'C': -g_C, 'H': g_H, 'G': g_G}
# PREDICTED: Gamma_4's upper block intertwines the two irreps.  WRONG (recorded).  The two blocks are INEQUIVALENT
# Cl(2,1)-modules -- the pseudoscalar's value +-1 is a representation invariant -- so no intertwiner exists.
Aint = sp.Matrix(2, 2, sp.symbols('a0:4', real=True))
sol_int = sp.solve([e for k in 'CHG' for e in list(Aint*rho_m[k] - rho_p[k]*Aint)], list(Aint), dict=True)
check("a6 [my prediction was WRONG] there is NO 2x2 A with A rho_-(x) = rho_+(x) A for all three generators: the two blocks are INEQUIVALENT irreps (the pseudoscalar's sign is an invariant)",
      len(sol_int) == 1 and all(v == 0 for v in sol_int[0].values()) and len(Aint.subs(sol_int[0]).free_symbols) == 0)
check("a6' what Gamma_4 actually does: conjugation Gamma_4 x Gamma_4^{-1} = -x on every generator -- the GRADE INVOLUTION, an outer automorphism of Cl(2,1) that Gamma_4 makes inner in the larger algebra; it carries I to -I and so the + rep to the - rep as ABSTRACT reps, not by a change of basis",
      all(zM(G4_space*g*G4_space.inv() + g) for g in gens.values()) and zM(G4_space*I4*G4_space.inv() + I4))
n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT [DERIVED]: a fourth generator anticommuting with the seat's three exists, is unique up to two scalars,")
print("  is off-diagonal, anticommutes with the pseudoscalar, and exchanges the two sheets.  Without it Cl(2,1) has a")
print("  2-dim commutant (two superselection sectors: the hbar/G ordering is a label).  With it the algebra is simple")
print("  (M_4(R)) and the sheets can superpose and interfere.  The fourth direction's JOB is derived; its SIGNATURE")
print("  (spacelike -> Cl(3,1), 3+1; timelike -> Cl(2,2), 2+2) is open and is the kill condition for the 3+1 reading.")
