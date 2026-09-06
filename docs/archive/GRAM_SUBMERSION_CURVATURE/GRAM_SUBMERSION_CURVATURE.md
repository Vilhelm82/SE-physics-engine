# Curvature of the Gram submersion

## Setup

$$E=(S^2)^3 \ni A=[a_1\,a_2\,a_3],\qquad \pi:A\mapsto G_{ij}=a_i\!\cdot\!a_j,\qquad B=\pi(E)=\text{elliptope}$$

$$\dim E=6,\qquad \dim(SO(3)\text{-orbit})=3,\qquad \dim B=3$$

Product round metric = restriction of the Euclidean metric on $(\mathbb R^3)^3$.

| object | definition |
|---|---|
| vertical | $\omega\mapsto(\omega\times a_1,\ \omega\times a_2,\ \omega\times a_3)$ |
| momentum map | $J(v)=\sum_i a_i\times v_i$ |
| locked inertia tensor | $\mathbb I=3\,\mathrm{Id}-M,\qquad M=\sum_i a_ia_i^{\mathsf T}=AA^{\mathsf T}$ |
| mechanical connection | $\mathcal A(v)=\mathbb I^{-1}J(v)$ |
| horizontal | $v_i\cdot a_i=0$ and $J(v)=0$ |

$\mathbb I$ from $\sum_i\lVert\omega\times a_i\rVert^2=3\lVert\omega\rVert^2-\omega^{\mathsf T}M\omega$.

Horizontal lift of a base direction $dG$: solve the $9\times9$ system

$$v_i\cdot a_i=0,\qquad \sum_i a_i\times v_i=0,\qquad a_i\cdot v_j+a_j\cdot v_i=dG_{ij}$$

---

## Theorem (closed-form curvature)

$$\boxed{\;F(X,Y)=2\,\mathbb I^{-1}\!\Big(\sum_i X_i\times Y_i\Big),\qquad \mathbb I=3\,\mathrm{Id}-\textstyle\sum_i a_ia_i^{\mathsf T}\;}$$

for $X,Y$ horizontal.

**Proof.** $\mathcal A(X)=\mathcal A(Y)=0$ horizontally, so
$F(X,Y)=-\mathcal A([X,Y])=-\mathbb I^{-1}J([X,Y])$.
From $dJ(X,Y)=X[J(Y)]-Y[J(X)]-J([X,Y])$ and $J(X)=J(Y)=0$: $J([X,Y])=-dJ(X,Y)$.
With $J=\sum_i a_i\times da_i$, i.e. $J^k=\sum_i\varepsilon^{kmn}a_i^m\,da_i^n$:

$$dJ^k(X,Y)=\sum_i\varepsilon^{kmn}\big(X_i^mY_i^n-Y_i^mX_i^n\big)=2\sum_i(X_i\times Y_i)^k. \qquad\square$$

**O'Neill.** $K_{\text{base}}(X,Y)=K_{\text{total}}(X,Y)+3\lVert A_XY\rVert^2$, $A_XY=\tfrac12[X,Y]^{\text{vert}}$.
Total space $(\mathbb{CP}^1)^3$ Kähler homogeneous; excess curvature carried entirely by $F$.

**Verification** vs. integrated holonomy, $\lVert F\rVert_{\text{closed}}/(\theta/h^2)$:

| $\gamma$ | ratio |
|---|---|
| $(0.30,0.20,0.15)$ | 0.99981 |
| $(0.10,-0.25,0.40)$ | 0.99980 |
| $(-0.50,-0.50,-0.35)$ | 0.99807 |
| $(0.55,0.55,0.55)$ | 0.99950 |

---

## Exponent law [proved]

$$\Delta=\det G,\qquad V=\det A,\qquad V^2=\Delta$$

$$\boxed{\;\lVert F\rVert\ \sim\ C(b)\,\Delta^{-m/2}=C(b)\,|V|^{-m}\;}$$

$m$ = local vanishing order of $\Delta$; $C(b)$ depends on the boundary point $b$ (NOT universal).

### Proof of the exponents

Lift system $M(\varepsilon)X=b$ with $M_0=M(0)$ at the degenerate frame; divergent part of every lift lies in $\ker M_0$.

**$m=1$ (coplanar).** Exact: $\dim\ker M_0=1$ and the kernel is purely out-of-plane,
$N=(z_1e_z,\,z_2e_z,\,z_3e_z)$. Hence any two lifts have singular parts $X_{-1}=\alpha N$, $Y_{-1}=\beta N$, and

$$\sum_i X_{-1,i}\times Y_{-1,i}=\alpha\beta\sum_i z_i^2\,(e_z\times e_z)=0$$

identically — per site, parallel vectors. The $\varepsilon^{-2}$ term of $S=\sum_iX_i\times Y_i$ vanishes; leading term $\varepsilon^{-1}$ (cross terms in-plane $\times$ out-of-plane); $\mathbb I$ regular at coplanar frames. So $\lVert F\rVert=O(\varepsilon^{-1})=O(\Delta^{-1/2})$. $\square$

**$m=2$ (node).** Exact: $\dim\ker M_0=4$, kernel entirely $\perp e_3$. Divergent lift parts are $\perp e_3$, so their cross products are $\parallel e_3$ — no cancellation — and land exactly on the dying eigendirection of the inertia tensor:

$$\mathbb I_{33}=\frac{12t^2}{(1+t^2)^2}=3\sin^2\varepsilon\quad\text{(exact)}.$$

Numerator $\varepsilon^{-2}$ along $e_3$, times $\mathbb I^{-1}\sim(3\varepsilon^2)^{-1}$ on $e_3$: $\lVert F\rVert=O(\varepsilon^{-4})=O(\Delta^{-1})$ since $\Delta\sim\varepsilon^4$. $\square$

**Mechanism of the doubling:** kernel geometry. Parallel kernel (1-dim, shared direction) $\Rightarrow$ cancellation; perpendicular kernel (4-dim) $\Rightarrow$ no cancellation *and* alignment with the degenerating inertia direction.

### Exact-arithmetic pin (rational circle $c=\tfrac{1-t^2}{1+t^2},\ s=\tfrac{2t}{1+t^2}$; no floats)

slopes of $\log\lVert F\rVert^2$ vs $\log\Delta$:

| stratum | $t$-refinement | slope | target |
|---|---|---|---|
| coplanar | $1/20\to1/40\to1/80\to1/160$ | $-1.00496,\ -1.00125,\ -1.00031$ | $-1$ |
| node | $1/10\to1/20\to1/40\to1/80$ | $-2.02231,\ -2.00544,\ -2.00135$ | $-2$ |

### The constant is point-dependent

$C^2(b)=\lim\lVert F\rVert^2\Delta^{m}$ (Richardson-extrapolated from exact values):

| boundary point | $C$ |
|---|---|
| coplanar, in-plane dirs $(1,0),(\tfrac35,\tfrac45),(-\tfrac5{13},-\tfrac{12}{13})$ | $2.00156$ |
| coplanar, in-plane dirs $(1,0),(\tfrac45,\tfrac35),(-\tfrac{12}{13},\tfrac5{13})$ | $2.61798$ |
| coplanar, angles $(0°,100°,215°)$ | $1.0022$ |
| node, transverse dirs $(1,0),(-\tfrac35,\tfrac45),(-\tfrac5{13},-\tfrac{12}{13})$ | $0.99825$ |

$C$ varies by a factor $>2.6$ across boundary points: **any universal-constant claim is refuted.**

### Closed form of $C(b)$ on the smooth stratum [proved]

With $e_1=\gamma_{12}+\gamma_{13}+\gamma_{23}$, $e_2=\gamma_{12}\gamma_{13}+\gamma_{12}\gamma_{23}+\gamma_{13}\gamma_{23}$, $e_3=\gamma_{12}\gamma_{13}\gamma_{23}$ of the **boundary** Gram entries (on $\Delta=0$, so $e_3=\tfrac{e_1^2-2e_2-1}{2}$):

$$\boxed{\;C^2(b)\;=\;\frac{P(e_1,e_2)}{\big(\det\mathbb I_\parallel\big)^{4}},\qquad
\det\mathbb I_\parallel \;=\; \sum_{i<j}\sigma_{ij}^2 \;=\; 2(1-e_3)\;}$$

$$P=4e_1^6-4e_1^5-24e_1^4e_2-17e_1^4+16e_1^3e_2+52e_1^3+52e_1^2e_2^2+68e_1^2e_2+54e_1^2$$
$$\qquad-16e_1e_2^2-104e_1e_2-48e_1-40e_2^3-116e_2^2-108e_2-9$$

$\mathbb I_\parallel$ = in-plane block of the inertia tensor at the coplanar frame; $\det\mathbb I_\parallel=\det\!\big(\sum_i u_iu_i^{\mathsf T}\big)=\sum\sigma_{ij}^2$ by Cauchy–Binet ($\operatorname{tr}=3$). **The inertia tensor controls the constant as well as the exponent.**

**Derivation route** (leading-order cascade, all exact): $\ker M_0=\langle(\,\sin(\theta_2{-}\theta_3)e_z,\ \sin(\theta_3{-}\theta_1)e_z,\ \sin(\theta_1{-}\theta_2)e_z\,)\rangle$ [verified, ratios exactly $5/4$ at the test config]; solvability fixes the singular coefficients $c_p$; the in-plane reduced system ($x_i=w_ia_i^\perp+r_ia_i$, $dG_{ij}=(w_i-w_j)\sigma_{ij}+\dots$, $\sum w_i=0$) yields $S^{(-1)}_{pq}=c_pT(X_q^0)-c_qT(X_p^0)$; on-shell reduction $s_b^2\to1-\gamma_{12}^2$, $s_bs_g\to\gamma_{23}-\gamma_{12}\gamma_{13}$; symmetrize; reduce mod $\Delta$.

**Verification (all exact rationals).** Against the independent $9\times9$ recipe: config A $C^2=\tfrac{28981090229}{7233948160}$ ✓, config B $C^2=\tfrac{333432430282250}{48648964964161}$ ✓, plus a random boundary point ✓. $S_3$-invariance on the locus: exact, all six permutations. Positivity: 35/35 rational boundary samples $\ge0$.

**Distinguished value.** At the trine ($\gamma_{ij}=-\tfrac12$, the symmetric wall):

$$C^2=\frac{16}{27},\qquad C=\frac{4}{3\sqrt3}.$$

**Remaining open (scoped).** The node stratum: a node is a singular boundary point, and the limit of $\lVert F\rVert^2\Delta^2$ there may depend on approach direction; the tested direction gave $C\approx0.99825$. General-direction closed form via the same cascade with the 4-dim kernel — open.

---

## Inertia tensor connection

$$\mathbb I=3\,\mathrm{Id}-M,\qquad M=\sum_i a_ia_i^{\mathsf T}=AA^{\mathsf T}$$

$M$ and $G=A^{\mathsf T}A$ share nonzero spectrum: $\mathrm{rank}\,M=\mathrm{rank}\,G$, $\det M=\det G=\Delta$, $\mathrm{tr}\,M=3$.

$$\mathrm{eig}(\mathbb I)=3-\mathrm{eig}(M),\qquad \sum\mathrm{eig}(\mathbb I)=6,\qquad \det\mathbb I=\prod_k(3-\mu_k)$$

| stratum | $\mathrm{rank}\,G$ | $\mathrm{eig}(M)$ | $\mathrm{eig}(\mathbb I)$ | $\det\mathbb I$ |
|---|---|---|---|---|
| generic | 3 | $\mu_1,\mu_2,\mu_3>0$ | $3-\mu_k$ | $>0$ |
| coplanar ($\Delta=0$, rank 2) | 2 | $\mu_1,\mu_2,0$ | $3-\mu_1,3-\mu_2,3$ | $>0$ |
| node (rank 1) | 1 | $3,0,0$ | $0,3,3$ | $0$ |
| trine (coplanar, $\gamma_{ij}=-\tfrac12$) | 2 | $\tfrac32,\tfrac32,0$ | $\tfrac32,\tfrac32,3$ | $\tfrac{27}{4}$ |

$$\{\det\mathbb I=0\}=\{\text{nodes}\}\subsetneq\{\Delta=0\}$$

Degeneracy loci are nested, not equal: $\mathbb I$ is blind to the coplanar stratum and sees only the four rank-one points. $\Delta$ sees both.

$$m=1\iff\det\mathbb I\neq0\ \text{on }\{\Delta=0\};\qquad m=2\iff\det\mathbb I=0$$

so $m$ is readable off $\mathbb I$ alone.

**Same integer, two readouts:**

$$\text{parity of }m\ \longrightarrow\ \text{sheet flip of }V\text{ (odd flips, even does not)}$$
$$\text{value of }m\ \longrightarrow\ \text{curvature exponent }\lVert F\rVert\sim|V|^{-m}$$

**Open:** $\mathbb I$ does not currently appear in the seated-root construction. $M=AA^{\mathsf T}$ is the same matrix as the Gram data in a transposed role. Test whether $\{\det\mathbb I=0\}$ coincides with an independently named stratum.

---

## Functions

```python
import numpy as np
from scipy.integrate import solve_ivp

def inertia(A):
    return 3*np.eye(3) - A @ A.T                      # I = 3 Id - sum a_i a_i^T

def horiz(A, dG):
    """Horizontal lift. dG keyed by (0,1),(0,2),(1,2). Returns 3x3, columns v_i."""
    M = np.zeros((9, 9)); b = np.zeros(9); r = 0
    for i in range(3):                                 # v_i . a_i = 0
        M[r, 3*i:3*i+3] = A[:, i]; r += 1
    for k in range(3):                                 # sum_i a_i x v_i = 0
        e = np.eye(3)[k]
        for i in range(3):
            M[r, 3*i:3*i+3] = np.array(
                [np.dot(e, np.cross(A[:, i], np.eye(3)[m])) for m in range(3)])
        r += 1
    for (i, j) in [(0, 1), (0, 2), (1, 2)]:            # a_i.v_j + a_j.v_i = dG_ij
        M[r, 3*j:3*j+3] += A[:, i]
        M[r, 3*i:3*i+3] += A[:, j]
        b[r] = dG[(i, j)]; r += 1
    v, *_ = np.linalg.lstsq(M, b, rcond=None)
    return v.reshape(3, 3).T

def F_closed(A, pa, pb):
    """Curvature 2-form component for base coordinate directions pa, pb in {0,1,2}
       indexing (g12, g13, g23).  F = 2 I^{-1} sum_i (X_i x Y_i)."""
    keys = [(0, 1), (0, 2), (1, 2)]
    da = {k: 0.0 for k in keys}; db = dict(da)
    da[keys[pa]] = 1.0; db[keys[pb]] = 1.0
    X = horiz(A, da); Y = horiz(A, db)
    s = sum(np.cross(X[:, i], Y[:, i]) for i in range(3))
    return 2 * np.linalg.solve(inertia(A), s)

def F_norm(A):
    return np.sqrt(sum(np.linalg.norm(F_closed(A, p, q))**2
                       for p, q in [(0, 1), (0, 2), (1, 2)]))

def frame_from(g12, g13, g23):
    a1 = np.array([1., 0, 0]); s = np.sqrt(1 - g12**2)
    a2 = np.array([g12, s, 0])
    x = g13; y = (g23 - g12*g13)/s
    return np.column_stack([a1, a2, np.array([x, y, np.sqrt(max(1 - x*x - y*y, 0.0))])])

def frame_near_node(eps, u=None):
    """Frame approaching the node gamma=(1,1,1); Delta ~ eps^4."""
    if u is None:
        u = [np.array([1., 0.]), np.array([-.5, np.sqrt(3)/2]), np.array([-.5, -np.sqrt(3)/2])]
    return np.column_stack([np.array([eps*ui[0], eps*ui[1], np.sqrt(1 - eps*eps*ui@ui)])
                            for ui in u])

# ---- holonomy (independent check of F_closed) ----

def _rhs(t, y, d):
    A = y.reshape(3, 3)
    return horiz(A, {(0, 1): d[0], (0, 2): d[1], (1, 2): d[2]}).reshape(-1)

def transport(A, d, T=1.0):
    s = solve_ivp(_rhs, [0, T], A.reshape(-1), args=(d,), rtol=1e-11, atol=1e-13)
    A = s.y[:, -1].reshape(3, 3)
    return A / np.linalg.norm(A, axis=0)

def hol_angle(A0, h, plane):
    """Rotation angle after a square loop of side h in the given coordinate plane."""
    A = A0.copy(); i, j = plane
    ds = [np.zeros(3) for _ in range(4)]
    ds[0][i] = h; ds[1][j] = h; ds[2][i] = -h; ds[3][j] = -h
    for d in ds:
        A = transport(A, d)
    R = A @ np.linalg.inv(A0)
    return np.arccos(np.clip((np.trace(R) - 1)/2, -1, 1))

def F_holonomy(A0, hb):
    """angle/h^2, Richardson-extrapolated h->0, combined over the three planes.
       Use hb <= 0.05*Delta near the branch locus."""
    comps = []
    for plane in [(0, 1), (0, 2), (1, 2)]:
        v = [hol_angle(A0, h, plane)/h**2 for h in (hb, hb/2, hb/4)]
        r1 = [2*v[k+1] - v[k] for k in range(2)]
        comps.append(2*r1[1] - r1[0])
    return np.sqrt(sum(c*c for c in comps))
```

### Exact-arithmetic functions (sympy, rationals only)

```python
import sympy as sp
from sympy import Rational as R

def cross_s(u,v):
    return sp.Matrix([u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]])

def build_M(a):
    """exact 9x9 horizontal-lift system for a frame a = [a1,a2,a3]"""
    M=sp.zeros(9,9); r=0
    for i in range(3):
        for c in range(3): M[r,3*i+c]=a[i][c]
        r+=1
    for k in range(3):
        for i in range(3):
            for c in range(3):
                M[r,3*i+c]=cross_s(a[i],sp.Matrix([1 if j==c else 0 for j in range(3)]))[k]
        r+=1
    for (i,j) in [(0,1),(0,2),(1,2)]:
        for c in range(3):
            M[r,3*j+c]+=a[i][c]; M[r,3*i+c]+=a[j][c]
        r+=1
    return M

def lift_exact(a,dG):
    M=build_M(a); b=sp.zeros(9,1)
    for idx,_ in enumerate([(0,1),(0,2),(1,2)]): b[6+idx]=dG[idx]
    X=M.solve(b)
    return [X[3*i:3*i+3,0] for i in range(3)]

def inertia_s(a):
    return 3*sp.eye(3)-sum((v*v.T for v in a),sp.zeros(3,3))

def F2_exact(a):
    """||F||^2, exact rational, summed over the three coordinate planes"""
    I=inertia_s(a)
    H=[lift_exact(a,[1,0,0]),lift_exact(a,[0,1,0]),lift_exact(a,[0,0,1])]
    tot=0
    for p,q in [(0,1),(0,2),(1,2)]:
        S=sum((cross_s(H[p][i],H[q][i]) for i in range(3)),sp.zeros(3,1))
        Fpq=2*I.solve(S)
        tot+=(Fpq.T*Fpq)[0,0]
    return sp.simplify(tot)

def cs(t):
    """rational point on the circle: exact unit vectors, no floats"""
    return (1-t**2)/(1+t**2), 2*t/(1+t**2)
```

Numerics: `solve_ivp` rtol `1e-11`, atol `1e-13`; loop side $h\le0.05\,\Delta$ near the branch locus (larger $h$ lets the loop sample a varying $\Delta$ and biases the exponent).
