# THM-TARGET N — Kerr at second order in spin

**Outcome:** the minimal continuation of THM-M is **killed before the Kerr quadrupole gate**.

Artifact: `thm_n_kerr_quadrupole.py`. Regression tests:
`test_thm_n_kerr_quadrupole.py` (4/4). The theorem script carries eight exact
checks by two independent angular-average routes.

## Question

Does the existing THM-M construction — boost each element's static presented
rank-two field and superpose the elements linearly — continue from the correct
first-order Lense–Thirring field to Kerr's second-order mass quadrupole

$$Q_{\rm Kerr}=-\frac{J^2}{Mc^2}?$$

This is the $\ell=2$ member of the Kerr field-multipole relation
$M_\ell+iS_\ell=M(ia)^\ell$ (comparison-stage lineage: Hansen, 1974,
doi:10.1063/1.1666501).

## Admissible inputs

- THM-K/H-5, to first order in $r_s$: each static element carries
  $h_{00}=\rho$ and $h_{ij}=\rho n_i n_j$, with
  $\rho=2G\,dm/(c^2s)$.
- BARE-1: the presented field transforms as a rank-two Lorentz tensor.
- K-6: linear superposition over the source.
- The same uniformly rotating thin spherical shell used by THM-M, with
  $J=(2/3)MR^2\Omega$.

Not admitted: the Kerr metric or quadrupole, Einstein's equations, a source
stress/tension law, retardation, nonlinear gravitational self-field terms, or
an after-the-result cancellation term.

## The exact expansion

Let $\epsilon=\Omega R/c$, $\delta=R/r$, $N_z=\cos\Theta$, and
$P_2=(3\cos^2\Theta-1)/2$. For a source direction $e$, put
$a=N\mathbin{\cdot}e$ and $b=\hat z\times e$. With
$S=r/|rN-Re|$, the literal continuation of THM-M gives

$$
\frac{h'_{00}}{2GM/(c^2r)}
=\langle S\rangle
+\epsilon^2\!\left[\langle b^2S\rangle
+\langle(b\mathbin{\cdot}N)^2S^3\rangle\right]
+O(\epsilon^4).
$$

Exact integration over the shell returns

$$
\frac{h'_{00}}{2GM/(c^2r)}
=1+\epsilon^2\left[
\frac89-\frac29P_2
-\frac{2}{15}\delta^2P_2+O(\delta^4)
\right]+O(\epsilon^4).
$$

There is an ambiguity invisible at first order: full Lorentz covariance can
also be read as transforming the field argument into each element's
instantaneous rest frame. The exact boost gives
$s'/s=1+(\beta\cdot n)^2/2+O(\beta^4)$, halving the rod term. Repeating the
integration changes the $1/r$ bracket to

$$\frac79-\frac19P_2,$$

while leaving the formal $-(2/15)\delta^2P_2$ coefficient unchanged. Thus both
admissible minimal readings fail the same gate.

## Verdict

1. **Asymptotic gate — killed.** A stationary Kerr exterior has a scalar mass
   monopole at order $1/r$ in asymptotically Cartesian mass-centred
   asymptotics. The construction instead produces a direction-dependent
   $P_2/r$ term: coefficient $-2/9$ under the literal THM-M continuation and
   $-1/9$ under the full-argument continuation. An isotropic mass
   renormalisation cannot absorb it.

2. **Formal quadrupole gate — also killed.** Ignoring the already-fatal term,
   the $1/r^3$ coefficient would be

   $$Q_{\rm model}=-\frac{2}{15}\frac{M\Omega^2R^4}{c^2}
   =\frac{3}{10}Q_{\rm Kerr}.$$

The result does **not** kill the seated-root model. It kills the proposition
that THM-M's first-order boost plus linear superposition is already a
second-order rotating field law. The first-order coefficient was genuinely
recovered; at second order, the missing source-stress/retardation/nonlinear
closure becomes visible. Which of those words belongs to the model is now a
derivation debt, not permission to import a cancellation.

## Consequence for the paper

The current scope statement — first order in $r_s$ and $J$, full rotating
solution not reached — survives exactly. Any stronger phrase such as “the Kerr
sector is derived” does not. The next valid Kerr attempt must first derive a
nonlinear source/field closure which predicts both cancellation of the
$P_2/r$ term and the remaining $Q=-J^2/(Mc^2)$ coefficient before comparison.
