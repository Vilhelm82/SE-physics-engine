#!/usr/bin/env python3
"""Cella full-tensor curvature numerators and certified weighted finite jets.

Run ``python -B cella_tensor_valuation.py`` for the exact replay.  Import
CurvatureNumerators for arbitrary symbolic coordinate metrics, or MonomialMetric
for finite sums of monomials with smooth/analytic units.  No intermediate
Laurent window is imposed: jet depth follows from the requested output weight,
and all polynomial numerator terms are assembled before output filtering.

Standard tools retained explicitly: coordinate Levi-Civita/Riemann definitions,
matrix adjugates, finite Taylor jets and exact SymPy arithmetic.  The manuscript
docs/2026-09-05-cella-tensor-valuation.md supplies the proofs and scope.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import cached_property, lru_cache
from itertools import product
from pathlib import Path
import json
import random
import sys
import time

import sympy as s


def clean(expr):
    reduced = s.cancel(s.sympify(expr))
    # Arbitrary transverse jets carry many independent Derivative atoms. Full
    # factorization over those generators is unnecessary for exact cancellation.
    return s.factor_terms(reduced) if reduced.has(s.Derivative) else s.factor(reduced)


class CurvatureNumerators:
    """R_ijkl = N_ijkl / det(g), with sphere-positive curvature convention.

    ``first[d][i,j]`` and ``second[d][e][i,j]`` optionally supply an exact
    coordinate two-jet at a point.  Otherwise derivatives are taken symbolically.
    Components always refer to the supplied coordinate basis.
    """

    def __init__(self, metric, coordinates, first=None, second=None):
        self.g = s.ImmutableMatrix(metric)
        self.coordinates = tuple(coordinates)
        self.n = len(self.coordinates)
        if self.g.shape != (self.n, self.n) or self.g != self.g.T:
            raise ValueError("A symmetric coordinate metric is required")
        self.first, self.second = first, second

    @cached_property
    def delta(self):
        return clean(self.g.det())

    @cached_property
    def adjugate(self):
        return self.g.adjugate().applyfunc(clean)

    @lru_cache(None)
    def derivative(self, i, j, *directions):
        if len(directions) == 1 and self.first is not None:
            return self.first[directions[0]][i, j]
        if len(directions) == 2 and self.second is not None:
            return self.second[directions[0]][directions[1]][i, j]
        return s.diff(self.g[i, j], *(self.coordinates[d] for d in directions))

    @lru_cache(None)
    def first_kind(self, a, b, c):
        return (self.derivative(a, b, c) + self.derivative(a, c, b)
                - self.derivative(b, c, a)) / 2

    @lru_cache(None)
    def riemann_numerator(self, i, j, k, l):
        if i == j or k == l:
            return s.S.Zero
        if i > j:
            return -self.riemann_numerator(j, i, k, l)
        if k > l:
            return -self.riemann_numerator(i, j, l, k)
        if (i, j) > (k, l):
            return self.riemann_numerator(k, l, i, j)
        linear = (self.derivative(i, l, j, k)
                  + self.derivative(j, k, i, l)
                  - self.derivative(i, k, j, l)
                  - self.derivative(j, l, i, k)) / 2
        quadratic = sum(
            self.adjugate[a, b]
            * (self.first_kind(a, i, l) * self.first_kind(b, j, k)
               - self.first_kind(a, i, k) * self.first_kind(b, j, l))
            for a in range(self.n) for b in range(self.n))
        return clean(self.delta * linear + quadratic)

    def riemann(self, i, j, k, l):
        return clean(self.riemann_numerator(i, j, k, l) / self.delta)

    @lru_cache(None)
    def covariant_riemann_numerator(self, indices, directions=()):
        """Numerator of iterated covariant derivatives, denominator delta^(1+r).

        ``directions[-1]`` is the last derivative applied.  Earlier derivative
        indices are themselves covariant tensor slots in the next derivative.
        A supplied numerical two-jet cannot determine these higher derivatives.
        """
        indices, directions = tuple(indices), tuple(directions)
        if len(indices) != 4:
            raise ValueError("Four lower curvature indices required")
        if not directions:
            return self.riemann_numerator(*indices)
        if self.first is not None:
            raise ValueError("Covariant curvature derivatives require higher metric jets")
        d = directions[-1]
        previous_directions = directions[:-1]
        previous = self.covariant_riemann_numerator(indices, previous_directions)
        previous_power = len(directions)
        out = (self.delta*s.diff(previous,self.coordinates[d])
               - previous_power*s.diff(self.delta,self.coordinates[d])*previous)
        slots = indices+previous_directions
        for position,index in enumerate(slots):
            for a in range(self.n):
                replaced = slots[:position]+(a,)+slots[position+1:]
                other = self.covariant_riemann_numerator(replaced[:4],replaced[4:])
                if other != 0:
                    out -= sum(self.adjugate[a,b]*self.first_kind(b,d,index)
                               for b in range(self.n))*other
        return clean(out)

    @lru_cache(None)
    def ricci_numerator(self, j, l):
        return clean(sum(self.adjugate[i, k]
                         * self.riemann_numerator(i, j, k, l)
                         for i in range(self.n) for k in range(self.n)))

    @cached_property
    def ricci(self):
        return s.Matrix(self.n, self.n, lambda j, l:
                        clean(self.ricci_numerator(j, l) / self.delta**2))

    @cached_property
    def scalar_numerator(self):
        return clean(sum(self.adjugate[j, l] * self.ricci_numerator(j, l)
                         for j in range(self.n) for l in range(self.n)))

    @cached_property
    def scalar(self):
        return clean(self.scalar_numerator / self.delta**3)

    @cached_property
    def nonzero_numerators(self):
        return {indices: value for indices in product(range(self.n), repeat=4)
                if (value := self.riemann_numerator(*indices)) != 0}

    def contraction_numerator(self, curvature_count, pairing, derivative_orders=None):
        """Complete contraction of k curvature tensors with optional derivatives.

        ``pairing`` pairs all 4*k tensor slots once.  For Riemann squared use
        k=2 and ((0,4),(1,5),(2,6),(3,7)); no signature positivity is assumed.
        With r total covariant derivatives, there are 4*k+r slots, and the
        denominator is det(g)^(3*k+3*r/2).  Every slot must be paired exactly once.
        """
        k = curvature_count
        derivative_orders = tuple(derivative_orders or (0,)*k)
        if len(derivative_orders) != k or any(r < 0 for r in derivative_orders):
            raise ValueError("One nonnegative derivative count per curvature factor")
        total_slots = 4*k+sum(derivative_orders)
        if k < 1 or sorted(a for pair in pairing for a in pair) != list(range(total_slots)):
            raise ValueError("Pairing must cover every curvature slot once")
        entries = []
        for order in derivative_orders:
            if order == 0:
                entries.append(tuple(self.nonzero_numerators.items()))
            else:
                entries.append(tuple((slots,value)
                    for slots in product(range(self.n),repeat=4+order)
                    if (value := self.covariant_riemann_numerator(slots[:4],slots[4:])) != 0))
        answer = s.S.Zero
        for factors in product(*entries):
            indices = tuple(index for entry, _ in factors for index in entry)
            coefficient = s.prod(value for _, value in factors)
            for a, b in pairing:
                coefficient *= self.adjugate[indices[a], indices[b]]
                if coefficient == 0:
                    break
            answer += coefficient
        return clean(answer)

    @cached_property
    def kretschmann(self):
        numerator = self.contraction_numerator(2, ((0, 4), (1, 5), (2, 6), (3, 7)))
        return clean(numerator / self.delta**6)

    @cached_property
    def scalar_gradient_squared(self):
        if self.first is not None:
            raise ValueError("Scalar curvature derivatives require higher metric jets")
        gradient_numerators = [self.delta*s.diff(self.scalar_numerator,coordinate)
            -3*self.scalar_numerator*s.diff(self.delta,coordinate)
            for coordinate in self.coordinates]
        return clean(sum(self.adjugate[a,b]*gradient_numerators[a]*gradient_numerators[b]
                         for a in range(self.n) for b in range(self.n))/self.delta**9)


def required_unit_degree(cutoff, degree, derivative_order, component_floor, weights):
    """Conservative proved depth; depends on the output weight, never a fixed window.

    A numerator with d metric factors and r total derivatives has lower bound
    b=d*mu-r*max(w).  Degree N Taylor units suffice when (N+1)*min(w)>cutoff-b;
    r extra degrees are retained to control differentiated Taylor remainders.
    """
    weights = tuple(map(s.Rational, weights))
    if not weights or min(weights) <= 0:
        raise ValueError("Positive rational weights are required by this replay")
    bound = degree * component_floor - derivative_order * max(weights)
    return max(0, int(s.floor((cutoff - bound) / min(weights)))) + derivative_order


def total_taylor(expr, normals, degree):
    """Taylor polynomial in normals; arbitrary transverse functions stay symbolic."""
    scale = s.Dummy("jet_scale")
    shifted = s.sympify(expr).subs({z: scale*z for z in normals}, simultaneous=True)
    return s.expand(s.series(shifted, scale, 0, degree+1).removeO().subs(scale, 1))


def weighted_coefficients(expr, normals, weights, ratios, cutoff=None, transverse=None):
    """Group an exact finite Laurent sum AFTER all collisions and ratio substitutions.

    Real exponents are allowed by the theorem.  The executable API uses rational
    exponents/weights and requires each coefficient to be independent of normals.
    A zero result means no terms in this supplied polynomial, not smooth flatness.
    """
    normals, weights, ratios = tuple(normals), tuple(map(s.Rational, weights)), tuple(ratios)
    if not (len(normals) == len(weights) == len(ratios)) or min(weights) <= 0:
        raise ValueError("Matching normals, positive weights and approach ratios required")
    out = {}
    expr = s.expand(s.sympify(expr).subs(transverse or {}))
    for term in s.Add.make_args(expr):
        powers = term.as_powers_dict()
        alpha = tuple(s.Rational(powers.get(z, 0)) for z in normals)
        coefficient = clean(term / s.prod(z**a for z, a in zip(normals, alpha)))
        if coefficient.has(*normals):
            raise ValueError("Finite Laurent polynomial required; Taylor-expand units first")
        weight = sum(w*a for w, a in zip(weights, alpha))
        if cutoff is None or weight <= cutoff:
            out[weight] = out.get(weight, 0) + coefficient*s.prod(
                s.sympify(rho)**a for rho, a in zip(ratios, alpha))
    return {weight: value for weight in sorted(out)
            if (value := clean(out[weight])) != 0}


def rational_leading(expr, normals, weights, ratios, transverse=None):
    """Exact leading order of a supplied rational Laurent expression along a path."""
    numerator, denominator = s.cancel(expr).as_numer_denom()
    nc = weighted_coefficients(numerator, normals, weights, ratios, transverse=transverse)
    dc = weighted_coefficients(denominator, normals, weights, ratios, transverse=transverse)
    if not dc:
        raise ZeroDivisionError("This path lies in the supplied denominator's zero set")
    if not nc:
        return {"status": "zero on this exact path", "order": None, "coefficient": s.S.Zero}
    nu, delta = min(nc), min(dc)
    return {"status": "finite leading coefficient", "order": nu-delta,
            "coefficient": clean(nc[nu]/dc[delta])}


@dataclass
class MonomialMetric:
    """g_ij = sum z^alpha h_ij,alpha, stored once for each i<=j.

    Coordinates start with the normal variables.  ``terms[(i,j)]`` contains
    ``(alpha, unit_expression)`` pairs, where alpha has normal_count entries.
    Missing components are exactly zero.  Units may include transverse functions.
    """
    coordinates: tuple
    normal_count: int
    terms: dict

    @property
    def normals(self):
        return self.coordinates[:self.normal_count]

    def metric(self, unit_degree=None):
        metric = s.zeros(len(self.coordinates))
        for (i, j), entries in self.terms.items():
            if i > j:
                raise ValueError("Store symmetric components once, with i<=j")
            value = s.S.Zero
            for alpha, unit in entries:
                if len(alpha) != self.normal_count:
                    raise ValueError("Exponent dimension differs from normal_count")
                h = unit if unit_degree is None else total_taylor(unit, self.normals, unit_degree)
                value += s.prod(z**s.Rational(a) for z, a in zip(self.normals, alpha))*h
            metric[i, j] = metric[j, i] = value
        return metric

    def numerator_jet(self, weights, cutoff, ratios, kind="scalar", component=None, transverse=None):
        """Certified coefficients through cutoff; no claim of zero beyond cutoff.

        For scalar, Ricci and lower Riemann the denominator powers are 3,2,1.
        Exact zeros beyond the requested cutoff are deliberately not inferred.
        """
        n = len(self.coordinates)
        mu = min(sum(s.Rational(w)*s.Rational(a) for w, a in zip(weights, alpha))
                 for entries in self.terms.values() for alpha, _ in entries)
        degree, derivatives, power = {
            "determinant": (n, 0, 0), "riemann": (n+1, 2, 1),
            "ricci": (2*n, 2, 2), "scalar": (3*n-1, 2, 3),
            "covariant_riemann": ((len(component or ())-3)*n+1,len(component or ())-2,len(component or ())-3)
        }[kind]
        if kind == "covariant_riemann" and len(component or ()) < 4:
            raise ValueError("Curvature component followed by its derivative indices required")
        depth = required_unit_degree(s.Rational(cutoff), degree, derivatives, mu, weights)
        tensor = CurvatureNumerators(self.metric(depth), self.coordinates)
        if kind == "determinant":
            numerator = tensor.delta
        elif kind == "riemann":
            numerator = tensor.riemann_numerator(*component)
        elif kind == "ricci":
            numerator = tensor.ricci_numerator(*component)
        elif kind == "covariant_riemann":
            numerator = tensor.covariant_riemann_numerator(tuple(component[:4]),tuple(component[4:]))
        else:
            numerator = tensor.scalar_numerator
        coefficients = weighted_coefficients(numerator, self.normals, weights, ratios,
                                             cutoff=cutoff, transverse=transverse)
        return {"coefficients": coefficients, "unit_taylor_degree": depth,
                "verified_through_weight": s.Rational(cutoff), "denominator_power": power,
                "status": "nonzero coefficient found" if coefficients else "no coefficient through cutoff"}

    def leading_through(self, weights, ratios, numerator_cutoff, determinant_cutoff,
                        kind="scalar", component=None, transverse=None):
        """Find a certified leading term if both requested finite searches reach it.

        An empty jet is returned as undecided beyond its stated bound.  Increasing
        the cutoffs eventually succeeds exactly under the theorem's finite-first-
        nonzero hypotheses.  Smooth flatness or an identically zero tensor needs
        separate data; this method never silently declares either from a zero jet.
        """
        numerator = self.numerator_jet(weights, numerator_cutoff, ratios, kind,
                                       component, transverse)
        determinant = self.numerator_jet(weights, determinant_cutoff, ratios,
                                         "determinant", transverse=transverse)
        if not numerator["coefficients"] or not determinant["coefficients"]:
            return {"status":"no leading term certified within these bounds",
                    "numerator":numerator,"determinant":determinant}
        nu, delta = min(numerator["coefficients"]), min(determinant["coefficients"])
        power = numerator["denominator_power"]
        return {"status":"finite leading coefficient", "order":nu-power*delta,
                "coefficient":clean(numerator["coefficients"][nu]
                                     / determinant["coefficients"][delta]**power),
                "numerator":numerator,"determinant":determinant}


def direct_two_jet_riemann(metric, first, second):
    """Independent rational Christoffel + differentiated-inverse oracle at a point."""
    n = metric.rows
    inverse = [[Fraction(value) for value in row] for row in metric.inv().tolist()]
    g = [[Fraction(value) for value in row] for row in metric.tolist()]
    dg = [[[Fraction(first[d][i, j]) for j in range(n)] for i in range(n)] for d in range(n)]
    ddg = [[[[Fraction(second[d][e][i, j]) for j in range(n)] for i in range(n)]
            for e in range(n)] for d in range(n)]
    dg_inverse = [[[-sum(inverse[a][p]*dg[d][p][q]*inverse[q][b]
                         for p in range(n) for q in range(n))
                    for b in range(n)] for a in range(n)] for d in range(n)]
    gamma = [[[sum(inverse[a][p]*(dg[b][p][c]+dg[c][p][b]-dg[p][b][c])/2
                   for p in range(n)) for c in range(n)] for b in range(n)] for a in range(n)]
    dgamma = [[[[sum(
        dg_inverse[d][a][p]*(dg[b][p][c]+dg[c][p][b]-dg[p][b][c])/2
        + inverse[a][p]*(ddg[d][b][p][c]+ddg[d][c][p][b]-ddg[d][p][b][c])/2
        for p in range(n)) for c in range(n)] for b in range(n)]
        for a in range(n)] for d in range(n)]
    return { (i,j,k,l): sum(g[i][a]*(
        dgamma[k][a][l][j]-dgamma[l][a][k][j]
        +sum(gamma[a][k][b]*gamma[b][l][j]-gamma[a][l][b]*gamma[b][k][j]
             for b in range(n))) for a in range(n))
        for i,j,k,l in product(range(n), repeat=4) }


def main():
    start = time.monotonic()
    checks = {}
    results = {}

    def check(name, actual, expected=0):
        residual = clean(actual-expected)
        if residual != 0:
            raise AssertionError((name, residual))
        checks[name] = {"verified": True}
        if "--progress" in sys.argv:
            print(f"PASS {name}",file=sys.stderr,flush=True)

    # All inertia possibilities in dimensions 2,3,4, with non-diagonal metrics,
    # arbitrary rational first derivatives, and symmetric second derivatives.
    rng = random.Random(20260905)
    oracle_components = 0
    for n in (2, 3, 4):
        coordinates = s.symbols(f"q0:{n}")
        for negatives in range(n+1):
            for sample in range(2):
                change = s.eye(n)
                for i in range(n):
                    for j in range(i+1,n):
                        change[i,j] = rng.choice((-2,-1,1,2))
                metric = change.T*s.diag(*([-1]*negatives+[1]*(n-negatives)))*change
                def symmetric():
                    out = s.zeros(n)
                    for i in range(n):
                        for j in range(i,n):
                            out[i,j] = out[j,i] = rng.randint(-3,3)
                    return out
                first = [symmetric() for _ in range(n)]
                second = [[None]*n for _ in range(n)]
                for d in range(n):
                    for e in range(d,n):
                        second[d][e] = second[e][d] = symmetric()
                tensor = CurvatureNumerators(metric, coordinates, first, second)
                oracle = direct_two_jet_riemann(metric, first, second)
                for indices, expected in oracle.items():
                    if tensor.riemann(*indices) != expected:
                        raise AssertionError((n, negatives, sample, indices))
                    i,j,k,l = indices
                    if clean(tensor.riemann(i,j,k,l)+tensor.riemann(i,k,l,j)
                             +tensor.riemann(i,l,j,k)) != 0:
                        raise AssertionError("algebraic Bianchi identity")
                    oracle_components += 1
                checks[f"full_tensor_n{n}_negative{negatives}_sample{sample}"] = {"verified": True}
    results["independent_rational_two_jet_components"] = oracle_components

    # Calibration and signed Cella laws.
    x, y, z = s.symbols("x y z", positive=True)
    sphere = CurvatureNumerators(s.diag(1,s.sin(x)**2),(x,y))
    check("sphere_positive_convention", s.trigsimp(sphere.scalar), 2)
    for sign in (1,-1):
        generic = CurvatureNumerators(s.diag(sign*2*x**2,3+5*x,7-2*x),(x,y,z))
        lead = rational_leading(generic.scalar,(x,),(1,),(1,))
        check(f"generic_signed_order_{sign}",lead["order"],-3)
        check(f"generic_signed_coefficient_{sign}",lead["coefficient"],s.Rational(29,42)/sign)
        parity = CurvatureNumerators(s.diag(sign*3*x**2,-2*x**-2,5*x**-2),(x,y,z))
        check(f"parity_signed_law_{sign}",parity.scalar,-s.Rational(14,3)/sign*x**-4)

    # A smooth nontrivial pullback of the diagonal parity metric.  The new chart
    # has off-diagonal components and changes the normal defining function.
    u, v, w = s.symbols("u v w", positive=True)
    old = (x,y,z)
    new = (u,v,w)
    mapping = s.Matrix([u*(1+v),v+u,w+u*v])
    jac = mapping.jacobian(new)
    source = s.diag(3*x**2,2*x**-2,-5*x**-2)
    original = CurvatureNumerators(source,old)
    pulled_metric = jac.T*source.subs(dict(zip(old,mapping)),simultaneous=True)*jac
    pulled = CurvatureNumerators(pulled_metric,new)
    check("non_diagonal_pullback_scalar",pulled.scalar,
          original.scalar.subs(dict(zip(old,mapping)),simultaneous=True))
    # Transform EVERY component, then compare the independently recomputed tensor.
    old_components = {indices:clean((numerator/original.delta).subs(
        dict(zip(old,mapping)),simultaneous=True))
        for indices,numerator in original.nonzero_numerators.items()}
    for i,j,k,l in product(range(3),repeat=4):
        transformed = sum(value*jac[a,i]*jac[b,j]*jac[c,k]*jac[d,l]
                          for (a,b,c,d),value in old_components.items())
        if clean(pulled.riemann(i,j,k,l)-transformed) != 0:
            raise AssertionError(("full pullback tensor",i,j,k,l))
    checks["non_diagonal_pullback_all_81_components"] = {"verified":True}
    check("pullback_divisor_order",rational_leading(pulled.scalar,(u,),(1,),(1,),{v:2,w:1})["order"],-4)
    check("pullback_divisor_coefficient",rational_leading(pulled.scalar,(u,),(1,),(1,),{v:2,w:1})["coefficient"],-s.Rational(14,243))

    # Off-diagonal components can change a leading law when they are new metric
    # data, rather than the pullback of the complete original metric.
    c = s.symbols("c",real=True)
    changed = CurvatureNumerators(s.Matrix([[3*x**2,c],[c,2*x**-2]]),(x,y))
    check("off_diagonal_changes_parity_residue",changed.scalar,-12/(6-c**2)*x**-4)
    results["independent_off_diagonal_residue"] = str(changed.scalar)
    changed_generic = CurvatureNumerators(s.Matrix([[3*x**2,c*x],[c*x,2+5*x]]),(x,y))
    changed_lead = rational_leading(changed_generic.scalar,(x,),(1,),(1,))
    check("off_diagonal_generic_order",changed_lead["order"],-3)
    check("off_diagonal_changes_generic_residue",changed_lead["coefficient"],5/(6-c**2))

    # Special ratio cancellation, then the next inward layer.  A scalar-flat
    # result cannot be inferred by deleting the two outer monomials.
    corner = CurvatureNumerators(s.diag(x**-2*y**-2,x**2,-y**2),(z,x,y))
    check("signed_corner_exact",corner.scalar,-6*x**-4+6*y**-4)
    balanced = rational_leading(corner.scalar,(x,y),(1,1),(1,1))
    check("balanced_front_cancels",balanced["coefficient"],0)
    changed_path = clean(corner.scalar.subs({x:u,y:u*(1+u)},simultaneous=True))
    next_layer = rational_leading(changed_path,(u,),(1,),(1,))
    check("curved_path_next_layer_order",next_layer["order"],-3)
    check("curved_path_next_layer_coefficient",next_layer["coefficient"],-24)

    # The determinant can also lose its first face; classify it before division.
    det_cancel = CurvatureNumerators(s.Matrix([[1,1],[1,1+x-y+x*x]]),(x,y))
    check("determinant_multimonomial",det_cancel.delta,x-y+x*x)
    dlead = rational_leading(det_cancel.delta,(x,y),(1,1),(1,1))
    check("determinant_face_cancellation_order",dlead["order"],2)
    actual = rational_leading(det_cancel.scalar,(x,y),(1,1),(1,1))
    check("determinant_face_cancellation_curvature_order",actual["order"],-4)
    check("determinant_face_cancellation_curvature_coefficient",actual["coefficient"],s.Rational(1,2))
    results["determinant_cancellation_curvature"] = {k:str(value) for k,value in actual.items()}

    # No arbitrary unit-jet window: a late coefficient requires a later request.
    late = 9
    h = (1+x**late)**2
    smooth_germ = MonomialMetric((x,y),1,{
        (0,0):[((0,),1+h)],(0,1):[((0,),h)],(1,1):[((0,),h)]})
    early = smooth_germ.numerator_jet((1,),4,(1,))
    late_jet = smooth_germ.numerator_jet((1,),late-2,(1,))
    if early["coefficients"]:
        raise AssertionError("A zero initial jet was misclassified")
    checks["finite_jet_early_answer_is_only_bounded"] = {"verified":True}
    check("finite_jet_late_numerator_order",min(late_jet["coefficients"]),late-2)
    check("finite_jet_late_numerator_coefficient",late_jet["coefficients"][late-2],-2*late*(late-1))
    exact_late = CurvatureNumerators(smooth_germ.metric(),(x,y))
    check("finite_jet_matches_exact_curvature",exact_late.scalar,
          -2*late*(late-1)*x**(late-2)/(1+x**late))
    late_answer = smooth_germ.leading_through((1,),(1,),late-2,0)
    check("finite_jet_quotient_order",late_answer["order"],late-2)
    check("finite_jet_quotient_coefficient",late_answer["coefficient"],-2*late*(late-1))
    results["automatic_jet_depths"] = {
        "early_cutoff":4,"early_unit_degree":early["unit_taylor_degree"],
        "later_cutoff":late-2,"later_unit_degree":late_jet["unit_taylor_degree"]}

    h3 = (1+x**3)**2
    differential = CurvatureNumerators(s.Matrix([[1+h3,h3],[h3,h3]]),(x,y))
    cov_x = differential.covariant_riemann_numerator((0,1,0,1),(0,))/differential.delta**2
    cov_y = differential.covariant_riemann_numerator((0,1,0,1),(1,))/differential.delta**2
    check("covariant_derivative_metric_compatibility_x",cov_x,
          s.diff(differential.scalar,x)*differential.delta/2)
    check("covariant_derivative_metric_compatibility_y",cov_y)
    cov_yy = differential.covariant_riemann_numerator((0,1,0,1),(1,1))/differential.delta**3
    check("second_covariant_derivative_includes_derivative_slot",cov_yy,
          -18*x*x*(1+x**3)*(1-2*x**3))
    check("differential_scalar_invariant",differential.scalar_gradient_squared,
          144*(1-2*x**3)**2/(1+x**3)**4)
    contracted = differential.contraction_numerator(2,
        ((0,2),(1,3),(5,7),(6,8),(4,9)),derivative_orders=(1,1))/differential.delta**9
    check("general_differential_contraction",contracted,differential.scalar_gradient_squared)
    differential_germ = MonomialMetric((x,y),1,{
        (0,0):[((0,),1+h3)],(0,1):[((0,),h3)],(1,1):[((0,),h3)]})
    dj = differential_germ.leading_through((1,),(1,),0,0,"covariant_riemann",(0,1,0,1,0))
    check("covariant_derivative_finite_jet_order",dj["order"],0)
    check("covariant_derivative_finite_jet_coefficient",dj["coefficient"],-6)

    # Original LEAD7's arbitrary transverse data, now through the same general
    # non-diagonal numerator API.  No fixed [-8,8] Laurent window is used.
    A2,A3,P0,P1,P2,Q0,Q1,Q2 = (s.Function(name)(y,z)
                              for name in ("A2","A3","P0","P1","P2","Q0","Q1","Q2"))
    transverse_generic = MonomialMetric((x,y,z),1,{
        (0,0):[((2,),A2+A3*x)],(1,1):[((0,),P0+P1*x+P2*x*x)],
        (2,2):[((0,),Q0+Q1*x+Q2*x*x)]})
    tg = transverse_generic.leading_through((1,),(1,),3,2)
    check("arbitrary_transverse_generic_order",tg["order"],-3)
    check("arbitrary_transverse_generic_coefficient",tg["coefficient"],(P1/P0+Q1/Q0)/A2)
    B,B4,C1,C10,C2,C20 = (s.Function(name)(y,z)
                         for name in ("B","B4","C1","C10","C2","C20"))
    transverse_parity = MonomialMetric((x,y,z),1,{
        (0,0):[((2,),B+B4*x*x)],(1,1):[((-2,),C1+C10*x*x)],
        (2,2):[((-2,),C2+C20*x*x)]})
    tp = transverse_parity.leading_through((1,),(1,),-10,-2)
    check("arbitrary_transverse_parity_order",tp["order"],-4)
    check("arbitrary_transverse_parity_coefficient",tp["coefficient"],-14/B)

    # Lorentzian null curvature: all scalar contractions vanish, while the
    # lower curvature tensor detects the two opposite tidal directions.
    uu, vv, xx, yy = s.symbols("u v X Y",real=True)
    wave_metric = s.Matrix([[xx**2-yy**2,1,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]])
    wave = CurvatureNumerators(wave_metric,(uu,vv,xx,yy))
    check("null_wave_nonzero_RuXuX",wave.riemann(0,2,0,2),-1)
    check("null_wave_nonzero_RuYuY",wave.riemann(0,3,0,3),1)
    check("null_wave_scalar_zero",wave.scalar)
    for entry in wave.ricci:
        if entry != 0:
            raise AssertionError("Trace-free null wave Ricci tensor")
    checks["null_wave_ricci_zero"] = {"verified":True}
    check("null_wave_kretschmann_zero",wave.kretschmann)

    # Native horizon member, kept genuinely off-diagonal.  A general F(v,r)
    # provides exact scalar and tidal-invariant formulas, including time dependence.
    av, ar, theta, phi = s.symbols("v r theta phi",real=True)
    F = s.Function("F")(av,ar)
    horizon_metric = s.Matrix([[-F,1,0,0],[1,0,0,0],
                               [0,0,ar**2,0],[0,0,0,ar**2*s.sin(theta)**2]])
    horizon = CurvatureNumerators(horizon_metric,(av,ar,theta,phi))
    expected_scalar = -s.diff(F,ar,2)-4*s.diff(F,ar)/ar+2*(1-F)/ar**2
    expected_k = s.diff(F,ar,2)**2+4*s.diff(F,ar)**2/ar**2+4*(1-F)**2/ar**4
    check("off_diagonal_horizon_determinant",horizon.delta,-ar**4*s.sin(theta)**2)
    check("off_diagonal_horizon_scalar",s.trigsimp(horizon.scalar-expected_scalar))
    check("off_diagonal_horizon_kretschmann",s.trigsimp(horizon.kretschmann-expected_k))
    check("off_diagonal_horizon_radial_tensor",horizon.riemann(0,1,0,1),s.diff(F,ar,2)/2)
    check("off_diagonal_horizon_time_dependent_tidal",horizon.riemann(0,2,0,2),
          ar*(F*s.diff(F,ar)-s.diff(F,av))/2)
    rho, L, r0 = s.symbols("rho L r0",positive=True)
    native_F = rho*(2-rho)/(1+2*rho-rho**2)
    substituted_scalar = expected_scalar.subs({s.diff(F,ar,2):s.diff(native_F,rho,2)/L**2,
        s.diff(F,ar):s.diff(native_F,rho)/L,F:native_F},simultaneous=True).subs(ar,r0+L*rho)
    substituted_k = expected_k.subs({s.diff(F,ar,2):s.diff(native_F,rho,2)/L**2,
        s.diff(F,ar):s.diff(native_F,rho)/L,F:native_F},simultaneous=True).subs(ar,r0+L*rho)
    check("native_horizon_finite_scalar",substituted_scalar.subs(rho,0),
          10/L**2-8/(L*r0)+2/r0**2)
    check("native_horizon_finite_kretschmann",substituted_k.subs(rho,0),
          100/L**4+16/(L**2*r0**2)+4/r0**4)
    results["horizon_scalar"] = str(expected_scalar)
    results["horizon_kretschmann"] = str(expected_k)
    results["native_horizon_scalar"] = str(clean(substituted_scalar.subs(rho,0)))
    results["native_horizon_kretschmann"] = str(clean(substituted_k.subs(rho,0)))

    report = {
        "scope":"Exact coordinate tensor and finite-jet extension for a supplied metric; native connection and load selection remain explicit model data.",
        "standard_tools":["coordinate Levi-Civita and curvature definitions","adjugate identity","finite Taylor jets","SymPy exact algebra"],
        "external_sources_used":[],
        "verified_count":len(checks),"checks":checks,"results":results,
        "elapsed_seconds":round(time.monotonic()-start,3)}
    destination = Path(__file__).resolve().parent / "docs/cella-tensor-valuation-checks.json"
    destination.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))


if __name__ == "__main__":
    main()
