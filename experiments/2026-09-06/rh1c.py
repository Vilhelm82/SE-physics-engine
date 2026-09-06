#!/usr/bin/env python3
"""RH-1c: signed timing response, local-fit stability and Bernoulli resources."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))  # repo root on sys.path (reorg 2026-09-06)
from functools import lru_cache
import numpy as np
from scipy.optimize import minimize_scalar
from rlq.rh1_common import (TAU,TRINE,I2,trine_loop,db,receipt,check)

DELTAS=np.array([-1e-3,-3e-4,-1e-4,0.,1e-4,3e-4,1e-3])


@lru_cache(None)
def leak(eps,delta,dt):
    before=I2.copy()
    for a in TRINE[:-1]:
        before=trine_loop(float(a),eps,delta)[2:,2:]@before
    e=trine_loop(float(TRINE[-1]),eps,delta,0.,dt)[:2,2:]
    return e@before, e, float(np.linalg.norm(before)**2/2)


def rate(eps,delta,dt):
    return float(np.linalg.norm(leak(eps,delta,dt)[0])**2/2)


def alignment(a,b):
    denom=np.linalg.norm(a)*np.linalg.norm(b)
    if denom<1e-12:
        return dict(complex_alignment=None,real_alignment=None,zero_tangent=True)
    inner=np.vdot(a,b)
    return dict(complex_alignment=float(abs(inner)/denom),
                real_alignment=float(inner.real/denom),zero_tangent=False)


def quadratic_window(eps,delta,width):
    xs=np.array([-1.,-.5,0.,.5,1.])
    ys=np.array([rate(eps,delta,float(x*width)) for x in xs])
    baseline=ys[2]
    coeff=np.polynomial.polynomial.polyfit(xs,ys-baseline,2)
    vertex=float(-coeff[1]/(2*coeff[2])*width) if coeff[2]>0 else None
    return dict(width=width,rates=ys,coefficients=coeff,
        vertex=vertex,curvature=float(2*coeff[2]/width**2),
        relative_fit_residual=float(np.linalg.norm(ys-baseline-np.polynomial.polynomial.polyval(xs,coeff))/max(np.linalg.norm(ys-baseline),1e-30)))


def run():
    rows, tangents, summaries, checks, resources=[],[],{},[],[]
    paired_resources=[]
    timing_convergence=[]
    widths=[.02*TAU,.05*TAU]
    for eps in (.01,0.):
        print(f'RH-1c: epsilon={eps}',flush=True)
        p0s=[]
        asym={w:[] for w in widths}
        for delta in DELTAS:
            delta=float(delta)
            p0=rate(eps,delta,0.)
            p0s.append(p0)
            fits=[quadratic_window(eps,delta,w) for w in widths+[.002*TAU]]
            dither=[]
            for w in widths:
                pm,pp=rate(eps,delta,-w),rate(eps,delta,w)
                asym[w].append(pp-pm)
                dither.append(dict(dt=w,p_minus=pm,p_plus=pp,asymmetry=pp-pm,
                    conditional_on_reaching_minus=pm/leak(eps,delta,-w)[2],
                    conditional_on_reaching_plus=pp/leak(eps,delta,w)[2]))
            # Q H(T) P = 0 makes p'(0)=0 exactly; a wide-window quadratic
            # vertex is not automatically the physical stationary point.
            f=leak(eps,delta,0.)[0]
            beta=-1j*np.array([[0.,1+eps],[1+eps,0.]])@f
            slope=float(np.vdot(f,beta).real)
            check(checks,f'eps={eps} delta={delta}: endpoint rate slope',slope,1e-16)
            candidates=[0.]+[f['vertex'] for f in fits if f['vertex'] is not None and abs(f['vertex'])<=.05*TAU]
            # Search each side separately: the stationary endpoint can lie
            # between very shallow off-centre extrema at zero gain.
            for bounds in [(-.05*TAU,0.),(0.,.05*TAU)]:
                opt=minimize_scalar(lambda t: rate(eps,delta,float(t)),bounds=bounds,
                                    method='bounded',options={'xatol':1e-8})
                if not opt.success:
                    raise RuntimeError(opt.message)
                candidates.extend([float(opt.x),*bounds])
            ps=[rate(eps,delta,float(t)) for t in candidates]
            ib=int(np.argmin(ps))
            rows.append(dict(eps=eps,delta=delta,p_nominal=p0,dithers=dither,local_quadratic_fits=fits,
                             tested_best_dt=candidates[ib],tested_best_rate=ps[ib],
                             timing_search_interval=[-.05*TAU,.05*TAU],
                             fraction_removed=max(0.,(p0-ps[ib])/p0) if p0>1e-28 else None,
                             endpoint_derivative=slope))
        selected=next(r for r in rows if r['eps']==eps and r['delta']==1e-3)
        fine_rates=[]
        for dt in (0.,selected['tested_best_dt']):
            before=I2.copy()
            for i,a in enumerate(TRINE):
                l=trine_loop(float(a),eps,1e-3,0.,dt if i==2 else 0.,max_step=.06)
                if i==2:
                    fine_rates.append(float(np.linalg.norm(l[:2,2:]@before)**2/2))
                before=l[2:,2:]@before
        coarse=np.array([selected['p_nominal'],selected['tested_best_rate']])
        # The scientifically small quantity is the removed depth. A common
        # ODE offset in both rates must not be mistaken for a timing benefit.
        check(checks,f'eps={eps}: minimum-depth ODE step halving',
              (coarse[0]-coarse[1])-(fine_rates[0]-fine_rates[1]),1e-18)
        timing_convergence.append(dict(eps=eps,delta=1e-3,dt=selected['tested_best_dt'],
            coarse_rates=coarse,fine_rates=fine_rates,
            coarse_improvement=float(coarse[0]-coarse[1]),
            fine_improvement=float(fine_rates[0]-fine_rates[1])))
        p0s=np.array(p0s)
        # All seven points determine baseline, linear and quadratic terms.
        x=DELTAS/1e-3
        coeff=np.polynomial.polynomial.polyfit(x,p0s,2)
        even_excess=(p0s[:3]+p0s[:3:-1])/2-p0s[3]
        exponent=float(np.polyfit(np.log(abs(DELTAS[:3])),np.log(even_excess),1)[0]) if np.all(even_excess>0) else None
        laws=[]
        for w,ys in asym.items():
            ys=np.array(ys)
            fit=np.polynomial.polynomial.polyfit(x,ys,2)
            odd=(ys[4:]-ys[2::-1])/2
            laws.append(dict(dt=w,intercept=float(fit[0]),linear_slope=float(fit[1]/1e-3),
                quadratic_coefficient=float(fit[2]/1e-6),odd_max=float(np.max(abs(odd))),
                values=ys,relative_fit_residual=float(np.linalg.norm(ys-np.polynomial.polynomial.polyval(x,fit))/max(np.linalg.norm(ys),1e-30))))
        tangent_records=[]
        for h in (1e-5,5e-6):
            a=(leak(eps,h,0.)[1]-leak(eps,-h,0.)[1])/(2*h)
            ac=(leak(eps,h,0.)[0]-leak(eps,-h,0.)[0])/(2*h)
            e=leak(eps,0.,0.)[1]
            f=leak(eps,0.,0.)[0]
            b=-1j*np.array([[0.,1+eps],[1+eps,0.]])@e
            bc=-1j*np.array([[0.,1+eps],[1+eps,0.]])@f
            # Matrix columns are amplitude vectors on |d>, |r>; their
            # Frobenius pairing/2 is the exact Haar-averaged pairing.
            entry=dict(eps=eps,step=h,alpha=a,beta=b,alpha_cumulative=ac,beta_cumulative=bc,
                alpha_haar_squared=float(np.linalg.norm(a)**2/2),
                beta_haar_squared=float(np.linalg.norm(b)**2/2),local=alignment(a,b),
                cumulative=alignment(ac,bc))
            tangent_records.append(entry)
            tangents.append(entry)
        check(checks,f'eps={eps}: detuning derivative convergence',
              np.linalg.norm(tangent_records[0]['alpha']-tangent_records[1]['alpha']),1e-7)
        for w in widths:
            for sign in (-1,1):
                p=rate(eps,0.,sign*w)
                deriv=(rate(eps,1e-5,sign*w)-rate(eps,-1e-5,sign*w))/(2e-5)
                fisher=deriv**2/(p*(1-p)) if 0<p<1 else 0.
                resources.append(dict(eps=eps,dt=sign*w,p=p,derivative_at_zero=deriv,
                    fisher_at_zero=fisher,shots_for_1e_minus5=(1/(1e-10*fisher) if fisher>1e-20 else None),
                    DB_dither=db(leak(eps,0.,sign*w)[1]),
                    interpretation='Formal local Bernoulli bound for known Haar input ensemble; not a resolved signed estimator if derivative is numerical zero.'))
            ps=np.array([rate(eps,0.,sign*w) for sign in (-1,1)])
            ddelta=np.array([(rate(eps,1e-5,sign*w)-rate(eps,-1e-5,sign*w))/(2e-5) for sign in (-1,1)])
            deps=np.array([(rate(eps+1e-5,0.,sign*w)-rate(eps-1e-5,0.,sign*w))/(2e-5) for sign in (-1,1)])
            scores=np.stack((ddelta,deps),axis=1)
            fi=scores.T@(scores/(2*ps[:,None]*(1-ps[:,None])))
            weighted=scores/np.sqrt(2*ps[:,None]*(1-ps[:,None]))
            nuisance=weighted[:,1]
            effective=float(np.linalg.norm(weighted[:,0]-nuisance*np.dot(nuisance,weighted[:,0])/np.dot(nuisance,nuisance))**2) if np.dot(nuisance,nuisance)>0 else float(fi[0,0])
            asym_fi=float((ddelta[1]-ddelta[0])**2/(2*np.sum(ps*(1-ps))))
            paired_resources.append(dict(eps=eps,dt=w,full_fisher=fi,
                known_gain_fisher=float(fi[0,0]),
                known_gain_shots=float(1e10/fi[0,0]),
                asymmetry_only_fisher=asym_fi,asymmetry_only_shots=float(1e10/asym_fi),
                unknown_gain_effective_fisher=effective,
                unknown_gain_shots=float(1e10/effective) if effective>1e-20 else None))
        summaries[str(eps)]=dict(p_nominal=p0s,
            seven_point_quadratic_coefficients=[float(coeff[0]),float(coeff[1]/1e-3),float(coeff[2]/1e-6)],
            quadratic_fit_relative_residual=float(np.linalg.norm(p0s-np.polynomial.polynomial.polyval(x,coeff))/max(np.linalg.norm(p0s-p0s[3]),1e-30)),
            even_excess_exponent=exponent,asymmetry_laws=laws,
            alignment=tangent_records[-1]['local'],
            fork=('c: signed response, no resolved minimum shift at 1% gain' if eps else
                  'c: signed response, only a shallow higher-order minimum shift at zero gain'),
            diagnostic_kill=False)
    return dict(summary=summaries,rows=rows,tangents=tangents,bernoulli_resources=resources,
        paired_resources=paired_resources,timing_search_convergence=timing_convergence,checks=checks,
        probability_convention='p_N is final first-click probability per original input, Haar averaged; also report conditioned-on-reaching rates.',
        fit_warning='At zero gain beta vanishes. A quadratic fit of a quartic timing minimum has no first-order window-shift interpretation.')


if __name__=='__main__':
    receipt('rh1c',run())
