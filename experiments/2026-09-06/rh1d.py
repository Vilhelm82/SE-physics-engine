#!/usr/bin/env python3
"""RH-1d: virtual frame audit and independent laboratory integration of a lock.

The literal positive-phase specification is run as written. The corrected
negative phase and common-mode rate are run separately and labelled.
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))  # repo root on sys.path (reorg 2026-09-06)
import numpy as np
from scipy.linalg import polar
from rlq.rh1_common import (I2,G,Z,TAU,TRINE,word_loops,trine_loop,compose_discard,
                        conditional_metrics,receipt,check)


def z_character(w):
    denom=np.linalg.norm(w-I2)
    return float(np.linalg.norm(w@Z-Z@w)/denom) if denom>1e-14 else None


def locked_loops(eps,dd,dr,lock,sign=-1,dt=0.):
    lab,rotating=[],[]
    for j,a in enumerate(TRINE):
        start=j*3*TAU
        end=start+3*TAU+(dt if j==2 else 0.)
        l=trine_loop(float(a),eps,dd,dr,dt if j==2 else 0.,tuple(lock),sign,start)
        r0=np.diag(np.exp(sign*1j*start*np.array(lock)))
        r1=np.diag(np.exp(sign*1j*end*np.array(lock)))
        lab.append(l)
        rotating.append(r1.conj().T@l@r0)
    return np.array(lab),np.array(rotating)


def last_rate(loops):
    return float(np.linalg.norm(compose_discard(loops)[1][-1])**2/2)


def run():
    virtual,physical,residuals,checks,summary=[],[],[],[],{}
    for kind in ('trine','five'):
        cases=[(0.,1e-4,0.)]
        if kind=='five': cases += [(0.,0.,1e-4),(0.,1e-4,1e-4),(0.,1e-4,-1e-4)]
        for case,loops in zip(cases,word_loops(kind,cases)):
            k,_,discard=compose_discard(loops)
            v,j=polar(k)
            unit=conditional_metrics(v)['conditional']
            nonunit=conditional_metrics(j,target=I2)['conditional']
            fraction=unit/discard['conditional']
            w=v@G.conj().T
            wp=np.exp(-1j*np.angle(np.trace(w)))*w
            raw,projective=z_character(w),z_character(wp)
            fork=('a' if fraction>.99 and raw is not None and raw<1e-2 else
                  'b' if fraction>.99 else 'c')
            check(checks,f'{kind} {case}: polar reconstruction',np.linalg.norm(v@j-k),2e-14)
            check(checks,f'{kind} {case}: virtual unitary',np.linalg.norm(v.conj().T@v-I2),3e-14)
            virtual.append(dict(word=kind,case=case,delta_com=(case[1]+case[2])/2,
                delta_diff=case[1]-case[2],conditional=discard['conditional'],unitary_infidelity=unit,
                unitary_fraction=fraction,positive_polar_infidelity=nonunit,
                residual_after_full_unitary_readout=nonunit,raw_Zc=raw,phase_removed_Zc=projective,
                global_phase=float(np.angle(np.trace(w))),fork_using_spec_metric=fork,
                V=v,J=j,W=w))
    delta=1e-4
    for mode in ('differential','common'):
        dd,dr=(delta/2,-delta/2) if mode=='differential' else (delta,delta)
        generator=np.array([0.,0.,.5,-.5 if mode=='differential' else .5])
        matched_hat=delta if mode=='differential' else 2*delta
        rows=[]
        for eps in (0.,.01):
            nominal=word_loops('trine',[(eps,0.,0.)])[0]
            nk,_,nm=compose_discard(nominal)
            for name,sign,hat in [('literal_spec',1,delta),
                                  ('negative_phase_same_rate',-1,delta),
                                  ('matched_negative_phase',-1,matched_hat)]:
                lock=tuple(hat*generator)
                lab,rotating=locked_loops(eps,dd,dr,lock,sign)
                # From U=R V: H_rot = H0 + D + sign*hat*generator.
                residual=(dd+sign*lock[2],dr+sign*lock[3])
                reference=word_loops('trine',[(eps,*residual)])[0]
                identity=float(np.linalg.norm(rotating-reference))
                check(checks,f'{mode} eps={eps} {name}: laboratory gauge identity',identity,1e-12)
                lk,_,lm=compose_discard(lab)
                rk,_,rm=compose_discard(rotating)
                row=dict(mode=mode,eps=eps,protocol=name,phase_sign=sign,delta_hat=hat,
                    site_residual=residual,gauge_defect=identity,
                    laboratory_error_from_nominal=float(np.linalg.norm(lk-nk)),
                    rotating_error_from_nominal=float(np.linalg.norm(rk-nk)),
                    laboratory=lm,rotating=rm,nominal=nm,
                    erasure_relative_change=(rm['erasure']-nm['erasure'])/nm['erasure'] if nm['erasure']>1e-20 else None)
                physical.append(row)
                rows.append(row)
                if name=='matched_negative_phase':
                    check(checks,f'{mode} eps={eps}: exact-lock erasure',rm['erasure']-nm['erasure'],1e-13)
                    dither=[]
                    for dt in (-.05*TAU,-.02*TAU,0.,.02*TAU,.05*TAU):
                        _,locked=locked_loops(eps,dd,dr,lock,sign,dt)
                        nominal_dt=nominal.copy()
                        nominal_dt[-1]=trine_loop(float(TRINE[-1]),eps,0.,0.,dt)
                        actual,expected=last_rate(locked),last_rate(nominal_dt)
                        check(checks,f'{mode} eps={eps} dt={dt}: window profile restoration',actual-expected,2e-14)
                        dither.append(dict(dt=dt,locked_rate=actual,nominal_rate=expected))
                    row['dither_profile']=dither
        # More than the specified +/-10% pair is needed to fit a power.
        # Use zero gain to isolate detuning; the gain floor is reported above.
        offsets=np.array([-.2,-.1,-.05,-.025,.025,.05,.1,.2])
        values=[]
        for fraction in offsets:
            hat=matched_hat*(1+fraction)
            _,ls=locked_loops(0.,dd,dr,tuple(hat*generator),-1)
            _,_,metrics=compose_discard(ls)
            values.append(metrics['conditional'])
            residuals.append(dict(mode=mode,delta=delta,delta_hat=hat,fraction=fraction,
                residual_rate=matched_hat-hat,conditional=metrics['conditional']))
        power=float(np.polyfit(np.log(abs(offsets*matched_hat)),np.log(values),1)[0])
        summary[mode]=dict(residual_power=power,
            literal_spec_returns_to_nominal=False,
            matched_negative_phase_restores_profile=True,
            required_phase='exp(-i delta_hat t (D_d +/- D_r)/2)',
            required_delta_hat=matched_hat)
    summary['virtual']=[{k:r[k] for k in ('word','case','unitary_fraction','raw_Zc','phase_removed_Zc','fork_using_spec_metric')}
                        for r in virtual]
    return dict(summary=summary,virtual=virtual,physical=physical,residual_sweep=residuals,checks=checks,
        model_extension='Complex Hermitian couplings R H0 R† replace the original real-symmetric control family. No hardware realization or stochastic tracking tested.',
        corrections=['Literal plus phase adds detuning in the rotating frame.',
            'Equality to nominal is in R† U, not in the fixed laboratory frame.',
            'For delta_d=delta_r=delta and generator (D_d+D_r)/2, exact cancellation needs delta_hat=2 delta.',
            'Raw Zc is global-phase dependent; phase-removed Zc is reported separately.',
            'Window restoration means the full epsilon-only rate profile, whose broad quadratic-fit vertex need not be zero.'])


if __name__=='__main__':
    receipt('rh1d',run())
