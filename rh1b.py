#!/usr/bin/env python3
"""RH-1b: fixed calibration, complete click enumeration, explicit replay loss."""
import numpy as np
from scipy.linalg import polar
from rh1_common import (G,I2,word_loops,foldback,pattern_kraus,channel_metrics,
    replay_branches,compose_discard,receipt,check)


def run():
    rows, calibrations, checks, summary = [], {}, [], {}
    for kind in ('trine','five'):
        eps_values=[.005,.01,.02]
        loops=word_loops(kind,[(eps,0.,0.) for eps in eps_values])
        cal=loops[1]
        cs=[foldback(l[2:,2:],l[:2,2:]) for l in cal]
        us=[polar(l[:2,2:])[0] for l in cal]
        calibrations[kind]=[dict(dump=j+1,U=u,V=polar(l[2:,2:])[0],C=c)
                            for j,(l,c,u) in enumerate(zip(cal,cs,us))]
        for eps,ls in zip(eps_values,loops):
            k,fs,discard=compose_discard(ls)
            branches=[(l[2:,2:],c@l[:2,2:]) for l,c in zip(ls,cs)]
            recovered=channel_metrics(pattern_kraus(branches))
            check(checks,f'{kind} eps={eps}: recovery TP',recovered['tp_defect'],3e-12)
            check(checks,f'{kind} eps={eps}: independent fidelity formula',
                  recovered['infidelity']-recovered['direct_infidelity'],2e-14)
            # Six Pauli eigenstates form a qubit 2-design: an independent
            # density-operator check of the complete recovered channel.
            states=np.array([[1,0],[0,1],[1,1],[1,-1],[1,1j],[1,-1j]],complex)
            states/=np.linalg.norm(states,axis=1)[:,None]
            fidelities=[]
            for psi in states:
                rho=np.outer(psi,psi.conj())
                for choices in branches:
                    rho=sum(b@rho@b.conj().T for b in choices)
                fidelities.append(np.vdot(G@psi,rho@(G@psi)).real)
            check(checks,f'{kind} eps={eps}: density evolution',
                  1-np.mean(fidelities)-recovered['infidelity'],3e-14)
            replay=channel_metrics(pattern_kraus([replay_branches(l[2:,2:],l[:2,2:],u)
                                                 for l,u in zip(ls,us)]))
            check(checks,f'{kind} eps={eps}: replay TP',replay['tp_defect'],4e-12)
            incomplete=pattern_kraus([(l[2:,2:],l[2:,2:]@u.conj().T@l[:2,2:])
                                     for l,u in zip(ls,us)])
            replay['without_second_leak_mean_trace']=channel_metrics(incomplete)['mean_trace']
            replay['difference_from_foldback']=replay['infidelity']-recovered['infidelity']
            # The purely gain-driven native loop has an exact spectator.
            # Its local error after polar recovery is dephasing with factor
            # sqrt(1-ell); it cannot be an error-free scaled unitary.
            local_distortion=[]
            for j,(l,c) in enumerate(zip(ls,cs),1):
                v=polar(l[2:,2:])[0]
                lm=channel_metrics([l[2:,2:],c@l[:2,2:]],target=v)
                local_distortion.append(dict(dump=j,infidelity=lm['infidelity']))
            rows.append(dict(word=kind,eps=eps,discard=discard,recovered=recovered,
                             replay=replay,per_dump=local_distortion,
                             ratio_to_conditional=recovered['infidelity']/discard['conditional'],
                             ratio_to_unconditional=recovered['infidelity']/discard['unconditional']))
        these=[r for r in rows if r['word']==kind]
        mid=these[1]
        fork='a' if mid['ratio_to_conditional']<=10 else (
             'b' if .1<=mid['ratio_to_unconditional']<=10 else 'between')
        excess=[r['recovered']['infidelity']-mid['recovered']['infidelity'] for r in (these[0],these[2])]
        exponent=(float(np.log(excess[1]/excess[0])/np.log(2)) if min(excess)>0 else None)
        raw_power=float(np.polyfit(np.log(eps_values),np.log([r['recovered']['infidelity'] for r in these]),1)[0])
        summary[kind]=dict(fork=fork,calibrated_infidelity=mid['recovered']['infidelity'],
            ratio_to_discard_conditional=mid['ratio_to_conditional'],
            ratio_to_discard_unconditional=mid['ratio_to_unconditional'],
            off_calibration_excess=excess,excess_exponent=exponent,
            excess_fit_status='not a positive excess law around calibration' if exponent is None else 'two-distance diagnostic only',
            raw_gain_exponent=raw_power)
    return dict(rows=rows,calibrations=calibrations,checks=checks,summary=summary,
        replay_convention='Calibrated U fixed; current S is the physically replayed loop. A second leak is explicitly replaced by I/2. The incomplete two-branch trace is also reported.',
        scope='Ideal coherent capture only. Foldback kill from RH-1a retained; runs quantify residual distortion.')


if __name__=='__main__':
    receipt('rh1b',run())
