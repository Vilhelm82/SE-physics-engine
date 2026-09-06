#!/usr/bin/env python3
"""RH-1a: prescribed per-dump and cumulative flag tests; no fitted controls."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))  # repo root on sys.path (reorg 2026-09-06)
import numpy as np
from scipy.linalg import polar
from rlq.rh1_common import (I2, word_loops, compose_discard, db, haar_quadrature,
                        receipt, check)

CASES = [(1e-2,0.,0.), (1e-3,0.,0.), (0.,1e-4,0.), (0.,0.,1e-4),
         (0.,1e-4,1e-4), (1e-2,1e-4,0.)]


def run():
    rows, families, checks, summary = [], [], [], {}
    states, weights = haar_quadrature()
    for kind in ('trine','five'):
        endpoints = word_loops(kind,CASES)
        for case, loops in zip(CASES,endpoints):
            _, fs, discard = compose_discard(loops)
            dumps = []
            for k,(l,f) in enumerate(zip(loops,fs),1):
                e=l[:2,2:]
                info=db(e)
                probabilities=np.sum(abs(states@e.T)**2,axis=1)
                mean=float(weights@probabilities)
                variance=float(weights@((probabilities-mean)**2))
                check(checks,f'{kind} {case} dump {k}: Haar mean',mean-info['mean'],1e-15)
                if info['DB'] is not None:
                    # Exact Haar variance/mean^2 = DB^2/3, unlike a finite
                    # sample's extrema, which depend on axis coverage.
                    check(checks,f'{kind} {case} dump {k}: Haar contrast',
                          np.sqrt(3*variance)/mean-info['DB'],1e-10)
                check(checks,f'{kind} {case} dump {k}: unitary',
                      np.linalg.norm(l.conj().T@l-np.eye(4)),3e-11)
                info.update(dump=k,cumulative=db(f),
                            haar_spread_over_mean=float(np.ptp(probabilities)/mean),
                            exact_spread_over_mean=2*info['DB'] if info['DB'] is not None else None)
                if case[1]==case[2]==0:
                    info['U']=polar(e)[0]
                    info['E']=e
                dumps.append(info)
            check(checks,f'{kind} {case}: total probability',discard['completeness_defect'],3e-11)
            rows.append(dict(word=kind,case=case,dumps=dumps,discard=discard))
        gainrows=[r for r in rows if r['word']==kind and r['case'][1:]==(0.,0.)]
        vals=[d['DB'] for r in gainrows for d in r['dumps']]
        summary[kind]=dict(fork='a' if max(vals)<1e-4 else 'b',
                           foldback_kill=all(v>1e-2 for v in vals),DB_gain_range=[min(vals),max(vals)])
        h=1e-5
        for step in (h,h/2):
            ls=word_loops(kind,[(step,0.,0.),(-step,0.,0.),(0.,step,0.),(0.,-step,0.)])
            for k in range(ls.shape[1]):
                a=(ls[0,k,:2,2:]-ls[1,k,:2,2:])/(2*step)
                b=(ls[2,k,:2,2:]-ls[3,k,:2,2:])/(2*step)
                families.append(dict(word=kind,dump=k+1,step=step,A=a,B=b,
                    Omega=float(np.linalg.norm(a.conj().T@b)/(np.linalg.norm(a)*np.linalg.norm(b)))))
    return dict(summary=summary,rows=rows,family_overlap=families,checks=checks,
                polar_note='Rank-deficient E uses a unitary polar completion. Null-space columns are arbitrary.')


if __name__=='__main__':
    receipt('rh1a',run())
