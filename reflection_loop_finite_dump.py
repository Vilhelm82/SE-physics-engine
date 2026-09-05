#!/usr/bin/env python3
"""Five native return loops with detuning correction and finite absorbing dumps.

The absorber is off during transport unless background_loss is requested.
No endpoint projection is made. No-click fidelity includes surviving bright
amplitude. All loops have the existing continuous zero-area entry/exit ramps.
"""
import argparse
from functools import lru_cache
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from reflection_loop_composite import Stage, return_time, coefficients, primitive
from reflection_loop_reference_echo import FramedStage, rotation, response_integrals
from reflection_loop_compression import noise_generators, packed_matrix

I = np.eye(4, dtype=complex)
P = np.diag([0., 0., 1., 1.])
Q = I-P
DD, DR = np.diag([0., 0., 1., 0.]), np.diag([0., 0., 0., 1.])
G = np.diag([-1., 1.])
STRETCH = 3*np.sqrt(3)/4


def controls(kind='five'):
    if kind == 'five':
        return (np.pi*np.array([1/6, 5/6, 7/6, 3/4, 1/4]),
                np.array([1, 1, 1, -1, -1]), np.array([1, 1, 1, STRETCH, STRETCH]))
    if kind == 'trine':
        return np.pi*np.array([1/3, 2/3, 1/3]), np.ones(3), np.ones(3)
    raise ValueError(kind)


def primitive_stages(n=1, slew=10.):
    if slew <= 0:
        raise ValueError('positive slew required')
    tau, b = return_time(n), 1/np.sqrt(2)
    if slew < np.pi/(2*tau):
        raise ValueError('slew is below the arc requirement')
    return [Stage(b/slew, q_start=0, q_end=-b),
            Stage((1+b)/slew, q_start=-b, q_end=1)] + [
            Stage(tau, arc=j) for j in range(3)] + [
            Stage((1+b)/slew, q_start=1, q_end=-b),
            Stage(b/slew, q_start=-b, q_end=0)]


def physical_loops(n=1, kind='five', slew=10.):
    angles, signs, stretches = controls(kind)
    base = primitive_stages(n, slew)
    return [[FramedStage(stage, rotation(a), stretch, sign < 0)
             for stage in (base if sign > 0 else reversed(base))]
            for a, sign, stretch in zip(angles, signs, stretches)]


def resources(n=1, kind='five', kappa=10., action=40., slew=10.):
    if kappa <= 0 or action < 0:
        raise ValueError('positive dump rate and nonnegative dump action required')
    angles, signs, stretches = controls(kind)
    td = action/kappa
    loop_time = sum(s.duration for s in primitive_stages(n, slew))*sum(stretches)
    return dict(n=n, kind=kind, loop_count=len(angles), loop_time=float(loop_time),
                dump_time=float(td), dump_rate=kappa, dump_action=action,
                bright_amplitude_factor=float(np.exp(-action/2)),
                total_time=float(loop_time+len(angles)*td),
                nominal_mean_transport_exposure=float(sum(stretches)*coefficients(n)['exposure']/2))


def loop_endpoints(cases, n=1, kind='five', background_loss=0., max_step=.06):
    """Batch finite laboratory ODE; inverse identity is used only at endpoints.

    cases contains (common gain error, d-site frequency, r-site frequency).
    For a physical inverse, additive frequencies change sign before adjoint;
    the real absorbing term keeps its sign. Background loss follows the
    instantaneous bright projector, including its continued ramp support.
    """
    angles, signs, stretches = controls(kind)
    rotations = np.array([rotation(a) for a in angles])
    errors = np.asarray(cases)
    eps = errors[:, 0, None, None, None]
    lab_noise = errors[:, 1, None, None]*DD+errors[:, 2, None, None]*DR
    noise = rotations.transpose(0, 2, 1)[None]@lab_noise[:, None]@rotations[None]
    noise *= (signs*stretches)[None, :, None, None]
    state = np.broadcast_to(I, (len(cases), len(angles), 4, 4)).copy()
    shape = state.shape
    for stage in primitive_stages(n):
        def rhs(t, value):
            generator = -1j*((1+eps)*stage.h(t)+noise)
            if background_loss:
                generator -= background_loss*stretches[None, :, None, None]*stage.loss_projector(t)/2
            return (generator@value.reshape(shape)).ravel()
        solved = solve_ivp(rhs, (0, stage.duration), state.ravel(), method='DOP853',
                           rtol=2e-13, atol=2e-15, max_step=max_step)
        assert solved.success, solved.message
        state = solved.y[:, -1].reshape(shape)
    result = []
    for j, (r, sign) in enumerate(zip(rotations, signs)):
        u = state[:, j]
        if sign < 0:
            u = u.conj().swapaxes(-1, -2)
        result.append(r@u@r.T)
    return np.array(result).transpose(1, 0, 2, 3)


def compose(endpoints, cases, kappa=10., action=40.):
    """Full no-click operator and summed flag effect, with no projection."""
    if kappa <= 0 or action < 0:
        raise ValueError('invalid dump parameters')
    td, r = action/kappa, np.exp(-action/2)
    errors = np.asarray(cases)
    dump = np.broadcast_to(np.diag([r, r, 1., 1.]), (len(cases), 4, 4)).astype(complex).copy()
    dump[:, 2, 2], dump[:, 3, 3] = np.exp(-1j*errors[:, 1:3]*td).T
    state = np.broadcast_to(I, dump.shape).copy()
    flag = np.zeros_like(state)
    for j in range(endpoints.shape[1]):
        u = endpoints[:, j]
        # This term is zero for lossless transport, and includes all its
        # absorbing flags when background_loss is nonzero.
        flag += state.conj().swapaxes(-1, -2)@(I-u.conj().swapaxes(-1, -2)@u)@state
        state = u@state
        flag += (1-r*r)*(state.conj().swapaxes(-1, -2)@Q@state)
        state = dump@state
    return state, flag


def metrics(k, flag=None):
    a, b = G@k[2:, 2:], k[:2, 2:]
    tf = a-np.trace(a)*np.eye(2)/2
    effect = k[:, 2:].conj().T@k[:, 2:]
    extrema = np.linalg.eigvalsh(effect)
    success = float(np.trace(effect).real/2)
    numerator = float(np.linalg.norm(b)**2/2+np.linalg.norm(tf)**2/3)
    loss = 1-success if flag is None else float(np.trace(flag[2:, 2:]).real/2)
    return dict(erasure=loss, success=success,
                conditional_infidelity=numerator/success,
                conditional_bright_population=float(np.linalg.norm(b)**2/(2*success)),
                erasure_spread=float(extrema[1]-extrema[0]),
                equal_input_conditional_bounds=[numerator/float(extrema[1]), numerator/float(extrema[0])])


def laboratory_word(case, n=1, kind='five', kappa=10., action=40., background_loss=0., max_step=.05):
    """Independent chronological ODE including the finite dumps themselves."""
    gain, d, r = case
    noise = d*DD+r*DR
    state = I.copy()
    td = action/kappa
    for loop in physical_loops(n, kind):
        for stage in loop:
            def rhs(t, value):
                generator = -1j*((1+gain)*stage.h(t)+noise)
                generator -= background_loss*stage.loss_projector(t)/2
                return (generator@value.reshape(4, 4)).ravel()
            ans = solve_ivp(rhs, (0, stage.duration), state.ravel(), method='DOP853',
                            rtol=2e-13, atol=2e-15, max_step=max_step)
            assert ans.success
            state = ans.y[:, -1].reshape(4, 4)
        generator = -1j*noise-kappa*Q/2
        ans = solve_ivp(lambda t, v: (generator@v.reshape(4, 4)).ravel(),
                        (0, td), state.ravel(), method='DOP853', rtol=2e-13,
                        atol=2e-15, max_step=min(max_step, .5/kappa))
        assert ans.success
        state = ans.y[:, -1].reshape(4, 4)
    return state


def static_response(n=1, kappa=10., action=40.):
    """Exact nominal dynamics and independent quadrature for all six noises."""
    noises = np.array(noise_generators(), complex)
    td, r = action/kappa, np.exp(-action/2)
    dump = P+r*Q
    dump_first = -1j*td*(dump@noises)
    zero, first = I.copy(), np.zeros_like(noises)
    for loop in physical_loops(n):
        local = response_integrals(loop, noises, nodes=96)
        u, one = local['endpoint'], -1j*local['endpoint']@local['moments']
        first, zero = u@first+one@zero, u@zero
        first, zero = dump@first+dump_first@zero, dump@zero
    return zero, first


@lru_cache(maxsize=None)
def joint_loop_jets(n=1, kind='five', order=4):
    """Reuse the native finite Fourier algebra for all joint gain/d-site orders."""
    from reflection_loop_short_correction import arc_mixed_jet, powers_through, multiply_jets
    from reflection_loop_composite import J4
    from math import factorial
    angles, signs, stretches = controls(kind)
    rotations = np.array([rotation(a) for a in angles])
    noises = rotations.transpose(0, 2, 1)@DD@rotations
    powers = powers_through(order)
    state = np.zeros((len(angles), len(powers), 4, 4), complex)
    state[:, 0] = I
    for stage in primitive_stages(n):
        if stage.arc >= 0:
            block = arc_mixed_jet(stage.arc, noises, n, order)
        else:
            block = np.zeros_like(state)
            area = (stage.q_start+stage.q_end)*stage.duration/2
            u = stage.u(stage.duration)
            block[:, 0] = u
            for k, (p, q) in enumerate(powers[1:], 1):
                if p and q:
                    continue
                block[:, k] = (u@np.linalg.matrix_power(J4, p)*(-1j*area)**p/factorial(p)
                               if p else u@noises*(-1j*stage.duration)**q/factorial(q))
        state = multiply_jets(block, state, powers)
    out = []
    for r, sign, stretch, row in zip(rotations, signs, stretches, state):
        if sign < 0:
            row = row.conj().swapaxes(-1, -2)
        out.append(r@row@r.T*np.array([(sign*stretch)**q for p, q in powers])[:, None, None])
    return np.array(out)


def joint_jet(n=1, kappa=10., action=40., order=4):
    from reflection_loop_short_correction import powers_through, multiply_jets
    from math import factorial
    powers = powers_through(order)
    td, r = action/kappa, np.exp(-action/2)
    dump = np.zeros((len(powers), 4, 4), complex)
    dump[0] = P+r*Q
    for j, (p, q) in enumerate(powers[1:], 1):
        if p == 0:
            dump[j] = (-1j*td)**q/factorial(q)*DD
    result = np.zeros_like(dump)
    result[0] = I
    for loop in joint_loop_jets(n, order=order):
        result = multiply_jets(dump, multiply_jets(loop, result, powers), powers)
    return dict(zip(powers, result))


def exact_gain_word(error, n=1, kappa=10., action=40., ideal=False):
    angles, signs, stretches = controls()
    f = primitive(error, n)
    dump = P if ideal else P+np.exp(-action/2)*Q
    out = I.copy()
    for a, sign in zip(angles, signs):
        rot = rotation(a)
        out = dump@rot@(f if sign > 0 else f.conj().T)@rot.T@out
    return out


def symbolic_checks():
    import sympy as s
    mu, duration, b = s.symbols('mu duration b', real=True)
    g, dd, dr = s.diag(-1, 1), s.diag(1, 0), s.diag(0, 1)
    y = s.Matrix([[0, -s.I], [s.I, 0]])
    def rot(a):
        return s.Matrix([[s.cos(a), -s.sin(a)], [s.sin(a), s.cos(a)]])
    def clean(a):
        return a.applyfunc(s.simplify)
    angles = [s.pi/6, 5*s.pi/6, 7*s.pi/6, 3*s.pi/4, s.pi/4]
    scales = [1, 1, 1, 3*s.sqrt(3)/4, 3*s.sqrt(3)/4]
    prefix, sums, dump, active, signed = s.eye(2), [s.zeros(2), s.zeros(2)], s.zeros(2), s.zeros(2), s.zeros(2)
    for j, (a, scale) in enumerate(zip(angles, scales)):
        r = rot(a); u = clean(r*g*r.T)
        moment = clean(prefix.T*r*dd*r.T*prefix)
        active += moment; signed += (1 if j < 3 else -1)*moment
        c, sn = s.cos(a), s.sin(a)
        for k, raw in enumerate([c*c*mu*dd+sn*sn*duration*dr-c*sn*b*y,
                                 sn*sn*mu*dd+c*c*duration*dr+c*sn*b*y]):
            a0 = clean(r*raw*r.T)
            if j >= 3:
                a0 = clean(u*a0*u)
            sums[k] += scale*prefix.T*a0*prefix
        prefix = clean(u*prefix)
        dump += prefix.T*dd*prefix
    assert clean(prefix-g) == s.zeros(2)
    for a in sums:
        assert clean(a-s.trace(a)*s.eye(2)/2) == s.zeros(2)
    assert clean(dump-s.Rational(5, 2)*s.eye(2)) == s.zeros(2)
    assert clean(active-s.Rational(5, 2)*s.eye(2)) == s.zeros(2)
    assert clean(signed-s.eye(2)/2) == s.zeros(2)
    z = s.symbols('z')
    pair = rot(s.pi/4)*s.diag(z,1)*rot(s.pi/4).T*rot(3*s.pi/4)*s.diag(z,1)*rot(3*s.pi/4).T
    assert clean(pair-z*s.eye(2)) == s.zeros(2)
    def projected_word(angles):
        out=s.eye(2)
        for a in angles:
            r=rot(a);out=clean(r*s.diag(z,1)*r.T*out)
        return out
    swap=s.Matrix([[0,1],[1,0]])
    assert clean(projected_word(angles[:3])-swap*projected_word([s.pi/3,2*s.pi/3,s.pi/3])*swap) == s.zeros(2)
    return dict(nominal_reflection=True, both_code_site_responses_scalar=True,
                dump_response_scalar=True, attenuation_and_phase_moments_scalar=True,
                projected_orthogonal_pair_exactly_scalar=True,
                ideal_conditional_gain_map_equivalent_to_trine=True)


def run_checks():
    checks = {}
    def close(name, value, target=0., tol=1e-9):
        err = float(np.max(abs(np.asarray(value)-target)))
        assert err <= tol, (name, err, tol)
        checks[name] = dict(error=err, tolerance=tol)
    report = dict(symbolic=symbolic_checks(), checks=checks, cases={}, dump_sweep=[], background_sweep=[])
    errors = [(0.,0.,0.), (.01,0.,0.), (.001,0.,0.), (0.,1e-4,0.), (0.,1e-3,0.),
              (.001,1e-4,0.), (.001,-1e-4,0.), (.01,1e-4,0.), (.001,1e-3,0.),
              (.001,0.,1e-4), (.001,1e-4,-1e-4)]
    for n in (1,2,3):
        endpoints = loop_endpoints(errors,n)
        full, flags = compose(endpoints,errors)
        close(f'n{n}_nominal_code_return',full[0][:,2:],np.r_[np.zeros((2,2)),G],2e-11)
        close(f'n{n}_probability_completeness',full.conj().swapaxes(-1,-2)@full+flags,I,2e-12)
        close(f'n{n}_loops_unitary',endpoints.conj().swapaxes(-1,-2)@endpoints,I,2e-11)
        close(f'n{n}_flag_effect_positive',max(0.,-float(np.linalg.eigvalsh(flags).min())),tol=2e-12)
        close(f'n{n}_gain_exact_endpoint',full[1],exact_gain_word(.01,n),3e-12)
        zero, first = static_response(n)
        for j, matrix in enumerate(first):
            a=G@matrix[2:,2:]
            close(f'n{n}_static_generator_{j}_logical',a-np.trace(a)*np.eye(2)/2,tol=2e-10)
        jet=joint_jet(n)
        close(f'n{n}_jet_static_site',jet[0,1],first[2],tol=3e-10)
        close(f'n{n}_jet_nominal',jet[0,0],zero,tol=2e-12)
        for p,q in ((1,0),(0,1)):
            a=G@jet[p,q][2:,2:]
            close(f'n{n}_first_logical_{p}_{q}',a-np.trace(a)*np.eye(2)/2,tol=3e-10)
        # Tiny finite-dump first-order leakage is retained in these predictions.
        for name,case,index in [('joint',errors[5],5),('detuning',errors[3],3)]:
            eps,d,_=case
            predicted=sum(eps**p*d**q*v for (p,q),v in jet.items())
            close(f'n{n}_{name}_fourth_order_endpoint',predicted,full[index],2e-9 if n==1 else 3e-7)
            close(f'n{n}_{name}_joint_prediction',metrics(predicted)['conditional_infidelity']/metrics(full[index])['conditional_infidelity'],1.,.015)
        close(f'n{n}_detuning_quartic_ratio',metrics(full[4])['conditional_infidelity']/metrics(full[3])['conditional_infidelity'],1e4,30.)
        row=resources(n)
        row['errors']=[dict(gain=e,delta_d=d,delta_r=r,**metrics(k,f)) for (e,d,r),k,f in zip(errors,full,flags)]
        row['joint_coefficients']={f'{p},{q}':packed_matrix(v) for (p,q),v in jet.items()}
        row['first_order_bright_gain_coefficient']=float(np.linalg.norm(jet[1,0][:2,2:])**2/2)
        row['first_order_bright_detuning_coefficient']=float(np.linalg.norm(jet[0,1][:2,2:])**2/2)
        report['cases'][f'five_n{n}']=row
        trine_errors=[(.01,0.,0.),(.001,1e-4,0.)]
        tr,fl=compose(loop_endpoints(trine_errors,n,'trine'),trine_errors)
        report['cases'][f'trine_n{n}']=dict(**resources(n,'trine'),errors=[dict(gain=e,delta_d=d,delta_r=r,**metrics(k,f)) for (e,d,r),k,f in zip(trine_errors,tr,fl)])
        # Gain-only ideal projections expose the endpoint-erasure n law.
        tiny=1e-5;b=coefficients(n)['b']
        ideal=exact_gain_word(tiny,n,ideal=True)
        close(f'n{n}_endpoint_erasure_coefficient',metrics(ideal)['erasure']/tiny**2,25*b*b/2,.005)
        gamma=1e-7
        damped,flag=compose(loop_endpoints([(0.,0.,0.)],n,background_loss=gamma),[(0.,0.,0.)])
        exposed=resources(n)['nominal_mean_transport_exposure']
        close(f'n{n}_nominal_background_loss_slope',metrics(damped[0],flag[0])['erasure']/gamma,exposed,2e-5)
        close(f'n{n}_background_loss_scalar_on_code',flag[0][2:,2:]/gamma,exposed*np.eye(2),4e-5)
        if n==1:
            fine=compose(loop_endpoints([errors[5]],n,max_step=.03),[errors[5]])[0][0]
            close('finite_integration_step_halving',fine,full[5],5e-12)
            direct=laboratory_word(errors[5])
            close('independent_chronological_finite_dumps',direct,full[5],5e-12)
            close('conditional_fidelity_convergence',metrics(direct)['conditional_infidelity']/metrics(full[5])['conditional_infidelity'],1.,2e-5)
        print(f'n={n}: finite dumps, six static responses, joint response and return sweep checked',flush=True)
    sweep_errors=[(.001,1e-4,0.),(.01,0.,0.)]
    endpoints=loop_endpoints(sweep_errors)
    for rate in (1.,10.,100.):
        for action in (8.,12.,16.,20.,24.,28.,32.,40.):
            matrices,flags=compose(endpoints,sweep_errors,rate,action)
            report['dump_sweep'].append(dict(**resources(kappa=rate,action=action),errors=[metrics(k,f) for k,f in zip(matrices,flags)]))
    for loss in (1e-6,1e-4,1e-3):
        for n in (1,2,3):
            endpoint=loop_endpoints(sweep_errors,n,background_loss=loss)
            matrices,flags=compose(endpoint,sweep_errors)
            close(f'background_{loss}_n{n}_instrument_completeness',matrices.conj().swapaxes(-1,-2)@matrices+flags,I,2e-12)
            report['background_sweep'].append(dict(background_rate=loss,**resources(n),errors=[metrics(k,f) for k,f in zip(matrices,flags)]))
    direct=laboratory_word(sweep_errors[0],background_loss=1e-4)
    matrix=compose(loop_endpoints([sweep_errors[0]],background_loss=1e-4),[sweep_errors[0]])[0][0]
    close('background_loss_inverse_and_direct_dump_ODE',direct,matrix,6e-12)
    # Physical continuity and limits include every inverse and zero-area ramp.
    loops=physical_loops()
    for j,loop in enumerate(loops):
        close(f'loop{j}_zero_endpoints',[loop[0].h(0),loop[-1].h(loop[-1].duration)],tol=2e-14)
        for k,(left,right) in enumerate(zip(loop,loop[1:])):
            close(f'loop{j}_join_{k}',left.h(left.duration),right.h(0),2e-14)
    report['check_count']=len(checks)
    report['scope']='Five finite native return loops. First-order logical correction for six static generators. No endpoint projection; residual bright amplitude is included. Dump-only and additional transport-loss models are reported separately. No global optimum or hardware claim.'
    return report


def plot(path,report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    ns=np.array([1,2,3])
    for kind,color in [('trine','#976b3b'),('five','#147d72')]:
        rows=[report['cases'][f'{kind}_n{n}'] for n in ns]
        axes[0,0].plot(ns,[100*r['errors'][1 if kind=='five' else 0]['erasure'] for r in rows],'o-',color=color,label=kind)
        axes[0,1].semilogy([r['total_time'] for r in rows],[r['errors'][5 if kind=='five' else 1]['conditional_infidelity'] for r in rows],'o-',color=color,label=kind)
    for rate,color in [(1.,'#976b3b'),(10.,'#147d72'),(100.,'#546da8')]:
        rows=[r for r in report['dump_sweep'] if r['dump_rate']==rate]
        axes[1,0].semilogy([r['total_time'] for r in rows],[r['errors'][0]['conditional_infidelity'] for r in rows],'o-',color=color,label=rf'$\kappa/a={rate:g}$')
    for n,color in zip(ns,['#147d72','#546da8','#976b3b']):
        rows=[r for r in report['background_sweep'] if r['n']==n]
        axes[1,1].loglog([r['background_rate'] for r in rows],[r['errors'][1]['erasure'] for r in rows],'o-',color=color,label=f'n={n}')
    axes[0,0].set(xlabel='Return index n',ylabel='Flag probability (%)',title=r'Endpoint flags at 1% gain, $\Delta=0$',xticks=ns)
    axes[0,1].set(xlabel=r'Total duration $aT$',ylabel='Conditional infidelity',title=r'$\epsilon=10^{-3},\ \Delta_d/a=10^{-4}$')
    axes[1,0].set(xlabel=r'Total duration $aT$',ylabel='Conditional infidelity',title='Finite dump action: 8 to 40; n=1')
    axes[1,1].set(xlabel=r'Additional transport loss rate $\gamma/a$',ylabel='Flag probability',title='Where higher n helps: background loss, 1% gain')
    for ax in axes.flat:
        ax.grid(alpha=.2);ax.legend(frameon=False)
    fig.suptitle('Detuning-corrected reflection with five finite dumps')
    fig.savefig(path,dpi=170)
    plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json');parser.add_argument('--plot')
    args=parser.parse_args();report=run_checks()
    if args.json:Path(args.json).write_text(json.dumps(report,indent=2)+'\n')
    else:print(json.dumps(report,indent=2))
    if args.plot:plot(args.plot,report)
