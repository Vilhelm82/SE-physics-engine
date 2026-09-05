#!/usr/bin/env python3
"""Circuit QEC tests of the finite-dump instrument with explicit boundary repair.

The tested instrument is Pauli randomized independently at each location.
Code states survive the boundary repair; all non-code states become I/2.
Syndrome CNOTs, preparation and measurement have separately declared noise.
Heralds are sampled in the circuit, then used to reweight the decoder per shot.
No detector-error-model sampling of heralded channels is used.
"""
import argparse
from collections import OrderedDict, defaultdict
from dataclasses import dataclass
import hashlib
import itertools
import json
from pathlib import Path
import time

import numpy as np
from scipy import sparse
from scipy.stats import beta
import stim
import pymatching
from ldpc import BpOsdDecoder

ROOT = Path(__file__).resolve().parent
PAULIS = np.array([np.eye(2), [[0, 1], [1, 0]], [[0, -1j], [1j, 0]],
                   [[1, 0], [0, -1]]], dtype=complex)
G = np.diag([-1., 1.])


def unpack(value):
    return np.array(value['real']) + 1j*np.array(value['imag'])


def operating_point(n=1, gain=.01, background=1e-6):
    path = ROOT/'docs/reflection-loop-finite-dump-checks.json'
    rows = json.loads(path.read_text())['operating_curves']
    return next(r for r in rows if r['n'] == n and r['gain'] == gain
                and np.isclose(r['background_rate'], background, rtol=1e-9, atol=1e-20))


def native_channel(row):
    """Exact Pauli twirl of the repaired, two-outcome instrument, not a fit to F.

    Phi_h(rho) = tr(E_h rho) I/2.
    Phi_0(rho) = A rho A^dagger + tr((E_u+B^dagger B) rho) I/2.
    w_P = |tr(P A)/2|^2 + tr(E_u+B^dagger B)/8.
    """
    k, eh, eu = map(unpack, (row['survival_code_map'], row['dump_effect'], row['background_effect']))
    a, b = G@k[2:], k[:2]
    lost = eu+b.conj().T@b
    h = float(np.trace(eh).real/2)
    r = float(np.trace(lost).real/2)
    weights = np.abs(np.einsum('pij,ji->p', PAULIS, a)/2)**2+r/4
    assert np.linalg.norm(k.conj().T@k+eh+eu-np.eye(2)) < 3e-11
    assert abs(weights.sum()+h-1) < 3e-11
    assert np.linalg.eigvalsh(lost).min() > -1e-13
    conditional = weights/(1-h)
    conditional[0] = 1-conditional[1:].sum()
    assert np.all(conditional >= 0)
    return dict(herald=h, pauli=conditional.tolist(), replaced_unheralded=r,
                repaired_conditional_average_infidelity=float(2*conditional[1:].sum()/3),
                original_conditional_average_infidelity=row['unheralded_residual'],
                n=row['n'], gain=row['gain'], background_rate=row['background_rate'],
                delta_d=row['delta_d'], delta_r=row['delta_r'], total_time=row['total_time'],
                dump_rate=row['dump_rate'], dump_action=row['dump_action'])


def verify_twirl(row):
    k, eh, eu = map(unpack, (row['survival_code_map'], row['dump_effect'], row['background_effect']))
    a, b = G@k[2:], k[:2]
    e = eu+b.conj().T@b
    ch = native_channel(row)
    errors = []
    for rho in PAULIS:
        direct = sum(p@(a@(p@rho@p)@a.conj().T+
                        np.trace(e@p@rho@p)*np.eye(2)/2)@p for p in PAULIS)/4
        expected = (1-ch['herald'])*sum(w*p@rho@p for w, p in zip(ch['pauli'], PAULIS))
        flag = sum(np.trace(eh@p@rho@p)*np.eye(2)/2 for p in PAULIS)/4
        errors += [np.linalg.norm(direct-expected),
                   np.linalg.norm(flag-ch['herald']*np.trace(rho)*np.eye(2)/2)]
    return float(max(errors))


def rref(a):
    a = np.array(a, dtype=np.uint8).copy()
    pivots, row = [], 0
    for col in range(a.shape[1]):
        candidates = np.flatnonzero(a[row:, col])
        if not len(candidates):
            continue
        chosen = row+int(candidates[0])
        a[[row, chosen]] = a[[chosen, row]]
        others = np.flatnonzero(a[:, col])
        others = others[others != row]
        a[others] ^= a[row]
        pivots.append(col)
        row += 1
        if row == a.shape[0]:
            break
    return a[:row], pivots


def kernel(a):
    reduced, pivots = rref(a)
    free = [i for i in range(a.shape[1]) if i not in pivots]
    out = np.zeros((len(free), a.shape[1]), dtype=np.uint8)
    for j, col in enumerate(free):
        out[j, col] = 1
        out[j, pivots] = reduced[:, col]
    return out


def quotient_basis(candidates, stabilizers):
    basis, pivots = rref(stabilizers)
    selected = []
    for v in candidates:
        residue = v.copy()
        for b, p in zip(basis, pivots):
            if residue[p]:
                residue ^= b
        if np.any(residue):
            selected.append(v.copy())
            basis, pivots = rref(np.vstack([basis, residue]))
    return np.array(selected, dtype=np.uint8)


def bicycle(ell=6, m=6):
    x = np.kron(np.roll(np.eye(ell, dtype=np.uint8), 1, axis=1), np.eye(m, dtype=np.uint8))
    y = np.kron(np.eye(ell, dtype=np.uint8), np.roll(np.eye(m, dtype=np.uint8), 1, axis=1))
    at = [np.linalg.matrix_power(x, 3), y, y@y]
    bt = [np.linalg.matrix_power(y, 3), x, x@x]
    a, b = np.bitwise_xor.reduce(at), np.bitwise_xor.reduce(bt)
    hx, hz = np.hstack([a, b]), np.hstack([b.T, a.T])
    lx, lz = quotient_basis(kernel(hz), hx), quotient_basis(kernel(hx), hz)
    assert not np.any(hx@hz.T % 2)
    assert len(lx) == len(lz) == 12
    assert len(rref(lx@lz.T % 2)[1]) == 12
    return hx, hz, lx, lz, at, bt


def distance_six_certificate(h, logical):
    """Exhaustive meet-in-the-middle proof: no weight <=5; exhibit weight 6.

    Every support of size <=5 splits into a pair and a triple. Equal syndrome
    and unequal logical parity detects a nontrivial logical XOR of supports.
    """
    n = h.shape[1]
    syndrome = [sum(int(v)<<j for j, v in enumerate(h[:, i])) for i in range(n)]
    parity = [sum(int(v)<<j for j, v in enumerate(logical[:, i])) for i in range(n)]
    low = defaultdict(list)
    for w in range(3):
        for support in itertools.combinations(range(n), w):
            s = l = mask = 0
            for q in support:
                s ^= syndrome[q]; l ^= parity[q]; mask ^= 1<<q
            for previous_l, previous_mask in low[s]:
                assert l == previous_l, ('logical of weight <=4', mask^previous_mask)
            low[s].append((l, mask))
    triples, witness, count = {}, None, 0
    for support in itertools.combinations(range(n), 3):
        s = l = mask = 0
        for q in support:
            s ^= syndrome[q]; l ^= parity[q]; mask ^= 1<<q
        for previous_l, previous_mask in low.get(s, []):
            assert l == previous_l, ('logical of weight <=5', mask^previous_mask)
        if s in triples and triples[s][0] != l:
            wmask = triples[s][1]^mask
            assert wmask.bit_count() == 6
            witness = [q for q in range(n) if wmask>>q & 1]
        else:
            triples[s] = (l, mask)
        count += 1
    assert witness is not None
    return dict(distance=6, triples_checked=count, witness=witness)


@dataclass
class CodeCircuit:
    name: str
    distance: int
    data: list
    rounds: int
    basis: str
    circuit: stim.Circuit
    keep_detectors: list
    logical_qubits: int


def surface_circuit(d, basis, p):
    source = stim.Circuit.generated(f'surface_code:rotated_memory_{basis.lower()}', distance=d,
        rounds=d, after_clifford_depolarization=p,
        before_round_data_depolarization=p,
        before_measure_flip_probability=p, after_reset_flip_probability=p).flattened()
    coords = source.get_final_qubit_coordinates()
    data = sorted(q for q, c in coords.items() if int(c[0])%2 and int(c[1])%2)
    xchecks, out, h_count = set(), stim.Circuit(), 0
    for ins in source:
        if ins.name == 'H':
            if h_count == 0:
                xchecks = {t.value for t in ins.targets_copy()}
            if h_count%2 == 0:
                out.append('TICK', tag='fd_layer')
            h_count += 1
        out.append(ins)
    assert h_count == 2*d and len(data) == d*d
    relevant = xchecks if basis == 'X' else set(coords)-set(data)-xchecks
    positions = {tuple(coords[q][:2]) for q in relevant}
    keep = [i for i, c in out.get_detector_coordinates().items() if tuple(c[:2]) in positions]
    return CodeCircuit(f'rotated_surface_d{d}', d, data, d, basis, out, keep, 1)


def bicycle_circuit(basis, p, rounds=6):
    hx, hz, lx, lz, at, bt = bicycle()
    n, half = hx.shape[1], hx.shape[0]
    data = list(range(n))
    ax, az = list(range(n,n+half)), list(range(n+half,2*n))
    c, measurement_count, previous = stim.Circuit(), 0, {}
    c.append('R' if basis == 'Z' else 'RX', data)
    terms_x = [*at, *bt]
    terms_z = [*(v.T for v in bt), *(v.T for v in at)]
    order_x, order_z = [None,1,4,3,5,0,2], [3,5,0,1,2,4,None]
    # Published CNOT ordering, with both ancilla families reset at the
    # cycle start and read at its end. Added waits have explicit idle noise.
    for t in range(rounds+1):
        noise = 0. if t == 0 else p
        if t:
            c.append('TICK', tag='fd_layer')
        c.append('RX',ax); c.append('R',az)
        if noise:
            c.append('Z_ERROR',ax,noise); c.append('X_ERROR',az,noise)
            c.append('DEPOLARIZE1',data,noise)
        for ox,oz in zip(order_x,order_z):
            pairs = []
            for label,order,terms,anc in [('X',ox,terms_x,ax),('Z',oz,terms_z,az)]:
                if order is None:
                    continue
                neighbours=np.argmax(terms[order],axis=1)+(half if order>=3 else 0)
                pairs += [q for ancilla,qdata in zip(anc,neighbours)
                          for q in ((ancilla,int(qdata)) if label=='X' else (int(qdata),ancilla))]
            assert len(pairs)==len(set(pairs)), 'CNOT layer uses a qubit twice'
            c.append('CX',pairs)
            if noise:
                c.append('DEPOLARIZE2',pairs,noise)
                idle=sorted(set(range(2*n))-set(pairs))
                if idle:
                    c.append('DEPOLARIZE1',idle,noise)
            c.append('TICK')
        if noise:
            c.append('Z_ERROR',ax,noise); c.append('X_ERROR',az,noise)
            c.append('DEPOLARIZE1',data,noise)
        for label,anc in [('X',ax),('Z',az)]:
            c.append('MX' if label=='X' else 'M',anc)
            records=np.arange(measurement_count,measurement_count+half)
            measurement_count+=half
            if label==basis:
                for j,rec in enumerate(records):
                    refs=[stim.target_rec(int(rec-measurement_count))]
                    if label in previous:
                        refs.append(stim.target_rec(int(previous[label][j]-measurement_count)))
                    c.append('DETECTOR',refs,[j,t])
            previous[label]=records
    # Ideal final data readout closes the decoding time boundary.
    c.append('M' if basis == 'Z' else 'MX', data)
    final_start = measurement_count
    measurement_count += n
    checks, logical = (hz, lz) if basis == 'Z' else (hx, lx)
    for j, row in enumerate(checks):
        refs = [stim.target_rec(int(final_start+q-measurement_count)) for q in np.flatnonzero(row)]
        refs.append(stim.target_rec(int(previous[basis][j]-measurement_count)))
        c.append('DETECTOR', refs, [j, rounds+1])
    for j, row in enumerate(logical):
        c.append('OBSERVABLE_INCLUDE', [stim.target_rec(int(final_start+q-measurement_count))
                 for q in np.flatnonzero(row)], j)
    c.detector_error_model(approximate_disjoint_errors=True)
    return CodeCircuit('bicycle_72_12_6', 6, data, rounds, basis, c,
                       list(range(c.num_detectors)), 12)


def insert_instrument(code, channel, mode='sample'):
    """Preserve original measurement references while inserting herald records."""
    out, flags, mapping = stim.Circuit(), [], {}
    old_m = new_m = slot = 0
    h, probs = channel['herald'], np.array(channel['pauli'])
    bit_p = [probs[1]+probs[2], probs[3]+probs[2]]
    if mode == 'hidden':
        bit_p = [(1-h)*p+h/2 for p in bit_p]
    for ins in code.circuit.flattened():
        if ins.tag == 'fd_layer':
            # G=-Z: the irrelevant common sign is omitted; Stim tracks the
            # known ideal logical/syndrome reference of the Z layer.
            out.append('Z', code.data)
            if mode == 'sample':
                out.append('PAULI_CHANNEL_1', code.data, probs[1:])
                out.append('HERALDED_ERASE', code.data, h)
                flags.extend(range(new_m, new_m+len(code.data)))
                new_m += len(code.data)
                slot += len(code.data)
            else:
                for q in code.data:
                    for label, prob in zip(('X', 'Z'), bit_p):
                        # Tiny positive model-only priors preserve all columns
                        # for reweighting, including exact noiseless limits.
                        out.append(label+'_ERROR', [q], max(float(prob), 1e-15), tag=f'fd:{slot}:{label}')
                    slot += 1
            continue
        targets = [stim.target_rec(mapping[old_m+t.value]-new_m) if t.is_measurement_record_target else t
                   for t in ins.targets_copy()]
        out.append(ins.name, targets, ins.gate_args_copy(), tag=ins.tag)
        if ins.name in ('M', 'MX', 'MY', 'MR', 'MRX', 'MRY'):
            for j in range(len(targets)):
                mapping[old_m+j] = new_m+j
            old_m += len(targets); new_m += len(targets)
    assert old_m == code.circuit.num_measurements
    assert slot == len(code.data)*code.rounds
    assert out.num_detectors == code.circuit.num_detectors
    return out, flags


@dataclass
class DecoderModel:
    h: sparse.csc_matrix
    logical: sparse.csc_matrix
    probabilities: np.ndarray
    slots: dict


def decoder_model(circuit, keep, logical_count):
    dem = circuit.detector_error_model(approximate_disjoint_errors=True).flattened()
    rowmap = {original: current for current, original in enumerate(keep)}
    grouped = OrderedDict()
    for ins in dem:
        if ins.type != 'error':
            continue
        detectors = tuple(sorted(rowmap[t.val] for t in ins.targets_copy()
                                 if t.is_relative_detector_id() and t.val in rowmap))
        logical = tuple(sorted(t.val for t in ins.targets_copy() if t.is_logical_observable_id()))
        if not detectors and not logical:
            continue
        key = (detectors, logical)
        entry = grouped.setdefault(key, [0., set()])
        p = ins.args_copy()[0]
        entry[0] = entry[0]+p-2*entry[0]*p
        if ins.tag.startswith('fd:'):
            entry[1].add(int(ins.tag.split(':')[1]))
    hr, hc, lr, lc, probabilities, slots = [], [], [], [], [], defaultdict(set)
    for col, ((ds, ls), (p, locations)) in enumerate(grouped.items()):
        hr.extend(ds); hc.extend([col]*len(ds))
        lr.extend(ls); lc.extend([col]*len(ls))
        probabilities.append(p)
        for slot in locations:
            slots[slot].add(col)
    h = sparse.csc_matrix((np.ones(len(hr), dtype=np.uint8), (hr, hc)), shape=(len(keep), len(grouped)))
    l = sparse.csc_matrix((np.ones(len(lr), dtype=np.uint8), (lr, lc)), shape=(logical_count, len(grouped)))
    return DecoderModel(h, l, np.array(probabilities), dict(slots))


class Decoder:
    def __init__(self, model, kind):
        self.model, self.kind, self.cache = model, kind, OrderedDict()
        self.base_weights = np.log((1-model.probabilities)/model.probabilities)
        if kind == 'matching':
            assert max(np.diff(model.h.indptr), default=0) <= 2
            self.base = self.make_matching(())
        else:
            self.base = BpOsdDecoder(model.h, error_channel=model.probabilities.tolist(),
                max_iter=30, bp_method='minimum_sum', ms_scaling_factor=.625,
                schedule='parallel', osd_method='osd_cs', osd_order=2)

    def make_matching(self, erased):
        weights = self.base_weights.copy()
        weights[list(erased)] = 0.
        return pymatching.Matching(self.model.h, weights=weights, faults_matrix=self.model.logical)

    def erased_columns(self, flags):
        return tuple(sorted(set().union(*(self.model.slots.get(int(s), set()) for s in np.flatnonzero(flags)))))

    def decode(self, syndromes, flags=None):
        k = self.model.logical.shape[0]
        output = np.zeros((len(syndromes), k), dtype=np.uint8)
        if self.kind == 'matching':
            groups = defaultdict(list)
            for i, s in enumerate(syndromes):
                if np.any(s):
                    groups[self.erased_columns(flags[i]) if flags is not None else ()].append(i)
            for erased, indices in groups.items():
                if not erased:
                    decoder = self.base
                elif erased in self.cache:
                    decoder = self.cache[erased]
                    self.cache.move_to_end(erased)
                else:
                    decoder = self.make_matching(erased)
                    self.cache[erased] = decoder
                    if len(self.cache) > 1024:
                        self.cache.popitem(last=False)
                output[indices] = decoder.decode_batch(syndromes[indices])
        else:
            previous = None
            for i, s in enumerate(syndromes):
                if not np.any(s):
                    continue
                erased = self.erased_columns(flags[i]) if flags is not None else ()
                if erased != previous:
                    p = self.model.probabilities.copy()
                    p[list(erased)] = .5
                    self.base.update_channel_probs(p)
                    previous = erased
                estimate = self.base.decode(s)
                assert np.array_equal((self.model.h@estimate)%2, s), 'decoder violated syndrome'
                output[i] = (self.model.logical@estimate)%2
        return output


def interval(failures, shots):
    if failures == 0:
        return [0., float(-np.expm1(np.log(.05)/shots))]
    return [float(beta.ppf(.025, failures, shots-failures+1)),
            float(beta.ppf(.975, failures+1, shots-failures)) if failures < shots else 1.]


def run_case(code, channel, shots, seed, output_dir=None, batch_size=1000):
    actual, flag_indices = insert_instrument(code, channel)
    conditional, _ = insert_instrument(code, channel, 'conditional')
    hidden, _ = insert_instrument(code, channel, 'hidden')
    model = decoder_model(conditional, code.keep_detectors, code.logical_qubits)
    hidden_model = decoder_model(hidden, code.keep_detectors, code.logical_qubits)
    kind = 'matching' if code.name.startswith('rotated') else 'bposd'
    informed, blind = Decoder(model, kind), Decoder(hidden_model, kind)
    sampler = actual.compile_sampler(seed=seed)
    converter = actual.compile_m2d_converter()
    failure_counts = np.zeros(2, dtype=int)
    logical_counts = np.zeros((2, code.logical_qubits), dtype=int)
    discordant = np.zeros(2, dtype=int)
    erasures = blocks_with_flag = 0
    started = time.monotonic()
    for offset in range(0, shots, batch_size):
        count = min(batch_size, shots-offset)
        measurements = sampler.sample(count)
        detectors, obs = converter.convert(measurements=measurements, separate_observables=True)
        flags = measurements[:, flag_indices]
        syndrome = np.ascontiguousarray(detectors[:, code.keep_detectors], dtype=np.uint8)
        pred = [informed.decode(syndrome, flags), blind.decode(syndrome)]
        failures = [np.any(p != obs, axis=1) for p in pred]
        for j in range(2):
            failure_counts[j] += int(np.sum(failures[j]))
            logical_counts[j] += np.sum(pred[j] != obs, axis=0)
        discordant += [np.sum(failures[0] & ~failures[1]), np.sum(failures[1] & ~failures[0])]
        erasures += int(flags.sum())
        blocks_with_flag += int(np.any(flags, axis=1).sum())
        if (offset+count)%50000 == 0 or offset+count == shots:
            print(f'{code.name} {code.basis}: {offset+count}/{shots}; failures visible/hidden '
                  f'{failure_counts.tolist()}; {time.monotonic()-started:.1f}s', flush=True)
    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)
        actual.to_file(output_dir/f'{code.name}_{code.basis}.stim')
    result = dict(code=code.name, distance=code.distance, data_qubits=len(code.data),
        logical_qubits=code.logical_qubits, syndrome_rounds=code.rounds, basis=code.basis,
        shots=shots, seed=seed, decoder=kind, channel=channel,
        detector_count=model.h.shape[0], error_mechanisms=model.h.shape[1],
        circuit_qubits=actual.num_qubits, native_gate_locations=len(flag_indices),
        native_parallel_layer_time=channel.get('total_time', 0)*code.rounds,
        sampled_erasure_events=erasures, sampled_blocks_with_flag=blocks_with_flag,
        paired_discordance_visible_only_hidden_only=discordant.tolist(), elapsed_seconds=time.monotonic()-started)
    for j, label in enumerate(('flags_visible', 'flags_hidden')):
        result[label] = dict(failures=int(failure_counts[j]), rate=float(failure_counts[j]/shots),
            interval_95=interval(int(failure_counts[j]), shots), logical_failures=logical_counts[j].tolist())
    return result


def checks():
    results = {}
    for n in (1, 2, 3):
        error = verify_twirl(operating_point(n=n))
        assert error < 1e-11
        results[f'n{n}_exact_instrument_twirl_error'] = error
    hx, hz, lx, lz, _, _ = bicycle()
    results['bicycle_X_distance'] = distance_six_certificate(hz, lz)
    results['bicycle_Z_distance'] = distance_six_certificate(hx, lx)
    results['bicycle_parameters'] = [72, 72-len(rref(hx)[1])-len(rref(hz)[1]), 6]
    channel = native_channel(operating_point())
    zero = dict(herald=0., pauli=[1.,0.,0.,0.])
    for factory in (lambda b,p: surface_circuit(3,b,p), bicycle_circuit):
        for basis in ('X','Z'):
            code = factory(basis, 0.)
            actual, flags = insert_instrument(code, zero)
            samples = actual.compile_detector_sampler(seed=8).sample(32, append_observables=True)
            assert not np.any(samples)
            results[f'{code.name}_{basis}_noiseless'] = True
            conditional, _ = insert_instrument(code, channel, 'conditional')
            dm = decoder_model(conditional, code.keep_detectors, code.logical_qubits)
            decoder = Decoder(dm, 'matching' if code.name.startswith('rotated') else 'bposd')
            # All individual model faults, including their known logical action.
            syndromes, expected = dm.h.T.toarray(), dm.logical.T.toarray()
            prediction = decoder.decode(syndromes)
            assert np.array_equal(prediction, expected), f'{code.name} {basis} single fault'
            results[f'{code.name}_{basis}_single_native_faults'] = len(syndromes)
            noisy = factory(basis, .001)
            model_circuit, _ = insert_instrument(noisy, channel, 'conditional')
            noisy_dm = decoder_model(model_circuit, noisy.keep_detectors, noisy.logical_qubits)
            noisy_decoder = Decoder(noisy_dm, 'matching' if code.name.startswith('rotated') else 'bposd')
            pred = noisy_decoder.decode(noisy_dm.h.T.toarray())
            assert np.array_equal(pred, noisy_dm.logical.T.toarray()), 'single circuit fault'
            results[f'{code.name}_{basis}_single_circuit_faults'] = noisy_dm.h.shape[1]
            # Erasure reweighting has a separate correction check; probability
            # conservation or an unweighted decoder check would not test it.
            if code.name.startswith('rotated'):
                cases = list(itertools.combinations(range(len(code.data)*code.rounds), 2))
            else:
                rng = np.random.default_rng(921)
                cases = [tuple(rng.choice(len(code.data),5,replace=False)) for _ in range(128)]
            erasure_trials = 0
            for erased_slots in cases:
                f = np.zeros(len(code.data)*code.rounds, dtype=bool)
                f[list(erased_slots)] = True
                columns = decoder.erased_columns(f)
                patterns = np.array(list(itertools.product((0,1),repeat=len(columns))),dtype=np.uint8)
                sy = np.asarray(patterns@dm.h[:,list(columns)].T.toarray()%2,dtype=np.uint8)
                expected = patterns@dm.logical[:,list(columns)].T.toarray()%2
                pred = decoder.decode(sy,np.broadcast_to(f,(len(sy),len(f))))
                assert np.array_equal(pred,expected), 'known erasure recovery'
                erasure_trials += len(sy)
            results[f'{code.name}_{basis}_known_erasure_patterns_corrected'] = erasure_trials
            # A maximally erased code must lose the initially prepared logical
            # information even though every location has a known flag.
            erased = dict(herald=1., pauli=[1.,0.,0.,0.])
            probe, _ = insert_instrument(code, erased)
            det, obs = probe.compile_detector_sampler(seed=11).sample(512, separate_observables=True)
            assert .3 < obs[:,0].mean() < .7
            results[f'{code.name}_{basis}_complete_erasure_randomizes_logical'] = float(obs[:,0].mean())
    return results


def plot(report, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    colors = {'baseline':'#666666','five':'#147d72','background':'#a26736'}
    for ax, basis in zip(axes, ('X','Z')):
        for scenario in sorted({r['scenario'] for r in report['runs']}):
            rows = sorted((r for r in report['runs'] if r['basis']==basis and r['code'].startswith('rotated')
                           and r['scenario']==scenario), key=lambda r:r['distance'])
            if not rows:
                continue
            for label, style in [('flags_visible','-'),('flags_hidden','--')]:
                if scenario=='baseline' and label=='flags_hidden':
                    continue
                x = [r['distance'] for r in rows]
                y = [r[label]['rate'] if r[label]['failures'] else r[label]['interval_95'][1] for r in rows]
                ax.plot(x,y,style,color=colors.get(scenario), label=f'{scenario}, '+('flags used' if style=='-' else 'flags hidden'))
                for r,xx,yy in zip(rows,x,y):
                    lo,hi=r[label]['interval_95']
                    if r[label]['failures']:
                        ax.errorbar(xx,yy,yerr=[[yy-lo],[hi-yy]],fmt='o',color=colors.get(scenario),capsize=3)
                    else:
                        ax.scatter(xx,yy,marker='v',color=colors.get(scenario))
        ax.set(xlabel='Code distance d (d noisy rounds)',ylabel='Logical failure probability per block',
               title=f'Logical {basis} memory',yscale='log',xticks=[3,5,7,9])
        ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.suptitle(f'Finite-dump gate in repeated circuit QEC; syndrome noise p={report["stack_noise"]:g}\n'
                 'Boundary replacement and independent Pauli randomization; triangles: 95% upper bounds')
    fig.savefig(path,dpi=180)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--shots', type=int, default=20000)
    parser.add_argument('--bb-shots', type=int, default=2000)
    parser.add_argument('--distances', type=int, nargs='+', default=[3,5,7])
    parser.add_argument('--noise', type=float, default=.001)
    parser.add_argument('--scenarios', nargs='+', default=['baseline','five'])
    parser.add_argument('--basis', nargs='+', default=['X','Z'])
    parser.add_argument('--seed', type=int, default=20260906)
    parser.add_argument('--json', type=Path, default=ROOT/'docs/qec-distance-stack-checks.json')
    parser.add_argument('--plot', type=Path)
    parser.add_argument('--circuits', type=Path)
    args=parser.parse_args()
    import importlib.metadata
    report=dict(scope='Circuit QEC of the repaired, independently Pauli-randomized finite-dump instrument; '
                'native gate on every data qubit before each noisy extraction round; independent reservoirs. '
                'Other circuit operations have the separately declared stack noise. No persistent-leakage or hardware-threshold claim.',
                stack_noise=args.noise,versions={n:importlib.metadata.version(n) for n in ('numpy','scipy','stim','pymatching','ldpc')},
                runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                source_sha256=hashlib.sha256((ROOT/'docs/reflection-loop-finite-dump-checks.json').read_bytes()).hexdigest(),runs=[])
    if args.check:
        report['checks']=checks()
        print(json.dumps(report['checks'],indent=2),flush=True)
    for scenario in args.scenarios:
        channel = (dict(herald=0.,pauli=[1.,0.,0.,0.],total_time=124.87350186269815)
                   if scenario=='baseline' else native_channel(operating_point(background=1e-3 if scenario=='background' else 1e-6)))
        codes=[surface_circuit(d,b,args.noise) for d in args.distances for b in args.basis]
        if args.bb_shots:
            codes += [bicycle_circuit(b,args.noise) for b in args.basis]
        for code in codes:
            shots=args.shots if code.name.startswith('rotated') else args.bb_shots
            if not shots:
                continue
            seed=args.seed+len(report['runs'])*1009
            result=run_case(code,channel,shots,seed,args.circuits/scenario if args.circuits else None)
            result['scenario']=scenario
            report['runs'].append(result)
            args.json.write_text(json.dumps(report,indent=2)+'\n')
    args.json.write_text(json.dumps(report,indent=2)+'\n')
    if args.plot and report['runs']:
        plot(report,args.plot)


if __name__=='__main__':
    main()
