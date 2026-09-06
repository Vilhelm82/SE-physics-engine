"""Independent forced-event checks of the circuit/decoder herald interface."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))  # repo root on sys.path (reorg 2026-09-06)
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
import itertools
import unittest

import numpy as np
import stim

from rlq.qec_distance_stack import (bicycle, bicycle_circuit, decoder_model, insert_instrument,
                                native_channel, operating_point, surface_circuit)


def force_one_erasure(circuit, selected_slot):
    result, slot = stim.Circuit(), 0
    for instruction in circuit:
        if instruction.name != 'HERALDED_ERASE':
            result.append(instruction)
            continue
        for target in instruction.targets_copy():
            result.append('HERALDED_ERASE', [target], float(slot == selected_slot))
            slot += 1
    return result


class HeraldLocations(unittest.TestCase):
    def test_parallel_bicycle_extraction_equals_sequential_on_all_states(self):
        hx, hz, _, _, _, _ = bicycle()
        n, m = hx.shape[1], hx.shape[0]
        parallel, sequential = stim.Circuit(), stim.Circuit()
        for instruction in bicycle_circuit('X', 0.).circuit:
            if instruction.name in ('M', 'MX'):
                break
            if instruction.name == 'CX':
                parallel.append(instruction)
        for j, row in enumerate(hx):
            for q in np.flatnonzero(row):
                sequential.append('CX', [n+j, int(q)])
        for j, row in enumerate(hz):
            for q in np.flatnonzero(row):
                sequential.append('CX', [int(q), n+m+j])
        self.assertEqual(stim.Tableau.from_circuit(parallel), stim.Tableau.from_circuit(sequential))

    def test_every_native_erasure_location_against_actual_circuit(self):
        native = native_channel(operating_point())
        zero = dict(herald=0., pauli=[1., 0., 0., 0.])
        count = 0
        for factory in (lambda b: surface_circuit(3, b, 0.), lambda b: bicycle_circuit(b, 0.)):
            for basis in ('X', 'Z'):
                code = factory(basis)
                actual, flag_records = insert_instrument(code, zero)
                surrogate, _ = insert_instrument(code, native, 'conditional')
                model = decoder_model(surrogate, code.keep_detectors, code.logical_qubits)
                converter = actual.compile_m2d_converter()
                for selected in range(len(flag_records)):
                    forced = force_one_erasure(actual, selected)
                    records = forced.compile_sampler(seed=500+selected).sample(16)
                    flags = records[:, flag_records]
                    expected_flag = np.zeros(len(flag_records), dtype=bool)
                    expected_flag[selected] = True
                    np.testing.assert_array_equal(flags, np.broadcast_to(expected_flag, flags.shape))
                    det, obs = converter.convert(measurements=records, separate_observables=True)
                    columns = sorted(model.slots.get(selected, ()))
                    patterns = np.array(list(itertools.product((0, 1), repeat=len(columns))), dtype=np.uint8)
                    h = model.h[:, columns].T.toarray()
                    l = model.logical[:, columns].T.toarray()
                    allowed = {row.tobytes() for row in np.hstack([patterns@h%2, patterns@l%2]).astype(np.uint8)}
                    observed = np.hstack([det[:, code.keep_detectors], obs]).astype(np.uint8)
                    self.assertTrue(all(row.tobytes() in allowed for row in observed),
                                    (code.name, basis, selected))
                    count += 1
        self.assertEqual(count, 918)
        print(f'Forced herald checks: {count} locations, {16*count} circuit shots', flush=True)


if __name__ == '__main__':
    unittest.main()
