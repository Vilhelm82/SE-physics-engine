"""Mechanical controls for RH-1; scientific forks are reported, not asserted."""
import importlib.util
import unittest

import numpy as np


class InstrumentControls(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(importlib.util.find_spec('rh1_common'),
                             'RH-1 instrument implementation is missing')
        import rh1_common
        self.m = rh1_common

    def test_isotropic_click_preserves_a_qubit(self):
        m = self.m
        s = np.sqrt(.8)*m.G
        e = np.sqrt(.2)*np.array([[0, 1j], [1, 0]])
        c = m.foldback(s, e)
        self.assertLess(m.db(e)['DB'], 1e-14)
        self.assertLess(m.channel_metrics([s, c@e])['infidelity'], 1e-28)
        self.assertLess(m.channel_metrics([s, c@e])['tp_defect'], 1e-14)

    def test_rank_one_click_cannot_restore_superposition(self):
        m = self.m
        s = np.diag([np.sqrt(.8), 1.])
        e = np.diag([np.sqrt(.2), 0.])
        c = m.foldback(s, e)
        result = m.channel_metrics([s, c@e], target=np.eye(2))
        self.assertEqual(m.db(e)['DB'], 1.)
        self.assertAlmostEqual(result['infidelity'], (1-np.sqrt(.8))/3, places=14)
        self.assertLess(result['tp_defect'], 1e-14)

    def test_conditional_metric_is_stable_below_machine_epsilon(self):
        m = self.m
        theta = 1e-9
        k = np.sqrt(.8)*m.G@np.diag([np.exp(1j*theta), np.exp(-1j*theta)])
        self.assertAlmostEqual(m.conditional_metrics(k)['conditional']/theta**2,
                               2/3, places=10)

    def test_replay_has_a_complete_second_leak_branch(self):
        m = self.m
        s = np.diag([np.sqrt(.8), 1.])
        e = np.diag([np.sqrt(.2), 0.])
        ks = m.replay_branches(s, e, m.polar(e)[0])
        self.assertLess(m.channel_metrics(ks, target=np.eye(2))['tp_defect'], 1e-14)

    def test_trine_matches_original_laboratory_harness(self):
        from split1_split_loop import loop_prop
        m=self.m
        alpha=np.pi/3
        old=loop_prop(m.TAU,gain=1.01,delta=1e-4,alpha=alpha,exact_return=True)
        new=m.trine_loop(alpha,.01,1e-4)
        self.assertLess(np.linalg.norm(old-new),3e-12)

    def test_five_gain_formula_matches_chronological_inverse_integration(self):
        from reflection_loop_finite_dump import loop_endpoints_with_effects
        m=self.m
        exact=m.word_loops('five',[(.01,0.,0.)])[0]
        direct,_=loop_endpoints_with_effects([(.01,0.,0.)],[0.])
        self.assertLess(np.linalg.norm(exact-direct[0]),3e-12)

    def test_per_dump_column_resolves_rank_and_zero_leakage(self):
        from reflection_loop_finite_dump import loop_data_blindness
        m=self.m
        self.assertEqual(loop_data_blindness([m.I4]),[None])
        ls=m.word_loops('five',[(.01,0.,0.)])[0]
        self.assertEqual(loop_data_blindness(ls),[1.]*5)


if __name__ == '__main__':
    unittest.main()
