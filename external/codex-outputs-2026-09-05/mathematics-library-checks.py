"""Read-only probes of the live Cella categorical implementation.

Run with python -B; results are printed, source files are never modified.
"""
import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path

sys.dont_write_bytecode = True
ROOT = Path('/home/williaml/Cella Framework')
sys.path.insert(0, str(ROOT / 'engine/src'))
from cella.continuation.r3_ac_fold import (
    ACFoldObject, hostile_ac_fold_parameters, ac_fold_route,
    reduce_ac_fold_morphism, compose_ac_fold_routes,
)
from cella.continuation.selected_skeleton import (
    identity_skeleton_morphism, compose_skeleton_morphisms as compose,
)
from cella.continuation.cce5 import corridor_word_endpoints, CCE5Error

x = ACFoldObject.build(hostile_ac_fold_parameters(), 'high')
r0 = reduce_ac_fold_morphism(ac_fold_route(x, 0))
r2 = reduce_ac_fold_morphism(ac_fold_route(x, 2))
a = ac_fold_route(x, 1)
b = ac_fold_route(a.target, 1)
checks = {
    'reduction_0_equals_2': r0 == r2,
    'left_identity': compose(identity_skeleton_morphism(r0.source), r0) == r0,
    'reduction_preserves_composition': (
        reduce_ac_fold_morphism(compose_ac_fold_routes(a, b))
        == compose(reduce_ac_fold_morphism(a), reduce_ac_fold_morphism(b))
    ),
    'associativity': compose(compose(r0, r0), r0) == compose(r0, compose(r0, r0)),
}
assert all(value is False for value in checks.values())
try:
    corridor_word_endpoints('')
except CCE5Error:
    empty_word_rejected = True
else:
    empty_word_rejected = False
assert empty_word_rejected

# Fine spacing 1, coarse spacing 2, ties to even grid index.
x_round = Fraction(3, 4)
fine = round(x_round)
coarse_after_fine = 2 * round(Fraction(fine, 2))
coarse_direct = 2 * round(x_round / 2)
assert fine == 1 and coarse_after_fine == coarse_direct == 0

sources = [
    ROOT / 'engine/src/cella/continuation/selected_skeleton.py',
    ROOT / 'engine/src/cella/continuation/r3_ac_fold.py',
    ROOT / 'engine/src/cella/continuation/cce5.py',
]
print(json.dumps({
    'categorical_law_results': checks,
    'cce5_empty_word_rejected': empty_word_rejected,
    'rounding_midpoint_preimage_without_failure': {
        'x': str(x_round), 'fine': fine,
        'coarse_after_fine': coarse_after_fine, 'coarse_direct': coarse_direct,
    },
    'interpretation': 'Four categorical failures reproduced; CCE5 identity encoding absent; rounding containment is strict in this example.',
    'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
}, indent=2))
