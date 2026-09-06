"""Reorganise the seated-root root. Run from the repo root. Uses git mv; rewrites imports and replay lines; writes docs/MOVED.md."""
import re, glob, os, subprocess, sys
from pathlib import Path
ROOT = Path('.').resolve()
def sh(*a): return subprocess.check_output(a, text=True)
def gitmv(src, dst):
    if not Path(src).exists():
        assert Path(dst).exists(), f'{src} missing and {dst} absent'; return          # already moved (idempotent rerun)
    Path(dst).parent.mkdir(parents=True, exist_ok=True)
    tracked = subprocess.run(['git', 'ls-files', '--error-unmatch', src], capture_output=True).returncode == 0
    if tracked: sh('git', 'mv', src, dst)
    else: Path(src).rename(dst); sh('git', 'add', dst)                                 # untracked: move and stage

QC_LIB = ['reflection_loop_dynamics', 'reflection_loop_composite', 'reflection_loop_reference_echo', 'reflection_loop_finite_dump',
          'reflection_loop_detuning_order', 'reflection_loop_detuning_certificate', 'reflection_loop_compression', 'reflection_loop_short_correction',
          'qec_distance_stack', 'rh1_common', 'rh2_precision', 'rh4_capture']
EXP = {'2026-09-05': ['split1_split_loop', 'qi_error_harness', 'claire_v0_thirteen'],
       '2026-09-06': ['rh1a', 'rh1b', 'rh1c', 'rh1d', 'rh1_cella_surface', 'rh1_gimbal_probe', 'rh1_passive_reference', 'rh2', 'rh2_a0_closed_form',
                      'rh2p_five', 'rh3_qec_accumulator', 'rh4_readings', 'rh_figures']}
TESTS = ['test_rh1', 'test_rh1_passive_reference', 'test_qec_distance_stack', 'test_thm_n_kerr_quadrupole']
KEEP = {'figures', 'build_html'}
allpy = sorted(p[:-3] for p in glob.glob('*.py'))
SUITES = [m for m in allpy if m not in QC_LIB and m not in sum(EXP.values(), []) and m not in TESTS and m not in KEEP]
moved = {}
for m in QC_LIB: gitmv(f'{m}.py', f'rlq/{m}.py'); moved[m] = f'rlq/{m}.py'
for d, ms in EXP.items():
    for m in ms: gitmv(f'{m}.py', f'experiments/{d}/{m}.py'); moved[m] = f'experiments/{d}/{m}.py'
for m in TESTS: gitmv(f'{m}.py', f'tests/{m}.py'); moved[m] = f'tests/{m}.py'
for m in SUITES: gitmv(f'{m}.py', f'suites/{m}.py'); moved[m] = f'suites/{m}.py'
# logs and receipts sitting in the root
for f in glob.glob('*.log') + glob.glob('*_results.json'):
    stem = f.rsplit('.', 1)[0]
    if stem in ('qi_error_harness', 'claire_v0_thirteen'): gitmv(f, f'experiments/2026-09-05/{f}')
    else: gitmv(f, f'suites/logs/{f}')
if Path('qec_stack_requirements.txt').exists(): gitmv('qec_stack_requirements.txt', 'rlq/requirements-qec.txt')
Path('tests/__init__.py').write_text(''); sh('git', 'add', 'tests/__init__.py')

BOOT = "import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[{depth}]))  # repo root on sys.path (reorg 2026-09-06)\n"
def rewrite_imports(path, libs, depth=None, sibling_dir=False):
    s = Path(path).read_text()
    for m in libs:
        s = re.sub(rf'^(\s*)from {m} import', rf'\1from rlq.{m} import', s, flags=re.M)
        s = re.sub(rf'^(\s*)import {m}\b(?! as)', rf'\1import rlq.{m} as {m}', s, flags=re.M)
        s = re.sub(rf'^(\s*)import {m} as (\w+)', rf'\1import rlq.{m} as \2', s, flags=re.M)
    if depth is not None and 'repo root on sys.path' not in s:
        # insert the bootstrap after the module docstring / shebang / __future__ block
        lines = s.splitlines(keepends=True); i = 0
        if lines and lines[0].startswith('#!'): i = 1
        if i < len(lines) and lines[i].lstrip().startswith(('"""', "'''")):
            q = lines[i].lstrip()[:3]
            if lines[i].count(q) >= 2 and len(lines[i].strip()) > 3: i += 1
            else:
                i += 1
                while i < len(lines) and q not in lines[i]: i += 1
                i += 1
        while i < len(lines) and lines[i].startswith('from __future__'): i += 1
        boot = BOOT.format(depth=depth)
        if sibling_dir: boot += "_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))\n"
        lines.insert(i, boot); s = ''.join(lines)
    Path(path).write_text(s)

for m in QC_LIB: rewrite_imports(f'rlq/{m}.py', QC_LIB, depth=1)
for f in ['rlq/instruments.py', 'rlq/channels.py', 'rlq/words.py', 'rlq/decoder.py', 'rlq/precision.py']: rewrite_imports(f, QC_LIB)
for d, ms in EXP.items():
    for m in ms: rewrite_imports(f'experiments/{d}/{m}.py', QC_LIB, depth=2)
for m in TESTS: rewrite_imports(f'tests/{m}.py', QC_LIB, depth=1, sibling_dir=True)
for m in ['test_thm_n_kerr_quadrupole']:   # imports a suite sibling
    s = Path(f'tests/{m}.py').read_text().replace("_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))", "_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]/'suites'))"); Path(f'tests/{m}.py').write_text(s)
for m in ['cella_normal_geometry']: rewrite_imports(f'suites/{m}.py', [], depth=1, sibling_dir=True)

# path strings inside the moved code
def sub(path, pairs):
    s = Path(path).read_text()
    for a, b in pairs: s = s.replace(a, b)
    Path(path).write_text(s)
sub('rlq/qec_distance_stack.py', [("ROOT = Path(__file__).resolve().parent", "ROOT = Path(__file__).resolve().parents[1]")])
sub('rlq/rh1_common.py', [("files = [f'{test}.py','rh1_common.py','split1_split_loop.py',\n             'reflection_loop_finite_dump.py','reflection_loop_dynamics.py',",
                           "files = [f'experiments/2026-09-06/{test}.py','rlq/rh1_common.py','experiments/2026-09-05/split1_split_loop.py',\n             'rlq/reflection_loop_finite_dump.py','rlq/reflection_loop_dynamics.py',")])
sub('rlq/rh2_precision.py', [("files = ['rh2_precision.py', 'rh1_common.py', 'reflection_loop_dynamics.py', 'reflection_loop_composite.py', 'reflection_loop_reference_echo.py']",
                              "files = ['rlq/rh2_precision.py', 'rlq/rh1_common.py', 'rlq/reflection_loop_dynamics.py', 'rlq/reflection_loop_composite.py', 'rlq/reflection_loop_reference_echo.py']")])
sub('experiments/2026-09-06/rh2p_five.py', [("files = ['rh2p_five.py', 'rh2_precision.py', 'rh1_common.py', 'reflection_loop_dynamics.py', 'reflection_loop_composite.py']",
                                              "files = ['experiments/2026-09-06/rh2p_five.py', 'rlq/rh2_precision.py', 'rlq/rh1_common.py', 'rlq/reflection_loop_dynamics.py', 'rlq/reflection_loop_composite.py']")])
sub('experiments/2026-09-06/rh2_a0_closed_form.py', [("files = ['rh2_a0_closed_form.py', 'rh2_precision.py', 'rh1_common.py', 'reflection_loop_dynamics.py']",
                                                       "files = ['experiments/2026-09-06/rh2_a0_closed_form.py', 'rlq/rh2_precision.py', 'rlq/rh1_common.py', 'rlq/reflection_loop_dynamics.py']")])
sub('experiments/2026-09-06/rh3_qec_accumulator.py', [("('rh3_qec_accumulator.py', 'qec_distance_stack.py', 'docs/reflection-loop-finite-dump-checks.json')",
                                                        "('experiments/2026-09-06/rh3_qec_accumulator.py', 'rlq/qec_distance_stack.py', 'docs/reflection-loop-finite-dump-checks.json')")])
sub('experiments/2026-09-06/rh4_readings.py', [("ROOT = Path(__file__).resolve().parent", "ROOT = Path(__file__).resolve().parents[2]"),
                                                ("('rh4_readings.py', 'rh4_capture.py', 'reflection_loop_finite_dump.py')", "('experiments/2026-09-06/rh4_readings.py', 'rlq/rh4_capture.py', 'rlq/reflection_loop_finite_dump.py')")])
sub('experiments/2026-09-06/rh_figures.py', [("ROOT = Path(__file__).resolve().parent/'docs'", "ROOT = Path(__file__).resolve().parents[2]/'docs'")])
sub('rlq/rh4_capture.py', [("('rh4_capture.py', 'reflection_loop_finite_dump.py', 'qec_distance_stack.py')", "('rlq/rh4_capture.py', 'rlq/reflection_loop_finite_dump.py', 'rlq/qec_distance_stack.py')")])
sub('experiments/2026-09-06/rh4_capture.py', [("modules=[ROOT/'reflection_loop_finite_dump.py', ROOT/'qec_distance_stack.py']", "modules=[ROOT/'rlq/reflection_loop_finite_dump.py', ROOT/'rlq/qec_distance_stack.py']")])

# replay lines in docs/*.md
for md in glob.glob('docs/*.md'):
    s = Path(md).read_text(); t = s
    for m, dst in moved.items():
        t = re.sub(rf'(?<![\w/]){m}\.py(?!\w)', dst, t)
    t = t.replace('python3 -m unittest test_rh1 -v', 'python3 -m unittest tests.test_rh1 -v')
    if t != s: Path(md).write_text(t)
lines = ['# Root reorganisation, 2026-09-06\n', 'Every file that moved, and where. Receipts committed before this date hash the old paths; git history has them.\n',
         'Run everything from the repository root.\n', '| from | to |\n|---|---|\n'] + [f'| {m}.py | {dst} |\n' for m, dst in moved.items()] + \
        ['| *.log, curv1_path1_results.json | suites/logs/ (qi_error_harness.log, claire_v0_thirteen.log -> experiments/2026-09-05/) |\n',
         '| qec_stack_requirements.txt | rlq/requirements-qec.txt |\n']
Path('docs/MOVED.md').write_text(''.join(lines)); sh('git', 'add', 'docs/MOVED.md')
print(f'moved {len(moved)} modules; root .py now:', sorted(p for p in glob.glob("*.py")))
