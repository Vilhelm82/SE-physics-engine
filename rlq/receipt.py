"""Receipt: the honest-accounting page as an object, with an artifact folder per run.

    r = Receipt('rh4_capture', 'capture model at the absorbing interface', driver=__file__)
    r.check('capture off == FD herald', value, 1e-9)        # must pass; the run exits nonzero otherwise
    r.fork('RH-4 register', 'a', 'ideal register = baseline')  # an OUTCOME; never pass/fail
    r.held_out('occupancy readout as an operation', ...)      # required; close() refuses an empty list
    r.predictions('docs/prereg/RH-4/DESIGNER-PREDICTIONS.md') # copied into the run folder, hash recorded
    r.record('grid', rows)                                    # any JSON-able data
    r.figure(fig, 'readings')                                 # saved as runs/<date>/<experiment>/readings.png
    sys.exit(r.close())                                       # writes receipt.json, run.log, appends runs/INDEX.md
"""
import hashlib, json, shutil, subprocess, sys, time, datetime
from pathlib import Path
import numpy as np
from . import ROOT


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def norm(value):
    if isinstance(value, np.ndarray):
        return dict(real=value.real.tolist(), imag=value.imag.tolist()) if np.iscomplexobj(value) else value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, Path): return str(value)
    raise TypeError(type(value).__name__)


class _Tee:
    def __init__(self, *streams): self.streams = streams
    def write(self, s):
        for st in self.streams: st.write(s)
    def flush(self):
        for st in self.streams: st.flush()


class Receipt:
    def __init__(self, experiment, description, driver=None, date=None, modules=(), runs_root=None):
        self.t0 = time.time()
        self.date = date or datetime.date.today().isoformat()
        base = (runs_root or ROOT/'runs')/self.date
        self.dir = base/experiment; k = 2
        while self.dir.exists(): self.dir = base/f'{experiment}-{k}'; k += 1
        self.dir.mkdir(parents=True)
        self.driver = Path(driver).resolve() if driver else None
        self.modules = [Path(m).resolve() for m in modules]
        self.data = dict(experiment=experiment, description=description, date=self.date, run_dir=str(self.dir.relative_to(ROOT) if self.dir.is_relative_to(ROOT) else self.dir),
                         checks=[], forks=[], held_out=[], records={}, figures=[], predictions=None, corrections=[])
        self._log = open(self.dir/'run.log', 'w'); self._stdout = sys.stdout; sys.stdout = _Tee(self._stdout, self._log)
        print(f"== {experiment}: {description}\n   run dir {self.data['run_dir']}", flush=True)

    # ---- the four kinds of statement
    def check(self, name, value, tol):
        ok = bool(np.isfinite(value) and abs(value) <= tol)
        self.data['checks'].append(dict(name=name, value=float(value), tolerance=float(tol), passed=ok))
        print(f"  [{'ok' if ok else 'FAIL'}] {name}: {value:.3e} (tol {tol:.0e})", flush=True); return ok
    def fork(self, name, outcome, detail=None):
        self.data['forks'].append(dict(name=name, outcome=outcome, detail=detail)); print(f"  fork {name}: ({outcome}) {detail or ''}", flush=True)
    def held_out(self, *items):
        self.data['held_out'] += [str(i) for i in items]
    def correction(self, text):
        self.data['corrections'].append(dict(time=datetime.datetime.now().isoformat(timespec='seconds'), text=text)); print(f"  correction: {text}", flush=True)
    def predictions(self, path):
        src = ROOT/path; dst = self.dir/'predictions.md'; shutil.copy(src, dst)
        self.data['predictions'] = dict(source=str(path), sha256=sha256(src), copied_to=str(dst))
    def record(self, key, value): self.data['records'][key] = value
    def figure(self, fig, name, dpi=170):
        out = self.dir/f'{name}.png'; fig.savefig(out, dpi=dpi); self.data['figures'].append(str(out)); print(f"  figure {out}", flush=True); return out

    # ---- provenance and closing
    def _provenance(self):
        git = lambda *a: subprocess.check_output(['git', *a], text=True, cwd=ROOT).strip()
        files = sorted(set([self.driver] if self.driver else []) | set(self.modules) | set((ROOT/'rlq').glob('*.py')))
        return dict(head=git('rev-parse', 'HEAD'), dirty=bool(git('status', '--porcelain')), runtime_s=round(time.time() - self.t0, 1),
                    python=sys.version.split()[0], sha256={str(f.relative_to(ROOT)): sha256(f) for f in files})
    def close(self):
        if not self.data['held_out']:
            raise RuntimeError('Receipt refused: held_out list is empty. Say what this run does not cover.')
        self.data['provenance'] = self._provenance()
        c = self.data['checks']; self.data['check_count'] = len(c); self.data['checks_passed'] = sum(x['passed'] for x in c)
        (self.dir/'receipt.json').write_text(json.dumps(self.data, default=norm, indent=2, allow_nan=True) + '\n')
        failed = [x['name'] for x in c if not x['passed']]
        print(f"\nRESULT: {self.data['checks_passed']}/{self.data['check_count']} checks passed; {len(self.data['forks'])} forks recorded; "
              f"{len(self.data['held_out'])} held out; {len(self.data['figures'])} figures.  runtime {self.data['provenance']['runtime_s']}s", flush=True)
        if failed: print('FAILED:', *failed, sep='\n  ', flush=True)
        sys.stdout = self._stdout; self._log.close()
        index = ROOT/'runs'/'INDEX.md'
        if not index.exists(): index.write_text('# Runs\n\n| date | experiment | checks | forks | figures | head | description |\n|---|---|---:|---:|---:|---|---|\n')
        with open(index, 'a') as f:
            f.write(f"| {self.date} | [{self.data['experiment']}]({self.data['run_dir']}/receipt.json) | {self.data['checks_passed']}/{self.data['check_count']} | "
                    f"{len(self.data['forks'])} | {len(self.data['figures'])} | {self.data['provenance']['head'][:7]}{'*' if self.data['provenance']['dirty'] else ''} | {self.data['description']} |\n")
        return 1 if failed else 0
