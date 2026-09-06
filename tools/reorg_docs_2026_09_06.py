"""Reorganise docs/. Run from the repo root. git mv + one path map applied to code, root md, and docs md."""
import re, glob, os, subprocess
from pathlib import Path
def sh(*a): return subprocess.check_output(a, text=True)
def mv(src, dst):
    if not Path(src).exists(): return
    Path(dst).parent.mkdir(parents=True, exist_ok=True)
    tracked = subprocess.run(['git', 'ls-files', '--error-unmatch', src], capture_output=True).returncode == 0
    (sh('git', 'mv', src, dst) if tracked else (Path(src).rename(dst), sh('git', 'add', dst)))

D = Path('docs'); moves = {}   # old relative-to-docs -> new relative-to-docs
STANDING = {'CONJECTURE-COSMOLOGY.md', 'EXTERNAL-MATHEMATICS-DEBT.md', 'LABELLED-MODEL.md', 'TESTING-SCHEDULE.md',
            'mathematical-foundation-audit.md', 'lorentzian-foundation-upgrade-programme-2026-09-05.md'}
for p in sorted(D.iterdir()):
    n = p.name
    if p.is_dir():
        if n == 'prereg': continue
        if n == 'seated_root_figures': moves[n] = f'figures/{n}'
        elif n == 'qec-distance-circuits': moves[n] = f'receipts/{n}'
        else: moves[n] = f'archive/{n}'
    elif n in STANDING or n == 'README.md': continue
    elif re.match(r'20\d\d-\d\d-\d\d', n): moves[n] = f'results/{n[:10]}/{n}'
    elif n.endswith(('.json', '.log')): moves[n] = f'receipts/{n}'
    elif n.endswith('.png'): moves[n] = f'figures/{n}'
    else: moves[n] = f'archive/{n}'
for old, new in moves.items(): mv(f'docs/{old}', f'docs/{new}')

# ---- path rewriting
def new_docs_path(old): return f'docs/{moves.get(old, old)}'
def rewrite_code_and_root(text):
    # any 'docs/<old>' reference (absolute, repo-relative, in code or prose) -> 'docs/<new>'
    def rep(m):
        rest = m.group(1); first = rest.split('/')[0]
        if first in moves and not (first == 'prereg'):
            tail = rest[len(first):]
            return 'docs/' + moves[first] + tail
        return m.group(0)
    text = text.replace('docs/chat snippets', 'docs/archive/chat snippets')            # the one name with a space
    return re.sub(r'docs/([A-Za-z0-9_.\-]+(?:/[A-Za-z0-9_.\-]+)*)', rep, text)
def rewrite_docs_md(text, new_rel):
    """Links inside a moved md: bare '(<old>)' or '(./<old>)' sibling references, and '(../<rootfile>)' repo-root links."""
    depth = new_rel.count('/')                     # results/2026-09-06/x.md -> 2
    up = '../'*depth
    text = rewrite_code_and_root(text)
    def sib(m):
        name = m.group(2)
        if name in moves: return m.group(1) + up + moves[name] + m.group(3)
        return m.group(0)
    text = re.sub(r'(\]\()(?:\./)?([A-Za-z0-9_.\-]+\.(?:md|json|log|png|py))(\))', sib, text)
    text = re.sub(r'(\]\()\.\./([A-Za-z0-9_.\-/]+)(\))', lambda m: m.group(1) + up + '../' + m.group(2) + m.group(3), text)
    return text
for f in glob.glob('rlq/*.py') + glob.glob('experiments/*/*.py') + glob.glob('tests/*.py') + glob.glob('*.md') + ['docs/MOVED.md']:
    p = Path(f); s = p.read_text(errors='replace'); t = rewrite_code_and_root(s)
    if t != s: p.write_text(t)
for f in glob.glob('docs/**/*.md', recursive=True):
    p = Path(f); rel = str(p.relative_to(D)); s = p.read_text(errors='replace')
    t = rewrite_docs_md(s, rel) if '/' in rel else rewrite_code_and_root(s)
    if t != s: p.write_text(t)
# ---- index
lines = ['# docs/\n', '\n', 'Standing documents live here; everything dated or generated has a home by kind.\n', '\n',
         '| where | what |\n|---|---|\n',
         '| `LABELLED-MODEL.md`, `CONJECTURE-COSMOLOGY.md`, `TESTING-SCHEDULE.md`, `EXTERNAL-MATHEMATICS-DEBT.md`, foundation audit/programme | standing state documents |\n',
         '| `results/<date>/` | dated specs, reports and results, by the date in the filename |\n',
         '| `receipts/` | every `*-checks.json`, run log, and the QEC circuit files; the machine-verified side of each result |\n',
         '| `figures/` | figures produced by runners; `figures/seated_root_figures/` the paper figures |\n',
         '| `prereg/` | sealed designer predictions per experiment |\n',
         '| `archive/` | Gram-submersion and inertia-node notes, chat snippets, the missing-window note, plan files |\n',
         '\n', f'Moved 2026-09-06 by `tools/reorg_docs_2026_09_06.py`; {len(moves)} entries. Full map:\n\n| from | to |\n|---|---|\n']
lines += [f'| docs/{o} | docs/{n} |\n' for o, n in moves.items()]
(D/'README.md').write_text(''.join(lines)); sh('git', 'add', 'docs/README.md')
print(f'moved {len(moves)} entries; docs root now:', sorted(p.name for p in D.iterdir()))
