#!/usr/bin/env python3
# =============================================================================
# build_html.py - render the canonical markdown paper to a SELF-CONTAINED html
#
# WHY THIS SHAPE: the .md stays canonical (git-diffable, toolchain-free,
# readable in eighty years). This script is a one-way RENDERER: never edit the
# html, edit the md and re-run. Images are base64-embedded so the output is a
# single portable file; math is rendered by MathJax (CDN - the only external
# dependency, and the md remains the offline fallback).
# Browser print-to-PDF gives the PDF; @media print rules below make it behave.
# =============================================================================
import base64, io, re, pathlib, markdown

SRC = 'PAPER-seated-root-v0.1.md'
OUT = 'PAPER-seated-root.html'
here = pathlib.Path('.')
text = io.open(SRC, encoding='utf-8').read()

# --- 1. protect math from the markdown parser (underscores become <em> etc.)
store = []
def stash(m):
    store.append(m.group(0)); return f'@@MATHTOKEN{len(store)-1}@@'
text = re.sub(r'\$\$.+?\$\$', stash, text, flags=re.S)
text = re.sub(r'(?<!\\)\$(?:[^$`\n]|\\\$)+?\$', stash, text)

# --- 2. markdown -> html
body = markdown.markdown(text, extensions=['tables', 'fenced_code', 'attr_list'])

# --- 3. restore math verbatim
def unstash(m): return store[int(m.group(1))]
body = re.sub(r'@@MATHTOKEN(\d+)@@', unstash, body)

# --- 4. embed images as data URIs -> one portable file
def embed(m):
    alt, src = m.group(1), m.group(2)
    p = here / src
    if not p.exists():
        return m.group(0)
    b64 = base64.b64encode(p.read_bytes()).decode()
    return (f'<figure><img alt="{alt}" src="data:image/png;base64,{b64}">'
            f'</figure>')
body = re.sub(r'<img alt="([^"]*)" src="([^"]+)"\s*/?>', embed, body)
print(f"figures embedded: {body.count('data:image/png;base64')}")

CSS = """
:root{--ink:#1a1a1a;--mid:#4a4a4a;--acc:#1f5c8b;--warn:#a33b20;--rule:#d8d8d4;
      --bg:#fdfdfb;--tint:#f4f6f8}
*{box-sizing:border-box}
html{-webkit-text-size-adjust:100%}
body{background:var(--bg);color:var(--ink);margin:0;
     font:16px/1.62 "Iowan Old Style","Palatino Linotype",Palatino,Georgia,serif;}
main{max-width:47rem;margin:0 auto;padding:4rem 1.5rem 6rem}
h1{font-size:1.95rem;line-height:1.22;margin:0 0 .3rem;letter-spacing:-.01em}
h2{font-size:1.15rem;font-weight:500;color:var(--mid);font-style:italic;
   margin:0 0 2.2rem;line-height:1.4}
h3{font-size:1.16rem;margin:3rem 0 .9rem;padding-bottom:.35rem;
   border-bottom:1px solid var(--rule);letter-spacing:.005em}
p{margin:0 0 1.05rem;hyphens:auto}
strong{font-weight:650}
em{color:var(--mid)}
a{color:var(--acc)}
hr{border:0;border-top:1px solid var(--rule);margin:2.6rem 0}
blockquote{margin:1.4rem 0;padding:.1rem 0 .1rem 1.1rem;
           border-left:3px solid var(--acc);color:var(--mid)}
code{font:.86em/1.5 ui-monospace,"SF Mono",Menlo,Consolas,monospace;
     background:var(--tint);padding:.1em .35em;border-radius:3px}
pre{background:var(--tint);padding:1rem;overflow-x:auto;border-radius:4px;
    border:1px solid var(--rule)}
pre code{background:none;padding:0}
table{border-collapse:collapse;width:100%;margin:1.4rem 0;font-size:.93rem}
th,td{border-bottom:1px solid var(--rule);padding:.5rem .6rem;text-align:left}
th{font-weight:650;border-bottom:2px solid var(--ink)}
figure{margin:2.2rem 0;text-align:center}
figure img{max-width:100%;height:auto;border:1px solid var(--rule);
           border-radius:3px;background:#fff}
.abstract{background:var(--tint);border:1px solid var(--rule);border-radius:4px;
          padding:1.3rem 1.5rem;font-size:.95rem;margin:0 0 2rem}
.abstract h3{margin-top:0;border:0;font-size:1rem;letter-spacing:.06em;
             text-transform:uppercase;color:var(--mid)}
.meta{color:var(--mid);font-size:.9rem;margin-bottom:2.4rem}
mjx-container{overflow-x:auto;overflow-y:hidden}
mjx-container[display="true"]{margin:1.2rem 0!important}
@media print{
  body{background:#fff;font-size:10.6pt}
  main{max-width:none;padding:0}
  h3{break-after:avoid-page}
  figure,table,pre,blockquote{break-inside:avoid}
  .abstract{background:none}
}
@media (max-width:640px){main{padding:2rem 1.1rem 4rem}body{font-size:15px}}
"""

HTML = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>The Seated Root</title>
<style>{CSS}</style>
<script>window.MathJax={{tex:{{inlineMath:[['$','$']],
  displayMath:[['$$','$$']],processEscapes:true}},
  options:{{skipHtmlTags:['script','noscript','style','textarea','pre','code']}}}};
</script>
<script id="MathJax-script" async
  src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
</head><body><main>
{body}
</main></body></html>"""

io.open(OUT, 'w', encoding='utf-8').write(HTML)
kb = pathlib.Path(OUT).stat().st_size / 1024
print(f"wrote {OUT}  ({kb:.0f} KB, self-contained)")
