"""
Generate a simple static HTML dashboard to browse results and plots under `micromol_results`.
"""
import os
from pathlib import Path

ROOT = Path('micromol_results')
OUT = ROOT / 'dashboard'
OUT.mkdir(parents=True, exist_ok=True)

html_lines = [
    '<!doctype html>',
    '<html>',
    '<head>',
    '<meta charset="utf-8"/>',
    '<title>MicroMol Results Dashboard</title>',
    '<style>body{font-family:Arial,sans-serif;padding:20px} .folder{margin-bottom:24px} img{max-width:400px;margin:8px 0;border:1px solid #ccc}</style>',
    '</head>',
    '<body>',
    '<h1>MicroMol Results Dashboard</h1>',
    '<p>Browse generated result folders and plots.</p>',
]

for folder in sorted(ROOT.iterdir()):
    if not folder.is_dir():
        continue
    # skip dashboard output itself
    if folder.name == 'dashboard':
        continue
    html_lines.append(f'<div class="folder"><h2>{folder.name}</h2>')
    # list csv and png files
    files = sorted(folder.rglob('*.png')) + sorted(folder.rglob('*.csv')) + sorted(folder.rglob('*.json'))
    if not files:
        html_lines.append('<p>No artifacts found.</p></div>')
        continue
    for f in files:
        rel = f.relative_to(OUT.parent)
        if f.suffix.lower() in ['.png', '.jpg', '.jpeg', '.gif']:
            html_lines.append(f'<div><strong>{f.name}</strong><br/><a href="{rel}"><img src="{rel}" alt="{f.name}"/></a></div>')
        else:
            html_lines.append(f'<div><a href="{rel}">{f.name}</a></div>')
    html_lines.append('</div>')

html_lines.append('</body></html>')

with open(OUT / 'index.html', 'w', encoding='utf-8') as f:
    f.write('\n'.join(html_lines))

print(f"Dashboard written to {OUT / 'index.html'}")
