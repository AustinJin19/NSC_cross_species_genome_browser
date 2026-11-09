import json
import math
import os
from pathlib import Path
import textwrap

os.environ.setdefault('MPLCONFIGDIR', '/tmp/regulatory-enrichment-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

ROOT = Path(__file__).resolve().parent
d = json.loads((ROOT / 'report_data.json').read_text())
groups = ['Shared', 'Corrected only', 'Uncorrected only']
lookup = {(r['Group'], r['Term_ID']): r for r in d['comparison']}
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})
fig = plt.figure(figsize=(13, 10), facecolor='white')
gs = fig.add_gridspec(2, 1, left=.455, right=.79, top=.89, bottom=.15,
                      height_ratios=[10, 2], hspace=.16)
cmap = plt.get_cmap('YlOrRd')
norm = Normalize(0, 3)
area = lambda n: 20 + 23 * n
for ai, library in enumerate(['GO:BP', 'KEGG']):
    ax = fig.add_subplot(gs[ai])
    ids = [t for lib, t in d['selected_terms'] if lib == library]
    labels = [next(r['Term'] for r in d['comparison'] if r['Term_ID'] == t) for t in ids]
    ax.set_xlim(-.5, 2.5)
    ax.set_ylim(len(ids)-.5, -.5)
    ax.set_yticks(range(len(ids)), ['\n'.join(textwrap.wrap(t, 49)) for t in labels])
    ax.set_xticks(range(3), ['Shared\n143 genes', 'Corrected only\n71 / 73 mapped', 'Uncorrected only\n114 genes'])
    ax.tick_params(axis='both', length=0, pad=9)
    ax.tick_params(axis='x', labelbottom=ai == 1)
    ax.grid(axis='both', color='#E7EBEE', linewidth=.6, zorder=0)
    for spine in ax.spines.values():
        spine.set_color('#28323B')
        spine.set_linewidth(.8)
    ax.set_title('GO biological processes' if library == 'GO:BP' else 'KEGG pathways',
                 loc='left', fontsize=11, fontweight='bold', pad=10)
    for y, tid in enumerate(ids):
        for x, group in enumerate(groups):
            r = lookup[group, tid]
            q, k = r['BH_FDR'], r['Overlap']
            if q is None:
                ax.text(x, y, '—', ha='center', va='center', color='#AFB7BD', fontsize=10)
            elif q < .05:
                ax.scatter(x, y, s=area(k), c=[cmap(norm(-math.log10(q)))],
                           edgecolors='#5D362C', linewidths=.5, zorder=3)
            else:
                ax.scatter(x, y, s=area(k), facecolors='white', edgecolors='#8D969D',
                           linewidths=1, zorder=3)

cax = fig.add_axes([.85, .65, .018, .20])
cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
cb.set_label('−log₁₀(BH FDR)', fontsize=10, labelpad=8)
cb.outline.set_linewidth(.6)
legax = fig.add_axes([.825, .31, .17, .29])
legax.axis('off')
handles = [legax.scatter([], [], s=area(n), facecolors='#C9CDD1', edgecolors='#68737D', linewidths=.7) for n in [2, 5, 9]]
legend = legax.legend(handles, ['2 genes', '5 genes', '9 genes'], title='Overlap count',
                      loc='upper left', frameon=False, labelspacing=1.2, borderpad=0, handletextpad=1)
legax.add_artist(legend)
h = [legax.scatter([], [], s=70, color=cmap(norm(2.5)), edgecolors='#5D362C', linewidths=.5),
     legax.scatter([], [], s=70, facecolors='white', edgecolors='#8D969D', linewidths=1)]
legax.legend(h, ['FDR < 0.05', 'FDR ≥ 0.05'], loc='lower left', frameon=False, borderpad=0, labelspacing=1)
fig.text(.03, .965, 'Pathway enrichment by body-mass correction status', fontsize=17, weight='bold', color='#173B55')
fig.text(.03, .935, 'Regulatory-activity hits pooled across TSS ±20–1000 kb windows', fontsize=11, color='#536977')
fig.text(.03, .052, 'Same background: 19,196 mapped tested genes. BH FDR is calculated separately for each group and annotation library.', fontsize=9, color='#536977')
fig.text(.03, .029, 'Representative terms shown; all 27 significant term–group results are in the workbook. A dash indicates no returned overlap term.', fontsize=9, color='#536977')
for ext in ['png', 'pdf', 'svg']:
    fig.savefig(ROOT / f'pathway_enrichment_comparison.{ext}', dpi=180, facecolor='white')
plt.close(fig)
print('Saved comparison plot (PNG, PDF, SVG).')
