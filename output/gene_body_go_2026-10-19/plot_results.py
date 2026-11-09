import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

here = Path(__file__).resolve().parent
out = here.parent.parent / 'outputs/go-gene-body-2026-10-19'
data = json.loads((here / 'report_data.json').read_text())
by_id = {r['go_id']: r for r in data['groups']['negative']['terms']}
rows = [by_id[gid] for gid in data['representative_ids']]
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'pdf.fonttype': 42, 'ps.fonttype': 42})
fig = plt.figure(figsize=(11.8, 6.7), facecolor='white')
ax = fig.add_axes([0.36, 0.25, 0.43, 0.53])
y = np.arange(len(rows))[::-1]
points = ax.scatter([r['fold_enrichment'] for r in rows], y,
                    s=[r['overlap'] * 1.2 for r in rows],
                    c=[-np.log10(r['fdr']) for r in rows], cmap='viridis',
                    vmin=0, vmax=11, edgecolors='white', linewidths=.7, zorder=3)
ax.set_yticks(y)
ax.set_yticklabels([r['term'] for r in rows])
ax.tick_params(axis='y', length=0, pad=10)
ax.set_xlabel('Fold enrichment versus the tested-gene background', labelpad=12)
ax.set_xlim(1.15, 2.0)
ax.set_ylim(-.7, len(rows)-.3)
ax.grid(axis='x', color='#E0E4E8', linewidth=.7)
for spine in ['top', 'right', 'left']:
    ax.spines[spine].set_visible(False)
ax.spines['bottom'].set_color('#8090A0')
for yi, r in zip(y, rows):
    ax.annotate(str(r['overlap']), (r['fold_enrichment'], yi), xytext=(14, 0),
                textcoords='offset points', va='center', fontsize=9, color='#223446')
cax = fig.add_axes([.855, .38, .018, .34])
cb = fig.colorbar(points, cax=cax)
cb.set_label('−log₁₀(GO-term BH FDR)', labelpad=10)
cb.outline.set_visible(False)
fig.text(.035, .95, 'GO Biological Process enrichment', fontsize=17, weight='bold', color='#123B60')
fig.text(.035, .90, 'Positive list: no significant terms     |     Negative list: 212 significant terms', fontsize=12)
fig.text(.035, .855, 'Selected terms from the negative list (1,936 input genes; 1,894 mapped)', fontsize=11, color='#41576A')
fig.text(.035, .13, 'Dot size and labels: contributing genes. Background: 20,200 tested symbols; 19,350 mapped genes.', fontsize=9)
fig.text(.035, .09, 'Inputs selected at nominal p < 0.05; no individual gene passes association FDR < 0.05. GO FDR is corrected per list.', fontsize=9)
fig.text(.035, .05, 'Gene body + 2 kb; log₁₀ activity and MLS; log₁₀ length and mass covariates. Human GO annotations (g:Profiler).', fontsize=9)
fig.savefig(out / 'GO_enrichment_summary.png', dpi=220)
fig.savefig(out / 'GO_enrichment_summary.pdf')
plt.close(fig)
