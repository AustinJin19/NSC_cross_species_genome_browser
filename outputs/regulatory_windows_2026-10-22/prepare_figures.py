"""Create publication-ready window comparison figures and supplemental report data."""
import csv
import gzip
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/nsc-regulatory-window-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
d = json.loads((ROOT/'report_data.json').read_text())
windows = [r['window_kb'] for r in d['summary']]
all_rows = []
for w in windows:
    with gzip.open(ROOT/f'tss_{w}kb/all_tested.tsv.gz', 'rt') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            for k in ['window_kb', 'n_species', 'df']:
                r[k] = int(r[k])
            for k in ['partial_r','p_value','MLS_coefficient','mass_coefficient_per_log10g','model_intercept','BH_FDR_window','BH_FDR_all_windows']:
                r[k] = float(r[k])
            all_rows.append(r)
d['joint_significant'] = [dict(r, direction='Positive' if r['partial_r']>0 else 'Negative')
                          for r in all_rows if r['BH_FDR_all_windows'] < .05]
(ROOT/'report_data.json').write_text(json.dumps(d, indent=2))

plt.rcParams.update({'font.family':'Arial','font.size':10,'pdf.fonttype':42,'ps.fonttype':42,
                     'axes.linewidth':.8,'savefig.facecolor':'white','axes.facecolor':'white'})
fig, axes = plt.subplots(1,2,figsize=(12.5,8.5),gridspec_kw={'width_ratios':[1,1.4]},layout='constrained')
ax=axes[0]
y=np.arange(len(windows))
positive=np.array([r['FDR_positive'] for r in d['summary']])
negative=np.array([r['FDR_negative'] for r in d['summary']])
ax.barh(y,positive,color='#3C5488',height=.6,label='Positive')
ax.barh(y,negative,left=positive,color='#E64B35',height=.6,label='Negative')
for i,(a,b) in enumerate(zip(positive,negative)):
    ax.text(a+b+2,i,f'{a} + / {b} −',va='center',fontsize=10)
ax.set(yticks=y,yticklabels=[f'±{w:,} kb' for w in windows],xlabel='Genes with window BH FDR < 0.05',ylabel='TSS window',xlim=(0,max(positive+negative)*1.33))
ax.invert_yaxis()
ax.legend(frameon=False,loc='upper right')
ax.grid(axis='x',color='#E5E5E5',linewidth=.6)
ax.set_axisbelow(True)

genes=[r['gene'] for r in d['comparison'] if r['significant_windows']>=3]
genes += [g for g in ['SIRT6','FBXO25'] if g not in genes]
lookup={(r['gene'],r['window_kb']):r for r in all_rows}
values=np.array([[lookup.get((g,w),{}).get('partial_r',np.nan) for w in windows] for g in genes])
flags=np.array([[lookup.get((g,w),{}).get('BH_FDR_window',1)<.05 for w in windows] for g in genes])
ax=axes[1]
im=ax.imshow(values,cmap='RdBu_r',vmin=-1,vmax=1,aspect='auto')
for i,j in zip(*np.where(flags)):
    ax.text(j,i,'●',ha='center',va='center',fontsize=7,color='black')
ax.set(yticks=np.arange(len(genes)),yticklabels=genes,xticks=np.arange(len(windows)),
       xticklabels=[f'±{w:,}' for w in windows],xlabel='TSS window (kb)')
ax.tick_params(axis='both',length=0)
ax.set_xticks(np.arange(-.5,len(windows),1),minor=True)
ax.set_yticks(np.arange(-.5,len(genes),1),minor=True)
ax.grid(which='minor',color='white',linewidth=.7)
ax.tick_params(which='minor',bottom=False,left=False)
cb=fig.colorbar(im,ax=ax,shrink=.5,pad=.03)
cb.set_label('Partial r')
ax.text(0,-.08,'●  Window BH FDR < 0.05',ha='left',va='top',transform=ax.transAxes,fontsize=10)
fig.savefig(ROOT/'regulatory_MLS_window_comparison.png',dpi=220,bbox_inches='tight')
fig.savefig(ROOT/'regulatory_MLS_window_comparison.pdf',bbox_inches='tight')
plt.close(fig)

(ROOT/'README.txt').write_text('''REGULATORY ACTIVITY AND MLS ACROSS TSS WINDOWS | 22 October 2026

Model: log10(species-median ATAC peak count) ~ log10(MLS in years) + log10(body mass in grams).
Windows are symmetric TSS flanks: +/-20, 50, 100, 200, 500 and 1,000 kb. Total widths are 40, 100, 200, 400, 1,000 and 2,000 kb.
Primary significance: Benjamini-Hochberg FDR < 0.05 separately across tested genes within each window, matching the earlier 71-gene analysis. Additional joint FDR pools all tested gene-window associations across the six windows. Joint BH thresholds can be higher or lower than per-window thresholds; joint FDR is a different correction family, not always more conservative for an individual window.

Count sample-specific peaks overlapping each clipped genomic interval; take the species median. Select the highest-median locus before filtering, without fallback to another locus. Exclude human, macaque, bowhead whale, missing/lower-bound MLS, missing mass, nonpositive medians and regions with <50% available genomic span. At least 5 species per regression. Region availability uses chromosome/scaffold bounds, not ATAC read coverage or internal assembly gaps. No gene-length covariate or pseudocount.

The 500 kb results exactly reproduce the prior 18,675 tested genes and 71 FDR-significant positive genes. All fitted associations were cross-checked with an independent full-model coefficient t-test; 2,346 locus-window examples were checked against the browser scorer. Inputs were not modified.

Files:
Regulatory_activity_MLS_six_windows.xlsx: summary; window-FDR significant genes; joint-FDR significant genes; all-window comparisons for the 216 genes significant in any window; species and methods.
tss_*kb/: all tested associations (gzip-compressed TSV); separate positive and negative lists for window FDR, joint FDR, and nominal p-value thresholds.
window_summary.tsv; significant_genes_all_windows.tsv; significant_gene_comparison.tsv; significant_gene_window_overlap.tsv: machine-readable summary outputs.
regulatory_MLS_window_comparison.pdf/.png: counts and partial-r heatmap; heatmap rows are genes significant in at least three windows plus SIRT6 and FBXO25; dots indicate window FDR < 0.05. All other genes remain in the workbook and TSV files.
species.tsv; coverage_by_window_species.tsv; manifest.json: traits, availability, source-file metadata and validation.

Species sets and selected loci can change across windows. Nested windows and nearby genes share peaks. Genes are matched by symbols, not verified one-to-one orthology. Models do not adjust for phylogeny, sequencing depth, batch, or global peak-count differences. These associations do not establish causality.
''')
print(json.dumps({'joint_significant_pairs':len(d['joint_significant']),'plot_genes':genes}))
