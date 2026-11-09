from pathlib import Path
import json, textwrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

out=Path(__file__).resolve().parent
d=json.loads((out/'locus_analysis.json').read_text())
rows=d['clusters']; n=len(rows)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none'})
fig,(ax,txt)=plt.subplots(1,2,figsize=(16,12),gridspec_kw={'width_ratios':[1.15,1.6],'wspace':0.04})
for i,r in enumerate(rows):
    y=n-i-1
    color='#D55E00' if 'SIRT6' in r['genes'] else ('#2166AC' if r['n_genes']>1 else '#888888')
    points=[(g['tss_1based']-r['first_tss_1based'])/1000 for g in d['genes'] if g.get('cluster')==r['cluster']]
    ax.plot([0,max(points)],[y,y],color=color,lw=3,alpha=.35,zorder=2)
    ax.scatter(points,[y]*len(points),color=color,s=27,zorder=3,edgecolors='white',linewidths=.3)
    label=r['genes'].replace('; ', ', ')
    txt.text(0,y,textwrap.fill(label,78),va='center',fontsize=9,color='#8B3500' if 'SIRT6' in label else '#222222')
    ax.axhline(y,color='#EFEFEF',lw=.6,zorder=0)
    txt.axhline(y,color='#EFEFEF',lw=.6,zorder=0)
labels=[f"{r['chrom']}: {r['first_tss_1based']/1e6:.3f}–{r['last_tss_1based']/1e6:.3f} Mb  ({r['n_genes']})" for r in rows]
ax.set_yticks(list(reversed(range(n))),labels)
ax.set_xlim(-30,1100); ax.set_ylim(-.7,n-.25)
ax.set_xlabel('TSS distance from the first gene in each group (kb)')
ax.xaxis.set_major_locator(MultipleLocator(250)); ax.grid(axis='x',color='#E5E5E5',lw=.65)
ax.tick_params(axis='y',length=0,pad=9)
for s in ['top','right','left']: ax.spines[s].set_visible(False)
ax.spines['bottom'].set_color('#AAAAAA')
txt.set_xlim(0,1); txt.set_ylim(ax.get_ylim()); txt.axis('off')
ax.set_title('hg38 TSS range (number of genes)',loc='left',fontsize=11,pad=16)
txt.set_title('Significant genes, in genomic order',loc='left',fontsize=11,pad=16)
fig.suptitle('Genomic locations of the 71 positive regulatory hits',x=.035,ha='left',y=.984,fontsize=17,fontweight='bold')
fig.text(.035,.951,'69 genes mapped to hg38: 63 in 18 multi-gene groups and 6 singletons. All original hits have BH FDR < 0.05.',fontsize=11)
fig.text(.035,.035,'Groups connect adjacent TSSs less than 1 Mb apart, so their ±500 kb windows overlap. Connected groups need not share one common interval.',fontsize=10)
fig.text(.035,.017,'MS4A20 and ZFP865 lack a unique human symbol/anchor mapping here; their native-species coordinates are retained in the workbook.',fontsize=10)
fig.subplots_adjust(left=.26,right=.99,top=.90,bottom=.09)
fig.savefig(out/'Regulatory_71_genomic_loci.png',dpi=180,facecolor='white')
fig.savefig(out/'Regulatory_71_genomic_loci.svg',facecolor='white')
