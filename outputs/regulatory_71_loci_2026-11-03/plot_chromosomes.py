"""Chromosome-scale map of the original 71 regulatory/MLS hits (hg38)."""
from pathlib import Path
import os, json, gzip, csv
OUT = Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR', str(OUT / '.matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
from matplotlib.font_manager import FontProperties

ROOT = OUT.parents[1]
d = json.loads((OUT / 'locus_analysis.json').read_text())
chroms = [f'chr{i}' for i in range(1, 23)]
sizes = {}
for line in (ROOT / 'pipelines/atac_super_enhancers/reference_sizes/Human.chrom.sizes').read_text().splitlines():
    c, n = line.split('\t')
    if c in chroms: sizes[c] = int(n)
bands = {c: [] for c in chroms}
with gzip.open(OUT / 'hg38.cytoBandIdeo.txt.gz', 'rt') as f:
    for line in f:
        c, start, end, name, stain = line.rstrip().split('\t')
        if c in bands: bands[c].append((int(start), int(end), name, stain))
genes = {c: sorted([g for g in d['genes'] if g['chrom'] == c], key=lambda g: g['tss_1based']) for c in chroms}
assert sum(map(len, genes.values())) == 69
for c in chroms:
    assert min(r[0] for r in bands[c]) == 0
    assert max(r[1] for r in bands[c]) == sizes[c]
    assert all(1 <= g['tss_1based'] <= sizes[c] for g in genes[c])
counts = [{'chromosome': c, 'length_bp': sizes[c], 'mapped_hits': len(genes[c]), 'genes': '; '.join(g['gene'] for g in genes[c])} for c in chroms]
with open(OUT / 'chromosome_hit_counts.tsv', 'w') as f:
    w = csv.DictWriter(f, fieldnames=list(counts[0]), delimiter='\t'); w.writeheader(); w.writerows(counts)

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'svg.fonttype': 'none', 'pdf.fonttype': 42, 'ps.fonttype': 42})
blue = '#0072B2'; dark = '#183D57'
fig, ax = plt.subplots(figsize=(16, 9))
fig.subplots_adjust(left=.06, right=.985, top=.865, bottom=.16)
stains = {'gneg':'#FFFFFF','gpos25':'#D6D6D6','gpos50':'#A3A3A3','gpos75':'#757575',
          'gpos100':'#404040','gvar':'#E4E4E4','stalk':'#CCCCCC','acen':'#C5B0AF'}
ys = {c: (len(chroms)-1-i)*1.5 for i, c in enumerate(chroms)}
bar_h = .18
ax.set_ylim(-.9, max(ys.values())+.9)
ax.set_xlim(-2,253)
fig.canvas.draw()
renderer=fig.canvas.get_renderer()
font=FontProperties(family='DejaVu Sans',size=8.8)
label_boxes=[]

def label_positions(items):
    """Minimize label displacement while preserving genomic order and spacing."""
    if not items: return []
    widths=[renderer.get_text_width_height_descent(g['gene'],font,False)[0]/ax.bbox.width*255 for g in items]
    offsets=[0.]
    for i in range(1,len(items)):
        offsets.append(offsets[-1]+(widths[i-1]+widths[i])/2+1.7)
    targets=[(g['tss_1based']-1)/1e6-offsets[i] for i,g in enumerate(items)]
    blocks=[]
    for i,t in enumerate(targets):
        blocks.append([t,1,[i]])
        while len(blocks)>1 and blocks[-2][0]>blocks[-1][0]:
            b=blocks.pop(); a=blocks.pop(); n=a[1]+b[1]
            blocks.append([(a[0]*a[1]+b[0]*b[1])/n,n,a[2]+b[2]])
    fitted=[0.]*len(items)
    for value,_,indices in blocks:
        for i in indices: fitted[i]=value+offsets[i]
    shift=max(0,widths[0]/2-fitted[0])
    fitted=[x+shift for x in fitted]
    shift=max(0,fitted[-1]+widths[-1]/2-250)
    fitted=[x-shift for x in fitted]
    assert fitted[0]-widths[0]/2>=-.001
    return fitted
for c in chroms:
    y = ys[c]; length = sizes[c] / 1e6
    for start, end, name, stain in bands[c]:
        start /= 1e6; end /= 1e6
        if stain == 'acen':
            points = [(start,y-bar_h/2),(start,y+bar_h/2),(end,y)] if name.startswith('p') else [(end,y-bar_h/2),(end,y+bar_h/2),(start,y)]
            ax.add_patch(Polygon(points, closed=True, facecolor=stains[stain], edgecolor='#888888', linewidth=.25, zorder=2))
        else:
            ax.add_patch(Rectangle((start,y-bar_h/2),end-start,bar_h,facecolor=stains[stain],edgecolor='none',zorder=2))
    # Separate arm outlines preserve the centromere notch.
    cents = [b for b in bands[c] if b[3]=='acen']
    c0, c1 = min(b[0] for b in cents)/1e6, max(b[1] for b in cents)/1e6
    for lo,hi in [(0,c0),(c1,length)]:
        ax.add_patch(Rectangle((lo,y-bar_h/2),hi-lo,bar_h,fill=False,edgecolor='#5C6166',linewidth=.55,zorder=3))
    # One label row per chromosome fits the widescreen layout without overlap.
    for side in [1]:
        items=genes[c]
        for g,tx in zip(items,label_positions(items)):
            x=(g['tss_1based']-1)/1e6
            ax.plot([x,x],[y+side*.09,y+side*.18],color=blue,lw=.9,zorder=4)
            ax.plot(x,y+side*.18,'o',ms=3.0,color=blue,zorder=6)
            ax.plot([x,x,tx,tx],[y+side*.18,y+side*.22,y+side*.31,y+side*.35],color=blue,lw=.6,alpha=.8,zorder=4)
            text=ax.text(tx,y+side*.38,g['gene'],ha='center',va='bottom' if side>0 else 'top',fontproperties=font,color=dark,zorder=7)
            label_boxes.append(text)
    ax.axhline(y-.75,color='#EAEDEF',lw=.55,zorder=0)

ax.set_yticks([ys[c] for c in chroms],[c.removeprefix('chr') for c in chroms])
ax.tick_params(axis='y',length=0,pad=10,labelsize=10)
ax.xaxis.set_major_locator(MultipleLocator(50));ax.xaxis.set_minor_locator(MultipleLocator(25))
ax.grid(axis='x',which='major',color='#D8DDE1',lw=.65,zorder=0)
ax.set_xlabel('Chromosome position (Mb)',labelpad=10)
ax.set_ylabel('Human chromosome',labelpad=10)
for side in ['top','right','left']:ax.spines[side].set_visible(False)
ax.spines['bottom'].set_color('#8D969C')
fig.text(.06,.967,'Chromosomal distribution of the 71 regulatory hits',ha='left',fontsize=18,fontweight='bold',color=dark)
fig.text(.06,.936,'69 mapped genes on hg38 autosomes. Original positive associations with MLS: BH FDR < 0.05.',ha='left',fontsize=10.5,color='#465763')
legend=[Line2D([0],[0],color=blue,lw=1,marker='o',ms=4,label='Gene TSS')]
fig.legend(handles=legend,loc='upper left',bbox_to_anchor=(.057,.922),frameon=False,fontsize=9.5,handlelength=1.2)
fig.text(.06,.065,'Labels are displaced for readability; connectors point to the actual TSS positions. X and Y are omitted.',fontsize=9.5,color='#465763')
fig.text(.06,.044,'TSS ±500 kb peak counts; log10(activity) vs log10(MLS), adjusted for log10(body mass). MS4A20 and ZFP865 are unmapped.',fontsize=9,color='#465763')
fig.text(.06,.023,'Coordinates: local GENCODE v50 / hg38. Cytobands: UCSC hg38 cytoBandIdeo.',fontsize=9,color='#65747D',url='https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/')
fig.canvas.draw()
boxes=[t.get_window_extent() for t in label_boxes]
assert len(boxes)==69
overlaps=[(label_boxes[i].get_text(),label_boxes[j].get_text()) for i in range(len(boxes)) for j in range(i+1,len(boxes)) if boxes[i].overlaps(boxes[j])]
assert not overlaps, overlaps
base=OUT/'Regulatory_71_chromosome_map'
fig.savefig(str(base)+'.png',dpi=200,facecolor='white')
fig.savefig(str(base)+'.svg',facecolor='white')
fig.savefig(str(base)+'.pdf',facecolor='white')
print(json.dumps({'mapped':sum(map(len,genes.values())), 'chr19_hits':len(genes['chr19']), 'output':str(base)},indent=2))
