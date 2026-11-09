from pathlib import Path
import csv
import gzip
import hashlib
import json
import math
import os

os.environ.setdefault('MPLCONFIGDIR', '/tmp/sirt6_deeploop_mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MPath
from matplotlib.colors import Normalize, LinearSegmentedColormap

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent.parent
SOURCE = ROOT / 'GSE167200_Liver.top300K.txt.gz'
ANNO = OUT / 'hg19_ncbiRefSeq_chr19_3500000_4700000.json'
transcripts = json.loads(ANNO.read_text())['ncbiRefSeq']
canonical = next(t for t in transcripts if t['name'] == 'NM_016539.4')
gene_start, gene_end = canonical['txStart'] + 1, canonical['txEnd']
tss = gene_end
assert canonical['strand'] == '-'

def parse(s):
    chrom, coords = s.split(':')
    start, end = map(int, coords.split('-'))
    return chrom, start, end

def overlap(a, start, end):
    return a[0] == 'chr19' and a[1] <= end and a[2] >= start

def gene_names(a):
    return sorted({t['name2'] for t in transcripts if t['chrom'] == a[0]
                   and t['txStart'] < a[2] and t['txEnd'] >= a[1]})

direct, nearby = [], []
windows = {k: {'either_anchor': 0, 'both_anchors': 0} for k in (200000, 500000)}
columns = set()
previous_score = math.inf
distance_matches = 0
total = 0
with gzip.open(SOURCE, 'rt') as handle:
    for rank, line in enumerate(handle, 1):
        total += 1
        fields = line.rstrip().split('\t')
        columns.add(len(fields))
        assert len(fields) == 4
        a, b = map(parse, fields[:2])
        score, gap = float(fields[2]), int(fields[3])
        assert score <= previous_score + 1e-12
        previous_score = score
        assert a[0] == b[0] and a[2] < b[1]
        distance_matches += gap == b[1] - a[2]
        for radius in windows:
            flags = [overlap(x, tss-radius, tss+radius) for x in (a,b)]
            windows[radius]['either_anchor'] += any(flags)
            windows[radius]['both_anchors'] += all(flags)
        if any(overlap(x,tss-500000,tss+500000) for x in (a,b)):
            nearby.append({'genomewide_rank':rank,'anchor1_hg19':fields[0],
                           'anchor2_hg19':fields[1],'DeepLoop_score':score,
                           'reported_anchor_gap_bp':gap})
        contains = [overlap(x,tss,tss) for x in (a,b)]
        if not any(contains):
            continue
        idx = contains.index(True)
        own, other = (a,b)[idx], (a,b)[1-idx]
        direct.append({
            'genomewide_rank': rank,
            'SIRT6_anchor_hg19': fields[idx],
            'partner_anchor_hg19': fields[1-idx],
            'DeepLoop_score': score,
            'reported_anchor_gap_bp': gap,
            'partner_midpoint_distance_from_TSS_kb': abs((other[1]+other[2])/2-tss)/1000,
            'partner_genomic_side': 'lower-coordinate' if other[2]<tss else 'higher-coordinate',
            'partner_overlapping_genes': ';'.join(gene_names(other)) or 'Intergenic',
            'SIRT6_anchor_overlapping_genes': ';'.join(gene_names(own)),
            'adjacent_anchor_bins': gap == 1,
        })
assert total == 300000 and columns == {4} and distance_matches == total
assert len(direct) == 14
assert len({r['partner_anchor_hg19'] for r in direct}) == 14
anchors = {r['SIRT6_anchor_hg19'] for r in direct}
assert anchors == {'chr19:4147554-4186620'}
for name, rows in [('SIRT6_TSS_anchor_contacts_hg19.tsv',direct),
                   ('SIRT6_plusminus500kb_any_anchor_hg19.tsv',nearby)]:
    with (OUT/name).open('w',newline='') as h:
        writer=csv.DictWriter(h, fieldnames=list(rows[0]),delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

shared = parse(next(iter(anchors)))
summary = {
    'source_file':SOURCE.name,
    'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    'genome_build':'hg19', 'organism':'Homo sapiens', 'tissue':'Liver',
    'source_study':'GSE167200 / DeepLoop; original liver Hi-C GSM1419086',
    'total_ranked_contact_pixels':total,
    'columns':['anchor1','anchor2','DeepLoop score','anchor2.start - anchor1.end (bp)'],
    'gene_reference':{'transcript':'NM_016539.4','start_1based':gene_start,'end_1based':gene_end,
                      'strand':'-','TSS_1based':tss},
    'SIRT6_TSS_anchor':next(iter(anchors)),
    'anchor_width_bp_inclusive':shared[2]-shared[1]+1,
    'anchor_gene_overlaps':gene_names(shared),
    'TSS_overlapping_contact_pixels':len(direct),
    'adjacent_bin_contacts':sum(r['adjacent_anchor_bins'] for r in direct),
    'other_contact_pixels':sum(not r['adjacent_anchor_bins'] for r in direct),
    'nearby_window_counts':windows,
    'interpretation':'Ranked DeepLoop contact pixels, not 14 independent or experimentally validated SIRT6-specific loops. Multiple adjacent partner bins may represent one broad contact feature.',
    'coordinate_handling':'Original anchor coordinates retained verbatim. Contiguous intervals have next.start = previous.end + 1, so displayed widths use inclusive coordinates. The overlap results are unchanged by a one-base boundary convention.',
    'sources':[
        'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167200',
        'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167200/suppl/GSE167200_README.txt',
        'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167200/suppl/GSE167200_Readme_39cool.txt',
        'https://github.com/JinLabBioinfo/DeepLoop#about-loopcalling',
        'https://api.genome.ucsc.edu/getData/track?genome=hg19;track=ncbiRefSeq;chrom=chr19;start=3500000;end=4700000',
    ],
}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2))

# Scientific figure: contacts use anchor centers, with the full anchor shaded.
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none'})
fig = plt.figure(figsize=(12.8,8.3))
gs = fig.add_gridspec(3,1,height_ratios=[3.4,1.8,2.0],hspace=.42)
ax, genes_ax, detail = [fig.add_subplot(gs[i]) for i in range(3)]
norm = Normalize(min(r['DeepLoop_score'] for r in direct),max(r['DeepLoop_score'] for r in direct))
cmap = LinearSegmentedColormap.from_list('loop_blues', [plt.get_cmap('Blues')(.35+i*.65/255) for i in range(256)])
xmin,xmax=3.72,4.24
own_mid=(shared[1]+shared[2])/2/1e6
for r in reversed(direct):
    partner=parse(r['partner_anchor_hg19'])
    x=(partner[1]+partner[2])/2/1e6
    height=(own_mid-x)*.75
    color=cmap(norm(r['DeepLoop_score']))
    vertices=[(x,0),(x,height),(own_mid,height),(own_mid,0)]
    ax.add_patch(PathPatch(MPath(vertices,[MPath.MOVETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4]),
                           facecolor='none',edgecolor=color,lw=1.8))
    ax.plot([partner[1]/1e6,partner[2]/1e6],[0,0],color=color,lw=4,solid_capstyle='butt')
ax.axvspan(shared[1]/1e6,shared[2]/1e6,color='#E8A46A',alpha=.25)
ax.axvline(tss/1e6,color='#B63F32',lw=1.2,ls='--')
ax.text(tss/1e6,.31,'SIRT6 TSS',color='#B63F32',ha='center',va='bottom')
ax.set(xlim=(xmin,xmax),ylim=(-.008,.34))
ax.set_title('Human liver DeepLoop contacts anchored over SIRT6 (hg19)',loc='left',fontweight='bold',pad=17)
ax.set_yticks([])
ax.tick_params(axis='x',labelbottom=False)
ax.spines[['top','right','left']].set_visible(False)
cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=ax,fraction=.018,pad=.015)
cb.set_label('DeepLoop score')

# Merge transcript extents by gene for an uncluttered neighborhood track.
models={}
for tr in transcripts:
    name=tr['name2']
    if tr['txEnd']<xmin*1e6 or tr['txStart']>xmax*1e6:
        continue
    if name not in models:
        models[name]=[tr['txStart']+1,tr['txEnd'],tr['strand']]
    else:
        models[name][0]=min(models[name][0],tr['txStart']+1)
        models[name][1]=max(models[name][1],tr['txEnd'])
lane_end=[]
for name,(s,e,strand) in sorted(models.items(),key=lambda item:item[1][0]):
    center=(s+e)/2/1e6
    left=min(s/1e6,center-len(name)*.0020)
    right=max(e/1e6,center+len(name)*.0020)
    lane=next((i for i,v in enumerate(lane_end) if v+.006<left),len(lane_end))
    if lane==len(lane_end):lane_end.append(right)
    else:lane_end[lane]=right
    color='#B63F32' if name=='SIRT6' else '#274E70'
    y=-lane
    genes_ax.plot([s/1e6,e/1e6],[y,y],lw=3,color=color)
    genes_ax.plot([center],[y],marker='>' if strand=='+' else '<',color=color,ms=5)
    genes_ax.text(center,y-.15,name,ha='center',va='top',fontsize=8,color=color)
genes_ax.axvspan(shared[1]/1e6,shared[2]/1e6,color='#E8A46A',alpha=.25)
genes_ax.set(xlim=(xmin,xmax),ylim=(-len(lane_end)+.3,.6),xlabel='Chromosome 19 position (Mb, hg19)')
genes_ax.set_yticks([])
genes_ax.spines[['top','right','left']].set_visible(False)

detail.axvspan(shared[1]/1e6,shared[2]/1e6,color='#E8A46A',alpha=.25)
detail.axvline(tss/1e6,color='#B63F32',lw=1.2,ls='--')
for i,name in enumerate(['CREB3L3','SIRT6','ANKRD24']):
    s,e,strand=models[name]
    y=2-i
    color='#B63F32' if name=='SIRT6' else '#274E70'
    detail.plot([s/1e6,e/1e6],[y,y],lw=3,color=color)
    detail.plot([(s+e)/2/1e6],[y],marker='>' if strand=='+' else '<',color=color,ms=7)
    detail.text(4.131,y,name,ha='left',va='center',fontsize=9,color=color)
detail.set(xlim=(4.13,4.23),ylim=(-.55,2.5),xlabel='Chromosome 19 position (Mb, hg19)')
detail.set_yticks([])
detail.spines[['top','right','left']].set_visible(False)
detail.set_title('Shared 39.1 kb anchor overlaps three genes',loc='left',fontsize=11,fontweight='bold')
fig.subplots_adjust(left=.06,right=.94,top=.94,bottom=.095)
arc_pos, gene_pos = ax.get_position(), genes_ax.get_position()
genes_ax.set_position([arc_pos.x0,gene_pos.y0,arc_pos.width,gene_pos.height])
fig.text(.06,.025,'14 ranked contact pixels; adjacent bins can belong to one broader contact feature. These contacts are not SIRT6-specific.',fontsize=9,color='#4B5563')
fig.savefig(OUT/'SIRT6_liver_contacts_hg19.png',dpi=180)
fig.savefig(OUT/'SIRT6_liver_contacts_hg19.svg')
plt.close(fig)
print(json.dumps({k:v for k,v in summary.items() if k not in ('sources','source_sha256')},indent=2))
