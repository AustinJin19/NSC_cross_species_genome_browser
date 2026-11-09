#!/usr/bin/env python3
"""Publication-sized plots of MLS correlation concordance and shared genes."""
import csv
import json
import os
import sys
from pathlib import Path

ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'cross_species_visual_tool/web_app').is_dir())
sys.path.insert(0, str(ROOT / 'cross_species_visual_tool/web_app/.deps'))
os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).resolve().parent / 'mpl_cache'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list

OUT = ROOT / 'outputs/activity_omics_comparison_2026-10-20'
ATAC = ['tss2kb','tss5kb','body','tss200kb','tss500kb']
ALL = ATAC + ['rna','protein']
LABELS = ['TSS ±2 kb','TSS ±5 kb','Gene body + 2 kb','TSS ±200 kb','TSS ±500 kb','Transcriptome','Proteome']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,
                     'axes.linewidth':.8,'savefig.facecolor':'white'})


def records(path):
    with path.open() as h:
        return list(csv.DictReader(h, delimiter='\t'))


def save(fig, name):
    fig.savefig(OUT / f'{name}.pdf', bbox_inches='tight')
    fig.savefig(OUT / f'{name}.png', dpi=240, bbox_inches='tight')
    plt.close(fig)


def matrix_plot(ax, values, rows, cols, cmap, vmin, vmax, fmt, title, bar_label):
    im=ax.imshow(values,cmap=cmap,vmin=vmin,vmax=vmax,aspect='auto')
    ax.set_xticks(range(len(cols)),cols); ax.set_yticks(range(len(rows)),rows)
    ax.tick_params(length=0,pad=8)
    ax.set_xticks(np.arange(-.5,len(cols),1),minor=True)
    ax.set_yticks(np.arange(-.5,len(rows),1),minor=True)
    ax.grid(which='minor',color='white',linewidth=1.2); ax.tick_params(which='minor',length=0)
    for i in range(len(rows)):
        for j in range(len(cols)):
            val=values[i,j]
            rgba=im.cmap(im.norm(val)); luminance=.2126*rgba[0]+.7152*rgba[1]+.0722*rgba[2]
            ax.text(j,i,format(val,fmt),ha='center',va='center',color='white' if luminance<.48 else '#172330',fontsize=11)
    ax.set_title(title,pad=12,fontsize=11)
    return im


def comparison_heatmaps():
    source=records(OUT/'comparison_summary.tsv')
    for variant in ('mass_adjusted','without_mass'):
        d={(r['method_a'],r['cohort']):r for r in source if r['variant']==variant}
        cols=[f'Transcriptome\n{int(d[(ATAC[0],"rna")]["tested_genes"]):,} genes',
              f'Proteome\n{int(d[(ATAC[0],"protein")]["tested_genes"]):,} genes']
        for threshold in ('nominal','FDR'):
            metrics=[('spearman_rho','Correlation of MLS r values','RdBu_r',-1,1,'.2f'),
                     (f'{threshold}_positive_both','Shared positive genes','Oranges',0,None,',.0f'),
                     (f'{threshold}_negative_both','Shared negative genes','Blues',0,None,',.0f')]
            fig,axes=plt.subplots(1,3,figsize=(11.7,4.5),gridspec_kw={'width_ratios':[1,1,1]})
            for k,(metric,title,cmap,vmin,vmax,fmt) in enumerate(metrics):
                arr=np.array([[float(d[(m,c)][metric]) for c in ('rna','protein')] for m in ATAC])
                if vmax is None:vmax=max(1,np.max(arr))
                if k:title+='\n'+('p < 0.05 in both' if threshold=='nominal' else 'BH q < 0.05 in both')
                im=matrix_plot(axes[k],arr,LABELS[:5] if k==0 else ['']*5,cols,cmap,vmin,vmax,fmt,title,'')
                bar=fig.colorbar(im,ax=axes[k],orientation='horizontal',fraction=.065,pad=.19)
                bar.set_label('Spearman ρ' if k==0 else 'Genes',fontsize=9)
                if k:bar.locator=matplotlib.ticker.MaxNLocator(nbins=4,integer=True);bar.update_ticks()
            fig.subplots_adjust(wspace=.25,left=.16,right=.99,bottom=.16,top=.87)
            save(fig,f'comparison_{variant}_{threshold}')


def joint_map():
    rows=records(OUT/'joint_correlation_map.tsv')
    lookup={(r['method_a'],r['method_b']):r for r in rows}
    matrix=np.array([[float(lookup[(a,b)]['spearman_rho']) for b in ALL] for a in ALL])
    fig,ax=plt.subplots(figsize=(8.6,7.3))
    im=matrix_plot(ax,matrix,LABELS,LABELS,'RdBu_r',-1,1,'.2f','', '')
    plt.setp(ax.get_xticklabels(),rotation=45,ha='right',rotation_mode='anchor')
    bar=fig.colorbar(im,ax=ax,fraction=.046,pad=.04);bar.set_label('Spearman ρ of gene–MLS correlation coefficients')
    ax.set_title(f"Shared cohort: {int(rows[0]['tested_genes']):,} genes",fontsize=10,pad=10)
    save(fig,'joint_correlation_map')


def gene_heatmap():
    fits={m:{r['gene']:r for r in records(OUT/'mass_adjusted/joint'/f'{m}_correlations.tsv')} for m in ALL}
    genes=sorted(set.intersection(*(set(v) for v in fits.values())))
    candidates=[]
    for g in genes:
        hits=[]
        for m in ATAC:
            for assay in ('rna','protein'):
                a,b=fits[m][g],fits[assay][g]
                if float(a['r'])*float(b['r'])>0 and max(float(a['p']),float(b['p']))<.05:
                    hits.append(max(float(a['p']),float(b['p'])))
        if hits:candidates.append((min(hits),g))
    selected=[g for _,g in sorted(candidates)[:50]]
    if not selected:return
    values=np.array([[float(fits[m][g]['r']) for m in ALL] for g in selected])
    order=leaves_list(linkage(values,method='average',metric='euclidean')) if len(selected)>1 else [0]
    values=values[order]; selected=[selected[i] for i in order]
    fig,ax=plt.subplots(figsize=(8.8,max(5,len(selected)*.235+1.6)))
    im=ax.imshow(values,cmap='RdBu_r',vmin=-1,vmax=1,aspect='auto')
    ax.set_xticks(range(7),LABELS);ax.set_yticks(range(len(selected)),selected)
    ax.tick_params(length=0,pad=5,labelsize=9)
    plt.setp(ax.get_xticklabels(),rotation=45,ha='right',rotation_mode='anchor')
    for i,g in enumerate(selected):
        for j,m in enumerate(ALL):
            r=fits[m][g]
            if float(r['p'])<.05:ax.text(j,i,'*',ha='center',va='center',fontsize=11,color='white' if abs(float(r['r']))>.62 else 'black')
    bar=fig.colorbar(im,ax=ax,fraction=.038,pad=.03);bar.set_label('MLS partial r')
    bar.ax.set_title('* p < 0.05',fontsize=9,pad=12)
    save(fig,'shared_gene_correlation_heatmap')
    with (OUT/'heatmap_genes.tsv').open('w') as h:
        w=csv.writer(h,delimiter='\t');w.writerow(['gene','n_species',*[f'{m}_r' for m in ALL],*[f'{m}_p' for m in ALL]])
        for g in selected:w.writerow([g,fits['body'][g]['n_species'],*[fits[m][g]['r'] for m in ALL],*[fits[m][g]['p'] for m in ALL]])


if __name__=='__main__':
    comparison_heatmaps();joint_map();gene_heatmap()
