#!/usr/bin/env python3
"""Direct across-gene ATAC/expression correlations separately within species."""
import csv
import gzip
import json
import math
import os
import tempfile
from collections import defaultdict
from pathlib import Path

from prepare import ROOT, OUT, METHODS, write, query
import numpy as np
from scipy.stats import false_discovery_control, pearsonr, spearmanr, rankdata
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'nsc_direct_omics_mpl'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

LABELS={'tss2kb':'TSS ±2 kb','tss5kb':'TSS ±5 kb','body':'Body + 2 kb','tss200kb':'TSS ±200 kb','tss500kb':'TSS ±500 kb'}


def read(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as handle:
        yield from csv.DictReader(handle,delimiter='\t')


def correlate(x,y):
    result=dict(n_genes=len(x),zero_activity_genes=int(np.sum(x==0)))
    if len(x)<5 or np.ptp(x)==0 or np.ptp(y)==0:
        return dict(**result,status='insufficient or constant data')
    pearson=pearsonr(x,y);spearman=spearmanr(x,y)
    np.testing.assert_allclose(float(pearson.statistic),np.corrcoef(x,y)[0,1],rtol=1e-10,atol=1e-12)
    np.testing.assert_allclose(float(spearman.statistic),np.corrcoef(rankdata(x),rankdata(y))[0,1],rtol=1e-10,atol=1e-12)
    result.update(status='ok',pearson_raw_r=float(pearson.statistic),pearson_raw_p=float(pearson.pvalue),
                  spearman_rho=float(spearman.statistic),spearman_p=float(spearman.pvalue))
    positive=x>0
    result['n_log10_genes']=int(positive.sum())
    if positive.sum()>=5 and np.ptp(x[positive])>0 and np.ptp(y[positive])>0:
        logged=pearsonr(np.log10(x[positive]),y[positive])
        result.update(pearson_log10_r=float(logged.statistic),pearson_log10_p=float(logged.pvalue))
    return result


def rna_protein_correlations(expr,species):
    """Compare the supplied expression scales directly, without another log."""
    rows=[]
    for spec in species:
        sid=spec['species']
        rna,protein=expr[(sid,'rna')],expr[(sid,'protein')]
        genes=sorted(set(rna)&set(protein))
        pairs=[dict(gene=g,rna_logTPM_median=rna[g][0],protein_normalized_median=protein[g][0],
                    RNA_samples=rna[g][1],protein_samples=protein[g][1]) for g in genes]
        write(OUT/'rna_protein_paired_values'/f'{sid}.tsv.gz',pairs,
              ['gene','rna_logTPM_median','protein_normalized_median','RNA_samples','protein_samples'])
        row=dict(species=sid,common_name=spec['common'],n_genes=len(genes),status='missing paired assay')
        if len(genes)>=5:
            x=np.array([rna[g][0] for g in genes]);y=np.array([protein[g][0] for g in genes])
            assert np.isfinite(x).all() and np.isfinite(y).all()
            if np.ptp(x)>0 and np.ptp(y)>0:
                pr,sr=pearsonr(x,y),spearmanr(x,y)
                np.testing.assert_allclose(pr.statistic,np.corrcoef(x,y)[0,1],atol=1e-12)
                np.testing.assert_allclose(sr.statistic,np.corrcoef(rankdata(x),rankdata(y))[0,1],atol=1e-12)
                row.update(status='ok',pearson_r=float(pr.statistic),pearson_p=float(pr.pvalue),
                           spearman_rho=float(sr.statistic),spearman_p=float(sr.pvalue))
            else:
                row['status']='constant expression profile'
        elif rna and protein:
            row['status']='fewer than five shared genes'
        rows.append(row)
    valid=[r for r in rows if r['status']=='ok']
    for prefix in ['pearson','spearman']:
        for row,q in zip(valid,false_discovery_control([r[prefix+'_p'] for r in valid])):
            row[prefix+'_q']=float(q)
    write(OUT/'transcriptome_proteome_correlations.tsv',rows,
          ['species','common_name','status','n_genes','pearson_r','pearson_p','pearson_q','spearman_rho','spearman_p','spearman_q'])
    summary=dict(species=len(valid),median_Pearson_r=float(np.median([r['pearson_r'] for r in valid])),
                 min_Pearson_r=min(r['pearson_r'] for r in valid),max_Pearson_r=max(r['pearson_r'] for r in valid),
                 median_Spearman_rho=float(np.median([r['spearman_rho'] for r in valid])),
                 min_genes=min(r['n_genes'] for r in valid),max_genes=max(r['n_genes'] for r in valid))
    return rows,summary


def heatmap(rows,species,metric,title,cohort='common_five_regions',rna_protein=None):
    lookup={(r['species'],r['assay'],r['method']):r for r in rows if r['cohort']==cohort}
    rna_protein={r['species']:r for r in (rna_protein or [])}
    fig=plt.figure(figsize=(17,12.8))
    axes=[fig.add_axes([.18,.17,.29,.76]),fig.add_axes([.51,.17,.29,.76])]
    rp_ax=fig.add_axes([.84,.17,.075,.76])
    cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#eeeeee')
    for ax,assay,assay_title in zip(axes,['rna','protein'],['Transcriptome: source logTPM','Proteome: source normalized log abundance']):
        vals=np.full((len(species),5),np.nan)
        labels=[]
        for i,s in enumerate(species):
            gene_counts=[]
            for j,m in enumerate(METHODS):
                row=lookup.get((s['species'],assay,m),{})
                vals[i,j]=row.get(metric,np.nan)
                gene_counts.append(row.get('n_log10_genes' if metric=='pearson_log10_r' else 'n_genes',0))
            label=s['common'].replace('Bramus lutescens','Transcaucasian mole vole')
            if s['species']=='ASM':label+=' †'
            low,high=min(gene_counts),max(gene_counts)
            count_label=f'{low:,}' if low==high else f'{low:,}–{high:,}'
            if cohort=='common_five_regions':
                labels.append(f'{label} ({count_label})' if high else label+' (NA)')
            else:
                labels.append(label if high else label+' (NA)')
        im=ax.imshow(vals,cmap=cmap,vmin=-1,vmax=1,aspect='auto',interpolation='nearest')
        ax.set_xticks(range(5),[LABELS[m] for m in METHODS],rotation=30,ha='right',fontsize=10)
        ax.set_yticks(range(len(species)),labels if assay=='rna' else ['']*len(species),fontsize=9)
        ax.tick_params(length=0,pad=5)
        ax.set_title(assay_title,fontsize=12,pad=12)
        for i in range(len(species)):
            for j in range(5):
                v=vals[i,j]
                content='—' if not np.isfinite(v) else f'{v:.2f}'
                if np.isfinite(v):
                    row=lookup[(species[i]['species'],assay,METHODS[j])]
                    count=row['n_log10_genes'] if metric=='pearson_log10_r' else row['n_genes']
                    content+=f"\n(n={count:,})"
                ax.text(j,i,content,ha='center',va='center',fontsize=8,
                        color='white' if np.isfinite(v) and abs(v)>.65 else '#222222')
        ax.set_xticks(np.arange(-.5,5,1),minor=True)
        ax.set_yticks(np.arange(-.5,len(species),1),minor=True)
        ax.grid(which='minor',color='white',lw=.7)
        ax.tick_params(which='minor',length=0)
    # RNA/protein share one gene set per species; ATAC windows do not apply.
    rp_metric='spearman_rho' if metric=='spearman_rho' else 'pearson_r'
    rp_vals=np.array([[rna_protein.get(s['species'],{}).get(rp_metric,np.nan)] for s in species])
    rp_ax.imshow(rp_vals,cmap=cmap,vmin=-1,vmax=1,aspect='auto',interpolation='nearest')
    rp_ax.set_title('RNA ↔ protein',fontsize=12,pad=12)
    rp_ax.set_xticks([0],['All shared genes'],rotation=30,ha='right',fontsize=10)
    rp_ax.set_yticks(range(len(species)),['']*len(species))
    rp_ax.tick_params(length=0)
    rp_ax.set_yticks(np.arange(-.5,len(species),1),minor=True)
    rp_ax.grid(which='minor',color='white',lw=.7);rp_ax.tick_params(which='minor',length=0)
    for i,s in enumerate(species):
        row=rna_protein.get(s['species'],{})
        value=row.get(rp_metric,np.nan)
        content='—' if not np.isfinite(value) else f"{value:.2f}\n(n={row['n_genes']:,})"
        rp_ax.text(0,i,content,ha='center',va='center',fontsize=8,
                   color='white' if np.isfinite(value) and abs(value)>.65 else '#222222')
    fig.suptitle(title,fontsize=16,y=.985)
    cax=fig.add_axes([.94,.25,.012,.50])
    fig.colorbar(im,cax=cax,label='Correlation',ticks=[-1,-.5,0,.5,1])
    cohort_note=('Same genes across the five regions within each species/assay.' if cohort=='common_five_regions' else
                 'Each cell uses all available genes matched for that region/assay; counts may differ between regions.')
    if metric=='pearson_log10_r':
        cohort_note='Shared five-region gene set; nonpositive activity is excluded separately by region, giving the displayed count ranges.'
    fig.text(.18,.025,'One observation per gene; median across samples within each species. Numbers in parentheses: matched genes.\n'+
             cohort_note+' Gray = no matching assay.\n'
             'RNA ↔ protein uses all shared expression genes, independently of ATAC availability; neither expression scale is logged again.\n'
             '† Acomys RNA species/common-name fields disagree; see METHODS.txt.',fontsize=9,ha='left',va='bottom')
    return fig


def plot_results(results,species,rna_protein=None):
    plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42})
    with PdfPages(OUT/'within_species_activity_RNA_protein.pdf') as pdf:
        for key,metric,title in [
            ('pearson_raw','pearson_raw_r','Within-species gene activity, transcriptome and proteome correlations (Pearson r)'),
            ('spearman','spearman_rho','Within-species gene activity, transcriptome and proteome correlations (Spearman ρ)'),
            ('pearson_log10','pearson_log10_r','Within-species correlations with log₁₀(ATAC activity) (Pearson r)')]:
            fig=heatmap(results,species,metric,title,rna_protein=rna_protein)
            fig.savefig(OUT/f'{key}_heatmap.png',dpi=220,bbox_inches='tight',facecolor='white')
            pdf.savefig(fig,bbox_inches='tight',facecolor='white');plt.close(fig)
    with PdfPages(OUT/'all_available_genes_activity_RNA_protein.pdf') as pdf:
        for key,metric,title in [
            ('pearson_raw','pearson_raw_r','Within-species activity and expression correlations: all available matched genes (Pearson r)'),
            ('spearman','spearman_rho','Within-species activity and expression correlations: all available matched genes (Spearman ρ)')]:
            fig=heatmap(results,species,metric,title,cohort='all_available_per_region',rna_protein=rna_protein)
            fig.savefig(OUT/f'all_available_{key}_heatmap.png',dpi=220,bbox_inches='tight',facecolor='white')
            pdf.savefig(fig,bbox_inches='tight',facecolor='white');plt.close(fig)


def all_available_summary(results):
    rows=[]
    for assay in ['rna','protein']:
        for method in METHODS:
            selected=[r for r in results if r['cohort']=='all_available_per_region' and r['assay']==assay and r['method']==method and r['status']=='ok']
            rows.append(dict(assay=assay,method=method,species=len(selected),
                             min_genes=min(r['n_genes'] for r in selected),max_genes=max(r['n_genes'] for r in selected),
                             median_Pearson_r=float(np.median([r['pearson_raw_r'] for r in selected])),
                             median_Spearman_rho=float(np.median([r['spearman_rho'] for r in selected]))))
    write(OUT/'all_available_method_summary.tsv',rows,list(rows[0]))
    return rows


def main():
    assert not query('SIRT6')['stale'], 'Rebuild expression index and prepared inputs first.'
    manifest=json.loads((OUT/'inputs.json').read_text())
    species=sorted(manifest['species'],key=lambda s:s['common'].lower())
    expr=defaultdict(dict)
    for r in read(OUT/'expression_species_medians.tsv.gz'):
        expr[(r['species'],r['assay'])][r['gene']]=(float(r['median']),int(r['n_samples']))
    rp_rows,rp_summary=rna_protein_correlations(expr,species)
    results=[];counts=[];exclusions=[]
    for spec in species:
        sid=spec['species']
        activity={m:{} for m in METHODS}
        stats=defaultdict(int)
        for r in read(Path(spec['activity_cache'])):
            method,gene=r['method'],r['gene']
            assert gene not in activity[method],(sid,method,gene)
            value=float(r['median_activity'])
            if r['region_eligible']!='1':
                stats['below_50pct_available']+=1;continue
            if not math.isfinite(value) or value<0:
                stats['nonfinite_or_negative']+=1;continue
            activity[method][gene]=r
        exclusions.append(dict(species=sid,**stats))
        union=set.union(*(set(activity[m]) for m in METHODS))
        paired=[]
        for gene in sorted(union & (set(expr[(sid,'rna')]) | set(expr[(sid,'protein')]))):
            record=dict(species=sid,gene=gene)
            for m in METHODS:
                r=activity[m].get(gene)
                if r:
                    record[m+'_CPM_bp']=float(r['median_activity'])
                    record[m+'_ATAC_samples']=int(r['n_samples'])
                    record[m+'_available_fraction']=float(r['available_fraction'])
                    record[m+'_locus']=r['locus']
            for assay in ['rna','protein']:
                val=expr[(sid,assay)].get(gene)
                if val:record[assay+'_median'],record[assay+'_samples']=val
            paired.append(record)
        fields=['species','gene']+[f'{m}_{suffix}' for m in METHODS for suffix in ['CPM_bp','ATAC_samples','available_fraction','locus']]+['rna_median','rna_samples','protein_median','protein_samples']
        write(OUT/'paired_values'/f'{sid}.tsv.gz',paired,fields)
        for assay in ['rna','protein']:
            available=expr[(sid,assay)]
            if not available:
                continue
            common=set(available).intersection(*(set(activity[m]) for m in METHODS))
            counts.append(dict(species=sid,common_name=spec['common'],assay=assay,
                               expression_genes=len(available),common_five_region_genes=len(common)))
            write(OUT/'gene_sets'/f'{sid}_{assay}_common_five_regions.tsv',
                  [dict(gene=g) for g in sorted(common)],['gene'])
            for cohort in ['common_five_regions','all_available_per_region']:
                for method in METHODS:
                    genes=sorted(common if cohort=='common_five_regions' else set(activity[method]) & set(available))
                    x=np.array([float(activity[method][g]['median_activity']) for g in genes])
                    y=np.array([available[g][0] for g in genes])
                    res=correlate(x,y)
                    results.append(dict(species=sid,common_name=spec['common'],assay=assay,method=method,cohort=cohort,**res))
        print(sid+': completed direct across-gene correlations',flush=True)
    # Correct the descriptive tests separately for each metric and declared cohort.
    for cohort in ['common_five_regions','all_available_per_region']:
        for prefix in ['pearson_raw','spearman','pearson_log10']:
            valid=[r for r in results if r['cohort']==cohort and prefix+'_p' in r]
            for r,q in zip(valid,false_discovery_control([r[prefix+'_p'] for r in valid])):
                r[prefix+'_q']=float(q)
    fields=['species','common_name','assay','method','cohort','status','n_genes','zero_activity_genes',
            'pearson_raw_r','pearson_raw_p','pearson_raw_q','spearman_rho','spearman_p','spearman_q',
            'n_log10_genes','pearson_log10_r','pearson_log10_p','pearson_log10_q']
    write(OUT/'within_species_correlations.tsv',results,fields)
    write(OUT/'matched_gene_counts.tsv',counts,list(counts[0]))
    write(OUT/'activity_exclusions.tsv',exclusions,['species','below_50pct_available','nonfinite_or_negative'])
    overview=[]
    for assay in ['rna','protein']:
        for method in METHODS:
            subset=[r for r in results if r['cohort']=='common_five_regions' and r['assay']==assay and r['method']==method and r['status']=='ok']
            record=dict(assay=assay,method=method,species=len(subset),min_genes=min(r['n_genes'] for r in subset),max_genes=max(r['n_genes'] for r in subset))
            for metric in ['pearson_raw_r','spearman_rho','pearson_log10_r']:
                valid=[r[metric] for r in subset if metric in r]
                record.update({metric+'_median':float(np.median(valid)),metric+'_min':min(valid),metric+'_max':max(valid)})
            overview.append(record)
    write(OUT/'method_summary.tsv',overview,list(overview[0]))
    all_available=all_available_summary(results)
    plot_results(results,species,rp_rows)
    summary=dict(species=len(species),RNA_species=sum(bool(expr[(s['species'],'rna')]) for s in species),
                 protein_species=sum(bool(expr[(s['species'],'protein')]) for s in species),methods=overview,
                 notes=manifest['notes'],all_available_methods=all_available,transcriptome_proteome=rp_summary)
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (OUT/'METHODS.txt').write_text(
        'DIRECT WITHIN-SPECIES GENE ACTIVITY / EXPRESSION COMPARISON\n\n'
        'Each observation is one gene, within one species. ATAC, RNA and protein are aggregated separately across samples by species medians. Samples from different assays are not assumed to be paired individuals. This is NOT an MLS analysis and not a correlation of MLS coefficients.\n\n'
        'ATAC: browser original CPM-normalized bigWigs, exact sum of signal times base-pair overlap. The activity unit is CPM·bp, not per-base CPM and not TMM. Five methods: TSS ±2 kb, ±5 kb, body + 2 kb upstream, TSS ±200 kb and ±500 kb. Fixed TSS windows have total requested widths 4, 10, 400 and 1000 kb. No additional log transform for the primary Pearson analysis. Body sums include introns and depend on gene length. No gene-length or body-mass covariates are used.\n'
        'Loci: existing browser annotations, including all gene biotypes. Select the highest median activity locus per symbol/species/method before applying the availability filter. Exclude when less than 50% of the requested window lies within chromosome/scaffold bounds; no fallback to a lower-activity locus. This is assembly-span availability, not coverage or N content. Retain zero activity for raw Pearson/Spearman. Drop nonfinite or negative values.\n'
        'Every reused cache was verified against current locus hash, bigWig paths, sizes and modification times. Missing or stale species were quantified with the existing exact bigWig integration code, including direct-query spot checks. inputs.json records cache paths and source fingerprints.\n\n'
        'RNA: species median of supplied log.tpm in liver_omics.sqlite (CrossSpecies_Count.xlsx and existing mouse liver supplement rnaLongTPM.xlsx); exclude BMR_Turk_Lu as in the browser. Protein: species median of source quant.norm normalized on log.quant from protLongDF.xlsx. These transformed source values are not logged again. Check source workbook sizes/timestamps against index before analysis.\n'
        'Species: all current browser species with matching liver RNA or protein, including bowhead whale and cave nectar bat. No MLS-value filtering. Species without a matching assay appear gray or are documented in species_coverage.tsv. Existing browser sample exclusions are preserved. Ellobius lutescens RNA is mapped to Bramus lutescens. Acomys RNA has a conflicting species/common-name field; follow the source species field and flag this row.\n\n'
        'Primary heatmaps: for each species and expression assay, intersect available genes across ALL five ATAC methods and the assay, so the five methods use exactly the same gene set. RNA and protein use their own matched gene sets; they need not be identical. Supplemental all_available_per_region results use every gene matched within that method/assay. Gene matching uses uppercase symbols; unique one-to-one orthology across species is not established or required for these within-species comparisons.\n'
        'Statistics: ordinary Pearson r of unlogged species-median CPM·bp against species-median source RNA log.tpm or protein quant.norm. Also report Spearman rho with average ranks for ties. A separate log10(activity) Pearson sensitivity drops nonpositive activity without adding a pseudocount and reports its own n_log10_genes. At least five genes and nonconstant profiles are required.\n'
        'RNA versus protein: independently intersect all available RNA and protein symbols within each species, without requiring ATAC availability or a particular window. Correlate median source log.tpm with median source quant.norm directly, with no additional logarithm or covariates. Both raw and rank correlations use the same genes; nonpositive transformed values are valid. transcriptome_proteome_correlations.tsv contains coefficients, counts and p/q; rna_protein_paired_values/ contains every gene pair. BH correction is separate for Pearson and Spearman across these species. This rightmost heatmap column remains the same when the ATAC window or transform changes.\n'
        'Two-sided p-values and BH q-values are exported, corrected separately by correlation metric and gene-matching cohort across species/method/assay tests. They are exploratory: neighboring/overlapping gene windows are not independent, and the primary interpretation should emphasize effect size. No phylogenetic regression, MLS, batch removal or other covariates. Correlations do not establish a causal regulatory link.\n'
        'The primary Pearson and Spearman coefficients are independently verified using numpy correlations of values or average ranks. Paired gene-level inputs, including sample counts, selected loci and available fractions, are saved in paired_values/.\n\n'
        'Reproduce: /usr/bin/python3 outputs/direct_activity_omics_2026-10-21/scripts/prepare.py\n'
        'Then: /usr/bin/python3 outputs/direct_activity_omics_2026-10-21/scripts/compare.py\n')
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    main()
