#!/usr/bin/env python3
"""Compare species-median TMM ATAC and liver RNA associations with lifespan."""
import csv
import gzip
import hashlib
import json
import math
import re
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

import run_pca_all_genes as base
from scipy.stats import false_discovery_control, pearsonr, spearmanr, t

np, plt = base.np, base.plt
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
WEB = ROOT / 'cross_species_visual_tool/web_app'
OUT = HERE / 'tmm_mls_transcriptome'
sys.path.insert(0, str(WEB))
from omics import query


def write(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def fit(x, y, mass=None):
    if len(x) < 5 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return None
    cov = np.ones((len(x), 1)) if mass is None else np.column_stack((np.ones(len(x)), mass))
    if np.linalg.matrix_rank(cov) != cov.shape[1]:
        return None
    q, _ = np.linalg.qr(cov, mode='reduced')
    rx, ry = x - q @ (q.T @ x), y - q @ (q.T @ y)
    xx, yy = rx @ rx, ry @ ry
    if xx < 1e-12 or yy < 1e-12:
        return None
    r = float(np.clip(rx @ ry / math.sqrt(xx * yy), -1, 1))
    df = len(x) - cov.shape[1] - 1
    p = float(2 * t.sf(abs(r) * math.sqrt(df / max(1-r*r, 1e-300)), df))
    return dict(r=r, p=p, df=df)


def add_q(rows, prefixes):
    for prefix in prefixes:
        for row, q in zip(rows, false_discovery_control([r[prefix+'p'] for r in rows])):
            row[prefix+'q'] = float(q)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    traits = json.loads((WEB / 'data/species_mls.json').read_text())
    meta = base.read_table(HERE / 'pca_scores.tsv')
    species = sorted({r['species'] for r in meta})
    common = dict(re.findall(r'^  ([A-Za-z_]+) = "([^"]+)"', (ROOT / 'plot_bundle_sample_qc.R').read_text(), re.M))
    eligible = [s for s in species if traits[s.replace('_', ' ')].get('years') is not None
                and traits[s.replace('_', ' ')].get('lower_bound_years') is None]
    assert len(species) == 32 and len(eligible) == 31
    with (HERE / 'logcpm_tmm_core.tsv').open() as handle:
        reader = csv.reader(handle, delimiter='\t')
        samples = next(reader)[1:]
        raw = list(reader)
    genes = [r[0].upper() for r in raw]
    values = np.array([r[1:] for r in raw], dtype=float)
    assert len(genes) == len(set(genes)) == 6459 and len(samples) == 75
    assert set(samples) == {r[''] for r in meta} and np.isfinite(values).all()
    sample_index = {s:i for i,s in enumerate(samples)}
    group = {s:[sample_index[r['']] for r in meta if r['species'] == s] for s in species}
    medians = np.column_stack([np.median(values[:,group[s]], axis=1) for s in species])
    atac = {g:dict(zip(species, medians[i])) for i,g in enumerate(genes)}
    write(OUT / 'species_median_tmm_logcpm.tsv',
          [dict(gene=g, **atac[g]) for g in genes], ['gene', *species])

    assert not query('SIRT6')['stale'], 'Liver RNA index requires rebuilding.'
    aliases = {'Ellobius lutescens': 'Bramus_lutescens'}
    records = defaultdict(list)
    seen = set()
    rna_samples = defaultdict(set)
    taxonomy_notes = set()
    with sqlite3.connect(f'file:{WEB}/data/liver_omics.sqlite?mode=ro', uri=True) as con:
        for gene, scientific, sample, value, common_field in con.execute(
                'SELECT symbol,species,sample,log,common FROM measurements WHERE assay="rna" AND sample!="BMR_Turk_Lu"'):
            sid = aliases.get(scientific, scientific.replace(' ', '_'))
            if sid not in species or gene.upper() not in atac:
                continue
            key = gene.upper(), sid, sample
            assert key not in seen, key
            seen.add(key)
            if value is None or not math.isfinite(value):
                continue
            records[(gene.upper(), sid)].append(value)
            rna_samples[sid].add(sample)
            if scientific == 'Acomys cahirinus' and common_field == 'Acomys russatus':
                taxonomy_notes.add('Acomys RNA species field is cahirinus but common field is russatus; an ASM-excluded sensitivity is included.')
        sources = json.loads(con.execute("SELECT value FROM metadata WHERE key='sources'").fetchone()[0])
    rna = defaultdict(dict)
    for (gene,sid), vals in records.items():
        rna[gene][sid] = float(np.median(vals))
    write(OUT / 'species_median_rna_source_logtpm.tsv',
          [dict(gene=g, **rna[g]) for g in sorted(rna)], ['gene', *species])
    inclusion = []
    for s in species:
        trait = traits[s.replace('_', ' ')]
        inclusion.append(dict(species=s, common_name=common[s], ATAC_samples=len(group[s]),
                              RNA_samples=len(rna_samples[s]), MLS_years=trait.get('years'),
                              MLS_lower_bound=trait.get('lower_bound_years'), body_mass_g=trait.get('body_mass_g'),
                              included_ATAC_MLS=s in eligible, included_matched_pool=s in eligible and bool(rna_samples[s]),
                              trait_source=trait.get('source'), trait_note=trait.get('note')))
    write(OUT / 'species_inclusion.tsv', inclusion, list(inclusion[0]))
    summaries = []
    all_atac_summaries = []
    validations = []
    plotted = {}
    for scale in ['MLS', 'log10_MLS']:
        xvals = {s:float(traits[s.replace('_', ' ')]['years']) for s in eligible}
        if scale == 'log10_MLS':
            xvals = {s:math.log10(v) for s,v in xvals.items()}
        all_atac = []
        x = np.array([xvals[s] for s in eligible])
        for gene in genes:
            result = fit(x, np.array([atac[gene][s] for s in eligible]))
            if result:
                all_atac.append(dict(gene=gene, n_species=len(eligible), **result, species=','.join(eligible)))
        add_q(all_atac, [''])
        all_atac_summaries.append(dict(MLS_scale=scale, tested_genes=len(all_atac), n_species=len(eligible),
                                       positive=sum(r['r'] > 0 for r in all_atac),
                                       negative=sum(r['r'] < 0 for r in all_atac),
                                       positive_p05=sum(r['r'] > 0 and r['p'] < .05 for r in all_atac),
                                       negative_p05=sum(r['r'] < 0 and r['p'] < .05 for r in all_atac),
                                       positive_FDR05=sum(r['r'] > 0 and r['q'] < .05 for r in all_atac),
                                       negative_FDR05=sum(r['r'] < 0 and r['q'] < .05 for r in all_atac)))
        write(OUT / scale / 'all_available_TMM_correlations.tsv', sorted(all_atac,key=lambda r:r['p']),
              ['gene','n_species','r','p','q','df','species'])
        for variant, mass_correct, exclude_asm in [('unadjusted',False,False),
                                                   ('body_mass_adjusted',True,False),
                                                   ('exclude_ASM',False,True)]:
            results = []
            skipped = []
            for gene in genes:
                ids = sorted(set(eligible) & set(rna[gene]))
                if exclude_asm:
                    ids = [s for s in ids if s != 'Acomys_cahirinus']
                x = np.array([xvals[s] for s in ids])
                y1 = np.array([atac[gene][s] for s in ids])
                y2 = np.array([rna[gene][s] for s in ids])
                mass = np.array([math.log10(traits[s.replace('_',' ')]['body_mass_g']) for s in ids]) if mass_correct else None
                a, b = fit(x,y1,mass), fit(x,y2,mass)
                if a is None or b is None:
                    skipped.append(dict(gene=gene,n_species=len(ids),reason='fewer than 5 species or constant/rank-deficient fit'))
                    continue
                direction = 'positive_both' if a['r'] > 0 and b['r'] > 0 else 'negative_both' if a['r'] < 0 and b['r'] < 0 else 'opposite'
                results.append(dict(gene=gene,n_species=len(ids),direction=direction,
                                    TMM_r=a['r'],TMM_p=a['p'],RNA_r=b['r'],RNA_p=b['p'],df=a['df'],species=','.join(ids)))
            add_q(results, ['TMM_', 'RNA_'])
            results.sort(key=lambda r:max(r['TMM_p'],r['RNA_p']))
            fields=['gene','n_species','direction','TMM_r','TMM_p','TMM_q','RNA_r','RNA_p','RNA_q','df','species']
            folder = OUT / scale / variant
            write(folder / 'all_gene_comparisons.tsv',results,fields)
            write(folder / 'excluded_genes.tsv',skipped,['gene','n_species','reason'])
            summary = dict(MLS_scale=scale,variant=variant,tested_genes=len(results),
                           min_species=min(r['n_species'] for r in results),max_species=max(r['n_species'] for r in results),
                           median_species=float(np.median([r['n_species'] for r in results])),
                           positive_both=sum(r['direction']=='positive_both' for r in results),
                           negative_both=sum(r['direction']=='negative_both' for r in results),
                           opposite=sum(r['direction']=='opposite' for r in results),
                           correlation_profile_pearson=float(pearsonr([r['TMM_r'] for r in results],[r['RNA_r'] for r in results]).statistic),
                           correlation_profile_spearman=float(spearmanr([r['TMM_r'] for r in results],[r['RNA_r'] for r in results]).statistic))
            for threshold,col in [('all',None),('nominal_p05','p'),('FDR05','q')]:
                selection = results if col is None else [r for r in results if r['TMM_'+col] < .05 and r['RNA_'+col] < .05]
                for direction in ['positive_both','negative_both','opposite']:
                    selected = [r for r in selection if r['direction']==direction]
                    write(folder / f'{threshold}_{direction}.tsv',selected,fields)
                    if col:
                        summary[f'{threshold}_{direction}'] = len(selected)
                if col:
                    summary[f'{threshold}_TMM'] = sum(r['TMM_'+col] < .05 for r in results)
                    summary[f'{threshold}_RNA'] = sum(r['RNA_'+col] < .05 for r in results)
            assert summary['positive_both'] + summary['negative_both'] + summary['opposite'] == len(results)
            for row in [results[0],results[len(results)//2],results[-1]]:
                ids=row['species'].split(',')
                x=np.array([xvals[s] for s in ids])
                design=np.column_stack([np.ones(len(ids)),x]+([
                    np.array([math.log10(traits[s.replace('_',' ')]['body_mass_g']) for s in ids])
                ] if mass_correct else []))
                for prefix,data in [('TMM_',atac),('RNA_',rna)]:
                    y=np.array([data[row['gene']][s] for s in ids])
                    beta=np.linalg.lstsq(design,y,rcond=None)[0]
                    se=math.sqrt(np.sum((y-design@beta)**2)/(len(ids)-design.shape[1])*np.linalg.inv(design.T@design)[1,1])
                    independent_p=float(2*t.sf(abs(beta[1]/se),len(ids)-design.shape[1]))
                    assert math.isclose(independent_p,row[prefix+'p'],rel_tol=1e-8,abs_tol=1e-12)
                    validations.append(dict(scale=scale,variant=variant,gene=row['gene'],assay=prefix[:-1],
                                            p_correlation=row[prefix+'p'],p_OLS=independent_p))
            summaries.append(summary)
            if variant=='unadjusted':
                plotted[scale]=(results,summary)
    write(OUT/'comparison_summary.tsv',summaries,list(summaries[0]))
    write(OUT/'all_available_TMM_summary.tsv',all_atac_summaries,list(all_atac_summaries[0]))
    write(OUT/'OLS_validation.tsv',validations,list(validations[0]))

    plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(12,5.4))
    for ax,scale in zip(axes,['MLS','log10_MLS']):
        rows,summary=plotted[scale]
        ax.scatter([r['TMM_r'] for r in rows],[r['RNA_r'] for r in rows],s=9,c='#8294a5',alpha=.23,edgecolors='none',rasterized=True)
        for col,color,label in [('p','#e07a22','Both p < 0.05'),('q','#005b96','Both FDR < 0.05')]:
            points=[r for r in rows if r['TMM_'+col]<.05 and r['RNA_'+col]<.05]
            if points:
                ax.scatter([r['TMM_r'] for r in points],[r['RNA_r'] for r in points],s=21,c=color,alpha=.9,label=label,edgecolors='white',linewidths=.2)
        ax.axhline(0,color='#777777',lw=.8,ls='--');ax.axvline(0,color='#777777',lw=.8,ls='--')
        ax.set(xlim=(-1,1),ylim=(-1,1),aspect='equal',xlabel='TMM ATAC–lifespan Pearson r',ylabel='Transcriptome–lifespan Pearson r')
        ax.set_title(('MLS (years)' if scale=='MLS' else 'log₁₀(MLS)')+f" · {len(rows):,} matched genes",fontsize=12)
        ax.text(.03,.96,f"Positive in both: {summary['positive_both']:,}\nBoth positive, p < 0.05: {summary['nominal_p05_positive_both']}\nBoth positive, FDR < 0.05: {summary['FDR05_positive_both']}",
                transform=ax.transAxes,va='top',fontsize=9)
        ax.legend(loc='lower right',fontsize=8,frameon=False)
    fig.tight_layout()
    fig.savefig(OUT/'TMM_transcriptome_MLS_comparison.png',dpi=220,bbox_inches='tight',facecolor='white')
    fig.savefig(OUT/'TMM_transcriptome_MLS_comparison.pdf',bbox_inches='tight',facecolor='white')
    plt.close(fig)
    provenance=dict(source_sha256=hashlib.sha256((HERE/'logcpm_tmm_core.tsv').read_bytes()).hexdigest(),
                    RNA_sources=sources,ATAC_genes=len(genes),ATAC_samples=len(samples),ATAC_species=len(species),
                    eligible_ATAC_MLS_species=len(eligible),matched_RNA_species_pool=sum(r['included_matched_pool'] for r in inclusion),
                    taxonomy_notes=sorted(taxonomy_notes),summaries=summaries,all_available_TMM=all_atac_summaries)
    (OUT/'summary.json').write_text(json.dumps(provenance,indent=2)+'\n')
    (OUT/'METHODS.txt').write_text(
        'Input ATAC: ../logcpm_tmm_core.tsv, 6,459 genes and 75 samples from 32 species. Values are supplied TMM log-CPM.\n'
        'Aggregate the already-logged values by median across samples within each species. No new log transform, TMM normalization, or batch removal. This is not a reconstructed raw-CPM median. Original TMM log pseudocount and TSS window construction were not supplied.\n'
        'RNA: current liver_omics.sqlite species median of supplied log.tpm, retaining original liver sources and mouse supplement. Exclude BMR_Turk_Lu as in the browser. No additional log transform.\n'
        'MLS: current project species_mls.json. Eonycteris spelaea has only a >20-year bound and is excluded from numeric correlations. Human, macaque and bowhead are absent from the TMM input. The BMR MLS remains the existing Spalax-complex proxy.\n'
        'For each gene, TMM and RNA correlations use EXACTLY the same eligible species. At least five species required. Zero-valued logged measurements are retained; absent RNA measurements are missing, not zero. Genes are matched by uppercase symbol, not verified one-to-one orthology.\n'
        'Ellobius lutescens RNA maps to Bramus lutescens. Acomys RNA has a conflicting common-name field; the primary analysis follows its species field and exclude_ASM provides sensitivity.\n'
        'Primary: ordinary Pearson correlation against raw MLS, plus log10 MLS. No body mass, gene-length or phylogenetic covariates. body_mass_adjusted is a separate partial-Pearson sensitivity using log10 grams for BOTH assays.\n'
        'Both correlation signs and two-sided nominal p<0.05 are reported. BH FDR is computed separately for each assay, scale and variant across the matched tested gene set. FDR-positive overlap requires q<0.05 in both assays AND positive r in both.\n'
        'All-available TMM correlations use 31 eligible species and their own 6,459-gene FDR family; use matched comparisons for RNA overlap, not the all-available TMM p/q.\n'
        'Correlation-profile Pearson/Spearman values describe agreement across genes; no inferential p is assigned because genes need not be independent. The comparison is of MLS associations, not direct ATAC-versus-RNA expression correlation.\n'
        'OLS_validation.tsv independently checks the correlation p-values against full linear-model coefficient tests.\n'
        'Reproduce from project root: /usr/bin/python3 tmm_pca/compare_tmm_mls_rna.py\n')
    print(json.dumps({'eligible_ATAC_species':len(eligible),'matched_RNA_species_pool':provenance['matched_RNA_species_pool'],
                      'summaries':summaries},indent=2))


if __name__=='__main__':
    main()
