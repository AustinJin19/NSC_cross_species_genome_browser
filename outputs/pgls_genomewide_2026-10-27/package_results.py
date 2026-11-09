"""Summarize PGLS hits, OLS retention, nominal hits, and joint multiplicity."""
import csv,gzip,hashlib,json,zipfile
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.stats import false_discovery_control
from analyze import METHODS,FIELDS,write
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1];APP=ROOT/'cross_species_visual_tool/web_app'
LABELS={'body':'Gene body + 2 kb','tss2kb':'TSS ±2 kb','tss5kb':'TSS ±5 kb','tss100kb':'TSS ±100 kb','tss200kb':'TSS ±200 kb','tss500kb':'TSS ±500 kb'}

def main():
    rows=[];summaries=[];method_rows={}
    for method in METHODS:
        summary=json.loads((OUT/'results'/method/'summary.json').read_text())
        with gzip.open(OUT/'results'/method/'all_tested.tsv.gz','rt') as f:records=list(csv.DictReader(f,delimiter='\t'))
        for r in records:
            for k in ['MLS_beta','MLS_SE','p_value','BH_FDR','lambda','GLS_R2','model_intercept','mass_beta','OLS_beta','OLS_p_value','OLS_BH_FDR','OLS_partial_r']:r[k]=float(r[k])
            for k in ['n_species','df']:r[k]=int(r[k])
            r['length_beta']=float(r['length_beta']) if r['length_beta'] else None
            r['retained_OLS_FDR_hit']=r['retained_OLS_FDR_hit']=='True'
        summary['family']='Regulatory peak count' if method.startswith('peaks') else 'Gene activity'
        summary['region']=f'TSS ±{method[5:-2]} kb' if method.startswith('peaks') else LABELS[method]
        summary['retained_nominal_same_direction']=sum(r['p_value']<.05 and r['OLS_p_value']<.05 and r['MLS_beta']*r['OLS_beta']>0 for r in records)
        summary['OLS_nominal_positive']=sum(r['OLS_p_value']<.05 and r['OLS_beta']>0 for r in records)
        summary['OLS_nominal_negative']=sum(r['OLS_p_value']<.05 and r['OLS_beta']<0 for r in records)
        # Report a clearly labeled sensitivity analysis with BH recomputed for n>=10.
        sufficient=[r for r in records if r['n_species']>=10]
        q10=false_discovery_control([r['p_value'] for r in sufficient])
        summary['n10_tested']=len(sufficient);summary['n10_FDR_positive']=int(sum(q<.05 and r['MLS_beta']>0 for q,r in zip(q10,sufficient)));summary['n10_FDR_negative']=int(sum(q<.05 and r['MLS_beta']<0 for q,r in zip(q10,sufficient)))
        method_rows[method]=records;rows.extend(records);summaries.append(summary)
    qjoint=false_discovery_control([r['p_value'] for r in rows])
    for r,q in zip(rows,qjoint):r['BH_FDR_all_12_methods']=float(q)
    for family in ['gene','regulatory']:
        selected=[r for r in rows if r['method'].startswith('peaks')==(family=='regulatory')]
        qfamily=false_discovery_control([r['p_value'] for r in selected])
        for r,q in zip(selected,qfamily):r['BH_FDR_within_family']=float(q)
    final_fields=FIELDS+['BH_FDR_within_family','BH_FDR_all_12_methods']
    hits=[r for r in rows if r['BH_FDR']<.05]
    nominal=[r for r in rows if r['p_value']<.05]
    for summary in summaries:
        rs=method_rows[summary['method']]
        summary['joint_FDR_total']=sum(r['BH_FDR_all_12_methods']<.05 for r in rs)
        summary['family_FDR_total']=sum(r['BH_FDR_within_family']<.05 for r in rs)
    write(OUT/'window_summary.tsv',summaries,list(summaries[0]))
    write(OUT/'PGLS_FDR05_genes_all_windows.tsv',hits,final_fields)
    write(OUT/'PGLS_nominal_p05_genes_all_windows.tsv.gz',nominal,final_fields)
    write(OUT/'all_tested_PGLS_and_OLS.tsv.gz',rows,final_fields)
    unique=defaultdict(list)
    for r in hits:unique[r['gene']].append(r)
    unique_rows=[dict(gene=g,significant_methods='; '.join(r['method'] for r in rs),
                      directions='; '.join(r['direction'] for r in rs),minimum_FDR=min(r['BH_FDR'] for r in rs),
                      also_significant_with_OLS=any(r['retained_OLS_FDR_hit'] for r in rs)) for g,rs in sorted(unique.items())]
    write(OUT/'PGLS_FDR05_unique_genes.tsv',unique_rows,['gene','significant_methods','directions','minimum_FDR','also_significant_with_OLS'])
    sirt6=[r for r in rows if r['gene']=='SIRT6'];write(OUT/'SIRT6_all_windows.tsv',sirt6,final_fields)
    report=dict(summary=summaries,FDR_hits=hits,unique_FDR_genes=unique_rows,SIRT6=sirt6,
                total_tests=len(rows),unique_significant_genes=len(unique),
                distinct_gene_activity_hits=len({r['gene'] for r in hits if not r['method'].startswith('peaks')}),
                distinct_regulatory_hits=len({r['gene'] for r in hits if r['method'].startswith('peaks')}),
                joint_FDR_gene_method_hits=sum(r['BH_FDR_all_12_methods']<.05 for r in rows))
    (OUT/'report.json').write_text(json.dumps(report,indent=2))
    metadata=dict(date='2026-10-27',model='log10(species median activity) ~ log10(MLS) + log10(body mass); body additionally log10(gene length kb)',
        phylogeny='TimeTree 6 shared root-to-MRCA branch covariance; original full-tree root retained when subsetting',
        lambda_estimation='Per-gene Gaussian profile maximum likelihood in [0,1]; same model as browser NSCPGLS',
        test='Two-sided coefficient t test conditional on estimated lambda (plug-in approximation); df=n-number_of_coefficients',
        primary_significance='BH FDR < 0.05 separately within each method, across all tested genes',
        additional_corrections='BH within each six-method family; BH pooled across all 12 gene-method scans',
        nominal_significance='Unadjusted p<0.05, reported separately from FDR',
        minimum_species=5,available_region_fraction_minimum=.5,maximum_species=30,
        species_excluded=['Human','Macaque','Balaena_mysticetus','missing/lower-bound MLS or unavailable body mass'],
        sample_QC='Current pre-TSS-filter website sample set; TSS<4 samples retained; legacy NMR1 exclusion preserved',
        selection='Highest species-median locus chosen separately for each method before coverage/log filters; no lower-locus fallback',
        comparison='OLS refitted on exactly the same species/genes and covariates within each method; different windows can have different eligible species',
        measurements='Gene activity: unweighted integrated bigWig CPM·bp. Regulatory activity: overlapping peak interval count; not CPM-in-peaks.',
        limitations=['Symbol matches are not verified one-to-one orthologs.','No batch, library-quality or peak-number normalization correction was added.','Lambda point estimates only; no confidence intervals or bootstrap.','PGLS adjusts residual relatedness, not genome/annotation quality.','The 1000 kb regulatory scan is the previous offline window; not an additional website option.'],
        files={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [APP/'data/phylogenetic_covariance.json',APP/'data/species_mls.json',APP/'data/sample_qc_exclusions.json',OUT/'pgls_batch.py',OUT/'analyze.py',OUT/'prepare_inputs.py']})
    (OUT/'manifest.json').write_text(json.dumps(metadata,indent=2))
    archive=OUT/'PGLS_MLS_all_windows_results.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for p in [OUT/'window_summary.tsv',OUT/'PGLS_FDR05_genes_all_windows.tsv',OUT/'PGLS_FDR05_unique_genes.tsv',OUT/'PGLS_nominal_p05_genes_all_windows.tsv.gz',OUT/'all_tested_PGLS_and_OLS.tsv.gz',OUT/'SIRT6_all_windows.tsv',OUT/'manifest.json',OUT/'input_validation.json',OUT/'species_inventory.json',OUT/'OLS_reproduction.json',OUT/'browser_validation.json',OUT/'R_caper_validation.tsv']:
            if p.exists():z.write(p,p.name)
        for p in (OUT/'results').rglob('*_genes.txt'):z.write(p,p.relative_to(OUT))
    print(json.dumps({k:v for k,v in report.items() if k not in ['FDR_hits','SIRT6']},indent=2))

if __name__=='__main__':main()
