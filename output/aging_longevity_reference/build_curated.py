from pathlib import Path
import csv
root=Path(__file__).resolve().parent
rows=[]
def add(genes,paper,evidence,direction,context,notes,url):
    for g in genes.split():rows.append(dict(gene_symbol=g,paper=paper,evidence_type=evidence,direction=direction,context=context,interpretation=notes,source=url))
nature='https://www.nature.com/articles/s41586-026-10542-3'
cell='https://doi.org/10.1016/j.cell.2023.05.002'
gwas='https://www.nature.com/articles/s41467-019-11558-2'
elife='https://doi.org/10.7554/eLife.39856'
add('GPNMB VSIG4 CDKN1A EDA2R','Tyshkovskiy 2026','Expression association','increases with age','Cross-species aging','Aging marker; not proof that higher expression shortens lifespan',nature)
add('NREP COL1A1 COL3A1','Tyshkovskiy 2026','Expression association','decreases with age','Cross-species aging','Aging marker; not proof that higher expression extends lifespan',nature)
add('CDKN1A LGALS3','Tyshkovskiy 2026','Protein association','associated with mortality/multimorbidity','UK Biobank','Association, not causal lifespan manipulation; see source for outcome-specific effects',nature)
add('IGF1 CUL4B','Tyshkovskiy 2023','Expression association','lower expression in longer-lived species','Comparative mammalian longevity','Does not establish the direction of ATAC accessibility',cell)
add('RPL28 MPG','Tyshkovskiy 2023','Expression association','higher expression in longer-lived species','Comparative mammalian longevity','Does not establish the direction of ATAC accessibility',cell)
add('APOE','Deelen 2019','Longevity GWAS','epsilon2 favorable; epsilon4 unfavorable','Human survival percentile','Allele-specific association, not gene expression direction',gwas)
add('FOXO3 CDKN2A CDKN2B','Deelen 2019','Prior longevity locus replication','not assigned','Human longevity','Locus-level evidence; genes listed separately do not imply independent causal signals',gwas)
add('GPR78','Deelen 2019','GWAS nearby-gene candidate','rs7676745 A associated with lower longevity odds','Human 90th-percentile longevity','Discovery association; not tested in de novo replication owing to assay availability; gene not proven causal',gwas)
add('KANSL1 CRHR1 ARL17A LRRC37A2','Deelen 2019','Genetically predicted expression association','higher predicted expression associated with longevity','Human tissue-specific transcriptome association','Not measured liver expression; correlated loci may share a signal',gwas)
add('ANKRD31 BLOC1S1','Deelen 2019','Genetically predicted expression association','lower predicted expression associated with longevity','Human tissue-specific transcriptome association','Not measured liver expression; not proof of causal gene',gwas)
add('ABO ZC3HC1 IGF2R','Timmers 2019','Parental lifespan GWAS locus','not assigned','Human parental lifespan','New replicated loci labeled by nearby gene; causal gene and expression direction unresolved',elife)
add('CDKN2B-AS1 ATXN2 BRAP FURIN FES ZW10 PSORS1C3 EBF1 FOXO3','Timmers 2019','Parental lifespan GWAS locus','not assigned','Human parental lifespan','Validated loci; ATXN2/BRAP and FURIN/FES are shared loci; genes are not independent causal discoveries',elife)
with (root/'curated_gene_reference.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
(root/'curated_gene_symbols.txt').write_text('\n'.join(sorted({r['gene_symbol'] for r in rows}))+'\n')
print(len(rows),'evidence rows;',len({r['gene_symbol'] for r in rows}),'unique genes')
