"""Choose human TSS anchors from supplied GENCODE transcript tags."""
import gzip
import json
import re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
DEST=ROOT/'cross_species_visual_tool/index/Human_anchors.json'

def build():
 candidates={}
 with gzip.open(ROOT/'reference_genomes/selected3/gencode.v50.primary_assembly.annotation.gtf.gz','rt') as handle:
  for line in handle:
   if line.startswith('#'):continue
   f=line.rstrip().split('\t')
   if f[2]!='transcript':continue
   attrs=dict(re.findall(r'(\w+) "([^"]*)"',f[8]))
   tags=set(re.findall(r'tag "([^"]*)"',f[8]))
   symbol=attrs.get('gene_name',attrs['gene_id']).upper()
   method=next((tag for tag in ['MANE_Select','Ensembl_canonical','appris_principal_1','basic'] if tag in tags),'longest transcript')
   rank=['MANE_Select','Ensembl_canonical','appris_principal_1','basic','longest transcript'].index(method)
   start,end=int(f[3]),int(f[4]);tx=attrs['transcript_id']
   key=(symbol,f[0],f[6]);score=(rank,-(end-start),tx)
   if key not in candidates or score<candidates[key][0]:
    candidates[key]=(score,dict(symbol=symbol,chrom=f[0],strand=f[6],start=start,end=end,tss=end if f[6]=='-' else start,transcript=tx,tssMethod=method))
 DEST.write_text(json.dumps([v[1] for v in candidates.values()],indent=2)+'\n')
 print(f'Wrote {len(candidates)} human anchors')
if __name__=='__main__':build()
