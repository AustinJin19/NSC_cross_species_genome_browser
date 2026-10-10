"""Build portable transcript indexes for the additional local ATAC assemblies."""
import gzip
import json
import re
from pathlib import Path
import pysam
import pyBigWig

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'cross_species_visual_tool/index'
SOURCES = {
 'Human': ('selected3', 'gencode.v50.primary_assembly.annotation.gtf.gz', 'hg38.fa.fai'),
 'Mouse_mm10': ('mouse_mm10', 'Mus_musculus_mm10.ucsc_ncbiRefSeq.gtf.gz', 'mm10.fa.fai'),
 'Macaque': ('macaque', 'Macaca_mulatta.gtf.gz', 'Macaca_mulatta.fa.fai'),
 'Rabbit': ('Rabbit', 'Sylvilagus_floridanus.gtf.gz', 'Sylvilagus_floridanus.fa.fai'),
 'ASM': ('Africa_spiny_mouse', 'Acomys_cahirinus.gz', 'Acomys_cahirinus.fa.fai'),
}

def build(key, folder, annotation, fai):
 directory = ROOT / 'reference_genomes' / folder
 sizes = {f[0]:int(f[1]) for line in (directory/fai).read_text().splitlines() if (f:=line.split('\t'))}
 tracks = sorted((ROOT / 'data/atac' / folder).glob('*.bw'))
 for track in tracks:
  with pyBigWig.open(str(track)) as bw:
   mismatches = [c for c,n in bw.chroms().items() if sizes.get(c)!=n]
   if mismatches: raise ValueError(f'{track.name}: FASTA/bigWig sequence mismatch: {mismatches[:5]}')
 with pyBigWig.open(str(tracks[0])) as bw: chroms = bw.chroms()
 aliases = {}
 # Read FASTA descriptions immediately preceding each indexed sequence offset.
 # This uses explicit scaffold names in the supplied FASTA, never inferred numbering.
 with (directory / fai.removesuffix('.fai')).open('rb') as fasta:
  for line in (directory/fai).read_text().splitlines():
   f=line.split('\t');offset=int(f[2]);fasta.seek(max(0,offset-4096))
   header=fasta.read(min(offset,4096)).decode().split('>')[-1]
   for alias in re.findall(r'\bscaffold_\d+\b',header): aliases[alias]=f[0]
 transcripts = {}; skipped = 0
 with gzip.open(directory/annotation,'rt') as handle:
  for line in handle:
   if line.startswith('#'): continue
   f=line.rstrip().split('\t')
   if len(f)!=9 or f[2] not in ('transcript','exon'): continue
   chrom,_,feature,start,end,_,strand,_,attrs=f
   chrom=aliases.get(chrom,chrom)
   a,b=int(start)-1,int(end)
   if chrom not in chroms or a<0 or b>chroms[chrom] or b<=a:
    skipped+=1; continue
   fields=dict(re.findall(r'(\w+) "([^"]*)"',attrs))
   tx=fields.get('transcript_id')
   if not tx: continue
   # bed2gtf encodes the original gene symbol in transcript_id: accession#symbol#id.
   symbol=fields.get('gene_name') or fields.get('gene') or (tx.split('#')[1] if '#' in tx else fields.get('gene_id'))
   if not symbol: continue
   record=transcripts.setdefault((chrom,strand,tx),dict(symbol=symbol,start=a,end=b,exons=set()))
   record['start']=min(record['start'],a);record['end']=max(record['end'],b)
   if feature=='exon':record['exons'].add((a,b))
 rows=[]
 for (chrom,strand,tx),r in transcripts.items():
  exons=sorted(r['exons'])
  if not exons: continue
  rows.append((chrom,r['start'],r['end'],tx,r['symbol'],strand,','.join(str(x[0]) for x in exons)+',',','.join(str(x[1]) for x in exons)+','))
 rows.sort(key=lambda r:(r[0],r[1],r[2],r[3]))
 plain=OUT/f'{key}_models.tsv'
 plain.write_text(''.join('\t'.join(map(str,r))+'\n' for r in rows))
 target=str(plain)+'.gz'
 pysam.tabix_compress(str(plain),target,force=True)
 pysam.tabix_index(target,preset='bed',force=True)
 plain.unlink()
 info=dict(source=str(directory/annotation),transcripts=len(rows),genes=len({r[4].upper() for r in rows}),unsupported_or_invalid_features=skipped,samples=[p.name for p in tracks])
 Path(target+'.json').write_text(json.dumps(info,indent=2)+'\n')
 if key == 'Human':
  from build_human_anchors import build as build_anchors
  build_anchors()
 print(key,info,flush=True)

if __name__=='__main__':
 for key,args in SOURCES.items():build(key,*args)
