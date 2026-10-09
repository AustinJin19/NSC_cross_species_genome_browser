"""Windowed retrieval of NMR TE repeat annotations."""
from pathlib import Path
import pysam
DB=Path(__file__).resolve().parent/'data'/'nmr_TE.bed.gz'
CLASSES=[('LINE','#377eb8'),('SINE','#008b8b'),('LTR','#984ea3'),('DNA','#d97727'),('Retroposon','#66733b'),('RC','#ad5379')]
def query(chrom,lo,hi):
    regions=[];available=DB.exists() and Path(str(DB)+'.tbi').exists()
    if available and hi>0:
        with pysam.TabixFile(str(DB)) as tb:
            if chrom in tb.contigs:
                for line in tb.fetch(chrom,max(0,lo-1),hi):
                    r=line.split('\t');regions.append(dict(start=int(r[1])+1,end=int(r[2]),name=r[3],repeatClass=r[3].split('/')[0],score=r[4],strand=r[5]))
    return dict(available=available,regions=regions,classes=[dict(name=n,color=c) for n,c in CLASSES],assembly='mHetGlaV3')
