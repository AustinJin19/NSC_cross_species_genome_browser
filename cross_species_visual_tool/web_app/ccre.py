"""ENCODE4 biosample-agnostic hg38 cCREs; indexed local bigBed queries."""
from pathlib import Path
import pyBigWig
DATA=Path(__file__).resolve().parents[2]/'ENCODE4/encodeCcreRegistry.hg38.bb'
def query(chrom,lo,hi):
    result=dict(available=DATA.exists(),assembly='hg38',source='ENCODE4 cCRE v4',regions=[])
    if not result['available']:return result
    with pyBigWig.open(str(DATA)) as bed:
        size=bed.chroms(chrom)
        if size is None or hi<1 or lo>size:return result
        for start,end,text in bed.entries(chrom,max(0,lo-1),min(size,hi)) or []:
            fields=text.split('\t')
            result['regions'].append(dict(start=start+1,end=end,name=fields[0],category=fields[6],color='rgb('+fields[5]+')'))
    return result
