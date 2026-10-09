"""Build a sorted tabix index of TE annotations on browser-supported NMR sequences."""
import gzip,json,subprocess,tempfile
from pathlib import Path
from collections import Counter
import pyBigWig,pysam
ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'mHetGlaV3.primary.filteredRepeats.bed.gz'
OUT=Path(__file__).resolve().parent/'data'/'nmr_TE.bed.gz'
CLASSES={'LINE','SINE','LTR','DNA','Retroposon','RC'}
def build():
    counts=Counter()
    with pyBigWig.open(str(ROOT/'cross_species_visual_tool/mole_rat_peaks/NMR3Liver.cpm.bw')) as bw:chroms=bw.chroms()
    with tempfile.TemporaryDirectory(dir=OUT.parent) as tmp:
        plain=Path(tmp)/'regions.bed';sorted_path=Path(tmp)/'sorted.bed'
        with gzip.open(SOURCE,'rt') as src,plain.open('w') as dst:
            for line in src:
                if not line.strip() or line.startswith('#'):continue
                r=line.split();counts['input_rows']+=1
                if r[3].split('/')[0] not in CLASSES:counts['non_TE_rows']+=1;continue
                c,a,b=r[0],int(r[1]),int(r[2])
                if c not in chroms:counts['TE_unsupported_sequence']+=1;continue
                if a<0 or b<=a or b>chroms[c]:counts['TE_invalid_coordinates']+=1;continue
                dst.write('\t'.join(r[:6])+'\n');counts['indexed_TE_rows']+=1
        with sorted_path.open('w') as dest:subprocess.run(['sort','-k1,1','-k2,2n','-k3,3n',str(plain)],stdout=dest,env={'LC_ALL':'C','PATH':'/usr/bin:/bin'},check=True)
        pysam.tabix_compress(str(sorted_path),str(OUT),force=True)
        pysam.tabix_index(str(OUT),preset='bed',force=True)
    meta=dict(counts,source=SOURCE.name,size=SOURCE.stat().st_size,mtime_ns=SOURCE.stat().st_mtime_ns,classes=sorted(CLASSES))
    OUT.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n');print(meta,flush=True)
if __name__=='__main__':build()
