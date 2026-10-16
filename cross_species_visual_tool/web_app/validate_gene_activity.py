"""Verify output sums independently against bigWig interval records."""
import csv
import gzip
import json
from pathlib import Path
import numpy as np
import pyBigWig
from calculate_gene_activity import ROOT

import argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,default=ROOT/'output/gene_activity_signac_regions')
out=parser.parse_args().output
manifest=json.loads((out/'manifest.json').read_text())
checked=0
for sample in manifest['samples']:
    path=Path(sample['path'])
    candidates=[]
    with gzip.open(out/f"{sample['species']}_gene_activity_long.tsv.gz",'rt') as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            if row['sample']==path.stem and row['status']=='ok' and float(row['signal_sum'])>0:
                candidates.append(row)
                if len(candidates)==5: break
    with pyBigWig.open(str(path)) as bw:
        for row in candidates:
            chrom=row['chrom']
            if chrom not in bw.chroms(): chrom=next(c for c in bw.chroms() if c.removeprefix('chr')==chrom.removeprefix('chr'))
            lo=max(0,int(row['region_start'])); hi=min(bw.chroms(chrom),int(row['region_end']))
            intervals=bw.intervals(chrom,lo,hi) or []
            total=sum((min(hi,b)-max(lo,a))*v for a,b,v in intervals if np.isfinite(v))
            assert np.isclose(total,float(row['signal_sum']),rtol=1e-6,atol=1e-6),(path,row['locus_id'],total,row['signal_sum'])
            checked+=1
print(f'Independent interval integration matched {checked} nonzero gene/sample scores across {len(manifest["samples"])} samples.')
