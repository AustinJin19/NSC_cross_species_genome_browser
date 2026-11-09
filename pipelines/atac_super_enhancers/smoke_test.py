#!/usr/bin/env python3
"""Exercise real ROSE stitching, BAM scoring, R cutoff and summaries on synthetic data."""
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import pysam
import workflow as w


def main():
    if not shutil.which('Rscript'):
        raise SystemExit('Rscript is needed for the smoke test')
    with tempfile.TemporaryDirectory(prefix='rose_smoke_',dir='/tmp') as tmp:
        root = Path(tmp)
        env = os.environ.copy()
        adapter = not shutil.which('samtools')
        if adapter:
            # Exercise actual SAMtools C code distributed with pysam, not fabricated output.
            cli = root/'bin/samtools'
            cli.parent.mkdir()
            cli.write_text(f'#!{sys.executable}\nimport pysam,sys\nif sys.argv[1]=="--version":\n print("samtools "+pysam.__samtools_version__+" (pysam adapter)"); sys.exit(0)\nout=getattr(pysam.samtools,sys.argv[1])(*sys.argv[2:])\nif isinstance(out,str): sys.stdout.write(out)\nelif isinstance(out,bytes): sys.stdout.buffer.write(out)\n')
            cli.chmod(0o755)
            env['PATH'] = str(cli.parent)+os.pathsep+env['PATH']
        sizes = {'NC_123.1':1000000,'chr2':1000000}
        (root/'chrom.sizes').write_text(''.join(f'{c}\t{n}\n' for c,n in sizes.items()))
        with gzip.open(root/'models.gz','wt') as f:
            for i,c in enumerate(sizes):
                f.write(f'{c}\t800000\t810000\ttx{i}\tGENE{i}\t+\t800000,\t810000,\n')
        records,peaks = [],[]
        for i in range(20):
            tid = 0 if i < 10 else 1
            c = list(sizes)[tid]
            start = 2000+(i%10)*30000
            peaks.append((c,start,start+1000))
            copies = (5+i*2) if i < 17 else [300,900,2200][i-17]
            for j in range(copies):
                for mate,flag,offset,nextoffset,tlen in [(1,99,100,250,200),(2,147,250,100,-200)]:
                    read = pysam.AlignedSegment()
                    read.query_name = f'p{i}_{j}'
                    read.query_sequence = 'A'*50
                    read.query_qualities = pysam.qualitystring_to_array('I'*50)
                    read.flag = flag
                    read.reference_id = tid
                    read.reference_start = start+offset+j%40
                    read.mapping_quality = 60
                    read.cigarstring = '50M'
                    read.next_reference_id = tid
                    read.next_reference_start = start+nextoffset+j%40
                    read.template_length = tlen
                    records.append(read)
        records.sort(key=lambda r:(r.reference_id,r.reference_start))
        header = dict(HD=dict(VN='1.6',SO='coordinate'),SQ=[dict(SN=c,LN=n) for c,n in sizes.items()])
        with pysam.AlignmentFile(root/'input.bam','wb',header=header) as out:
            for r in records:
                out.write(r)
        pysam.index(str(root/'input.bam'))
        (root/'peaks.bed').write_text(''.join(f'{c}\t{a}\t{b}\n' for c,a,b in peaks))
        rows = [dict(species='Test',common_name='Synthetic test',scientific_name='Synthetic test',
                     sample=s,bam='input.bam',peaks='peaks.bed',models='models.gz',chrom_sizes='chrom.sizes',
                     exclude_bed='',exclude_chroms='',tss_enrichment='10',qc_flag='') for s in ('test1','test2')]
        w.write_table(root/'samples.tsv',rows,w.FIELDS)
        command = [sys.executable,str(w.HERE/'workflow.py')]
        try:
            for i in range(2):
                subprocess.run(command+['run-sample','--manifest',str(root/'samples.tsv'),'--project-root',str(root),
                               '--output',str(root/'results'),'--index',str(i)],env=env,check=True)
            subprocess.run(command+['summarize','--manifest',str(root/'samples.tsv'),'--output',str(root/'results')],env=env,check=True)
            result = json.loads((root/'results/Test/test1/summary.json').read_text())
            assert result['candidate_se_count']>0
            table = w.read_table(root/'results/Test/test1/all_regions.tsv')
            assert {r['chrom'] for r in table} == set(sizes)
            summary = w.read_table(root/'results/species_summary/species_counts.tsv')[0]
            assert int(summary['shared_ge2_regions'])>0
            assert summary['shared_ge2_bp'] == summary['union_bp']
            # Identical inputs can be resumed without rerunning computation.
            subprocess.run(command+['run-sample','--manifest',str(root/'samples.tsv'),'--project-root',str(root),
                           '--output',str(root/'results'),'--index','0'],env=env,check=True)
            print(json.dumps(dict(status='passed',candidate_se_count=result['candidate_se_count'],
                                  samtools_cli='pysam SAMtools adapter' if adapter else 'samtools',
                                  native_and_chr_prefixed_contigs=True,replicate_summary=True,resume=True),indent=2))
        except Exception:
            for log in (root/'results').glob('*/*/run.log'):
                print(str(log)+'\n'+log.read_text()[-8000:],file=sys.stderr)
            raise


if __name__ == '__main__':
    main()
