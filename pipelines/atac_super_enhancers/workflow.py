#!/usr/bin/env python3
"""Portable, per-sample ATAC super-enhancer calling with pinned ROSE algorithms."""
import argparse
import ast
import bisect
import collections
import contextlib
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
VENDOR = HERE / 'vendor' / 'ROSE'
FIELDS = ['species', 'common_name', 'scientific_name', 'sample', 'bam', 'peaks',
          'models', 'chrom_sizes', 'exclude_bed', 'exclude_chroms',
          'tss_enrichment', 'qc_flag']


def read_table(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_table(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def dump(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def text_open(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def resolve(path, root):
    p = Path(path)
    return p if p.is_absolute() else root / p


def identifier(s):
    if not re.fullmatch(r'[A-Za-z0-9_-]+', s):
        raise ValueError('Unsafe species/sample identifier: ' + s)
    return s


def manifest(path):
    rows = read_table(path)
    keys = []
    for r in rows:
        keys.append((identifier(r['species']), identifier(r['sample'])))
    if not rows or len(keys) != len(set(keys)):
        raise ValueError('Manifest is empty or has duplicate species/sample rows')
    return rows


def export_manifest(args):
    """Read the configured browser catalog without importing/running its server."""
    import pyBigWig
    root = args.project_root.resolve()
    web = root / 'cross_species_visual_tool/web_app'
    literals = {}
    for node in ast.parse((web / 'server.py').read_text()).body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id in ('SPECS', 'MODEL_INDEX'):
                    literals[target.id] = ast.literal_eval(node.value)
    specs = literals['SPECS'] + json.loads((web / 'data/additional_species.json').read_text())
    policy = json.loads((web / 'data/sample_qc_exclusions.json').read_text())
    excluded = set(policy.get('always_exclude', []))
    if policy.get('enabled', True):
        excluded.update(r['sample'] for r in policy['excluded_samples'])
    qc = {r['sample']: r for r in read_table(root / 'bundle_sample_qc_summary.tsv')}
    rows = []
    refdir = args.output.resolve().parent / 'reference_sizes'
    refdir.mkdir(parents=True, exist_ok=True)
    for sid, name, latin, _, folder, prefix in specs:
        sizes = None
        models = root / 'cross_species_visual_tool/index' / (literals['MODEL_INDEX'].get(sid, sid) + '_models.tsv.gz')
        if not models.is_file():
            raise FileNotFoundError(models)
        for bwpath in sorted((root / 'cross_species_visual_tool' / folder).resolve().glob('*.bw')):
            if not bwpath.name.lower().startswith(prefix.lower()):
                continue
            sample = bwpath.name.removesuffix('.bw').removesuffix('.cpm').removesuffix('_scaled')
            if sample in excluded:
                continue
            peak = bwpath.with_name(sample + '_peaks.narrowPeak')
            if sample in ('ASM3Liver', 'ASM4Liver'):
                peak = root / 'data/atac/liver_atac_bundle/Acomys_cahirinus' / peak.name
            if not peak.exists():
                peak = Path(str(peak) + '.gz')
            if not peak.is_file():
                raise FileNotFoundError(peak)
            with pyBigWig.open(str(bwpath)) as bw:
                current = bw.chroms()
            if sizes is not None and current != sizes:
                raise ValueError('Inconsistent reference sequence dictionaries: ' + sid)
            sizes = current
            value = qc.get(sample, {}).get('tss_enrichment', '')
            if value in ('-', 'NA', 'nan'):
                value = ''
            low = value and float(value) < 4
            sizepath = refdir / (sid + '.chrom.sizes')
            rows.append(dict(species=sid, common_name=name, scientific_name=latin, sample=sample,
                             bam='', peaks=str(peak.relative_to(root)), models=str(models.relative_to(root)),
                             chrom_sizes=str(sizepath.relative_to(root)), exclude_bed='',
                             exclude_chroms='chrM;MT;M;chrMT;chrY', tss_enrichment=value,
                             qc_flag='TSS_below_4_review' if low else ('QC_not_in_summary' if not value else '')))
        if sizes is not None:
            (refdir / (sid + '.chrom.sizes')).write_text(''.join(f'{c}\t{n}\n' for c, n in sorted(sizes.items())))
    write_table(args.output, rows, FIELDS)
    print(f'Exported {len(rows)} samples / {len({r["species"] for r in rows})} species to {args.output}')


def discover(args):
    rows = manifest(args.manifest)
    files = sorted(p.resolve() for p in args.bam_root.rglob('*.bam'))
    unresolved = []
    for row in rows:
        # Prefix must end at a delimiter; do not confuse NMR1 with NMR10.
        pattern = re.compile(re.escape(row['sample']) + r'(?:[._-]|$)', re.I)
        matches = [p for p in files if pattern.match(p.name)]
        if len(matches) == 1:
            row['bam'] = str(matches[0])
        else:
            row['bam'] = ''
            unresolved.append(dict(sample=row['sample'], matches=[str(p) for p in matches]))
    write_table(args.output, rows, list(rows[0]))
    dump(str(args.output) + '.unresolved.json', unresolved)
    print(f'Matched {len(rows)-len(unresolved)}/{len(rows)} BAMs. Review {args.output}; ambiguous matches remain blank.')


def sizes_from(path):
    with Path(path).open() as f:
        sizes = {c: int(n) for c, n, *_ in (line.split() for line in f if line.strip())}
    if not sizes or any(not re.fullmatch(r'[A-Za-z0-9_.-]+', c) or n < 1 for c,n in sizes.items()):
        raise ValueError('Invalid or shell-unsafe reference sequence name/length')
    return sizes


def convert_models(source, output, sizes, excluded):
    """UCSC genePredExt layout, with synthetic unique transcript IDs for ROSE."""
    header = ['bin','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd',
              'exonCount','exonStarts','exonEnds','score','name2','cdsStartStat','cdsEndStat','exonFrames']
    genes = collections.defaultdict(list)
    mapping = []
    with text_open(source) as f, Path(output).open('w') as out:
        out.write('\t'.join(header) + '\n')
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            c, a, b, tx, gene, strand, xs, xe, *_ = line.rstrip('\n').split('\t')
            if c in excluded:
                continue
            a, b = int(a), int(b)
            if c not in sizes or not (0 <= a < b <= sizes[c]) or strand not in ('+', '-'):
                raise ValueError(f'Invalid annotation interval: {c}:{a}-{b} ({source})')
            starts = [int(x) for x in xs.rstrip(',').split(',')]
            ends = [int(x) for x in xe.rstrip(',').split(',')]
            if len(starts) != len(ends) or any(not (a <= x < y <= b) for x, y in zip(starts, ends)):
                raise ValueError('Invalid exon coordinates: ' + tx)
            rid = f'TX{len(mapping)+1:09d}'
            # Coding bounds are not inferred; ROSE only uses transcript/TSS coordinates and names.
            record = [0, rid, c, strand, a, b, a, a, len(starts), xs.rstrip(',')+',',
                      xe.rstrip(',')+',', 0, gene, 'unk', 'unk', ','.join(['-1']*len(starts))+',']
            out.write('\t'.join(map(str, record)) + '\n')
            mapping.append(dict(rose_id=rid, transcript=tx, gene=gene))
            # Same end-coordinate TSS convention as pinned ROSE for minus-strand transcripts.
            genes[c].append((a if strand == '+' else b, gene, tx))
    if not mapping:
        raise ValueError('No usable transcripts')
    write_table(str(output) + '.transcript_ids.tsv', mapping, ['rose_id', 'transcript', 'gene'])
    for c in genes:
        genes[c].sort()
    return genes, len(mapping)


def bed_records(path):
    with text_open(path) as f:
        for lineno, line in enumerate(f, 1):
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            p = line.split()
            try:
                c, a, b = p[0], int(p[1]), int(p[2])
                if a < 0 or b <= a:
                    raise ValueError()
            except (ValueError, IndexError):
                raise ValueError(f'Invalid BED interval at {path}:{lineno}')
            yield c, a, b


def interval_index(records):
    groups = collections.defaultdict(list)
    for c, a, b in records:
        groups[c].append((a, b))
    result = {}
    for c, spans in groups.items():
        spans.sort()
        starts, maxends, end = [], [], -1
        for a, b in spans:
            end = max(end, b)
            starts.append(a)
            maxends.append(end)
        result[c] = (starts, maxends)
    return result


def overlaps(index, c, a, b):
    starts, ends = index.get(c, ([], []))
    i = bisect.bisect_left(starts, b) - 1
    return i >= 0 and ends[i] > a


def prepare_peaks(source, output, sizes, excluded, exclude_bed=None):
    mask = [] if not exclude_bed else list(bed_records(exclude_bed))
    for c, a, b in mask:
        if c not in sizes or b > sizes[c]:
            raise ValueError(f'Exclusion BED/reference mismatch: {c}:{a}-{b}')
    masks = interval_index(mask)
    counts = collections.Counter()
    peaks = set()
    for c, a, b in bed_records(source):
        counts['input'] += 1
        if c not in sizes or b > sizes[c]:
            raise ValueError(f'Peak/reference mismatch: {c}:{a}-{b} in {source}')
        if c in excluded:
            counts['excluded_contig'] += 1
        elif overlaps(masks, c, a, b):
            counts['excluded_mask'] += 1
        elif (c, a, b) in peaks:
            counts['duplicate_interval'] += 1
        else:
            peaks.add((c, a, b))
    if not peaks:
        raise ValueError('No candidate peaks remain')
    with Path(output).open('w') as f:
        for i, (c, a, b) in enumerate(sorted(peaks), 1):
            f.write(f'{c}\t{a}\t{b}\tPEAK{i:09d}\t0\n')
    counts['retained'] = len(peaks)
    return dict(counts)


def check_bam(path, sizes, min_mapq=30, scan=100000):
    import pysam
    pysam.quickcheck(str(path))
    with pysam.AlignmentFile(str(path), 'rb') as bam:
        if bam.header.to_dict().get('HD', {}).get('SO') != 'coordinate':
            raise ValueError('BAM must be coordinate sorted: ' + str(path))
        if not bam.has_index():
            raise ValueError('BAM needs an index (.bai or .csi): ' + str(path))
        actual = dict(zip(bam.references, bam.lengths))
        bad = [c for c, n in sizes.items() if actual.get(c) != n]
        if bad:
            raise ValueError('BAM differs from browser reference sequence dictionary: ' + ', '.join(bad[:8]))
        counts = collections.Counter()
        for i, read in enumerate(bam.fetch(until_eof=True)):
            if i >= scan:
                break
            counts['reads_checked'] += 1
            if read.flag & (4 | 256 | 512 | 1024 | 2048):
                counts['flagged_reads'] += 1
            if read.mapping_quality < min_mapq:
                counts['low_mapq'] += 1
            if not read.is_paired or not read.is_proper_pair:
                counts['not_proper_pair'] += 1
        if not counts['reads_checked']:
            raise ValueError('Empty BAM')
        if counts['flagged_reads'] or counts['low_mapq'] or counts['not_proper_pair']:
            raise ValueError('BAM is not a filtered paired-end ATAC BAM: ' + str(dict(counts)))
    return dict(counts, scan_limit=scan, min_mapq=min_mapq)


def verify_vendor():
    record = json.loads((HERE / 'vendor/source.json').read_text())
    for rel, digest in record['files_sha256'].items():
        if sha(VENDOR / rel) != digest:
            raise ValueError('Vendored ROSE was modified: ' + rel)
    return record


def load_rose():
    verify_vendor()
    sys.path[:0] = [str(VENDOR / 'lib'), str(VENDOR / 'bin')]
    import ROSE_utils as utils
    import ROSE_main as main
    import ROSE_bamToGFF as mapper
    # Upstream adds/removes chr from some scaffold names. Keep exact BAM/peak names instead.
    utils.checkChrStatus = lambda _: 1
    original = utils.makeSearchLocus
    def clipped_search(locus, up, down):
        value = original(locus, up, down)
        return utils.Locus(value.chr(), max(1, value.start()), max(1, value.end()), value.sense(), value.ID())
    utils.makeSearchLocus = clipped_search
    return utils, main, mapper


def annotate_candidates(rows, genes):
    coordinates = {c: [r[0] for r in records] for c, records in genes.items()}
    result = []
    for row in rows:
        c, a, b = row['CHROM'], int(row['START']), int(row['STOP'])
        if b <= a:
            raise ValueError('ROSE produced a non-positive interval')
        candidates = genes.get(c, [])
        positions = coordinates.get(c, [])
        near = bisect.bisect_left(positions, (a+b)/2)
        closest = min(candidates[max(0,near-1):near+1], key=lambda x: abs(x[0]-(a+b)/2), default=None)
        left = bisect.bisect_left(positions, a-50000)
        right = bisect.bisect_right(positions, b+50000)
        result.append(dict(chrom=c, start=a, end=b, region_id=row['REGION_ID'],
                           signal_rpm_bp=float(row[list(row)[6]]),
                           rank=int(row['stitchedPeakRank']), is_candidate_se=int(row['isSuper']),
                           constituent_peaks=int(row['NUM_LOCI']),
                           nearest_gene=closest[1] if closest else '',
                           nearest_transcript=closest[2] if closest else '',
                           nearest_tss_distance_bp=max(a-closest[0], closest[0]-b, 0) if closest else '',
                           genes_tss_within_50kb=';'.join(sorted({g[1] for g in candidates[left:right]}))))
    return result


def run_sample(args):
    rows = manifest(args.manifest)
    if args.index is not None:
        if not 0 <= args.index < len(rows):
            raise ValueError('Sample index outside manifest')
        row = rows[args.index]
    else:
        matches = [r for r in rows if r['sample'] == args.sample]
        if len(matches) != 1:
            raise ValueError('Select one unique sample')
        row = matches[0]
    root = args.project_root.resolve()
    if not row.get('bam'):
        raise ValueError('BAM path is blank for ' + row['sample'])
    paths = {k: resolve(row[k], root).resolve() for k in ('bam','peaks','models','chrom_sizes')}
    if row.get('exclude_bed'):
        paths['exclude_bed'] = resolve(row['exclude_bed'], root).resolve()
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    for program in ('samtools','Rscript'):
        if not shutil.which(program):
            raise RuntimeError('Missing executable: ' + program)
    settings = dict(stitch_bp=args.stitch, tss_exclusion_bp=args.tss,
                    min_mapq=args.min_mapq, min_tss_qc=args.min_tss_qc,
                    read_extension_bp=200, floor=1, normalization='ROSE_RPM',
                    exclude_chroms=row.get('exclude_chroms','').split(';'))
    if args.min_tss_qc is not None:
        if not row.get('tss_enrichment'):
            raise ValueError('TSS QC unavailable; cannot apply requested TSS filter')
        if float(row['tss_enrichment']) < args.min_tss_qc:
            raise ValueError('Sample is below requested TSS QC threshold; remove it from the run manifest')
    settings['exclude_chroms'] = sorted(set(settings['exclude_chroms']) | {'chrY'})
    signature = dict(settings=settings, inputs={k: dict(path=str(p), size=p.stat().st_size,
                   mtime_ns=p.stat().st_mtime_ns, sha256=sha(p) if k != 'bam' else None) for k,p in paths.items()},
                   rose=verify_vendor(), workflow_sha256=sha(__file__))
    job = args.output.resolve() / row['species'] / row['sample']
    if re.search(r'[^A-Za-z0-9_./-]', str(job)):
        raise ValueError('ROSE working/output path must contain only letters, numbers, _, -, /, and .')
    if (job/'SUCCESS.json').exists():
        previous = json.loads((job/'SUCCESS.json').read_text())
        if previous['signature'] == signature and all(sha(job/p)==h for p,h in previous['outputs_sha256'].items()):
            print('Already complete: ' + str(job))
            return
        raise ValueError('Inputs/settings changed; use a new output directory')
    if job.exists():
        raise ValueError('Incomplete output exists; inspect its log and use a new output directory: ' + str(job))
    job.mkdir(parents=True)
    dump(job/'inputs.json', dict(sample=row, signature=signature))
    import pysam
    versions = dict(python=sys.version, pysam=pysam.__version__)
    for tool in ('samtools','Rscript'):
        info = subprocess.run([tool,'--version'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,check=True)
        versions[tool] = info.stdout.decode(errors='replace').splitlines()[0]
    dump(job/'runtime_versions.json', versions)
    start = time.time()
    with (job/'run.log').open('w', buffering=1) as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        try:
            sizes = sizes_from(paths['chrom_sizes'])
            qc = check_bam(paths['bam'], sizes, args.min_mapq)
            dump(job/'bam_validation.json', qc)
            excluded = set(settings['exclude_chroms'])
            genes, n_tx = convert_models(paths['models'], job/'annotation_refseq.ucsc', sizes, excluded)
            peak_counts = prepare_peaks(paths['peaks'], job/'peaks.bed', sizes, excluded, paths.get('exclude_bed'))
            # Safe, short symlink name for legacy upstream shell calls. No input BAM is changed.
            bam = job/'signal.bam'
            bam.symlink_to(paths['bam'])
            for candidate in (Path(str(paths['bam'])+'.bai'), paths['bam'].with_suffix('.bai'),
                              Path(str(paths['bam'])+'.csi'), paths['bam'].with_suffix('.csi')):
                if candidate.is_file():
                    Path(str(bam)+candidate.suffix).symlink_to(candidate)
                    break
            utils, main, mapper = load_rose()
            utils.bedToGFF(str(job/'peaks.bed'), str(job/'peaks.gff'))
            stitched, debug = main.regionStitching(str(job/'peaks.gff'), args.stitch, args.tss,
                                                   str(job/'annotation_refseq.ucsc'), args.tss > 0)
            n_stitched = len(stitched.getLoci())
            if n_stitched < 3:
                raise ValueError('Fewer than three stitched regions remain; cannot rank SEs')
            utils.unParseTable(utils.locusCollectionToGFF(stitched), str(job/'stitched.gff'), '\t')
            utils.unParseTable(debug, str(job/'stitching_debug.tsv'), '\t')
            for stem in ('stitched','peaks'):
                mapped = mapper.mapBamToGFF(str(bam), str(job/(stem+'.gff')), 'both', 200, 1, True, 1)
                values = [float(r[2]) for r in mapped[1:]]
                if len(mapped) < 2 or not all(math.isfinite(v) for v in values):
                    raise ValueError('Missing or nonfinite mapped signal')
                utils.unParseTable(mapped, str(job/(stem+'_signal.bam_MAPPED.gff')), '\t')
            main.mapCollection(stitched, utils.gffToLocusCollection(str(job/'peaks.gff')),
                               [str(bam)], str(job)+'/', str(job/'REGION_MAP.txt'), 'stitched')
            with (job/'REGION_MAP.txt').open() as handle:
                mapped_rows = list(csv.DictReader(handle, delimiter='\t'))
            signals = [float(r['signal.bam']) for r in mapped_rows]
            if not signals or max(signals) <= min(signals):
                raise ValueError('No signal variation; ROSE cutoff is undefined')
            subprocess.run(['Rscript', str(VENDOR/'bin/ROSE_callSuper.R'), str(job)+'/',
                            str(job/'REGION_MAP.txt'), 'ATAC', 'NONE'], check=True, stdout=log, stderr=log)
            with (job/'ATAC_AllStitched.table.txt').open() as f:
                called = list(csv.DictReader((line for line in f if not line.startswith('#')), delimiter='\t'))
            if len(called) != n_stitched:
                raise ValueError('ROSE result count differs from stitched region count')
            result = annotate_candidates(called, genes)
            write_table(job/'all_regions.tsv', result, list(result[0]))
            se = [r for r in result if r['is_candidate_se']]
            write_table(job/'candidate_se.tsv', se, list(result[0]))
            with (job/'candidate_se.bed').open('w') as f:
                for r in sorted(se, key=lambda x:(x['chrom'],x['start'],x['end'])):
                    f.write(f'{r["chrom"]}\t{r["start"]}\t{r["end"]}\t{r["region_id"]}\n')
            report = dict(species=row['species'], sample=row['sample'], common_name=row['common_name'],
                          candidate_se_count=len(se), stitched_regions=n_stitched, annotation_transcripts=n_tx,
                          peaks=peak_counts, qc_flag=row.get('qc_flag',''),
                          blacklist_applied=bool(paths.get('exclude_bed')), elapsed_seconds=time.time()-start)
            dump(job/'summary.json', report)
            outputs = ['all_regions.tsv','candidate_se.tsv','candidate_se.bed','summary.json',
                       'runtime_versions.json','ATAC_Plot_points.png']
            dump(job/'SUCCESS.json', dict(signature=signature, outputs_sha256={p:sha(job/p) for p in outputs}))
        except Exception:
            import traceback
            traceback.print_exc()
            raise
    print(f'{row["sample"]}: {len(se)} candidate SEs; {job}')


def support_segments(sample_beds):
    """Exact half-open segments supported by distinct samples, not chain overlaps."""
    events = collections.defaultdict(lambda: collections.defaultdict(list))
    for sample, spans in sample_beds.items():
        for c,a,b in spans:
            events[c][a].append((sample,1))
            events[c][b].append((sample,-1))
    output = []
    for c, points in sorted(events.items()):
        active = collections.Counter()
        previous = None
        for position, changes in sorted(points.items()):
            names = tuple(sorted(s for s,n in active.items() if n > 0))
            if previous is not None and position > previous and names:
                if output and output[-1][0] == c and output[-1][2] == previous and output[-1][3] == names:
                    output[-1] = (c, output[-1][1], position, names)
                else:
                    output.append((c,previous,position,names))
            for sample, delta in changes:
                active[sample] += delta
            previous = position
    return output


def merged_support(segments, minimum):
    merged = []
    for c,a,b,names in segments:
        if len(names) < minimum:
            continue
        if merged and merged[-1][0] == c and merged[-1][2] == a:
            merged[-1] = (c,merged[-1][1],b)
        else:
            merged.append((c,a,b))
    return merged


def summarize(args):
    rows = manifest(args.manifest)
    samples, grouped, missing, signatures = [], collections.defaultdict(dict), [], []
    for row in rows:
        job = args.output / row['species'] / row['sample']
        if not (job/'SUCCESS.json').is_file():
            missing.append(row['sample'])
            continue
        success = json.loads((job/'SUCCESS.json').read_text())
        for p,h in success['outputs_sha256'].items():
            if sha(job/p) != h:
                raise ValueError('Output changed after completion: ' + str(job/p))
        signatures.append(success['signature'])
        samples.append(json.loads((job/'summary.json').read_text()))
        grouped[row['species']][row['sample']] = list(bed_records(job/'candidate_se.bed'))
    if missing:
        raise ValueError('Some samples are incomplete. Subset the manifest explicitly for a partial summary: '+', '.join(missing))
    common = {(s['settings']['stitch_bp'],s['settings']['tss_exclusion_bp'],s['settings']['min_mapq'],s['settings']['min_tss_qc'],
               s['rose']['commit'],s['workflow_sha256']) for s in signatures}
    if len(common) != 1:
        raise ValueError('Cannot combine runs with different analysis settings/versions')
    destination = args.output/'species_summary'
    destination.mkdir(parents=True, exist_ok=True)
    sample_rows = [{k:v for k,v in r.items() if k != 'peaks'} for r in samples]
    write_table(destination/'sample_counts.tsv', sample_rows, list(sample_rows[0]))
    species_rows = []
    for sid, beds in sorted(grouped.items()):
        segments = support_segments(beds)
        write_table(destination/(sid+'_support_segments.tsv'),
                    [dict(chrom=c,start=a,end=b,n_samples=len(names),samples=';'.join(names)) for c,a,b,names in segments],
                    ['chrom','start','end','n_samples','samples'])
        record = dict(species=sid, n_samples=len(beds))
        for label, minimum in [('union',1),('shared_ge2',2)]:
            merged = merged_support(segments, minimum)
            with (destination/(sid+'_candidate_se_'+label+'.bed')).open('w') as f:
                for i,(c,a,b) in enumerate(merged,1):
                    f.write(f'{c}\t{a}\t{b}\t{sid}_{label}_{i:06d}\n')
            record[label+'_regions'] = len(merged)
            record[label+'_bp'] = sum(b-a for _,a,b in merged)
        record['replicate_support_status'] = 'available' if len(beds)>1 else 'single_sample_only'
        species_rows.append(record)
    write_table(destination/'species_counts.tsv', species_rows, list(species_rows[0]))
    print('Wrote per-sample and species summaries to '+str(destination))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    ex = sub.add_parser('export-manifest', help='Snapshot current browser species, peaks, annotations and sequence sizes')
    ex.add_argument('--project-root',type=Path,required=True)
    ex.add_argument('--output',type=Path,required=True)
    ex.set_defaults(func=export_manifest)
    discover_parser = sub.add_parser('discover-bams', help='Suggest unambiguous sample-name-prefix BAM matches')
    discover_parser.add_argument('--manifest',type=Path,required=True)
    discover_parser.add_argument('--bam-root',type=Path,required=True)
    discover_parser.add_argument('--output',type=Path,required=True)
    discover_parser.set_defaults(func=discover)
    run = sub.add_parser('run-sample')
    run.add_argument('--manifest',type=Path,required=True)
    run.add_argument('--project-root',type=Path,required=True)
    run.add_argument('--output',type=Path,required=True)
    select = run.add_mutually_exclusive_group(required=True)
    select.add_argument('--sample')
    select.add_argument('--index',type=int,help='Zero-based manifest data-row index')
    run.add_argument('--stitch',type=int,default=12500)
    run.add_argument('--tss',type=int,default=2500)
    run.add_argument('--min-mapq',type=int,default=30)
    run.add_argument('--min-tss-qc',type=float)
    run.set_defaults(func=run_sample)
    summ = sub.add_parser('summarize')
    summ.add_argument('--manifest',type=Path,required=True)
    summ.add_argument('--output',type=Path,required=True)
    summ.set_defaults(func=summarize)
    args = parser.parse_args()
    if hasattr(args,'stitch') and (args.stitch < 0 or args.tss < 0 or args.min_mapq < 0):
        parser.error('Distances and mapping-quality thresholds must be nonnegative')
    args.func(args)


if __name__ == '__main__':
    main()
