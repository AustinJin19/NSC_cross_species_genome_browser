#!/usr/bin/env python3
"""Local cross-species ATAC browser. Run: python3 server.py"""
import argparse
import sys
import secrets
import threading
import bisect
import csv
import io
import math
import re
import gzip
import json
from collections import defaultdict, OrderedDict
from functools import lru_cache
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse, parse_qs

import numpy as np
import pyBigWig
import pysam
import chipseq
import ccre
import repeats
import microc
import abc_tracks
import human_hic

ROOT = Path(__file__).resolve().parent.parent
STATIC = Path(__file__).resolve().parent / 'www'
# Optional project-local package install; ordinary pip installs also work.
sys.path.insert(0, str(Path(__file__).resolve().parent / '.deps'))
SPECS = [
    ('Human', 'Human (hg38)', 'Homo sapiens', '#5275b5', '../data/atac/selected3', ''),
    ('NMR', 'Naked mole-rat', 'Heterocephalus glaber', '#477763', '../data/atac/mole_rat_peaks', 'NMR'),
    ('BMR', 'Blind mole-rat', 'Nannospalax galili', '#927598', '../data/atac/mole_rat_peaks', 'BMR'),
    ('DMR', 'Damaraland mole-rat', 'Fukomys damarensis', '#bd9158', '../data/atac/mole_rat_peaks', 'DMR'),
    ('Mouse', 'Mouse (mm10)', 'Mus musculus', '#6389a1', '../data/atac/mouse_mm10', 'wtmice'),
    ('Rat', 'Rat', 'Rattus norvegicus', '#90945b', '../data/atac/rat', 'Rat'),
    ('Macaque', 'Rhesus macaque', 'Macaca mulatta', '#b66579', '../data/atac/macaque', 'SRR'),
    ('Rabbit', 'Eastern cottontail', 'Sylvilagus floridanus', '#997342', '../data/atac/Rabbit', 'CTR'),
    ('ASM', 'African spiny mouse', 'Acomys cahirinus', '#369a98', '../data/atac/Africa_spiny_mouse', 'ASM'),
]
MODEL_INDEX = {'Mouse': 'Mouse_mm10'}
EXCLUDED_SAMPLES = {'NMR1Liver'}
CATALOG, SAMPLES, GENES = [], {}, {}
PDF_DOWNLOADS = OrderedDict()
PDF_LOCK = threading.Lock()


def load_catalog():
    """Read existing portable BED-like transcript indexes; no R/RDS dependency."""
    CATALOG.clear(); SAMPLES.clear(); GENES.clear()
    for sid, name, latin, color, folder, prefix in SPECS:
        paths = sorted(p for p in (ROOT / folder).glob('*.bw')
                       if p.name.lower().startswith(prefix.lower())
                       and p.name.removesuffix('.bw').removesuffix('.cpm').removesuffix('_scaled') not in EXCLUDED_SAMPLES)
        SAMPLES[sid] = {p.name.removesuffix('.bw').removesuffix('.cpm').removesuffix('_scaled'): p for p in paths}
        if not paths:
            raise RuntimeError(f'No bigWig tracks for {sid}')
        CATALOG.append(dict(id=sid, name=name, latin=latin, color=color,
                            samples=list(SAMPLES[sid]), default=next(iter(SAMPLES[sid]))))
        grouped = {}
        with gzip.open(ROOT / 'index' / f'{MODEL_INDEX.get(sid, sid)}_models.tsv.gz', 'rt') as handle:
            for line in handle:
                chrom, start, end, tx, symbol, strand, *_ = line.rstrip().split('\t')
                key = (symbol.upper(), chrom, strand)
                start, end = int(start) + 1, int(end)
                if key in grouped:
                    g = grouped[key]
                    g['start'], g['end'] = min(g['start'], start), max(g['end'], end)
                else:
                    grouped[key] = dict(chrom=chrom, start=start, end=end, strand=strand)
        genes = defaultdict(list)
        for (symbol, chrom, strand), g in grouped.items():
            g['tss'] = g['end'] if strand == '-' else g['start']
            genes[symbol].append(g)
        if sid == 'Human':
            anchors = json.loads((ROOT / 'index/Human_anchors.json').read_text())
            by_locus = {(a['symbol'], a['chrom'], a['strand']): a for a in anchors}
            for symbol, loci in genes.items():
                for g in loci:
                    anchor = by_locus.get((symbol, g['chrom'], g['strand']))
                    if anchor:
                        g.update(anchor)
        GENES[sid] = genes


def resolve_chrom(chrom, names):
    if chrom in names:
        return chrom
    candidates = [n for n in names if n.removeprefix('chr') == chrom.removeprefix('chr')]
    return candidates[0] if len(candidates) == 1 else None


@lru_cache(maxsize=32)
def peak_index(path):
    stem = Path(path).name.removesuffix('.bw').removesuffix('.cpm').removesuffix('_scaled')
    path = Path(path).with_name(stem + '_peaks.narrowPeak')
    if not path.exists():
        path = Path(str(path) + '.gz')
    if not path.exists():
        return None
    groups = defaultdict(list)
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as handle:
        for line in handle:
            fields = line.split('\t')
            groups[fields[0]].append((int(fields[1]) + 1, int(fields[2])))
    result = {}
    for chrom, rows in groups.items():
        rows.sort()
        starts = [r[0] for r in rows]
        prefix_end = np.maximum.accumulate([r[1] for r in rows]).tolist()
        result[chrom] = (rows, starts, prefix_end)
    return result


def get_panel(sid, gene, window, pan, sample, bins=700):
    hits = GENES[sid].get(gene, [])
    if not hits:
        return dict(id=sid, error=f'{gene} is absent from this annotation')
    hit = max(hits, key=lambda g: g['end'] - g['start'])
    if sample not in SAMPLES[sid]:
        sample = next(iter(SAMPLES[sid]))
    path = SAMPLES[sid][sample]
    strand = hit['strand']
    sgn = -1 if strand == '-' else 1
    centre = hit['tss'] + sgn * pan
    lo, hi = math.floor(centre - window), math.floor(centre + window)
    width = hi - lo + 1
    # Read exact per-bin sums instead of allocating an array for every base.
    # Uncovered bases within the chromosome are zero; outside it is missing.
    count = min(bins, width)
    edges = [i * width // count for i in range(count + 1)]
    signal = []
    with pyBigWig.open(str(path)) as bw:
        chrom = resolve_chrom(hit['chrom'], bw.chroms())
        if chrom is None:
            return dict(id=sid, error='Annotation sequence is absent from this bigWig')
        chrom_size = bw.chroms(chrom)
        for a, b in zip(edges[:-1], edges[1:]):
            left, right = (hi-b+1, hi-a) if strand == '-' else (lo+a, lo+b-1)
            start, end = max(1, left), min(chrom_size, right)
            value = None
            if start <= end:
                total = bw.stats(chrom, start-1, end, type='sum', exact=True)[0]
                value = float(total or 0) / (end-start+1)
            signal.append(dict(rel=float((a+b-1)/2 + (centre-hi if strand == '-' else lo-centre)), value=value))
    models = []
    with pysam.TabixFile(str(ROOT / 'index' / f'{MODEL_INDEX.get(sid, sid)}_models.tsv.gz')) as tb:
        if hi > 0 and lo <= chrom_size and hit['chrom'] in tb.contigs:
            for line in tb.fetch(hit['chrom'], max(0, lo-1), min(chrom_size, hi)):
                f = line.split('\t')
                models.append(dict(start=int(f[1])+1, end=int(f[2]), tx=f[3], gname=f[4], strand=f[5], ex_s=f[6], ex_e=f[7]))
    selected = {}
    for model in sorted(models, key=lambda m: (m['tx'] == hit.get('transcript'), m['end']-m['start']), reverse=True):
        selected.setdefault(model['gname'], model)
    peak_warning = 'Peak coordinates do not match this assembly' if sid == 'ASM' else None
    peaks, peak_data = [], None if peak_warning else peak_index(path)
    if peak_data is not None:
        pc = resolve_chrom(hit['chrom'], peak_data)
        if pc:
            rows, starts, prefix_end = peak_data[pc]
            first, last = bisect.bisect_left(prefix_end, lo), bisect.bisect_right(starts, hi)
            peaks = [dict(start=a, end=b) for a,b in rows[first:last] if b >= lo]
    return dict(id=sid, chrom=hit['chrom'], **{'from':lo, 'to':hi}, strand=strand,
                tss=hit['tss'], tssTranscript=hit.get('transcript'), tssMethod=hit.get('tssMethod'), centre=centre, sample=sample, signalUnit='scaled' if sid == 'Human' else 'CPM', ambiguous=len(hits)>1,
                signal=signal, hic=human_hic.query(hit['chrom'],lo,hi) if sid=='Human' else None, models=sorted(selected.values(), key=lambda m: (m['gname'].upper() != gene, -(m['end']-m['start'])))[:5], modelCount=len(selected), peaks=peaks,
                abc=abc_tracks.query(hit['chrom'],lo,hi) if sid=='NMR' else None, microc=microc.query(sample,hit['chrom'],lo,hi) if sid=='NMR' else None, peakAvailable=peak_data is not None, peakWarning=peak_warning, ccre=ccre.query(hit['chrom'],lo,hi) if sid=='Human' else None, chipTracks=chipseq.query(hit['chrom'],lo,hi,sid), teRepeats=repeats.query(hit['chrom'],lo,min(hi,chrom_size)) if sid=='NMR' else None)


def browse(req):
    gene = str(req.get('gene', 'HAS2')).strip().upper()
    if not gene or len(gene) > 100:
        raise ValueError('Enter a gene symbol')
    window = float(req.get('window', 5000))
    if not math.isfinite(window) or window < 0.5 or (window * 2) % 1:
        raise ValueError('Window span must be a positive whole number of base pairs')
    pan = int(req.get('pan', 0))
    if abs(pan) > 2_000_000_000:
        raise ValueError('Pan exceeds supported coordinate range')
    bins = int(req.get('bins', 700))
    if not 10 <= bins <= 5000:
        raise ValueError('Bin count must be between 10 and 5000')
    panels = []
    panel_specs = req.get('panels')
    if panel_specs is None:
        panel_specs = [dict(key=sid, species=sid, sample=req.get('samples', {}).get(sid))
                       for sid in dict.fromkeys(req.get('species', list(GENES)))]
    if not isinstance(panel_specs, list) or len(panel_specs) > 40:
        raise ValueError('Use at most 40 panels')
    seen = set()
    for spec in panel_specs:
        sid = spec['species']
        key = str(spec['key'])
        if key in seen or len(key) > 100:
            raise ValueError('Panel identifiers must be unique and at most 100 characters')
        seen.add(key)
        if sid not in GENES:
            raise ValueError('Unknown species')
        sample = spec.get('sample')
        if sample is not None and sample not in SAMPLES[sid]:
            raise ValueError('Unknown sample for this species')
        panel = get_panel(sid, gene, window, pan, sample, bins)
        panel['panelKey'] = key
        panels.append(panel)
    return dict(gene=gene, window=window, pan=pan, bins=bins, panels=panels, annotation=annotation_counts(gene))


def annotation_counts(gene):
    species = [sid for sid in GENES if gene.upper() in GENES[sid]]
    return dict(gene=gene.upper(), species=len(species), samples=sum(len(SAMPLES[sid]) for sid in species),
                totalSpecies=len(GENES), totalSamples=sum(len(s) for s in SAMPLES.values()))


def export_csv(data):
    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow(['panel_id','species','sample','gene','chromosome','strand','relative_to_tss_bp','genomic_bin_center_1based','mean_signal','signal_unit'])
    for p in data['panels']:
        if 'error' in p:
            continue
        for v in p['signal']:
            writer.writerow([p['panelKey'], p['id'], p['sample'], data['gene'], p['chrom'], p['strand'],
                             v['rel']+data['pan'], p['centre']+(-1 if p['strand']=='-' else 1)*v['rel'], v['value'], p['signalUnit']])
    return output.getvalue()


class Handler(BaseHTTPRequestHandler):
    def send(self, status, data, content_type='application/json', filename=None):
        if content_type == 'application/json':
            data = json.dumps(data, allow_nan=False).encode()
        self.send_response(status)
        self.send_header('Content-Type', content_type)
        self.send_header('Content-Length', str(len(data)))
        if filename:
            self.send_header('Content-Disposition', f'attachment; filename="{filename}"')
        self.send_header('Cache-Control', 'no-store')
        self.send_header('X-Content-Type-Options', 'nosniff')
        self.end_headers()
        try:
            self.wfile.write(data)
        except (BrokenPipeError, ConnectionResetError):
            pass

    def do_GET(self):
        url = urlparse(self.path)
        if url.path in ('/api/omics','/api/omics-export'):
            try:
                import omics
                gene = parse_qs(url.query).get('gene', ['HAS2'])[0]
                if url.path == '/api/omics-export':
                    name = re.sub(r'[^A-Za-z0-9_-]', '_', gene) + '_liver_omics.csv'
                    return self.send(200, omics.export_csv(gene).encode(), 'text/csv; charset=utf-8', name)
                return self.send(200, omics.query(gene))
            except (ValueError, OSError) as exc:
                return self.send(400, {'error':str(exc)})
        if url.path.startswith('/api/download-pdf/'):
            with PDF_LOCK:
                item = PDF_DOWNLOADS.get(url.path.rsplit('/', 1)[-1])
            if item is None:
                return self.send(404, {'error':'PDF expired; export it again'})
            return self.send(200, item[0], 'application/pdf', item[1])
        if url.path == '/api/export':
            try:
                raw = parse_qs(url.query).get('view', ['{}'])[0]
                if len(raw) > 20000:
                    raise ValueError('Request too large')
                data = browse(json.loads(raw))
                filename = re.sub(r'[^A-Za-z0-9_-]', '_', data['gene']) + '_signal.csv'
                return self.send(200, export_csv(data).encode(), 'text/csv; charset=utf-8', filename)
            except (ValueError, TypeError, KeyError, AttributeError) as exc:
                return self.send(400, {'error':str(exc)})
        if url.path == '/api/catalog':
            return self.send(200, CATALOG)
        if url.path == '/api/genes':
            q = parse_qs(url.query).get('q', [''])[0].upper()
            keys = set().union(*(g.keys() for g in GENES.values()))
            return self.send(200, dict(choices=[annotation_counts(k) for k in sorted(k for k in keys if k.startswith(q))[:30]], match=annotation_counts(q)))
        files = {'/statistics.js':('statistics.js','text/javascript'), '/omics.js':('omics.js','text/javascript'), '/ur-logo.svg':('ur-logo.svg','image/svg+xml'), '/':('index.html','text/html; charset=utf-8'), '/style.css':('style.css','text/css'), '/browser.js':('browser.js','text/javascript')}
        if url.path not in files:
            return self.send(404, {'error':'Not found'})
        name, mime = files[url.path]
        return self.send(200, (STATIC / name).read_bytes(), mime)

    def do_POST(self):
        if self.path == '/api/export-pdf':
            try:
                from pdf_export import export_pdf
                size = int(self.headers.get('Content-Length', 0))
                if not 0 < size <= 20_000_000:
                    raise ValueError('Invalid PDF request size')
                payload = json.loads(self.rfile.read(size))
                pdf = export_pdf(payload)
                filename = re.sub(r'[^A-Za-z0-9_-]', '_', str(payload.get('gene','region'))) + '_NSC.pdf'
                token = secrets.token_urlsafe(24)
                with PDF_LOCK:
                    PDF_DOWNLOADS[token] = (pdf, filename)
                    while len(PDF_DOWNLOADS) > 8:
                        PDF_DOWNLOADS.popitem(last=False)
                return self.send(200, {'url':'/api/download-pdf/' + token})
            except ImportError:
                return self.send(500, {'error':'Install the PDF dependencies from requirements.txt'})
            except Exception as exc:
                return self.send(400, {'error':'PDF export failed: ' + str(exc)})
        if self.path != '/api/browse':
            return self.send(404, {'error':'Not found'})
        try:
            size = int(self.headers.get('Content-Length', 0))
            if not 0 < size <= 20000:
                raise ValueError('Invalid request size')
            req = json.loads(self.rfile.read(size))
            if not isinstance(req, dict):
                raise ValueError('Expected an object')
            self.send(200, browse(req))
        except (ValueError, TypeError, KeyError, AttributeError) as exc:
            self.send(400, {'error':str(exc)})
        except Exception as exc:
            print(f'Data query failed: {exc}', flush=True)
            self.send(500, {'error':'Unable to read this interval. Check the server log.'})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--port', type=int, default=8765)
    args = parser.parse_args()
    load_catalog()
    print(f'NSC: {len(CATALOG)} species ready at http://127.0.0.1:{args.port}', flush=True)
    ThreadingHTTPServer(('127.0.0.1', args.port), Handler).serve_forever()
