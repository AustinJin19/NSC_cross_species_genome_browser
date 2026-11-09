import csv
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import atac_super_enhancers as se
from import_atac_super_enhancers import supported_spans


class AtacSuperEnhancers(unittest.TestCase):
    def test_coordinates_overlaps_and_unavailable_status(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(se, 'DATA', Path(tmp)):
            root = Path(tmp)
            info = dict(file='Mouse.bed', regions=3, sharedBp=120,
                        nSamples=2, samples=['a', 'b'], supportStatus='available')
            (root / 'metadata.json').write_text(json.dumps({'species': {'Mouse': info}}))
            (root / 'Mouse.bed').write_text('chr1\t100\t200\tfirst\nchr1\t110\t120\tinside\nchr1\t200\t210\tnext\n')
            self.assertEqual(se.query('Mouse', 'chr1', 100, 100)['regions'], [])
            first = se.query('Mouse', 'chr1', 101, 101)['regions'][0]
            self.assertEqual((first['start'], first['end']), (101, 200))
            self.assertEqual([r['name'] for r in se.query('Mouse', 'chr1', 200, 200)['regions']], ['first'])
            self.assertEqual([r['name'] for r in se.query('Mouse', 'chr1', 201, 201)['regions']], ['next'])
            self.assertEqual(se.query('Mouse', 'chr2', 1, 500)['regions'], [])
            self.assertEqual(se.query('Mouse', 'chr1', -10, -1)['regions'], [])
            self.assertEqual(se.query('Mouse', 'chr1', 10, 1)['regions'], [])
            self.assertEqual(se.query('NMR', 'chr1', 1, 500)['status'], 'not_provided')
            self.assertEqual(se.query('Mouse', 'chr1', 1, 500, 'union')['status'], 'not_provided')
            (root / 'Mouse.bed').unlink()
            self.assertEqual(se.query('Mouse', 'chr1', 1, 500)['status'], 'file_unavailable')

    def test_single_sample_union_is_available_without_shared_fallback(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(se, 'DATA', Path(tmp)):
            root = Path(tmp)
            info = dict(file='shared.bed', regions=0, sharedBp=0, nSamples=1,
                        samples=['one'], supportStatus='single_sample_only',
                        union=dict(file='union.bed', regions=1, coveredBp=100, supportStatus='available'))
            (root / 'metadata.json').write_text(json.dumps({'species': {'Bat': info}}))
            (root / 'shared.bed').write_text('')
            (root / 'union.bed').write_text('chr1\t100\t200\tunion1\n')
            shared = se.query('Bat', 'chr1', 1, 300)
            union = se.query('Bat', 'chr1', 1, 300, 'union')
            self.assertEqual(shared['status'], 'single_sample_only')
            self.assertEqual(shared['regions'], [])
            self.assertEqual(union['status'], 'available')
            self.assertEqual(union['minimumSupport'], 1)
            self.assertEqual(union['coveredBp'], 100)
            self.assertEqual(union['regions'], [dict(start=101, end=200, name='union1')])
            self.assertEqual(se.query('Bat', 'chr1', 100, 100, 'union')['regions'], [])
            self.assertEqual(se.query('Bat', 'chr1', 201, 201, 'union')['regions'], [])
            (root / 'union.bed').unlink()
            self.assertEqual(se.query('Bat', 'chr1', 1, 300, 'union')['status'], 'file_unavailable')
            with self.assertRaises(ValueError):
                se.query('Bat', 'chr1', 1, 300, 'typo')

    def test_shared_bases_not_transitive_union(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'support.tsv'
            path.write_text('chrom\tstart\tend\tn_samples\tsamples\n'
                            'chr1\t0\t10\t1\ta\nchr1\t10\t20\t2\ta;b\n'
                            'chr1\t20\t30\t2\tb;c\nchr1\t30\t40\t1\tc\n'
                            'chr1\t40\t50\t2\ta;c\n')
            self.assertEqual(supported_spans(path, {'a', 'b', 'c'}), [('chr1', 10, 30), ('chr1', 40, 50)])
            self.assertEqual(supported_spans(path, {'a', 'b', 'c'}, minimum=1), [('chr1', 0, 50)])
            with self.assertRaises(ValueError):
                supported_spans(path, {'a', 'b'})

    @unittest.skipUnless((se.DATA / 'metadata.json').exists(), 'ATAC SE data not installed')
    def test_installed_data_counts_hashes_and_browser_integration(self):
        import server
        server.load_catalog()
        meta = json.loads((se.DATA / 'metadata.json').read_text())
        self.assertEqual(len(meta['species']), 31)
        self.assertEqual(sum(s['regions'] for s in meta['species'].values()), 48377)
        for sid, info in meta['species'].items():
            for selected in (info, info['union']):
                path = se.DATA / selected['file']
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), selected['sha256'])
                with path.open() as handle:
                    rows = list(csv.reader(handle, delimiter='\t'))
                self.assertEqual(len(rows), selected['regions'])
                self.assertEqual(sum(int(r[2])-int(r[1]) for r in rows), selected.get('coveredBp', selected.get('sharedBp')))
            result = se.query(sid, 'chr1', 1, 1)
            self.assertTrue(result['available'])
            union = se.query(sid, 'chr1', 1, 1, 'union')
            self.assertTrue(union['available'])
            self.assertEqual(union['status'], 'available')
            self.assertGreater(union['total'], 0)
            if info['nSamples'] == 1:
                self.assertEqual(result['status'], 'single_sample_only')
                self.assertEqual(result['regions'], [])
        for sid in ['Human', 'Mouse', 'NMR', 'Macaque']:
            p = server.get_panel(sid, 'SIRT6', 500000, 0, None, 10)
            self.assertNotIn('error', p)
            track = p['atacSuperEnhancers']
            self.assertTrue(track['regions'], sid)
            self.assertTrue(all(r['start'] <= p['to'] and r['end'] >= p['from'] for r in track['regions']))
            union = p['atacSuperEnhancerUnion']
            self.assertTrue(union['regions'], sid)
            self.assertEqual(union['minimumSupport'], 1)
            self.assertTrue(all(any(u['start'] <= r['start'] and u['end'] >= r['end'] for u in union['regions'])
                                for r in track['regions']))
        samples = list(server.SAMPLES['NMR'])[:2]
        tracks = [server.get_panel('NMR', 'SIRT6', 500000, 0, sample, 10)['atacSuperEnhancers'] for sample in samples]
        self.assertEqual(tracks[0], tracks[1])
        self.assertEqual(se.query('Balaena_mysticetus', 'chr1', 1, 100)['status'], 'not_provided')


if __name__ == '__main__':
    unittest.main()
