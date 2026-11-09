import json
import sqlite3
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pyBigWig

import activity
from signal_coverage import annotate_signal, positive_bases


def track(path, intervals, size=1000):
    with pyBigWig.open(str(path), 'w') as bw:
        bw.addHeader([('chr1', size)])
        if intervals:
            bw.addEntries(['chr1'] * len(intervals), [a for a, b, v in intervals],
                          ends=[b for a, b, v in intervals], values=[float(v) for a, b, v in intervals])


class PositiveSignalCoverage(unittest.TestCase):
    def test_exact_overlap_zeros_negatives_missing_and_full_denominator(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            track(root/'a.bw', [(0, 500, 2), (500, 700, 0), (700, 800, -1), (800, 900, float('nan'))])
            for start, end, bp, eligible in [(0, 1000, 500, True), (-1, 1000, 500, False),
                                            (100, 600, 400, True), (500, 1000, 0, False),
                                            (1000, 2000, 0, False)]:
                row = dict(samples={'a': 100.}, median=100.)
                annotate_signal(row, '1', start, end, root)
                self.assertEqual(row['positive_signal_bp'], bp)
                self.assertEqual(row['positive_signal_fraction'], bp/(end-start))
                self.assertEqual(row['signal_eligible'], eligible)
                self.assertEqual(row['median'], 100.)  # retain the original score
            annotate_signal(row, 'missing', 0, 1000, root)
            self.assertIsNone(row['positive_signal_fraction'])
            self.assertFalse(row['signal_eligible'])

    def test_species_median_uses_measured_samples_and_unknown_is_not_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            track(root/'a.bw', [(0, 400, 1)])
            track(root/'b.bw', [(0, 600, 1)])
            row = dict(samples={'a': 400., 'b': 600., 'unmeasured': None})
            annotate_signal(row, 'chr1', 0, 1000, root)
            self.assertEqual(row['positive_signal_fraction'], .5)
            self.assertTrue(row['signal_eligible'])
            self.assertFalse(row['window_details']['a']['signal_eligible'])
            row['samples']['unmeasured'] = 10.
            annotate_signal(row, 'chr1', 0, 1000, root)
            self.assertIsNone(row['positive_signal_fraction'])
            self.assertFalse(row['signal_eligible'])

    def test_live_query_uses_genomic_span_without_positive_signal_filter(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            track(root/'sample.bw', [(1000, 4000, 2)], size=10000)
            db = root/'activity.sqlite'
            (root/'species_mls.json').write_text('{}')
            with sqlite3.connect(db) as con:
                con.execute('CREATE TABLE activity (species TEXT,symbol TEXT,symbol_key TEXT,locus TEXT,median REAL,unit TEXT,samples TEXT)')
                con.execute('INSERT INTO activity VALUES (?,?,?,?,?,?,?)',
                            ('test', 'G', 'G', 'G|chr1:3000-4000:+', 6000., 'CPM_signal_bp', json.dumps({'sample': 6000.})))
            specs = [('test', 'Test', 'Test species', '#333', str(root), 'sample')]
            with patch.object(activity, 'DB', db), patch('server.SPECS', specs), \
                 patch('gene_names.candidates', return_value=['G']), \
                 patch('signal_coverage.annotate_signal', side_effect=AssertionError('Live queries must not apply the historical signal filter')), \
                 patch('gene_names.match_records', side_effect=lambda gene, records: (records, {})):
                for strategy in ['body', 'tss2kb', 'tss5kb', 'tss100kb', 'tss200kb', 'tss500kb']:
                    result = activity.query('G', strategy)
                    row = result['records'][0]
                    self.assertNotIn('positive_signal_fraction', row)
                    self.assertNotIn('signal_eligible', row)
                    self.assertEqual(row['region_eligible'], row['available_fraction'] >= .5)
                    self.assertEqual(row['region_eligible'], strategy in ('body', 'tss2kb', 'tss5kb'))
                    self.assertEqual(row['median'], 6000.)
                    self.assertNotIn('minimum_positive_signal_fraction', result)
                    self.assertEqual(result['minimum_available_fraction'], .5)


if __name__ == '__main__':
    unittest.main()
