import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pyBigWig

import activity
import peak_activity
import regulatory_strength
from region_availability import annotate_region, body_region, coverage


class RegionAvailabilityTest(unittest.TestCase):
    def test_boundary_geometry_and_strand(self):
        for start, end, expected in [(-501, 499, False), (-500, 500, True),
                                     (-499, 501, True), (500, 1500, True),
                                     (501, 1501, False), (1000, 2000, False)]:
            self.assertEqual(coverage('1', start, end, {'chr1': 1000})['region_eligible'], expected)
        self.assertFalse(coverage('missing', 0, 100, {'chr1': 1000})['region_eligible'])
        self.assertFalse(coverage('chr1', 0, 0, {'chr1': 1000})['region_eligible'])
        self.assertEqual(body_region('G|chr1:100-200:+'), ('chr1', -1900, 200))
        self.assertEqual(body_region('G|chr1:100-200:-'), ('chr1', 100, 2200))

    def test_all_window_methods_use_genomic_span_not_signal_or_peak_coverage(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with pyBigWig.open(str(root/'sample.bw'), 'w') as bw:
                bw.addHeader([('chr1', 1000)])
                bw.addEntries(['chr1'], [400], ends=[410], values=[2.])
            (root/'sample_peaks.narrowPeak').write_text('chr1\t400\t410\n')
            specs = [('test', '', '', '', str(root), '')]
            row = dict(species='test', locus='G|chr1:400-601:+',
                       samples={'sample': None}, unit='CPM_signal_bp')
            methods = [lambda rows, flank: activity.tss_records(rows, 1, flank),
                       regulatory_strength.score_records, peak_activity.score_records]
            with patch('server.SPECS', specs), patch('calculate_gene_activity.SPECS', specs):
                for method in methods:
                    for strand in ('+', '-'):
                        locus = dict(row, locus='G|chr1:400-601:'+strand)
                        half = method([locus], 1000)[0]
                        below = method([locus], 1001)[0]
                        self.assertEqual(half['available_fraction'], .5)
                        self.assertTrue(half['region_eligible'])
                        self.assertFalse(below['region_eligible'])
                        self.assertEqual(below['available_bp'], 1000)
                        self.assertIsNotNone(below['median'])
                signal = peak_activity.score_records([row], 1000)[0]
                self.assertEqual(signal['median'], 20.)
                self.assertEqual(signal['window_details']['sample']['peak_bp'], 10)
                self.assertTrue(signal['region_eligible'])  # 10 peak bp is not the denominator
                zero = peak_activity.score_records([dict(row, locus='G|chr1:450-601:+')], 10)[0]
                self.assertEqual(zero['median'], 0.)
                self.assertTrue(zero['region_eligible'])

    def test_measured_sample_bounds_and_missing_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for sample, size in [('a', 1000), ('b', 400)]:
                with pyBigWig.open(str(root/(sample+'.bw')), 'w') as bw:
                    bw.addHeader([('chr1', size)])
                    bw.addEntries(['chr1'], [0], ends=[size], values=[0.])
            row = dict(samples={'a': 0., 'b': 0.}, median=0.)
            annotate_region(row, 'chr1', 0, 1000, root)
            self.assertFalse(row['region_eligible'])
            self.assertEqual(row['available_fraction'], .4)
            row['samples']['b'] = None
            annotate_region(row, 'chr1', 0, 1000, root)
            self.assertTrue(row['region_eligible'])  # missing sample is not a zero measurement
            row['samples']['missing'] = 5.
            annotate_region(row, 'chr1', 0, 1000, root)
            self.assertFalse(row['region_eligible'])
            self.assertIsNone(row['available_fraction'])

    def test_muskrat_sirt6_geometry(self):
        sizes = {'PVIU01017863.1': 40717}
        self.assertEqual(coverage('PVIU01017863.1', 14537-50000, 14537+50000, sizes)['available_fraction'], .40717)
        self.assertTrue(coverage('PVIU01017863.1', 14537-20000, 14537+20000, sizes)['region_eligible'])
        self.assertEqual(coverage(*body_region('SIRT6|PVIU01017863.1:14537-18958:+'), sizes)['available_fraction'], 1.)


if __name__ == '__main__':
    unittest.main()
