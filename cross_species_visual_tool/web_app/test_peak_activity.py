import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pyBigWig

import peak_activity


class PeakActivityTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.bw = self.root / 'sample.bw'
        with pyBigWig.open(str(self.bw), 'w') as bw:
            bw.addHeader([('chr1', 1000), ('chr2', 100)])
            bw.addEntries(['chr1'] * 3, [0, 100, 300], ends=[100, 200, 1000], values=[2., 3., 1.])
        self.bed = self.root / 'sample_peaks.narrowPeak'
        self.bed.write_text('chr1\t20\t80\nchr1\t40\t60\nchr1\t70\t120\nchr1\t120\t150\nchr1\t200\t220\n')

    def score(self, start, end, chrom='chr1'):
        bs, ps = self.bw.stat(), self.bed.stat()
        return peak_activity.window_score(str(self.bw), bs.st_size, bs.st_mtime_ns,
                                          str(self.bed), ps.st_size, ps.st_mtime_ns, chrom, start, end)

    def test_clipped_peak_union_integral(self):
        value, details = self.score(50, 125)
        self.assertEqual(value, 175.)  # 50 bases * 2 + 25 bases * 3
        self.assertEqual(details['peak_bp'], 75)
        self.assertEqual(details['covered_bp'], 75)
        self.assertEqual(self.score(0, 1000)[0], 310.)  # nonpeak signal excluded

    def test_empty_signal_and_empty_peaks_are_zero(self):
        value, details = self.score(200, 220)
        self.assertEqual(value, 0.)
        self.assertEqual(details['peak_bp'], 20)
        self.assertEqual(details['covered_bp'], 0)
        self.assertEqual(self.score(220, 300)[0], 0.)
        self.assertEqual(self.score(0, 100, 'chr2')[0], 0.)
        self.bed.write_text('')
        self.assertEqual(self.score(0, 100)[0], 0.)

    def test_clipping_alias_and_incompatible_chromosomes(self):
        value, details = self.score(-50, 50, '1')
        self.assertEqual(value, 60.)
        self.assertEqual(details['window_bp'], 50)
        self.assertEqual(details['status'], 'clipped')
        self.assertIsNone(self.score(0, 50, 'chrAbsent')[0])
        self.assertIsNone(self.score(1001, 1100)[0])
        self.bed.write_text('differentAssembly\t20\t80\n')
        value, details = self.score(0, 100)
        self.assertIsNone(value)
        self.assertEqual(details['status'], 'peak_assembly_mismatch')

    def test_union_matches_basewise_reference(self):
        rng = np.random.default_rng(42)
        intervals = []
        expected = np.zeros(1000, dtype=bool)
        for _ in range(80):
            a = int(rng.integers(0, 980)); b = min(1000, a + int(rng.integers(1, 40)))
            intervals.append((a, b)); expected[a:b] = True
        self.bed.write_text(''.join(f'chr1\t{a}\t{b}\n' for a, b in intervals))
        signal = np.r_[np.full(100, 2.), np.full(100, 3.), np.zeros(100), np.ones(700)]
        for lo, hi in [(0, 1000), (23, 186), (198, 357), (900, 999)]:
            value, details = self.score(lo, hi)
            self.assertAlmostEqual(value, float(np.sum(signal[lo:hi] * expected[lo:hi])))
            self.assertEqual(details['peak_bp'], int(np.sum(expected[lo:hi])))

    def test_strand_sample_median_units_and_missing_files(self):
        row = dict(species='test', locus='GENE|chr1:100-201:+', samples={'sample':None}, unit='CPM_signal_bp')
        with patch('server.SPECS', [('test', '', '', '', str(self.root), '')]):
            plus = peak_activity.score_records([row], 50)[0]
            minus = peak_activity.score_records([dict(row, locus='GENE|chr1:100-201:-')], 50)[0]
            self.assertEqual(plus['median'], 250.)
            self.assertEqual(minus['median'], 0.)
            self.assertEqual(plus['unit'], 'CPM_signal_bp')
            self.assertEqual(plus['region_start'], 50)
            self.assertEqual(minus['region_start'], 150)
            human = peak_activity.score_records([dict(row, unit='source_scaled_signal_bp')], 50)[0]
            self.assertEqual(human['unit'], 'source_scaled_signal_bp')
            self.bed.unlink()
            missing = peak_activity.score_records([row], 50)[0]
            self.assertIsNone(missing['median'])
            self.assertEqual(missing['window_details']['sample']['status'], 'peak_file_missing')


if __name__ == '__main__':
    unittest.main()
