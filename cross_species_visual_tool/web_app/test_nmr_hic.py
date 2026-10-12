"""Integration checks against the optional local NMR hypothalamus mcool matrix."""
import json
import unittest
import numpy as np
import nmr_hic


@unittest.skipUnless(nmr_hic.PATH.exists(), 'Local NMR Hi-C file not installed')
class NMRHiC(unittest.TestCase):
    def test_values_and_bin_origin_match_reader(self):
        import hictkpy
        lo, hi = 32011275, 32811275
        result = nmr_hic.query('chr9', lo, hi)
        self.assertTrue(result['available'], result)
        self.assertEqual(result['assembly'], 'mHetGlaV3')
        self.assertEqual(result['resolution'], 10000)
        self.assertEqual(result['normalization'], 'weight')
        self.assertEqual(result['origin'], 32010000)
        direct = hictkpy.File(str(nmr_hic.PATH), result['resolution'])
        m = direct.fetch(f'chr9:{lo-1}-{hi}', normalization=result['normalization']).to_numpy()
        for i,j,value in result['cells']:
            if np.isfinite(m[i,j]):
                self.assertAlmostEqual(value, float(m[i,j]))
            else:
                self.assertIsNone(value)
        json.dumps(result, allow_nan=False)

    def test_chromosome_bounds_and_missing_chromosome(self):
        result = nmr_hic.query('chr1', -5000, 5000)
        self.assertTrue(result['available'])
        self.assertEqual(result['origin'], 0)
        self.assertFalse(nmr_hic.query('chr1', 300000000, 300001000)['available'])
        self.assertFalse(nmr_hic.query('chrMissing', 1, 1000)['available'])

    def test_whole_chromosome_is_bounded(self):
        result = nmr_hic.query('chr1', 1, 134629806)
        self.assertTrue(result['available'])
        self.assertLessEqual(result['bins'], 161)
        self.assertLessEqual(len(result['cells']), 13041)
