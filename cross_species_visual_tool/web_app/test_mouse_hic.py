"""Integration checks against the optional local mouse liver ZT10 mcool matrix."""
import json
import unittest
import numpy as np
import mouse_hic


@unittest.skipUnless(mouse_hic.PATH.exists(), 'Local mouse Hi-C file not installed')
class MouseHiC(unittest.TestCase):
    def test_values_and_bin_origin_match_reader(self):
        import hictkpy
        lo, hi = 32011275, 32811275
        result = mouse_hic.query('chr9', lo, hi)
        self.assertTrue(result['available'], result)
        self.assertEqual(result['assembly'], 'mm10')
        self.assertEqual(result['resolution'], 10000)
        self.assertEqual(result['normalization'], 'weight')
        self.assertEqual(result['origin'], 32010000)
        direct = hictkpy.File(str(mouse_hic.PATH), result['resolution'])
        m = direct.fetch(f'chr9:{lo-1}-{hi}', normalization=result['normalization']).to_numpy()
        for i,j,value in result['cells']:
            if np.isfinite(m[i,j]):
                self.assertAlmostEqual(value, float(m[i,j]))
            else:
                self.assertIsNone(value)
        json.dumps(result, allow_nan=False)

    def test_chromosome_bounds_and_missing_chromosome(self):
        result = mouse_hic.query('chr1', -5000, 5000)
        self.assertTrue(result['available'])
        self.assertEqual(result['origin'], 0)
        self.assertFalse(mouse_hic.query('chr1', 300000000, 300001000)['available'])
        self.assertFalse(mouse_hic.query('chrMissing', 1, 1000)['available'])

    def test_whole_chromosome_is_bounded(self):
        result = mouse_hic.query('chr19', 1, 61431566)
        self.assertTrue(result['available'])
        self.assertLessEqual(result['bins'], 161)
        self.assertLessEqual(len(result['cells']), 13041)
