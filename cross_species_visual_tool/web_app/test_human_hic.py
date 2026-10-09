"""Integration checks against the optional local 4DN matrix."""
import json
import unittest
import numpy as np
import human_hic


@unittest.skipUnless(human_hic.PATH.exists(), 'Local human Hi-C file not installed')
class HumanHiC(unittest.TestCase):
    def test_values_and_bin_origin_match_reader(self):
        import hictkpy
        lo, hi = 102430563, 102530563
        result = human_hic.query('chr12', lo, hi)
        self.assertTrue(result['available'], result)
        self.assertEqual(result['origin'], 102430000)
        direct = hictkpy.File(str(human_hic.PATH), result['resolution'])
        m = direct.fetch(f'12:{lo-1}-{hi}', normalization=result['normalization']).to_numpy()
        for i,j,value in result['cells']:
            if np.isfinite(m[i,j]):
                self.assertAlmostEqual(value, float(m[i,j]))
            else:
                self.assertIsNone(value)
        json.dumps(result, allow_nan=False)

    def test_chromosome_bounds_and_missing_chromosome(self):
        result = human_hic.query('chr1', -5000, 5000)
        self.assertTrue(result['available'])
        self.assertEqual(result['origin'], 0)
        self.assertFalse(human_hic.query('chr1', 300000000, 300001000)['available'])
        self.assertFalse(human_hic.query('chrMissing', 1, 1000)['available'])

    def test_whole_chromosome_is_bounded(self):
        result = human_hic.query('chr1', 1, 248956422)
        self.assertTrue(result['available'])
        self.assertLessEqual(result['bins'], 161)
        self.assertLessEqual(len(result['cells']), 13041)
