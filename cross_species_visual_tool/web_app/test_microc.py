import unittest
from unittest.mock import patch
import microc

class LoopTests(unittest.TestCase):
    def test_anchor_overlap_and_boundaries(self):
        row=dict(chrom1='chr1',start1=101,end1=200,chrom2='chr1',start2=501,end2=600,name='x',score=2)
        with patch.object(microc,'read_loops',return_value=[row]):
            self.assertEqual(microc.query('x','chr1',100,100)['total'],0)
            self.assertEqual(microc.query('x','chr1',200,200)['total'],1)
            self.assertEqual(microc.query('x','chr1',201,500)['total'],0)
            self.assertEqual(microc.query('x','chr1',501,501)['total'],1)
            self.assertEqual(microc.query('x','chr1',1,700)['total'],1)
            self.assertEqual(microc.query('x','chr2',1,700)['total'],0)
    def test_samples_and_limit(self):
        for sample in ('NMR3Liver','NMR4Liver','NMR8Liver'):
            rows=microc.read_loops(sample)
            self.assertTrue(rows)
            self.assertTrue(all(r['name'].startswith(sample+'_') for r in rows))
            hit=microc.query(sample,'chr3',1,200000000,limit=2)
            self.assertGreater(hit['total'],2)
            self.assertEqual(len(hit['regions']),2)
            self.assertGreaterEqual(hit['regions'][0]['score'],hit['regions'][1]['score'])
        self.assertFalse(microc.query('NMR1Liver','chr3',1,100)['available'])

if __name__=='__main__': unittest.main()
