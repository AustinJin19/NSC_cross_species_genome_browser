import gzip
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pyBigWig
from calculate_gene_activity import gene_regions, integrate, quantify

class GeneActivity(unittest.TestCase):
    def test_integral(self):
        values,cov=integrate([(10,20,2.),(25,30,4.)],np.array([0,12,20,26]),np.array([40,28,25,30]))
        np.testing.assert_allclose(values,[40,28,0,16])
        np.testing.assert_allclose(cov,[15,11,0,4])
    def test_regions(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'models.gz'
            with gzip.open(p,'wt') as f:
                f.write('chr1\t3000\t4000\ta\tPLUS\t+\nchr1\t3500\t5000\tb\tPLUS\t+\nchr1\t8000\t9000\tc\tPLUS\t+\nchr1\t200\t1000\td\tMINUS\t-\n')
            rows=gene_regions(p)
            plus=[r for r in rows if r['gene_symbol']=='PLUS']
            self.assertEqual(len(plus),2)
            self.assertEqual((plus[0]['region_start'],plus[0]['region_end']),(1000,5000))
            self.assertEqual(plus[0]['transcript_count'],2)
            minus=next(r for r in rows if r['gene_symbol']=='MINUS')
            self.assertEqual((minus['region_start'],minus['region_end']),(200,3000))
    def test_chunk_boundary(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'track.bw'
            with pyBigWig.open(str(p),'w') as bw:
                bw.addHeader([('chr1',10000000)])
                bw.addEntries(['chr1'],[4999900],ends=[5000100],values=[2.])
            rows=[dict(chrom='chr1',region_start=4999800,region_end=5000200)]
            fast=quantify(p,rows,chunked=True)
            self.assertEqual(fast[0][0],400)
            self.assertEqual(fast[1][0],200)
    def test_tss_window(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'models.gz'
            with gzip.open(p,'wt') as f:
                f.write('chr1\t300000\t400000\ta\tPLUS\t+\nchr1\t300000\t400000\tb\tMINUS\t-\n')
            rows=gene_regions(p,tss_flank=200000)
            by_name={r['gene_symbol']:r for r in rows}
            self.assertEqual((by_name['PLUS']['region_start'],by_name['PLUS']['region_end']),(100000,500000))
            self.assertEqual((by_name['MINUS']['region_start'],by_name['MINUS']['region_end']),(199999,599999))
    def test_bigwig(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'track.bw'
            with pyBigWig.open(str(p),'w') as bw:
                bw.addHeader([('chr1',100)])
                bw.addEntries(['chr1','chr1'],[10,25],ends=[20,30],values=[2.,4.])
            rows=[dict(chrom=c,region_start=a,region_end=b) for c,a,b in [('chr1',-10,110),('chr1',12,28),('missing',0,10),('chr1',40,50)]]
            v,c,w,status=quantify(p,rows)
            fast=quantify(p,rows,chunked=True)
            np.testing.assert_allclose(fast[0],v,equal_nan=True)
            np.testing.assert_allclose(fast[1],c,equal_nan=True)
            np.testing.assert_allclose(v[[0,1,3]],[40,28,0])
            self.assertTrue(np.isnan(v[2]))
            self.assertEqual(status,['clipped','ok','chromosome_missing','ok'])
            with pyBigWig.open(str(p)) as bw:
                self.assertAlmostEqual(v[1],bw.stats('chr1',12,28,type='sum',exact=True)[0])
