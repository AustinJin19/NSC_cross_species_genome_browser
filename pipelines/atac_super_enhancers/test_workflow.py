import contextlib
import gzip
import io
import tempfile
import unittest
from pathlib import Path

import workflow as w


class PipelineTests(unittest.TestCase):
    def test_support_counts_samples_not_duplicate_calls(self):
        beds = {'a':[('c',0,10),('c',5,15)], 'b':[('c',8,12)], 'c':[('c',12,20)]}
        segments = w.support_segments(beds)
        self.assertEqual(w.merged_support(segments,2), [('c',8,15)])
        self.assertFalse(any(len(names)>2 for _,_,_,names in segments))
        self.assertEqual(w.merged_support(segments,1), [('c',0,20)])

    def test_adjacent_samples_do_not_share_bases(self):
        seg = w.support_segments({'a':[('c',0,10)],'b':[('c',10,20)]})
        self.assertEqual(w.merged_support(seg,2), [])

    def test_exclusion_boundaries_and_nested_intervals(self):
        idx = w.interval_index([('c',2,100),('c',10,12)])
        self.assertTrue(w.overlaps(idx,'c',90,101))
        self.assertFalse(w.overlaps(idx,'c',100,101))
        self.assertFalse(w.overlaps(idx,'c',0,2))

    def test_peak_reference_mismatch_is_error(self):
        with tempfile.TemporaryDirectory() as d:
            src, out = Path(d)/'a.bed', Path(d)/'b.bed'
            src.write_text('wrong_assembly\t2\t10\n')
            with self.assertRaisesRegex(ValueError,'reference mismatch'):
                w.prepare_peaks(src,out,{'c':100},set())

    def test_duplicate_and_blacklisted_peaks(self):
        with tempfile.TemporaryDirectory() as d:
            src,out,mask = (Path(d)/x for x in ('a.bed','b.bed','mask.bed'))
            src.write_text('c\t0\t10\nc\t0\t10\nc\t15\t30\nchrM\t1\t20\n')
            mask.write_text('c\t20\t25\n')
            counts = w.prepare_peaks(src,out,{'c':100,'chrM':100},{'chrM'},mask)
            self.assertEqual(counts,dict(input=4,duplicate_interval=1,excluded_mask=1,excluded_contig=1,retained=1))
            self.assertEqual(list(w.bed_records(out)),[('c',0,10)])

    def test_annotation_header_first_transcript_and_strand(self):
        with tempfile.TemporaryDirectory() as d:
            src,out = Path(d)/'m.gz',Path(d)/'refseq.ucsc'
            with gzip.open(src,'wt') as f:
                f.write('scaffold_1\t10\t100\ttx1\tGENE1\t+\t10,80,\t20,100,\n')
                f.write('scaffold_1\t200\t300\ttx1\tGENE2\t-\t200,\t300,\n')
            genes,n = w.convert_models(src,out,{'scaffold_1':1000},set())
            self.assertEqual(n,2)
            self.assertEqual([r[0] for r in genes['scaffold_1']],[10,300])
            utils,_,_ = w.load_rose()
            start = utils.makeStartDict(str(out))
            self.assertEqual(len(start),2)
            self.assertEqual(start['TX000000001']['name'],'GENE1')

    def test_native_contig_round_trip_and_stitching(self):
        with tempfile.TemporaryDirectory() as d, contextlib.redirect_stdout(io.StringIO()):
            p = Path(d)
            (p/'a.bed').write_text('NC_123.1\t10000\t11000\ta\t0\nNC_123.1\t12000\t13000\tb\t0\nNC_123.1\t80000\t81000\tc\t0\n')
            (p/'m.tsv').write_text('NC_123.1\t50000\t51000\ttx\tG\t+\t50000\t51000\n')
            w.convert_models(p/'m.tsv',p/'refseq.ucsc',{'NC_123.1':100000},set())
            utils,main,_ = w.load_rose()
            utils.bedToGFF(str(p/'a.bed'),str(p/'a.gff'))
            result,_ = main.regionStitching(str(p/'a.gff'),12500,2500,str(p/'refseq.ucsc'),True)
            self.assertEqual(sorted((l.chr(),l.start(),l.end()) for l in result.getLoci()),
                             [('NC_123.1',10000,13000),('NC_123.1',80000,81000)])

    def test_vendor_integrity(self):
        self.assertEqual(w.verify_vendor()['commit'],'d845b98e8e8fa657b3d3e9af91653a8eea7449b4')


if __name__ == '__main__':
    unittest.main()
