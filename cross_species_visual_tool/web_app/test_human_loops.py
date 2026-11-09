import gzip
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from build_human_liver_loops import parse_anchor, valid_mapping
import human_loops


class LiverLoops(unittest.TestCase):
    def test_coordinate_conversion_and_mapping_failures(self):
        self.assertEqual(parse_anchor('chr19:4147554-4186620'), ('chr19', 4147553, 4186620))
        sizes = {'chr19': 60_000_000}
        anchor = ('chr19', 4147556, 4186623, '+')
        self.assertEqual(valid_mapping([anchor], sizes, 39067), (anchor, 'mapped'))
        self.assertEqual(valid_mapping([], sizes, 39067)[1], 'unmapped')
        self.assertEqual(valid_mapping([anchor, anchor], sizes, 39067)[1], 'ambiguous')
        self.assertEqual(valid_mapping([('chr19', -1, 10, '+')], sizes, 11)[1], 'invalid_target_bounds')
        self.assertEqual(valid_mapping([('chr19', 1, 1000, '+')], sizes, 10)[1], 'extreme_length_change')

    def test_both_anchors_boundaries_deduplication_and_limit(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'test.bedpe.gz'
            with gzip.open(path, 'wt') as handle:
                handle.write('chr19\t100\t110\tchr19\t200\t210\tlower\t2\t.\t.\told1\told2\n')
                handle.write('chr19\t105\t115\tchr19\t205\t215\thigher\t3\t.\t.\told3\told4\n')
            with patch.object(human_loops, 'LOOPS', path):
                self.assertEqual(human_loops.query('chr19', 100, 100)['total'], 0)
                self.assertEqual(human_loops.query('chr19', 101, 101)['total'], 1)
                self.assertEqual(human_loops.query('chr19', 110, 110)['total'], 2)
                self.assertEqual(human_loops.query('chr19', 210, 210)['total'], 2)
                self.assertEqual(human_loops.query('chr19', 216, 216)['total'], 0)
                result = human_loops.query('chr19', 1, 1000, limit=1)
                self.assertEqual(result['total'], 2)  # not four despite both anchors in view
                self.assertEqual([r['name'] for r in result['regions']], ['higher'])
                self.assertEqual(human_loops.query('chr1', 1, 1000)['regions'], [])
            with patch.object(human_loops, 'LOOPS', Path(tmp) / 'absent.gz'):
                self.assertFalse(human_loops.query('chr19', 1, 1000)['available'])

    @unittest.skipUnless(human_loops.LOOPS.exists(), 'Converted data not installed')
    def test_real_sirt6_tss_contacts(self):
        result = human_loops.query('chr19', 4182563, 4182563)
        self.assertEqual(result['total'], 14)
        self.assertEqual(result['assembly'], 'hg38')
        self.assertEqual(result['regions'][0]['name'], 'GSE167200_Liver_013899')
        self.assertAlmostEqual(result['regions'][0]['score'], 6.7270731925964355)
        for row in result['regions']:
            i = 1 if row['start1'] <= 4182563 <= row['end1'] else 2
            self.assertEqual((row[f'start{i}'], row[f'end{i}']), (4147557, 4186623))


if __name__ == '__main__':
    unittest.main()
