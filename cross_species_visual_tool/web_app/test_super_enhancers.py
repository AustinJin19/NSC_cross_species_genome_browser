import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import super_enhancers as se


class SuperEnhancers(unittest.TestCase):
    def test_bed_boundaries_long_spans_and_candidate_subset(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(se, 'DATA', Path(tmp)):
            root = Path(tmp)
            (root / 'human_hg38_liver_super_enhancers.bed').write_text(
                'chr1\t100\t200\tfirst\nchr1\t110\t120\tinside\nchr1\t200\t210\tnext\n')
            (root / 'human_hg38_candidate_liver_specific_super_enhancers.bed').write_text(
                'chr1\t100\t200\tfirst\n')
            (root / 'metadata.json').write_text(json.dumps({'species': {'Human': {
                'regions': {'first': {'nLiverDatasets': 2, 'liverDatasets': ['a', 'b']}}}}}))
            self.assertEqual(se.query('Human', 'chr1', 100, 100)['regions'], [])
            first = se.query('Human', 'chr1', 101, 101)['regions'][0]
            self.assertEqual((first['start'], first['end']), (101, 200))
            self.assertTrue(first['candidateLiverSpecific'])
            self.assertEqual(first['nLiverDatasets'], 2)
            self.assertEqual([r['name'] for r in se.query('Human', 'chr1', 200, 200)['regions']], ['first'])
            self.assertEqual([r['name'] for r in se.query('Human', 'chr1', 201, 201)['regions']], ['next'])
            self.assertEqual(se.query('Human', 'chr2', 1, 500)['regions'], [])
            self.assertEqual(se.query('Human', 'chr1', -50, -1)['regions'], [])
            self.assertEqual(se.query('Human', 'chr1', 10, 1)['regions'], [])

    def test_unavailable_files_and_other_species(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(se, 'DATA', Path(tmp)):
            self.assertIsNone(se.query('NMR', 'chr1', 1, 500))
            self.assertFalse(se.query('Human', 'chr1', 1, 500)['available'])
            (Path(tmp) / 'mouse_mm10_liver_super_enhancers.bed').write_text('chr1\t0\t10\tmouse\n')
            result = se.query('Mouse', 'chr1', 1, 1)
            self.assertTrue(result['available'])
            self.assertFalse(result['candidateAvailable'])
            self.assertEqual(result['regions'][0]['start'], 1)
            self.assertEqual(result['assembly'], 'mm10')

    @unittest.skipUnless((se.DATA / 'metadata.json').exists(), 'dbSUPER data not installed')
    def test_installed_data_and_browser_species(self):
        import server
        server.load_catalog()
        for sid, count, candidates in [('Human', 2252, 198), ('Mouse', 857, 241)]:
            result = se.query(sid, 'chr1', 1, 300_000_000)
            self.assertEqual((result['total'], result['candidateTotal']), (count, candidates))
            self.assertTrue(result['regions'])
            p = server.get_panel(sid, 'ERBB2', 250000, 0, None, 10)
            self.assertNotIn('error', p)
            self.assertEqual(p['superEnhancers']['assembly'], 'hg38' if sid == 'Human' else 'mm10')
            self.assertTrue(all(r['start'] <= p['to'] and r['end'] >= p['from']
                                for r in p['superEnhancers']['regions']))
        self.assertIsNone(server.get_panel('NMR', 'ERBB2', 5000, 0, None, 10)['superEnhancers'])


if __name__ == '__main__':
    unittest.main()
