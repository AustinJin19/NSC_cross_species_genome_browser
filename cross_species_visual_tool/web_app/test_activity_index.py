import csv, gzip, json, tempfile, unittest
from pathlib import Path
from unittest.mock import patch
import activity
class ActivityIndexTest(unittest.TestCase):
    def test_matching_and_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            (root/'manifest.json').write_text(json.dumps({'species':[{'species':'NMR'}]}))
            (root/'species_mls.json').write_text('{}')
            with gzip.open(root/'NMR_gene_activity_matrix.tsv.gz','wt') as f:
                w=csv.writer(f,delimiter='\t')
                w.writerow(['locus_id','gene_symbol','sample','species_median','valid_samples','unit'])
                w.writerows([['a','Gene',0,0,1,'CPM_signal_bp'],['b','GENE',2,2,1,'CPM_signal_bp'],['c','Missing','NA','NA',0,'CPM_signal_bp']])
            db=root/'activity.sqlite';activity.build(root,db)
            with patch.object(activity,'DB',db):
                rows=activity.query('gene')['records']
                self.assertEqual(len(rows),2)
                self.assertTrue(all(r['ambiguous'] for r in rows))
                self.assertEqual(rows[0]['samples']['sample'],0)
                self.assertIsNone(activity.query('Missing')['records'][0]['median'])
                self.assertEqual(activity.query('absent')['records'],[])
                with self.assertRaises(ValueError): activity.query('')
