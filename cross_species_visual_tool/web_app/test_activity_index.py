import csv, gzip, json, sqlite3, tempfile, unittest
from pathlib import Path
from unittest.mock import patch
import activity
import pyBigWig
class ActivityIndexTest(unittest.TestCase):
    def test_matching_and_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            (root/'manifest.json').write_text(json.dumps({'species':[{'species':'NMR'}]}))
            (root/'species_mls.json').write_text('{}')
            with gzip.open(root/'NMR_gene_activity_matrix.tsv.gz','wt') as f:
                w=csv.writer(f,delimiter='\t')
                w.writerow(['locus_id','gene_symbol','sample','species_median','valid_samples','unit'])
                w.writerows([['Gene|chr1:2000-2100:+','Gene',0,0,1,'CPM_signal_bp'],['GENE|chr1:0-100:+','GENE',2,2,1,'CPM_signal_bp'],['Missing|chr1:500-600:+','Missing','NA','NA',0,'CPM_signal_bp']])
            with pyBigWig.open(str(root/'sample.bw'),'w') as bw:
                bw.addHeader([('chr1',10000)])
                bw.addEntries(['chr1'],[0],ends=[10000],values=[1.])
            db=root/'activity.sqlite';activity.build(root,db)
            with patch.object(activity,'DB',db), patch('server.SPECS',[('NMR','','','',str(root),'')]):
                rows=activity.query('gene')['records']
                self.assertEqual(len(rows),2)
                self.assertTrue(all(r['ambiguous'] for r in rows))
                rows.sort(key=lambda r:r['median'])
                self.assertEqual(rows[0]['samples']['sample'],0)
                self.assertTrue(rows[0]['region_eligible'])
                self.assertFalse(rows[1]['region_eligible'])
                self.assertEqual(rows[1]['median'],2)  # audit score preserved
                self.assertAlmostEqual(rows[1]['available_fraction'],100/2100)
                self.assertIsNone(activity.query('Missing')['records'][0]['median'])
                self.assertEqual(activity.query('absent')['records'],[])
                with self.assertRaises(ValueError): activity.query('')

    def test_replaced_samples_are_recomputed_for_every_strategy(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'species_mls.json').write_text('{}')
            for sample, signal in [('SRR1', 2.), ('SRR2', 4.)]:
                with pyBigWig.open(str(root/f'{sample}.cpm.bw'), 'w') as bw:
                    bw.addHeader([('chr1', 10000)])
                    bw.addEntries(['chr1'], [0], ends=[10000], values=[signal])
                (root/f'{sample}_peaks.narrowPeak').write_text('chr1\t1990\t2010\n')
            db = root/'activity.sqlite'
            old_samples = json.dumps({'122_scaled': 999999})
            with sqlite3.connect(db) as con:
                con.execute('CREATE TABLE activity (species TEXT,symbol TEXT,symbol_key TEXT,locus TEXT,median REAL,unit TEXT,samples TEXT)')
                con.execute('INSERT INTO activity VALUES (?,?,?,?,?,?,?)',
                            ('Human','G','G','G|chr1:2000-2100:+',999999,'source_scaled_signal_bp',old_samples))
            with patch.object(activity, 'DB', db), patch('server.SPECS', [('Human','','','',str(root),'SRR')]):
                strategies = ['body','tss2kb','tss5kb','tss100kb','tss200kb','tss500kb']
                strategies += [f'{kind}{kb}kb' for kind in ['peaks','peak_signal'] for kb in [20,50,100,200,500]]
                for strategy in strategies:
                    with self.subTest(strategy=strategy):
                        row = activity.query('G', strategy)['records'][0]
                        self.assertEqual(set(row['samples']), {'SRR1.cpm','SRR2.cpm'})
                        self.assertEqual(row['unit'], 'peak_count' if strategy.startswith('peaks') else 'CPM_signal_bp')
                        if strategy == 'body': expected = 6300.
                        elif strategy.startswith('peaks'): expected = 1.
                        elif strategy.startswith('peak_signal'): expected = 60.
                        else: expected = min(10000, 2000 + int(strategy[3:-2])*1000) * 3.
                        self.assertAlmostEqual(row['median'], expected)
            with sqlite3.connect(db) as con:
                self.assertEqual(con.execute('SELECT samples FROM activity').fetchone()[0], old_samples)
