import gzip
import json
import sqlite3
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pyBigWig
import activity
import server
import sample_qc
from sample_qc import QC_POLICY, exclusions_for_policy, filter_records, sample_id


class SampleQC(unittest.TestCase):
    def setUp(self):
        self.excluded = exclusions_for_policy(dict(QC_POLICY, enabled=True))
        patcher = patch.object(sample_qc, 'EXCLUDED_SAMPLES', self.excluded)
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_filter_can_be_disabled_to_restore_previous_samples(self):
        restored = exclusions_for_policy(dict(QC_POLICY, enabled=False))
        self.assertEqual(restored, {'NMR1Liver'})
        with patch.object(sample_qc, 'EXCLUDED_SAMPLES', restored):
            rows = filter_records([
                dict(samples={'CH2Liver.cpm':1000, 'keep':10, 'keep2':30}, median=20),
                dict(samples={'BHW14B11Liver.cpm':200}, median=200),
            ])
            self.assertEqual(rows[0]['median'], 30)
            self.assertIn('CH2Liver.cpm', rows[0]['samples'])
            self.assertEqual(rows[1]['samples'], {'BHW14B11Liver.cpm':200})

    def test_identifiers_and_median_recalculation(self):
        for name in ['CH2Liver', 'CH2Liver.cpm', 'CH2Liver.cpm.bw', '/tmp/CH2Liver_scaled.bw']:
            self.assertEqual(sample_id(name), 'CH2Liver')
        original = dict(samples={'CH2Liver.cpm': 1000, 'keep': 10, 'keep2': 30, 'missing': None}, median=30)
        row = filter_records([original])[0]
        self.assertEqual(row['median'], 20)
        self.assertNotIn('CH2Liver.cpm', row['samples'])
        self.assertIn('CH2Liver.cpm', original['samples'])
        self.assertEqual(filter_records([dict(samples={'CH2Liver': 20}, median=20)]), [])
        self.assertTrue({'AGP1Liver','ASM3Liver','ASM4Liver'}.isdisjoint(self.excluded))
        self.assertEqual(len(self.excluded), 12)

    def test_all_catalog_samples_excluded(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder/'BHW14B11Liver.cpm.bw').touch()
            with patch.object(server, 'SPECS', [('Balaena_mysticetus','','','',tmp,'BHW')]), \
                 patch.object(server, 'EXCLUDED_SAMPLES', self.excluded), \
                 patch.object(server, 'CATALOG', []), patch.object(server, 'SAMPLES', {}), patch.object(server, 'GENES', {}):
                server.load_catalog()
                self.assertEqual(server.CATALOG, [])
                self.assertEqual(server.SAMPLES, {})
                self.assertEqual(server.GENES, {})

    def test_all_activity_strategies_and_locus_selection(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);db=root/'activity.sqlite'
            (root/'species_mls.json').write_text('{}')
            for sample,signal,peaks in [('keep',2.,[(1990,2000)]),('keep2',4.,[(1990,2000),(2005,2010)])]:
                with pyBigWig.open(str(root/(sample+'.bw')),'w') as bw:
                    bw.addHeader([('chr1',10000)])
                    bw.addEntries(['chr1'],[0],ends=[10000],values=[signal])
                (root/(sample+'_peaks.narrowPeak')).write_text(''.join(f'chr1\t{a}\t{b}\n' for a,b in peaks))
            with sqlite3.connect(db) as con:
                con.execute('CREATE TABLE activity (species TEXT,symbol TEXT,symbol_key TEXT,locus TEXT,median REAL,unit TEXT,samples TEXT)')
                con.executemany('INSERT INTO activity VALUES (?,?,?,?,?,?,?)',[
                    ('test','G','G','G|chr1:2000-2100:+',1000,'CPM_signal_bp',json.dumps({'CH2Liver.cpm':1000,'keep':10,'keep2':30})),
                    ('test','G','G','G|chr1:4000-4100:+',1,'CPM_signal_bp',json.dumps({'CH2Liver.cpm':1,'keep':40,'keep2':60})),
                    ('test','Gone','GONE','Gone|chr1:2000-2100:+',200,'CPM_signal_bp',json.dumps({'BHW14B11Liver.cpm':200})),
                ])
            with patch.object(activity,'DB',db), patch.object(server,'SPECS',[('test','','','',tmp,'')]):
                rows=activity.query('G')['records']
                self.assertEqual([r['median'] for r in rows],[20,50])
                self.assertTrue(max(rows,key=lambda r:r['median'])['locus'].endswith('4000-4100:+'))
                # A retired sample set is replaced by current catalog samples,
                # even when every historical sample was excluded by QC.
                refreshed=activity.query('Gone')['records'][0]
                self.assertEqual(set(refreshed['samples']), {'keep','keep2'})
                self.assertEqual(refreshed['median'], 6300.)
                for strategy in ['body','tss2kb','tss5kb','tss100kb','tss200kb','tss500kb',
                                 'peaks20kb','peaks50kb','peaks100kb','peaks200kb','peaks500kb',
                                 'peak_signal20kb','peak_signal50kb','peak_signal100kb','peak_signal200kb','peak_signal500kb']:
                    result=activity.query('G',strategy)['records']
                    for row in result:
                        self.assertEqual(set(row['samples']),{'keep','keep2'})
                        self.assertEqual(set(row['window_details']),{'keep','keep2'})
                    if strategy.startswith('peaks'):
                        self.assertEqual(result[0]['median'],1.5)
                    if strategy.startswith('peak_signal'):
                        self.assertEqual(result[0]['median'],40.)


if __name__=='__main__':
    unittest.main()
