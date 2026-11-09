"""Integration checks against the installed human CPM tracks and hg38 models."""
import sqlite3
import statistics
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pyBigWig

import activity
import server


class HumanATAC(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        server.load_catalog()
        cls.ids = {'SRR28763537', 'SRR28763586', 'SRR28763842', 'SRR28763939'}

    def test_catalog_and_base_resolution_signal(self):
        catalog = next(s for s in server.CATALOG if s['id'] == 'Human')
        self.assertEqual(set(catalog['samples']), self.ids)
        self.assertEqual(catalog['default'], 'SRR28763586')
        self.assertEqual(server.default_sample('Human'), 'SRR28763586')
        for sample, path in server.SAMPLES['Human'].items():
            self.assertEqual(path.parent.name, 'Homo_sapiens')
            self.assertTrue(server.peak_path(path).is_file())
            panel = server.get_panel('Human', 'SIRT6', 50, 0, sample, 101)
            self.assertEqual(panel['signalUnit'], 'CPM')
            self.assertEqual(panel['tssMethod'], 'MANE_Select')
            self.assertTrue(panel['peakAvailable'])
            self.assertTrue(panel['peaks'])
            with pyBigWig.open(str(path)) as bw:
                expected = np.nan_to_num(bw.values(panel['chrom'], panel['from']-1, panel['to']))
            if panel['strand'] == '-': expected = expected[::-1]
            np.testing.assert_allclose([p['value'] for p in panel['signal']], expected)

    def test_live_activity_matches_new_bigwigs_and_peaks(self):
        # Isolate real human records to keep this integration check fast.
        with sqlite3.connect(activity.DB) as con:
            rows = con.execute("SELECT * FROM activity WHERE species='Human' AND symbol_key='SIRT6'").fetchall()
        self.assertTrue(rows)
        with tempfile.TemporaryDirectory() as tmp:
            db = Path(tmp)/'activity.sqlite'
            (Path(tmp)/'species_mls.json').write_text('{}')
            with sqlite3.connect(db) as con:
                con.execute('CREATE TABLE activity (species TEXT,symbol TEXT,symbol_key TEXT,locus TEXT,median REAL,unit TEXT,samples TEXT)')
                con.executemany('INSERT INTO activity VALUES (?,?,?,?,?,?,?)', rows)
            with patch.object(activity, 'DB', db):
                for strategy in ['body', 'tss2kb', 'peaks50kb', 'peak_signal50kb']:
                    for row in activity.query('SIRT6', strategy)['records']:
                        self.assertEqual(set(row['samples']), {s+'.cpm' for s in self.ids})
                        chrom = row['locus'].rsplit('|',1)[1].rsplit(':',2)[0]
                        lo, hi = row['region_start'], row['region_end']
                        expected_values = []
                        for sample, value in row['samples'].items():
                            path = server.SAMPLES['Human'][sample.removesuffix('.cpm')]
                            with pyBigWig.open(str(path)) as bw:
                                signal = np.nan_to_num(bw.values(chrom, lo, hi))
                            if strategy.startswith('peak'):
                                mask = np.zeros(hi-lo, dtype=bool)
                                count = 0
                                with server.peak_path(path).open() as f:
                                    for line in f:
                                        c, a, b, *_ = line.split('\t')
                                        a, b = int(a), int(b)
                                        if c == chrom and a < hi and b > lo:
                                            count += 1
                                            mask[max(a,lo)-lo:min(b,hi)-lo] = True
                                expected = count if strategy.startswith('peaks') else float(signal[mask].sum())
                            else:
                                expected = float(signal.sum())
                            # bigWig stores float32; prefix-sum integration and
                            # basewise summation can differ at rounding precision.
                            np.testing.assert_allclose(value, expected, rtol=1e-7, atol=1e-6)
                            expected_values.append(expected)
                        np.testing.assert_allclose(row['median'], statistics.median(expected_values), rtol=1e-7, atol=1e-6)
                        self.assertEqual(row['unit'], 'peak_count' if strategy.startswith('peaks') else 'CPM_signal_bp')
                        self.assertTrue(row['region_eligible'])


if __name__ == '__main__':
    unittest.main()
