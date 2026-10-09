"""Integration checks against the bundled real tracks (no web server needed)."""
import unittest
import csv
import io
import numpy as np
import pyBigWig
import server

class GenomicQueries(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        server.load_catalog()

    def test_all_species_and_peak_overlap(self):
        result = server.browse({'gene':'has2'})
        self.assertEqual(len(result['panels']), 9)
        for p in result['panels']:
            self.assertNotIn('error', p)
            self.assertTrue(any((s['value'] or 0) > 0 for s in p['signal']))
            self.assertTrue(all(peak['start'] <= p['to'] and peak['end'] >= p['from'] for peak in p['peaks']))

    def test_added_assemblies_and_sample_coordinates(self):
        expected = {'Human': ('selected3', 'chr12', 3), 'Mouse': ('mouse_mm10', 'chr10', 3),
                    'Macaque': ('macaque', 'NC_133416.1', 4),
                    'Rabbit': ('Rabbit', 'OX463231.1', 2),
                    'ASM': ('Africa_spiny_mouse', 'JAULSH010000016.1', 2)}
        for sid, (folder, chrom, count) in expected.items():
            self.assertEqual(len(server.SAMPLES[sid]), count)
            for sample, path in server.SAMPLES[sid].items():
                self.assertEqual(path.resolve().parent.name, folder)
                p = server.get_panel(sid, 'IGF1', 50, 0, sample, 101)
                self.assertEqual(p['chrom'], chrom)
                with pyBigWig.open(str(path)) as bw:
                    raw = np.nan_to_num(np.array(bw.values(chrom, p['from']-1, p['to'])))
                if p['strand'] == '-': raw = raw[::-1]
                np.testing.assert_allclose([s['value'] for s in p['signal']], raw)
                self.assertEqual(p['models'][0]['gname'].upper(), 'IGF1')
                self.assertEqual(p['peakAvailable'], sid != 'ASM')
                self.assertEqual(p['signalUnit'], 'scaled' if sid == 'Human' else 'CPM')
                if sid == 'ASM':
                    self.assertIn('assembly', p['peakWarning'])
                    self.assertEqual(p['peaks'], [])

    def test_human_mane_tss_and_displayed_transcript(self):
        p = server.get_panel('Human','ALDH1A2',5000,0,'122',101)
        self.assertEqual(p['strand'], '-')
        self.assertEqual(p['tss'], 58065711)
        self.assertEqual(p['centre'], 58065711)
        self.assertEqual(p['tssTranscript'], 'ENST00000249750.9')
        self.assertEqual(p['tssMethod'], 'MANE_Select')
        self.assertEqual(p['models'][0]['tx'], p['tssTranscript'])
        self.assertEqual(p['models'][0]['end'], p['tss'])
        for g in server.GENES['Human'].values():
            for locus in g:
                self.assertEqual(locus['tss'], locus['end'] if locus['strand']=='-' else locus['start'])

    def test_base_resolution_matches_bigwig_both_strands(self):
        for sid in ['NMR', 'Mouse']:
            p = server.get_panel(sid, 'HAS2', 50, 0, None)
            with pyBigWig.open(str(server.SAMPLES[sid][p['sample']])) as bw:
                chrom = server.resolve_chrom(p['chrom'], bw.chroms())
                raw = np.nan_to_num(np.array(bw.values(chrom, p['from']-1, p['to'])))
            if p['strand'] == '-':
                raw = raw[::-1]
            np.testing.assert_allclose([s['value'] for s in p['signal']], raw)
            self.assertEqual(p['signal'][50]['rel'], 0)

    def test_pan_follows_transcription(self):
        for sid in ['NMR', 'Mouse']:
            p = server.get_panel(sid, 'HAS2', 50, 1000, None)
            self.assertEqual(p['centre']-p['tss'], -1000 if p['strand']=='-' else 1000)

    def test_missing_gene_and_empty_selection(self):
        self.assertEqual(server.browse({'species':[]})['panels'], [])
        self.assertTrue(all('error' in p for p in server.browse({'gene':'NOT_A_REAL_GENE_123'})['panels']))

    def test_off_chromosome_is_missing_not_zero(self):
        p = server.get_panel('NMR','HAS2',50,-100_000_000,None)
        self.assertTrue(all(s['value'] is None for s in p['signal']))
        self.assertEqual(p['peaks'],[])

    def test_wide_window_keeps_selected_gene(self):
        p = server.get_panel('Rat', 'TP53', 500000, 0, None)
        self.assertEqual(p['models'][0]['gname'].upper(), 'TP53')
        self.assertGreater(p['modelCount'], 5)

    def test_csv_preserves_coordinates_and_samples(self):
        result = server.browse({'species':['Mouse'], 'window':50, 'pan':1000})
        p = result['panels'][0]
        rows = list(csv.DictReader(io.StringIO(server.export_csv(result))))
        self.assertEqual(len(rows), 101)
        self.assertEqual(rows[50]['sample'], p['sample'])
        self.assertEqual(float(rows[50]['genomic_bin_center_1based']), p['centre'])
        self.assertEqual(float(rows[50]['relative_to_tss_bp']), 1000)

    def test_multiple_panels_keep_bins_samples_and_export_ids(self):
        specs = [dict(key='a',species='NMR',sample='NMR3Liver'),
                 dict(key='b',species='NMR',sample='NMR4Liver')]
        result = server.browse({'panels':specs, 'bins':123, 'window':6172.5})
        self.assertEqual([p['panelKey'] for p in result['panels']], ['a','b'])
        self.assertEqual([p['sample'] for p in result['panels']], ['NMR3Liver','NMR4Liver'])
        for p in result['panels']:
            self.assertEqual(len(p['signal']),123)
            self.assertEqual(p['to']-p['from'],12345)
        rows = list(csv.DictReader(io.StringIO(server.export_csv(result))))
        self.assertEqual(len(rows),246)
        self.assertEqual({r['panel_id'] for r in rows},{'a','b'})

    def test_annotation_counts_do_not_depend_on_visible_panels(self):
        counts = server.browse({'panels':[]})['annotation']
        self.assertEqual(counts['species'],9)
        self.assertEqual(counts['samples'],26)
        self.assertEqual(server.annotation_counts('not_a_gene')['samples'],0)
        tp53 = server.annotation_counts('TP53')
        self.assertEqual(tp53['species'],5)
        self.assertEqual(tp53['samples'],16)

    def test_invalid_bins_and_duplicate_keys(self):
        for bins in [0,9,5001]:
            with self.assertRaises(ValueError):
                server.browse({'bins':bins})
        with self.assertRaises(ValueError):
            server.browse({'panels':[dict(key='a',species='NMR'),dict(key='a',species='NMR')]})
        for window in [0,-1,0.75,float('nan'),float('inf')]:
            with self.assertRaises(ValueError):
                server.browse({'window':window})

    def test_unrestricted_window_and_exact_bins(self):
        for span in [1, 98, 2_000_000, 1_000_000_000]:
            p = server.browse({'gene':'IGF1','species':['Mouse'],'window':span/2,'bins':10})['panels'][0]
            self.assertEqual(p['to']-p['from'],span)
            self.assertEqual(len(p['signal']),min(10,span+1))
        p = server.get_panel('Mouse','IGF1',5000,0,None,73)
        with pyBigWig.open(str(server.SAMPLES['Mouse'][p['sample']])) as bw:
            values=np.nan_to_num(np.array(bw.values(p['chrom'],p['from']-1,p['to'])))
        if p['strand']=='-':values=values[::-1]
        edges=np.linspace(0,len(values),74,dtype=int)
        np.testing.assert_allclose([r['value'] for r in p['signal']], [values[a:b].mean() for a,b in zip(edges[:-1],edges[1:])],rtol=1e-6)

    def test_excluded_sample_not_available(self):
        nmr = next(c for c in server.CATALOG if c['id'] == 'NMR')
        self.assertEqual(nmr['default'], 'NMR3Liver')
        self.assertNotIn('NMR1Liver', nmr['samples'])
        with self.assertRaises(ValueError):
            server.browse({'panels':[dict(key='a',species='NMR',sample='NMR1Liver')]})

    def test_invalid_species_rejected(self):
        with self.assertRaises(ValueError):
            server.browse({'species':['wrong']})

    def test_replicates_are_distinct(self):
        a = server.get_panel('NMR','HAS2',5000,0,'NMR3Liver')
        b = server.get_panel('NMR','HAS2',5000,0,'NMR4Liver')
        self.assertNotEqual(a['signal'], b['signal'])

if __name__ == '__main__':
    unittest.main()
