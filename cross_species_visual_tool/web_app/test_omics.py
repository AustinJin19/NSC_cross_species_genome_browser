import csv
import io
import unittest
import omics

class LiverOmics(unittest.TestCase):
    def test_source_rna_value_and_tissue(self):
        data=omics.query('0610009b22rik')
        self.assertEqual(data['tissue'],'Liver')
        self.assertFalse(any(r['sample'] in ('BMR_Turk_H','BMR_Turk_Lu') for r in data['records']))
        self.assertFalse(any(r['sample']=='BMR_Turk_Lu' for r in omics.query('IGF1')['records']))
        self.assertNotIn('BMR_Turk_Lu',omics.export_csv('IGF1'))

    def test_source_protein_value(self):
        row=next(r for r in omics.query('0610009B22Rik')['records'] if r['assay']=='protein' and r['sample']=='Tombline_DM-1_21-129.raw')
        self.assertEqual(row['source_row'],11)
        self.assertEqual(row['raw'],6079856.5)
        self.assertAlmostEqual(row['normalized'],6.84447464141982)

    def test_missing_protein_not_zero(self):
        data=omics.query('HAS2')
        self.assertEqual(len([r for r in data['records'] if r['assay']=='rna']),4)
        self.assertEqual([r for r in data['records'] if r['assay']=='protein'],[])
        self.assertEqual(omics.query('NOT_A_REAL_GENE_123')['records'],[])

    def test_catalog_and_export(self):
        data=omics.query('ALB')
        for assay,n in [('rna',37),('protein',33)]:
            self.assertEqual(len([r for r in data['catalog'] if r['assay']==assay]),n)
            self.assertTrue(any(r['assay']==assay for r in data['records']))
        rows=list(csv.DictReader(io.StringIO(omics.export_csv('ALB'))))
        self.assertEqual(len(rows),len(data['records']))
        self.assertTrue(all(r['tissue']=='Liver' for r in rows))
        self.assertEqual(data['stale'],[])

    def test_lifespan_sources_and_export(self):
        import math
        data=omics.query('ALB')
        nmr=data['lifespans']['Heterocephalus glaber']
        self.assertEqual(nmr['years'],41)
        self.assertAlmostEqual(nmr['log10_years'],math.log10(41))
        self.assertIn('MLS_corrected',nmr['source'])
        self.assertIsNone(data['lifespans']['Eonycteris spelaea']['years'])
        self.assertIn('Anage.zip',data['lifespans']['Sus scrofa']['source'])
        rows=list(csv.DictReader(io.StringIO(omics.export_csv('ALB'))))
        for row in rows:
            expected=data['lifespans'][row['species']]['years']
            self.assertEqual(row['MLS_years'],'' if expected is None else str(expected))

    def test_igf1_includes_naked_mole_rat(self):
        from statistics import median
        data=omics.query('igf1')
        rows=[r for r in data['records'] if r['assay']=='rna' and r['species']=='Heterocephalus glaber']
        self.assertEqual(len(rows),5)
        self.assertEqual({r['sample'] for r in rows},{'SRR17216157','SRR17216166','SRR17216173','SRR17216182','SRR17216183'})
        self.assertAlmostEqual(median(r['log'] for r in rows),2.28681316777468)
        self.assertAlmostEqual(median(r['raw'] for r in rows),192.55891)
        self.assertTrue(all(r['normalized'] is None for r in rows))
        self.assertEqual(data['lifespans']['Heterocephalus glaber']['years'],41)

    def test_cross_species_bmr_liver(self):
        data=omics.query('IGF1')
        rows=[r for r in data['records'] if r['assay']=='rna' and r['species']=='Nannospalax galili']
        self.assertEqual(len(rows),1)
        row=rows[0]
        self.assertEqual(row['sample'],'BMR_Turk_Liv')
        self.assertEqual(row['source_row'],222603)
        self.assertAlmostEqual(row['raw'],4.110711)
        self.assertAlmostEqual(row['log'],0.70848132320648)
        self.assertIsNone(row['normalized'])
        self.assertEqual(data['sources']['rna']['file'],'CrossSpecies_Count.xlsx')
        self.assertNotIn('BMR_Turk_Lu',omics.export_csv('IGF1'))

    def test_mouse_supplement_provenance(self):
        data=omics.query('IGF1')
        rows=[r for r in data['records'] if r['assay']=='rna' and r['species']=='Mus musculus']
        self.assertEqual(len(rows),5)
        self.assertTrue(all(r['workbook']=='rnaLongTPM.xlsx' and r['normalized'] is None for r in rows))
        row=next(r for r in rows if r['sample']=='SRR17216343')
        self.assertAlmostEqual(row['raw'],893.091634)
        self.assertAlmostEqual(row['log'],2.95138203121333)
        exported=list(csv.DictReader(io.StringIO(omics.export_csv('IGF1'))))
        self.assertEqual(next(r for r in exported if r['sample']=='SRR17216343')['workbook'],'rnaLongTPM.xlsx')
        self.assertEqual(next(r for r in exported if r['sample']=='BMR_Turk_Liv')['workbook'],'CrossSpecies_Count.xlsx')
