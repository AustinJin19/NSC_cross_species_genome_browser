"""Check name resolution, paralog boundaries and real-panel consistency."""
import csv
import io
import sqlite3
import unittest
from unittest.mock import patch
import gene_names


class ResolutionRules(unittest.TestCase):
    def test_approved_symbol_is_not_an_alias(self):
        resolver = gene_names.GeneNames(dict(approved=['A','B'], aliases={'A':'B','OLD':'B'}))
        self.assertEqual(resolver.canonical('a'), 'A')
        self.assertEqual(resolver.canonical('old'), 'B')

    def test_ambiguous_aliases_are_not_combined(self):
        resolver = gene_names.GeneNames(dict(approved=['A'], aliases={'OLD1':'A','OLD2':'A'}))
        self.assertEqual(resolver.resolve('A','X',{'OLD1','OLD2'})['method'], 'ambiguous')
        self.assertEqual(resolver.resolve('A','X',{'A','OLD1'})['symbol'], 'A')

    def test_one_to_many_is_not_resolved_by_missing_annotation(self):
        resolver = gene_names.GeneNames(dict(approved=['A'], orthologs={'Mouse':{'A':[
            dict(symbol='B',support=['Ensembl'],reciprocal_unique=False),
            dict(symbol='C',support=['Ensembl'],reciprocal_unique=False)]}}))
        self.assertEqual(resolver.resolve('A','Mouse',{'B'})['method'], 'ambiguous')
        self.assertIsNone(resolver.resolve('A','Mouse',{'B'})['symbol'])

    def test_hbe1_does_not_include_other_globins(self):
        resolver=gene_names.resolver()
        m=resolver.resolve('HBE1','Mouse',{'HBB-Y','HBB-BH1','HBB-BS'})
        self.assertEqual(m['symbol'],'HBB-Y')
        self.assertEqual(m['method'],'ortholog')
        self.assertEqual(resolver.canonical('Hbb-y'),'HBE1')
        self.assertIsNone(resolver.resolve('HBE1','Mouse',{'HBB-BH1','HBB-BS'})['symbol'])

    def test_projected_mapping_is_scoped_and_labelled(self):
        self.assertEqual(gene_names.resolve('HBE1','Rabbit',{'HBB-Y'})['method'],'projected_ortholog')
        self.assertEqual(gene_names.resolve('HBE1','ASM',{'HBB-Y'})['method'],'projected_ortholog')
        self.assertIsNone(gene_names.resolve('HBE1','BMR',{'HBB-Y'})['symbol'])
        self.assertIsNone(gene_names.resolve('HBE1','Rabbit',{'HBB-Y'},False)['symbol'])


class RealDataMatching(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import server
        server.load_catalog()

    def test_counts_and_both_search_directions(self):
        import server
        a=server.annotation_counts('HBE1');b=server.annotation_counts('HBB-Y')
        self.assertEqual(a['species'],b['species'])
        exact=sum('HBE1' in g for g in server.GENES.values())
        self.assertEqual(a['species'],exact+3)
        self.assertEqual(a['mappedSpecies'],3)
        for query in ['HBE1','HBB-Y']:
            data=server.browse(dict(gene=query,species=['Mouse','Rabbit','ASM'],window=50,bins=10))
            self.assertEqual(data['gene'],'HBE1')
            for p in data['panels']:
                self.assertNotIn('error',p)
                self.assertEqual(p['matchedGene'],'Hbb-y')
                self.assertEqual(p['models'][0]['gname'].upper(),'HBB-Y')
            exported=list(csv.DictReader(io.StringIO(server.export_csv(data))))
            self.assertTrue(all(r['annotation_gene']=='Hbb-y' for r in exported))

    def test_activity_and_peak_strategies_use_same_gene(self):
        import activity
        for strategy in ['body','tss5kb','peaks20kb','peak_signal20kb']:
            with self.subTest(strategy=strategy):
                a=activity.query('HBE1',strategy);b=activity.query('Hbb-y',strategy)
                self.assertEqual([(r['locus'],r['median']) for r in a['records']],
                                 [(r['locus'],r['median']) for r in b['records']])
                for sid in ['Mouse','Rabbit','ASM']:
                    rows=[r for r in a['records'] if r['species']==sid]
                    self.assertTrue(rows)
                    self.assertTrue(all(r['symbol_key']=='HBB-Y' for r in rows))
                    self.assertTrue(all(r['median'] is not None for r in rows))

    def test_unrelated_activity_genes_preserve_loci(self):
        import activity
        for symbol in ['SIRT6','HAS2']:
            with sqlite3.connect(activity.DB) as con:
                old=set(con.execute('SELECT species,locus FROM activity WHERE symbol_key=? UNION SELECT species,locus FROM activity_additions WHERE symbol_key=?',(symbol,symbol)))
            current={(r['species'],r['locus']) for r in activity.query(symbol,'tss5kb')['records']}
            self.assertEqual(current,old)

    def test_omics_source_values_and_symbols_preserved(self):
        import omics
        with omics.connect() as con:
            expected=[dict(r) for r in con.execute('SELECT * FROM measurements WHERE symbol="SIRT6" AND sample!="BMR_Turk_Lu"')]
        result=omics.query('SIRT6')
        actual=[{k:v for k,v in r.items() if k!='gene_match'} for r in result['records']]
        order=lambda r:(r['assay'],r['species'],r['sample'],r['source_row'])
        self.assertEqual(sorted(actual,key=order),sorted(expected,key=order))
        a=omics.query('HBE1');b=omics.query('Hbb-y')
        self.assertEqual(a['records'],b['records'] if a['requested_gene']==b['requested_gene'] else
                         [dict(r,gene_match=dict(r['gene_match'],query='HBE1')) for r in b['records']])


if __name__=='__main__':
    unittest.main()
