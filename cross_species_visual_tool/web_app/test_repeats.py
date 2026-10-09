import unittest
import repeats
import server
class TERegions(unittest.TestCase):
 def test_bed_boundaries_and_classes(self):
  self.assertEqual(repeats.query('chr1',1651,1651)['regions'],[])
  r=repeats.query('chr1',1652,1666)['regions']
  self.assertEqual([(v['start'],v['end'],v['name']) for v in r],[(1652,1666,'LTR/ERVL-MaLR')])
  r=repeats.query('chr1',1667,1667)['regions']
  self.assertEqual(r[0]['name'],'SINE/Alu')
  self.assertEqual(repeats.query('chr1',3264,3305)['regions'],[])
  self.assertEqual(repeats.query('chr1___fragment_1',1,10000)['regions'],[])
 def test_nmr_only_and_negative_strand(self):
  server.load_catalog()
  p=server.get_panel('NMR','IGF1R',50000,0,'NMR3Liver',50)
  self.assertEqual(p['strand'],'-')
  self.assertTrue(p['teRepeats']['available'])
  self.assertTrue(p['teRepeats']['regions'])
  self.assertTrue(all(r['start']<=p['to'] and r['end']>=p['from'] for r in p['teRepeats']['regions']))
  self.assertIsNone(server.get_panel('Mouse','IGF1',5000,0,'WtMice11Liver',50)['teRepeats'])
