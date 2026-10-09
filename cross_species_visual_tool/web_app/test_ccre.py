import unittest
import ccre
import server
class CCRE(unittest.TestCase):
 def test_boundaries(self):
  def has(pos):return any(r['name']=='EH38E1638326' for r in ccre.query('chr12',pos,pos)['regions'])
  self.assertFalse(has(102430773));self.assertTrue(has(102430774))
  self.assertTrue(has(102431042));self.assertFalse(has(102431043))
  self.assertEqual(ccre.query('missing',1,10)['regions'],[])
  self.assertEqual(ccre.query('chr12',-100,-1)['regions'],[])
 def test_human_only(self):
  server.load_catalog()
  for p in server.browse({'gene':'IGF1','window':50000})['panels']:
   if p['id']=='Human':
    self.assertTrue(p['ccre']['available']); self.assertGreater(len(p['ccre']['regions']),0)
    self.assertTrue(all(r['start']<=p['to'] and r['end']>=p['from'] for r in p['ccre']['regions']))
   else:self.assertIsNone(p['ccre'])
