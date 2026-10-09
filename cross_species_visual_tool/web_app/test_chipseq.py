import unittest
import server
import chipseq
class ChipRegions(unittest.TestCase):
 def test_bed_boundary(self):
  # First H3K4me3 interval is BED chr2:53131305-53132813.
  t=lambda lo,hi:next(t for t in chipseq.query('chr2',lo,hi) if t['name']=='H3K4me3')['regions']
  self.assertFalse(any(r['name']=='hetGlaH3K4me31' for r in t(53131305,53131305)))
  r=next(r for r in t(53131306,53131306) if r['name']=='hetGlaH3K4me31')
  self.assertEqual((r['start'],r['end']),(53131306,53132813))
  self.assertTrue(any(r['name']=='hetGlaH3K4me31' for r in t(53132813,53132813)))
  self.assertFalse(any(r['name']=='hetGlaH3K4me31' for r in t(53132814,53132814)))
 def test_species_tracks(self):
  server.load_catalog()
  for gene in ['IGF1','IGF1R']:
   p=server.get_panel('NMR',gene,50000,0,'NMR3Liver',50)
   self.assertEqual(len(p['chipTracks']),2)
   for t in p['chipTracks']:
    self.assertTrue(t['available'])
    for r in t['regions']:self.assertTrue(r['start']<=p['to'] and r['end']>=p['from'])
  for sid,assembly in chipseq.ASSEMBLIES.items():
   p=server.get_panel(sid,'IGF1',50000,0,None,50)
   self.assertEqual(len(p['chipTracks']),2)
   for track in p['chipTracks']:
    self.assertTrue(track['available'])
    self.assertEqual(track['assembly'],assembly)
    self.assertTrue(all(r['start']<=p['to'] and r['end']>=p['from'] for r in track['regions']))
   for name,color,_ in chipseq.TRACKS:
    mark='H3K27Ac' if name=='H3K27ac' else name
    row=(chipseq.LIFTED/f'{sid}_{mark}_{assembly}.bed').read_text().splitlines()[0].split()
    chrom,start,end,ident=row[:4]; start,end=int(start),int(end)
    def hits(pos):return next(t for t in chipseq.query(chrom,pos,pos,sid) if t['name']==name)['regions']
    self.assertFalse(any(r['name']==ident for r in hits(start)))
    self.assertTrue(any(r['name']==ident for r in hits(start+1)))
    self.assertTrue(any(r['name']==ident for r in hits(end)))
    self.assertFalse(any(r['name']==ident for r in hits(end+1)))
  for sid in ['BMR','DMR','Rabbit','ASM']:
   self.assertEqual(server.get_panel(sid,'IGF1',5000,0,None,50)['chipTracks'],[])
