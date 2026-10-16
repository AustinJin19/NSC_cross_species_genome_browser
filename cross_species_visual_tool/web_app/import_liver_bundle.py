"""Organize and validate the supplied liver bundle; build additional browser indexes.
Run with PYTHONPATH=cross_species_visual_tool/web_app/.deps python3 cross_species_visual_tool/web_app/import_liver_bundle.py
Existing browser species are preserved; raw source manifest remains unchanged.
"""
import csv, gzip, json, re, shutil
from pathlib import Path
from collections import defaultdict
import pyBigWig, pysam
ROOT=Path(__file__).resolve().parents[2]
BUNDLE=ROOT/'data/liver_atac_bundle'
INDEX=ROOT/'cross_species_visual_tool/index'
EXISTING={'Homo_sapiens','Heterocephalus_glaber','Nannospalax_galili','Fukomys_damarensis','Mus_musculus','Rattus_norvegicus','Macaca_mulatta','Sylvilagus_floridanus','Acomys_cahirinus'}
NAMES={'Hydrochoerus_hydrochaeris':'Capybara','Eonycteris_spelaea':'Cave nectar bat','Myocastor_coypus':'Nutria','Balaena_mysticetus':'Bowhead whale','Artibeus_jamaicensis':'Jamaican fruit bat','Castor_canadensis':'North American beaver','Cricetulus_griseus':'Chinese hamster','Chinchilla_lanigera':'Chinchilla','Octodon_degus':'Degu','Eptesicus_fuscus':'Big brown bat','Sciurus_niger':'Fox squirrel','Meriones_unguiculatus':'Mongolian gerbil','Mesocricetus_auratus':'Golden hamster','Molossus_molossus':'Pallas’s mastiff bat','Erethizon_dorsatum':'North American porcupine','Tamiasciurus_hudsonicus':'American red squirrel','Marmota_monax':'Woodchuck','Peromyscus_maniculatus':'Deer mouse','Microtus_pennsylvanicus':'Meadow vole','Cavia_porcellus':'Guinea pig','Sciurus_carolinensis':'Eastern gray squirrel'}
COLORS=['#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#DC0000','#7E6148','#B09C85']
def local_track_path(source, target, overrides):
    """Keep consolidated paths and deliberately removed samples on re-import."""
    relative = str(target.relative_to(ROOT))
    if relative in overrides:
        retained = overrides[relative]
        if retained is None:
            return '-'
        path = (ROOT / retained).resolve()
        if ROOT.resolve() not in path.parents or not path.is_file():
            raise ValueError('Missing or invalid retained bundle input: ' + retained)
        return retained
    move(source, target)
    return relative if target.exists() else '-'

def move(src,dst):
    dst.parent.mkdir(parents=True,exist_ok=True)
    if src.exists():
        if dst.exists():raise RuntimeError('Refusing overwrite: '+str(dst))
        shutil.move(str(src),str(dst))
def build(sp,ref,trackdir):
    sizes={f[0]:int(f[1]) for l in (ref/(sp+'.fa.gz.fai')).read_text().splitlines() if (f:=l.split('\t'))}
    tracks=sorted(trackdir.glob('*.bw')); chroms={}
    for p in tracks:
        with pyBigWig.open(str(p)) as bw:
            cs=bw.chroms()
            bad=[c for c,n in cs.items() if sizes.get(c)!=n]
            if bad:raise ValueError(f'{p.name}: sequence length/name mismatch: {bad[:5]}')
            chroms.update(cs)
    # Read indexed compressed FASTA directly, without creating expanded genome copies.
    with pysam.FastaFile(str(ref/(sp+'.fa.gz'))) as fa:
        c=next(iter(chroms));assert len(fa.fetch(c,0,min(20,chroms[c])))==min(20,chroms[c])
    for p in trackdir.glob('*narrowPeak*'):
        with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as handle:
            for l in handle:
                if not l.strip() or l.startswith(('#','track','browser')):continue
                f=l.split();a,b=int(f[1]),int(f[2])
                if f[0] not in sizes or not 0<=a<b<=sizes[f[0]]:raise ValueError(f'Invalid peak coordinates: {p.name}: {l[:80]}')
    txs={};skipped=0
    with gzip.open(ref/(sp+'.gtf.gz'),'rt') as f:
        for l in f:
            if l.startswith('#'):continue
            v=l.rstrip().split('\t')
            if len(v)!=9 or v[2] not in ('transcript','exon'):continue
            c,_,feature,a,b,_,strand,_,attrs=v;a,b=int(a)-1,int(b)
            if c not in chroms or not 0<=a<b<=chroms[c]:skipped+=1;continue
            attr=dict(re.findall(r'(\w+) "([^"]*)"',attrs));tx=attr.get('transcript_id')
            if not tx:continue
            # ORF annotation labels contain a homology name plus a source-species suffix.
            label=attr.get('annotation_label', '').rsplit('_',1)[0]
            symbol=attr.get('gene_name') or attr.get('gene') or label or (tx.split('#')[1] if '#' in tx else attr.get('gene_id'))
            if not symbol:continue
            r=txs.setdefault((c,strand,tx),[a,b,symbol,set()]);r[0]=min(a,r[0]);r[1]=max(b,r[1])
            if feature=='exon':r[3].add((a,b))
    rows=[]
    for (c,strand,tx),(a,b,symbol,exons) in txs.items():
        if exons:
            exons=sorted(exons);rows.append((c,a,b,tx,symbol,strand,','.join(str(e[0]) for e in exons)+',',','.join(str(e[1]) for e in exons)+','))
    if not rows:raise ValueError('No valid annotated transcripts for '+sp)
    plain=INDEX/(sp+'_models.tsv');plain.write_text(''.join('\t'.join(map(str,r))+'\n' for r in sorted(rows)))
    target=str(plain)+'.gz';pysam.tabix_compress(str(plain),target,force=True);pysam.tabix_index(target,preset='bed',force=True);plain.unlink()
    return {'transcripts':len(rows),'genes':len({r[4].upper() for r in rows}),'unsupported_features':skipped,'samples':len(tracks)}
def main():
    records=list(csv.DictReader((BUNDLE/'manifest.tsv').open(),delimiter='\t'));groups=defaultdict(list)
    override_path=BUNDLE/'path_overrides.json'
    overrides=json.loads(override_path.read_text()) if override_path.exists() else {}
    local=[]
    for r in records:groups[r['species']].append(r)
    specs=[];report=[]
    for i,(sp,rs) in enumerate(groups.items()):
        ref=ROOT/'reference_genomes/liver_atac_bundle'/sp;tracks=ROOT/'data/atac/liver_atac_bundle'/sp
        for p in list((BUNDLE/'genome').glob(sp+'.*'))+list((BUNDLE/'annotation').glob(sp+'.*')):move(p,ref/p.name)
        for r in rs:
            item=dict(r)
            for key,folder in [('bigwig','bigwig'),('narrowPeak','peaks')]:
                if r[key]!='-':
                    name=Path(r[key]).name
                    item[key]=local_track_path(BUNDLE/folder/name,tracks/name,overrides)
            item['fasta']=str((ref/(sp+'.fa.gz')).relative_to(ROOT));item['gtf']=str((ref/(sp+'.gtf.gz')).relative_to(ROOT));local.append(item)
        status={'species':sp}
        if sp in EXISTING:status['status']='existing browser sources retained; see manifest.local.tsv for consolidated paths'
        elif not list(tracks.glob('*.bw')):status['status']='not added: no bigWig signal supplied'
        else:
            try:
                status.update(build(sp,ref,tracks));status['status']='added'
                specs.append([sp,NAMES.get(sp,sp.replace('_',' ')),sp.replace('_',' '),COLORS[i%len(COLORS)],'../'+str(tracks.relative_to(ROOT)),''])
            except Exception as e:status['status']='not added: '+str(e)
        report.append(status);print(status,flush=True)
    with (BUNDLE/'manifest.local.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(local)
    (BUNDLE/'import_report.json').write_text(json.dumps(report,indent=2)+'\n')
    (ROOT/'cross_species_visual_tool/web_app/data/additional_species.json').write_text(json.dumps(specs,indent=2)+'\n')
if __name__=='__main__':main()
