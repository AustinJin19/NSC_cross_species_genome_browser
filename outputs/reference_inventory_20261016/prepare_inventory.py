import ast, csv, gzip, hashlib, json, pathlib, re

ROOT = pathlib.Path('/Users/zjin27/Downloads/CS_Genome_Browser')
OUT = ROOT/'outputs/reference_inventory_20261016'
evidence = json.loads((OUT/'local_evidence.json').read_text())
byfolder = {r['folder']:r for r in evidence}
public = {}
for p in (OUT/'ncbi').glob('*.json'):
    for r in json.loads(p.read_text()).get('reports',[]):
        if 'accession' in r:
            public[r['accession']] = r

original = [
 ('Human','Human','Homo sapiens','selected3','GCA_000001405.15'),
 ('Mouse','House mouse','Mus musculus','mouse_mm10','GCA_000001635.2'),
 ('NMR','Naked mole-rat','Heterocephalus glaber','Heterocephalus_glaber','GCA_964261345.1'),
 ('DMR','Damaraland mole-rat','Fukomys damarensis','Fukomys_damarensis','GCF_012274545.1'),
 ('BMR','Blind mole-rat','Nannospalax galili','Nannospalax_galili','GCF_000622305.1'),
 ('Rat','Norway rat','Rattus norvegicus','rat','GCF_036323735.1'),
 ('Macaque','Rhesus macaque','Macaca mulatta','macaque','GCF_049350105.2'),
 ('Rabbit','Eastern cottontail','Sylvilagus floridanus','Rabbit','GCA_949820135.1'),
 ('ASM','African spiny mouse','Acomys cahirinus','Africa_spiny_mouse','GCA_030555615.1'),
]
accessions = {
 'Hydrochoerus_hydrochaeris':'GCA_015741225.1',
 'Eonycteris_spelaea':'GCA_062009622.1',
 'Myocastor_coypus':'', 'Balaena_mysticetus':'',
 'Artibeus_jamaicensis':'GCF_021234435.1',
 'Castor_canadensis':'GCF_047511655.1',
 'Cricetulus_griseus':'GCF_003668045.3',
 'Chinchilla_lanigera':'GCF_000276665.1',
 'Octodon_degus':'GCF_000260255.1',
 'Bramus_lutescens':'GCA_001685075.1',
 'Eptesicus_fuscus':'GCF_027574615.1',
 'Sciurus_niger':'GCA_020740815.1',
 'Meriones_unguiculatus':'GCF_030254825.1',
 'Mesocricetus_auratus':'GCF_017639785.1',
 'Molossus_molossus':'GCF_014108415.1',
 'Erethizon_dorsatum':'GCA_028451465.1',
 'Tamiasciurus_hudsonicus':'GCA_057657825.1',
 'Marmota_monax':'GCF_021218885.2',
 'Peromyscus_maniculatus':'GCF_049852395.1',
 'Microtus_pennsylvanicus':'GCF_037038515.1',
 'Cavia_porcellus':'GCF_034190915.1',
 'Sciurus_carolinensis':'GCF_902686445.1',
 'Cricetomys_gambianus':'GCA_004027575.1',
 'Scalopus_aquaticus':'GCA_004024925.1',
 'Ondatra_zibethicus':'GCA_004026605.1',
 'Tamias_striatus':'',
}
species = [dict(id=i,common=c,scientific=s,folder='reference_genomes/'+f,accession=a,status='Active') for i,c,s,f,a in original]
for i,c,s,*_ in json.loads((ROOT/'cross_species_visual_tool/web_app/data/additional_species.json').read_text()):
    species.append(dict(id=i,common=c,scientific=s,folder='reference_genomes/liver_atac_bundle/'+i,accession=accessions[i],status='Active'))
for i,c in [('Cricetomys_gambianus','Gambian pouched rat'),('Scalopus_aquaticus','Eastern mole'),('Ondatra_zibethicus','Muskrat'),('Tamias_striatus','Eastern chipmunk')]:
    species.append(dict(id=i,common=c,scientific=i.replace('_',' '),folder='reference_genomes/liver_atac_bundle/'+i,accession=accessions[i],status='Not active: no bigWig'))

bedurl='https://github.com/alejandrogzi/bed2gtf'
togaurl='https://github.com/hillerlab/TOGA'
rows=[]
for s in species:
    ev=byfolder[s['folder']]; fa=ev['fasta'][0]
    an=next((a for a in ev['annotations'] if 'hetGla2_' not in a['path']),ev['annotations'][0])
    aid=s['accession']; rec=public.get(aid,{})
    name=rec.get('assembly_info',{}).get('assembly_name','Unresolved')
    genome_url=f'https://www.ncbi.nlm.nih.gov/datasets/genome/{aid}/' if aid else ''
    ann_url=rec.get('annotation_info',{}).get('report_url','')
    ann='NCBI RefSeq annotation (Gnomon / RefSeq features)'
    release='Not retained in local GTF'
    scope='Assembly record; public annotation report for this assembly. Local annotation release not verified.'
    notes='Assembly identified from local sequence accessions / FASTA build name and NCBI metadata. No full sequence checksum comparison.'
    first=an['features'][0]
    if '\tbed2gtf\t' in first:
        ann='Projected gene models converted to GTF with bed2gtf'
        ann_url=bedurl
        scope='Assembly record; annotation link is conversion software, not the local GTF download.'
        notes+=' Original prediction release and exact public annotation download were not recorded.'
    if s['id']=='Human':
        name='GRCh38 / hg38 (primary assembly)'
        ann='GENCODE comprehensive primary-assembly annotation'
        release='GENCODE v50 (Ensembl 116); header date 2026-04-08'
        genome_url='https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/'
        ann_url='https://www.gencodegenes.org/human/release_50.html'
        scope='hg38 assembly resource; exact GENCODE release with primary-assembly GTF download.'
        notes='Local FASTA has 194 sequences, consistent with the GRCh38 primary set. GENCODE v50 is a GRCh38.p14 release; local primary FASTA excludes patch/alternate contigs. Original FASTA download URL is not retained.'
    if s['id']=='Mouse':
        name='GRCm38 / mm10'
        ann='UCSC NCBI RefSeq gene models for mm10'
        release='ncbiRefSeq.2021-04-23 (GTF feature source)'
        genome_url='https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/'
        ann_url='https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.ncbiRefSeq.gtf.gz'
        scope='UCSC assembly resource and matching dated GTF download.'
        notes='Active browser uses Mouse_mm10_models.tsv.gz. Bundle mouse GTF is a different HAVANA/GENCODE annotation. mouse_legacy contains GRCm39 and is not used by the current Python browser.'
    if s['id']=='NMR':
        name='mHetGlaV3 (local chr-name distribution)'
        ann='Ensembl-derived gene models transferred with Liftoff; rtracklayer GTF export'
        release='Local export: rtracklayer 1.62.0, 2024-12-15; Ensembl release not recorded'
        genome_url='https://zenodo.org/records/14497788'
        ann_url='https://github.com/agshumate/Liftoff'
        scope='Assembly-associated Zenodo dataset; annotation link is transfer software. Exact local GTF distribution unverified.'
        notes='mHetGlaV3 accession supported by local project configuration and public dataset. Local names are chr1, chr2, etc. Exact Zenodo file version was not checksum-matched. Active GTF has ENSHGL identifiers, not the separate hetGla2 RefSeq liftOver file.'
    if s['id'] in ['DMR','BMR']:
        ann='NCBI RefSeq annotation'
        release='NCBI Annotation Release 102 (local GTF header)'
        scope='Assembly record and exact annotation release report.'
    if s['id']=='ASM':
        release='bed2gtf 1.9.3, 2026-03-29; upstream annotation release not recorded'
        notes='Assembly identified by JAULSH01 WGS project, isolate DEM-2023a and full indexed length. GTF uses scaffold_N aliases; browser builder maps these to JAULSH accessions using FASTA descriptions. This is Acomys cahirinus.'
    if s['id']=='Rabbit':
        release='bed2gtf 1.9.3, 2026-04-01; upstream annotation release not recorded'
        notes='Assembly supported by OX463223.1 / CATKVR01 and exact indexed assembly size. Original active GTF contains mouse ENSMUS projections; the retained bundle GTF contains human ENSG projections and is not the active index source.'
    if s['id']=='Hydrochoerus_hydrochaeris':
        ann='Local TOGA2-run projected models converted with bed2gtf'
        notes='Manifest labels assembly/run as toga2; this is not an assembly accession. Public HiC_scaffold_1 and _2 names and lengths match the local FASTA. Original scaffold names are preserved; full-file checksum not checked.'
    if s['id']=='Eonycteris_spelaea':
        ann='Local TOGA2-run projected models converted with bed2gtf'
        notes='All 305 local sequence names and lengths match the NCBI sequence report for GCA_062009622.1 using original names s1, s2, etc. Sequence bases were not checksum-compared. Manifest label toga2 is a pipeline label.'
    if s['id']=='Myocastor_coypus':
        name='Local Hi-C scaffold assembly (HiC_scaffold_N)'
        ann='Local TOGA2-run projected models converted with bed2gtf'
        notes='Manifest says toga2, not an assembly name. Local FASTA is 2,632,872,920 bp, 2,956 sequences. Exact public assembly accession/download was not recovered; do not substitute the different MyoCoy_v1_BIUU assembly.'
        scope='No verified assembly URL; annotation link is conversion software only.'
    if s['id']=='Balaena_mysticetus':
        name='HLbalMyst1A (local manifest label)'
        ann='Custom metassembly ORF gene models, BED12-to-GTF'
        release='HLbalMyst1A.metassembly.orf.gtf (upstream filename); version not recorded'
        ann_url=''
        scope='No exact public assembly or annotation URL recovered.'
        notes='Manifest explicitly points to HLbalMyst1A.fa and HLbalMyst1A.metassembly.orf.gtf. Local sequence names are HAP1_SUPER_N; transcript IDs include xORF and annotation_label. Do not equate with the older public HLbalMys1 assembly.'
    if s['id']=='Bramus_lutescens':
        notes='Bramus lutescens is represented in older resources as Ellobius lutescens. LOEQ01 identifies ASM168507v1. Local annotation uses projected human ENSG identifiers; exact prediction release is not retained.'
    if s['id']=='Peromyscus_maniculatus':
        notes='Public assembly organism is Peromyscus maniculatus bairdii. Local browser label is Peromyscus maniculatus. Local NC_134852.1 chromosome matches HU_Pman_BW_mat_3.1.'
    if s['id']=='Tamias_striatus':
        name='Local contig assembly (ptg identifiers)'
        ann='MAKER gene annotation (TST identifiers)'
        ann_url='https://www.yandell-lab.org/software/maker.html'
        notes='FASTA starts with ptg000001l; local assembly has 393 sequences. GTF source is maker and includes TST gene identifiers. Exact public assembly accession and annotation release were not recovered.'
        scope='No verified assembly URL; annotation link is MAKER software, not the local dataset.'
    if s['id']=='Erethizon_dorsatum':
        ann='NCBI GenBank submitted Gnomon gene models (O8D03 locus tags)'
        ann_url=genome_url
        scope='NCBI assembly page provides submitted annotation downloads. Local annotation release not verified.'
    if s['id']=='Tamiasciurus_hudsonicus':
        notes='Local CM172909.1 / JBXXPO01 identify mTamHud3.hap1. Indexed size is 5.72 Gb. Local GTF is bed2gtf with human ENSG identifiers; it is not assumed to be the submitted annotation currently shown by NCBI.'
    rows.append({**s,'assembly':name,'annotation':ann,'release':release,'assembly_url':genome_url,'annotation_url':ann_url,'link_scope':scope,'notes':notes,'fasta':fa['path'],'gtf':an['path'],'index':f'cross_species_visual_tool/index/{"Mouse_mm10" if s["id"]=="Mouse" else s["id"]}_models.tsv.gz' if s['status']=='Active' else '', 'sequence_count':fa['sequences'],'total_bp':fa['total_bp'],'public_annotation_release':rec.get('annotation_info',{}).get('name','')})

# File-level provenance also includes retained alternative references.
active_by_gtf={r['gtf']:r for r in rows}
scientific_by_folder={r['folder']:r['scientific'] for r in rows}
files=[]
for ev in evidence:
    for an in ev['annotations']:
        fa=ev['fasta'][0]
        match=active_by_gtf.get(an['path'])
        sci=match['scientific'] if match else scientific_by_folder.get(ev['folder'],ev['folder'].split('/')[-1].replace('_',' '))
        state=match['status'] if match else 'Retained alternative; not active'
        if ev['folder'].endswith('mouse_legacy'):
            sci='Mus musculus';state='Legacy R-browser reference (GRCm39)'
        elif 'hetGla2_' in an['path']:
            sci='Heterocephalus glaber';state='Alternative RefSeq liftOver; not active'
        elif match is None and '/liver_atac_bundle/' in ev['folder']:
            state='Retained bundle copy / alternative; original index remains active'
        digest=hashlib.sha256((ROOT/an['path']).read_bytes()).hexdigest()
        files.append({'scientific':sci,'status':state,'fasta':fa['path'],'gtf':an['path'],'index':match['index'] if match else '', 'sequence_count':fa['sequences'],'total_bp':fa['total_bp'],'header':fa['header'],'annotation_header':'\n'.join(an['comments']) if an['comments'] else 'No leading GTF comments retained','feature_source':an['features'][0].split('\t')[1], 'example_transcript':re.search(r'transcript_id "([^"]+)"',an['features'][0]).group(1) if re.search(r'transcript_id "([^"]+)"',an['features'][0]) else 'First feature has no transcript ID','gtf_sha256':digest})
files.sort(key=lambda r:(r['scientific'],r['status']!='Active',r['gtf']))
assert len(rows)==35 and sum(r['status']=='Active' for r in rows)==31
assert len(files)==44
assert len({r['scientific'] for r in rows})==35
for r in rows:
    assert (ROOT/r['fasta']).exists() and (ROOT/r['gtf']).exists()
    if r['index']: assert (ROOT/r['index']).exists()
(OUT/'inventory.json').write_text(json.dumps({'as_of':'2026-10-16','root':str(ROOT),'species':rows,'files':files},indent=2))
print(json.dumps({'species':len(rows),'active':sum(r['status']=='Active' for r in rows),'file_pairs':len(files),'unresolved_assemblies':[r['scientific'] for r in rows if not r['accession']],'missing_annotation_urls':[r['scientific'] for r in rows if not r['annotation_url']]},indent=2))
