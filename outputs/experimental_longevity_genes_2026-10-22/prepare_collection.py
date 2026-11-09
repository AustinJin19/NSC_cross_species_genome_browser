"""Curated primary-study evidence plus a separately labelled GenAge screen.

No expression/GWAS correlation is treated as a lifespan-extension experiment.
Run with the bundled Python runtime; no workbook authoring occurs here.
"""
from pathlib import Path
import collections
import csv
import hashlib
import io
import json
import zipfile

ROOT = Path(__file__).resolve().parent
DATE = "2026-10-22"
ROWS = []
SPECIES = {
    "Mouse": "Mus musculus", "Worm": "Caenorhabditis elegans",
    "Fruit fly": "Drosophila melanogaster", "Mosquito": "Anopheles stephensi",
}

def add(gene, animal, direction, intervention, setting, result, context,
        mechanism, limitation, year, title, url, other=""):
    ROWS.append({
        "Record ID": f"EXP{len(ROWS)+1:03d}",
        "Gene in tested organism": gene,
        "Animal": animal,
        "Scientific name": SPECIES[animal],
        "Beneficial perturbation direction": direction,
        "Intervention": intervention,
        "Tissue and timing": setting,
        "Reported lifespan result": result,
        "Sex / experimental context": context,
        "Mechanism category": mechanism,
        "Interpretation / limits": limitation,
        "Publication year": year,
        "Primary paper": title,
        "Primary source URL": url,
        "Additional primary source URL": other,
    })

# MOUSE: symbols use mouse capitalization; no automatic human orthology claim.
add("Il11", "Mouse", "Decrease", "Germline knockout", "Systemic; lifelong",
    "Pooled median lifespan increased 24.9% (120.9 to 151 weeks).",
    "Males and females; pooled percentage is not a sex-specific estimate.",
    "Inflammatory signalling", "Il11ra1 healthspan findings are not substituted for this Il11 lifespan experiment.", 2024,
    "Inhibition of IL-11 signalling extends mammalian healthspan and lifespan",
    "https://www.nature.com/articles/s41586-024-07701-9")
add("Il11", "Mouse", "Decrease", "Neutralizing anti-IL-11 antibody", "Adult treatment beginning at 75 weeks",
    "Median lifespan increased 22.5% in males and 25% in females.",
    "Both sexes; treatment of otherwise normally ageing mice.",
    "Inflammatory signalling", "Targeted protein inhibition; not a knockout or evidence that lower liver ATAC necessarily has the same effect.", 2024,
    "Inhibition of IL-11 signalling extends mammalian healthspan and lifespan",
    "https://www.nature.com/articles/s41586-024-07701-9")
add("Sirt6", "Mouse", "Increase", "Transgenic overexpression", "Systemic; lifelong",
    "Lifespan extended in both male and female transgenic mice.",
    "C57BL/6JOlaHsd; earlier mixed-background study reported a male-specific benefit.",
    "Chromatin / metabolic homeostasis", "Strain and sex matter; activity of this enzyme is not measured directly by ATAC accessibility.", 2021,
    "Restoration of energy homeostasis by SIRT6 extends healthy lifespan",
    "https://www.nature.com/articles/s41467-021-23545-7")
add("Has2", "Mouse", "Increase", "Overexpression of naked-mole-rat Has2 transgene", "Mouse expressing the naked-mole-rat gene",
    "Longer lifespan and improved health measures were reported.",
    "Lifespan experiment conducted in mice, not naked mole-rats.",
    "Extracellular matrix / hyaluronan", "The donor-species transgene and high-molecular-mass hyaluronan context matter; not a generic native Has2 expression result.", 2023,
    "Increased hyaluronan by naked mole-rat Has2 improves healthspan in mice",
    "https://www.nature.com/articles/s41586-023-06463-0")
add("Becn1", "Mouse", "Activate / alter function", "F121A knock-in", "Systemic; lifelong",
    "Pooled median lifespan increased from 26 to 29 months; extension in both sexes.",
    "C57BL/6; normal-ageing comparison included.",
    "Autophagy", "Activating mutation weakens BCL2 binding; this is not Becn1 deletion.", 2018,
    "Disruption of the beclin 1-BCL2 autophagy regulatory complex promotes longevity in mice",
    "https://www.nature.com/articles/s41586-018-0162-7")
add("Cat", "Mouse", "Increase", "Mitochondria-targeted human catalase transgene (MCAT)", "Mitochondrial targeting; lifelong",
    "Median and maximum lifespan increased by averages of about 5 and 5.5 months, respectively.",
    "Mitochondrial construct compared with other subcellular targeting strategies.",
    "Mitochondrial redox control", "Subcellular localization is central; cannot generalize to any increase in catalase.", 2005,
    "Extension of murine life span by overexpression of catalase targeted to mitochondria",
    "https://pubmed.ncbi.nlm.nih.gov/15879174/")
add("Tert", "Mouse", "Increase", "AAV-mediated telomerase gene delivery", "Treatment at 1 or 2 years of age",
    "Median lifespan increased 24% for the 1-year cohort and 13% for the 2-year cohort.",
    "Adult and old mice; non-progeroid lifespan experiment.",
    "Telomere maintenance", "Catalytically inactive TERT did not reproduce the benefit; vector and age of delivery matter.", 2012,
    "Telomerase gene therapy in adult and old mice delays aging and increases longevity without increasing cancer",
    "https://pubmed.ncbi.nlm.nih.gov/22585399/",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC3494070/")
add("Tert", "Mouse", "Increase", "Rosa26 insertion of an EF1a-Tert expression cassette", "Stable transgene knock-in",
    "Longer lifespan was reported in the engineered mice.",
    "Gain-of-expression knock-in; distinct from a point mutation in endogenous Tert.",
    "Telomere maintenance", "Primary abstract checked; effect size and sex-specific survival were not extracted here. Independent replication not assessed.", 2024,
    "Telomerase reverse transcriptase gene knock-in unleashes enhanced longevity and accelerated damage repair in mice",
    "https://pubmed.ncbi.nlm.nih.gov/39660787/")
add("Fgf21", "Mouse", "Increase", "Transgenic overexpression", "Chronic endocrine overexpression",
    "Lifespan extended in both sexes.",
    "Chow-fed mice; lifelong transgene study.",
    "Endocrine / nutrient sensing", "Developmental and reproductive consequences matter; does not establish that any adult FGF21 treatment extends life.", 2012,
    "The starvation hormone, fibroblast growth factor-21, extends lifespan in mice",
    "https://pubmed.ncbi.nlm.nih.gov/23066506/",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC3466591/")
add("Irs1", "Mouse", "Decrease", "Global knockout", "Systemic; lifelong",
    "Lifespan extension reported initially in females and reproduced in a later cohort including males.",
    "Sex/cohort effects; 2008 discovery and 2011 follow-up.",
    "Insulin / IGF signalling", "This is Irs1, not Irs2; deleting different insulin-signalling genes can have different consequences.", 2008,
    "Evidence for lifespan extension and delayed age-related biomarkers in insulin receptor substrate 1 null mice",
    "https://doi.org/10.1096/fj.07-9261com",
    "https://doi.org/10.1371/journal.pone.0016144")
add("Rps6kb1", "Mouse", "Decrease", "S6K1 knockout", "Systemic; lifelong",
    "Lifespan extended in females; male survival did not show the same benefit.",
    "Sex-specific result.",
    "mTOR / translation", "S6K1 is encoded by Rps6kb1; do not conflate with Rps6kb2 or all mTOR pathway genes.", 2009,
    "Ribosomal protein S6 kinase 1 signaling regulates mammalian life span",
    "https://pubmed.ncbi.nlm.nih.gov/19797661/")
add("Pappa", "Mouse", "Decrease", "Germline knockout", "Systemic; lifelong",
    "Mean and maximum lifespan increased approximately 30-40% in the report.",
    "Congenital loss; reduced growth/local IGF availability.",
    "Local IGF signalling", "Developmental effects are part of this model; adult knockdown is listed separately.", 2007,
    "Loss of pregnancy-associated plasma protein A extends lifespan in mice",
    "https://pubmed.ncbi.nlm.nih.gov/17681037/")
add("Pappa", "Mouse", "Decrease", "Inducible gene knockdown", "Adult induction",
    "Adult induction increased female lifespan.",
    "Female mice; adult intervention.",
    "Local IGF signalling", "Separate evidence from congenital Pappa-null dwarfism; no numerical effect transcribed here.", 2017,
    "Inducible knockdown of pregnancy-associated plasma protein-A gene expression in adult female mice extends life span",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC5506424/")
add("Ghr", "Mouse", "Decrease", "Growth hormone receptor knockout", "Systemic; lifelong",
    "Ghr deletion increased lifespan; receptor antagonism did not reproduce the same survival effect in the comparison.",
    "Growth-hormone-resistant dwarf mice.",
    "GH / IGF endocrine axis", "Congenital growth and endocrine changes; not equivalent to a small change in adult hepatic accessibility.", 2003,
    "Deletion, but not antagonism, of the mouse growth hormone receptor ... increased life span [title shortened]",
    "https://pubmed.ncbi.nlm.nih.gov/12933651/")
add("Kl", "Mouse", "Increase", "Klotho transgenic overexpression", "Systemic transgene; lifelong",
    "Lifespan extension reported in two independent transgenic lines.",
    "Male and female mice.",
    "Endocrine / mineral metabolism", "Reproductive trade-offs and transgene context should accompany interpretation.", 2005,
    "Suppression of aging in mice by the hormone Klotho",
    "https://pubmed.ncbi.nlm.nih.gov/16123266/")
add("Insr", "Mouse", "Decrease", "Fat-specific insulin receptor knockout (FIRKO)", "Adipose-targeted deletion",
    "Longer lifespan was reported in FIRKO mice.",
    "Tissue-targeted insulin receptor loss.",
    "Insulin signalling / adipose metabolism", "Not a whole-body Insr knockout; results should not be extrapolated to liver-specific Insr suppression.", 2003,
    "Extended longevity in mice lacking the insulin receptor in adipose tissue",
    "https://pubmed.ncbi.nlm.nih.gov/12543978/")
add("Myc", "Mouse", "Decrease", "Heterozygous deletion / haploinsufficiency", "Systemic partial reduction; lifelong",
    "Median lifespan increased 10.7% in males and 20.9% in females (15.1% pooled).",
    "Both sexes.",
    "Growth / ribosome biogenesis", "Partial loss, not a complete knockout; growth and cancer effects accompany the intervention.", 2015,
    "Reduced Expression of MYC Increases Longevity and Enhances Healthspan",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC4624921/")
add("Atg5", "Mouse", "Increase", "Ubiquitous transgenic overexpression", "Systemic; lifelong",
    "Median lifespan increased 17.2% in the report; both sexes showed extension.",
    "Transgenic mouse study.",
    "Autophagy", "A particular overexpression model; not evidence that loss of autophagy genes is beneficial.", 2013,
    "Overexpression of Atg5 in mice activates autophagy and extends lifespan",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC3753544/")
add("Mtor", "Mouse", "Decrease", "Hypomorphic allele", "Systemic; approximately 25% of wild-type mTOR protein",
    "Median lifespan increased about 22% in males and 19% in females.",
    "Both sexes; lifelong partial reduction.",
    "Nutrient sensing / mTOR", "Not complete knockout. Some health outcomes improved while susceptibility to infection increased.", 2013,
    "Increased mammalian lifespan and a segmental and tissue-specific slowing of aging following genetic reduction of mTOR expression",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC3784301/")
add("Pten", "Mouse", "Increase", "Extra genomic copies / modest overexpression", "Systemic; lifelong",
    "Longer lifespan reported, with a benefit beyond the reduction in cancer deaths.",
    "Pten transgenic mice.",
    "PI3K-AKT signalling / metabolism", "Gene dosage and tissue context matter; no effect size extracted here.", 2012,
    "Pten Positively Regulates Brown Adipose Function, Energy Expenditure, and Longevity",
    "https://doi.org/10.1016/j.cmet.2012.02.001")
add("Adcy5", "Mouse", "Decrease", "Adenylyl cyclase 5 knockout", "Systemic; lifelong",
    "Median lifespan increased approximately 30%.",
    "AC5-null mice.",
    "cAMP signalling / stress resistance", "Primary report supports the direction; a systematic independent-replication review was not performed.", 2007,
    "Type 5 adenylyl cyclase disruption increases longevity and protects against stress",
    "https://pubmed.ncbi.nlm.nih.gov/17662940/")
add("Sirt1", "Mouse", "Increase", "Brain-specific transgenic overexpression (BRASTO)", "Brain; hypothalamic DMH/LH mechanism",
    "Lifespan extension reported with this brain-specific transgenic strategy.",
    "Both sexes in the reported study.",
    "Neuroendocrine regulation", "Do not substitute whole-body overexpression: the lifespan benefit is construct/tissue dependent.", 2013,
    "Sirt1 extends life span and delays aging in mice through the regulation of Nk2 homeobox 1 in the DMH and LH",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC3794712/")
add("Igf1r", "Mouse", "Decrease", "Heterozygous knockout", "Systemic partial receptor reduction",
    "Original study: female mean lifespan increased 33%; the male increase was not statistically significant.",
    "Sex- and genetic-background-dependent.",
    "Insulin / IGF signalling", "Context-sensitive evidence; follow-up results do not justify treating the original magnitude as universal. Not a full knockout.", 2003,
    "IGF-1 receptor regulates lifespan and resistance to oxidative stress in mice",
    "https://pubmed.ncbi.nlm.nih.gov/12483226/",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC3223158/")
add("Bub1b", "Mouse", "Increase", "BubR1 transgenic overexpression", "Systemic; sustained high expression",
    "Median and maximum lifespan extended.",
    "Normal-ageing transgenic cohort; separate from rescue of BubR1-deficient progeroid mice.",
    "Chromosome segregation / genome stability", "BubR1 is the protein encoded by Bub1b. Sex-specific numerical effects were not extracted here.", 2013,
    "Increased expression of BubR1 protects against aneuploidy and cancer and extends healthy lifespan",
    "https://pubmed.ncbi.nlm.nih.gov/23242215/")
add("Trpv1", "Mouse", "Decrease", "Knockout", "Systemic; lifelong",
    "Median lifespan increased 11.9% in males and 15.9% in females.",
    "C57BL/6; ad libitum feeding.",
    "Sensory / neuroendocrine signalling", "Knockout lifespan evidence is distinct from short-term health effects of downstream inhibitors.", 2014,
    "TRPV1 Pain Receptors Regulate Longevity and Metabolism by Neuropeptide Signaling",
    "https://pubmed.ncbi.nlm.nih.gov/24855942/")

# NEMATODE: retain organism-specific symbols, not speculative mammalian mappings.
add("daf-2", "Worm", "Decrease", "Reduced-function mutation", "Insulin/IGF-like receptor signalling",
    "Adult lifespan more than doubled in the classic experiment.",
    "Active fertile adults under the tested temperature conditions.",
    "Insulin-like signalling", "Requires daf-16; this does not mean daf-16 loss extends life. Development and temperature matter.", 1993,
    "A C. elegans mutant that lives twice as long as wild type",
    "https://www.nature.com/articles/366461a0")
add("daf-16", "Worm", "Increase", "DAF-16 transgenic overexpression", "Integrated transgene",
    "Increased lifespan reported for the DAF-16-overexpressing transgenic animals.",
    "Construct-dependent; developmental/reproductive effects accompany expression changes.",
    "FOXO transcriptional regulation", "Direct overexpression evidence is separate from daf-16 being necessary for daf-2 longevity. No one-to-one human FOXO mapping assigned.", 2001,
    "daf-16 integrates developmental and environmental inputs to mediate aging in the nematode Caenorhabditis elegans",
    "https://pubmed.ncbi.nlm.nih.gov/11747825/")
add("hsf-1", "Worm", "Increase", "Transgenic overexpression", "Lifelong expression in the tested lines",
    "Lifespan increased in independent overexpression lines.",
    "Adult nematode lifespan; heat-shock survival alone is not the inclusion endpoint.",
    "Proteostasis / stress response", "Expression level and transgene construction matter.", 2003,
    "Regulation of aging and age-related disease by DAF-16 and heat-shock factor",
    "https://pubmed.ncbi.nlm.nih.gov/12750521/")
add("pat-10", "Worm", "Increase", "Transgenic overexpression", "Increased cytoskeletal support",
    "Overexpression increased lifespan and thermotolerance.",
    "Direct pat-10 intervention, not only an association with HSF-1.",
    "Cytoskeletal maintenance", "Stress resistance and lifespan were separately assessed; no numerical longevity estimate extracted here.", 2014,
    "HSF-1-mediated cytoskeletal integrity determines thermotolerance and life span",
    "https://pubmed.ncbi.nlm.nih.gov/25324391/")
add("skn-1", "Worm", "Activate / alter function", "Constitutively active SKN-1B/C S393A transgene", "Increased intestinal nuclear SKN-1",
    "Longer lifespan reported in the activating-transgene experiments.",
    "Benefit can occur independently of daf-16.",
    "Stress response / detoxification", "Activation/localization mutation; not equivalent to any change in skn-1 RNA abundance.", 2008,
    "Direct inhibition of the longevity-promoting factor SKN-1 by insulin-like signaling in C. elegans",
    "https://pubmed.ncbi.nlm.nih.gov/18358814/")
add("rsks-1", "Worm", "Decrease", "RNAi / loss of function", "Reduced S6 kinase signalling",
    "RNAi increased mean adult lifespan from 16.4 to 20.0 days in the displayed comparison.",
    "Control and intervention from the same experiment.",
    "Translation / nutrient sensing", "Not comparable numerically with a different laboratory's lifespan baseline.", 2007,
    "Inhibition of mRNA translation extends lifespan in Caenorhabditis elegans",
    "https://pubmed.ncbi.nlm.nih.gov/17266680/")
add("ifg-1", "Worm", "Decrease", "RNAi knockdown", "Reduced translation initiation",
    "RNAi increased mean adult lifespan from 16.4 to 22.3 days in the displayed comparison.",
    "Partial knockdown of eIF4G homolog.",
    "Translation", "Complete disruption of essential translation is not implied by an RNAi longevity result.", 2007,
    "Inhibition of mRNA translation extends lifespan in Caenorhabditis elegans",
    "https://pubmed.ncbi.nlm.nih.gov/17266680/")
add("isp-1", "Worm", "Decrease", "qm150 partial-function mutation", "Mitochondrial complex III perturbation",
    "Lifespan extension accompanied reduced respiratory activity.",
    "Specific mitochondrial mutant.",
    "Mitochondrial respiration / stress signalling", "Partial impairment; severe or complete respiratory-chain loss is not equivalent.", 2001,
    "Mitochondrial Electron Transport Is a Key Determinant of Life Span in Caenorhabditis elegans",
    "https://doi.org/10.1016/S1534-5807(01)00071-5")
add("hlh-30", "Worm", "Increase", "Transgenic overexpression", "Autophagy/lysosome transcriptional programme",
    "Overexpression extended lifespan.",
    "Direct overexpression result, in addition to requirement in other longevity models.",
    "Autophagy / lysosomes", "Tissue-specific constructs need their own validation; mammalian TFEB lifespan effects cannot be assumed.", 2013,
    "The TFEB orthologue HLH-30 regulates autophagy and modulates longevity in Caenorhabditis elegans",
    "https://pubmed.ncbi.nlm.nih.gov/23925298/")
add("xbp-1", "Worm", "Activate / alter function", "Neuronal expression of spliced active XBP-1s", "Neuronal rab-3 promoter",
    "Neuronal XBP-1s expression extended organismal lifespan.",
    "Cell-nonautonomous ER stress response.",
    "ER proteostasis", "Spliced active isoform and neuronal origin are essential distinctions from generic xbp-1 expression.", 2013,
    "XBP-1 Is a Cell-Nonautonomous Regulator of Stress Resistance and Longevity",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC4771415/")
add("lipl-4", "Worm", "Increase", "Constitutive overexpression", "Intestinal lysosomal lipase",
    "LIPL-4 overexpression promoted longevity.",
    "Lysosome-to-nucleus lipid signalling experiment.",
    "Lysosomal lipid metabolism", "Downstream genes required for the response are not automatically direct lifespan-extension interventions.", 2015,
    "Lysosomal signaling molecules regulate longevity in Caenorhabditis elegans",
    "https://doi.org/10.1126/science.1258857",
    "https://escholarship.org/uc/item/30n9p7k4")

# FLY
add("chico", "Fruit fly", "Decrease", "Loss-of-function mutation / heterozygosity", "Insulin receptor substrate",
    "Homozygous females: median lifespan up to 48% longer; heterozygotes also showed extension.",
    "Homozygous male mutants did not show the same benefit and could be shorter-lived.",
    "Insulin-like signalling", "Sex and gene dosage must be retained with the gene label.", 2001,
    "Extension of life-span by loss of CHICO, a Drosophila insulin receptor substrate protein",
    "https://pubmed.ncbi.nlm.nih.gov/11292874/")
add("foxo", "Fruit fly", "Increase", "Adult fat-body dFOXO expression", "Adult pericerebral fat body",
    "Lifespan extended following targeted adult expression.",
    "Tissue-, driver- and diet-dependent.",
    "FOXO / endocrine signalling", "Not a result for every tissue, and not direct proof for any particular human FOXO paralog.", 2004,
    "Drosophila dFOXO controls lifespan and regulates insulin signalling in brain and fat body",
    "https://pubmed.ncbi.nlm.nih.gov/15175753/")
add("park", "Fruit fly", "Increase", "Parkin overexpression", "Adult-ageing transgenic experiments",
    "Overexpression extended lifespan and reduced proteotoxicity.",
    "Lifespan extension in the tested flies, not just rescue of a Parkinsonian model.",
    "Mitochondrial quality control / proteostasis", "No universal mammalian PRKN longevity effect is inferred.", 2013,
    "Parkin overexpression during aging reduces proteotoxicity, alters mitochondrial dynamics, and extends lifespan",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC3666724/")
add("Thor", "Fruit fly", "Activate / alter function", "Activated d4E-BP transgene", "Translational repression under rich nutrition",
    "Activated 4E-BP increased lifespan under a rich diet.",
    "Mutation and nutritional context are relevant.",
    "Translation / mitochondrial metabolism", "Ordinary wild-type 4E-BP overexpression did not give the same effect under high nutrition.", 2009,
    "4E-BP extends lifespan upon dietary restriction by enhancing mitochondrial activity in Drosophila",
    "https://pubmed.ncbi.nlm.nih.gov/19804760/")
add("Pten", "Fruit fly", "Increase", "Muscle-specific overexpression", "Muscle",
    "Median lifespan increased from 55 to 66 days; maximum from 76 to 88 days in the reported comparison.",
    "Muscle-directed transgene experiment.",
    "PI3K-AKT / proteostasis", "Tissue-specific perturbation; no numerical effect is transferred to mice or humans.", 2010,
    "FOXO/4E-BP signaling in Drosophila muscles regulates organism-wide proteostasis during aging",
    "https://pubmed.ncbi.nlm.nih.gov/21111239/")
add("Atg8a", "Fruit fly", "Increase", "Enhanced neuronal expression", "Ageing brain; pan-neural driver",
    "Average adult lifespan increased 56% in the report.",
    "Neural autophagy experiment.",
    "Autophagy / proteostasis", "Do not generalize to arbitrary tissues or to mammalian LC3/GABARAP family members without orthology and perturbation evidence.", 2008,
    "Promoting basal levels of autophagy in the nervous system enhances longevity and oxidant resistance in adult Drosophila",
    "https://pubmed.ncbi.nlm.nih.gov/18059160/")

add("PTEN", "Mosquito", "Increase", "Transgenic overexpression", "Midgut epithelium",
    "Midgut PTEN overexpression increased mosquito lifespan.",
    "Anopheles stephensi; tissue-specific IIS inhibition.",
    "PI3K-AKT / insulin-like signalling", "Longevity is separate from resistance to malaria infection. Sex-specific effect size not extracted here.", 2013,
    "Overexpression of phosphatase and tensin homolog improves fitness and decreases Plasmodium falciparum development in Anopheles stephensi",
    "https://pubmed.ncbi.nlm.nih.gov/23774695/",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC3788859/")

def write_tsv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

with zipfile.ZipFile(ROOT / "sources/genage_models.zip") as archive:
    content = archive.read("genage_models.csv")
    (ROOT / "sources/genage_models.csv").write_bytes(content)
    (ROOT / "sources/release.html").write_bytes(archive.read("release.html"))
    genage = list(csv.DictReader(io.StringIO(content.decode("utf-8-sig"))))
fungi = {"Saccharomyces cerevisiae", "Schizosaccharomyces pombe", "Podospora anserina"}
screen = []
for r in genage:
    if r["organism"] in fungi or "Increase" not in r["lifespan effect"]:
        continue
    screen.append({
        "GenAge ID": r["GenAge ID"],
        "Gene symbol": r["symbol"],
        "Scientific name": r["organism"],
        "Gene name": r["name"],
        "Entrez Gene ID": r["entrez gene id"],
        "Lifespan effect (GenAge)": r["lifespan effect"],
        "Longevity influence (GenAge)": r["longevity influence"],
        "Evidence status": "Database screening record; experiment and replication not individually reviewed",
        "Source": "GenAge Build 21 (2023-08-28)",
        "Source URL": "https://genomics.senescence.info/genes/models.html",
        "NCBI Gene URL": "https://www.ncbi.nlm.nih.gov/gene/" + r["entrez gene id"] if r["entrez gene id"] else "",
    })
screen.sort(key=lambda r: (r["Scientific name"], r["Gene symbol"].lower(), r["GenAge ID"]))
unique = sorted({(r["Scientific name"], r["Gene in tested organism"]) for r in ROWS})
assert len(screen) == 852
assert len({r["Record ID"] for r in ROWS}) == len(ROWS)
assert all(r["Primary source URL"].startswith("https://") for r in ROWS)
stats = {
    "curation_date": DATE,
    "primary_intervention_records": len(ROWS),
    "unique_gene_species_pairs": len(unique),
    "unique_genes_by_species": dict(collections.Counter(a for a, _ in unique)),
    "interventions_by_species": dict(collections.Counter(r["Scientific name"] for r in ROWS)),
    "direction_counts": dict(collections.Counter(r["Beneficial perturbation direction"] for r in ROWS)),
    "screening_records": len(screen),
    "screening_by_species": dict(collections.Counter(r["Scientific name"] for r in screen)),
    "screening_mixed_effects": sum(r["Lifespan effect (GenAge)"] == "Increase and Decrease" for r in screen),
    "genage_download_url": "https://genomics.senescence.info/genes/models_genes.zip",
    "genage_archive_sha256": hashlib.sha256((ROOT / "sources/genage_models.zip").read_bytes()).hexdigest(),
}
write_tsv(ROOT / "primary_lifespan_experiments.tsv", ROWS, list(ROWS[0]))
write_tsv(ROOT / "genage_animal_screening.tsv", screen, list(screen[0]))
gene_rows = []
for organism, gene in unique:
    experiments = [r for r in ROWS if r["Scientific name"] == organism and r["Gene in tested organism"] == gene]
    gene_rows.append({"Scientific name": organism, "Gene": gene,
                      "Beneficial direction": "; ".join(sorted({r["Beneficial perturbation direction"] for r in experiments})),
                      "Record IDs": "; ".join(r["Record ID"] for r in experiments)})
write_tsv(ROOT / "curated_gene_list.tsv", gene_rows, list(gene_rows[0]))
mouse_symbols = sorted(gene.upper() for species, gene in unique if species == "Mus musculus")
(ROOT / "mouse_gene_symbols_uppercase.txt").write_text("\n".join(mouse_symbols) + "\n")
(ROOT / "collection.json").write_text(json.dumps({"experiments": ROWS, "screen": screen, "stats": stats}, indent=2, ensure_ascii=False))
(ROOT / "manifest.json").write_text(json.dumps(stats, indent=2))
(ROOT / "README.txt").write_text(f"""EXPERIMENTAL ANIMAL LIFESPAN GENE COLLECTION | {DATE}

The primary table has {len(ROWS)} selected intervention records covering {len(unique)} gene-species pairs. It records direct animal lifespan experiments in mice, C. elegans, Drosophila and mosquitoes. Increased lifespan refers to the specific intervention, controls and context in the linked paper. This is a curated starting collection, not a systematic review or a universal claim of reproducibility.

INCLUSION: gene deletion, partial loss, RNAi, knock-in/activating mutation, overexpression, gene delivery or direct targeted protein inhibition, with a reported increase in whole-animal lifespan relative to a control. Disease/progeria rescue alone, stress survival alone, cell-culture survival, human GWAS and expression correlations are not sufficient for the primary table.

INTERPRETATION: directions refer to the perturbed gene/protein, not the net direction of a downstream pathway. Loss of an inhibitor can activate its pathway. Gene dosage, tissue, timing, sex and background matter. Knock-in is not synonymous with activation; the exact allele or cassette is stated. Percentages are study-specific, not pooled estimates. Numerical effects, sample sizes and sex-stratified details were not exhaustively extracted; qualitative entries do not imply absence of an effect size in the paper. Replication and publication-bias assessment is incomplete. Positive and null findings are not balanced as in a systematic review.

GENAGE SCREEN: {len(screen)} animal organism-gene records from GenAge Build 21 (2023-08-28), downloaded {DATE}. The filter retains records whose original lifespan-effect field contains 'Increase' and excludes the three fungal species in the file. It retains {stats['screening_mixed_effects']} records labelled 'Increase and Decrease'. These are screening candidates, not individually paper-verified interventions. The bulk export does not provide enough detail to assign an exact KO/OE direction safely. Original 'longevity influence' labels are retained without conversion to intervention instructions. Newer manually curated papers need not appear in this older database build.

Source: https://genomics.senescence.info/genes/models.html
Download: https://genomics.senescence.info/genes/models_genes.zip

GENE MATCHING: keep organism identity. Worm/fly symbols were not automatically converted into human symbols; homology can be one-to-many. The uppercase mouse list is only case-normalized for symbol matching and does not establish orthology in other mammals. Has2 in the mouse experiment is a naked-mole-rat transgene. ATAC accessibility is not identical to gene expression, protein activity or an experimental perturbation; correlation directions alone cannot reproduce these causal experiments.

FILES: primary_lifespan_experiments.tsv (paper-linked intervention records), curated_gene_list.tsv (one gene per species), mouse_gene_symbols_uppercase.txt (mouse matching list), genage_animal_screening.tsv (broad database screen), sources/ (unaltered downloaded source and release information), manifest.json (counts and source checksum). The Excel workbook presents the main evidence and the separately labelled screening set.
""", encoding="utf-8")
print(json.dumps(stats, indent=2))
