#!/usr/bin/env bash
# Version-pinned Multiz470 references for 29 species in the NSC browser.
# Source URLs and provider checksums verified 2026-10-27. Bash 3+; macOS/Linux.
set -euo pipefail

usage() {
  cat <<'HELP'
Usage: bash download_multiz470_references.sh [options]

  --out-dir DIR     Destination (default: ./multiz470_reference_genomes)
  --only IDS        Comma-separated Multiz IDs, e.g. hg38,mm10,HLacoCah1
  --dry-run         Print selected URLs and destinations without downloading
  --list            Print the embedded download manifest and exit
  --convert-2bit    Also convert the two custom 2bit references to FASTA.gz
                   Requires UCSC twoBitToFa on PATH before starting
  -h, --help       Show this help

Downloads 32 files (~25.13 GB) covering all 29 overlapping species by default.
Includes both available versions for mouse, fruit bat and woodchuck.
The extra rat build HLratNor7 is suppressed at NCBI and is explicitly skipped;
rn6 is included. No newer assembly is silently substituted.

Files are grouped by scientific species name, then Multiz assembly ID.
HTTPS is used for transfer; the companion ftp_paths.txt lists FTP equivalents.
Interrupted transfers retain .part files and resume on the next run. Existing
completed files are verified and skipped. Checksum failures are reported and
files are left in place for inspection; they are never silently overwritten.
30 FASTA.gz files have publisher MD5 checksums. The two Hiller 2bit files have
byte-size and magic-byte checks only (no published checksum was found).
Different NCBI/Multiz sequence naming conventions may require alias mapping.
HELP
}

manifest() {
  # ID, species directory, filename, format, bytes, MD5, HTTPS URL, status
  cat <<'MULTIZ_MANIFEST'
hg38	Homo_sapiens	hg38.fa.gz	fasta.gz	983659424	1c9dcaddfa41027f17cd8f7a82c7293b	https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz	available
hetGla2	Heterocephalus_glaber	hetGla2.fa.gz	fasta.gz	774035182	98c398f2f24b80b839ca743b218af9e8	https://hgdownload.soe.ucsc.edu/goldenPath/hetGla2/bigZips/hetGla2.fa.gz	available
nanGal1	Nannospalax_galili	GCF_000622305.1_S.galili_v1.0_genomic.fna.gz	fasta.gz	936159309	c0c548a3e7b45317aef5d2d2b50a867b	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/622/305/GCF_000622305.1_S.galili_v1.0/GCF_000622305.1_S.galili_v1.0_genomic.fna.gz	available
HLfukDam2	Fukomys_damarensis	GCA_012274545.1_DMR_v1.0_HiC_genomic.fna.gz	fasta.gz	731905891	c2f0185a43c603b377591fffc701e655	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/274/545/GCA_012274545.1_DMR_v1.0_HiC/GCA_012274545.1_DMR_v1.0_HiC_genomic.fna.gz	available
mm10	Mus_musculus	mm10.fa.gz	fasta.gz	870141360	db005b65828db31735f384e4c5787be5	https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz	available
mm39	Mus_musculus	mm39.fa.gz	fasta.gz	870543764	41ace1b7157b98b393746aef5d1287a8	https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz	available
rn6	Rattus_norvegicus	rn6.fa.gz	fasta.gz	896422133	3a7f116d0ffce00ffd19f8885b3d28c3	https://hgdownload.soe.ucsc.edu/goldenPath/rn6/bigZips/rn6.fa.gz	available
HLratNor7	Rattus_norvegicus	-	-	0	-	-	unavailable
rheMac10	Macaca_mulatta	rheMac10.fa.gz	fasta.gz	941070922	704053b75131c042f09ef4bae0b1a5ab	https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/bigZips/rheMac10.fa.gz	available
HLacoCah1	Acomys_cahirinus	GCA_004027535.1_AcoCah_v1_BIUU_genomic.fna.gz	fasta.gz	734019227	8a0ab73711ef592824cc00dd232efed8	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/027/535/GCA_004027535.1_AcoCah_v1_BIUU/GCA_004027535.1_AcoCah_v1_BIUU_genomic.fna.gz	available
HLhydHyd1	Hydrochoerus_hydrochaeris	GCA_004027455.1_HydHyd_v1_BIUU_genomic.fna.gz	fasta.gz	873091633	a3d060ef8e5539d6a5df70d125cd09bf	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/027/455/GCA_004027455.1_HydHyd_v1_BIUU/GCA_004027455.1_HydHyd_v1_BIUU_genomic.fna.gz	available
HLeonSpe1	Eonycteris_spelaea	GCA_003508835.1_Espe.v1_genomic.fna.gz	fasta.gz	626523388	7893d9f5e2a91e89628b494a8329d764	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/508/835/GCA_003508835.1_Espe.v1/GCA_003508835.1_Espe.v1_genomic.fna.gz	available
HLmyoCoy1	Myocastor_coypus	GCA_004027025.1_MyoCoy_v1_BIUU_genomic.fna.gz	fasta.gz	894555059	397c773584a850e0dc720386408aa018	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/027/025/GCA_004027025.1_MyoCoy_v1_BIUU/GCA_004027025.1_MyoCoy_v1_BIUU_genomic.fna.gz	available
HLbalMys1	Balaena_mysticetus	HLbalMys1.2bit	2bit	623609214	-	https://bds.mpi-cbg.de/hillerlab/120MammalAlignment/Human120way/data/Genome2bitFiles/HLbalMys1.2bit	available
HLartJam1	Artibeus_jamaicensis	GCA_004027435.1_ArtJam_v1_BIUU_genomic.fna.gz	fasta.gz	781823167	3232caaef3318a6684b4d69ccc0d8040	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/027/435/GCA_004027435.1_ArtJam_v1_BIUU/GCA_004027435.1_ArtJam_v1_BIUU_genomic.fna.gz	available
HLartJam2	Artibeus_jamaicensis	GCA_014825515.1_WHU_Ajam_v2_genomic.fna.gz	fasta.gz	695411006	3e3de2b6b287d58b82948f09db618d54	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/825/515/GCA_014825515.1_WHU_Ajam_v2/GCA_014825515.1_WHU_Ajam_v2_genomic.fna.gz	available
HLcasCan3	Castor_canadensis	HLcasCan3.2bit	2bit	653105646	-	https://genome.senckenberg.de/download/TOGA/MammalianDNAZooAssemblies/HLcasCan3.2bit	available
HLcriGri3	Cricetulus_griseus	GCA_003668045.1_CriGri-PICR_genomic.fna.gz	fasta.gz	750620924	d05d38abc1ec09c7b5870ae437c89df7	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/668/045/GCA_003668045.1_CriGri-PICR/GCA_003668045.1_CriGri-PICR_genomic.fna.gz	available
chiLan1	Chinchilla_lanigera	GCF_000276665.1_ChiLan1.0_genomic.fna.gz	fasta.gz	735935420	a46fb100cff8db5338411b4c8deec622	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/276/665/GCF_000276665.1_ChiLan1.0/GCF_000276665.1_ChiLan1.0_genomic.fna.gz	available
octDeg1	Octodon_degus	GCF_000260255.1_OctDeg1.0_genomic.fna.gz	fasta.gz	822465315	ee08f13da9110b851cebd37352aad134	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/260/255/GCF_000260255.1_OctDeg1.0/GCF_000260255.1_OctDeg1.0_genomic.fna.gz	available
HLellLut1	Bramus_lutescens	GCA_001685075.1_ASM168507v1_genomic.fna.gz	fasta.gz	729499326	d48d7e6c6d23e03b7f0aa20e107f1df4	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/685/075/GCA_001685075.1_ASM168507v1/GCA_001685075.1_ASM168507v1_genomic.fna.gz	available
eptFus1	Eptesicus_fuscus	GCF_000308155.1_EptFus1.0_genomic.fna.gz	fasta.gz	588077879	4c43d3ebd71dfde198d3f3d9a7072931	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/308/155/GCF_000308155.1_EptFus1.0/GCF_000308155.1_EptFus1.0_genomic.fna.gz	available
HLmerUng1	Meriones_unguiculatus	GCA_002204375.1_MunDraft-v1.0_genomic.fna.gz	fasta.gz	776659510	818a32f838b1734536ac817223fa57c7	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/204/375/GCA_002204375.1_MunDraft-v1.0/GCA_002204375.1_MunDraft-v1.0_genomic.fna.gz	available
mesAur1	Mesocricetus_auratus	GCF_000349665.1_MesAur1.0_genomic.fna.gz	fasta.gz	679322068	483d7350bbb9d7cab2f7022e59949545	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/349/665/GCF_000349665.1_MesAur1.0/GCF_000349665.1_MesAur1.0_genomic.fna.gz	available
HLmolMol2	Molossus_molossus	GCA_014108415.1_mMolMol1.p_genomic.fna.gz	fasta.gz	717441598	f7f2b0a63e3633b22c891bdc797716c0	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/108/415/GCA_014108415.1_mMolMol1.p/GCA_014108415.1_mMolMol1.p_genomic.fna.gz	available
HLereDor1	Erethizon_dorsatum	GCA_006547115.1_GSC_porc_1.0_genomic.fna.gz	fasta.gz	765037104	c9737ff987112f1032fda0954c670766	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/006/547/115/GCA_006547115.1_GSC_porc_1.0/GCA_006547115.1_GSC_porc_1.0_genomic.fna.gz	available
HLmarMon1	Marmota_monax	GCA_901343595.1_MONAX5_genomic.fna.gz	fasta.gz	807713633	3aaa14dc8b4604f7d01b93eb7b08ee0a	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/901/343/595/GCA_901343595.1_MONAX5/GCA_901343595.1_MONAX5_genomic.fna.gz	available
HLmarMon2	Marmota_monax	GCA_014533835.1_Woodchuck_1.0_genomic.fna.gz	fasta.gz	861680456	355ef68442be4a8bf5c5b07eb3cf3e6f	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/533/835/GCA_014533835.1_Woodchuck_1.0/GCA_014533835.1_Woodchuck_1.0_genomic.fna.gz	available
HLperManBai2	Peromyscus_maniculatus	GCA_003704035.1_HU_Pman_2.1_genomic.fna.gz	fasta.gz	766251034	95b19b0156317935a821ecf90047c211	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/704/035/GCA_003704035.1_HU_Pman_2.1/GCA_003704035.1_HU_Pman_2.1_genomic.fna.gz	available
cavPor3	Cavia_porcellus	cavPor3.fa.gz	fasta.gz	875388603	5ad43ad5ec4e05cadbcb20dfb564633d	https://hgdownload.soe.ucsc.edu/goldenPath/cavPor3/bigZips/cavPor3.fa.gz	available
HLsciCar1	Sciurus_carolinensis	GCA_902686445.1_mSciCar1.1_genomic.fna.gz	fasta.gz	893683921	ba5bad8144041f099e58a7998bc12dec	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/902/686/445/GCA_902686445.1_mSciCar1.1/GCA_902686445.1_mSciCar1.1_genomic.fna.gz	available
HLscaAqu1	Scalopus_aquaticus	GCA_004024925.1_ScaAqu_v1_BIUU_genomic.fna.gz	fasta.gz	652033225	7b084496fb86149208e86f1dab0e689a	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/024/925/GCA_004024925.1_ScaAqu_v1_BIUU/GCA_004024925.1_ScaAqu_v1_BIUU_genomic.fna.gz	available
HLondZib1	Ondatra_zibethicus	GCA_004026605.1_OndZib_v1_BIUU_genomic.fna.gz	fasta.gz	819474788	d135b338cb0874094fb8705a569f1227	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/026/605/GCA_004026605.1_OndZib_v1_BIUU/GCA_004026605.1_OndZib_v1_BIUU_genomic.fna.gz	available
MULTIZ_MANIFEST
}

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
out_dir='./multiz470_reference_genomes'
only=''
dry_run=0
convert=0
while [ "$#" -gt 0 ]; do
  case "$1" in
    --out-dir) [ "$#" -ge 2 ] && [ -n "$2" ] || die '--out-dir needs a directory'; out_dir=$2; shift 2 ;;
    --only) [ "$#" -ge 2 ] && [ -n "$2" ] || die '--only needs one or more IDs'; only=$2; shift 2 ;;
    --dry-run) dry_run=1; shift ;;
    --list) manifest; exit 0 ;;
    --convert-2bit) convert=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1" ;;
  esac
done

if [ -n "$only" ]; then
  case "$only" in ,*|*,|*,,*|*[!a-zA-Z0-9_,]*) die 'Use comma-separated assembly IDs without spaces' ;; esac
  IFS=',' read -r -a requested <<< "$only"
  for id in "${requested[@]}"; do
    manifest | awk -F '\t' -v id="$id" '$1==id {found=1} END {exit !found}' || die "Unknown Multiz ID: $id (use --list)"
  done
fi

if [ "$dry_run" -eq 0 ]; then
  for tool in curl gzip awk wc od tr; do command -v "$tool" >/dev/null 2>&1 || die "Missing command: $tool"; done
  if command -v md5sum >/dev/null 2>&1; then md5_tool=md5sum
  elif command -v md5 >/dev/null 2>&1; then md5_tool=md5
  else die 'Install md5sum (Linux coreutils) or use macOS md5'; fi
  if [ "$convert" -eq 1 ]; then
    command -v twoBitToFa >/dev/null 2>&1 || die 'Missing twoBitToFa; install UCSC twoBitToFa or omit --convert-2bit'
  fi
  mkdir -p "$out_dir"
fi

verify_file() {
  local path=$1 format=$2 expected_bytes=$3 expected_md5=$4 size digest magic
  [ -f "$path" ] || return 1
  size=$(wc -c < "$path" | tr -d '[:space:]')
  [ "$size" = "$expected_bytes" ] || return 1
  if [ "$expected_md5" != '-' ]; then
    if [ "$md5_tool" = md5sum ]; then digest=$(md5sum "$path" | awk '{print $1}')
    else digest=$(md5 -q "$path"); fi
    [ "$digest" = "$expected_md5" ] || return 1
  fi
  if [ "$format" = 'fasta.gz' ]; then
    gzip -t "$path" || return 1
  elif [ "$format" = '2bit' ]; then
    magic=$(od -An -N4 -tx1 "$path" | tr -d '[:space:]')
    case "$magic" in 1a412743|4327411a) ;; *) return 1 ;; esac
  else return 1; fi
}

get_reference() {
  local target=$1 url=$2 format=$3 size=$4 digest=$5 partial
  partial="${target}.part"
  if [ -e "$target" ]; then
    verify_file "$target" "$format" "$size" "$digest" || {
      printf 'Existing file fails verification; inspect or rename it: %s\n' "$target" >&2; return 1;
    }
    printf 'Verified existing: %s\n' "$target"
    return 0
  fi
  # A completed .part may remain if interrupted just before the final rename.
  if ! verify_file "$partial" "$format" "$size" "$digest"; then
    if ! curl --fail --location --retry 5 --retry-delay 5 --connect-timeout 30 \
      --speed-limit 1024 --speed-time 120 --continue-at - --output "$partial" "$url"; then
      printf 'Transfer failed; partial file retained: %s\n' "$partial" >&2; return 1
    fi
  fi
  verify_file "$partial" "$format" "$size" "$digest" || {
    printf 'Downloaded file fails verification: %s\n' "$partial" >&2; return 1;
  }
  mv "$partial" "$target" || return 1
  printf 'Downloaded and verified: %s\n' "$target"
}

convert_reference() {
  local source=$1 fasta=$2
  if [ -e "$fasta" ]; then
    gzip -t "$fasta" || return 1
    printf 'Existing converted FASTA: %s\n' "$fasta"
    return 0
  fi
  # pipefail prevents a failed conversion from being accepted as an empty gzip.
  twoBitToFa "$source" stdout | gzip -c > "${fasta}.part" || return 1
  gzip -t "${fasta}.part" || return 1
  mv "${fasta}.part" "$fasta" || return 1
}

available=0
unavailable=0
failed=0
while IFS=$'\t' read -r id species filename format size digest url status; do
  if [ -n "$only" ]; then
    case ",$only," in *",$id,"*) ;; *) continue ;; esac
  fi
  if [ "$status" != 'available' ]; then
    printf 'UNAVAILABLE: %s (suppressed NCBI release; no FASTA; rat rn6 is available)\n' "$id" >&2
    unavailable=$((unavailable + 1))
    continue
  fi
  available=$((available + 1))
  dest="$out_dir/$species/$id"
  if [ "$dry_run" -eq 1 ]; then
    printf '%s\t%s\t%s/%s\n' "$id" "$url" "$dest" "$filename"
    continue
  fi
  mkdir -p "$dest"
  if ! get_reference "$dest/$filename" "$url" "$format" "$size" "$digest"; then
    failed=$((failed + 1)); continue
  fi
  if [ "$convert" -eq 1 ] && [ "$format" = '2bit' ]; then
    if ! convert_reference "$dest/$filename" "$dest/$id.fa.gz"; then
      printf 'FASTA conversion failed: %s\n' "$id" >&2
      failed=$((failed + 1))
    fi
  fi
done < <(manifest)

printf 'Selected available: %s; unavailable: %s; failed: %s\n' "$available" "$unavailable" "$failed" >&2
[ "$failed" -eq 0 ] || exit 1
[ "$available" -gt 0 ] || exit 2
