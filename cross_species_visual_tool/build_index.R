# =====================================================================
# build_index.R  —  run once per data set, before starting the app
#
#   Rscript build_index.R config.tsv index/
#
# For each species this produces:
#   <species>_genes.rds     gene table for search and TSS lookup (in memory)
#   <species>_models.tsv.gz transcript models, bgzipped + tabix indexed
#   <species>_meta.rds      file paths, sequence-name maps, chromosome sizes
#
# The GTF is read exactly once, here. The app never touches it.
#
#   BiocManager::install(c("rtracklayer","GenomicRanges","Rsamtools"))
#   install.packages("data.table")
# =====================================================================

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(Rsamtools)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[1] else "config.tsv"
index_dir   <- if (length(args) >= 2) args[2] else "index"

if (!file.exists(config_path)) stop("Config not found: ", config_path)
dir.create(index_dir, showWarnings = FALSE, recursive = TRUE)


# ---------------------------------------------------------------------
# Sequence name reconciliation
#
# The single most common reason a browser like this shows empty panels is
# that the GTF says "1", the bigWig says "chr1", and the FASTA header says
# "1 dna:chromosome chromosome:GRCh38:1:...". Resolve it once, here, and
# store the mapping rather than guessing at every query.
# ---------------------------------------------------------------------

strip_prefix <- function(x) sub("^chr", "", x)

build_seqmap <- function(from, to) {
  # Returns a named vector: names are GTF-style, values are target-style.
  direct <- intersect(from, to)
  map <- stats::setNames(direct, direct)

  missing <- setdiff(from, direct)
  if (length(missing)) {
    key <- stats::setNames(to, strip_prefix(to))
    hit <- key[strip_prefix(missing)]
    ok  <- !is.na(hit)
    if (any(ok)) map <- c(map, stats::setNames(as.character(hit[ok]), missing[ok]))
  }
  map
}

report_seqmap <- function(label, from, map) {
  unmapped <- setdiff(from, names(map))
  if (!length(unmapped)) {
    message(sprintf("    %-8s all %d sequences matched", label, length(from)))
  } else {
    message(sprintf("    %-8s %d of %d matched; unmatched: %s",
                    label, length(map), length(from),
                    paste(utils::head(unmapped, 6), collapse = ", ")))
  }
}


# ---------------------------------------------------------------------
# GTF -> gene table + transcript models
# ---------------------------------------------------------------------

parse_gtf <- function(gtf_path) {
  message("    reading GTF (this is the slow step)")

  # Exons alone are enough. Deriving transcript and gene bounds from them
  # avoids depending on feature types that vary between GENCODE, Ensembl
  # and NCBI GTFs.
  ex <- rtracklayer::import(gtf_path, format = "gtf", feature.type = "exon")
  if (!length(ex)) stop("No exon features found in ", gtf_path)

  mc <- mcols(ex)
  pick <- function(...) {
    for (nm in c(...)) if (nm %in% names(mc)) return(as.character(mc[[nm]]))
    NULL
  }

  gid <- pick("gene_id", "geneID", "gene")
  if (is.null(gid)) stop("No gene_id attribute in ", gtf_path)
  tid <- pick("transcript_id", "transcriptID")
  if (is.null(tid)) tid <- gid          # gene-level GTF; treat as one transcript
  gname <- pick("gene_name", "gene_symbol", "Name")
  if (is.null(gname)) gname <- gid
  gname[is.na(gname) | !nzchar(gname)] <- gid[is.na(gname) | !nzchar(gname)]
  biotype <- pick("gene_biotype", "gene_type")
  if (is.null(biotype)) biotype <- rep(NA_character_, length(ex))

  dt <- data.table(
    chrom  = as.character(seqnames(ex)),
    start  = start(ex),
    end    = end(ex),
    strand = as.character(strand(ex)),
    tx     = tid,
    gid    = gid,
    gname  = gname,
    btype  = biotype
  )
  dt <- dt[strand %in% c("+", "-")]
  message(sprintf("    %s exons across %s transcripts",
                  format(nrow(dt), big.mark = ","),
                  format(uniqueN(dt$tx), big.mark = ",")))

  # ---- transcripts ----
  txs <- dt[order(chrom, start), .(
    chrom  = chrom[1],
    start  = min(start),
    end    = max(end),
    strand = strand[1],
    gid    = gid[1],
    gname  = gname[1],
    btype  = btype[1],
    ex_starts = paste(start - 1L, collapse = ","),   # 0-based, matching BED
    ex_ends   = paste(end,        collapse = ",")
  ), by = tx]
  txs[, tss := fifelse(strand == "-", end, start)]

  # ---- genes ----
  genes <- txs[, .(
    chrom  = chrom[1],
    start  = min(start),
    end    = max(end),
    strand = strand[1],
    gname  = gname[1],
    btype  = btype[1],
    n_tx   = .N,
    # Canonical TSS: the most upstream transcript start.
    tss    = if (strand[1] == "-") max(end) else min(start),
    tss_all = paste(sort(unique(tss)), collapse = ",")
  ), by = gid]
  genes[, key := toupper(trimws(gname))]

  list(genes = genes, txs = txs)
}


write_models <- function(txs, out_prefix) {
  # A plain tab-delimited file, sorted and bgzipped, indexed on cols 1/2/3.
  # Tabix keeps range queries off the heap, which matters once you have
  # five species open at once.
  path_raw <- paste0(out_prefix, "_models.tsv")
  setorder(txs, chrom, start)

  out <- txs[, .(chrom, start = start - 1L, end, tx, gname, strand,
                 ex_starts, ex_ends, btype)]
  data.table::fwrite(out, path_raw, sep = "\t", col.names = FALSE,
                     quote = FALSE, na = ".")

  bgz <- Rsamtools::bgzip(path_raw, dest = paste0(path_raw, ".gz"), overwrite = TRUE)
  unlink(path_raw)
  Rsamtools::indexTabix(bgz, seq = 1, start = 2, end = 3, skip = 0, zeroBased = TRUE)
  bgz
}


# ---------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------

cfg <- data.table::fread(config_path, sep = "\t", header = TRUE,
                         colClasses = "character", data.table = FALSE)
required <- c("species", "fasta", "gtf", "bigwig")
missing_cols <- setdiff(required, names(cfg))
if (length(missing_cols)) {
  stop("config is missing column(s): ", paste(missing_cols, collapse = ", "))
}
if (!"colour" %in% names(cfg)) {
  cfg$colour <- rep(c("#2E5E4E", "#7A4B72", "#B5622B", "#31586E",
                      "#6B6438", "#8C3B3B"), length.out = nrow(cfg))
}

manifest <- list()

for (i in seq_len(nrow(cfg))) {
  sp <- cfg$species[i]
  message("\n", sp)

  for (f in c(cfg$fasta[i], cfg$gtf[i], cfg$bigwig[i])) {
    if (!file.exists(f)) stop("  missing file: ", f)
  }

  # ---- FASTA ----
  fa_path <- normalizePath(cfg$fasta[i])
  if (!file.exists(paste0(fa_path, ".fai"))) {
    message("    indexing FASTA")
    Rsamtools::indexFa(fa_path)
  }
  # scanFaIndex() returns a GRanges; names() on a plain GRanges is NULL
  # (that's per-range metadata, not seqnames), so pull sequence names via
  # seqnames() instead. Confirmed against the actually-installed Rsamtools.
  fa_seqs <- as.character(GenomicRanges::seqnames(Rsamtools::scanFaIndex(fa_path)))
  # FASTA headers often carry a description after the first space; the
  # .fai keeps only the first token, which is what we want.

  # ---- bigWig ----
  bw_path <- normalizePath(cfg$bigwig[i])
  bw_si   <- seqinfo(BigWigFile(bw_path))
  bw_seqs <- seqnames(bw_si)

  # ---- GTF ----
  parsed <- parse_gtf(cfg$gtf[i])
  gtf_seqs <- unique(parsed$genes$chrom)

  message("    reconciling sequence names")
  map_bw <- build_seqmap(gtf_seqs, bw_seqs)
  map_fa <- build_seqmap(gtf_seqs, fa_seqs)
  report_seqmap("bigWig", gtf_seqs, map_bw)
  report_seqmap("FASTA",  gtf_seqs, map_fa)

  usable <- intersect(names(map_bw), names(map_fa))
  dropped <- setdiff(gtf_seqs, usable)
  if (length(dropped)) {
    message(sprintf("    dropping %d sequence(s) absent from bigWig or FASTA",
                    length(dropped)))
  }
  genes <- parsed$genes[chrom %in% usable]
  txs   <- parsed$txs[chrom %in% usable]

  prefix <- file.path(index_dir, gsub("[^A-Za-z0-9_.-]", "_", sp))
  saveRDS(genes, paste0(prefix, "_genes.rds"), compress = "xz")
  models_path <- write_models(txs, prefix)

  meta <- list(
    species    = sp,
    colour     = cfg$colour[i],
    fasta      = fa_path,
    bigwig     = bw_path,
    models     = normalizePath(models_path),
    map_bw     = map_bw,
    map_fa     = map_fa,
    chrom_size = stats::setNames(as.integer(seqlengths(bw_si)[map_bw[usable]]), usable),
    n_genes    = nrow(genes)
  )
  saveRDS(meta, paste0(prefix, "_meta.rds"))

  message(sprintf("    %s genes indexed", format(nrow(genes), big.mark = ",")))
  manifest[[sp]] <- basename(prefix)
}

saveRDS(manifest, file.path(index_dir, "manifest.rds"))
message("\nIndex written to ", normalizePath(index_dir))
message("Start the app with: shiny::runApp('app.R')")
