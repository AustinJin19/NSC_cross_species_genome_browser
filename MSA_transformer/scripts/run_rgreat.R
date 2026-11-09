#!/usr/bin/env Rscript
# Local GREAT analysis. Input BED coordinates are converted once to 1-based GRanges.
options(stringsAsFactors = FALSE, timeout = 300)
suppressPackageStartupMessages({
  library(rGREAT)
  library(GenomicRanges)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args)) args[1] else 'regions_SINE_cCRE_both_annotated.tsv'
outdir <- if (length(args) > 1) args[2] else 'outputs/rgreat_hg38'
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(x, name) write.table(x, file.path(outdir, name), sep = '\t',
                                         quote = FALSE, row.names = FALSE, na = 'NA')
x <- read.delim(input, check.names = FALSE)
stopifnot(all(c('chrom','start','end','Length') %in% names(x)),
          all(x$start >= 0), all(x$end > x$start), all(x$Length == x$end - x$start))
key <- paste(x$chrom, x$start, x$end, sep = ':')
stopifnot(!anyDuplicated(key))
gr <- GRanges(x$chrom, IRanges(x$start + 1L, x$end), strand = '*')
names(gr) <- sprintf('region_%04d', seq_len(nrow(x)))
gr$region_row <- seq_len(nrow(x))
gr_mid <- gr
start(gr_mid) <- end(gr_mid) <- floor((start(gr) + end(gr))/2)
midkey <- paste(seqnames(gr_mid), start(gr_mid), sep = ':')
stopifnot(!anyDuplicated(midkey))
results <- list()
tables <- list()
for (ontology in c('BP','MF','CC')) {
  message('Running GO:', ontology)
  cache <- file.path(outdir, paste0('GREAT_GO_', ontology, '.rds'))
  if (length(args) > 2 && args[3] == '--resume' && file.exists(cache)) {
    res <- readRDS(cache)
    stopifnot(identical(ranges(res@gr), ranges(gr)),
              identical(as.character(seqnames(res@gr)), as.character(seqnames(gr))))
  } else {
    res <- great(gr, paste0('GO:', ontology),
               tss_source = 'TxDb.Hsapiens.UCSC.hg38.knownGene',
               mode = 'basalPlusExt', extend_from = 'TSS',
               basal_upstream = 5000, basal_downstream = 1000, extension = 1000000,
               min_gene_set_size = 5, background = NULL, exclude = 'gap', cores = 1)
  }
  print(res)
  # Retain all tested terms; do not recalculate FDR after filtering by hit count.
  tb <- getEnrichmentTable(res, min_region_hits = 0)
  tb$ontology <- ontology
  tables[[ontology]] <- tb
  results[[ontology]] <- res
  saveRDS(res, file.path(outdir, paste0('GREAT_GO_', ontology, '.rds')))
}
combined <- do.call(rbind, tables)
rownames(combined) <- NULL
combined$p_adjust_all_GO <- p.adjust(combined$p_value, method = 'BH')
combined$p_adjust_hyper_all_GO <- p.adjust(combined$p_value_hyper, method = 'BH')
combined <- combined[order(combined$p_adjust_all_GO, combined$p_value), ]
write_tsv(combined, 'GREAT_GO_all_terms.tsv')
sig <- subset(combined, p_adjust_all_GO < 0.05 & fold_enrichment > 1)
write_tsv(sig, 'GREAT_GO_significant_FDR05.tsv')
for (ontology in names(tables)) {
  write_tsv(combined[combined$ontology == ontology, ], paste0('GREAT_GO_', ontology, '.tsv'))
}

# The midpoint setting aligns association assignment with local GREAT enrichment.
res <- results[['BP']]
assoc <- getRegionGeneAssociations(res, by_middle_points = TRUE, use_symbols = FALSE)
ak <- paste(seqnames(assoc), start(assoc), sep = ':')
region_rows <- match(ak, midkey)
stopifnot(!anyNA(region_rows))
ng <- lengths(assoc$annotated_genes)
links <- data.frame(region_row = rep(region_rows, ng),
                    entrez_id = as.character(unlist(assoc$annotated_genes)),
                    signed_midpoint_to_TSS_bp = as.integer(unlist(assoc$dist_to_TSS)))
links <- unique(links)
symbol_map <- mapIds(org.Hs.eg.db, keys = unique(links$entrez_id), column = 'SYMBOL',
                     keytype = 'ENTREZID', multiVals = 'first')
links$gene_symbol <- unname(symbol_map[links$entrez_id])
links$gene_symbol[is.na(links$gene_symbol)] <- links$entrez_id[is.na(links$gene_symbol)]
links <- links[order(links$region_row, abs(links$signed_midpoint_to_TSS_bp), links$entrez_id), ]
links$chrom <- x$chrom[links$region_row]
links$start <- x$start[links$region_row]
links$end <- x$end[links$region_row]
links <- links[,c('region_row','chrom','start','end','entrez_id','gene_symbol','signed_midpoint_to_TSS_bp')]
write_tsv(links, 'GREAT_region_gene_links.tsv')

# Export the annotation domains and identify the nearest TSS in this same gene set.
domain_gr <- res@extended_tss
names(domain_gr) <- NULL
domains <- unique(as.data.frame(domain_gr)[, c('seqnames','start','end','gene_id','tss_position','tss_strand')])
colnames(domains)[1:3] <- c('chrom','domain_start_1based','domain_end_1based')
write_tsv(domains, 'GREAT_regulatory_domains.tsv')
tss <- unique(domains[,c('chrom','gene_id','tss_position','tss_strand')])
all_symbol_map <- mapIds(org.Hs.eg.db, keys = unique(as.character(tss$gene_id)),
                         column = 'SYMBOL', keytype = 'ENTREZID', multiVals = 'first')
tss$gene_symbol <- unname(all_symbol_map[as.character(tss$gene_id)])
tss$gene_symbol[is.na(tss$gene_symbol)] <- tss$gene_id[is.na(tss$gene_symbol)]
annot <- x
annot$region_row <- seq_len(nrow(x))
annot$great_gene_count <- tabulate(links$region_row, nbins = nrow(x))
annot$great_gene_symbols <- 'NA'
annot$great_entrez_ids <- 'NA'
annot$great_signed_midpoint_to_TSS_bp <- 'NA'
for (i in unique(links$region_row)) {
  z <- links[links$region_row == i, ]
  annot$great_gene_symbols[i] <- paste(z$gene_symbol, collapse = ';')
  annot$great_entrez_ids[i] <- paste(z$entrez_id, collapse = ';')
  annot$great_signed_midpoint_to_TSS_bp[i] <- paste(z$signed_midpoint_to_TSS_bp, collapse = ';')
}
annot$nearest_TSS_gene_symbols <- NA_character_
annot$nearest_TSS_entrez_ids <- NA_character_
annot$nearest_TSS_distance_bp <- NA_integer_
for (i in seq_len(nrow(x))) {
  z <- tss[tss$chrom == x$chrom[i], ]
  # Distance from nearest included region base to 1-based TSS (0 if overlapping).
  d <- pmax(x$start[i] + 1L - z$tss_position, z$tss_position - x$end[i], 0L)
  best <- min(d)
  nearest <- z[d == best, ]
  nearest <- nearest[order(nearest$gene_id), ]
  annot$nearest_TSS_gene_symbols[i] <- paste(nearest$gene_symbol, collapse = ';')
  annot$nearest_TSS_entrez_ids[i] <- paste(nearest$gene_id, collapse = ';')
  annot$nearest_TSS_distance_bp[i] <- best
}
eligible <- countOverlaps(gr_mid, res@background, ignore.strand = TRUE) > 0
annot$great_status <- ifelse(!eligible, 'midpoint_excluded_from_background',
                            ifelse(annot$great_gene_count == 0, 'no_regulatory_domain_overlap', 'assigned'))
write_tsv(annot, 'regions_SINE_cCRE_GREAT_annotated_hg38.tsv')

# Preserve evidence supporting the enrichment results and package/database versions.
capture.output(sessionInfo(), file = file.path(outdir, 'sessionInfo.txt'))
capture.output(AnnotationDbi::metadata(org.Hs.eg.db), file = file.path(outdir, 'OrgDb_metadata.txt'))
capture.output(AnnotationDbi::metadata(GO.db::GO.db), file = file.path(outdir, 'GO_metadata.txt'))
txdb <- getExportedValue('TxDb.Hsapiens.UCSC.hg38.knownGene','TxDb.Hsapiens.UCSC.hg38.knownGene')
capture.output(AnnotationDbi::metadata(txdb), file = file.path(outdir, 'TxDb_metadata.txt'))
capture.output(res@param, file = file.path(outdir, 'GREAT_parameters.txt'))
saveRDS(res@background, file.path(outdir, 'background_hg38_excluding_gaps.rds'))
write_tsv(data.frame(file = normalizePath(input), md5 = unname(tools::md5sum(input))), 'input_checksum.tsv')

# Verification: independent regulatory-domain overlaps must match rGREAT links.
ov <- findOverlaps(gr_mid[eligible], res@extended_tss, ignore.strand = TRUE)
expected <- unique(paste(which(eligible)[queryHits(ov)],
                         as.character(res@extended_tss$gene_id[subjectHits(ov)]), sep = ':'))
actual <- unique(paste(links$region_row, links$entrez_id, sep = ':'))
stopifnot(setequal(expected, actual), nrow(annot) == nrow(x),
          identical(annot[,names(x)], x), all(annot$nearest_TSS_distance_bp >= 0),
          sum(eligible) == res@n_total,
          isTRUE(all.equal(combined$p_adjust_all_GO, p.adjust(combined$p_value, 'BH'))))
saved <- read.delim(file.path(outdir, 'regions_SINE_cCRE_GREAT_annotated_hg38.tsv'),
                    check.names = FALSE, na.strings = character())
stopifnot(nrow(saved) == nrow(x), identical(saved[,names(x)], x))

summary <- data.frame(metric = c('input_regions','regions_in_background','assigned_regions',
  'unassigned_regions','unique_assigned_genes','region_gene_links','GO_terms_tested',
  'GO_terms_FDR05_all_GO','GO_BP_FDR05_all_GO','GO_MF_FDR05_all_GO','GO_CC_FDR05_all_GO'),
  value = c(nrow(x), sum(eligible), sum(annot$great_gene_count > 0),
  sum(annot$great_gene_count == 0), length(unique(links$entrez_id)), nrow(links),
  nrow(combined), nrow(sig), sum(sig$ontology == 'BP'), sum(sig$ontology == 'MF'), sum(sig$ontology == 'CC')))
write_tsv(summary, 'analysis_summary.tsv')
print(summary)
print(head(combined[,c('ontology','id','description','observed_region_hits','fold_enrichment','p_adjust_all_GO')], 12))

png(file.path(outdir, 'GREAT_GO_top_terms.png'), width = 2100, height = 1600, res = 180)
par(mfrow = c(3,1), mar = c(4.8,22,2.8,2), mgp = c(2.2,0.7,0))
for (ont in c('BP','MF','CC')) {
  z <- head(sig[sig$ontology == ont, ], 8)
  if (nrow(z) == 0) {
    plot.new(); title(main = paste('GO', ont, ': no terms at global FDR < 0.05'))
    next
  }
  z <- z[nrow(z):1, ]
  labs <- vapply(z$description, function(s) paste(strwrap(s, width = 48), collapse = '\n'), '')
  barplot(-log10(pmax(z$p_adjust_all_GO, .Machine$double.xmin)),
          names.arg = labs, horiz = TRUE, las = 1, col = '#287F8E', border = NA,
          cex.names = 0.75, cex.lab = 0.85, main = paste('GO', ont, '| whole genome excluding gaps'),
          xlab = '-log10(BH FDR across all three GO categories)')
  abline(v = -log10(0.05), lty = 2, col = 'gray45')
}
dev.off()
message('Completed local rGREAT analysis and verification.')
