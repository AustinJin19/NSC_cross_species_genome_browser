#!/usr/bin/env Rscript
# Post-analysis selection: preserve all original multiple-testing corrections.
base <- 'outputs/rgreat_hg38'
out <- file.path(base, 'focused_terms')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
d <- read.delim(file.path(base, 'GREAT_GO_all_terms.tsv'), check.names = FALSE)
keep <- with(d, gene_set_size >= 10 & gene_set_size <= 300 &
                fold_enrichment >= 2 & observed_region_hits >= 5 & p_adjust_all_GO < 0.05)
focused <- d[keep, ]
focused <- focused[order(focused$p_adjust_all_GO, focused$p_value, -focused$fold_enrichment), ]
strict <- focused[focused$fold_enrichment >= 3, ]
write_tsv <- function(x, name) write.table(x, file.path(out, name), sep = '\t',
                                         quote = FALSE, row.names = FALSE, na = 'NA')
write_tsv(focused, 'GREAT_specific_terms_size10-300_FE2.tsv')
write_tsv(strict, 'GREAT_specific_terms_size10-300_FE3.tsv')
bp <- focused[focused$ontology == 'BP', ]
write_tsv(bp, 'GREAT_specific_GO_BP.tsv')

# Display-only redundancy reduction. Greedily select by original global FDR;
# group terms with gene-set Jaccard >= 0.7 under an earlier representative.
suppressPackageStartupMessages(library(rGREAT))
res <- readRDS(file.path(base, 'GREAT_GO_BP.rds'))
sets <- lapply(res@gene_sets, unique)
representatives <- character()
mapping <- data.frame(term_id = character(), representative_id = character(), jaccard = numeric())
for (id in bp$id) {
  scores <- vapply(representatives, function(r) {
    length(intersect(sets[[id]], sets[[r]])) / length(union(sets[[id]], sets[[r]]))
  }, 0)
  if (length(scores) && max(scores) >= 0.7) {
    rep_id <- representatives[which.max(scores)]
    score <- max(scores)
  } else {
    representatives <- c(representatives, id)
    rep_id <- id
    score <- 1
  }
  mapping <- rbind(mapping, data.frame(term_id = id, representative_id = rep_id, jaccard = score))
}
reps <- bp[match(representatives, bp$id), ]
write_tsv(reps, 'GREAT_specific_GO_BP_representatives.tsv')
write_tsv(mapping, 'GO_BP_redundancy_mapping.tsv')
plot_data <- head(reps, 15)
write_tsv(plot_data, 'GREAT_specific_GO_BP_plot_data.tsv')
png(file.path(out, 'GREAT_specific_GO_BP.png'), width = 2400, height = 1500, res = 180)
par(mar = c(5, 28, 5, 3), mgp = c(2.5, 0.7, 0))
z <- plot_data[nrow(plot_data):1, ]
labels <- vapply(z$description, function(s) paste(strwrap(s, width = 49), collapse = '\n'), '')
b <- barplot(z$fold_enrichment, names.arg = labels, horiz = TRUE, las = 1,
             cex.names = 0.78, col = '#287F8E', border = NA,
             xlim = c(0, max(z$fold_enrichment) * 1.22),
             main = 'Specific enriched biological processes',
             xlab = 'Region-based fold enrichment vs whole-genome background')
text(z$fold_enrichment, b, labels = paste0('  ', z$gene_set_size, ' genes'), pos = 4, cex = 0.73)
mtext('10-300 genes | fold enrichment >=2 | >=5 region hits | global GO FDR <0.05', side = 3, line = 1, cex = 0.8)
mtext('Top 15 representatives ranked by original FDR; gene-set Jaccard >=0.70 grouped for display',
      side = 1, line = 3.6, cex = 0.65)
dev.off()

source_index <- match(focused$id, d$id)
stopifnot(!anyNA(source_index), identical(focused$p_adjust_all_GO, d$p_adjust_all_GO[source_index]),
          all(strict$fold_enrichment >= 3), nrow(mapping) == nrow(bp),
          all(mapping$representative_id %in% reps$id))
saved <- read.delim(file.path(out, 'GREAT_specific_terms_size10-300_FE2.tsv'))
stopifnot(nrow(saved) == sum(keep), identical(saved$id, focused$id))
summary <- data.frame(selection = c('Fold >=2, all GO', 'Fold >=2, BP', 'Fold >=2, MF',
                                    'Fold >=2, CC', 'Fold >=3, all GO', 'Fold >=3, BP',
                                    'Fold >=3, MF', 'Fold >=3, CC', 'BP representatives'),
                      terms = c(nrow(focused), sum(focused$ontology == 'BP'),
                                sum(focused$ontology == 'MF'), sum(focused$ontology == 'CC'),
                                nrow(strict), sum(strict$ontology == 'BP'),
                                sum(strict$ontology == 'MF'), sum(strict$ontology == 'CC'), nrow(reps)))
write_tsv(summary, 'filter_summary.tsv')
print(summary)
print(head(reps[,c('description','gene_set_size','observed_region_hits','fold_enrichment','p_adjust_all_GO')], 15), row.names = FALSE)
