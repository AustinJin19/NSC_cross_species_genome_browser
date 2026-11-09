#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(rGREAT); library(GenomicRanges)})
all_sine <- '--all-sine' %in% commandArgs(trailingOnly=TRUE)
prefix <- if (all_sine) 'all_SINE_enhancer' else 'Alu_enhancer'
group_label <- if (all_sine) 'Enhancer regions (all SINE types)' else 'Alu + enhancer regions'
out <- paste0('outputs/rgreat_',prefix,'_by_level')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(d, name) write.table(d, file.path(out,name), sep='\t', quote=FALSE, row.names=FALSE, na='NA')
x <- read.delim('regions_SINE_cCRE_both_annotated.tsv', check.names=FALSE)
x$original_region_row <- seq_len(nrow(x))
has_label <- function(v, labels) vapply(strsplit(v, ',', fixed=TRUE), function(z) any(trimws(z) %in% labels), logical(1))
sel <- has_label(x$ccre_type, c('Distal enhancer','Proximal enhancer'))
if (!all_sine) sel <- sel & has_label(x$sine_family, 'Alu')
x <- x[sel, ]
stopifnot(nrow(x)>0, all(x$end-x$start==x$Length), !anyDuplicated(paste(x$chrom,x$start,x$end)))
write_tsv(x, paste0(prefix,'_selected_regions.tsv'))
levels <- sort(unique(x$level))
all_tables <- list()
objects <- list()
summary <- list()
resume <- '--resume' %in% commandArgs(trailingOnly=TRUE)
for (lev in levels) {
  z <- x[x$level==lev, ]
  levname <- gsub(' ', '_', lev)
  dir.create(file.path(out,levname), showWarnings=FALSE)
  write_tsv(z, file.path(levname,'input_regions.tsv'))
  gr <- GRanges(z$chrom, IRanges(z$start+1L,z$end))
  gr$original_region_row <- z$original_region_row
  names(gr) <- as.character(z$original_region_row)
  results <- list()
  for (ont in c('BP','MF','CC')) {
    message('Running ',lev,' GO:',ont)
    cache <- file.path(out,levname,paste0('GREAT_GO_',ont,'.rds'))
    if (resume && file.exists(cache)) {
      res <- readRDS(cache)
      stopifnot(identical(ranges(res@gr),ranges(gr)), identical(as.character(seqnames(res@gr)),as.character(seqnames(gr))))
    } else {
      res <- great(gr,paste0('GO:',ont),'TxDb.Hsapiens.UCSC.hg38.knownGene',
                   mode='basalPlusExt',extend_from='TSS',basal_upstream=5000,
                   basal_downstream=1000,extension=1000000,min_gene_set_size=5,
                   background=NULL,exclude='gap',cores=1)
      saveRDS(res,cache)
    }
    results[[ont]] <- res
    tab <- getEnrichmentTable(res,min_region_hits=0)
    tab$ontology <- ont
    tab$level <- lev
    all_tables[[paste(lev,ont)]] <- tab
  }
  objects[[lev]] <- results[['BP']]
  assoc <- getRegionGeneAssociations(results[['BP']],by_middle_points=TRUE,use_symbols=FALSE)
  ids <- as.character(unlist(assoc$annotated_genes))
  summary[[lev]] <- data.frame(level=lev,input_regions=nrow(z),
    regions_in_background=results[['BP']]@n_total,
    assigned_regions=length(assoc),unique_assigned_genes=length(unique(ids)))
}
d <- do.call(rbind,all_tables); rownames(d) <- NULL
d$p_adjust_all_GO_within_level <- NA_real_
for (lev in levels) {
  ix <- d$level==lev
  d$p_adjust_all_GO_within_level[ix] <- p.adjust(d$p_value[ix],'BH')
}
d$p_adjust_all_levels_GO <- p.adjust(d$p_value,'BH')
d <- d[order(d$p_adjust_all_levels_GO,d$p_value), ]
write_tsv(d,'GREAT_all_levels_all_GO_terms.tsv')
criteria <- with(d,gene_set_size>=10 & gene_set_size<=300 & fold_enrichment>=2 & observed_region_hits>=5)
focused <- d[criteria & d$p_adjust_all_levels_GO<0.05, ]
write_tsv(focused,'GREAT_all_levels_focused_terms.tsv')
write_tsv(d[criteria & d$p_adjust_all_GO_within_level<0.05, ],'GREAT_focused_terms_within_level_FDR.tsv')

# Full-input gene links are deterministic under the same domains/background.
# Check against each new run before reusing IDs and nearest-TSS metadata.
annotations <- read.delim('outputs/rgreat_hg38/regions_SINE_cCRE_GREAT_annotated_hg38.tsv',check.names=FALSE,na.strings=character())
links <- read.delim('outputs/rgreat_hg38/GREAT_region_gene_links.tsv')
annotations <- annotations[match(x$original_region_row,annotations$region_row), ]
stopifnot(identical(annotations$chrom,x$chrom),identical(annotations$start,x$start),identical(annotations$end,x$end))
write_tsv(annotations,paste0(prefix,'_regions_GREAT_genes.tsv'))
links <- links[links$region_row %in% x$original_region_row, ]
links$level <- x$level[match(links$region_row,x$original_region_row)]
plot_tables <- list(); mappings <- list()
for (lev in levels) {
  levname <- gsub(' ','_',lev)
  z <- x[x$level==lev, ]
  res <- objects[[lev]]
  assoc <- getRegionGeneAssociations(res,by_middle_points=TRUE,use_symbols=FALSE)
  mids <- floor((z$start+1L+z$end)/2)
  matchrow <- match(paste(seqnames(assoc),start(assoc)),paste(z$chrom,mids))
  expected <- unique(paste(rep(z$original_region_row[matchrow],lengths(assoc$annotated_genes)),as.character(unlist(assoc$annotated_genes))))
  lk <- links[links$level==lev, ]
  stopifnot(setequal(expected,unique(paste(lk$region_row,lk$entrez_id))))
  write_tsv(lk,file.path(levname,'region_gene_links.tsv'))
  write_tsv(d[d$level==lev, ],file.path(levname,'GREAT_all_GO_terms.tsv'))
  f <- focused[focused$level==lev, ]
  write_tsv(f,file.path(levname,'GREAT_focused_terms.tsv'))
  write_tsv(f[f$fold_enrichment>=3, ],file.path(levname,'GREAT_focused_terms_FE3.tsv'))
  bp <- f[f$ontology=='BP', ]
  sets <- lapply(res@gene_sets,unique)
  reps <- character()
  map <- data.frame(level=character(),term_id=character(),representative_id=character(),jaccard=numeric())
  for (id in bp$id) {
    scores <- vapply(reps,function(r) length(intersect(sets[[id]],sets[[r]]))/length(union(sets[[id]],sets[[r]])),0)
    if (length(scores) && max(scores)>=0.7) {rid <- reps[which.max(scores)]; score <- max(scores)}
    else {reps <- c(reps,id); rid <- id; score <- 1}
    map <- rbind(map,data.frame(level=lev,term_id=id,representative_id=rid,jaccard=score))
  }
  representatives <- bp[match(reps,bp$id), ]
  write_tsv(representatives,file.path(levname,'GREAT_BP_representatives.tsv'))
  mappings[[lev]] <- map
  plot_tables[[lev]] <- head(representatives,7)
  summary[[lev]]$focused_BP <- sum(f$ontology=='BP')
  summary[[lev]]$focused_MF <- sum(f$ontology=='MF')
  summary[[lev]]$focused_CC <- sum(f$ontology=='CC')
  summary[[lev]]$focused_all_GO <- nrow(f)
  summary[[lev]]$focused_FE3_all_GO <- sum(f$fold_enrichment>=3)
}
write_tsv(do.call(rbind,mappings),'BP_redundancy_mapping.tsv')
plot_data <- do.call(rbind,plot_tables)
write_tsv(plot_data,'GREAT_BP_by_level_plot_data.tsv')
summary <- do.call(rbind,summary); rownames(summary)<-NULL
write_tsv(summary,'analysis_summary_by_level.tsv')
common_max <- max(plot_data$fold_enrichment,2)*1.38
draw_panel <- function(lev) {
  z <- plot_tables[[lev]]
  n <- summary$input_regions[summary$level==lev]
  if (!nrow(z)) {plot.new(); title(main=paste0(lev,' (n=',n,'): no terms pass the filters'));return(invisible(NULL))}
  z <- z[nrow(z):1, ]
  labs <- vapply(z$description,function(s) paste(strwrap(s,width=42),collapse='\n'),'')
  col <- c('Level 0'='#287F8E','Level 1'='#5D6CB0','Level 2'='#B87533','Level 3'='#377B56')[lev]
  b <- barplot(z$fold_enrichment,names.arg=labs,horiz=TRUE,las=1,col=col,border=NA,
    cex.names=0.76,xlim=c(0,common_max),main=paste0(lev,' | ',group_label,': ',n),
    xlab='Fold enrichment vs whole genome',cex.lab=0.9)
  text(z$fold_enrichment,b,labels=paste0('  adj-P=',formatC(z$p_adjust_all_levels_GO,format='e',digits=2)),pos=4,cex=0.65)
  abline(v=2,lty=3,col='gray60')
}
png(file.path(out,paste0('GREAT_',prefix,'_BP_by_level.png')),width=3200,height=2400,res=180)
par(mfrow=c(2,2),mar=c(4.5,19,3.5,1.5),oma=c(3,0,4,0),mgp=c(2.5,0.7,0))
for (lev in levels) draw_panel(lev)
mtext(paste0(group_label,': GREAT Biological Process enrichment by level'),outer=TRUE,side=3,line=1.8,cex=1.3)
mtext('10-300 genes; fold >=2; >=5 region hits; FDR <0.05 across all GO terms and all four levels',outer=TRUE,side=3,line=0.1,cex=0.9)
mtext('Top 7 representatives per level; Jaccard >=0.70 grouped. Labels: BH-adjusted P across all levels and GO categories.',outer=TRUE,side=1,line=1,cex=0.8)
dev.off()
for (lev in levels) {
  png(file.path(out,gsub(' ','_',lev),'GREAT_BP_focused.png'),width=1800,height=1200,res=160)
  par(mar=c(5,23,5,2),mgp=c(2.5,0.7,0))
  draw_panel(lev)
  mtext('10-300 genes; fold >=2; >=5 hits; FDR <0.05 across all levels and GO categories',side=3,line=0.5,cex=0.7)
  dev.off()
}
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
stopifnot(sum(summary$input_regions)==nrow(x),all(summary$regions_in_background==summary$input_regions),
          !anyDuplicated(paste(d$level,d$ontology,d$id)),
          isTRUE(all.equal(d$p_adjust_all_levels_GO,p.adjust(d$p_value,'BH'))),
          all(focused$gene_set_size>=10 & focused$gene_set_size<=300))
print(summary)
print(plot_data[,c('level','description','gene_set_size','observed_region_hits','fold_enrichment','p_adjust_all_levels_GO')],row.names=FALSE)
message('All per-level runs, exports, and validation completed.')
