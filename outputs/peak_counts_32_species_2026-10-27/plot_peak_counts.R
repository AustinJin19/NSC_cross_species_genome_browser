#!/usr/bin/env Rscript
# Horizontal, matched-species box plots of sample-level narrowPeak counts.
# Run: Rscript plot_peak_counts.R [sample_peak_counts.tsv] [output_prefix]
# Peak number: valid narrowPeak records, including multiple summits per interval.
# Distinct peak number: unique chromosome/start/end tuples within each sample.
# Excludes human, macaque and bowhead whale. No additional QC threshold is used.
# Missing NMR1Liver is omitted, not converted to zero. Boxes use ggplot2 quartiles
# and 1.5*IQR whiskers; every available sample is shown as a dot. n=1 species
# have a collapsed box. Both panels share the y order and the numeric x scale.

if (!requireNamespace("ggplot2", quietly=TRUE)) stop("Install ggplot2 first.")
library(ggplot2)
args <- commandArgs(trailingOnly=TRUE)
script_arg <- grep("^--file=", commandArgs(), value=TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
input <- if(length(args)>=1) args[1] else file.path(script_dir,"sample_peak_counts.tsv")
prefix <- if(length(args)>=2) args[2] else file.path(script_dir,"peak_counts_boxplots")
d <- read.delim(input, stringsAsFactors=FALSE, check.names=FALSE, quote="", comment.char="", na.strings=c("","NA"))
stopifnot(all(c("species","common_name","sample","peak_count","unique_peak_regions") %in% names(d)))
d <- d[!d$species %in% c("Homo_sapiens","Macaca_mulatta","Balaena_mysticetus"),]
missing <- d[is.na(d$peak_count) | is.na(d$unique_peak_regions),]
d <- d[!is.na(d$peak_count) & !is.na(d$unique_peak_regions),]
stopifnot(!anyDuplicated(d$sample), all(d$peak_count>=d$unique_peak_regions), all(d$unique_peak_regions>=0))
med <- tapply(d$peak_count,d$species,median)
common <- setNames(d$common_name[!duplicated(d$species)],d$species[!duplicated(d$species)])
species_order <- names(med)[order(-med,common[names(med)])]
n <- table(d$species)
species_labels <- setNames(paste0(common[species_order],"  (n=",n[species_order],")"),species_order)
metric_names <- c("Peak number","Distinct peak number")
long <- rbind(
  data.frame(species=d$species,sample=d$sample,metric=metric_names[1],value=d$peak_count),
  data.frame(species=d$species,sample=d$sample,metric=metric_names[2],value=d$unique_peak_regions)
)
long$species <- factor(long$species,levels=rev(species_order))
long$metric <- factor(long$metric,levels=metric_names)
colors <- setNames(c("#3C5488","#00A087"),metric_names)
xmax <- ceiling(max(long$value)/100000)*100000

peak_plot <- ggplot(long,aes(x=value,y=species)) +
  geom_boxplot(aes(fill=metric,colour=metric),orientation="y",width=0.56,
               alpha=0.22,linewidth=0.45,outlier.shape=NA) +
  geom_point(aes(colour=metric),position=position_jitter(width=0,height=0.085,seed=42),
             size=1.85,alpha=0.9) +
  facet_grid(cols=vars(metric),scales="fixed",drop=FALSE) +
  scale_y_discrete(labels=species_labels,drop=FALSE) +
  scale_x_continuous(limits=c(0,xmax),breaks=seq(0,xmax,by=200000),
                     labels=function(x) format(x,big.mark=",",scientific=FALSE,trim=TRUE),
                     expand=expansion(mult=c(0.015,0.015))) +
  scale_fill_manual(values=colors,guide="none") +
  scale_colour_manual(values=colors,guide="none") +
  labs(x=NULL,y=NULL,title=NULL,subtitle=NULL,caption=NULL) +
  theme_minimal(base_size=11,base_family="sans") +
  theme(
    plot.background=element_rect(fill="white",colour=NA),
    panel.background=element_rect(fill="white",colour=NA),
    panel.border=element_rect(fill=NA,colour="black",linewidth=0.6),
    strip.text=element_text(face="bold",size=11,margin=margin(8,5,8,5)),
    strip.background=element_rect(fill="white",colour="black",linewidth=0.6),
    axis.text.y=element_text(size=10,colour="#303846"),
    axis.text.x=element_text(size=10,colour="#303846"),
    panel.grid.minor=element_blank(),
    panel.grid.major.y=element_line(colour="#E8ECF0",linewidth=0.3),
    panel.grid.major.x=element_line(colour="#E8ECF0",linewidth=0.3),
    panel.spacing.x=grid::unit(1.1,"lines"),
    plot.margin=margin(8,10,8,8)
  )

# Verify every sample appears once per facet and every species has two boxes.
built <- ggplot_build(peak_plot)
stopifnot(nrow(built$data[[1]])==2*length(species_order), nrow(built$data[[2]])==2*nrow(d))
stopifnot(all(table(long$metric)==nrow(d)),!anyNA(built$data[[2]]$x),!anyNA(built$data[[2]]$y))
dir.create(dirname(prefix),recursive=TRUE,showWarnings=FALSE)
height <- 0.30*length(species_order)+0.7
ggsave(paste0(prefix,".pdf"),peak_plot,device=grDevices::pdf,width=11.8,height=height,units="in",bg="white",useDingbats=FALSE)
ggsave(paste0(prefix,".png"),peak_plot,device="png",width=11.8,height=height,units="in",dpi=300,bg="white")
write.table(data.frame(order=seq_along(species_order),species=species_order,common_name=common[species_order],n_samples=as.integer(n[species_order]),median_peak_entries=as.numeric(med[species_order])),
            paste0(prefix,"_species_order.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
message("Plotted ",nrow(d)," samples across ",length(species_order)," species in both panels.")
message("Missing samples omitted: ",paste(missing$sample,collapse=", "))
message("Verified ",nrow(built$data[[1]])," boxes and ",nrow(built$data[[2]])," sample points.")
message("Saved PDF and PNG: ",prefix)
if(interactive()) print(peak_plot)
