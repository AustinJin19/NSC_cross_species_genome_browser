#!/usr/bin/env Rscript
# Run: Rscript plot_peak_overlaps.R [species_peak_overlaps.tsv] [output_prefix]
# Shared segments use >=1 bp overlap and maximal contiguous coverage by the
# indicated number of distinct samples. One-sample species are unavailable.
library(ggplot2)
args <- commandArgs(trailingOnly=TRUE)
script_arg <- grep("^--file=",commandArgs(),value=TRUE)
script_dir <- if(length(script_arg)) dirname(normalizePath(sub("^--file=","",script_arg[1]))) else getwd()
input <- if(length(args)>=1) args[1] else file.path(script_dir,"species_peak_overlaps.tsv")
prefix <- if(length(args)>=2) args[2] else file.path(script_dir,"shared_peak_regions_barplot")
d <- read.delim(input,stringsAsFactors=FALSE,check.names=FALSE,quote="",comment.char="",na.strings=c("","NA"))
stopifnot(nrow(d)==32, sum(d$samples_available)==80, all(is.na(d$shared_at_least_two)==(d$samples_available<2)))
d <- d[order(is.na(d$shared_at_least_two),-d$shared_at_least_two,d$common_name,na.last=TRUE),]
species_order <- d$species
species_labels <- setNames(paste0(d$common_name,"  (n=",d$samples_available,")"),d$species)
metrics <- c("Shared in at least 2 samples","Shared in all available samples")
long <- rbind(data.frame(species=d$species,metric=metrics[1],value=d$shared_at_least_two),
              data.frame(species=d$species,metric=metrics[2],value=d$shared_all_samples))
long$species <- factor(long$species,levels=rev(species_order))
long$metric <- factor(long$metric,levels=metrics)
fmt <- function(x) format(x,big.mark=",",scientific=FALSE,trim=TRUE)
known <- long[!is.na(long$value),]
unknown <- long[is.na(long$value),]
known$label <- fmt(known$value)
xmax <- ceiling(max(known$value)*1.2/50000)*50000
plot <- ggplot(long,aes(y=species)) +
  geom_col(data=known,aes(x=value,fill=metric),width=0.68) +
  geom_text(data=known,aes(x=value,label=label),hjust=-0.12,size=2.8,colour="#25323B") +
  geom_text(data=unknown,aes(x=xmax*0.02,label="Not assessable (one sample)"),hjust=0,size=2.8,colour="#6C747B") +
  facet_grid(cols=vars(metric),scales="fixed",drop=FALSE) +
  scale_y_discrete(labels=species_labels,drop=FALSE) +
  scale_x_continuous(limits=c(0,xmax),breaks=seq(0,xmax,by=100000),labels=fmt,expand=expansion(mult=c(0.01,0.015))) +
  scale_fill_manual(values=setNames(c("#3C5488","#00A087"),metrics),guide="none") +
  labs(x="Number of shared peak regions",y=NULL,title=NULL,caption=NULL) +
  theme_minimal(base_size=11,base_family="sans") +
  theme(plot.background=element_rect(fill="white",colour=NA),
        panel.background=element_rect(fill="white",colour=NA),
        panel.border=element_rect(fill=NA,colour="black",linewidth=0.6),
        strip.text=element_text(face="bold",size=10.5,margin=margin(8,5,8,5)),
        strip.background=element_rect(fill="white",colour="black",linewidth=0.6),
        axis.text.y=element_text(size=10,colour="#303846"),axis.text.x=element_text(size=10,colour="#303846"),
        axis.title.x=element_text(size=11,margin=margin(t=8)),
        panel.grid.minor=element_blank(),panel.grid.major=element_line(colour="#E8ECF0",linewidth=0.3),
        panel.spacing.x=grid::unit(1.1,"lines"),plot.margin=margin(8,12,8,8))
built <- ggplot_build(plot)
stopifnot(nrow(built$data[[1]])==2*sum(d$samples_available>=2), nrow(built$data[[3]])==2*sum(d$samples_available<2))
ggsave(paste0(prefix,".pdf"),plot,device=grDevices::pdf,width=12,height=10.6,units="in",bg="white",useDingbats=FALSE)
ggsave(paste0(prefix,".png"),plot,device="png",width=12,height=10.6,units="in",dpi=300,bg="white")
message("Saved shared-region bar plots for ",nrow(d)," species; ",sum(d$samples_available<2)," have only one sample.")
if(interactive()) print(plot)
