#!/usr/bin/env Rscript
# Single-panel summary of contiguous peak-covered regions shared by >=2 samples.
# Input counts were calculated by count_peak_overlaps.py (>=1 shared bp).
# Exact duplicate entries within a sample never increase sample support.
library(ggplot2)
args <- commandArgs(trailingOnly=TRUE)
script_arg <- grep("^--file=",commandArgs(),value=TRUE)
script_dir <- if(length(script_arg)) dirname(normalizePath(sub("^--file=","",script_arg[1]))) else getwd()
input <- if(length(args)>=1) args[1] else file.path(script_dir,"species_peak_overlaps.tsv")
prefix <- if(length(args)>=2) args[2] else file.path(script_dir,"shared_peaks_at_least_two")
d <- read.delim(input,stringsAsFactors=FALSE,check.names=FALSE,quote="",comment.char="",na.strings=c("","NA"))
stopifnot(nrow(d)==32,sum(d$samples_available)==80,
          all(is.na(d$shared_at_least_two)==(d$samples_available<2)))
d <- d[order(is.na(d$shared_at_least_two),-d$shared_at_least_two,d$common_name,na.last=TRUE),]
labels <- setNames(paste0(d$common_name,"  (n=",d$samples_available,")"),d$species)
d$species <- factor(d$species,levels=rev(d$species))
known <- d[!is.na(d$shared_at_least_two),]
unknown <- d[is.na(d$shared_at_least_two),]
fmt <- function(x) format(x,big.mark=",",scientific=FALSE,trim=TRUE)
known$label <- fmt(known$shared_at_least_two)
xmax <- ceiling(max(known$shared_at_least_two)*1.12/25000)*25000
p <- ggplot(d,aes(y=species)) +
  geom_col(data=known,aes(x=shared_at_least_two),width=0.68,fill="#3C5488") +
  geom_text(data=known,aes(x=shared_at_least_two,label=label),hjust=-0.15,size=3.1,colour="#25323B") +
  geom_text(data=unknown,aes(x=xmax*0.02,label="Not assessable (one sample)"),hjust=0,size=3.1,colour="#6C747B") +
  scale_y_discrete(labels=labels,drop=FALSE) +
  scale_x_continuous(limits=c(0,xmax),breaks=seq(0,xmax,by=50000),labels=fmt,
                     expand=expansion(mult=c(0.005,0.01))) +
  labs(x="Peak regions shared in at least 2 samples",y=NULL,title=NULL,caption=NULL) +
  theme_minimal(base_size=11,base_family="sans") +
  theme(plot.background=element_rect(fill="white",colour=NA),
        panel.background=element_rect(fill="white",colour=NA),
        panel.border=element_rect(fill=NA,colour="black",linewidth=0.6),
        axis.text.y=element_text(size=10.5,colour="#303846"),
        axis.text.x=element_text(size=10,colour="#303846"),
        axis.title.x=element_text(size=11,margin=margin(t=9)),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_line(colour="#E8ECF0",linewidth=0.3),
        plot.margin=margin(8,12,8,8))
built <- ggplot_build(p)
stopifnot(nrow(built$data[[1]])==28,nrow(built$data[[3]])==4,
          all(built$data[[1]]$x==known$shared_at_least_two))
ggsave(paste0(prefix,".pdf"),p,device=grDevices::pdf,width=9.8,height=10.3,
       units="in",bg="white",useDingbats=FALSE)
ggsave(paste0(prefix,".png"),p,device="png",width=9.8,height=10.3,
       units="in",dpi=300,bg="white")
tab <- data.frame(Common_name=d$common_name,Scientific_name=gsub("_"," ",as.character(d$species)),
                  Available_samples=d$samples_available,Shared_peak_regions_at_least_two_samples=d$shared_at_least_two,
                  Status=ifelse(d$samples_available>=2,"Counted","Not assessable: one sample"))
write.table(tab,paste0(prefix,".tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
message("Saved single-panel PDF/PNG and table: 28 counts, 4 not assessable.")
if(interactive()) print(p)
