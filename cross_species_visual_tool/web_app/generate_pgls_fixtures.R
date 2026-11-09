# Regenerate independent numerical references. Requires ape, caper, jsonlite.
library(ape)
library(caper)
library(jsonlite)
arg <- grep("^--file=",commandArgs(),value=TRUE)
here <- dirname(normalizePath(sub("^--file=","",arg[1])))
tree <- read.tree(file.path(here,"data/timetree_pgls.nwk"))
payload <- fromJSON(file.path(here,"data/phylogenetic_covariance.json"))
names <- gsub(" ","_",payload$species)
C <- vcv.phylo(tree)[names,names]
stopifnot(max(abs(C-payload$covariance))<1e-7)
C <- C/mean(diag(C))
set.seed(47)
n <- length(names)
x <- rnorm(n,1,0.4)
logLength <- rnorm(n,1.5,0.5)
logMass <- rnorm(n,2.5,0.9)+0.3*x
error <- as.numeric(t(chol(0.65*C+0.35*diag(n)))%*%rnorm(n))
y <- 2+0.8*x+0.6*logLength-0.1*logMass+0.7*error
dat <- data.frame(species=names,x=x,y=y,logLength=logLength,logMass=logMass)
configs <- list(
 list(name="BM_unadjusted",lambda=1,adjustLength=FALSE,adjustMass=FALSE),
 list(name="lambda_zero_length",lambda=0,adjustLength=TRUE,adjustMass=FALSE),
 list(name="fixed_mass",lambda=0.55,adjustLength=FALSE,adjustMass=TRUE),
 list(name="BM_both",lambda=1,adjustLength=TRUE,adjustMass=TRUE),
 list(name="ML_unadjusted",lambda="ML",adjustLength=FALSE,adjustMass=FALSE),
 list(name="ML_length",lambda="ML",adjustLength=TRUE,adjustMass=FALSE),
 list(name="ML_mass",lambda="ML",adjustLength=FALSE,adjustMass=TRUE),
 list(name="ML_both",lambda="ML",adjustLength=TRUE,adjustMass=TRUE))
cases <- lapply(seq_along(configs),function(i){
 cfg <- configs[[i]]
 selected <- if(i%%2==0) seq_len(n) else which(!names %in% c("Homo_sapiens","Macaca_mulatta","Balaena_mysticetus","Cricetomys_gambianus"))
 d <- dat[rev(selected),]
 form <- reformulate(c("x",if(cfg$adjustLength) "logLength",if(cfg$adjustMass) "logMass"),response="y")
 cd <- comparative.data(tree,d,names.col="species",vcv=TRUE,na.omit=FALSE)
 f <- pgls(form,cd,lambda=cfg$lambda,bounds=list(lambda=c(0,1)))
 s <- summary(f); co <- s$coefficients
 expected <- list(slope=unname(co["x",1]),slopeSE=unname(co["x",2]),t=unname(co["x",3]),
                  p=unname(co["x",4]),lambda=unname(f$param["lambda"]),df=unname(s$df[2]),
                  modelIntercept=unname(co["(Intercept)",1]),rSquared=unname(s$r.squared))
 expected$intercept <- expected$modelIntercept
 if(cfg$adjustLength){expected$lengthSlope<-unname(co["logLength",1]);expected$intercept<-expected$intercept+expected$lengthSlope*mean(d$logLength)}
 if(cfg$adjustMass){expected$massSlope<-unname(co["logMass",1]);expected$intercept<-expected$intercept+expected$massSlope*mean(d$logMass)}
 points <- data.frame(species=gsub("_"," ",d$species),x=d$x,y=d$y,lengthKb=10^d$logLength,logMass=d$logMass)
 list(name=cfg$name,options=cfg[c("lambda","adjustLength","adjustMass")],points=points,expected=expected)
})
write_json(list(reference=paste0("caper ",packageVersion("caper"),"; ape ",packageVersion("ape")),
                covariance_max_error=max(abs(vcv.phylo(tree)[names,names]-payload$covariance)),cases=cases),
           file.path(here,"test_pgls_fixtures.json"),pretty=TRUE,auto_unbox=TRUE,digits=16)
cat("Saved",length(cases),"independent caper PGLS fixtures.\n")
