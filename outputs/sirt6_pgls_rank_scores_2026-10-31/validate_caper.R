library(ape)
library(caper)
library(jsonlite)
arg <- grep('^--file=',commandArgs(),value=TRUE)
here <- dirname(normalizePath(sub('^--file=','',arg[1])))
tree <- read.tree(file.path(here,'../../cross_species_visual_tool/web_app/data/timetree_pgls.nwk'))
cases <- fromJSON(file.path(here,'validation_cases.json'),simplifyVector=FALSE)
results <- lapply(cases,function(case){
  d <- do.call(rbind,lapply(case$points,as.data.frame))
  C <- do.call(rbind,lapply(case$covariance,unlist)); dimnames(C) <- list(d$species,d$species)
  cd <- comparative.data(tree,d,names.col='species',vcv=TRUE,na.omit=FALSE)
  cd$vcv <- C[rownames(cd$data),rownames(cd$data)]
  expected <- unlist(case$expected[c('beta','se','p','lambda_ml','gls_r2','intercept')])
  values <- function(fit){
    s <- summary(fit); co <- s$coefficients
    c(co['x',1],co['x',2],co['x',4],unname(fit$param['lambda']),s$r.squared,co['(Intercept)',1])
  }
  fixed <- pgls(y~x,cd,lambda=case$expected$lambda_ml,bounds=list(lambda=c(0,1)))
  fixed_error <- max(abs(values(fixed)-expected)); stopifnot(fixed_error<2e-7)
  fitted <- tryCatch(pgls(y~x,cd,lambda='ML',bounds=list(lambda=c(0,1))),error=function(e)NULL)
  optimizer <- 'caper ML'
  if(is.null(fitted) || max(abs(values(fitted)-expected))>2e-5){
    objective <- function(lam)as.numeric(logLik(pgls(y~x,cd,lambda=lam,bounds=list(lambda=c(0,1)))))
    grid <- seq(0,1,length.out=51); scores <- vapply(grid,objective,numeric(1))
    best <- which.max(scores); best_lambda <- grid[best]; best_score <- scores[best]
    for(i in seq_along(grid)){
      if((i>1 && scores[i]<scores[i-1])||(i<length(grid) && scores[i]<scores[i+1]))next
      candidate <- optimize(objective,c(max(0,grid[i]-.02),min(1,grid[i]+.02)),maximum=TRUE,tol=1e-10)
      if(candidate$objective>best_score){best_lambda<-candidate$maximum;best_score<-candidate$objective}
    }
    fitted <- pgls(y~x,cd,lambda=best_lambda,bounds=list(lambda=c(0,1)))
    optimizer <- 'R optimize of caper likelihood across grid maxima'
  }
  ml_error <- max(abs(values(fitted)-expected)); stopifnot(ml_error<2e-5)
  data.frame(assay=case$assay,gene=case$gene,n=nrow(d),
             max_fixed_lambda_error=fixed_error,max_ML_error=ml_error,
             lambda_python=case$expected$lambda_ml,lambda_R=unname(fitted$param['lambda']),
             p_python=case$expected$p,p_R=summary(fitted)$coefficients['x',4],optimizer=optimizer)
})
results <- do.call(rbind,results)
write.table(results,file.path(here,'R_caper_validation.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
cat(nrow(results),'fits validated. Max fixed-lambda error:',max(results$max_fixed_lambda_error),
    '; max independently optimized error:',max(results$max_ML_error),'\n')
