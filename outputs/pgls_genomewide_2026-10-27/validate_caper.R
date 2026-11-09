library(ape)
library(caper)
library(jsonlite)
arg<-grep('^--file=',commandArgs(),value=TRUE)
here<-dirname(normalizePath(sub('^--file=','',arg[1])))
tree<-read.tree(file.path(here,'../../cross_species_visual_tool/web_app/data/timetree_pgls.nwk'))
cases<-fromJSON(file.path(here,'validation_cases.json'),simplifyVector=FALSE)
results<-lapply(cases,function(c){
 d<-do.call(rbind,lapply(c$points,as.data.frame));d$species<-gsub(' ','_',d$species);d$logLength<-log10(d$lengthKb)
 C<-do.call(rbind,lapply(c$covariance,unlist));dimnames(C)<-list(d$species,d$species)
 form<-reformulate(c('x',if(c$method=='body')'logLength','logMass'),response='y')
 cd<-comparative.data(tree,d,names.col='species',vcv=TRUE,na.omit=FALSE)
 # Preserve the original full-tree root, exactly as the browser does. Pruning
 # to a pure subclade can otherwise remove a shared stem from caper's matrix.
 cd$vcv<-C[rownames(cd$data),rownames(cd$data)]
 scalar_fit<-function(){
  objective<-function(lam)as.numeric(logLik(pgls(form,cd,lambda=lam,bounds=list(lambda=c(0,1)))))
  grid<-seq(0,1,length.out=21);scores<-vapply(grid,objective,numeric(1));best<-which.max(scores)
  best_lambda<-grid[best];best_score<-scores[best]
  for(i in seq_along(grid)){
   if((i>1 && scores[i]<scores[i-1])||(i<length(grid) && scores[i]<scores[i+1]))next
   candidate<-optimize(objective,c(max(0,grid[i]-.05),min(1,grid[i]+.05)),maximum=TRUE,tol=1e-9)
   if(candidate$objective>best_score){best_lambda<-candidate$maximum;best_score<-candidate$objective}
  }
  pgls(form,cd,lambda=best_lambda,bounds=list(lambda=c(0,1)))
 }
 fit<-tryCatch(pgls(form,cd,lambda='ML',bounds=list(lambda=c(0,1))),error=function(e)e)
 optimizer_note<-'caper default ML optimizer'
 initial_lambda<-if(inherits(fit,'error'))NA_real_ else unname(fit$param['lambda'])
 initial_likelihood<-if(inherits(fit,'error'))NA_real_ else as.numeric(logLik(fit))
 if(inherits(fit,'error')){
  optimizer_note<-paste('R optimize of caper likelihood; default optimizer failed:',conditionMessage(fit))
  fit<-scalar_fit()
 }
 s<-summary(fit);co<-s$coefficients
 expected<-unlist(c$expected[c('MLS_beta','MLS_SE','p_value','lambda','GLS_R2','model_intercept')])
 observed<-c(co['x',1],co['x',2],co['x',4],unname(fit$param['lambda']),s$r.squared,co['(Intercept)',1])
 err<-max(abs(expected-observed))
 if(err>2e-5){
  fit<-scalar_fit();s<-summary(fit);co<-s$coefficients
  observed<-c(co['x',1],co['x',2],co['x',4],unname(fit$param['lambda']),s$r.squared,co['(Intercept)',1])
  err<-max(abs(expected-observed))
  optimizer_note<-'R optimize of caper likelihood across local maxima; default optimizer reached a lower likelihood'
  stopifnot(as.numeric(logLik(fit))>=initial_likelihood-1e-8)
 }
 stopifnot(err<2e-5)
 # Always verify the GLS calculation independently at the browser's lambda,
 # even if caper's optimizer reaches a different local likelihood maximum.
 fixed<-pgls(form,cd,lambda=c$expected$lambda,bounds=list(lambda=c(0,1)))
 sf<-summary(fixed);cf<-sf$coefficients
 fixed_values<-c(cf['x',1],cf['x',2],cf['x',4],c$expected$lambda,sf$r.squared,cf['(Intercept)',1])
 fixed_error<-max(abs(expected-fixed_values));stopifnot(fixed_error<2e-5)
 data.frame(method=c$method,gene=c$gene,n=nrow(d),lambda_scan=c$expected$lambda,lambda_caper=unname(fit$param['lambda']),
            max_difference_ML=err,max_difference_fixed_lambda=fixed_error,MLS_p_scan=c$expected$p_value,MLS_p_caper=co['x',4],
            initial_caper_lambda=initial_lambda,likelihood_improvement=as.numeric(logLik(fit))-initial_likelihood,optimizer=optimizer_note,row.names=NULL)
})
r<-do.call(rbind,results)
write.table(r,file.path(here,'R_caper_validation.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
print(r[r$max_difference_ML>2e-5,],row.names=FALSE)
cat(nrow(r),'fits: max fixed-lambda error',max(r$max_difference_fixed_lambda),'; max ML error',max(r$max_difference_ML),'\n')
