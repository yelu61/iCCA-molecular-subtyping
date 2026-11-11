suppressPackageStartupMessages({library(data.table);library(GSVA)})
run <- '2-TranscriptomicCharacterization/runs/fig2_intermediate_20261026_v1'
oldrun <- '2-TranscriptomicCharacterization/runs/fig2_program_scores_20261026_v7'
v7 <- '8-Visualization/refinements/fig2_sfig2_transcriptomic_20261026_v7'
inputs <- c('governance/data/subtype_master_table.csv','1-Clustering/data/combined_expr_after_correction.rds',file.path(oldrun,'program_scores.rds'),'2-TranscriptomicCharacterization/publication/models/step2_wgcna_network.rds',file.path(v7,'tables/score_main_display_plan.csv'),file.path(v7,'tables/wgcna_main_effects.csv'),'6-ClinicalSurvival/publication/tables/step6_survival_population_by_study_and_state.csv')
protected <- unique(c(inputs,list.files(file.path(v7,'figures'),full.names=TRUE)))
before <- tools::md5sum(protected)
meta <- fread(inputs[1]); meta[,Sample:=as.character(Sample)]
states <- c('C1','C2','C3','C1_C2','C1_C3','C2_C3')
meta[,State:=factor(FinalSubtype,levels=states)];setorder(meta,State,Study,Sample)
stopifnot(nrow(meta)==586L,!anyDuplicated(meta$Sample),identical(as.integer(table(meta$State)),c(112L,223L,162L,41L,39L,9L)))
e<-new.env();load(inputs[2],e);expr<-e$corrected_expr_limma[,meta$Sample,drop=FALSE]
stopifnot(nrow(expr)==10531L,all(is.finite(expr)))
old<-readRDS(inputs[3]);sets<-as.data.table(old$membership)
eligible<-old$coverage[Eligible==TRUE,SetID]
sets<-sets[SetID%in%eligible & VariableMeasured==TRUE]
gs<-split(sets$Gene,sets$SetID)
stopifnot(length(gs)==156L,all(unlist(gs)%in%rownames(expr)))
param<-GSVA::gsvaParam(expr,gs,minSize=5,maxSize=Inf,kcdf='Gaussian',tau=1,maxDiff=TRUE,absRanking=FALSE)
scores<-GSVA::gsva(param,verbose=FALSE,BPPARAM=BiocParallel::SerialParam())
z<-t(scale(t(scores)))
stopifnot(identical(colnames(scores),meta$Sample),all(is.finite(z)),nrow(scores)==156)
long<-melt(as.data.table(scores,keep.rownames='SetID'),id.vars='SetID',variable.name='Sample',value.name='GSVA');long[,Sample:=as.character(Sample)]
long<-merge(long,meta[,.(Sample,Study,FinalSubtype,SubtypeClass)],by='Sample')
long<-merge(long,as.data.table(old$coverage)[,.(SetID,Collection,Signature,RequestedGenes,MeasuredGenes,Coverage)],by='SetID')
fwrite(long,file.path(run,'gsva_all_states_long.csv'));fwrite(sets,file.path(run,'scored_membership.csv'))
# Quantify numerical drift caused by population-dependent scoring; do not relabel as replication.
ids<-colnames(old$scores)
audit<-rbindlist(lapply(rownames(old$scores),function(s)data.table(SetID=s,MainN=length(ids),Pearson=cor(scores[s,ids],old$scores[s,ids]),Spearman=cor(scores[s,ids],old$scores[s,ids],method='spearman'),MaxAbsoluteDelta=max(abs(scores[s,ids]-old$scores[s,ids])))))
fwrite(audit,file.path(run,'gsva_main_population_sensitivity.csv'))
model<-readRDS(inputs[4]);loads<-as.data.table(model$reference_eigengene_loadings)
mods<-sort(unique(loads$Module)); raw<-matrix(NA_real_,586,length(mods),dimnames=list(meta$Sample,mods))
for(mod in mods){x<-loads[Module==mod];a<-sweep(expr[x$Gene,,drop=FALSE],1,x$ReferenceCenter,'-');a<-sweep(a,1,x$ReferenceScale,'/');raw[,mod]<-as.numeric(crossprod(x$PC1Loading,a))}
fixed<-raw;cal<-list()
for(st in unique(meta$Study)){
 allids<-meta[Study==st,Sample];mainids<-meta[Study==st & SubtypeClass=='Main',Sample]
 mu<-colMeans(raw[mainids,,drop=FALSE]);sdv<-apply(raw[mainids,,drop=FALSE],2,sd);stopifnot(all(is.finite(sdv)),all(sdv>0))
 fixed[allids,]<-sweep(sweep(raw[allids,,drop=FALSE],2,mu,'-'),2,sdv,'/')
 cal[[st]]<-data.table(Study=st,Module=mods,MainReferenceN=length(mainids),Center=mu,SD=sdv)
}
ref<-model$fixed_loading_projected_study_standardized_scores
err<-max(abs(fixed[rownames(ref),colnames(ref)]-ref));stopifnot(err<1e-10,all(is.finite(fixed)))
mlong<-melt(as.data.table(t(fixed),keep.rownames='Module'),id.vars='Module',variable.name='Sample',value.name='Score');mlong[,Sample:=as.character(Sample)]
mlong<-merge(mlong,meta[,.(Sample,Study,FinalSubtype,SubtypeClass)],by='Sample')
fwrite(mlong,file.path(run,'wgcna_all_states_long.csv'));fwrite(rbindlist(cal),file.path(run,'wgcna_fixed_study_calibration.csv'))
fwrite(loads,file.path(run,'wgcna_frozen_loadings.csv'))
fwrite(mlong[,.(N=.N,Median=median(Score),Q25=quantile(Score,.25),Q75=quantile(Score,.75)),by=.(Module,FinalSubtype)],file.path(run,'wgcna_state_summary.csv'))
fwrite(meta[,.(N=.N),by=.(Study,FinalSubtype,SubtypeClass)],file.path(run,'state_study_counts.csv'))
cover<-fread(inputs[7]);fwrite(cover[,.(TotalN=sum(TotalN),ValidOSN=sum(ValidOSN),Events=sum(Events)),by=.(FinalSubtype,SubtypeClass)],file.path(run,'os_coverage.csv'))
fwrite(meta[,.(Sample,Study,FinalSubtype,SubtypeClass)],file.path(run,'sample_order.csv'))
after<-tools::md5sum(protected);stopifnot(identical(before,after))
fwrite(data.table(Path=protected,BeforeMD5=unname(before),AfterMD5=unname(after)),file.path(run,'protected_input_audit.csv'))
saveRDS(list(meta=meta,scores=scores,z=z,coverage=old$coverage,modules=fixed),file.path(run,'intermediate_characterization.rds'))
writeLines(c('Status: EXECUTED_REVIEW_DRAFT',paste('WGCNA max absolute difference for saved 497:',format(err,digits=16)),paste('GSVA 497-sample Spearman range:',paste(range(audit$Spearman),collapse=' to ')),'No new biological labels or inferential P values. All six states retained; intermediate does not mean temporal transition.'),file.path(run,'RESULTS.md'))
capture.output(sessionInfo(),file=file.path(run,'sessionInfo.txt'))
cat('PASS: GSVA',dim(scores),'WGCNA',dim(fixed),'main score error',err,'\n')
