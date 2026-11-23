options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(GSVA);library(jsonlite);library(digest)})
cfg<-fromJSON('5-ClassifierValidation/config/external_programs_20261103.json');run<-file.path('5-ClassifierValidation/runs',cfg$run_id);par<-'5-ClassifierValidation/runs/bulk-rnaseq__20261030_step5_consensus_v1'
inputs<-c('5-ClassifierValidation/config/external_programs_20261103.json','5-ClassifierValidation/scripts/run_external_programs_20261103.R',paste0(par,'/tables/external_sample_predictions.csv'),paste0(par,'/tables/final_model_features.csv'),paste0(par,'/models/candidate_models.rds'),'2-TranscriptomicCharacterization/runs/fig2_program_scores_20261026_v7/gene_set_membership.csv','1-Clustering/data/combined_expr_after_correction.rds','0-Data/qc/0-Data-rnaseq-microarray1.Rdata','0-Data/qc/0-Data-rnaseq-microarray2.Rdata')
fwrite(data.table(Path=inputs,SHA256=vapply(inputs,function(p)digest(file=p,algo='sha256'),character(1))),file.path(run,'provenance/inputs.csv'))
modelhash<-digest(file=inputs[5],algo='sha256');pred<-fread(inputs[3])[Model=='Consensus_SVM'&Track=='Microarray'];stopifnot(nrow(pred)==506,all(pred$TotalFeatures==114),!anyDuplicated(pred[,.(Cohort,Sample)]))
feat<-fread(inputs[4])[Model=='Consensus_SVM',Gene];stopifnot(length(feat)==114)
members<-fread(inputs[6])[Collection=='Hallmark'&Signature%in%cfg$programs&VariableMeasured==TRUE,.(Program=Signature,Gene)];sets<-split(members$Gene,members$Program);stopifnot(length(sets)==9)
fwrite(members,file.path(run,'tables/frozen_membership.csv'));fwrite(data.table(Program=names(sets),FrozenN=lengths(sets),Overlap114=vapply(sets,function(g)sum(g%in%feat),integer(1)),ExcludedN=vapply(sets,function(g)sum(!g%in%feat),integer(1))),file.path(run,'tables/classifier_overlap.csv'))
e<-new.env();load(inputs[7],e);universe<-rownames(e$corrected_expr_limma);stopifnot(length(universe)==10531);fwrite(data.table(Gene=universe),file.path(run,'tables/frozen_background.csv'));rm(e)
raw<-list();for(p in inputs[8:9]){e<-new.env();load(p,e);for(n in intersect(ls(e),unique(pred$Cohort)))raw[[n]]<-get(n,e)}
cohorts<-sort(unique(pred$Cohort));stopifnot(length(cohorts)==6,all(cohorts%in%names(raw)));fwrite(pred,file.path(run,'tables/projection_labels.csv'))
scores<-list();coverage<-list();effects<-list();used<-list();audit<-list();drawsall<-list();set.seed(cfg$seed)
for(co in cohorts){cat('Scoring',co,'\n');flush.console();z<-pred[Cohort==co];setorder(z,Sample);x<-as.matrix(raw[[co]]);storage.mode(x)<-'double';stopifnot(all(z$Sample%in%colnames(x)));x<-x[,z$Sample,drop=FALSE]
 if(anyDuplicated(rownames(x))){nn<-table(rownames(x));x<-rowsum(x,rownames(x),reorder=FALSE)/as.numeric(nn[rownames(rowsum(x,rownames(x),reorder=FALSE))])}
 stopifnot(all(is.finite(x)));bg<-intersect(universe,rownames(x));x<-x[bg,,drop=FALSE];variable<-apply(x,1,sd)>0;x<-x[variable,,drop=FALSE]
 strata<-split(seq_len(nrow(z)),z$Prediction);draws<-replicate(cfg$bootstrap_reps,unlist(lapply(strata,function(ii)ii[sample.int(length(ii),length(ii),TRUE)])),simplify=FALSE);drawsall[[co]]<-draws
 for(sc in cfg$scenarios){xx<-if(sc=='Full')x else x[!rownames(x)%in%feat,,drop=FALSE];ss<-lapply(sets,function(g)if(sc=='Full')g else setdiff(g,feat));cv<-rbindlist(lapply(names(ss),function(p){g<-intersect(ss[[p]],rownames(xx));data.table(Cohort=co,Scenario=sc,Program=p,OriginalN=length(sets[[p]]),ScenarioN=length(ss[[p]]),Overlap114=sum(sets[[p]]%in%feat),MeasuredN=length(g),Coverage=length(g)/length(ss[[p]]),OriginalCoverage=length(g)/length(sets[[p]]),Eligible=length(g)>=cfg$min_genes&&length(g)/length(ss[[p]])>=cfg$min_coverage,BackgroundN=nrow(xx))}));coverage[[length(coverage)+1]]<-cv
 used[[length(used)+1]]<-rbindlist(lapply(names(ss),function(p)data.table(Cohort=co,Scenario=sc,Program=p,Gene=intersect(ss[[p]],rownames(xx)))))
 eligible<-cv[Eligible==TRUE,Program];stopifnot(length(eligible)>0);if(sc=='Exclude114')stopifnot(!any(rownames(xx)%in%feat),!any(unlist(ss)%in%feat))
 param<-gsvaParam(xx,ss[eligible],minSize=cfg$min_genes,maxSize=Inf,kcdf='Gaussian',tau=1,maxDiff=TRUE,absRanking=FALSE);gs<-gsva(param,verbose=FALSE,BPPARAM=BiocParallel::SerialParam());stopifnot(all(is.finite(gs)),identical(colnames(gs),z$Sample));saveRDS(gs,file.path(run,'models',paste0(co,'_',sc,'_scores.rds')))
 long<-as.data.table(as.table(gs));setnames(long,c('Program','Sample','Score'));long<-merge(long,z[,.(Sample,Prediction)],by='Sample');long[,`:=`(Cohort=co,Scenario=sc)];scores[[length(scores)+1]]<-long
 for(p in cfg$programs)for(cl in c('C1','C2','C3')){n1<-sum(z$Prediction==cl);n0<-nrow(z)-n1;cvp<-cv[Program==p];eligible<-isTRUE(cvp$Eligible)&&n1>=2&&n0>=2
  rr<-data.table(Cohort=co,Scenario=sc,Program=p,Class=cl,N=nrow(z),N1=n1,N0=n0,Eligible=eligible,Reason=if(!cvp$Eligible)'Coverage below fixed rule' else if(n1<2||n0<2)'Insufficient group n' else 'Estimated',Mean1=NA_real_,Mean0=NA_real_,G=NA_real_,Lower95=NA_real_,Upper95=NA_real_,P=NA_real_)
  if(eligible){v<-as.numeric(gs[p,]);a<-z$Prediction==cl;hedges<-function(ii){vv<-v[ii];aa<-a[ii];u<-vv[aa];b<-vv[!aa];sp<-sqrt(((length(u)-1)*var(u)+(length(b)-1)*var(b))/(length(vv)-2));(1-3/(4*length(vv)-9))*(mean(u)-mean(b))/sp};bs<-vapply(draws,hedges,numeric(1));stopifnot(all(is.finite(bs)));rr[,`:=`(Mean1=mean(v[a]),Mean0=mean(v[!a]),G=hedges(seq_along(v)),Lower95=quantile(bs,.025),Upper95=quantile(bs,.975),P=t.test(v[a],v[!a],var.equal=FALSE)$p.value)]}
  effects[[length(effects)+1]]<-rr
 }
 audit[[length(audit)+1]]<-data.table(Cohort=co,Scenario=sc,N=nrow(z),InputGenes=nrow(raw[[co]]),BackgroundN=nrow(xx),EligiblePrograms=sum(cv$Eligible),Nonfinite=0)
 }
}
fwrite(rbindlist(scores),file.path(run,'tables/scores_long.csv'));fwrite(rbindlist(coverage),file.path(run,'tables/coverage.csv'));fwrite(rbindlist(used),file.path(run,'tables/measured_membership.csv'));fwrite(rbindlist(audit),file.path(run,'tables/block_audit.csv'))
ef<-rbindlist(effects);ef[,Q:=p.adjust(P,'BH'),by=Scenario];fwrite(ef,file.path(run,'tables/program_effects.csv'));saveRDS(drawsall,file.path(run,'models/bootstrap_draws.rds'))
p1<-pred[grepl('Jusakual1',Cohort),Sample];p2<-pred[grepl('Jusakual2',Cohort),Sample];write_json(list(platform_sample_id_overlap=intersect(p1,p2),patient_linkage='Not established;platform blocks kept separate;no pooled cross-platform inference',n=506),file.path(run,'review/sample_unit.json'),pretty=TRUE,auto_unbox=TRUE)
stopifnot(nrow(ef)==324,nrow(rbindlist(coverage))==108,digest(file=inputs[5],algo='sha256')==modelhash)
write_json(list(passed=TRUE,rows_effects=nrow(ef),estimable=sum(ef$Eligible),n=506,model_unchanged=TRUE),file.path(run,'review/execution.json'),pretty=TRUE,auto_unbox=TRUE);capture.output(sessionInfo(),file=file.path(run,'provenance/sessionInfo.txt'));cat('DONE\n')
