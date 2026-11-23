#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(digest);library(e1071);library(pROC)})
source('5-ClassifierValidation/scripts/svm_rfe_helpers_20261031.R')
stage<-if(length(commandArgs(TRUE)))commandArgs(TRUE)[1] else 'run'
cfg<-fromJSON("5-ClassifierValidation/config/consensus_revision_20261030.json")
size_cfg<-fromJSON("5-ClassifierValidation/config/svm_rfe_20261031.json")
parent<-file.path("5-ClassifierValidation/runs",size_cfg$parent_run)
out<-file.path("5-ClassifierValidation/runs",size_cfg$run_id)
for(d in c("tables","checkpoints","provenance"))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
identity<-readLines(file.path(parent,"provenance/identity.txt"))
status<-function(state,detail){tmp<-file.path(out,paste0("status.",Sys.getpid(),".tmp"));write_json(list(state=state,detail=detail,time=format(Sys.time(),tz="UTC",usetz=TRUE)),tmp,auto_unbox=TRUE,pretty=TRUE);stopifnot(file.rename(tmp,file.path(out,"status.json")));cat(format(Sys.time(),"%H:%M:%S"),state,detail,"\n");flush.console()}
options(error=function(){status("FAILED",geterrmessage());quit(save="no",status=1)})
classes<-c("C1","C2","C3")
meta<-as.data.frame(fread(file.path(parent,"tables/population.csv"),colClasses=c(Sample="character")));meta$Subtype<-factor(meta$Subtype,levels=classes)
e<-new.env(parent=emptyenv());load("1-Clustering/data/combined_expr_before_correction.rds",e);X<-t(e$combined_expr[,meta$Sample,drop=FALSE]);storage.mode(X)<-"double";rm(e)
functions<-c("seed_for","atomic_save","correct_split","scale_split","ba","align_prob","fit_model","predict_prob")
for(expr in parse(file.path(parent,"provenance/run_step5_consensus_revision.R")))if(is.call(expr)&&identical(expr[[1]],as.name("<-"))&&as.character(expr[[2]])[1]%in%functions)eval(expr)
stopifnot(all(vapply(functions,exists,logical(1))))
source_identity<-digest(list(size_cfg,identity,digest(file='5-ClassifierValidation/scripts/svm_rfe_helpers_20261031.R',algo='sha256'),digest(file="5-ClassifierValidation/scripts/run_step5_svm_rfe_20261031.R",algo="sha256")),algo="sha256")
identity_file<-file.path(out,"provenance/identity.txt")
if(file.exists(identity_file))stopifnot(identical(readLines(identity_file),source_identity)) else writeLines(source_identity,identity_file)
for(p in c("5-ClassifierValidation/scripts/svm_rfe_helpers_20261031.R","5-ClassifierValidation/config/svm_rfe_20261031.json","5-ClassifierValidation/scripts/run_step5_svm_rfe_20261031.R"))if(!file.exists(file.path(out,"provenance",basename(p))))file.copy(p,file.path(out,"provenance",basename(p)))
cache_get<-function(indices){key<-digest(list(identity,sort(meta$Sample[indices])),algo="sha256");path<-file.path(parent,"cache",paste0("selection_",key,".rds"));stopifnot(file.exists(path));s<-readRDS(path);stopifnot(identical(s$training_samples,sort(meta$Sample[indices])));d<-as.data.table(s$support)[Methods>=2];setorder(d,-Methods,FRank,Gene);stopifnot(nrow(d)>=max(size_cfg$marker_counts));list(genes=d$Gene,path=path)}
parent_inputs<-fread(file.path(parent,'provenance/input_manifest.csv'))
checked<-parent_inputs[Role %in% c('expression','master','folds')]
stopifnot(all(vapply(checked$Path,function(p)digest(file=p,algo='sha256'),character(1))==checked$SHA256))
membership<-fread(file.path(parent,"tables/inner_membership.csv"),colClasses=c(Sample="character"))
files<-sort(list.files(file.path(parent,"checkpoints"),pattern="^nested_.*rds$",full.names=TRUE));stopifnot(length(files)==25,nrow(meta)==497,all(is.finite(X)))
audit<-list()
for(f in files){z<-readRDS(f);tr<-match(z$train,meta$Sample);tag<-sub(".rds$","",basename(f));cache_get(tr);m<-membership[Context==tag & Role=="Inner validation"];inner<-m$InnerFold[match(z$train,m$Sample)];stopifnot(!anyNA(inner));for(j in 1:4){s<-cache_get(tr[inner!=j]);audit[[length(audit)+1]]<-data.frame(Context=tag,InnerFold=j,TrainingN=sum(inner!=j),Cache=s$path,SHA256=digest(file=s$path,algo="sha256"))}}
fwrite(rbindlist(audit),file.path(out,"provenance/selection_cache_manifest.csv"));status("PREFLIGHT_PASS","25 outer splits and 100 inner cached selections verified; all support K=80")
run_split<-function(f){
 tag<-sub(".rds$","",basename(f));dest<-file.path(out,"checkpoints",basename(f));if(file.exists(dest))return(readRDS(dest))
 status("RUNNING",tag);z<-readRDS(f);train<-match(z$train,meta$Sample);test<-match(z$test,meta$Sample);stopifnot(!length(intersect(train,test)))
 m<-membership[Context==tag & Role=="Inner validation"];inner<-m$InnerFold[match(z$train,m$Sample)]
 preds<-array(NA_character_,dim=c(length(train),length(size_cfg$marker_counts),length(size_cfg$svm_cost)))
 pathlogs<-list();maxerr<-0
 for(j in 1:4){itr<-train[inner!=j];ite<-train[inner==j];s<-cache_get(itr);corrected<-correct_split(itr,ite)
  sc<-scale_split(corrected$train[,s$genes,drop=FALSE],corrected$test[,s$genes,drop=FALSE])
  for(ci in seq_along(size_cfg$svm_cost)){
   path<-rfe_path(sc$train,meta$Subtype[itr],size_cfg$svm_cost[ci],size_cfg$marker_counts)
   maxerr<-max(maxerr,path$max_decision_error)
   drops<-copy(path$drops);drops[,`:=`(Context=tag,Partition=paste0('Inner',j),Cost=size_cfg$svm_cost[ci])];pathlogs[[length(pathlogs)+1]]<-drops
   for(ki in seq_along(size_cfg$marker_counts)){K<-size_cfg$marker_counts[ki];genes<-path$sets[[as.character(K)]]
    model<-fit_model('SVM',sc$train[,genes,drop=FALSE],meta$Subtype[itr],size_cfg$svm_cost[ci],seed_for(paste(tag,j,K,ci,'rfe')))
    p<-predict_prob('SVM',model,sc$test[,genes,drop=FALSE]);preds[inner==j,ki,ci]<-classes[max.col(p,ties.method='first')]
   }
  }
 }

 stopifnot(!anyNA(preds));s<-cache_get(train);corrected<-correct_split(train,test);res<-list();tuning<-list();features<-list()
 sc<-scale_split(corrected$train[,s$genes,drop=FALSE],corrected$test[,s$genes,drop=FALSE]);outerpaths<-list()
 for(ki in seq_along(size_cfg$marker_counts)){K<-size_cfg$marker_counts[ki]
  score<-vapply(seq_along(size_cfg$svm_cost),function(ci)ba(meta$Subtype[train],preds[,ki,ci]),numeric(1));best<-which.max(score);C<-size_cfg$svm_cost[best];key<-as.character(C)
  if(is.null(outerpaths[[key]])){
   outerpaths[[key]]<-rfe_path(sc$train,meta$Subtype[train],C,size_cfg$marker_counts)
   maxerr<-max(maxerr,outerpaths[[key]]$max_decision_error)
   drops<-copy(outerpaths[[key]]$drops);drops[,`:=`(Context=tag,Partition='OuterTrain',Cost=C)];pathlogs[[length(pathlogs)+1]]<-drops
  }
  genes<-outerpaths[[key]]$sets[[as.character(K)]]
  model<-fit_model('SVM',sc$train[,genes,drop=FALSE],meta$Subtype[train],C,seed_for(paste(tag,K,'rfe-refit')));p<-predict_prob('SVM',model,sc$test[,genes,drop=FALSE])
  res[[ki]]<-data.frame(Context=tag,Repeat=z$predictions$Repeat[1],Fold=z$predictions$Fold[1],K=K,Sample=meta$Sample[test],Study=meta$Study[test],Truth=as.character(meta$Subtype[test]),Prediction=classes[max.col(p,ties.method='first')],ProbabilityC1=p[,1],ProbabilityC2=p[,2],ProbabilityC3=p[,3],Cost=C)
  tuning[[ki]]<-data.frame(Context=tag,K=K,Cost=size_cfg$svm_cost,InnerBalancedAccuracy=score,Selected=seq_along(score)==best)
  features[[ki]]<-data.frame(Context=tag,K=K,Cost=C,FeatureOrder=seq_along(genes),Gene=genes)
 }

 ans<-list(predictions=rbindlist(res),tuning=rbindlist(tuning),features=rbindlist(features),paths=rbindlist(pathlogs),max_decision_error=maxerr,train=z$train,test=z$test);atomic_save(ans,dest);cat("Completed",tag,"\n");flush.console();ans
}
if(stage=='smoke'){t<-system.time(z<-run_split(files[1]));write_json(list(passed=TRUE,seconds=unname(t['elapsed']),max_decision_error=z$max_decision_error),file.path(out,'review/smoke.json'),auto_unbox=TRUE,pretty=TRUE);print(t);quit(save='no')}
all<-parallel::mclapply(files,run_split,mc.cores=size_cfg$workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
if(any(vapply(all,inherits,logical(1),"try-error")))stop("A marker-count worker failed")
fwrite(rbindlist(lapply(all,`[[`,'paths')),file.path(out,'tables/rfe_deletion_paths.csv'));stopifnot(max(vapply(all,`[[`,numeric(1),'max_decision_error'))<1e-8)
raw<-rbindlist(lapply(all,`[[`,"predictions"));fwrite(raw,file.path(out,"tables/predictions_by_repeat.csv"));fwrite(rbindlist(lapply(all,`[[`,"tuning")),file.path(out,"tables/tuning.csv"));fwrite(rbindlist(lapply(all,`[[`,"features")),file.path(out,"tables/selected_features.csv"))
oof<-raw[,{stopifnot(.N==5,uniqueN(Repeat)==5);p<-colMeans(.SD);list(Study=Study[1],Truth=Truth[1],ProbabilityC1=p[1],ProbabilityC2=p[2],ProbabilityC3=p[3],Prediction=classes[which.max(p)])},by=.(K,Sample),.SDcols=paste0("Probability",classes)];fwrite(oof,file.path(out,"tables/oof_predictions.csv"))
metric<-function(z){rec<-vapply(classes,function(cl)mean(z$Prediction[z$Truth==cl]==cl),numeric(1));prec<-vapply(classes,function(cl)if(sum(z$Prediction==cl))sum(z$Prediction==cl&z$Truth==cl)/sum(z$Prediction==cl) else 0,numeric(1));f1<-ifelse(rec+prec>0,2*rec*prec/(rec+prec),0);auc<-vapply(classes,function(cl)as.numeric(pROC::auc(pROC::roc(as.integer(z$Truth==cl),z[[paste0("Probability",cl)]],levels=c(0,1),direction="<",quiet=TRUE))),numeric(1));c(MacroAUROC=mean(auc),BalancedAccuracy=mean(rec),MacroF1=mean(f1),setNames(rec,paste0(classes,"_Recall")))}
set.seed(size_cfg$seed);base<-oof[K==min(K)];setorder(base,Sample);groups<-split(seq_len(nrow(base)),interaction(base$Study,base$Truth,drop=TRUE));draws<-replicate(size_cfg$bootstrap_reps,unlist(lapply(groups,function(ii)ii[sample.int(length(ii),length(ii),TRUE)])),simplify=FALSE)
summary<-rbindlist(lapply(size_cfg$marker_counts,function(k){z<-oof[K==k];z<-z[match(base$Sample,Sample)];v<-metric(z);b<-t(vapply(draws,function(ix)metric(z[ix]),v));data.table(K=k,Metric=names(v),Estimate=as.numeric(v),Lower95=apply(b,2,quantile,.025),Upper95=apply(b,2,quantile,.975))}));fwrite(summary,file.path(out,"tables/marker_size_performance.csv"))
stopifnot(nrow(raw)==497*5*length(size_cfg$marker_counts))
stopifnot(nrow(oof)==497*length(size_cfg$marker_counts),all(abs(rowSums(oof[,paste0("Probability",classes),with=FALSE])-1)<1e-10))
write_json(list(passed=TRUE,outer_fits=25,marker_counts=size_cfg$marker_counts,n=497,predictions_per_K=497*5,scope=size_cfg$scope),file.path(out,"verification.json"),auto_unbox=TRUE,pretty=TRUE)
full<-fread(file.path(parent,'tables/nested_oof_predictions.csv'))[Model=='Consensus_SVM'];full<-full[match(base$Sample,Sample)]
fb<-vapply(draws,function(ix)metric(full[ix])['BalancedAccuracy'],numeric(1))
delta<-rbindlist(lapply(size_cfg$marker_counts,function(k){z<-oof[K==k];z<-z[match(base$Sample,Sample)];v<-metric(z)['BalancedAccuracy']-metric(full)['BalancedAccuracy'];b<-vapply(draws,function(ix)metric(z[ix])['BalancedAccuracy'],numeric(1))-fb;data.table(K=k,Estimate=unname(v),Lower95=quantile(b,.025),Upper95=quantile(b,.975))}))
fwrite(delta,file.path(out,'tables/paired_BA_difference.csv'))
bundle<-readRDS(file.path(parent,'models/candidate_models.rds'))$Consensus_SVM
w<-svm_pair_weights(bundle$model);d<-data.table(Gene=rownames(w),MeanSquaredWeight=rowMeans(w^2));d[,RelativeWeight:=MeanSquaredWeight/sum(MeanSquaredWeight)];setorder(d,-MeanSquaredWeight,Gene);d[,Rank:=seq_len(.N)]
fwrite(d,file.path(out,'tables/full_model_weight_ranking.csv'))
wl<-as.data.table(as.table(w));setnames(wl,c('Gene','Boundary','Weight'));fwrite(wl,file.path(out,'tables/full_model_pairwise_weights.csv'))
zz<-readRDS('8-Visualization/refinements/fig5_marker_story_20261030_v2/source_data/display_data.rds')$z[,bundle$genes,drop=FALSE]
ref<-attr(predict(bundle$model,zz,decision.values=TRUE),'decision.values');calc<-sweep(zz%*%w,2,bundle$model$rho,'-');err<-max(abs(ref-calc));stopifnot(err<1e-8)
write_json(list(passed=TRUE,genes=nrow(w),boundaries=ncol(w),max_decision_error=err),file.path(out,'review/full_model_weight_reconstruction.json'),auto_unbox=TRUE,pretty=TRUE)
status('COMPUTED_REVIEW','SVM-RFE sensitivity computed; no K selected; primary model unchanged')
