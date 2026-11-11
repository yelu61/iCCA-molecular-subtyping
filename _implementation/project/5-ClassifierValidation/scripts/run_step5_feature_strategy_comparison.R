#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(digest);library(e1071);library(pROC)})
scfg<-fromJSON("5-ClassifierValidation/config/feature_strategy_comparison_20261030.json")
parent<-file.path("5-ClassifierValidation/runs",scfg$parent_run)
cfg<-fromJSON(file.path(parent,"provenance/consensus_revision_20261030.json"))
out<-file.path("5-ClassifierValidation/runs",scfg$run_id)
for(d in c("tables","checkpoints","provenance"))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
classes<-c("C1","C2","C3");identity<-readLines(file.path(parent,"provenance/identity.txt"))
source_file<-"5-ClassifierValidation/scripts/run_step5_feature_strategy_comparison.R"
config_file<-"5-ClassifierValidation/config/feature_strategy_comparison_20261030.json"
run_identity<-digest(list(scfg,identity,digest(file=source_file,algo="sha256")),algo="sha256")
idf<-file.path(out,"provenance/identity.txt")
if(file.exists(idf))stopifnot(identical(readLines(idf),run_identity)) else writeLines(run_identity,idf)
for(p in c(source_file,config_file))if(!file.exists(file.path(out,"provenance",basename(p))))stopifnot(file.copy(p,file.path(out,"provenance",basename(p))))
status<-function(state,detail){tmp<-file.path(out,paste0("status.",Sys.getpid(),".tmp"));write_json(list(state=state,detail=detail,time=format(Sys.time(),tz="UTC",usetz=TRUE)),tmp,auto_unbox=TRUE,pretty=TRUE);stopifnot(file.rename(tmp,file.path(out,"status.json")));cat(format(Sys.time(),"%H:%M:%S"),state,detail,"\n");flush.console()}
options(error=function(){status("FAILED",geterrmessage());quit(save="no",status=1)})
meta<-as.data.frame(fread(file.path(parent,"tables/population.csv"),colClasses=c(Sample="character")));meta$Subtype<-factor(meta$Subtype,levels=classes)
env<-new.env(parent=emptyenv());load("1-Clustering/data/combined_expr_before_correction.rds",env);X<-t(env$combined_expr[,meta$Sample,drop=FALSE]);storage.mode(X)<-"double";rm(env)
functions<-c("seed_for","atomic_save","correct_split","scale_split","ba","align_prob","fit_model","predict_prob")
for(e in parse(file.path(parent,"provenance/run_step5_consensus_revision.R")))if(is.call(e)&&identical(e[[1]],as.name("<-"))&&as.character(e[[2]])[1]%in%functions)eval(e)
stopifnot(nrow(meta)==497,all(is.finite(X)),all(vapply(functions,exists,logical(1))))
input_manifest<-fread(file.path(parent,"provenance/input_manifest.csv"))
stopifnot(all(vapply(seq_len(nrow(input_manifest)),function(i)identical(digest(file=input_manifest$Path[i],algo="sha256"),input_manifest$SHA256[i]),logical(1))))
fwrite(input_manifest,file.path(out,"provenance/input_manifest.csv"))
membership<-fread(file.path(parent,"tables/inner_membership.csv"),colClasses=c(Sample="character"))
specs<-as.data.table(expand.grid(Budget=scfg$budgets,Strategy=scfg$strategies,stringsAsFactors=FALSE));specs[,Spec:=paste(Budget,Strategy,sep="_")]
cache_get<-function(indices){
 key<-digest(list(identity,sort(meta$Sample[indices])),algo="sha256");p<-file.path(parent,"cache",paste0("selection_",key,".rds"));s<-readRDS(p)
 stopifnot(identical(s$training_samples,sort(meta$Sample[indices])))
 d<-as.data.table(s$support);setorder(d,FRank,Gene)
 sets<-list(Limma=d[Limma==TRUE,Gene],RF=d[RF==TRUE,Gene],LASSO=d[LASSO==TRUE,Gene],Consensus2=d[Methods>=2,Gene],Core3=d[Methods==3,Gene])
 stopifnot(identical(sets$Consensus2,s$genes),min(lengths(sets))>=2)
 list(sets=sets,path=p,matched=length(sets$Core3))
}
genes_for<-function(s,strategy,budget){g<-s$sets[[strategy]];if(budget=="CoreMatched")head(g,s$matched) else g}
files<-sort(list.files(file.path(parent,"checkpoints"),pattern="^nested_.*rds$",full.names=TRUE));stopifnot(length(files)==25)
audit<-list()
for(f in files){z<-readRDS(f);tag<-sub(".rds$","",basename(f));tr<-match(z$train,meta$Sample);m<-membership[Context==tag & Role=="Inner validation"];inner<-m$InnerFold[match(z$train,m$Sample)];stopifnot(!anyNA(inner),!length(intersect(z$train,z$test)))
 for(j in 0:4){ii<-if(j==0)tr else tr[inner!=j];s<-cache_get(ii);audit[[length(audit)+1]]<-data.table(Context=tag,InnerFold=j,TrainingN=length(ii),MatchedN=s$matched,Cache=s$path,SHA256=digest(file=s$path,algo="sha256"))}}
fwrite(rbindlist(audit),file.path(out,"provenance/selection_cache_manifest.csv"));fwrite(specs,file.path(out,"tables/specifications.csv"))
status("PREFLIGHT_PASS","125 training-bound selectors; five strategies at native and partition-specific matched counts; fixed parent splits")
run_split<-function(f){
 tag<-sub(".rds$","",basename(f));dest<-file.path(out,"checkpoints",basename(f));if(file.exists(dest))return(readRDS(dest))
 z<-readRDS(f);tr<-match(z$train,meta$Sample);te<-match(z$test,meta$Sample);m<-membership[Context==tag & Role=="Inner validation"];inner<-m$InnerFold[match(z$train,m$Sample)]
 hard<-array(NA_character_,dim=c(length(tr),nrow(specs),length(scfg$svm_cost)));inner_features<-list()
 for(j in 1:4){itr<-tr[inner!=j];ival<-tr[inner==j];s<-cache_get(itr);corrected<-correct_split(itr,ival);fits<-new.env(parent=emptyenv())
  for(si in seq_len(nrow(specs))){sp<-specs[si];g<-genes_for(s,sp$Strategy,sp$Budget);sc<-scale_split(corrected$train[,g,drop=FALSE],corrected$test[,g,drop=FALSE]);inner_features[[length(inner_features)+1]]<-data.table(Context=tag,InnerFold=j,Spec=sp$Spec,Rank=seq_along(g),Gene=g)
   for(ci in seq_along(scfg$svm_cost)){key<-digest(list(g,scfg$svm_cost[ci]),algo="sha256")
    if(!exists(key,fits,inherits=FALSE)){
     seed<-if(sp$Spec=="Native_Consensus2")seed_for(paste(tag,j,"Consensus_SVM",ci)) else seed_for(paste(tag,j,key,"strategy"))
     fit<-fit_model("SVM",sc$train,meta$Subtype[itr],scfg$svm_cost[ci],seed);p<-predict_prob("SVM",fit,sc$test);assign(key,classes[max.col(p,ties.method="first")],fits)
    };hard[inner==j,si,ci]<-get(key,fits,inherits=FALSE)
   }
  }
 }
 stopifnot(!anyNA(hard));s<-cache_get(tr);corrected<-correct_split(tr,te);pred<-tuning<-features<-list();fits<-new.env(parent=emptyenv())
 for(si in seq_len(nrow(specs))){sp<-specs[si];g<-genes_for(s,sp$Strategy,sp$Budget);scores<-vapply(seq_along(scfg$svm_cost),function(ci)ba(meta$Subtype[tr],hard[,si,ci]),numeric(1));best<-which.max(scores);cost<-scfg$svm_cost[best]
  sc<-scale_split(corrected$train[,g,drop=FALSE],corrected$test[,g,drop=FALSE]);key<-digest(list(g,cost),algo="sha256")
  if(!exists(key,fits,inherits=FALSE)){
   seed<-if(sp$Spec=="Native_Consensus2")seed_for(paste(tag,"Consensus_SVM","refit")) else seed_for(paste(tag,key,"strategy-refit"))
   fit<-fit_model("SVM",sc$train,meta$Subtype[tr],cost,seed);assign(key,predict_prob("SVM",fit,sc$test),fits)
  };p<-get(key,fits,inherits=FALSE)
  pred[[si]]<-data.table(Context=tag,Repeat=z$predictions$Repeat[1],Fold=z$predictions$Fold[1],Budget=sp$Budget,Strategy=sp$Strategy,Spec=sp$Spec,Sample=meta$Sample[te],Study=meta$Study[te],Truth=as.character(meta$Subtype[te]),Prediction=classes[max.col(p,ties.method="first")],ProbabilityC1=p[,1],ProbabilityC2=p[,2],ProbabilityC3=p[,3],Features=length(g),Cost=cost)
  tuning[[si]]<-data.table(Context=tag,Spec=sp$Spec,Cost=scfg$svm_cost,InnerBalancedAccuracy=scores,Selected=seq_along(scores)==best)
  features[[si]]<-data.table(Context=tag,Spec=sp$Spec,Budget=sp$Budget,Strategy=sp$Strategy,Rank=seq_along(g),Gene=g)
 }
 ans<-list(predictions=rbindlist(pred),tuning=rbindlist(tuning),features=rbindlist(features),inner_features=rbindlist(inner_features),inner_hard=hard,train=z$train,test=z$test)
 atomic_save(ans,dest);cat("Completed",tag,"\n");flush.console();ans
}
all<-parallel::mclapply(files,run_split,mc.cores=scfg$workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
stopifnot(!any(vapply(all,inherits,logical(1),"try-error")))
for(n in c("predictions","tuning","features","inner_features"))fwrite(rbindlist(lapply(all,`[[`,n)),file.path(out,"tables",paste0(n,".csv")))
raw<-rbindlist(lapply(all,`[[`,"predictions"))
oof<-raw[,{stopifnot(.N==5,uniqueN(Repeat)==5);p<-colMeans(.SD);list(Study=Study[1],Truth=Truth[1],ProbabilityC1=p[1],ProbabilityC2=p[2],ProbabilityC3=p[3],Prediction=classes[which.max(p)])},by=.(Spec,Budget,Strategy,Sample),.SDcols=paste0("Probability",classes)]
fwrite(oof,file.path(out,"tables/oof_predictions.csv"))
parent_oof<-fread(file.path(parent,"tables/nested_oof_predictions.csv"),colClasses=c(Sample="character"))[Model=="Consensus_SVM"]
primary<-oof[Spec=="Native_Consensus2"];primary<-primary[match(parent_oof$Sample,Sample)]
primary_error<-max(abs(as.matrix(primary[,paste0("Probability",classes),with=FALSE])-as.matrix(parent_oof[,paste0("Probability",classes),with=FALSE])))
stopifnot(primary_error<1e-10)
metric<-function(z){rec<-vapply(classes,function(c)mean(z$Prediction[z$Truth==c]==c),numeric(1));f1<-vapply(classes,function(c)2*sum(z$Truth==c & z$Prediction==c)/(sum(z$Truth==c)+sum(z$Prediction==c)),numeric(1));aucs<-vapply(classes,function(c)as.numeric(pROC::auc(pROC::roc(as.integer(z$Truth==c),z[[paste0("Probability",c)]],levels=c(0,1),direction="<",quiet=TRUE))),numeric(1));c(Accuracy=mean(z$Truth==z$Prediction),BalancedAccuracy=mean(rec),MacroF1=mean(f1),MacroAUROC=mean(aucs))}
base<-oof[Spec=="Native_Consensus2"];setorder(base,Sample);set.seed(scfg$seed);groups<-split(seq_len(nrow(base)),interaction(base$Study,base$Truth,drop=TRUE));draws<-replicate(scfg$bootstrap_reps,unlist(lapply(groups,function(ii)ii[sample.int(length(ii),length(ii),TRUE)])),simplify=FALSE)
est<-boots<-list();perf<-list()
for(si in seq_len(nrow(specs))){sp<-specs[si];z<-oof[Spec==sp$Spec];z<-z[match(base$Sample,Sample)];v<-metric(z);b<-t(vapply(draws,function(ix)metric(z[ix]),v));est[[sp$Spec]]<-v;boots[[sp$Spec]]<-b;perf[[si]]<-data.table(Spec=sp$Spec,Budget=sp$Budget,Strategy=sp$Strategy,N=497,Metric=names(v),Estimate=as.numeric(v),Lower95=apply(b,2,quantile,.025),Upper95=apply(b,2,quantile,.975))}
fwrite(rbindlist(perf),file.path(out,"tables/performance_ci.csv"))
delta<-list();for(budget in scfg$budgets)for(strategy in setdiff(scfg$strategies,"Consensus2")){ref<-paste(budget,"Consensus2",sep="_");cmp<-paste(budget,strategy,sep="_");v<-est[[ref]]-est[[cmp]];b<-boots[[ref]]-boots[[cmp]];delta[[length(delta)+1]]<-data.table(Budget=budget,Primary=ref,Comparator=cmp,Metric=names(v),Difference=as.numeric(v),Lower95=apply(b,2,quantile,.025),Upper95=apply(b,2,quantile,.975))}
fwrite(rbindlist(delta),file.path(out,"tables/paired_differences.csv"))
features<-rbindlist(lapply(all,`[[`,"features"));sizes<-features[,.(Features=.N),by=.(Context,Spec,Budget,Strategy)];fwrite(sizes,file.path(out,"tables/outer_feature_counts.csv"))
fwrite(sizes[,.(MedianFeatures=median(Features),MinFeatures=min(Features),MaxFeatures=max(Features)),by=.(Spec,Budget,Strategy)],file.path(out,"tables/feature_count_summary.csv"))
stopifnot(nrow(raw)==497*5*10,nrow(oof)==497*10,all(abs(rowSums(oof[,paste0("Probability",classes),with=FALSE])-1)<1e-10))
write_json(list(passed=TRUE,n=497,specifications=10,outer_splits=25,inner_splits=100,primary_probability_error=primary_error,scope=scfg$scope),file.path(out,"verification.json"),auto_unbox=TRUE,pretty=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,"provenance/sessionInfo.txt"))
status("COMPUTED_REVIEW","Five screening strategies at native and count-aligned budgets; primary predictions reproduce parent; no new selected model")
