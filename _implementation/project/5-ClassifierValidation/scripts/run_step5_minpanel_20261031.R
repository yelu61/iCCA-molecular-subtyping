options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(digest);library(e1071);library(pROC);library(limma);library(randomForest);library(glmnet)})
source('5-ClassifierValidation/scripts/svm_rfe_helpers_20261031.R')
mc<-fromJSON('5-ClassifierValidation/config/minpanel_20261031.json');cfg<-fromJSON('5-ClassifierValidation/config/consensus_revision_20261030.json')
parent<-file.path('5-ClassifierValidation/runs',mc$parent_run);oldrfe<-file.path('5-ClassifierValidation/runs',mc$rfe_run);run<-file.path('5-ClassifierValidation/runs',mc$run_id)
classes<-c('C1','C2','C3');identity<-readLines(file.path(parent,'provenance/identity.txt'))
status<-function(state,detail){cat(format(Sys.time(),'%H:%M:%S'),state,detail,'\n');flush.console();write_json(list(state=state,detail=detail),file.path(run,'status.json'),auto_unbox=TRUE,pretty=TRUE)}
options(error=function(){status('FAILED',geterrmessage());quit(save='no',status=1)})
meta<-as.data.frame(fread(file.path(parent,'tables/population.csv'),colClasses=c(Sample='character')));meta$Subtype<-factor(meta$Subtype,levels=classes)
env<-new.env(parent=emptyenv());load('1-Clustering/data/combined_expr_before_correction.rds',env);X<-t(env$combined_expr[,meta$Sample,drop=FALSE]);storage.mode(X)<-'double';rm(env)
fun<-c('seed_for','atomic_save','make_folds','correct_split','limma_fit','scale_split','select_features','ba','align_prob','fit_model','predict_prob')
for(e in parse(file.path(parent,'provenance/run_step5_consensus_revision.R')))if(is.call(e)&&identical(e[[1]],as.name('<-'))&&as.character(e[[2]])[1]%in%fun)eval(e)
stopifnot(all(vapply(fun,exists,logical(1))),nrow(meta)==497,all(is.finite(X)))
source_id<-digest(list(mc,identity,digest(file='5-ClassifierValidation/scripts/run_step5_minpanel_20261031.R',algo='sha256'),digest(file='5-ClassifierValidation/scripts/svm_rfe_helpers_20261031.R',algo='sha256')),algo='sha256')
ip<-file.path(run,'provenance/identity.txt');if(file.exists(ip))stopifnot(identical(readLines(ip),source_id)) else writeLines(source_id,ip)
for(p in c('5-ClassifierValidation/config/minpanel_20261031.json','5-ClassifierValidation/scripts/run_step5_minpanel_20261031.R','5-ClassifierValidation/scripts/svm_rfe_helpers_20261031.R'))file.copy(p,file.path(run,'provenance',basename(p)),overwrite=FALSE)
inputs<-fread(file.path(parent,'provenance/input_manifest.csv'))[Role%in%c('expression','master','folds')];stopifnot(all(vapply(inputs$Path,function(p)digest(file=p,algo='sha256'),character(1))==inputs$SHA256))
fwrite(inputs,file.path(run,'provenance/input_manifest.csv'));model_sha<-digest(file=file.path(parent,'models/candidate_models.rds'),algo='sha256');writeLines(model_sha,file.path(run,'provenance/primary_model_sha256.txt'))
cache<-function(ii){key<-digest(list(identity,sort(meta$Sample[ii])),algo='sha256');p<-file.path(parent,'cache',paste0('selection_',key,'.rds'));s<-if(file.exists(p))readRDS(p) else select_features(ii);stopifnot(identical(s$training_samples,sort(meta$Sample[ii])));d<-as.data.table(s$support)[Methods>=2];setorder(d,-Methods,FRank,Gene);list(genes=d$Gene,parent_genes=s$genes,training_samples=s$training_samples)}
paths<-fread(file.path(oldrfe,'tables/rfe_deletion_paths.csv'));member<-fread(file.path(parent,'tables/inner_membership.csv'),colClasses=c(Sample='character'))
recover_path<-function(genes,tag,partition,C,tr,y){
 d<-paths[Context==tag & Partition==partition & Cost==C];setorder(d,-NBefore)
 if(nrow(d)){
  stopifnot(nrow(d)==length(genes)-2,identical(d$NBefore,seq.int(length(genes),3L)),!anyDuplicated(d$Dropped),all(d$Dropped%in%genes))
  sets<-list();g<-genes;sets[[as.character(length(g))]]<-g
  for(i in seq_len(nrow(d))){stopifnot(length(g)==d$NBefore[i],d$Dropped[i]%in%g);g<-g[g!=d$Dropped[i]];sets[[as.character(length(g))]]<-g}
  list(sets=sets,drops=d,max_decision_error=0,reused=TRUE)
 }else{p<-rfe_path(tr,y,C,2:length(genes));p$reused<-FALSE;p}
}
recalls<-function(y,p)setNames(vapply(classes,function(cl)mean(p[y==cl]==cl),numeric(1)),paste0(classes,'Recall'))
choose_run<-function(tag,train,test=NULL,inner,reference=NULL){
 dest<-file.path(run,'checkpoints',paste0(tag,'.rds'));if(file.exists(dest))return(readRDS(dest))
 status('RUNNING',tag);sel<-cache(train);splits<-list();basepred<-matrix(NA_character_,length(train),length(mc$costs));audit<-list();path_audit<-list()
 for(j in 1:4){itr<-train[inner!=j];ite<-train[inner==j];s<-cache(itr);z<-correct_split(itr,ite);sc<-scale_split(z$train[,s$genes,drop=FALSE],z$test[,s$genes,drop=FALSE]);plist<-list()
  for(ci in seq_along(mc$costs)){
   plist[[ci]]<-recover_path(s$genes,tag,paste0('Inner',j),mc$costs[ci],sc$train,meta$Subtype[itr])
   path_audit[[length(path_audit)+1]]<-data.table(Context=tag,Partition=paste0('Inner',j),Cost=mc$costs[ci],InitialN=length(s$genes),Reused=plist[[ci]]$reused,DecisionError=plist[[ci]]$max_decision_error)
   if(is.null(reference)){
    order<-match(s$parent_genes,s$genes);fit<-fit_model('SVM',sc$train[,order,drop=FALSE],meta$Subtype[itr],mc$costs[ci],seed_for(paste(tag,j,'Consensus_SVM',ci)))
    pp<-predict_prob('SVM',fit,sc$test[,order,drop=FALSE]);basepred[inner==j,ci]<-classes[max.col(pp,ties.method='first')]
   }
  }
  splits[[j]]<-list(train=itr,test=ite,sc=sc,paths=plist)
  audit[[j]]<-data.table(Context=tag,InnerFold=j,Sample=meta$Sample[train],Role=ifelse(inner==j,'Inner validation','Inner training'))
 }
 if(is.null(reference)){scores<-apply(basepred,2,function(p)ba(meta$Subtype[train],p));reference<-data.table(Parameter=mc$costs,InnerBalancedAccuracy=scores);reference[,Selected:=seq_len(.N)==which.max(scores)]}
 baseline<-max(reference$InnerBalancedAccuracy);baseline_cost<-reference$Parameter[which.max(reference$InnerBalancedAccuracy)]
 cap<-min(length(sel$genes),vapply(splits,function(s)ncol(s$sc$train),integer(1)));rows<-list();prediction_logs<-list();choice<-NULL
 for(K in seq.int(mc$minimum_k,cap)){
  preds<-matrix(NA_character_,length(train),length(mc$costs))
  for(j in 1:4){s<-splits[[j]];for(ci in seq_along(mc$costs)){g<-s$paths[[ci]]$sets[[as.character(K)]]
   fit<-fit_model('SVM',s$sc$train[,g,drop=FALSE],meta$Subtype[s$train],mc$costs[ci],seed_for(paste(tag,j,K,ci,'rfe')))
   pp<-predict_prob('SVM',fit,s$sc$test[,g,drop=FALSE]);preds[inner==j,ci]<-classes[max.col(pp,ties.method='first')]
  }}
  stopifnot(!anyNA(preds));rec<-t(apply(preds,2,function(p)recalls(meta$Subtype[train],p)))
  rr<-data.table(Context=tag,K=K,Cost=mc$costs,InnerBA=rowMeans(rec),C1Recall=rec[,1],C2Recall=rec[,2],C3Recall=rec[,3],FullInnerBA=baseline,FeasibleMaxK=cap)
  rr[,Eligible:=InnerBA>=baseline-mc$max_ba_loss-1e-12 & pmin(C1Recall,C2Recall,C3Recall)>=mc$min_class_recall-1e-12];rows[[length(rows)+1]]<-rr
  prediction_logs[[length(prediction_logs)+1]]<-rbindlist(lapply(seq_along(mc$costs),function(ci)data.table(Context=tag,K=K,Cost=mc$costs[ci],Sample=meta$Sample[train],InnerFold=inner,Truth=as.character(meta$Subtype[train]),Prediction=preds[,ci])))
  eligible<-rr[Eligible==TRUE];if(nrow(eligible)){setorder(eligible,-InnerBA,Cost);choice<-eligible[1];break}
 }
 fallback<-is.null(choice);K<-if(fallback)length(sel$genes) else choice$K;C<-if(fallback)baseline_cost else choice$Cost
 z<-correct_split(train,if(is.null(test))train else test);sc<-scale_split(z$train[,sel$genes,drop=FALSE],z$test[,sel$genes,drop=FALSE])
 if(fallback){g<-sel$parent_genes;fitseed<-seed_for(paste(tag,'Consensus_SVM','refit'))}else{
  pp<-recover_path(sel$genes,tag,'OuterTrain',C,sc$train,meta$Subtype[train]);g<-pp$sets[[as.character(K)]];fitseed<-seed_for(paste(tag,K,'rfe-refit'))
  path_audit[[length(path_audit)+1]]<-data.table(Context=tag,Partition='OuterTrain',Cost=C,InitialN=length(sel$genes),Reused=pp$reused,DecisionError=pp$max_decision_error)
 }
 fit<-fit_model('SVM',sc$train[,g,drop=FALSE],meta$Subtype[train],C,fitseed)
 result<-list(choice=data.table(Context=tag,K=K,Cost=C,Fallback=fallback,InitialN=length(sel$genes),FeasibleMaxK=cap,FullInnerBA=baseline,SelectedInnerBA=if(fallback)baseline else choice$InnerBA),tuning=rbindlist(rows),inner_predictions=rbindlist(prediction_logs),features=data.table(Context=tag,Rank=seq_along(g),Gene=g),train=meta$Sample[train],test=meta$Sample[test],membership=rbindlist(audit),path_audit=rbindlist(path_audit))
 if(!is.null(test)){
  p<-predict_prob('SVM',fit,sc$test[,g,drop=FALSE]);result$predictions<-data.table(Context=tag,Sample=meta$Sample[test],Study=meta$Study[test],Truth=as.character(meta$Subtype[test]),Prediction=classes[max.col(p,ties.method='first')],ProbabilityC1=p[,1],ProbabilityC2=p[,2],ProbabilityC3=p[,3],K=K,Cost=C,Fallback=fallback)
 }else{
  bundle<-list(model=fit,genes=g,center=sc$center[g],scale=sc$scale[g],study_offsets=z$offsets,global=z$global,classes=classes,kind='SVM',parameter=C,training_samples=meta$Sample[train],identity=source_id,status='EXPLORATORY_CANDIDATE_NOT_PRIMARY')
  atomic_save(bundle,file.path(run,'models/compact_candidate.rds'));w<-svm_pair_weights(fit);weights<-data.table(Gene=rownames(w),MeanSquaredWeight=rowMeans(w^2));weights[,RelativeWeight:=MeanSquaredWeight/sum(MeanSquaredWeight)];setorder(weights,-MeanSquaredWeight,Gene);weights[,Rank:=seq_len(.N)];fwrite(weights,file.path(run,'tables/full_candidate_weights.csv'))
 }
 atomic_save(result,dest);cat('Selected',tag,'K',K,'C',C,'fallback',fallback,'\n');flush.console();result
}
files<-sort(list.files(file.path(parent,'checkpoints'),pattern='^nested_.*rds$',full.names=TRUE));stopifnot(length(files)==25)
work<-function(f){z<-readRDS(f);tag<-sub('[.]rds$','',basename(f));train<-match(z$train,meta$Sample);test<-match(z$test,meta$Sample);m<-member[Context==tag & Role=='Inner validation'];inner<-m$InnerFold[match(z$train,m$Sample)];stopifnot(!anyNA(inner),!length(intersect(train,test)));ref<-as.data.table(z$tuning)[Model=='Consensus_SVM'];ans<-choose_run(tag,train,test,inner,ref);ans}
stage<-if(length(commandArgs(TRUE)))commandArgs(TRUE)[1] else 'run'
if(stage=='smoke'){tt<-system.time(work(files[1]));write_json(list(passed=TRUE,seconds=unname(tt['elapsed'])),file.path(run,'review/smoke.json'),auto_unbox=TRUE);print(tt);quit(save='no')}
all<-parallel::mclapply(files,work,mc.cores=mc$workers,mc.preschedule=FALSE,mc.set.seed=FALSE);stopifnot(!any(vapply(all,inherits,logical(1),'try-error')))
for(n in c('choice','tuning','inner_predictions','features','membership','path_audit'))fwrite(rbindlist(lapply(all,`[[`,n)),file.path(run,'tables',paste0(n,'.csv')))
raw<-rbindlist(lapply(all,`[[`,'predictions'));raw[,Repeat:=as.integer(sub('nested_R([0-9]+)_.*','\\1',Context))];fwrite(raw,file.path(run,'tables/predictions_by_repeat.csv'))
oof<-raw[,{stopifnot(.N==5,uniqueN(Repeat)==5);p<-colMeans(.SD);list(Study=Study[1],Truth=Truth[1],ProbabilityC1=p[1],ProbabilityC2=p[2],ProbabilityC3=p[3],Prediction=classes[which.max(p)])},by=Sample,.SDcols=paste0('Probability',classes)];setorder(oof,Sample)
stopifnot(nrow(oof)==497,nrow(raw)==2485,all(abs(rowSums(oof[,paste0('Probability',classes),with=FALSE])-1)<1e-10));fwrite(oof,file.path(run,'tables/oof_predictions.csv'))
metric<-function(z){rec<-recalls(z$Truth,z$Prediction);prec<-vapply(classes,function(cl)if(sum(z$Prediction==cl))mean(z$Truth[z$Prediction==cl]==cl) else 0,numeric(1));f1<-ifelse(rec+prec>0,2*rec*prec/(rec+prec),0);auc<-vapply(classes,function(cl)as.numeric(pROC::auc(pROC::roc(as.integer(z$Truth==cl),z[[paste0('Probability',cl)]],levels=c(0,1),direction='<',quiet=TRUE))),numeric(1));c(Accuracy=mean(z$Truth==z$Prediction),BalancedAccuracy=mean(rec),MacroF1=mean(f1),MacroAUROC=mean(auc),rec)}
full<-fread(file.path(parent,'tables/nested_oof_predictions.csv'),colClasses=c(Sample='character'))[Model=='Consensus_SVM'];full<-full[match(oof$Sample,Sample)];stopifnot(identical(full$Truth,oof$Truth))
set.seed(mc$seed);gr<-split(seq_len(nrow(oof)),interaction(oof$Study,oof$Truth,drop=TRUE));draws<-replicate(mc$bootstrap_reps,unlist(lapply(gr,function(ii)ii[sample.int(length(ii),length(ii),TRUE)])),simplify=FALSE)
vals<-list(Full=metric(full),Compact=metric(oof));bs<-lapply(list(Full=full,Compact=oof),function(z)t(vapply(draws,function(ii)metric(z[ii]),vals$Full)))
perf<-rbindlist(lapply(names(vals),function(n)data.table(Model=n,Metric=names(vals[[n]]),Estimate=as.numeric(vals[[n]]),Lower95=apply(bs[[n]],2,quantile,.025),Upper95=apply(bs[[n]],2,quantile,.975))));fwrite(perf,file.path(run,'tables/performance_ci.csv'))
delta<-bs$Compact-bs$Full;fwrite(data.table(Metric=names(vals$Full),Estimate=as.numeric(vals$Compact-vals$Full),Lower95=apply(delta,2,quantile,.025),Upper95=apply(delta,2,quantile,.975)),file.path(run,'tables/paired_differences.csv'))
status('FULL_CANDIDATE','Selecting a separate full-development candidate with four training-only inner folds')
inner<-make_folds(meta$Subtype,meta$Study,4,mc$full_inner_seed);final<-choose_run('full_candidate',seq_len(nrow(meta)),inner=inner)
for(n in c('choice','tuning','inner_predictions','features','membership','path_audit'))fwrite(final[[n]],file.path(run,'tables',paste0('full_',n,'.csv')))
stopifnot(digest(file=file.path(parent,'models/candidate_models.rds'),algo='sha256')==model_sha)
write_json(list(passed=TRUE,outer_splits=25,n=497,user_ba_loss=mc$max_ba_loss,user_recall_floor=mc$min_class_recall,full_candidate_K=final$choice$K,primary_model_unchanged=TRUE),file.path(run,'verification.json'),auto_unbox=TRUE,pretty=TRUE)
capture.output(sessionInfo(),file=file.path(run,'provenance/sessionInfo.txt'))
status('COMPUTED_REVIEW','Adaptive selection evaluated; separate compact candidate saved;114-gene primary unchanged')
