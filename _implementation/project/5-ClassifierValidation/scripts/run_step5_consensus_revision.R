#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE, warn=1)
suppressPackageStartupMessages({library(limma);library(randomForest);library(glmnet);library(e1071);library(jsonlite);library(digest);library(pROC)})
args <- commandArgs(TRUE)
stage <- if(length(args)) args[1] else "preflight"
stopifnot(stage %in% c("preflight","smoke","run","finalize"))
cfg_path <- "5-ClassifierValidation/config/consensus_revision_20261030.json"
cfg <- fromJSON(cfg_path)
run <- file.path("5-ClassifierValidation/runs",cfg$run_id)
for(d in c("tables","cache","checkpoints","models","provenance","logs")) dir.create(file.path(run,d),recursive=TRUE,showWarnings=FALSE)
classes <- c("C1","C2","C3")
write_table <- function(x,n) write.csv(x,file.path(run,"tables",paste0(n,".csv")),row.names=FALSE,na="")
status <- function(state,detail) {write_json(list(state=state,detail=detail,time=format(Sys.time(),tz="UTC",usetz=TRUE)),file.path(run,"status.json"),auto_unbox=TRUE,pretty=TRUE);cat(format(Sys.time(),"%H:%M:%S"),state,detail,"\n");flush.console()}
options(error=function(){status("FAILED",geterrmessage());traceback(10);quit(save="no",status=1,runLast=FALSE)})
seed_for <- function(x) strtoi(substr(digest(x,algo="xxhash32"),1,7),16L)
atomic_save <- function(x,path) {tmp<-paste0(path,".",Sys.getpid(),".tmp");saveRDS(x,tmp,compress=FALSE);stopifnot(file.rename(tmp,path))}
inputs <- c(expression="1-Clustering/data/combined_expr_before_correction.rds",master="governance/data/subtype_master_table.csv",folds="5-ClassifierValidation/publication/tables/step5_nested_cv_fold_definitions.csv",microarray1="0-Data/qc/0-Data-rnaseq-microarray1.Rdata",microarray2="0-Data/qc/0-Data-rnaseq-microarray2.Rdata")
stopifnot(all(file.exists(inputs)))
hashes <- vapply(inputs,function(p)digest(file=p,algo="sha256"),character(1))
identity <- digest(list(cfg,hashes),algo="sha256")
identity_file<-file.path(run,"provenance","identity.txt")
if(file.exists(identity_file)) stopifnot(identical(readLines(identity_file),identity)) else writeLines(identity,identity_file)
e<-new.env(parent=emptyenv());load(inputs[["expression"]],e);expr<-e$combined_expr
master<-read.csv(inputs[["master"]],check.names=FALSE,colClasses=c(Sample="character"))
meta<-master[master$SubtypeClass=="Main",c("Sample","Study","FinalSubtype")];names(meta)[3]<-"Subtype"
stopifnot(nrow(meta)==497,!anyDuplicated(meta$Sample),!anyDuplicated(rownames(expr)),!anyDuplicated(colnames(expr)),all(meta$Sample %in% colnames(expr)))
X<-t(as.matrix(expr[,meta$Sample,drop=FALSE]));storage.mode(X)<-"double";rm(expr,e)
stopifnot(all(is.finite(X)),all(meta$Subtype %in% classes));meta$Subtype<-factor(meta$Subtype,levels=classes)
folds<-read.csv(inputs[["folds"]],colClasses=c(Sample="character"))
stopifnot(setequal(unique(folds$Sample),meta$Sample),all(folds$Subtype==as.character(meta$Subtype[match(folds$Sample,meta$Sample)])))
make_folds<-function(y,study,k,seed) {set.seed(seed);out<-integer(length(y));for(z in levels(interaction(study,y,drop=TRUE))){ii<-which(interaction(study,y,drop=TRUE)==z);out[ii]<-sample(rep(seq_len(k),length.out=length(ii)))};out}
correct_split<-function(train,test) {
  stopifnot(length(train)>0,length(test)>0)
  tr<-X[train,,drop=FALSE];te<-X[test,,drop=FALSE];mu<-colMeans(tr)
  offsets<-lapply(split(seq_along(train),meta$Study[train]),function(ii)colMeans(tr[ii,,drop=FALSE])-mu)
  for(st in unique(meta$Study[train])){ii<-which(meta$Study[train]==st);tr[ii,]<-sweep(tr[ii,,drop=FALSE],2,offsets[[st]],"-")}
  unseen<-setdiff(unique(meta$Study[test]),names(offsets))
  for(st in unique(meta$Study[test])){ii<-which(meta$Study[test]==st);off<-if(st %in% names(offsets))offsets[[st]] else colMeans(te[ii,,drop=FALSE])-mu;te[ii,]<-sweep(te[ii,,drop=FALSE],2,off,"-")}
  list(train=tr,test=te,offsets=offsets,global=mu,unseen=unseen)
}
limma_fit<-function(x,y,study) {
  y<-factor(y,levels=classes);study<-droplevels(factor(study));design<-model.matrix(~0+y+study)
  colnames(design)<-make.names(colnames(design));stopifnot(qr(design)$rank==ncol(design))
  con<-makeContrasts(C1_vs_C2=yC1-yC2,C1_vs_C3=yC1-yC3,C2_vs_C3=yC2-yC3,levels=design)
  eBayes(contrasts.fit(lmFit(t(x),design),con),trend=TRUE,robust=TRUE)
}
scale_split<-function(tr,te){mu<-colMeans(tr);sdv<-apply(tr,2,sd);sdv[!is.finite(sdv)|sdv<1e-8]<-1;list(train=sweep(sweep(tr,2,mu,"-"),2,sdv,"/"),test=sweep(sweep(te,2,mu,"-"),2,sdv,"/"),center=mu,scale=sdv)}
select_features<-function(train) {
  key<-digest(list(identity,sort(meta$Sample[train])),algo="sha256");path<-file.path(run,"cache",paste0("selection_",key,".rds"))
  if(file.exists(path))return(readRDS(path))
  t0<-proc.time()[3];z<-correct_split(train,train);y<-meta$Subtype[train];x<-z$train
  fit<-limma_fit(x,y,meta$Study[train]);ranked<-rownames(topTable(fit,coef=1:3,number=Inf,sort.by="F"))
  lim<-unique(unlist(lapply(1:3,function(j){d<-topTable(fit,coef=j,number=Inf,sort.by="none");d<-d[is.finite(d$adj.P.Val)&d$adj.P.Val<cfg$limma_fdr,,drop=FALSE];up<-d[d$logFC>0,,drop=FALSE];dn<-d[d$logFC<0,,drop=FALSE];c(head(rownames(up)[order(-up$logFC,rownames(up))],cfg$limma_top_per_direction),head(rownames(dn)[order(dn$logFC,rownames(dn))],cfg$limma_top_per_direction))})))
  set.seed(seed_for(paste(key,"RF")));rf<-randomForest(x=x,y=y,ntree=cfg$rf_trees,importance=TRUE,keep.forest=FALSE)
  imp<-importance(rf,type=1);stopifnot(identical(colnames(imp),"MeanDecreaseAccuracy"));rf_genes<-head(rownames(imp)[order(-imp[,1],rownames(imp))],cfg$rf_top)
  lf<-make_folds(y,meta$Study[train],cfg$lasso_folds,seed_for(paste(key,"LASSO")))
  set.seed(seed_for(paste(key,"GLMNET")))
  lasso<-cv.glmnet(x,y,family="multinomial",alpha=cfg$lasso_alpha,type.multinomial="grouped",foldid=lf,type.measure="deviance",standardize=TRUE,parallel=FALSE)
  cf<-coef(lasso,s=cfg$lasso_lambda);las<-unique(unlist(lapply(cf,function(m)rownames(m)[as.vector(m)!=0 & rownames(m)!="(Intercept)"])))
  support<-data.frame(Gene=colnames(x),Limma=colnames(x)%in%lim,RF=colnames(x)%in%rf_genes,LASSO=colnames(x)%in%las)
  support$Methods<-rowSums(support[,c("Limma","RF","LASSO")]);support$FRank<-match(support$Gene,ranked)
  genes<-support$Gene[support$Methods>=cfg$consensus_min_methods];genes<-ranked[ranked %in% genes]
  if(length(genes)<cfg$minimum_features)stop("Consensus selection has fewer than minimum features; no fallback: ",length(genes))
  ans<-list(genes=genes,baseline=head(ranked,cfg$baseline_features),support=support,core=support$Gene[support$Methods==3],training_samples=sort(meta$Sample[train]),lasso_lambda=lasso$lambda.min,counts=c(Limma=length(lim),RF=length(rf_genes),LASSO=length(las),Consensus=length(genes),Core=length(which(support$Methods==3))),seconds=unname(proc.time()[3]-t0))
  atomic_save(ans,path);cat("selection",length(train),"samples",paste(names(ans$counts),ans$counts,collapse=";"),sprintf("%.1fs",ans$seconds),"\n");flush.console();ans
}
ba<-function(y,p)mean(vapply(classes,function(c)mean(p[y==c]==c),numeric(1)))
align_prob<-function(p){a<-matrix(0,nrow(p),3,dimnames=list(rownames(p),classes));a[,colnames(p)]<-p;stopifnot(all(is.finite(a)),all(rowSums(a)>0));a/rowSums(a)}
fit_model<-function(kind,tr,y,param,seed){set.seed(seed);if(kind=="SVM")return(svm(x=tr,y=factor(y,levels=classes),kernel="linear",cost=param,scale=FALSE,probability=TRUE));if(kind=="RF")return(randomForest(x=tr,y=factor(y,levels=classes),ntree=cfg$rf_trees,mtry=max(1,min(ncol(tr),round(sqrt(ncol(tr))*param)))));list(train=tr,y=factor(y,levels=classes),k=as.integer(param))}
predict_prob<-function(kind,model,te){if(kind=="SVM"){p<-predict(model,te,probability=TRUE);return(align_prob(attr(p,"probabilities")))};if(kind=="RF")return(align_prob(predict(model,te,type="prob")));dist<-outer(rowSums(te^2),rowSums(model$train^2),"+")-2*tcrossprod(te,model$train);p<-t(vapply(seq_len(nrow(te)),function(i){ix<-head(order(dist[i,],seq_len(ncol(dist))),model$k);as.numeric(table(factor(model$y[ix],levels=classes)))/length(ix)},numeric(3)));colnames(p)<-classes;rownames(p)<-rownames(te);p}
specs<-list(Consensus_SVM=list(kind="SVM",selection="genes",grid=cfg$svm_cost),Limma123_SVM=list(kind="SVM",selection="baseline",grid=cfg$svm_cost),Consensus_RF=list(kind="RF",selection="genes",grid=cfg$rf_mtry_multiplier),Consensus_kNN=list(kind="kNN",selection="genes",grid=cfg$knn_k))
tune<-function(train,inner,models,tag){
  allpred<-lapply(models,function(n)matrix(NA_character_,length(train),length(specs[[n]]$grid)));names(allpred)<-models
  todo<-sort(unique(inner));splits<-lapply(todo,function(k){list(train=train[inner!=k],test=train[inner==k],local=which(inner==k))})
  sels<-parallel::mclapply(splits,function(s)select_features(s$train),mc.cores=cfg$workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
  if(any(vapply(sels,inherits,logical(1),"try-error")))stop("Selection worker failed")
  audit<-list();ti<-list()
  for(j in seq_along(splits)){s<-splits[[j]];stopifnot(!length(intersect(s$train,s$test)));z<-correct_split(s$train,s$test)
    audit[[j]]<-data.frame(Context=tag,InnerFold=todo[j],Sample=meta$Sample[train],Role=ifelse(inner==todo[j],"Inner validation","Inner training"))
    for(n in models){sp<-specs[[n]];genes<-sels[[j]][[sp$selection]];scaled<-scale_split(z$train[,genes,drop=FALSE],z$test[,genes,drop=FALSE]);for(g in seq_along(sp$grid)){m<-fit_model(sp$kind,scaled$train,meta$Subtype[s$train],sp$grid[g],seed_for(paste(tag,j,n,g)));p<-predict_prob(sp$kind,m,scaled$test);allpred[[n]][s$local,g]<-classes[max.col(p,ties.method="first")]}}
  }
  chosen<-list();rows<-list()
  for(n in models){stopifnot(!anyNA(allpred[[n]]));scores<-apply(allpred[[n]],2,function(p)ba(meta$Subtype[train],p));chosen[[n]]<-specs[[n]]$grid[which.max(scores)];rows[[n]]<-data.frame(Context=tag,Model=n,Parameter=specs[[n]]$grid,InnerBalancedAccuracy=scores,Selected=seq_along(scores)==which.max(scores))}
  list(chosen=chosen,tuning=do.call(rbind,rows),membership=do.call(rbind,audit))
}
fit_split<-function(train,test,inner,models,tag,analysis,repeat_id=NA_integer_,fold=NA_integer_){
  stopifnot(!length(intersect(train,test)));path<-file.path(run,"checkpoints",paste0(tag,".rds"));if(file.exists(path))return(readRDS(path))
  status("RUNNING",tag);sel<-select_features(train);tu<-tune(train,inner,models,tag);z<-correct_split(train,test);pred<-list();bundles<-list()
  for(n in models){sp<-specs[[n]];genes<-sel[[sp$selection]];sc<-scale_split(z$train[,genes,drop=FALSE],z$test[,genes,drop=FALSE]);m<-fit_model(sp$kind,sc$train,meta$Subtype[train],tu$chosen[[n]],seed_for(paste(tag,n,"refit")));p<-predict_prob(sp$kind,m,sc$test);pred[[n]]<-data.frame(Analysis=analysis,Context=tag,Repeat=repeat_id,Fold=fold,Model=n,Sample=meta$Sample[test],Study=meta$Study[test],Truth=as.character(meta$Subtype[test]),Prediction=classes[max.col(p,ties.method="first")],ProbabilityC1=p[,1],ProbabilityC2=p[,2],ProbabilityC3=p[,3],Features=length(genes),Parameter=tu$chosen[[n]],UnseenStudyCentering=meta$Study[test]%in%z$unseen)}
  ans<-list(predictions=do.call(rbind,pred),selection=sel,tuning=tu$tuning,membership=tu$membership,train=meta$Sample[train],test=meta$Sample[test]);atomic_save(ans,path);ans
}
provenance<-data.frame(Role=names(inputs),Path=unname(inputs),SHA256=unname(hashes),Bytes=file.info(inputs)$size)
write.csv(provenance,file.path(run,"provenance","input_manifest.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(run,"provenance","sessionInfo.txt"))
write_table(meta,"population")
if(stage=="preflight"){
  required<-c("limma","randomForest","glmnet","e1071","pROC","jsonlite","digest")
  write_table(data.frame(Package=required,Version=vapply(required,function(p)as.character(packageVersion(p)),character(1))),"environment")
  checks<-data.frame(Check=c("canonical_main_497","finite_expression","labels_match_locked_folds","25_outer_splits","one_row_per_sample_study_namespace"),Pass=c(nrow(meta)==497,all(is.finite(X)),all(folds$Subtype==as.character(meta$Subtype[match(folds$Sample,meta$Sample)])),length(unique(paste(folds$Repeat,folds$Fold)))==25,!anyDuplicated(paste(meta$Study,meta$Sample))))
  write_table(checks,"preflight_checks");stopifnot(all(checks$Pass));status("PREFLIGHT_PASS",paste(nrow(X),"samples",ncol(X),"genes"));quit(save="no")
}
if(stage=="smoke"){
  f<-folds[folds$Repeat==1 & folds$Fold==1,];train<-match(f$Sample[f$Role=="Outer train"],meta$Sample);test<-match(f$Sample[f$Role=="Outer test"],meta$Sample)
  s<-select_features(train);stopifnot(all(s$genes %in% colnames(X)),identical(s$training_samples,sort(meta$Sample[train])),all(s$support$Methods[s$support$Gene%in%s$genes]>=2))
  z<-correct_split(train,test);sc<-scale_split(z$train[,s$genes,drop=FALSE],z$test[,s$genes,drop=FALSE]);for(kind in c("SVM","RF","kNN")){m<-fit_model(kind,sc$train,meta$Subtype[train],if(kind=="kNN")5 else 1,42);p<-predict_prob(kind,m,sc$test);stopifnot(identical(colnames(p),classes),max(abs(rowSums(p)-1))<1e-10,all(p>=0),nrow(p)==length(test))}
  original_test<-X[test,,drop=FALSE];X[test,]<-X[test,,drop=FALSE]+100
  changed<-correct_split(train,test);stopifnot(identical(z$train,changed$train),identical(z$offsets,changed$offsets))
  X[test,]<-original_test;stopifnot(!length(intersect(s$training_samples,meta$Sample[test])))
  write_table(data.frame(Check=c("consensus_support_at_least_two","SVM_RF_kNN_probability_contracts","test_expression_does_not_change_training_correction","training_test_membership_disjoint"),Pass=TRUE),"smoke_checks")
  write_table(data.frame(Method=names(s$counts),N=unname(s$counts)),"smoke_feature_counts");status("SMOKE_PASS",sprintf("selection %.1fs; three prediction interfaces pass",s$seconds));quit(save="no")
}
if(stage=="run"){
  for(source_path in c("5-ClassifierValidation/scripts/run_step5_consensus_revision.R","5-ClassifierValidation/scripts/finalize_step5_consensus_revision.R","5-ClassifierValidation/tests/verify_consensus_revision.py",cfg_path)){
    destination<-file.path(run,"provenance",basename(source_path));if(file.exists(destination))stopifnot(identical(digest(file=source_path,algo="sha256"),digest(file=destination,algo="sha256"))) else stopifnot(file.copy(source_path,destination))
  }
  nested<-list();for(rep in 1:cfg$outer_repeats)for(k in 1:cfg$outer_folds){f<-folds[folds$Repeat==rep & folds$Fold==k,];train<-match(f$Sample[f$Role=="Outer train"],meta$Sample);test<-match(f$Sample[f$Role=="Outer test"],meta$Sample);inner<-make_folds(meta$Subtype[train],meta$Study[train],cfg$inner_folds,cfg$seed+rep*100+k);tag<-sprintf("nested_R%d_F%d",rep,k);nested[[tag]]<-fit_split(train,test,inner,names(specs),tag,"Repeated nested CV",rep,k)}
  write_table(do.call(rbind,lapply(nested,`[[`,"predictions")),"nested_predictions_by_repeat");write_table(do.call(rbind,lapply(nested,`[[`,"tuning")),"nested_tuning");write_table(do.call(rbind,lapply(nested,`[[`,"membership")),"inner_membership")
  loso<-list();for(st in sort(unique(meta$Study))){train<-which(meta$Study!=st);test<-which(meta$Study==st);inner<-as.integer(factor(meta$Study[train]));tag<-paste0("loso_",gsub("[^A-Za-z0-9]","_",st));loso[[tag]]<-fit_split(train,test,inner,c("Consensus_SVM","Limma123_SVM"),tag,"LOSO cohort-transductive sensitivity")}
  write_table(do.call(rbind,lapply(loso,`[[`,"predictions")),"loso_predictions");write_table(do.call(rbind,lapply(loso,`[[`,"tuning")),"loso_tuning");write_table(do.call(rbind,lapply(loso,`[[`,"membership")),"loso_inner_membership")
  status("TRAINING_COMPLETE","Nested CV and LOSO complete; final models/projections next")
}
if(stage %in% c("run","finalize")) source("5-ClassifierValidation/scripts/finalize_step5_consensus_revision.R")
