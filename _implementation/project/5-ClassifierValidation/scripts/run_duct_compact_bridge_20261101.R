options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(e1071);library(jsonlite);library(digest)})
run<-'5-ClassifierValidation/review/20261101_duct_compact_bridge_v1';dir.create(file.path(run,'tables'),recursive=TRUE,showWarnings=FALSE)
parent<-'5-ClassifierValidation/runs/bulk-rnaseq__20261030_step5_consensus_v1';compact<-'5-ClassifierValidation/runs/bulk-rnaseq__20261031_step5_minpanel_v1'
inputs<-character();rd<-function(p){inputs<<-unique(c(inputs,p));fread(p,colClasses=c(Sample='character'))};wr<-function(x,n)fwrite(x,file.path(run,'tables',paste0(n,'.csv')))
meta<-as.data.frame(rd(file.path(parent,'tables/population.csv')));classes<-c('C1','C2','C3');meta$Subtype<-factor(meta$Subtype,levels=classes)
master<-rd('governance/data/subtype_master_table.csv')[SubtypeClass=='Main'];stopifnot(nrow(meta)==497,setequal(meta$Sample,master$Sample),all(as.character(meta$Subtype)==master$FinalSubtype[match(meta$Sample,master$Sample)]))
ep<-'1-Clustering/data/combined_expr_before_correction.rds';inputs<-c(inputs,ep);env<-new.env();load(ep,env)
panels<-list(Anchor3=c('S100P','CRP','CDH2'),Anchor4=c('S100P','CRP','CDH2','MUC6'));X<-t(env$combined_expr[unique(unlist(panels)),meta$Sample,drop=FALSE]);storage.mode(X)<-'double';stopifnot(all(is.finite(X)));rm(env)
hp<-file.path(parent,'provenance/run_step5_consensus_revision.R');inputs<-c(inputs,hp)
helpers<-c('seed_for','correct_split','scale_split','ba','align_prob','fit_model','predict_prob')
for(e in parse(hp))if(is.call(e)&&identical(e[[1]],as.name('<-'))&&as.character(e[[2]])[1]%in%helpers)eval(e)
stopifnot(all(vapply(helpers,exists,logical(1))))
member<-rd(file.path(compact,'tables/membership.csv'));oldpred<-rd(file.path(compact,'tables/predictions_by_repeat.csv'));costs<-c(.01,.1,1,10)
pred<-list();tuning<-list();audit<-list();scales<-list()
for(tag in sort(unique(oldpred$Context))){
 message('Duct reference, same outer split: ',tag);cp<-file.path(parent,'checkpoints',paste0(tag,'.rds'));inputs<-c(inputs,cp);ck<-readRDS(cp)
 train<-match(ck$train,meta$Sample);test<-match(ck$test,meta$Sample);stopifnot(!anyNA(c(train,test)),!length(intersect(train,test)),setequal(meta$Sample[test],oldpred[Context==tag]$Sample))
 im<-member[Context==tag&Role=='Inner validation'];inner<-im$InnerFold[match(meta$Sample[train],im$Sample)];stopifnot(!anyNA(inner),!anyDuplicated(im$Sample),length(unique(inner))==4)
 for(model in names(panels)){
  genes<-panels[[model]];ip<-matrix(NA_character_,length(train),length(costs))
  for(j in 1:4){itr<-train[inner!=j];ite<-train[inner==j];z<-correct_split(itr,ite);sc<-scale_split(z$train[,genes,drop=FALSE],z$test[,genes,drop=FALSE])
   stopifnot(!length(intersect(itr,ite)),max(abs(colMeans(sc$train)))<1e-8)
   for(ci in seq_along(costs)){m<-fit_model('SVM',sc$train,meta$Subtype[itr],costs[ci],seed_for(paste(tag,model,j,ci,'anchor')));p<-predict_prob('SVM',m,sc$test);ip[inner==j,ci]<-classes[max.col(p,ties.method='first')]}
   audit[[length(audit)+1]]<-data.table(Context=tag,Model=model,InnerFold=j,TrainN=length(itr),ValidationN=length(ite),TrainingTestOverlap=length(intersect(itr,ite)),OuterTestInInnerTraining=length(intersect(test,itr)))
  }
  stopifnot(!anyNA(ip));b<-apply(ip,2,function(p)ba(meta$Subtype[train],p));cstar<-costs[which.max(b)]
  tuning[[length(tuning)+1]]<-data.table(Context=tag,Model=model,Cost=costs,InnerBA=b,Selected=costs==cstar)
  z<-correct_split(train,test);sc<-scale_split(z$train[,genes,drop=FALSE],z$test[,genes,drop=FALSE]);m<-fit_model('SVM',sc$train,meta$Subtype[train],cstar,seed_for(paste(tag,model,'refit')));p<-predict_prob('SVM',m,sc$test)
  pred[[length(pred)+1]]<-data.table(Context=tag,Model=model,Sample=meta$Sample[test],Study=meta$Study[test],Truth=as.character(meta$Subtype[test]),Prediction=classes[max.col(p,ties.method='first')],ProbabilityC1=p[,1],ProbabilityC2=p[,2],ProbabilityC3=p[,3],Cost=cstar,Features=length(genes))
  scales[[length(scales)+1]]<-data.table(Context=tag,Model=model,Gene=genes,TrainingCenter=sc$center,TrainingSD=sc$scale)
 }
}
pred<-rbindlist(pred);wr(pred,'anchor_predictions_by_repeat');wr(rbindlist(tuning),'anchor_inner_tuning');wr(rbindlist(audit),'split_audit');wr(rbindlist(scales),'outer_training_scaling')
oof<-pred[,{stopifnot(.N==5);p<-colMeans(.SD);.(Study=Study[1],Truth=Truth[1],Prediction=classes[which.max(p)],ProbabilityC1=p[1],ProbabilityC2=p[2],ProbabilityC3=p[3])},by=.(Model,Sample),.SDcols=paste0('Probability',classes)]
cpoof<-rd(file.path(compact,'tables/oof_predictions.csv'));cpoof[,Model:='CompactWorkflow'];oof<-rbind(oof,cpoof,use.names=TRUE)
full<-rd(file.path(parent,'tables/nested_oof_predictions.csv'))[Model=='Consensus_SVM'];full[,Model:='FullWorkflow'];oof<-rbind(oof,full[,names(oof),with=FALSE],use.names=TRUE);setorder(oof,Model,Sample)
stopifnot(all(table(oof$Model)==497),all(abs(rowSums(as.matrix(oof[,paste0('Probability',classes),with=FALSE]))-1)<1e-10))
wr(oof,'all_oof_predictions')
auc<-function(truth,p,cl){a<-truth==cl;n1<-sum(a);n0<-sum(!a);if(min(n1,n0)==0)return(NA_real_);(sum(rank(p)[a])-n1*(n1+1)/2)/(n1*n0)}
metric<-function(d){rec<-vapply(classes,function(cl)mean(d$Prediction[d$Truth==cl]==cl),numeric(1));prec<-vapply(classes,function(cl)if(any(d$Prediction==cl))mean(d$Truth[d$Prediction==cl]==cl)else 0,numeric(1));au<-vapply(classes,function(cl)auc(d$Truth,d[[paste0('Probability',cl)]],cl),numeric(1));Y<-outer(d$Truth,classes,'==')+0;P<-as.matrix(d[,paste0('Probability',classes),with=FALSE]);c(BalancedAccuracy=mean(rec),setNames(rec,paste0(classes,'Recall')),Accuracy=mean(d$Prediction==d$Truth),MacroF1=mean(ifelse(rec+prec>0,2*rec*prec/(rec+prec),0)),setNames(au,paste0(classes,'AUROC')),MacroAUROC=mean(au),MulticlassBrier=mean(rowSums((P-Y)^2)))}
ref<-oof[Model=='Anchor3'];gr<-split(seq_len(nrow(ref)),interaction(ref$Study,ref$Truth,drop=TRUE));set.seed(20261101);draws<-replicate(1000,unlist(lapply(gr,function(ii)ii[sample.int(length(ii),length(ii),replace=TRUE)])),simplify=FALSE)
models<-unique(oof$Model);vals<-list();bs<-list();perf<-list();dif<-list()
for(model in models){d<-oof[Model==model];stopifnot(identical(d$Sample,ref$Sample),identical(d$Truth,ref$Truth));vals[[model]]<-metric(d);bs[[model]]<-t(vapply(draws,function(ii)metric(d[ii]),vals[[model]]));perf[[model]]<-data.table(Model=model,Metric=names(vals[[model]]),Estimate=as.numeric(vals[[model]]),Lower95=apply(bs[[model]],2,quantile,.025),Upper95=apply(bs[[model]],2,quantile,.975))}
for(base in c('Anchor3','Anchor4','FullWorkflow')){delta<-bs$CompactWorkflow-bs[[base]];dif[[base]]<-data.table(Comparison=paste0('CompactWorkflow-',base),Metric=names(vals$CompactWorkflow),Estimate=as.numeric(vals$CompactWorkflow-vals[[base]]),Lower95=apply(delta,2,quantile,.025),Upper95=apply(delta,2,quantile,.975))}
wr(rbindlist(perf),'performance_ci');wr(rbindlist(dif),'paired_differences')
conf<-oof[,as.data.table(table(factor(Truth,levels=classes),factor(Prediction,levels=classes))),by=Model];setnames(conf,c('V1','V2','N'),c('Truth','Prediction','Count'));conf[,`:=`(RowN=sum(Count),RowFraction=Count/sum(Count)),by=.(Model,Truth)];wr(conf,'confusion_matrices')
perstudy<-oof[,{d<-copy(.SD);m<-metric(d);.(Metric=names(m),Value=as.numeric(m),N=.N,C1N=sum(Truth=='C1'),C2N=sum(Truth=='C2'),C3N=sum(Truth=='C3'))},by=.(Model,Study)];wr(perstudy,'study_specific_performance')
paired<-merge(oof[Model=='Anchor3',.(Sample,Study,Truth,AnchorPrediction=Prediction)],oof[Model=='CompactWorkflow',.(Sample,CompactPrediction=Prediction)],by='Sample');paired[,Transition:=fifelse(AnchorPrediction==Truth,fifelse(CompactPrediction==Truth,'Both correct','Lost correct classification'),fifelse(CompactPrediction==Truth,'Corrected classification','Both incorrect'))];wr(paired,'patient_classification_transitions');wr(paired[,.(N=.N),by=.(Truth,Transition)],'patient_transition_counts')
weights<-rd(file.path(compact,'tables/full_candidate_weights.csv'));mp<-file.path(compact,'models/compact_candidate.rds');inputs<-c(inputs,mp);m<-readRDS(mp);stopifnot(setequal(weights$Gene,m$genes),length(m$genes)==17)
wr(data.table(Path=unique(inputs),SHA256=vapply(unique(inputs),function(p)digest(file=p,algo='sha256'),character(1))),'input_hashes')
write_json(list(Status='COMPUTED_REVIEW',Patients=497,OuterSplits=25,InnerSplits=4,PrimaryAnchorGenes=panels$Anchor3,SensitivityAnchorGenes=panels$Anchor4,Fixed17Validation=FALSE,CompactOOF='Existing adaptive training-fold selection; final separate candidate17',Anchor3BA=vals$Anchor3['BalancedAccuracy'],Anchor4BA=vals$Anchor4['BalancedAccuracy'],CompactBA=vals$CompactWorkflow['BalancedAccuracy'],PairedBAGain=vals$CompactWorkflow['BalancedAccuracy']-vals$Anchor3['BalancedAccuracy']),file.path(run,'summary.json'),pretty=TRUE,auto_unbox=TRUE)
print(rbindlist(perf)[Metric%chin%c('BalancedAccuracy','C1Recall','C2Recall','C3Recall')])
