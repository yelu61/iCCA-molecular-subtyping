options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(e1071);library(jsonlite);library(digest)})
run <- '5-ClassifierValidation/review/20261101_duct_compact_joint_v1'
parent <- '5-ClassifierValidation/runs/bulk-rnaseq__20261030_step5_consensus_v1'
compact <- '5-ClassifierValidation/runs/bulk-rnaseq__20261031_step5_minpanel_v1'
bridge <- '5-ClassifierValidation/review/20261101_duct_compact_bridge_v1'
rfe <- '5-ClassifierValidation/runs/bulk-rnaseq__20261031_step5_svm_rfe_v1'
for(d in c('tables','checkpoints','provenance','review','logs'))dir.create(file.path(run,d),recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(run,'COMPLETED')))stop('Completed immutable run already exists')
status<-function(s,detail){message(format(Sys.time(),'%H:%M:%S'),' ',s,' ',detail);write_json(list(state=s,detail=detail),file.path(run,'status.json'),pretty=TRUE,auto_unbox=TRUE)}
options(error=function(){status('FAILED',geterrmessage());quit(save='no',status=1)})
inputs<-character();hashes<-character()
track<-function(p){if(!p%in%inputs){inputs<<-c(inputs,p);hashes<<-c(hashes,setNames(digest(file=p,algo='sha256'),p))};p}
rd<-function(p){track(p);cols<-names(fread(p,nrows=0));fread(p,colClasses=if('Sample'%in%cols)c(Sample='character')else NULL)}
rr<-function(p)readRDS(track(p))
wr<-function(x,n)fwrite(x,file.path(run,'tables',paste0(n,'.csv')))
classes<-c('C1','C2','C3');anchors<-c('S100P','CRP','CDH2');costs<-c(.01,.1,1,10)
meta<-as.data.frame(rd(file.path(parent,'tables/population.csv')));meta$Subtype<-factor(meta$Subtype,levels=classes)
master<-rd('governance/data/subtype_master_table.csv')[SubtypeClass=='Main']
stopifnot(nrow(meta)==497,!anyDuplicated(meta$Sample),setequal(meta$Sample,master$Sample),all(as.character(meta$Subtype)==master$FinalSubtype[match(meta$Sample,master$Sample)]))
member<-rd(file.path(compact,'tables/membership.csv'));choice<-rd(file.path(compact,'tables/choice.csv'));features<-rd(file.path(compact,'tables/features.csv'))
paths<-rd(file.path(rfe,'tables/rfe_deletion_paths.csv'));oldpred<-rd(file.path(compact,'tables/predictions_by_repeat.csv'))
identity<-readLines(track(file.path(parent,'provenance/identity.txt')))
helpers<-c('seed_for','correct_split','scale_split','ba','align_prob','fit_model','predict_prob')
hp<-track(file.path(parent,'provenance/run_step5_consensus_revision.R'))
for(e in parse(hp))if(is.call(e)&&identical(e[[1]],as.name('<-'))&&as.character(e[[2]])[1]%in%helpers)eval(e)
stopifnot(all(vapply(helpers,exists,logical(1))))
track(file.path(parent,'models/candidate_models.rds'));track(file.path(compact,'models/compact_candidate.rds'))
track('5-ClassifierValidation/scripts/run_duct_compact_joint_20261101.R');track(file.path(run,'DESIGN.md'))
cache<-function(ii){key<-digest(list(identity,sort(meta$Sample[ii])),algo='sha256');p<-file.path(parent,'cache',paste0('selection_',key,'.rds'));s<-rr(p);stopifnot(identical(s$training_samples,sort(meta$Sample[ii])));d<-as.data.table(s$support)[Methods>=2];setorder(d,-Methods,FRank,Gene);list(genes=d$Gene,path=p,training=s$training_samples)}
get_genes<-function(tag,part,ii,K,C){s<-cache(ii);d<-paths[Context==tag&Partition==part&Cost==C];setorder(d,-NBefore)
 if(nrow(d)==0 && part=='OuterTrain'){
  ck<-rr(file.path(compact,'checkpoints',paste0(tag,'.rds')));stopifnot(setequal(ck$train,meta$Sample[ii]),ck$choice$K==K,ck$choice$Cost==C);g<-ck$features[order(ck$features$Rank)]$Gene;stopifnot(length(g)==K,all(g%in%s$genes));return(list(compact=g,union=unique(c(g,anchors)),cache=s$path))
 }
 if(nrow(d)!=length(s$genes)-2)stop(paste('Path incomplete:',tag,part,C,nrow(d),'expected',length(s$genes)-2));stopifnot(identical(d$NBefore,seq.int(length(s$genes),3L)),!anyDuplicated(d$Dropped),all(d$Dropped%in%s$genes));g<-s$genes
 for(i in which(d$NBefore>K)){stopifnot(length(g)==d$NBefore[i],d$Dropped[i]%in%g);g<-g[g!=d$Dropped[i]]}
 stopifnot(length(g)==K);list(compact=g,union=unique(c(g,anchors)),cache=s$path)
}
plans<-list();audit<-list();gene_audit<-list();counts<-list()
for(tag in sort(choice$Context)){
 ck<-rr(file.path(parent,'checkpoints',paste0(tag,'.rds')));tr<-match(ck$train,meta$Sample);te<-match(ck$test,meta$Sample)
 stopifnot(!anyNA(c(tr,te)),!length(intersect(tr,te)),setequal(meta$Sample[te],oldpred[Context==tag]$Sample))
 ch<-choice[Context==tag];stopifnot(nrow(ch)==1,!ch$Fallback)
 im<-member[Context==tag&Role=='Inner validation'];inner<-im$InnerFold[match(meta$Sample[tr],im$Sample)]
 stopifnot(!anyNA(inner),!anyDuplicated(im$Sample),setequal(inner,1:4))
 inner_sets<-list()
 for(j in 1:4){itr<-tr[inner!=j];ite<-tr[inner==j];g<-get_genes(tag,paste0('Inner',j),itr,ch$K,ch$Cost);inner_sets[[j]]<-c(list(train=itr,test=ite),g)
  audit[[length(audit)+1]]<-data.table(Context=tag,Partition=paste0('Inner',j),TrainN=length(itr),ValidationN=length(ite),TrainValidationOverlap=length(intersect(itr,ite)),OuterTestInTraining=length(intersect(te,itr)),OuterTestInValidation=length(intersect(te,ite)),SelectorCache=g$cache,SelectorMembershipExact=TRUE)
  gene_audit[[length(gene_audit)+1]]<-data.table(Context=tag,Partition=paste0('Inner',j),Gene=g$union,CompactMember=g$union%in%g$compact,AnchorMember=g$union%in%anchors)
  counts[[length(counts)+1]]<-data.table(Context=tag,Partition=paste0('Inner',j),CompactK=ch$K,RFEcost=ch$Cost,AnchorsAlreadyPresent=sum(anchors%in%g$compact),AnchorsAdded=length(setdiff(anchors,g$compact)),UnionK=length(g$union))
 }
 g<-get_genes(tag,'OuterTrain',tr,ch$K,ch$Cost);stopifnot(identical(g$compact,features[Context==tag][order(Rank)]$Gene))
 gene_audit[[length(gene_audit)+1]]<-data.table(Context=tag,Partition='OuterTrain',Gene=g$union,CompactMember=g$union%in%g$compact,AnchorMember=g$union%in%anchors)
 counts[[length(counts)+1]]<-data.table(Context=tag,Partition='OuterTrain',CompactK=ch$K,RFEcost=ch$Cost,AnchorsAlreadyPresent=sum(anchors%in%g$compact),AnchorsAdded=length(setdiff(anchors,g$compact)),UnionK=length(g$union))
 plans[[tag]]<-list(train=tr,test=te,inner=inner,splits=inner_sets,genes=g$union,compact=g$compact,K=ch$K,RFEcost=ch$Cost)
}
gene_audit<-rbindlist(gene_audit);counts<-rbindlist(counts);wr(gene_audit,'selected_features');wr(counts,'feature_counts');wr(rbindlist(audit),'split_audit');wr(member,'membership')
env<-new.env();ep<-track('1-Clustering/data/combined_expr_before_correction.rds');load(ep,env)
universe<-sort(unique(gene_audit$Gene));X<-t(env$combined_expr[universe,meta$Sample,drop=FALSE]);storage.mode(X)<-'double';rm(env)
stopifnot(all(is.finite(X)),all(anchors%in%colnames(X)));wr(data.table(Gene=universe),'evaluated_expression_genes')
all<-list()
for(tag in names(plans)){
 status('RUNNING',tag);p<-plans[[tag]];tr<-p$train;te<-p$test;ip<-array(NA_real_,c(length(tr),3,length(costs)));tlogs<-list();scales<-list()
 for(j in 1:4){s<-p$splits[[j]];z<-correct_split(s$train,s$test);sc<-scale_split(z$train[,s$union,drop=FALSE],z$test[,s$union,drop=FALSE]);stopifnot(max(abs(colMeans(sc$train)))<1e-8)
  for(ci in seq_along(costs)){fit<-fit_model('SVM',sc$train,meta$Subtype[s$train],costs[ci],seed_for(paste(tag,j,ci,'joint-anchor-tuning')));prob<-predict_prob('SVM',fit,sc$test);ip[p$inner==j,,ci]<-prob
   tlogs[[length(tlogs)+1]]<-data.table(Context=tag,InnerFold=j,Cost=costs[ci],Sample=meta$Sample[s$test],Truth=as.character(meta$Subtype[s$test]),Prediction=classes[max.col(prob,ties.method='first')],ProbabilityC1=prob[,1],ProbabilityC2=prob[,2],ProbabilityC3=prob[,3])
  }
  scales[[length(scales)+1]]<-data.table(Context=tag,Partition=paste0('Inner',j),Gene=s$union,TrainingCenter=sc$center,TrainingSD=sc$scale)
 }
 stopifnot(!anyNA(ip));b<-vapply(seq_along(costs),function(ci)ba(meta$Subtype[tr],classes[max.col(ip[,,ci],ties.method='first')]),numeric(1));cs<-costs[which.max(b)]
 tuning<-data.table(Context=tag,Cost=costs,InnerBA=b,Selected=costs==cs,InheritedCompactK=p$K,InheritedRFEcost=p$RFEcost)
 z<-correct_split(tr,te);sc<-scale_split(z$train[,p$genes,drop=FALSE],z$test[,p$genes,drop=FALSE]);fit<-fit_model('SVM',sc$train,meta$Subtype[tr],cs,seed_for(paste(tag,'joint-anchor-refit')));prob<-predict_prob('SVM',fit,sc$test)
 scales[[length(scales)+1]]<-data.table(Context=tag,Partition='OuterTrain',Gene=p$genes,TrainingCenter=sc$center,TrainingSD=sc$scale)
 pred<-data.table(Context=tag,Model='JointWorkflow',Sample=meta$Sample[te],Study=meta$Study[te],Truth=as.character(meta$Subtype[te]),Prediction=classes[max.col(prob,ties.method='first')],ProbabilityC1=prob[,1],ProbabilityC2=prob[,2],ProbabilityC3=prob[,3],Cost=cs,Features=length(p$genes),CompactK=p$K,RFEcost=p$RFEcost,AnchorsAdded=length(setdiff(anchors,p$compact)))
 ans<-list(predictions=pred,tuning=tuning,inner_predictions=rbindlist(tlogs),scaling=rbindlist(scales),train=meta$Sample[tr],test=meta$Sample[te],features=p$genes,fit=fit,offsets=z$offsets,center=sc$center,scale=sc$scale)
 saveRDS(ans,file.path(run,'checkpoints',paste0(tag,'.rds')));all[[tag]]<-ans
}
for(n in c('predictions','tuning','inner_predictions','scaling'))wr(rbindlist(lapply(all,`[[`,n)),if(n=='predictions')'joint_predictions_by_repeat' else n)
raw<-rbindlist(lapply(all,`[[`,'predictions'));raw[,Repeat:=as.integer(sub('nested_R([0-9]+)_.*','\\1',Context))]
joof<-raw[,{stopifnot(.N==5,uniqueN(Repeat)==5);p<-colMeans(.SD);.(Study=Study[1],Truth=Truth[1],Prediction=classes[which.max(p)],ProbabilityC1=p[1],ProbabilityC2=p[2],ProbabilityC3=p[3])},by=.(Model,Sample),.SDcols=paste0('Probability',classes)]
old<-rd(file.path(bridge,'tables/all_oof_predictions.csv'))[Model%chin%c('Anchor3','CompactWorkflow','FullWorkflow')];oof<-rbind(old,joof,use.names=TRUE);setorder(oof,Model,Sample)
stopifnot(all(table(oof$Model)==497),all(abs(rowSums(as.matrix(oof[,paste0('Probability',classes),with=FALSE]))-1)<1e-10));wr(oof,'all_oof_predictions')
auc<-function(truth,p,cl){a<-truth==cl;n1<-sum(a);n0<-sum(!a);if(min(n1,n0)==0)return(NA_real_);(sum(rank(p)[a])-n1*(n1+1)/2)/(n1*n0)}
metric<-function(d){rec<-vapply(classes,function(cl)mean(d$Prediction[d$Truth==cl]==cl),numeric(1));prec<-vapply(classes,function(cl)if(any(d$Prediction==cl))mean(d$Truth[d$Prediction==cl]==cl)else NA_real_,numeric(1));au<-vapply(classes,function(cl)auc(d$Truth,d[[paste0('Probability',cl)]],cl),numeric(1));Y<-outer(d$Truth,classes,'==')+0;P<-as.matrix(d[,paste0('Probability',classes),with=FALSE]);c(BalancedAccuracy=mean(rec),setNames(rec,paste0(classes,'Recall')),C3Sensitivity=unname(rec[3]),C3Specificity=mean(d$Prediction[d$Truth!='C3']!='C3'),C3PPV=unname(prec[3]),Accuracy=mean(d$Prediction==d$Truth),MacroF1=mean(ifelse(rec+prec>0,2*rec*prec/(rec+prec),0)),setNames(au,paste0(classes,'AUROC')),MacroAUROC=mean(au),MulticlassBrier=mean(rowSums((P-Y)^2)))}
ref<-oof[Model=='Anchor3'];groups<-split(seq_len(nrow(ref)),interaction(ref$Study,ref$Truth,drop=TRUE));set.seed(20261101);draws<-replicate(1000,unlist(lapply(groups,function(ii)ii[sample.int(length(ii),length(ii),replace=TRUE)])),simplify=FALSE)
models<-c('Anchor3','CompactWorkflow','JointWorkflow','FullWorkflow');vals<-list();bs<-list();perf<-list();dif<-list()
status('SUMMARIZING','Paired patient bootstrap with fixed OOF predictions')
for(model in models){d<-oof[Model==model];stopifnot(identical(d$Sample,ref$Sample),identical(d$Truth,ref$Truth));vals[[model]]<-metric(d);bs[[model]]<-t(vapply(draws,function(ii)metric(d[ii]),vals[[model]]));perf[[model]]<-data.table(Model=model,Metric=names(vals[[model]]),Estimate=as.numeric(vals[[model]]),Lower95=apply(bs[[model]],2,quantile,.025),Upper95=apply(bs[[model]],2,quantile,.975))}
for(pair in list(c('JointWorkflow','Anchor3'),c('JointWorkflow','CompactWorkflow'),c('JointWorkflow','FullWorkflow'),c('CompactWorkflow','Anchor3'))){a<-pair[1];b<-pair[2];delta<-bs[[a]]-bs[[b]];dif[[paste(a,b)]]<-data.table(Comparison=paste0(a,'-',b),Metric=names(vals[[a]]),Estimate=as.numeric(vals[[a]]-vals[[b]]),Lower95=apply(delta,2,quantile,.025),Upper95=apply(delta,2,quantile,.975))}
wr(rbindlist(perf),'performance_ci');wr(rbindlist(dif),'paired_differences');saveRDS(list(draws=draws,samples=ref$Sample,values=vals,bootstrap=bs),file.path(run,'review/paired_bootstrap.rds'))
conf<-oof[,as.data.table(table(factor(Truth,levels=classes),factor(Prediction,levels=classes))),by=Model];setnames(conf,c('V1','V2','N'),c('Truth','Prediction','Count'));conf[,`:=`(RowN=sum(Count),RowFraction=Count/sum(Count)),by=.(Model,Truth)];wr(conf,'confusion_matrices')
transition<-function(before,after)fifelse(before,fifelse(after,'Both correct','New error'),fifelse(after,'Corrected','Both incorrect'))
trans<-list();binary<-list()
for(pair in list(c('Anchor3','JointWorkflow'),c('CompactWorkflow','JointWorkflow'),c('Anchor3','CompactWorkflow'))){b<-oof[Model==pair[1]];a<-oof[Model==pair[2]];stopifnot(identical(a$Sample,b$Sample));p<-data.table(Comparison=paste0(pair[2],'-',pair[1]),Sample=b$Sample,Study=b$Study,Truth=b$Truth,Before=b$Prediction,After=a$Prediction);p[,Transition:=transition(Before==Truth,After==Truth)];trans[[length(trans)+1]]<-copy(p);p[,`:=`(ReferenceGroup=ifelse(Truth=='C3','C3','Non-C3'),BeforeC3=Before=='C3',AfterC3=After=='C3')];p[,Transition:=transition(BeforeC3==(Truth=='C3'),AfterC3==(Truth=='C3'))];binary[[length(binary)+1]]<-p}
trans<-rbindlist(trans);binary<-rbindlist(binary);wr(trans,'patient_classification_transitions');tc<-trans[,.(N=.N),by=.(Comparison,Truth,Transition)];tc<-merge(CJ(Comparison=unique(trans$Comparison),Truth=classes,Transition=c('Both correct','Corrected','New error','Both incorrect')),tc,all.x=TRUE);tc[is.na(N),N:=0L];wr(tc,'patient_transition_counts');wr(binary,'patient_C3_binary_transitions');bc<-binary[,.(N=.N),by=.(Comparison,ReferenceGroup,Transition)];bc<-merge(CJ(Comparison=unique(binary$Comparison),ReferenceGroup=c('C3','Non-C3'),Transition=c('Both correct','Corrected','New error','Both incorrect')),bc,all.x=TRUE);bc[is.na(N),N:=0L];wr(bc,'C3_binary_transition_counts')
perstudy<-oof[,{m<-metric(.SD);.(Metric=names(m),Estimate=as.numeric(m),N=.N)},by=.(Model,Study)];wr(perstudy,'study_specific_performance')
integrity<-data.table(Path=inputs,BeforeSHA256=unname(hashes),AfterSHA256=vapply(inputs,function(p)digest(file=p,algo='sha256'),character(1)));integrity[,Unchanged:=BeforeSHA256==AfterSHA256];stopifnot(all(integrity$Unchanged));wr(integrity,'input_hashes')
outer<-counts[Partition=='OuterTrain'];summary<-list(Status='COMPUTED_REVIEW_NOT_RELEASED',N=497,OuterSplits=25,InnerFolds=4,RawHeldoutPredictions=nrow(raw),JointSVMFits=425,CompactK=list(min=min(outer$CompactK),median=median(outer$CompactK),max=max(outer$CompactK)),UnionK=list(min=min(outer$UnionK),median=median(outer$UnionK),max=max(outer$UnionK)),AnchorsAdded=as.list(setNames(as.integer(table(outer$AnchorsAdded)),names(table(outer$AnchorsAdded)))),Fixed17Validation=FALSE,NewKSearch=FALSE,ComparatorsUnchanged=TRUE,InputHashesUnchanged=all(integrity$Unchanged),PrimaryModelUnchanged=TRUE)
write_json(summary,file.path(run,'summary.json'),pretty=TRUE,auto_unbox=TRUE)
file.copy('5-ClassifierValidation/scripts/run_duct_compact_joint_20261101.R',file.path(run,'provenance/run_duct_compact_joint_20261101.R'),overwrite=TRUE)
capture.output(sessionInfo(),file=file.path(run,'provenance/sessionInfo.txt'))
status('COMPUTED_REVIEW','Joint OOF comparison complete; independent verification still required')
print(rbindlist(perf)[Metric%chin%c('BalancedAccuracy','C3Sensitivity','C3Specificity','C3PPV')]);print(rbindlist(dif)[Comparison=='JointWorkflow-CompactWorkflow'&Metric%chin%c('BalancedAccuracy','C3Sensitivity','C3Specificity')])
