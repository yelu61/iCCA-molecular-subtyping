# Sourced by run_step5_consensus_revision.R after candidate evaluations finish.
status("FINALIZING","Summarizing cross-fitted predictions")
pred<-read.csv(file.path(run,"tables","nested_predictions_by_repeat.csv"),colClasses=c(Sample="character"))
loso_pred<-read.csv(file.path(run,"tables","loso_predictions.csv"),colClasses=c(Sample="character"))
aggregate_predictions<-function(d){do.call(rbind,lapply(split(d,paste(d$Model,d$Sample,sep="::")),function(z){stopifnot(nrow(z)==5,length(unique(z$Repeat))==5);p<-colMeans(z[,paste0("Probability",classes)]);data.frame(Model=z$Model[1],Sample=z$Sample[1],Study=z$Study[1],Truth=z$Truth[1],Prediction=classes[which.max(p)],ProbabilityC1=p[1],ProbabilityC2=p[2],ProbabilityC3=p[3])}))}
oof<-aggregate_predictions(pred);rownames(oof)<-NULL;write_table(oof,"nested_oof_predictions")
metrics<-function(d){
  y<-d$Truth;p<-d$Prediction;prob<-as.matrix(d[,paste0("Probability",classes)])
  recall<-vapply(classes,function(c)mean(p[y==c]==c),numeric(1));precision<-vapply(classes,function(c)if(sum(p==c))sum(p==c & y==c)/sum(p==c) else 0,numeric(1));f1<-ifelse(recall+precision>0,2*recall*precision/(recall+precision),0)
  onehot<-outer(y,classes,`==`)*1
  out<-c(Accuracy=mean(y==p),BalancedAccuracy=mean(recall),MacroF1=mean(f1),MulticlassBrier=mean(rowSums((prob-onehot)^2)))
  for(j in seq_along(classes)){obs<-as.integer(y==classes[j]);score<-prob[,j];ord<-order(-score,seq_along(score));ends<-c(which(diff(score[ord])!=0),length(score));tp<-cumsum(obs[ord])[ends];ap<-if(sum(obs))sum(diff(c(0,tp))/sum(obs)*(tp/ends)) else NA_real_;auc<-if(length(unique(obs))==2)as.numeric(pROC::auc(pROC::roc(obs,score,levels=c(0,1),direction="<",quiet=TRUE))) else NA_real_;out<-c(out,setNames(c(recall[j],precision[j],auc,ap),paste0(classes[j],c("_Sensitivity","_Precision","_AUROC","_AUPRC"))))};out
}
summarize_models<-function(d,label,bootstrap=TRUE){
  models<-sort(unique(d$Model));base<-d[d$Model==models[1],];base<-base[order(base$Sample),];n<-nrow(base)
  ds<-lapply(models,function(nm){z<-d[d$Model==nm,];z<-z[match(base$Sample,z$Sample),];stopifnot(!anyNA(z$Sample),identical(z$Truth,base$Truth));z});names(ds)<-models
  estimates<-lapply(ds,metrics);boots<-NULL
  if(bootstrap){set.seed(cfg$seed);groups<-split(seq_len(n),interaction(base$Study,base$Truth,drop=TRUE));draws<-replicate(cfg$bootstrap_reps,unlist(lapply(groups,function(ii)ii[sample.int(length(ii),length(ii),replace=TRUE)])),simplify=FALSE);boots<-lapply(ds,function(z)t(vapply(draws,function(ix)metrics(z[ix,]),estimates[[1]])))}
  rows<-lapply(models,function(nm){v<-estimates[[nm]];lo<-hi<-rep(NA_real_,length(v));if(bootstrap){lo<-apply(boots[[nm]],2,quantile,.025,na.rm=TRUE);hi<-apply(boots[[nm]],2,quantile,.975,na.rm=TRUE)};data.frame(Analysis=label,Model=nm,N=n,Metric=names(v),Estimate=as.numeric(v),Lower95=lo,Upper95=hi)})
  deltas<-list();if(bootstrap && cfg$primary_model%in%models)for(nm in setdiff(models,cfg$primary_model)){delta<-boots[[cfg$primary_model]]-boots[[nm]];v<-estimates[[cfg$primary_model]]-estimates[[nm]];deltas[[nm]]<-data.frame(Analysis=label,Primary=cfg$primary_model,Comparator=nm,Metric=names(v),Difference=as.numeric(v),Lower95=apply(delta,2,quantile,.025,na.rm=TRUE),Upper95=apply(delta,2,quantile,.975,na.rm=TRUE))}
  list(performance=do.call(rbind,rows),paired=if(length(deltas))do.call(rbind,deltas) else NULL)
}
perf<-summarize_models(oof,"Repeated nested CV");write_table(perf$performance,"nested_performance_ci");write_table(perf$paired,"nested_paired_model_differences")
lp<-summarize_models(loso_pred,"LOSO cohort-transductive sensitivity");write_table(lp$performance,"loso_performance_ci");write_table(lp$paired,"loso_paired_model_differences")
repeat_scores<-do.call(rbind,lapply(split(pred,paste(pred$Model,pred$Repeat)),function(z){v<-metrics(z);data.frame(Model=z$Model[1],Repeat=z$Repeat[1],Metric=names(v),Estimate=as.numeric(v))}));write_table(repeat_scores,"per_repeat_performance")
study_scores<-do.call(rbind,lapply(split(loso_pred,paste(loso_pred$Model,loso_pred$Study)),function(z){v<-metrics(z);data.frame(Model=z$Model[1],Study=z$Study[1],N=nrow(z),Metric=names(v),Estimate=as.numeric(v))}));write_table(study_scores,"loso_per_study_performance")
conf<-do.call(rbind,lapply(split(oof,oof$Model),function(z){d<-as.data.frame(table(Truth=factor(z$Truth,levels=classes),Prediction=factor(z$Prediction,levels=classes)));d$Model<-z$Model[1];d$RowFraction<-d$Freq/as.numeric(table(factor(z$Truth,levels=classes)))[match(d$Truth,classes)];d}));write_table(conf,"confusion_matrix")
cal<-list();bins<-list();curves<-list();ci<-0
for(nm in unique(oof$Model))for(cl in classes){ci<-ci+1;z<-oof[oof$Model==nm,];obs<-as.integer(z$Truth==cl);p<-z[[paste0("Probability",cl)]];clipped<-pmin(1-1e-6,pmax(1e-6,p));logit<-qlogis(clipped);warnings<-character();fits<-withCallingHandlers(list(intercept=glm(obs~1,offset=logit,family=binomial()),joint=glm(obs~logit,family=binomial())),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")});cal[[ci]]<-data.frame(Model=nm,Class=cl,N=nrow(z),Brier=mean((p-obs)^2),CalibrationIntercept=coef(fits$intercept)[1],CalibrationSlope=coef(fits$joint)[2],JointIntercept=coef(fits$joint)[1],Converged=all(vapply(fits,function(m)isTRUE(m$converged),logical(1))),ClippedProbabilities=sum(p!=clipped),Warnings=paste(unique(warnings),collapse=";"))
  bin<-cut(p,breaks=seq(0,1,.1),include.lowest=TRUE,right=TRUE);bins[[ci]]<-do.call(rbind,lapply(levels(bin),function(b){ix<-which(bin==b);data.frame(Model=nm,Class=cl,Bin=b,N=length(ix),Predicted=if(length(ix))mean(p[ix]) else NA_real_,Observed=if(length(ix))mean(obs[ix]) else NA_real_)}))
  roc<-pROC::roc(obs,p,levels=c(0,1),direction="<",quiet=TRUE);curves[[ci]]<-data.frame(Model=nm,Class=cl,Threshold=roc$thresholds,Sensitivity=roc$sensitivities,Specificity=roc$specificities)
}
write_table(do.call(rbind,cal),"calibration_metrics");write_table(do.call(rbind,bins),"calibration_bins");write_table(do.call(rbind,curves),"roc_coordinates")
files<-list.files(file.path(run,"checkpoints"),pattern="^nested_.*rds$",full.names=TRUE);stopifnot(length(files)==25);ck<-lapply(files,readRDS)
sizes<-do.call(rbind,lapply(seq_along(ck),function(i)data.frame(Context=sub(".rds$","",basename(files[i])),Method=names(ck[[i]]$selection$counts),N=unname(ck[[i]]$selection$counts))))
features<-do.call(rbind,lapply(seq_along(ck),function(i){s<-ck[[i]]$selection;d<-s$support;d$Context<-sub(".rds$","",basename(files[i]));d$Consensus<-d$Methods>=cfg$consensus_min_methods;d$Limma123<-d$Gene%in%s$baseline;d}));write_table(sizes,"outer_feature_sizes");write_table(features,"outer_feature_membership")
freq<-aggregate(features[,c("Limma","RF","LASSO","Consensus","Limma123")],list(Gene=features$Gene),sum);for(nm in c("Limma","RF","LASSO","Consensus","Limma123"))freq[[paste0(nm,"Fraction")]]<-freq[[nm]]/25;write_table(freq,"feature_selection_stability")
jac<-do.call(rbind,lapply(seq_along(ck),function(i)do.call(rbind,lapply(seq_len(i-1),function(j){a<-ck[[i]]$selection$genes;b<-ck[[j]]$selection$genes;data.frame(Context1=basename(files[i]),Context2=basename(files[j]),Shared=length(intersect(a,b)),Union=length(union(a,b)),Jaccard=length(intersect(a,b))/length(union(a,b)))}))));write_table(jac,"outer_consensus_jaccard")
status("FINAL_MODEL_FIT","Training full-development candidate models; not replacing locked models")
final_path<-file.path(run,"models","candidate_models.rds")
if(file.exists(final_path)){final_models<-readRDS(final_path)}else{
  all_idx<-seq_len(nrow(meta));sel<-select_features(all_idx);tp<-file.path(run,"checkpoints","final_tuning.rds")
  if(file.exists(tp))tu<-readRDS(tp) else {tu<-tune(all_idx,as.integer(factor(meta$Study)),names(specs),"full_development_inner_study");atomic_save(tu,tp)}
  write_table(tu$tuning,"final_model_tuning");write_table(tu$membership,"final_inner_membership");write_table(sel$support,"final_feature_method_support")
  z<-correct_split(all_idx,all_idx);final_models<-list()
  for(nm in names(specs)){sp<-specs[[nm]];genes<-sel[[sp$selection]];sc<-scale_split(z$train[,genes,drop=FALSE],z$test[,genes,drop=FALSE]);m<-fit_model(sp$kind,sc$train,meta$Subtype,tu$chosen[[nm]],seed_for(paste("final",nm)));final_models[[nm]]<-list(model=m,genes=genes,center=sc$center,scale=sc$scale,study_offsets=z$offsets,global=z$global,classes=classes,kind=sp$kind,parameter=tu$chosen[[nm]],training_samples=meta$Sample,identity=identity)}
  atomic_save(final_models,final_path)
}
write_table(do.call(rbind,lapply(names(final_models),function(nm){m<-final_models[[nm]];data.frame(Model=nm,Rank=seq_along(m$genes),Gene=m$genes,Parameter=m$parameter)})),"final_model_features")
status("EXTERNAL_PROJECTION","Projecting each microarray study/platform independently")
external<-list();ext_qc<-list();ext_index<-0
project_matrix<-function(x,cohort,track="Microarray"){
  storage.mode(x)<-"double";stopifnot(!is.null(rownames(x)),!is.null(colnames(x)),!anyDuplicated(colnames(x)))
  duplicate_rows<-sum(duplicated(rownames(x)));if(duplicate_rows){counts<-table(rownames(x));x<-rowsum(x,rownames(x),reorder=FALSE);x<-x/as.numeric(counts[rownames(x)])}
  for(nm in names(final_models)){
    m<-final_models[[nm]];genes<-m$genes;available<-intersect(genes,rownames(x));stopifnot(length(available)>0)
    sub<-x[available,,drop=FALSE];nonfinite<-sum(!is.finite(sub));sub[!is.finite(sub)]<-NA_real_;mu<-rowMeans(sub,na.rm=TRUE);sdv<-apply(sub,1,sd,na.rm=TRUE);badmu<-!is.finite(mu);badsd<-!is.finite(sdv)|sdv<1e-8;mu[badmu]<-0;sdv[badsd]<-1
    std<-sweep(sweep(sub,1,mu,"-"),1,sdv,"/");std[!is.finite(std)]<-0;inp<-matrix(0,ncol(x),length(genes),dimnames=list(colnames(x),genes));inp[,available]<-t(std)
    p<-predict_prob(m$kind,m$model,inp);ext_index<<-ext_index+1;external[[ext_index]]<<-data.frame(Track=track,Cohort=cohort,Model=nm,Sample=colnames(x),Prediction=classes[max.col(p,ties.method="first")],ProbabilityC1=p[,1],ProbabilityC2=p[,2],ProbabilityC3=p[,3],MaximumProbability=apply(p,1,max),Margin=apply(p,1,function(a)diff(sort(a,decreasing=TRUE)[2:1])),MatchedFeatures=length(available),TotalFeatures=length(genes),MissingRate=1-length(available)/length(genes),ProjectionOnly=TRUE,IndependentTruthAvailable=FALSE)
    ext_qc[[ext_index]]<<-data.frame(Track=track,Cohort=cohort,Model=nm,N=ncol(x),MatchedFeatures=length(available),TotalFeatures=length(genes),MissingGenes=paste(setdiff(genes,available),collapse=";"),DuplicateRowsCollapsed=duplicate_rows,NonfiniteMeasurements=nonfinite,ZeroOrUndefinedSDGenes=sum(badsd),NoFiniteMeasurementsGenes=sum(badmu),Normalization="within-cohort gene z-score (transductive)")
  }
}
for(path in unname(inputs[c("microarray1","microarray2")])){env<-new.env(parent=emptyenv());load(path,env);for(n in sort(ls(env,pattern="^log2dat_"))){x<-get(n,env);if(is.matrix(x)||is.data.frame(x))project_matrix(as.matrix(x),n)};rm(env)}
ext<-do.call(rbind,external);write_table(ext,"external_sample_predictions");write_table(do.call(rbind,ext_qc),"external_coverage_qc")
composition<-do.call(rbind,lapply(split(ext,paste(ext$Model,ext$Cohort)),function(z)data.frame(Model=z$Model[1],Cohort=z$Cohort[1],N=nrow(z),C1=sum(z$Prediction=="C1"),C2=sum(z$Prediction=="C2"),C3=sum(z$Prediction=="C3"),MatchedFeatures=z$MatchedFeatures[1],TotalFeatures=z$TotalFeatures[1],MeanMaximumProbability=mean(z$MaximumProbability),MedianMargin=median(z$Margin))))
write_table(composition,"external_cohort_summary")
oldproj<-read.csv("5-ClassifierValidation/publication/tables/step5_external_projection_sample_predictions.csv",colClasses=c(Sample="character"));oldproj<-oldproj[oldproj$Platform!="DepMap/CCLE",]
oldkey<-paste(oldproj$Platform,oldproj$Sample,sep="::");newkey<-paste(ext$Cohort[ext$Model==cfg$primary_model],ext$Sample[ext$Model==cfg$primary_model],sep="::")
checks<-data.frame(Check=c("four_models_497_oof_each","five_predictions_per_sample_model","25_completed_outer_splits","7_studies_loso_both_models","external_sample_inventory_matches_locked","input_hashes_unchanged","normalized_probabilities","primary_final_consensus_support"),Pass=c(all(table(oof$Model)==497),all(table(paste(pred$Model,pred$Sample))==5),length(files)==25,all(table(loso_pred$Model)==497)&&length(unique(loso_pred$Study))==7,setequal(oldkey,newkey)&&!anyDuplicated(newkey),identical(hashes,vapply(inputs,function(p)digest(file=p,algo="sha256"),character(1))),all(abs(rowSums(oof[,paste0("Probability",classes)])-1)<1e-10),all(read.csv(file.path(run,"tables","final_feature_method_support.csv"))$Methods[match(final_models[[cfg$primary_model]]$genes,read.csv(file.path(run,"tables","final_feature_method_support.csv"))$Gene)]>=2)))
write_table(checks,"technical_checks");stopifnot(all(checks$Pass))
primary<-perf$performance[perf$performance$Model==cfg$primary_model & perf$performance$Metric %in% c("Accuracy","BalancedAccuracy","MacroF1","MulticlassBrier"),]
summary<-c("# Step5 consensus candidate results","","Status: REVIEW. New candidate analysis; locked publication models and figures are unchanged.","",paste("Final primary consensus features:",length(final_models[[cfg$primary_model]]$genes)),"",apply(primary,1,function(r)sprintf("- %s: %.4f (95%% conditional patient-bootstrap interval %.4f–%.4f)",r[["Metric"]],as.numeric(r[["Estimate"]]),as.numeric(r[["Lower95"]]),as.numeric(r[["Upper95"]]))),"",paste("Outer consensus size range:",paste(range(sizes$N[sizes$Method=="Consensus"]),collapse="–")),paste("External microarray samples:",length(newkey),"across",length(unique(ext$Cohort)),"study/platform objects."),"","Interpretation: same-resource label reproducibility; no independent biological-truth/diagnostic validation. Patient-bootstrap intervals condition on predictions, not full pipeline retraining. Calibration results require review. LOSO and external cohort centering are transductive. Comparators do not replace the prespecified primary based on outer/external ranking. Gene/probe mapping and cohort preprocessing are inherited from protected platform inputs. DepMap/PRISM and downstream clinical inference are not refitted here.","","See DESIGN.md and tables/ for complete methods, OOF predictions, paired comparisons, feature stability, calibration, source checks and projection inventories.")
writeLines(summary,file.path(run,"Results.md"))
verification<-system2("python3","5-ClassifierValidation/tests/verify_consensus_revision.py",stdout=TRUE,stderr=TRUE)
writeLines(verification,file.path(run,"logs","independent_verification.log"))
if(!is.null(attr(verification,"status")) && attr(verification,"status")!=0)stop("Independent artifact verification failed; see logs/independent_verification.log")
status("COMPUTED_REVIEW",paste("All scoped technical checks pass; final candidate features",length(final_models[[cfg$primary_model]]$genes)))
