# Read-only extension of the fixed classifier: recover RF evidence and join external outcomes.
suppressPackageStartupMessages({library(data.table);library(randomForest);library(digest);library(jsonlite);library(survival)})
run <- '5-ClassifierValidation/runs/bulk-rnaseq__20261030_step5_consensus_v1'
out <- '8-Visualization/refinements/fig5_optimized_20261031_v5'
tab <- function(n)fread(file.path(run,'tables',paste0(n,'.csv')))
save_table <- function(x,n)fwrite(x,file.path(out,'source_data',paste0(n,'.csv')))
cfg <- fromJSON(file.path(run,'provenance/consensus_revision_20261030.json'))
identity <- readLines(file.path(run,'provenance/identity.txt'))
meta <- tab('population');meta[,Subtype:=factor(Subtype,levels=c('C1','C2','C3'))]
stopifnot(nrow(meta)==497,!anyDuplicated(meta$Sample))
rf_path <- file.path(out,'source_data/RF_importance_all_genes.csv')
if(!file.exists(rf_path)){
 env<-new.env();load('1-Clustering/data/combined_expr_before_correction.rds',env)
 x<-t(as.matrix(env$combined_expr[,meta$Sample,drop=FALSE]));storage.mode(x)<-'double';mu<-colMeans(x)
 offsets<-lapply(split(seq_len(nrow(meta)),meta$Study),function(ii)colMeans(x[ii,,drop=FALSE])-mu)
 for(st in unique(meta$Study)){ii<-which(meta$Study==st);x[ii,]<-sweep(x[ii,,drop=FALSE],2,offsets[[st]],'-')}
 key<-digest(list(identity,sort(meta$Sample)),algo='sha256')
 seed<-strtoi(substr(digest(paste(key,'RF'),algo='xxhash32'),1,7),16L)
 set.seed(seed);cat('Recovering full-development RF importance; exact original seed',seed,'\n');flush.console()
 fit<-randomForest(x=x,y=meta$Subtype,ntree=cfg$rf_trees,importance=TRUE,keep.forest=FALSE)
 imp<-importance(fit,type=1);stopifnot(identical(colnames(imp),'MeanDecreaseAccuracy'))
 d<-data.table(Gene=rownames(imp),RFImportance=as.numeric(imp[,1]));setorder(d,-RFImportance,Gene);d[,RFRank:=seq_len(.N)]
 support<-tab('final_feature_method_support');stopifnot(setequal(d[RFRank<=cfg$rf_top,Gene],support[RF==TRUE,Gene]))
 save_table(d,'RF_importance_all_genes')
 write_json(list(passed=TRUE,seed=seed,selector_key=key,rf_trees=cfg$rf_trees,top200_exact_membership=TRUE,importance_definition='randomForest importance(type=1, scale=TRUE): scaled permutation MeanDecreaseAccuracy',model_refit=FALSE),file.path(out,'review/RF_reproduction.json'),auto_unbox=TRUE,pretty=TRUE)
}
evidence<-fread('8-Visualization/refinements/fig5_marker_story_20261030_v2/source_data/marker_priority.csv')
evidence<-merge(evidence,fread(rf_path),by='Gene');setorder(evidence,RFRank)
means<-fread('8-Visualization/refinements/fig5_marker_logic_20261030_v3/source_data/complete114_evidence.csv')[,c('Gene','C1MeanTrainingZ','C2MeanTrainingZ','C3MeanTrainingZ'),with=FALSE]
evidence<-merge(evidence,means,by='Gene',all.x=TRUE);setorder(evidence,RFRank)
save_table(evidence,'feature_evidence_all_genes')
panel<-evidence[FinalSelected==TRUE];panel[,PanelRFRank:=seq_len(.N)]
stopifnot(nrow(panel)==114,all(panel$Methods>=2));save_table(panel,'feature_evidence_114')

# Use current frozen projection labels. Outcome fields never affect labels or feature selection.
pred<-tab('external_sample_predictions')[Model=='Consensus_SVM' & Track=='Microarray']
pred[,Study:=fifelse(grepl('2017Jusakual',Cohort),'2017Jusakual',fifelse(grepl('2020Job',Cohort),'2020Job',NA_character_))]
pred<-pred[!is.na(Study)]
clinical<-fread('0-Data/qc/sample_information.csv',na.strings=c('','NA'))
clinical<-clinical[Study%in%c('2017Jusakual','2020Job'),.(Sample=SampleID,Study,Time=as.numeric(`Overall survival (months)`),Event=as.numeric(`Vital state (1=Dead)`))]
stopifnot(!anyDuplicated(clinical[,.(Study,Sample)]),!anyDuplicated(pred[,.(Study,Sample)]))
d<-merge(pred,clinical,by=c('Study','Sample'),all.x=TRUE)
d[,Eligible:=is.finite(Time)&Time>0&Event%in%c(0,1)]
d[,Exclusion:=fifelse(Eligible,'Included',fifelse(!is.finite(Time),'Missing time',fifelse(Time<=0,'Nonpositive time','Missing/invalid event')))]
save_table(d,'external_survival_population')
counts<-d[,.(Projected=.N,Eligible=sum(Eligible),Excluded=sum(!Eligible),Events=sum(Event[Eligible])),by=Study];save_table(counts,'external_survival_counts')
curves<-list();risks<-list();tests<-list();classcounts<-list()
for(st in unique(d$Study)){
 z<-d[Study==st & Eligible==TRUE];z[,Prediction:=factor(Prediction,levels=c('C1','C2','C3'))]
 fit<-survfit(Surv(Time,Event)~Prediction,data=z,conf.type='log-log')
 ss<-summary(fit,censored=TRUE)
 curve<-data.table(Study=st,Class=sub('Prediction=','',as.character(ss$strata)),Time=ss$time,Survival=ss$surv,Lower95=ss$lower,Upper95=ss$upper,Risk=ss$n.risk,Events=ss$n.event,Censored=ss$n.censor)
 curve<-rbind(data.table(Study=st,Class=c('C1','C2','C3'),Time=0,Survival=1,Lower95=1,Upper95=1,Risk=as.numeric(table(z$Prediction)),Events=0,Censored=0),curve)
 setorder(curve,Class,Time);curves[[st]]<-curve
 tick_step<-if(max(z$Time)>100)24 else 12
 ticks<-seq(0,ceiling(max(z$Time)/tick_step)*tick_step,by=tick_step)
 risks[[st]]<-rbindlist(lapply(c('C1','C2','C3'),function(cl)data.table(Study=st,Class=cl,Time=ticks,AtRisk=vapply(ticks,function(t)sum(z$Time[z$Prediction==cl]>=t),integer(1)))))
 lr<-survdiff(Surv(Time,Event)~Prediction,data=z);tests[[st]]<-data.table(Study=st,N=nrow(z),Events=sum(z$Event),ChiSquare=lr$chisq,DF=length(lr$n)-1,P=pchisq(lr$chisq,length(lr$n)-1,lower.tail=FALSE))
 classcounts[[st]]<-z[,.(N=.N,Events=sum(Event)),by=.(Study,Prediction)]
}
tt<-rbindlist(tests);tt[,P_Holm:=p.adjust(P,method='holm')]
save_table(tt,'external_survival_tests');save_table(rbindlist(curves),'external_survival_curves');save_table(rbindlist(risks),'external_survival_at_risk');save_table(rbindlist(classcounts),'external_survival_class_counts')
write_json(list(passed=TRUE,projected_records=nrow(d),eligible=sum(d$Eligible),excluded=sum(!d$Eligible),label_source='Frozen Consensus_SVM external_sample_predictions',outcomes_used_for_selection=FALSE),file.path(out,'review/survival_join.json'),auto_unbox=TRUE,pretty=TRUE)
print(counts);print(tt)
