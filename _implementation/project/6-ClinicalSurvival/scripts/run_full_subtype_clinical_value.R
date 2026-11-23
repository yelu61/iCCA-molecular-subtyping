# Full eligible-cohort extension; structurally adapted from the verified observed-duct runner.
# See run-local PROTOCOL.md; the source runner is preserved unchanged.
options(stringsAsFactors=FALSE, warn=1)
suppressPackageStartupMessages({library(survival);library(data.table);library(jsonlite)})
root<-normalizePath(getwd(),mustWork=TRUE)
out<-file.path(root,'6-ClinicalSurvival/review/20261030_full_subtype_clinical_value_v1')
for(d in c('tables','models','review','figures'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
paths<-c('governance/data/subtype_master_table.csv',
 '3-DuctType/publication/tables/step3_all_samples_observed_and_predicted_duct.csv',
 '6-ClinicalSurvival/publication/tables/step6_incremental_prognostic_oof_predictions.csv')
protected<-unique(c(paths,list.files('6-ClinicalSurvival/publication/tables',full.names=TRUE),
 Sys.glob('8-Visualization/figures/main/Figure6*.pdf')))
before<-data.table(Path=protected,MD5=unname(tools::md5sum(protected)))
baseline<-file.path(out,'review/protected_before.csv')
if(!file.exists(baseline))fwrite(before,baseline)
stopifnot(file.exists(file.path(out,'PROTOCOL.md')))
master<-fread(paths[1]);duct<-fread(paths[2]);old<-fread(paths[3])[Model=='Clinical']
stopifnot(!anyDuplicated(master$Sample),!anyDuplicated(duct$Sample),!anyDuplicated(old$Sample),nrow(old)==331)
num<-function(x)suppressWarnings(as.numeric(x))
dat<-copy(master[SubtypeClass=='Main'])
dat<-merge(dat,duct[,.(Sample,DuctObserved)],by='Sample',all.x=TRUE,sort=FALSE)
dat[,`:=`(Age=num(Age_num),Male=fifelse(Sex=='Male',1,fifelse(Sex=='Female',0,NA_real_)),
 StageAdv=fifelse(Stage %chin% c('III','IV'),1,fifelse(Stage %chin% c('I','II'),0,NA_real_)),
 Time=num(`Overall survival (months)`),Event=num(`Vital state (1=Dead)`),
 Large=fifelse(DuctObserved=='large duct type',1,fifelse(DuctObserved=='small duct type',0,NA_real_)),
 C3=as.numeric(FinalSubtype=='C3'), C2=as.numeric(FinalSubtype=='C2'))]
dat[,`:=`(ValidOS=is.finite(Time)&Time>0&Event %in% c(0,1),
 CompleteClinical=complete.cases(Age,Male,StageAdv,Study),ObservedDuct=!is.na(Large))]
fwrite(dat[,.(Sample,Study,FinalSubtype,ValidOS,CompleteClinical,ObservedDuct)],file.path(out,'tables/population_audit.csv'))
z<-dat[ValidOS & CompleteClinical];setorder(z,Study,Sample)
stopifnot(nrow(z)==331,sum(z$Event)==153,uniqueN(z$Study)==4,
 setequal(z$Sample,old$Sample))
check<-merge(z[,.(Sample,Time,Event,Study)],old[,.(Sample,Time,Event,Study)],by='Sample')
stopifnot(max(abs(check$Time.x-check$Time.y))<1e-7,all(check$Event.x==check$Event.y),all(check$Study.x==check$Study.y))
z[,Study:=factor(Study)];z[,Event:=as.integer(Event)]
fwrite(z[,.(Sample,Study,FinalSubtype,DuctObserved,Age,Male,StageAdv,Time,Event,Large,C3,C2)],file.path(out,'tables/analysis_population.csv'))
fwrite(z[,.(N=.N,Events=sum(Event)),by=.(Study,DuctObserved,FinalSubtype)],file.path(out,'tables/patient_event_counts.csv'))
follow<-rbindlist(lapply(c(12,24),function(t)z[,.(N=.N,Cases=sum(Event==1&Time<=t),Controls=sum(Time>t),
 CensoredBefore=sum(Event==0&Time<=t)),by=Study][,Horizon:=t]))
fwrite(follow,file.path(out,'tables/horizon_support.csv'))

specs<-list(M0=c('Age','Male','StageAdv'),MC3=c('Age','Male','StageAdv','C3'),
 MAll=c('Age','Male','StageAdv','C2','C3'))
families<-c('RidgePH','RidgeTimeVarying')
warning_log<-list()
capture<-function(expr,context)withCallingHandlers(expr,warning=function(w){
 warning_log[[length(warning_log)+1L]]<<-data.table(Context=context,Warning=conditionMessage(w));invokeRestart('muffleWarning')})
make_folds<-function(dd,seed){set.seed(seed);ans<-integer(nrow(dd));groups<-interaction(dd$Study,dd$Event,drop=TRUE)
 for(g in levels(groups)){ii<-which(groups==g);ans[ii]<-sample(rep(1:5,length.out=length(ii)))};ans}
folds<-rbindlist(lapply(1:5,function(r)data.table(Repeat=r,Sample=z$Sample,Study=z$Study,Event=z$Event,Fold=make_folds(z,20261030+r))))
fwrite(folds,file.path(out,'tables/fold_definitions.csv'))

standardize<-function(train,new,vars){
 mu<-vapply(train[,..vars],mean,numeric(1));ss<-vapply(train[,..vars],sd,numeric(1))
 stopifnot(all(is.finite(ss)),all(ss>1e-8))
 for(v in vars){train[,(paste0('Z',v)) := (get(v)-mu[v])/ss[v]];new[,(paste0('Z',v)) := (get(v)-mu[v])/ss[v]]}
 list(train=train,new=new,center=mu,scale=ss)
}
expand_training<-function(train,interaction_vars){
 # Match survival's ordinary-Cox near-tie convention before making intervals.
 # This prevents effectively zero intervals from floating-point month values.
 train<-copy(train)
 fixed_time<-survival::aeqSurv(Surv(train$Time,train$Event))
 train[,Time:=unclass(fixed_time)[,1L]]
 event_times<-sort(unique(train[Event==1,Time]));pieces<-lapply(seq_len(nrow(train)),function(i){
  r<-train[i];ends<-sort(unique(c(event_times[event_times<=r$Time],r$Time)))
  a<-r[rep(1L,length(ends))];a[,`:=`(Start=c(0,head(ends,-1)),Stop=ends,Event=as.integer(ends==r$Time & r$Event==1))];a})
 long<-rbindlist(pieces)
 for(v in interaction_vars)long[,(paste0(v,'Time')):=get(v)*log1p(Stop/12)]
 long
}
fit_model<-function(train,new,model,family,context){
 vars<-specs[[model]];st<-standardize(copy(train),copy(new),vars);tr<-st$train;te<-st$new;v<-paste0('Z',vars)
 tv<-family=='RidgeTimeVarying';interaction_vars<-if(tv)intersect(v,c('ZAge','ZC2','ZC3')) else character()
 data<-if(tv)expand_training(tr,interaction_vars) else tr
 time_vars<-if(length(interaction_vars))paste0(interaction_vars,'Time') else character()
 design_vars<-c(v,time_vars)
 stopifnot(!any(design_vars %chin% c('Time','Event','Start','Stop','Sample')))
 rhs<-if(family=='OrdinaryPH')paste(design_vars,collapse=' + ') else
  paste0('ridge(',paste(design_vars,collapse=','),', theta=1, scale=FALSE)')
 lhs<-if(tv)'Surv(Start, Stop, Event)' else 'Surv(Time, Event)'
 form<-as.formula(paste(lhs,'~',rhs,'+ strata(Study)'))
 fit<-capture(coxph(form,data=data,ties='efron',x=TRUE,y=TRUE,model=TRUE,
  control=coxph.control(iter.max=100,eps=1e-9)),context)
 beta<-coef(fit);stopifnot(length(beta)==length(design_vars),all(is.finite(beta)))
 names(beta)<-design_vars
 bh<-capture(basehaz(fit,centered=FALSE),context)
 predict_risk<-function(newdata,times){
  ans<-matrix(NA_real_,nrow(newdata),length(times),dimnames=list(NULL,paste0('Risk',times)))
  for(i in seq_len(nrow(newdata))){
   b<-bh[sub('^Study=','',as.character(strata))==as.character(newdata$Study[i])]
   stopifnot(nrow(b)>0)
   dh<-diff(c(0,b$hazard));stopifnot(min(dh)>-1e-9);dh<-pmax(dh,0)
   linear<-sum(as.numeric(newdata[i,..v])*beta[v]);lp<-rep(linear,nrow(b))
   for(k in interaction_vars)lp<-lp+as.numeric(newdata[[k]][i])*log1p(b$time/12)*beta[paste0(k,'Time')]
   increments<-dh*exp(lp);stopifnot(all(is.finite(increments)))
   for(j in seq_along(times))ans[i,j]<- -expm1(-sum(increments[b$time<=times[j]]))
  };ans
 }
 bh<-as.data.table(bh)
 risk<-predict_risk(te,c(12,24));stopifnot(all(is.finite(risk)),all(risk>=0 & risk<=1))
 list(fit=fit,risk=risk,beta=beta,bh=bh,new=te,train=tr,center=st$center,scale=st$scale,interaction_vars=interaction_vars)
}

predictions<-list();full_fits<-list();parity<-list();diag<-list()
for(family in families){
 message('Fitting ',family)
 for(model in names(specs)){
  full<-fit_model(z,z,model,family,paste(family,model,'full'))
  full_fits[[paste(family,model)]]<-full[c('fit','center','scale','beta','interaction_vars')]
  if(family!='RidgeTimeVarying'){
   # Independent native survfit check of manual baseline-integration prediction.
   for(st in levels(z$Study)){
    ii<-which(as.character(z$Study)==st)[1];sf<-capture(survfit(full$fit,newdata=full$new[ii],se.fit=FALSE),paste(family,model,'native parity'))
    native<-1-as.numeric(summary(sf,times=c(12,24),extend=TRUE)$surv)
    err<-max(abs(native-full$risk[ii,]));stopifnot(err<1e-7)
    parity[[length(parity)+1L]]<-data.table(Family=family,Model=model,Study=st,MaxError=err)
   }
  }
  if(family=='OrdinaryPH'){
   ph<-capture(cox.zph(full$fit),paste(model,'PH diagnostic'))
   pp<-as.data.table(ph$table,keep.rownames='Term');pp[,Model:=model];diag[[model]]<-pp
  }
  for(r in 1:5)for(f in 1:5){ids<-folds[Repeat==r & Fold==f,Sample]
   train<-z[!Sample %chin% ids];te<-z[Sample %chin% ids]
   stopifnot(!any(train$Sample %chin% te$Sample),nrow(train)+nrow(te)==331,all(levels(z$Study) %chin% as.character(train$Study)))
   fit<-fit_model(train,te,model,family,paste(family,model,r,f,sep='/'))
   predictions[[length(predictions)+1L]]<-cbind(data.table(Family=family,Model=model,Repeat=r,Fold=f),
    te[,.(Sample,Study,Time,Event)],as.data.table(fit$risk))
  }
 }
}
pred<-rbindlist(predictions);setorder(pred,Family,Model,Repeat,Sample)
stopifnot(nrow(pred)==2*3*5*331,all(pred[, .N,by=.(Family,Model,Repeat,Sample)]$N==1))
fwrite(pred,file.path(out,'tables/predictions_by_repeat.csv'))
oof<-pred[,.(Study=Study[1],Time=Time[1],Event=Event[1],Risk12=mean(Risk12),Risk24=mean(Risk24)),by=.(Family,Model,Sample)]
setorder(oof,Family,Model,Sample);fwrite(oof,file.path(out,'tables/oof_predictions.csv'))
fwrite(rbindlist(parity),file.path(out,'review/native_prediction_parity.csv'))
fwrite(rbindlist(diag),file.path(out,'tables/ph_diagnostics.csv'))
saveRDS(full_fits,file.path(out,'models/full_data_fits.rds'))

g_eval<-function(z,times,left=FALSE){
 ans<-numeric(nrow(z));stopifnot(length(times)==nrow(z))
 for(st in unique(as.character(z$Study))){ii<-which(as.character(z$Study)==st);sf<-survfit(Surv(Time,1-Event)~1,data=z[ii])
  t<-times[ii]-if(left)1e-8 else 0;ind<-findInterval(t,sf$time);ans[ii]<-c(1,sf$surv)[ind+1L]}
 ans
}
weights_at<-function(z,t){gcase<-g_eval(z,z$Time,TRUE);gt<-g_eval(z,rep(t,nrow(z)))
 list(case=as.numeric(z$Event==1 & z$Time<=t)/pmax(gcase,1e-4),control=as.numeric(z$Time>t)/pmax(gt,1e-4),
      minG=min(c(gcase[z$Event==1 & z$Time<=t],gt[z$Time>t]))) }
metric_set<-function(z,weights=NULL){
 if(is.null(weights))weights<-lapply(c(12,24),function(t)weights_at(z,t))
 auc<-brier<-numeric(2);within<-NA_real_
 for(j in 1:2){risk<-z[[paste0('Risk',c(12,24)[j])]];w<-weights[[j]];ii<-which(w$case>0);kk<-which(w$control>0)
  pair<-outer(w$case[ii],w$control[kk]);order<-outer(risk[ii],risk[kk],'>')+.5*outer(risk[ii],risk[kk],'==')
  auc[j]<-sum(pair*order)/sum(pair);brier[j]<-mean(w$case*(1-risk)^2+w$control*risk^2)
  if(j==2){same<-outer(as.character(z$Study[ii]),as.character(z$Study[kk]),'==');within<-sum(pair*order*same)/sum(pair*same)}
 }
 c(CIndex12=concordance(Surv(Time,Event)~Risk12,data=z,reverse=TRUE)$concordance,
   CIndex12WithinStudy=concordance(Surv(Time,Event)~Risk12+strata(Study),data=z,reverse=TRUE)$concordance,
   AUC12=auc[1],AUC24=auc[2],AUC24WithinStudy=within,Brier12=brier[1],Brier24=brier[2])
}
contrasts<-data.table(Comparison=c('MC3-M0','MAll-M0','MAll-MC3'),
 First=c('MC3','MAll','MAll'),Reference=c('M0','M0','MC3'))
rows<-list();repeat_rows<-list();cal_rows<-list();weights_rows<-list();boot_n<-1000L
for(family in families){
 message('Evaluating paired OOF predictions: ',family)
 zz<-lapply(names(specs),function(m)copy(oof[Family==family & Model==m]));names(zz)<-names(specs)
 stopifnot(all(vapply(zz,function(a)identical(a$Sample,zz[[1]]$Sample),logical(1))))
 ww<-lapply(c(12,24),function(t)weights_at(zz[[1]],t));pts<-sapply(zz,metric_set,weights=ww)
 stopifnot(all(is.finite(pts)))
 for(j in 1:2)weights_rows[[length(weights_rows)+1L]]<-data.table(Family=family,Horizon=c(12,24)[j],MinCensoringSurvival=ww[[j]]$minG,
   MaxWeight=max(c(ww[[j]]$case,ww[[j]]$control)),FloorUsed=ww[[j]]$minG<1e-4)
 set.seed(20261030);boot<-array(NA_real_,dim=c(nrow(pts),ncol(pts),boot_n))
 groups<-split(seq_len(nrow(zz[[1]])),zz[[1]]$Study)
 for(b in seq_len(boot_n)){
  idx<-unlist(lapply(groups,function(ii)sample(ii,length(ii),replace=TRUE)),use.names=FALSE)
  base<-zz[[1]][idx];w<-lapply(c(12,24),function(t)weights_at(base,t))
  for(m in seq_along(zz))boot[,m,b]<-metric_set(zz[[m]][idx],w)
 }
 saveRDS(boot,file.path(out,'review',paste0('conditional_bootstrap_',family,'.rds')))
 for(m in seq_along(zz))for(k in seq_len(nrow(pts))){
  vv<-boot[k,m,];qs<-quantile(vv,c(.025,.975),na.rm=TRUE)
  rows[[length(rows)+1L]]<-data.table(Family=family,Comparison=names(zz)[m],Metric=rownames(pts)[k],
   Estimate=pts[k,m],Lower95=qs[1],Upper95=qs[2],Lower975=NA_real_,Upper975=NA_real_,BootstrapValid=sum(is.finite(vv)),N=331,Events=153)
 }
 for(j in seq_len(nrow(contrasts)))for(k in seq_len(nrow(pts))){
  a<-match(contrasts$First[j],names(zz));b<-match(contrasts$Reference[j],names(zz));vv<-boot[k,a,]-boot[k,b,]
  qs<-quantile(vv,c(.025,.975,.0125,.9875),na.rm=TRUE)
  rows[[length(rows)+1L]]<-data.table(Family=family,Comparison=contrasts$Comparison[j],Metric=rownames(pts)[k],Estimate=pts[k,a]-pts[k,b],
   Lower95=qs[1],Upper95=qs[2],Lower975=qs[3],Upper975=qs[4],BootstrapValid=sum(is.finite(vv)),N=331,Events=153)
 }
 for(r in 1:5){p<-pred[Family==family & Repeat==r];setorder(p,Model,Sample)
  for(m in names(specs)){rr<-metric_set(p[Model==m]);repeat_rows[[length(repeat_rows)+1L]]<-data.table(Family=family,Repeat=r,Model=m,Metric=names(rr),Estimate=as.numeric(rr))}}
 for(m in names(zz))for(j in 1:2){t<-c(12,24)[j];a<-zz[[m]];risk<-pmin(pmax(a[[paste0('Risk',t)]],1e-6),1-1e-6)
  y<-as.integer(a$Time<=t & a$Event==1);w<-ww[[j]]$case+ww[[j]]$control;keep<-w>0
  fit<-capture(glm(y[keep]~qlogis(risk[keep]),family=quasibinomial(),weights=w[keep]),paste(family,m,t,'calibration'))
  cal_rows[[length(cal_rows)+1L]]<-data.table(Family=family,Model=m,Horizon=t,Intercept=coef(fit)[1],Slope=coef(fit)[2],N=331)
 }
}
perf<-rbindlist(rows);fwrite(perf,file.path(out,'tables/model_performance_and_contrasts.csv'))
fwrite(perf[Comparison %chin% c('MC3-M0','MAll-M0') & Metric=='AUC24'],file.path(out,'tables/primary_comparisons.csv'))
fwrite(rbindlist(repeat_rows),file.path(out,'tables/per_repeat_performance.csv'))
fwrite(rbindlist(cal_rows),file.path(out,'tables/calibration.csv'))
fwrite(rbindlist(weights_rows),file.path(out,'review/ipcw_weight_diagnostics.csv'))
warnings<-if(length(warning_log))rbindlist(warning_log) else data.table(Context=character(),Warning=character())
fwrite(warnings,file.path(out,'review/warnings.csv'))
stopifnot(identical(before$MD5,unname(tools::md5sum(before$Path))))
prior<-fread(baseline);stopifnot(identical(prior$MD5,unname(tools::md5sum(prior$Path))))
writeLines(capture.output(sessionInfo()),file.path(out,'review/sessionInfo.txt'))
write_json(list(Status='COMPUTED_REVIEW',N=331,Events=153,Studies=levels(z$Study),Primary='RidgePH',
 PrimaryEndpoint='AUC24',PrimaryContrasts=c('MC3-M0','MAll-M0'),Repeats=5,Folds=5,ConditionalPairedBootstraps=boot_n,
 SourceMD5=as.list(tools::md5sum(c(paths,'6-ClinicalSurvival/scripts/run_full_subtype_clinical_value.R',file.path(out,'PROTOCOL.md')))),
 Warnings=nrow(warnings),ProtectedInputsUnchanged=TRUE),file.path(out,'review/run_summary.json'),pretty=TRUE,auto_unbox=TRUE)
print(perf[Comparison %chin% c('MC3-M0','MAll-M0') & Metric %chin% c('AUC24','CIndex12','Brier24')])
