# Exploratory association extension; protocol is fixed in the new review bundle.
options(stringsAsFactors=FALSE, warn=1)
suppressPackageStartupMessages({library(data.table);library(survival);library(jsonlite)})
out <- '6-ClinicalSurvival/review/20261101_duct_c3_stratified_v1'
owner <- '6-ClinicalSurvival/review/20261030_c3_vs_observed_duct_v1'
stopifnot(file.exists(file.path(out,'PROTOCOL.md')))
for(d in c('tables','models','review','figures'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
protected <- c(list.files(owner,recursive=TRUE,full.names=TRUE),
 'governance/data/subtype_master_table.csv',
 'governance/data/subtype_label_color_dictionary.csv','0-Data/qc/study_color_dictionary.csv')
protected <- protected[!file.info(protected)$isdir]
before <- data.table(Path=protected,MD5=unname(tools::md5sum(protected)))
fwrite(before,file.path(out,'review/protected_before.csv'))
z <- fread(file.path(owner,'tables/analysis_population.csv'))
master <- fread('governance/data/subtype_master_table.csv')
stopifnot(nrow(z)==79,sum(z$Event)==45,!anyDuplicated(z$Sample),
 all(z$Sample %in% master[SubtypeClass=='Main',Sample]),
 all(z$FinalSubtype==master$FinalSubtype[match(z$Sample,master$Sample)]),
 all(z$Event %in% 0:1),all(z$Time>0),!anyNA(z))
z[,`:=`(Study=factor(Study),Age10=Age/10,
 Duct=fifelse(Large==1,'Large duct','Small duct'),
 Molecular=fifelse(C3==1,'C3','Non-C3'))]
z[,Group:=paste(Duct,Molecular,sep=' / ')]
fwrite(z,file.path(out,'tables/analysis_population.csv'))
counts <- z[,.(N=.N,Events=sum(Event),Censored=sum(1-Event),MaxFollowup=max(Time)),by=.(Duct,Molecular)]
fwrite(counts,file.path(out,'tables/group_counts.csv'))
fwrite(z[,.(N=.N,Events=sum(Event)),by=.(Study,Duct,Molecular)],file.path(out,'tables/study_group_counts.csv'))
stopifnot(identical(sort(counts$N),c(7L,11L,15L,46L)))
support <- rbindlist(lapply(c(12,24),function(t) z[,.(N=.N,DeathsBy=sum(Event==1 & Time<=t),
 CensoredBy=sum(Event==0 & Time<=t),AtRisk=sum(Time>=t)),by=.(Study,Duct,Molecular)][,Month:=t]))
fwrite(support,file.path(out,'tables/horizon_support.csv'))
small <- z[Large==0];stopifnot(nrow(small)==61,sum(small$Event)==33,sum(small$C3)==15)

warning_log <- list()
capture <- function(expr,context)withCallingHandlers(expr,warning=function(w){
 warning_log[[length(warning_log)+1L]] <<- data.table(Context=context,Warning=conditionMessage(w));invokeRestart('muffleWarning')})
ph <- capture(coxph(Surv(Time,Event)~C3+Age10+Male+StageAdv+strata(Study),
 data=small,ties='efron',x=TRUE,y=TRUE,model=TRUE,
 control=coxph.control(iter.max=100,eps=1e-9)),'small ordinary PH')
phdiag <- capture(cox.zph(ph),'small PH diagnostics')
fwrite(as.data.table(phdiag$table,keep.rownames='Term'),file.path(out,'tables/ph_diagnostics.csv'))

# Split each patient's observed history at every event time. The covariate
# functions are evaluated at the interval endpoint, i.e. each event risk set.
event_times <- sort(unique(small[Event==1,Time]))
long <- rbindlist(lapply(seq_len(nrow(small)),function(i){
 r<-small[i];ends<-sort(unique(c(event_times[event_times<=r$Time],r$Time)))
 ans<-r[rep(1L,length(ends))]
 ans[,`:=`(Start=c(0,head(ends,-1)),Stop=ends,Event=as.integer(ends==r$Time & r$Event==1))]
 ans[,`:=`(C3Time=C3*log1p(Stop/12),AgeTime=Age10*log1p(Stop/12))];ans
}))
stopifnot(sum(long$Event)==33,uniqueN(long$Sample)==61,all(long$Stop>long$Start))
tv <- capture(coxph(Surv(Start,Stop,Event)~C3+Age10+Male+StageAdv+C3Time+AgeTime+strata(Study),
 data=long,ties='efron',x=TRUE,y=TRUE,model=TRUE,
 control=coxph.control(iter.max=100,eps=1e-9)),'small time-varying')

# Independent survival::tt implementation verifies coefficient/covariance
# equality for the time-varying parameterization.
tv_native <- capture(coxph(Surv(Time,Event)~C3+Age10+Male+StageAdv+tt(C3)+tt(Age10)+strata(Study),
 data=small,ties='efron',tt=function(x,t,...)x*log1p(t/12),
 control=coxph.control(iter.max=100,eps=1e-9)),'small native tt parity')
stopifnot(max(abs(coef(tv)-coef(tv_native)))<1e-7,
 max(abs(vcov(tv)-vcov(tv_native)))<1e-7)

coef_rows <- function(fit,family){b<-coef(fit);se<-sqrt(diag(vcov(fit)))
 data.table(Family=family,Term=names(b),LogHR=unname(b),SE=unname(se),HR=exp(b),
 Lower95=exp(b-1.96*se),Upper95=exp(b+1.96*se),P=2*pnorm(-abs(b/se)),N=61,Events=33)}
allcoef <- rbindlist(list(coef_rows(ph,'PH'),coef_rows(tv,'TimeVarying')))
fwrite(allcoef,file.path(out,'tables/small_duct_all_coefficients.csv'))
effect <- function(fit,family,month=NA_real_){L<-setNames(rep(0,length(coef(fit))),names(coef(fit)));L['C3']<-1
 if(family=='TimeVarying')L['C3Time']<-log1p(month/12)
 b<-sum(L*coef(fit));se<-sqrt(drop(t(L)%*%vcov(fit)%*%L))
 data.table(Family=family,Month=month,LogHR=b,SE=se,HR=exp(b),Lower95=exp(b-1.96*se),
 Upper95=exp(b+1.96*se),P=2*pnorm(-abs(b/se)),N=61,Events=33,C3N=15,NonC3N=46)}
effects <- rbindlist(list(effect(ph,'PH'),effect(tv,'TimeVarying',12),effect(tv,'TimeVarying',24)))
fwrite(effects,file.path(out,'tables/small_duct_c3_effects.csv'))
checks <- rbindlist(lapply(list(PH=ph,TimeVarying=tv),function(f){
 data.table(FiniteCoefficients=all(is.finite(coef(f))),FiniteCovariance=all(is.finite(vcov(f))),
 PositiveVariance=all(eigen(vcov(f),symmetric=TRUE,only.values=TRUE)$values>0),
 Iterations=f$iter,ConditionNumber=kappa(vcov(f)),Converged=f$iter<100)}),idcol='Family')
stopifnot(all(checks$FiniteCoefficients),all(checks$FiniteCovariance),all(checks$PositiveVariance),all(checks$Converged))
fwrite(checks,file.path(out,'review/model_diagnostics.csv'))

sf <- survfit(Surv(Time,Event)~Group,data=z,conf.type='log-log')
ss <- summary(sf,censored=TRUE)
curves <- data.table(Group=sub('^Group=','',as.character(ss$strata)),Time=ss$time,
 AtRisk=ss$n.risk,Deaths=ss$n.event,Censors=ss$n.censor,Survival=ss$surv,Lower95=ss$lower,Upper95=ss$upper)
curves <- merge(curves,z[,.(Group,Duct,Molecular)][!duplicated(Group)],by='Group')
fwrite(curves,file.path(out,'tables/km_curves.csv'))
max_time <- ceiling(max(z$Time)/24)*24
risk <- rbindlist(lapply(seq(0,max_time,by=24),function(t)z[,.(AtRisk=sum(Time>=t)),by=.(Duct,Molecular)][,Month:=t]))
fwrite(risk,file.path(out,'tables/km_risk_table.csv'))
s24<-summary(sf,times=c(12,24),extend=FALSE)
fwrite(data.table(Group=sub('^Group=','',as.character(s24$strata)),Month=s24$time,AtRisk=s24$n.risk,
 Survival=s24$surv,Lower95=s24$lower,Upper95=s24$upper),file.path(out,'tables/km_at_12_24.csv'))
saveRDS(list(PH=ph,TimeVarying=tv,TimeVaryingNative=tv_native,KM=sf),file.path(out,'models/fits.rds'))
warn <- if(length(warning_log))rbindlist(warning_log) else data.table(Context=character(),Warning=character())
fwrite(warn,file.path(out,'review/warnings.csv'))
stopifnot(identical(before$MD5,unname(tools::md5sum(before$Path))))
verification<-data.table(Check=c('Exact source population','Canonical labels','Four groups exhaustive',
 'Small-duct denominator','Event splitting preserves patients and events','Native tt coefficient parity',
 'Native tt covariance parity','Model fits finite and converged','Protected files unchanged'),
 Passed=TRUE,Detail=c('79 patients; 45 deaths','Main-only unchanged','15/46/11/7',
 '61 patients; 33 deaths; 15 C3; 46 non-C3','61 patients; 33 deaths',
 sprintf('maximum difference %.3g',max(abs(coef(tv)-coef(tv_native)))),
 sprintf('maximum difference %.3g',max(abs(vcov(tv)-vcov(tv_native)))),
 'PH and time-varying models',sprintf('%s files',nrow(before))))
fwrite(verification,file.path(out,'review/verification.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'review/sessionInfo.txt'))
write_json(list(Status='COMPUTED_EXPLORATORY_REVIEW',N=79,Events=45,SmallN=61,SmallEvents=33,
 WarningCount=nrow(warn),ProtectedUnchanged=TRUE,MaxPlotTime=max_time,
 SourceMD5=as.list(tools::md5sum(c(file.path(owner,'tables/analysis_population.csv'),file.path(out,'PROTOCOL.md'),
 '6-ClinicalSurvival/scripts/run_duct_c3_stratified_20261101.R')))),file.path(out,'review/run_summary.json'),pretty=TRUE,auto_unbox=TRUE)
print(counts);print(effects);print(phdiag$table);print(checks)
