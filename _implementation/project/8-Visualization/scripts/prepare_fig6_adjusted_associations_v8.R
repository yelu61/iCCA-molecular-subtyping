# Extract the complete covariate display from the unchanged Step6 time-varying model.
# This is a parity-gated reproduction, not new model selection or CV.
suppressPackageStartupMessages({library(data.table); library(survival)})
out <- '8-Visualization/refinements/fig6_ABC_20261031_v8'
master <- fread('governance/data/subtype_master_table.csv')
owner <- '6-ClinicalSurvival/publication/tables'
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
d <- copy(master[SubtypeClass=='Main'])
d[, os_months := num(`Overall survival (months)`)]
d[, os_event := num(`Vital state (1=Dead)`)]
d[!os_event %in% c(0,1), os_event := NA_real_]
d[, age_years := num(Age_num)]
d[, sex_clean := factor(fifelse(Sex %chin% c('Female','Male'),Sex,NA_character_),levels=c('Female','Male'))]
d[, stage_group := factor(fifelse(Stage %chin% c('I','II'),'I/II',fifelse(Stage %chin% c('III','IV'),'III/IV',NA_character_)),levels=c('I/II','III/IV'))]
d[, FinalSubtype := factor(FinalSubtype,levels=c('C1','C2','C3'))]
vars <- c('os_months','os_event','FinalSubtype','age_years','sex_clean','stage_group','Study')
d <- d[!is.na(os_months) & os_months>0 & os_event %in% c(0,1)]
d <- d[complete.cases(d[,..vars])]
pop <- fread('6-ClinicalSurvival/review/20261030_full_subtype_clinical_value_v1/tables/analysis_population.csv')
stopifnot(nrow(d)==331,sum(d$os_event)==153,!anyDuplicated(d$Sample),setequal(d$Sample,pop$Sample))
p <- pop[match(d$Sample,Sample)]
stopifnot(max(abs(d$os_months-p$Time))<1e-10,all(d$os_event==p$Event),all(d$age_years==p$Age),
 all(as.character(d$FinalSubtype)==p$FinalSubtype),all(as.integer(d$sex_clean=='Male')==p$Male),
 all(as.integer(d$stage_group=='III/IV')==p$StageAdv),all(d$Study==p$Study))
d[, `:=`(C2_indicator=as.integer(FinalSubtype=='C2'),C3_indicator=as.integer(FinalSubtype=='C3'))]
warnings_seen <- character()
fit <- withCallingHandlers(coxph(Surv(os_months,os_event) ~ C2_indicator+C3_indicator+age_years+sex_clean+stage_group+
 tt(C2_indicator)+tt(C3_indicator)+tt(age_years)+strata(Study),data=d,ties='efron',x=TRUE,
 tt=function(x,t,...) x*log(pmax(t,0.1))),warning=function(w){warnings_seen<<-c(warnings_seen,conditionMessage(w))})
stopifnot(length(warnings_seen)==0)
b <- coef(fit); V <- vcov(fit)
cv <- function(key,t){
 v <- setNames(rep(0,length(b)),names(b))
 if(key %in% c('C2 vs C1','C3 vs C1','C3 vs C2')){
  if(key %in% c('C2 vs C1','C3 vs C2')) {v['C2_indicator']<-if(key=='C2 vs C1')1 else -1;v['tt(C2_indicator)']<-v['C2_indicator']*log(t)}
  if(key %in% c('C3 vs C1','C3 vs C2')) {v['C3_indicator']<-1;v['tt(C3_indicator)']<-log(t)}
 } else if(key=='Age (per 10 years)') {v['age_years']<-10;v['tt(age_years)']<-10*log(t)}
 else if(key=='Male vs female') v['sex_cleanMale']<-1
 else if(key=='Stage III/IV vs I/II') v['stage_groupIII/IV']<-1
 else stop('Unknown contrast')
 eta<-sum(v*b);se<-sqrt(as.numeric(t(v)%*%V%*%v))
 data.table(Contrast=key,TimeMonths=t,HR=exp(eta),Lower95=exp(eta-1.96*se),Upper95=exp(eta+1.96*se),
 PValue=2*pnorm(-abs(eta/se)),N=nrow(d),Events=sum(d$os_event))
}
all <- rbindlist(lapply(c('Age (per 10 years)','Male vs female','Stage III/IV vs I/II','C2 vs C1','C3 vs C1','C3 vs C2'),
 function(k)rbindlist(lapply(c(3,6,12,24,36,60),function(t)cv(k,t)))))
prior <- fread(file.path(owner,'step6_cox_time_varying_subtype_effects.csv'))
parity <- merge(prior,all,by=c('Contrast','TimeMonths'),suffixes=c('_owner','_reproduced'))
stopifnot(nrow(parity)==18)
errs <- sapply(c('HR','Lower95','Upper95','PValue'),function(k) max(abs(parity[[paste0(k,'_owner')]]-parity[[paste0(k,'_reproduced')]])))
stopifnot(max(errs)<1e-9)
fwrite(parity,file.path(out,'review/time_varying_owner_parity.csv'))
fwrite(data.table(Field=names(errs),MaxAbsoluteError=errs),file.path(out,'review/time_varying_parity_summary.csv'))
fwrite(all,file.path(out,'tables/A_adjusted_HR_all_horizons.csv'))
shown <- all[TimeMonths %in% c(12,24) & Contrast!='C3 vs C2']
fwrite(shown,file.path(out,'tables/A_adjusted_HR_display.csv'))
fwrite(data.table(Term=names(b),Coefficient=as.numeric(b)),file.path(out,'tables/A_time_varying_coefficients.csv'))
fwrite(as.data.table(V,keep.rownames='Term'),file.path(out,'tables/A_time_varying_covariance.csv'))
file.copy(file.path(owner,'step6_cox_main_states_extended_sensitivity.csv'),file.path(out,'tables/A_extended_239_context.csv'),overwrite=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'review/association_sessionInfo.txt'))
message('Existing model reproduced; 18 owner subtype contrasts pass. All covariate HRs extracted.')
