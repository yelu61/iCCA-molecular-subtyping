#!/usr/bin/env Rscript
# Implements runs/20261027_tme_reaudit_v2/ANALYSIS_SPEC.md; no visual selection.
suppressPackageStartupMessages({library(data.table);library(sandwich)})
run<-'2-TranscriptomicCharacterization/runs/20261027_tme_reaudit_v2'
out<-file.path(run,'tables');dir.create(out,recursive=TRUE,showWarnings=FALSE)
without_ciber<-'--without-cibersort'%in%commandArgs(TRUE)
meta<-fread(file.path(run,'inputs/meta_all.csv'))[SubtypeClass=='Main'];meta[,Sample:=as.character(Sample)]
primary_studies<-setdiff(unique(meta$Study),'2023Deng')
stopifnot(nrow(meta)==497L,nrow(meta[Study%in%primary_studies])==433L)
read_scores<-function(path){
 d<-fread(path);d[,Sample:=as.character(Sample)]
 # Metadata is validated against the canonical Main table, not trusted from score files.
 d<-d[Sample%in%meta$Sample]
 if('Study'%in%names(d))stopifnot(all(d$Study==meta$Study[match(d$Sample,meta$Sample)]))
 if('Subtype'%in%names(d))stopifnot(all(d$Subtype==meta$Subtype[match(d$Sample,meta$Sample)]))
 d[,`:=`(Study=meta$Study[match(Sample,meta$Sample)],Subtype=meta$Subtype[match(Sample,meta$Sample)])]
 if(!'Eligible'%in%names(d))d[,Eligible:=TRUE]
 if('ModelEligible'%in%names(d))d[,Eligible:=Eligible & ModelEligible]
 if(!'Coverage'%in%names(d))d[,Coverage:=1]
 if(!'View'%in%names(d))d[,View:='perstudy_native_allRNA']
 d
}
xcov<-fread(file.path(run,'coverage/xCell_signature_coverage.csv'))
xcover<-xcov[,.(Coverage=median(Fraction)),by=.(View,Study,Score=CellType)]
xglobal<-xcov[,.(Coverage=median(Fraction)),by=.(View,Study)]
xglobal<-xglobal[,.(Score=c('ImmuneScore','StromaScore','MicroenvironmentScore'),Coverage=Coverage),by=.(View,Study)]
xcover<-rbind(xcover,xglobal)
read_xcell<-function(view){
 d<-read_scores(file.path(run,'raw',paste0('xCell_',view,'Main_scores_long.csv')))
 d[,Coverage:=NULL];d[,CoverageStudy:=if(grepl('^joint_',view))'joint' else Study]
 cv<-xcover[View==view];setnames(cv,'Study','CoverageStudy')
 d<-merge(d,cv[,.(CoverageStudy,Score,Coverage)],by=c('CoverageStudy','Score'),all.x=TRUE,sort=FALSE)
 stopifnot(!anyNA(d$Coverage));d[,CoverageStudy:=NULL];d
}
native<-read_scores(file.path(run,'raw/native_programs_scores_long.csv'))
epic<-read_scores(file.path(run,'raw/EPIC_scores_long.csv'))
epic[,Coverage:=1] # Audited global native reference; do not invent cell-specific cutoffs.
xcprimary<-read_xcell('joint_common_six_allRNA')
base<-rbindlist(list(native,epic,xcprimary),fill=TRUE)
ciber<-NULL
if(!without_ciber){
 cpath<-file.path(run,'raw/CIBERSORT_scores_long.csv');stopifnot(file.exists(cpath))
 ciber<-read_scores(cpath);stopifnot(uniqueN(ciber$Score)==22L,uniqueN(ciber$Sample)==497L)
 cats<-list('Total T cells'=c('T cells CD8','T cells CD4 naive','T cells CD4 memory resting','T cells CD4 memory activated','T cells follicular helper','T cells regulatory (Tregs)','T cells gamma delta'),
 'Total NK cells'=c('NK cells resting','NK cells activated'),'B / plasma'=c('B cells naive','B cells memory','Plasma cells'),
 'Mono / macrophages'=c('Monocytes','Macrophages M0','Macrophages M1','Macrophages M2'),
 'Dendritic total'=c('Dendritic cells resting','Dendritic cells activated'),
 'Granulocyte / mast'=c('Mast cells resting','Mast cells activated','Eosinophils','Neutrophils'))
 stopifnot(length(unlist(cats))==22L,!anyDuplicated(unlist(cats)),setequal(unlist(cats),unique(ciber$Score)))
 derived<-rbindlist(lapply(names(cats),function(n)ciber[Score%in%cats[[n]],.(Value=sum(Value),Eligible=all(Eligible),Coverage=1,Score=n,Method='CIBERSORT',View='derived_disjoint_composition'),by=.(Sample,Study,Subtype)]))
 stopifnot(all(abs(derived[,sum(Value),by=Sample]$V1-1)<1e-8))
 fwrite(rbindlist(lapply(names(cats),function(n)data.table(Category=n,NativeScore=cats[[n]]))),file.path(out,'CIBERSORT_category_membership.csv'))
 ciber<-rbindlist(list(ciber,derived),fill=TRUE);base<-rbindlist(list(base,ciber),fill=TRUE)
}
stopifnot(!anyDuplicated(base[,.(Sample,Method,Score)]))
# Eligibility concerns the primary modeling population, not biological validation.
prepare<-function(d,coverage_cut=.8){
 d<-copy(d)
 d[,InputEligible:=!is.na(Eligible)&Eligible&is.finite(Value)]
 d[,CoverageEligible:=is.finite(Coverage)&Coverage>=coverage_cut]
 d[Method%in%c('EPIC','CIBERSORT'),CoverageEligible:=TRUE]
 d[Method=='Prespecified program',CoverageEligible:=Coverage==1 & InputEligible]
 d[,PreScaleEligible:=InputEligible&CoverageEligible]
 d[,`:=`(StudyN=sum(PreScaleEligible),StudyMean=mean(Value[PreScaleEligible]),StudySD=sd(Value[PreScaleEligible])),by=.(Method,Score,Study)]
 d[,VariableStudy:=is.finite(StudySD)&StudySD>1e-12]
 d[,ModelEligible:=PreScaleEligible & (Method=='CIBERSORT'|VariableStudy)]
 d[,Unit:=ifelse(Method=='CIBERSORT','percentage points','within-study SD')]
 d[,ModelValue:=NA_real_]
 d[ModelEligible & Method!='CIBERSORT',ModelValue:=(Value-StudyMean)/StudySD]
 d[ModelEligible & Method=='CIBERSORT',ModelValue:=100*Value]
 d[,EligibilityReason:=fcase(!InputEligible,'Nonfinite or computationally unavailable',!CoverageEligible,'Below prespecified coverage',Method!='CIBERSORT'&!VariableStudy,'NE: no within-study variance',default='Eligible')]
 d
}
primary<-prepare(base[Study%in%primary_studies])
fwrite(primary,file.path(out,'primary_scores.csv'))
audit<-primary[,.(N=.N,EligibleN=sum(ModelEligible),MinCoverage=min(Coverage),StudySD=StudySD[1],Reason=paste(unique(EligibilityReason),collapse=';')),by=.(Method,Score,Study)]
fwrite(audit,file.path(out,'eligibility_audit.csv'))
parity<-list();parity_n<-0L
fit_one<-function(d,cl,min_per_group=1L,check=FALSE){
 d<-d[ModelEligible & is.finite(ModelValue)]
 n1<-sum(d$Subtype==cl);n0<-sum(d$Subtype!=cl);nst<-uniqueN(d$Study)
 empty<-list(Effect=NA_real_,SE=NA_real_,Lower=NA_real_,Upper=NA_real_,P=NA_real_,N=nrow(d),NTarget=n1,NOthers=n0,Studies=nst,Estimable=FALSE,Reason='Insufficient groups or observations')
 if(min(n1,n0)<min_per_group)return(empty)
 dat<-data.frame(y=d$ModelValue,g=as.numeric(d$Subtype==cl),Study=factor(d$Study))
 X<-if(nst>1)model.matrix(~g+Study,dat) else model.matrix(~g,dat)
 if(qr(X)$rank<ncol(X)||nrow(X)<=ncol(X)){empty$Reason<-'Rank-deficient or no residual degrees of freedom';return(empty)}
 inv<-solve(crossprod(X));b<-inv%*%crossprod(X,dat$y);res<-as.numeric(dat$y-X%*%b);h<-rowSums((X%*%inv)*X)
 if(any(1-h<1e-10)){empty$Reason<-'Unit leverage prevents HC3 covariance';return(empty)}
 V<-inv%*%crossprod(X,X*as.numeric((res/(1-h))^2))%*%inv
 se<-sqrt(max(0,V[2,2]));df<-nrow(X)-ncol(X);eff<-as.numeric(b[2]);crit<-qt(.975,df)
 if(!is.finite(se)||se<=1e-14){empty$Reason<-'Zero or nonfinite robust SE';return(empty)}
 if(check){
   fit<-if(nst>1)lm(y~g+Study,dat) else lm(y~g,dat)
   sv<-sandwich::vcovHC(fit,type='HC3');e2<-unname(coef(fit)['g']);s2<-sqrt(sv['g','g'])
   stopifnot(abs(eff-e2)<1e-8,abs(se-s2)<1e-8)
   parity_n<<-parity_n+1L;parity[[parity_n]]<<-data.table(Method=d$Method[1],Score=d$Score[1],Subtype=cl,EffectDelta=eff-e2,SEDelta=se-s2,N=nrow(d))
 }
 list(Effect=eff,SE=se,Lower=eff-crit*se,Upper=eff+crit*se,P=2*pt(-abs(eff/se),df),N=nrow(d),NTarget=n1,NOthers=n0,Studies=nst,Estimable=TRUE,Reason='Estimated')
}
contrast<-function(d,check=FALSE,min_per_group=1L){
 ans<-d[,rbindlist(lapply(c('C1','C2','C3'),function(cl)c(list(Subtype=cl,Unit=Unit[1]),fit_one(.SD,cl,min_per_group,check)))),by=.(Method,Score)]
 ans[,FDR:=p.adjust(P,method='BH'),by=Method]
 ans[,`:=`(FamilyTests=sum(is.finite(P)),FamilySpecified=.N),by=Method]
 setcolorder(ans,c('Method','Score','Subtype','Effect','SE','Lower','Upper','P','FDR','N','Studies','Unit'))
 ans
}
pc<-contrast(primary,check=TRUE);fwrite(pc,file.path(out,'primary_contrasts.csv'))
if(!without_ciber){
  derived_pc<-pc[Method=='CIBERSORT' & Score%in%names(cats)]
  closure<-derived_pc[,.(CategoryCount=.N,Complete=sum(Estimable),EffectSum=sum(Effect),SameN=uniqueN(N)==1L,N=N[1]),by=Subtype]
  stopifnot(all(closure$CategoryCount==6L),all(closure$Complete==6L),all(closure$SameN),all(abs(closure$EffectSum)<1e-8))
  cp<-primary[Method=='CIBERSORT' & ModelEligible==TRUE]
  stopifnot(max(abs(cp$ModelValue-100*cp$Value))<1e-10,all(cp$Unit=='percentage points'))
  closure[,`:=`(DisjointExhaustive=TRUE,PerSampleFractionClosure=TRUE,RawPercentagePointScale=TRUE)]
  fwrite(closure,file.path(out,'CIBERSORT_composition_closure_QA.csv'))
}
fwrite(rbindlist(parity),file.path(out,'HC3_lm_sandwich_parity.csv'))
ps<-rbindlist(lapply(primary_studies,function(st){d<-contrast(primary[Study==st],min_per_group=3L);d[,Study:=st];d}))
fwrite(ps,file.path(out,'per_study_contrasts.csv'))
lo<-rbindlist(lapply(primary_studies,function(st){d<-contrast(primary[Study!=st]);d[,LeftOut:=st];d}))
fwrite(lo,file.path(out,'loso_contrasts.csv'))
equal<-ps[,.(EstimableStudies=sum(Estimable),PositiveStudies=sum(Effect>0,na.rm=TRUE),NegativeStudies=sum(Effect<0,na.rm=TRUE),EqualStudyMean=if(any(Estimable))mean(Effect[Estimable]) else NA_real_,MinEffect=if(any(Estimable))min(Effect[Estimable]) else NA_real_,MaxEffect=if(any(Estimable))max(Effect[Estimable]) else NA_real_),by=.(Method,Score,Subtype,Unit)]
fwrite(equal,file.path(out,'equal_study_summary.csv'))
ss<-list();sk<-0L
add_sens<-function(d,name){sk<<-sk+1L;q<-contrast(d);q[,Sensitivity:=name];ss[[sk]]<<-q}
add_sens(prepare(base[Study%in%primary_studies],0),'six_study_relaxed_coverage')
add_sens(prepare(base[Method!='xCell']),'seven_study_native_coverage80')
add_sens(prepare(base[Method!='xCell'],0),'seven_study_native_relaxed_coverage')
for(v in c('joint_common_allRNA','perstudy_common_allRNA','perstudy_native_allRNA')){
  xd<-read_xcell(v)
  add_sens(prepare(xd[Study%in%primary_studies]),paste0('xCell_',v,'_six_study_coverage80'))
  add_sens(prepare(xd),paste0('xCell_',v,'_seven_study_coverage80'))
  add_sens(prepare(xd,0),paste0('xCell_',v,'_seven_study_relaxed_coverage'))
}
ef<-fread(file.path(run,'raw/EPIC_fit_diagnostics.csv'));ef[,Sample:=as.character(Sample)]
add_sens(primary[Method=='EPIC' & Sample%in%ef[is.finite(pearsonR)&pearsonR>=.3,Sample]],'EPIC_Pearson_ge0.3_six_study_fixed_scale')
if(!without_ciber){
 fp<-file.path(run,'raw/CIBERSORT_fit_diagnostics.csv')
 if(!file.exists(fp))fp<-file.path(run,'raw/CIBERSORT_fit_all.csv')
 if(!file.exists(fp))stop('CIBERSORT completed score inputs require fit diagnostics')
 cf<-fread(fp);cf[,Sample:=as.character(Sample)];cf<-cf[Sample%in%meta$Sample]
 cf[,Pass:=get('P-value')<.05]
 cfcounts<-cf[,.(FullN=.N,FitPassedN=sum(Pass)),by=Study]
 cfcounts[,`:=`(PrimaryStudy=Study%in%primary_studies,FitSubsetStudyEstimable=FitPassedN>=2L)]
 cfcounts[,Reason:=ifelse(FitSubsetStudyEstimable,'At least two retained samples','Singleton or empty fit-selected study; HC3 study intercept not estimable')]
 fwrite(cfcounts,file.path(out,'CIBERSORT_fit_subset_study_counts.csv'))
 add_sens(primary[Method=='CIBERSORT' & Sample%in%cf[Pass==TRUE,Sample]],'CIBERSORT_fit_Plt0.05_six_study')
 add_sens(prepare(ciber[Sample%in%cf[Pass==TRUE,Sample]]),'CIBERSORT_fit_Plt0.05_seven_study')
 # Preserve the unmodified fit-subset NE result above. This explicitly named
 # additional sensitivity excludes only singleton studies, which create h=1.
 enough<-cfcounts[FitSubsetStudyEstimable==TRUE,Study]
 add_sens(primary[Method=='CIBERSORT' & Study%in%enough & Sample%in%cf[Pass==TRUE,Sample]],'CIBERSORT_fit_Plt0.05_estimable_studies_six')
 add_sens(prepare(ciber[Study%in%enough & Sample%in%cf[Pass==TRUE,Sample]]),'CIBERSORT_fit_Plt0.05_estimable_studies_seven')
 fwrite(cf[,.(N=.N,Passed=sum(Pass),Rate=mean(Pass),MedianP=median(get('P-value'))),by=.(Study,Subtype)],file.path(out,'CIBERSORT_fit_selection_by_study_subtype.csv'))
 sel<-cf[Study%in%primary_studies]
 fit0<-glm(Pass~factor(Study),data=sel,family=binomial());fit1<-glm(Pass~factor(Study)+Subtype,data=sel,family=binomial())
 capture.output(anova(fit0,fit1,test='LRT'),file=file.path(out,'CIBERSORT_fit_selection_LRT.txt'))
 fwrite(as.data.table(summary(fit1)$coefficients,keep.rownames='Term'),file.path(out,'CIBERSORT_fit_selection_logistic.csv'))
}
sens<-rbindlist(ss,fill=TRUE);fwrite(sens,file.path(out,'sensitivity_contrasts.csv'))
selection_path<-'8-Visualization/refinements/fig2_transcriptomic_20261027_v12_tme/tables/E_selection_map.csv'
if(file.exists(selection_path)){
  selection<-fread(selection_path)
  ls<-lo[,.(LOSOEstimable=sum(Estimable),LOSOMin=if(any(Estimable))min(Effect[Estimable]) else NA_real_,LOSOMax=if(any(Estimable))max(Effect[Estimable]) else NA_real_),by=.(Method,Score,Subtype)]
  selected<-merge(selection,pc,by=c('Method','Score'));selected<-merge(selected,ls,by=c('Method','Score','Subtype'))
  selected[,LOSOSignReversal:=(Effect>0 & LOSOMin<0)|(Effect<0 & LOSOMax>0)]
  fwrite(selected,file.path(out,'selected_LOSO_robustness.csv'))
  sscomp<-merge(selected[,.(Method,Score,Subtype,PrimaryEffect=Effect,PrimaryFDR=FDR,PrimaryN=N)],sens,by=c('Method','Score','Subtype'))
  sscomp[,SignReversal:=sign(Effect)!=sign(PrimaryEffect)]
  fwrite(sscomp,file.path(out,'selected_sensitivity_comparison.csv'))
}
xc<-merge(pc[Method=='xCell',.(Method,Score,Subtype,PrimaryEffect=Effect,PrimaryFDR=FDR,PrimaryN=N)],sens[Method=='xCell'],by=c('Method','Score','Subtype'))
xc[,`:=`(EffectDelta=Effect-PrimaryEffect,SignChanged=sign(Effect)!=sign(PrimaryEffect))];fwrite(xc,file.path(out,'xCell_view_effect_differences.csv'))
# Old-v1 comparison is diagnostic only: sources, scales and populations changed.
old<-fread('2-TranscriptomicCharacterization/runs/20261026_tme_expansion_v1/tables/native_adjusted_contrasts.csv');setnames(old,'Feature','Score')
old[Method=='MCP-counter',Score:=sub('_MCPcounter$','',gsub('_',' ',Score))]
old[Method=='MCP-counter',Score:=sub(' MCPcounter$','',Score)]
old[Method=='Immune GSVA',Score:=sub('^Immune 28::','',Score)]
old[Method=='ESTIMATE',Score:=sub('_ESTIMATE$','',Score)]
chg<-merge(old[,.(Method,Score,Subtype,OldEffect=Effect,OldFDR=FDR,OldN=N)],pc,by=c('Method','Score','Subtype'))
chg[,ComparisonNote:=fifelse(Method=='CIBERSORT','Old standardized composition versus new percentage points: magnitude not directly comparable','Input/calibration/coverage and sample population differ; not a replication test')]
chg[,SignChanged:=sign(OldEffect)!=sign(Effect)];fwrite(chg,file.path(out,'old_v1_effect_comparison.csv'))
input_paths<-c(file.path(run,'ANALYSIS_SPEC.md'),file.path(run,'inputs/meta_all.csv'),file.path(run,'raw/native_programs_scores_long.csv'),file.path(run,'raw/EPIC_scores_long.csv'),file.path(run,'raw/EPIC_fit_diagnostics.csv'),file.path(run,'raw',paste0('xCell_',c('joint_common_six_allRNA','joint_common_allRNA','perstudy_common_allRNA','perstudy_native_allRNA'),'Main_scores_long.csv')),file.path(run,'coverage/xCell_signature_coverage.csv'),'2-TranscriptomicCharacterization/runs/20261026_tme_expansion_v1/tables/native_adjusted_contrasts.csv','2-TranscriptomicCharacterization/scripts/summarize_tme_reaudit_20261027.R')
if(!without_ciber)input_paths<-c(input_paths,file.path(run,'raw/CIBERSORT_scores_long.csv'),fp)
fwrite(data.table(Path=input_paths,MD5=unname(tools::md5sum(input_paths))),file.path(out,'statistics_input_hashes.csv'))
capture.output(sessionInfo(),file=file.path(run,'logs/statistics_sessionInfo.txt'))
writeLines(c(paste('Status:',if(without_ciber)'PRIMARY_OTHER_METHODS_COMPUTED_CIBERSORT_PENDING' else 'ALL_METHODS_COMPUTED_REVIEW'),
 paste('Primary samples:',uniqueN(primary$Sample),'contrasts:',nrow(pc),'estimable:',sum(pc$Estimable)),
 paste('HC3 vs lm/sandwich parity checks:',parity_n),
 'Computational QA does not validate biological fractions or source independence.'),file.path(out,'statistics_status.txt'))
cat('COMPLETE:',nrow(pc),'primary contrasts;',nrow(sens),'sensitivity contrasts;',parity_n,'HC3 parity checks\n')
