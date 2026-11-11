# Display-scale revision only: never rerun or overwrite v2 deconvolution.
suppressPackageStartupMessages({library(data.table);library(sandwich)})
old<-'2-TranscriptomicCharacterization/runs/20261027_tme_reaudit_v2'
run<-'2-TranscriptomicCharacterization/runs/20261027_tme_display_v3'
out<-file.path(run,'tables');dir.create(out,recursive=TRUE,showWarnings=FALSE)
d<-fread(file.path(old,'tables/primary_scores.csv'));d[,Sample:=as.character(Sample)]
for(e in parse('2-TranscriptomicCharacterization/scripts/summarize_tme_reaudit_20261027.R')){
 if(is.call(e)&&identical(e[[1]],as.name('<-'))&&as.character(e[[2]])%in%c('fit_one','contrast'))eval(e)
}
parity<-list();parity_n<-0L
members<-list('B cells total'=c('B cells naive','B cells memory'),'Macrophages total'=c('Macrophages M0','Macrophages M1','Macrophages M2'))
extra<-rbindlist(lapply(names(members),function(n)d[Method=='CIBERSORT' & Score%in%members[[n]],.(Value=sum(Value),Score=n,PreScaleEligible=all(PreScaleEligible)),by=.(Sample,Method,Study,Subtype)]))
d<-rbindlist(list(d,extra),fill=TRUE)
d[Method=='CIBERSORT',`:=`(StudyMean=mean(Value),StudySD=sd(Value)),by=.(Score,Study)]
d[Method=='CIBERSORT',ModelEligible:=PreScaleEligible & is.finite(Value) & is.finite(StudySD) & StudySD>1e-12]
d[Method=='CIBERSORT',ModelValue:=ifelse(ModelEligible,(Value-StudyMean)/StudySD,NA_real_)]
d[,Unit:='within-study SD']
d[Method=='ESTIMATE' & Score=='ESTIMATEScore',`:=`(Score='PurityProxy',Value=-Value,ModelValue=-ModelValue,StudyMean=-StudyMean)]
stopifnot(uniqueN(d$Sample)==433,!anyDuplicated(d[,.(Sample,Method,Score)]),uniqueN(d[,.(Method,Score)])==149)
pc<-contrast(d,check=TRUE);stopifnot(nrow(pc)==447)
orig<-fread(file.path(old,'tables/primary_contrasts.csv'))
pp<-merge(pc[Method=='ESTIMATE' & Score=='PurityProxy'],orig[Method=='ESTIMATE' & Score=='ESTIMATEScore'],by='Subtype')
stopifnot(max(abs(pp$Effect.x+pp$Effect.y))<1e-10,max(abs(pp$FDR.x-pp$FDR.y))<1e-10)
fwrite(d,file.path(out,'primary_scores.csv'));fwrite(pc,file.path(out,'primary_contrasts.csv'));fwrite(rbindlist(parity),file.path(out,'HC3_parity.csv'))
studies<-sort(unique(d$Study))
ps<-rbindlist(lapply(studies,function(st){z<-contrast(d[Study==st],min_per_group=3L);z[,Study:=st];z}))
lo<-rbindlist(lapply(studies,function(st){z<-contrast(d[Study!=st]);z[,LeftOut:=st];z}))
fwrite(ps,file.path(out,'per_study_contrasts.csv'));fwrite(lo,file.path(out,'loso_contrasts.csv'))
ss<-fread(file.path(old,'tables/sensitivity_contrasts.csv'))[Method!='CIBERSORT']
ss[Method=='ESTIMATE' & Score=='ESTIMATEScore',`:=`(Score='PurityProxy',Effect=-Effect,Lower=-Upper,Upper=-Lower)]
cf<-fread(file.path(old,'raw/CIBERSORT_fit_all.csv'));cf[,Sample:=as.character(Sample)]
ns<-cf[Study%in%studies & get('P-value')<.05,.N,by=Study][N>=2,Study]
z<-contrast(d[Method=='CIBERSORT' & Study%in%ns & Sample%in%cf[get('P-value')<.05,Sample]])
z[,Sensitivity:='CIBERSORT_fit_Plt0.05_estimable_studies_six_fixed_SD'];ss<-rbindlist(list(ss,z),fill=TRUE)
fwrite(ss,file.path(out,'sensitivity_contrasts.csv'))
fwrite(rbindlist(lapply(names(members),function(n)data.table(Score=n,NativeScore=members[[n]],Aggregation='Sum of non-overlapping native fractions'))),file.path(out,'new_CIBERSORT_membership.csv'))
fwrite(d[,.(N=.N,EligibleN=sum(ModelEligible),StudySD=StudySD[1]),by=.(Method,Score,Study)],file.path(out,'eligibility_audit.csv'))
sources<-c(file.path(old,'tables/primary_scores.csv'),file.path(old,'tables/primary_contrasts.csv'),file.path(old,'tables/sensitivity_contrasts.csv'),file.path(old,'raw/CIBERSORT_fit_all.csv'),'2-TranscriptomicCharacterization/scripts/summarize_tme_reaudit_20261027.R','2-TranscriptomicCharacterization/scripts/prepare_tme_grouped_v13.R')
fwrite(data.table(Path=sources,MD5=unname(tools::md5sum(sources))),file.path(run,'source_hashes.csv'))
writeLines(c('COMPUTED_REVIEW: display-scale revision, all deconvolution values unchanged.',paste(nrow(pc),'primary contrasts;',sum(pc$Estimable),'estimable.'),'PurityProxy is -ESTIMATEScore, replacing that test; no cosine conversion or new independent hypothesis.','CIBERSORT now within-study SD, not percentage points; model re-fitted and BH across22native+6existing+2new sums x3contrasts.','Anatomical grouping in E is visual, not an average of incompatible enrichment scores. Native readouts retained.','Zero-variance study/readouts excluded from SD inference, with actual N recorded. Original percentage-point results remain in v2.'),file.path(run,'CONTRACT.md'))
cat('PASS:',nrow(pc),'contrasts;',sum(pc$Estimable),'estimable\n')
