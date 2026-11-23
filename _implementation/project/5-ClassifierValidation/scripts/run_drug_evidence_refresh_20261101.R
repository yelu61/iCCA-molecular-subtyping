#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(sandwich);library(e1071);library(jsonlite);library(digest)})
out<-'5-ClassifierValidation/review/20261101_drug_evidence_refresh_v1'
dir.create(file.path(out,'tables'),recursive=TRUE,showWarnings=FALSE)
inputs<-character(); rd<-function(p){inputs<<-unique(c(inputs,p));fread(p)}
wr<-function(d,n)fwrite(d,file.path(out,'tables',paste0(n,'.csv')))
cross_join<-function(a,b){a<-copy(a);b<-copy(b);a[,JoinKey:=1L];b[,JoinKey:=1L];z<-merge(a,b,by='JoinKey',allow.cartesian=TRUE);z[,JoinKey:=NULL];z}
checks<-list();ck<-function(name,test){checks[[name]]<<-data.table(Check=name,Pass=isTRUE(test));stopifnot(isTRUE(test))}
owner<-function(n)rd(paste0('5-ClassifierValidation/publication/tables/step5_oncopredict_',n,'.csv'))
prev<-'8-Visualization/refinements/fig6_application_suite_20261031_v7/tables'
master<-rd('governance/data/subtype_master_table.csv')[SubtypeClass=='Main',.(Sample,Study,FinalSubtype)]
quality<-owner('gdsc2_cv_quality');batch<-owner('batch_method_sensitivity');drugs<-quality$Drug
raw<-merge(master,owner('gdsc2_studywise_predictions'),by='Sample');cent<-merge(master,owner('gdsc2_sample_centered_predictions'),by='Sample')
ck('497_main_unique',nrow(master)==497 && !anyDuplicated(master$Sample));ck('198_drug_models',length(drugs)==198)
source_effect<-owner('study_adjusted_effects');prior_screen<-rd(file.path(prev,'all_subtype_candidate_screen.csv'))
chosen<-prior_screen[Selected==TRUE];prior<-rd(file.path(prev,'candidate_selection_provenance.csv'));ids<-prior$Drug
classes<-c('C1','C2','C3');contr<-list('C2-C1'=c(SubtypeC2=1),'C3-C1'=c(SubtypeC3=1),'C3-C2'=c(SubtypeC3=1,SubtypeC2=-1))
fit_rows<-function(y,meta,drug,typ,extra=NULL){
 d<-data.frame(y=y,Subtype=factor(meta$FinalSubtype,levels=classes),Study=factor(meta$Study));if(!is.null(extra))d<-cbind(d,extra)
 f<-lm(reformulate(c('Subtype','Study',names(extra)),response='y'),d);b<-coef(f);stopifnot(all(is.finite(b)))
 v<-vcovHC(f,type='HC3');rbindlist(lapply(names(contr),function(co){L<-setNames(rep(0,length(b)),names(b));L[names(contr[[co]])]<-contr[[co]];est<-sum(L*b);se<-sqrt(drop(t(L)%*%v%*%L));data.table(Drug=drug,ScoreType=typ,Contrast=co,N=nrow(d),Estimate=est,SE_HC3=se,CI_low=est-1.96*se,CI_high=est+1.96*se,P=2*pnorm(abs(est/se),lower.tail=FALSE),DesignRank=f$rank,DesignColumns=length(b),DesignCondition=kappa(model.matrix(f)))}))
}
target_rows<-function(e){rbindlist(lapply(classes,function(st){
 pairs<-switch(st,C1=c('C2-C1','C3-C1'),C2=c('C2-C1','C3-C2'),C3=c('C3-C1','C3-C2'));signs<-switch(st,C1=c(-1,-1),C2=c(1,-1),C3=c(1,1))
 rbindlist(lapply(1:2,function(j){a<-copy(e[Contrast==pairs[j]]);s<-signs[j];a[,`:=`(Subtype=st,Other=setdiff(strsplit(pairs[j],'-',fixed=TRUE)[[1]],st),TargetEstimate=s*Estimate,TargetLower=if(s==1)CI_low else -CI_high,TargetUpper=if(s==1)CI_high else -CI_low)];a}))
}))}
screen_one<-function(e,d){
 l<-melt(d,id.vars=c('Sample','Study','FinalSubtype'),measure.vars=drugs,variable.name='Drug',value.name='Value',variable.factor=FALSE)
 l[,Residual:=Value-mean(Value),by=.(Drug,Study)];sdtab<-l[,.(WithinStudySD=sd(Residual)),by=Drug]
 ds<-rbindlist(lapply(classes,function(st){a<-l[,.(N_target=sum(FinalSubtype==st),N_rest=sum(FinalSubtype!=st),Difference=mean(Value[FinalSubtype==st])-mean(Value[FinalSubtype!=st])),by=.(Drug,Study)][N_target>=2&N_rest>=2];a[,Subtype:=st];a}))
 dir<-ds[,.(Studies=.N,StudiesLower=sum(Difference<0),FractionLower=mean(Difference<0)),by=.(Drug,Subtype)]
 ct<-target_rows(e[ScoreType=='SampleCentered']);s<-ct[,.(BothBelowZero=all(TargetUpper<0),BothFDR=all(FDR_global<.05),WeakerDifference=max(TargetEstimate)),by=.(Drug,Subtype)]
 for(t in list(quality,batch,sdtab))s<-merge(s,t,by='Drug');s<-merge(s,dir,by=c('Drug','Subtype'),all.x=TRUE)
 s[,WeakerStandardizedDifference:=WeakerDifference/WithinStudySD]
 s[,Eligible:=BothBelowZero & BothFDR & Pearson_r>=.3 & Primary_vs_pooled_standardize_Spearman>=.6 & Primary_vs_pooled_eb_Spearman>=.6 & Studies>=3 & FractionLower>=.7]
 s[is.na(Eligible),Eligible:=FALSE];s[,`:=`(Compound=toupper(sub('_[0-9]+$','',Drug)),Selected=FALSE,EligibleRank=NA_integer_)]
 for(st in classes){a<-s[Subtype==st&Eligible==TRUE][order(WeakerStandardizedDifference,Drug)];a<-a[!duplicated(Compound)];s[Subtype==st,EligibleRank:=match(Drug,a$Drug)];s[Subtype==st&Drug%chin%head(a$Drug,3),Selected:=TRUE]}
 list(screen=s,directions=ds,sd=sdtab,targets=ct)
}
all_effects<-list();all_screens<-list();all_targets<-list()
for(omit in c('None',sort(unique(master$Study)))){
 message('Downstream association exclusion: ',omit)
 a<-raw[Study!=omit];d<-cent[Study!=omit]
 e<-rbindlist(lapply(drugs,function(dr)rbind(fit_rows(a[[dr]],a,dr,'Raw'),fit_rows(d[[dr]],d,dr,'SampleCentered'))))
 e[,FDR_global:=p.adjust(P,'BH')];e[,FDR_within_score_contrast:=p.adjust(P,'BH'),by=.(ScoreType,Contrast)]
 s<-screen_one(e,d)
 if(omit=='None'){
   test<-merge(e,source_effect,by=c('Drug','ScoreType','Contrast'),suffixes=c('_new','_source'));wr(test,'baseline_effect_reproduction')
   ck('baseline_effects_reproduced',max(abs(test$Estimate_new-test$Estimate_source))<1e-10 && max(abs(test$SE_HC3_new-test$SE_HC3_source))<1e-10 && max(abs(test$FDR_global_new-test$FDR_global_source))<1e-9)
   ss<-merge(s$screen,prior_screen,by=c('Drug','Subtype'),suffixes=c('_new','_source'))
   ck('baseline_all_screen_flags_reproduced',all(ss$Selected_new==ss$Selected_source)&all(ss$Eligible_new==ss$Eligible_source));wr(s$sd,'fixed_display_sd')
 }
 e[,OmittedStudy:=omit];s$screen[,OmittedStudy:=omit];s$targets[,OmittedStudy:=omit]
 all_effects[[omit]]<-e;all_screens[[omit]]<-s$screen;all_targets[[omit]]<-s$targets
}
effects<-rbindlist(all_effects);screens<-rbindlist(all_screens);targets<-rbindlist(all_targets)
wr(effects,'all_study_exclusion_effects');wr(screens,'all_study_exclusion_screens');wr(targets,'all_study_exclusion_target_contrasts')
ret<-merge(screens[OmittedStudy!='None'],chosen[,.(Drug,Subtype)],by=c('Drug','Subtype'))
wr(ret,'original_candidate_exclusion_details')
wr(ret[,.(Omissions=.N,EligibleRetained=sum(Eligible),TopThreeRetained=sum(Selected),BothNegativeRetained=sum(BothBelowZero),BestEligibleRank=if(all(is.na(EligibleRank)))NA_integer_ else min(EligibleRank,na.rm=TRUE),WorstEligibleRank=if(all(is.na(EligibleRank)))NA_integer_ else max(EligibleRank,na.rm=TRUE)),by=.(Drug,Subtype)],'original_candidate_retention')

# Pair-specific per-study uncertainty for the nine originally selected candidates.
pair<-rbindlist(lapply(seq_len(nrow(chosen)),function(i){dr<-chosen$Drug[i];st<-chosen$Subtype[i]
 rbindlist(lapply(setdiff(classes,st),function(other)rbindlist(lapply(sort(unique(master$Study)),function(study){d<-cent[Study==study];a<-d[[dr]][d$FinalSubtype==st];b<-d[[dr]][d$FinalSubtype==other];n1<-length(a);n0<-length(b);eligible<-min(n1,n0)>=2
  if(eligible){v1<-var(a)/n1;v0<-var(b)/n0;se<-sqrt(v1+v0);df<-(v1+v0)^2/(v1^2/(n1-1)+v0^2/(n0-1));est<-mean(a)-mean(b);lo<-est-qt(.975,df)*se;hi<-est+qt(.975,df)*se}else{est<-se<-df<-lo<-hi<-NA_real_}
  data.table(Drug=dr,Subtype=st,Other=other,Study=study,N_target=n1,N_other=n0,Eligible=eligible,Estimate=est,SE=se,DF=df,CI_low=lo,CI_high=hi)
 }))))
}))
wr(pair,'candidate_pairwise_study_effects')
met<-pair[Eligible==TRUE,{
 w<-1/SE^2;mu<-sum(w*Estimate)/sum(w);Q<-sum(w*(Estimate-mu)^2);k<-.N;df<-k-1;C<-sum(w)-sum(w^2)/sum(w);tau<-max(0,(Q-df)/C);ww<-1/(SE^2+tau);mm<-sum(ww*Estimate)/sum(ww);q<-sum(ww*(Estimate-mm)^2)/df;se<-sqrt(max(1,q)/sum(ww));
 .(Studies=k,N_target=sum(N_target),N_other=sum(N_other),Estimate=mm,SE=se,CI_low=mm-qt(.975,df)*se,CI_high=mm+qt(.975,df)*se,Tau2=tau,I2_percent=if(Q>0)max(0,100*(Q-df)/Q)else 0,P=2*pt(abs(mm/se),df,lower.tail=FALSE),Method='DL tau2; modified Hartung-Knapp t interval')
},by=.(Drug,Subtype,Other)]
met[,FDR_18_postselection:=p.adjust(P,'BH')];wr(met,'candidate_pairwise_meta_analysis')

# Fixed current Figure5 full model; reproduce owner's external-cohort transform.
mp<-'5-ClassifierValidation/runs/bulk-rnaseq__20261030_step5_consensus_v1/models/candidate_models.rds';inputs<-c(inputs,mp);mh<-digest(file=mp,algo='sha256');m<-readRDS(mp)$Consensus_SVM
ck('current_model_114',length(m$genes)==114)
cellmeta<-rd('5-ClassifierValidation/Model.csv')[OncotreeSubtype=='Intrahepatic Cholangiocarcinoma']
cp<-'5-ClassifierValidation/OmicsExpressionTPMLogp1HumanAllGenes.csv';inputs<-c(inputs,cp);hd<-names(fread(cp,nrows=0));symbols<-sub(' \\(.*$','',hd);cols<-hd[symbols%chin%m$genes]
xx<-fread(cp,select=c('ModelID','is_default_entry',cols))[ModelID%chin%cellmeta$ModelID & is_default_entry==TRUE];ck('28_unique_eligible_lines',nrow(xx)==28&&!anyDuplicated(xx$ModelID))
mat<-t(as.matrix(xx[,..cols]));rownames(mat)<-sub(' \\(.*$','',cols);colnames(mat)<-xx$ModelID
ck('all_114_features_present',setequal(rownames(mat),m$genes));mat<-mat[m$genes,,drop=FALSE]
standardize<-function(x){mu<-rowMeans(x,na.rm=TRUE);sdv<-apply(x,1,sd,na.rm=TRUE);mu[!is.finite(mu)]<-0;sdv[!is.finite(sdv)|sdv<1e-8]<-1;z<-sweep(sweep(x,1,mu,'-'),1,sdv,'/');z[!is.finite(z)]<-0;t(z)}
predict_p<-function(z)attr(predict(m$model,z,probability=TRUE),'probabilities')[,classes,drop=FALSE]
inp<-standardize(mat);p<-predict_p(inp)
project<-data.table(ModelID=rownames(p),Prediction=classes[max.col(p,ties.method='first')],ProbabilityC1=p[,1],ProbabilityC2=p[,2],ProbabilityC3=p[,3],MaximumProbability=apply(p,1,max),Margin=apply(p,1,function(x)diff(sort(x,decreasing=TRUE)[2:1])),MatchedFeatures=114L,TotalFeatures=114L,ProjectionOnly=TRUE,IndependentTruthAvailable=FALSE)
project<-merge(project,cellmeta[,.(ModelID,CellLineName,PatientID)],by='ModelID')
old<-rd('5-ClassifierValidation/publication/tables/step5_external_projection_sample_predictions.csv')[Platform=='DepMap/CCLE']
project<-merge(project,old[,.(ModelID=Sample,Old123Prediction=Prediction)],by='ModelID');project[,LabelChanged:=Prediction!=Old123Prediction]
ck('probabilities_normalized',max(abs(rowSums(p)-1))<1e-10)
leave<-rbindlist(lapply(seq_len(ncol(mat)),function(j){z<-standardize(mat[,-j,drop=FALSE]);pp<-predict_p(z);data.table(ExcludedLine=colnames(mat)[j],ModelID=rownames(pp),Prediction=classes[max.col(pp,ties.method='first')],ProbabilityC1=pp[,1],ProbabilityC2=pp[,2],ProbabilityC3=pp[,3])}))
leave<-merge(leave,project[,.(ModelID,FullPrediction=Prediction)],by='ModelID');leave[,Agrees:=Prediction==FullPrediction];wr(leave,'cell_line_recentering_sensitivity')
project<-merge(project,leave[,.(RecenteringAgreement=mean(Agrees)),by=ModelID],by='ModelID')
ep<-'1-Clustering/data/combined_expr_before_correction.rds';inputs<-c(inputs,ep);env<-new.env();load(ep,env);expr<-env$combined_expr
train<-t(expr[m$genes,master$Sample,drop=FALSE]);for(st in unique(master$Study)){ii<-which(master$Study==st);train[ii,]<-sweep(train[ii,,drop=FALSE],2,m$study_offsets[[st]][m$genes],'-')}
train<-sweep(sweep(train,2,m$center[m$genes],'-'),2,m$scale[m$genes],'/')
distance<-function(a,b)sqrt(pmax(outer(rowSums(a*a),rowSums(b*b),'+')-2*tcrossprod(a,b),0))
dt<-distance(train,train);diag(dt)<-Inf;nn<-apply(dt,1,min);dl<-distance(inp,train);line_nn<-apply(dl,1,min)
domain<-data.table(ModelID=rownames(inp),NearestDevelopmentDistance=line_nn,DevelopmentNNPercentile=vapply(line_nn,function(v)mean(nn<=v),numeric(1)),Development95th=unname(quantile(nn,.95)))
project<-merge(project,domain,by='ModelID');wr(data.table(Sample=master$Sample,Subtype=master$FinalSubtype,NearestOtherDevelopmentDistance=nn),'development_domain_reference')
pc<-prcomp(train,center=TRUE,scale.=FALSE);pt<-predict(pc,inp)
wr(rbind(data.table(Unit=rownames(train),Domain='Development tumor',Subtype=master$FinalSubtype,PC1=pc$x[,1],PC2=pc$x[,2]),data.table(Unit=rownames(inp),Domain='Projected cell line',Subtype=project$Prediction[match(rownames(inp),project$ModelID)],PC1=pt[,1],PC2=pt[,2])),'model_input_domain_pca')
wr(data.table(Gene=m$genes,MissingMeasurements=rowSums(!is.finite(mat)),CellLineSD=apply(mat,1,sd),DevelopmentSD=apply(train,2,sd)),'cell_line_feature_coverage')

# Rebuild assay-resolved line-by-drug cells with new labels, without reinterpreting NA.
pr<-rd('5-ClassifierValidation/Drug_sensitivity_AUC_(PRISM_Repurposing_Secondary_Screen).csv');setnames(pr,1,'ModelID');overlap<-intersect(pr$ModelID,project$ModelID)
ck('prism_overlap_three',setequal(overlap,c('ACH-000209','ACH-000808','ACH-000976')))
pr<-melt(pr[ModelID%chin%overlap],id.vars='ModelID',variable.name='Assay',value.name='AUC',variable.factor=FALSE);pr[,PRISMCompound:=sub(' \\(BRD:.*$','',Assay)]
pr<-merge(pr,project[,.(ModelID,CellLineName,Prediction,MaximumProbability,Margin)],by='ModelID');ck('2397_unchanged_prism_observations',sum(is.finite(pr$AUC))==2397)
project[,PRISMAvailable:=ModelID%chin%overlap];wr(project,'cell_line_current114_predictions');wr(pr,'PRISM_all_assay_cells')
cross<-rd(file.path(out,'tables/chemical_crosswalk_all_links.csv'))[Status=='Accepted name/alias link']
linked<-merge(cross,pr,by=c('Assay','PRISMCompound'),allow.cartesian=TRUE);wr(linked,'PRISM_all_mapped_model_cells')
shared<-pr[Assay%chin%cross$Assay];wr(shared,'PRISM_unique_shared_assay_cells')
grid<-cross_join(prior[,.(Drug,SelectedFor)],project[PRISMAvailable==TRUE,.(ModelID,CellLineName,Prediction)])
cand<-merge(grid,linked[,.(Drug,ModelID,Assay,PRISMCompound,AUC,MatchType,Source)],by=c('Drug','ModelID'),all.x=TRUE,allow.cartesian=TRUE)
cand[,ObservationStatus:=fifelse(is.na(Assay),'No accepted assay identity',fifelse(is.finite(AUC),'Measured','Matched assay; missing measurement'))];wr(cand,'PRISM_candidate_model_cells')
wr(unique(cand[,.(PRISMCompound,Assay,ModelID,CellLineName,Prediction,AUC,ObservationStatus)]),'PRISM_candidate_unique_cells')

# Limited composition / proliferation-proxy sensitivities on identical matched patients.
tme<-rd('2-TranscriptomicCharacterization/runs/20261027_tme_reaudit_v2/tables/primary_scores.csv')[Method=='ESTIMATE'&Score%chin%c('ImmuneScore','StromalScore')&ModelEligible==TRUE]
tme<-dcast(tme,Sample~Score,value.var='ModelValue');cov<-merge(master,tme,by='Sample',all.x=TRUE)
cov[,MKI67:=as.numeric(expr['MKI67',Sample])];cov[,MKI67_z:=as.numeric(scale(MKI67)),by=Study]
sens<-rbindlist(lapply(list(Composition=c('ImmuneScore','StromalScore'),ProliferationProxy='MKI67_z'),function(cols){d<-merge(cent,cov[,c('Sample',cols),with=FALSE],by='Sample');d<-d[complete.cases(d[,..cols])]
 rbindlist(lapply(ids,function(dr){a<-fit_rows(d[[dr]],d,dr,'MatchedBaseline');b<-fit_rows(d[[dr]],d,dr,'Adjusted',as.data.frame(d[,..cols]));rbind(a,b)}))[,Covariates:=paste(cols,collapse='+')]
}))
sens[,FDR_60:=p.adjust(P,'BH'),by=.(Covariates,ScoreType)];wr(sens,'candidate_phenotype_sensitivity');wr(cov,'phenotype_covariates')

# Mechanism-to-gene map; generic cytotoxic mechanisms are not forced to one target.
ann<-rd(file.path(out,'tables/drug_target_resource_annotations.csv'))[Drug%chin%ids]
gene_map<-list(AZD1208=c('PIM1','PIM2','PIM3'),AZD7762=c('CHEK1','CHEK2'),Dabrafenib='BRAF',Ibrutinib='BTK',Irinotecan='TOP1',JAK1_8709='JAK1',LY2109761=c('TGFBR1','TGFBR2'),PD0325901=c('MAP2K1','MAP2K2'),'PRIMA-1MET'='TP53',Selumetinib=c('MAP2K1','MAP2K2'),'Sepantronium bromide'='BIRC5',Trametinib=c('MAP2K1','MAP2K2'),Ulixertinib=c('MAPK1','MAPK3'))
gm<-rbindlist(lapply(seq_len(nrow(ann)),function(i){a<-ann[i];g<-gene_map[[a$ModelCompound]];if(is.null(g))g<-NA_character_;a<-a[rep(1,length(g))];a[,Gene:=g];a[,EvidenceRole:=fifelse(is.na(Gene),'Non-gene-specific cytotoxic mechanism',fifelse(ModelCompound%chin%c('PRIMA-1MET','Sepantronium bromide'),'Mechanism-associated gene; not a direct-binding claim','Annotated pharmacologic target'))];a[,AnnotationSource:='GDSC release 8.2; standard gene-symbol mapping'];a[ModelCompound=='LY2109761',AnnotationSource:='Corrected from TGFB1 ligand to receptors; PMID 16452224; IUPHAR ligand 8075'];a}))
wr(gm,'candidate_mechanism_gene_map')
rna<-rd('1-Clustering/publication/tables/step1_main_only_study_adjusted_deg.csv')[,.(Gene,Subtype,RNAEffect=logFC_subtype_minus_others,RNASE=SE,RNAP=PValue,RNAFDR=FDR,RNAN=NTotal)]
protein<-rd('4-Multiomic/publication/tables/step4_mass_spec_study_adjusted_deg_main_subtypes.csv')[,.(Gene=Protein,Subtype,ProteinEffect=StandardizedDifference,ProteinP=PValue,ProteinFDR=FDR)]
genomic<-rd('4-Multiomic/publication/tables/step4_driver_snv_indel_study_stratified_tests.csv')[,.(Gene,Subtype,GenomicOR=CMHOddsRatio,GenomicLow=CMHLower95,GenomicHigh=CMHUpper95,GenomicFDR=CMHFDR,GenomicTargetN=TargetN,GenomicOtherN=OtherN)]
evidence<-cross_join(gm,data.table(Subtype=classes))
for(tb in list(rna,protein,genomic))evidence<-merge(evidence,tb,by=c('Gene','Subtype'),all.x=TRUE,allow.cartesian=TRUE)
evidence[,`:=`(RNAAvailable=is.finite(RNAEffect),ProteinAvailable=is.finite(ProteinEffect),GenomicAvailable=is.finite(GenomicOR),DependencyMeasured=FALSE,DrugEfficacyValidated=FALSE)]
wr(evidence,'candidate_target_multiomic_evidence')
wr(rd('4-Multiomic/publication/tables/step4_assay_specific_inclusion.csv'),'multiomic_assay_denominators')
cnv<-rd('4-Multiomic/publication/tables/step4_driver_cnv_frequency_main_subtypes_2022dong.csv');wr(cnv[Gene%chin%gm$Gene],'candidate_target_CNV_descriptive')
rppa<-rd('4-Multiomic/publication/tables/step4_tcga_rppa_subtype_associations.csv');wr(rppa[grepl('MEK|MAPK|ERK|CHK|CHK1|P53|BRAF|PIM|BTK|SMAD|JAK',Target,ignore.case=TRUE)],'candidate_pathway_RPPA_source_context')
# Broad RPPA text matches are context only, not an automatic phosphosite/target crosswalk.
tide<-rd('2-TranscriptomicCharacterization/publication/tables/step2_tide_study_specific_effects.csv');wr(tide,'TIDE_study_effects')
wr(rd('2-TranscriptomicCharacterization/publication/tables/step2_tide_meta_analysis.csv'),'TIDE_meta_analysis')
ck('model_unchanged',identical(mh,digest(file=mp,algo='sha256')))
wr(data.table(Path=unique(inputs),SHA256=vapply(unique(inputs),function(p)digest(file=p,algo='sha256'),character(1))),'analysis_input_hashes')
wr(rbindlist(checks),'technical_checks')
write_json(list(Status='COMPUTED_REVIEW',MainPatients=nrow(master),DrugModels=length(drugs),AssociationExclusions=length(all_effects)-1,OriginalCandidates=nrow(chosen),CurrentCellLines=nrow(project),CurrentProjectionCounts=as.list(table(project$Prediction)),LabelsChanged=sum(project$LabelChanged),PRISMLines=length(overlap),PRISMObserved=sum(is.finite(pr$AUC)),SharedAssays=uniqueN(shared$Assay),MappedCandidateObservedModelCells=sum(is.finite(cand$AUC)),MappedCandidateUniqueObservedAssayCells=nrow(unique(cand[is.finite(AUC),.(ModelID,Assay)])),ModelSHA256=mh),file.path(out,'summary.json'),pretty=TRUE,auto_unbox=TRUE)
message('Computed review artifacts complete: ',out)
