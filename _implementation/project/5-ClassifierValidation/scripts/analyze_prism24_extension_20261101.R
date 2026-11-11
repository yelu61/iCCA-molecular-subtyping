suppressPackageStartupMessages({library(data.table);library(jsonlite);library(digest)})
run<-'5-ClassifierValidation/review/20261101_drug_evidence_refresh_v1';ext<-file.path(run,'external_metadata');out<-file.path(run,'tables')
rd<-function(n)fread(file.path(out,paste0(n,'.csv')));wr<-function(x,n)fwrite(x,file.path(out,paste0(n,'.csv')))
proj<-rd('cell_line_current114_predictions');cm<-fread(file.path(ext,'PRISM_Repurposing_24Q2_Cell_Line_Info.csv'));tr<-fread(file.path(ext,'treatment_info.csv'));x<-fread(file.path(ext,'PRISM24Q2_response_icca19_allQC.csv'))
stopifnot(!anyDuplicated(cm$row_id),!anyDuplicated(tr$profile_id),!anyDuplicated(x[,.(row_id,profile_id)]),nrow(x)==91078)
x<-merge(x,cm[,.(row_id,ModelID=depmap_id,CellCulture=culture,CellScreen=screen,pool_id)],by='row_id');x<-merge(x,tr,by='profile_id');stopifnot(all(x$CellCulture==x$culture),all(x$CellScreen==x$screen),all(x$perturbation_type=='trt_cp'),uniqueN(x$ModelID)==19)
wr(x[,.(Records=.N,Passed=sum(PASS),FiniteLFC=sum(is.finite(LFC)),FiniteLFC_cb=sum(is.finite(LFC_cb)),Assays=uniqueN(broad_id)),by=.(ModelID,screen)],'PRISM24_line_QC_coverage')
pass<-x[PASS==TRUE & is.finite(LFC_cb)];collapsed<-pass[,.(NTechnical=.N,NReplicates=uniqueN(rep),MedianLFC_cb=median(LFC_cb),MedianLFC=median(LFC),MinLFC_cb=min(LFC_cb),MaxLFC_cb=max(LFC_cb),Name=unique(name)[1]),by=.(ModelID,row_id,pool_id,culture,screen,compound_plate,broad_id,dose)]
collapsed[,Eligible:=NReplicates>=2];wr(collapsed,'PRISM24_all_compound_line_summaries')
norm<-function(s)gsub('[[:space:]._-]','',toupper(s))
ann<-rd('drug_target_resource_annotations');md<-fread('oncoPredict_data/GLDS/GDSCv2/screened_compunds_rel_8.2.csv');allquality<-fread('5-ClassifierValidation/publication/tables/step5_oncopredict_gdsc2_cv_quality.csv')
allann<-data.table(Drug=allquality$Drug,ModelCompound=sub('_[0-9]+$','',allquality$Drug),DRUG_ID=as.integer(sub('^.*_','',allquality$Drug)));allann<-merge(allann,md,by='DRUG_ID')
cp<-unique(tr[,.(broad_id,name,perturbation_type)])
links<-rbindlist(lapply(seq_len(nrow(allann)),function(i){a<-allann[i];aliases<-unique(c(a$ModelCompound,a$DRUG_NAME,trimws(strsplit(a$SYNONYMS,',',fixed=TRUE)[[1]])))
 if(a$ModelCompound=='Selumetinib')aliases<-c('Selumetinib','AZD6244','ARRY-142886');if(a$ModelCompound=='5-Fluorouracil')aliases<-c(aliases,'Fluorouracil')
 hit<-cp[norm(name)%chin%norm(aliases)];if(!nrow(hit))return(NULL);hit[,`:=`(Drug=a$Drug,ModelCompound=a$ModelCompound,IdentitySource='Exact/punctuation-normalized name or GDSC release 8.2 synonym')];hit
}))
oldcross<-rd('chemical_crosswalk_all_links')[Status=='Accepted name/alias link'];oldcross[,Stem:=sub('^.*BRD:(BRD-[A-Z][0-9]+).*$','\\1',Assay)]
cp[,Stem:=sub('^(BRD-[A-Z][0-9]+).*$','\\1',broad_id)];stem<-merge(oldcross[,.(Drug,ModelCompound,Stem)],cp,by='Stem',allow.cartesian=TRUE)
stem<-stem[,.(broad_id,name,perturbation_type,Drug,ModelCompound,IdentitySource='Same Broad compound stem as accepted legacy assay; assay variants retained')]
links<-unique(rbind(links,stem));links<-links[,.(IdentitySource=paste(unique(IdentitySource),collapse='; ')),by=.(Drug,ModelCompound,broad_id,name,perturbation_type)]
wr(links,'PRISM24_model_identity_links')
experimental<-links[perturbation_type=='trt_cp'];matched<-merge(collapsed,experimental,by='broad_id',allow.cartesian=TRUE)
matched<-merge(matched,proj[,.(ModelID,CellLineName,PatientID,Prediction,MaximumProbability,Margin)],by='ModelID')
wr(matched,'PRISM24_all_matched_model_summaries');wr(merge(pass,experimental[,.(broad_id,Drug,ModelCompound)],by='broad_id',allow.cartesian=TRUE),'PRISM24_matched_technical_records')
chosen<-rd('original_candidate_retention');prior<-fread('8-Visualization/refinements/fig6_application_suite_20261031_v7/tables/candidate_selection_provenance.csv')
coverage<-merge(proj,matched[Eligible==TRUE,.(NewSharedModels=uniqueN(Drug),NewNamedCompounds=uniqueN(ModelCompound),NewOriginalCandidates=uniqueN(Drug[Drug%chin%chosen$Drug]),NewCandidateOrPrior=uniqueN(Drug[Drug%chin%prior$Drug])),by=ModelID],by='ModelID',all.x=TRUE)
coverage[,InNewScreenMetadata:=ModelID%chin%cm$depmap_id];coverage[,InNewScreenResponse:=ModelID%chin%x$ModelID]
for(k in c('NewSharedModels','NewNamedCompounds','NewOriginalCandidates','NewCandidateOrPrior'))coverage[is.na(get(k)),(k):=0L]
legacy<-rd('PRISM_candidate_model_cells')[is.finite(AUC)];coverage<-merge(coverage,legacy[,.(LegacyCandidateModelCount=uniqueN(Drug),LegacyNamedCompoundCount=uniqueN(sub('_[0-9]+$','',Drug))),by=ModelID],by='ModelID',all.x=TRUE)
coverage[is.na(LegacyCandidateModelCount),`:=`(LegacyCandidateModelCount=0L,LegacyNamedCompoundCount=0L)];wr(coverage,'external_resource_line_coverage')
summary<-list(ReleaseDOI='10.6084/m9.figshare.25917643.v1',EligibleMetadataLines=sum(coverage$InNewScreenMetadata),ObservedLines=uniqueN(x$ModelID),DistinctLocalPatientIDs=uniqueN(proj[ModelID%chin%x$ModelID]$PatientID),ProjectedSubtypeCounts=as.list(table(proj[ModelID%chin%x$ModelID]$Prediction)),AllSelectedTechnicalRecords=nrow(x),FailedPASSRecords=sum(!x$PASS),SharedExperimentalModels=uniqueN(experimental$Drug),SharedExperimentalNamedCompounds=uniqueN(experimental$ModelCompound),EligibleMatchedConditionCells=nrow(matched[Eligible==TRUE]),MatchedOriginalNineCandidates=uniqueN(experimental[Drug%chin%chosen$Drug]$Drug),MatchedCandidateOrPriorModels=uniqueN(experimental[Drug%chin%prior$Drug]$Drug),Endpoint='2.5 micromolar / 5 days; median ComBat-corrected log2 viability (LFC_cb), exact source dose retained',Conclusion='More eligible lines but zero original/prior candidate coverage; does not validate the candidate panel',SourceHash=digest(file=file.path(ext,'PRISM24Q2_response_icca19_allQC.csv'),algo='sha256'))
write_json(summary,file.path(run,'PRISM24_extension_summary.json'),pretty=TRUE,auto_unbox=TRUE);print(summary)
