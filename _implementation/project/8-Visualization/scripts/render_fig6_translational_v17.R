suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(digest)})
source('8-Visualization/scripts/fig_style.R')
# Quantitative grid: patient heterogeneity complements cross-study effect intervals.
# Same frozen scores, nine previously selected candidates, every Main patient.
out<-'8-Visualization/refinements/fig6_translational_20261101_v17'
for(d in c('figures','supplement','review','tables','qa'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
pal<-load_subtype_palette()$colors[c('C1','C2','C3')];invisible(load_study_palette())
src<-'8-Visualization/refinements/fig6_drug_evidence_20261101_v14'
inputs<-c(file.path(src,'tables/10_candidate_patient_heatmap.csv'),file.path(src,'tables/27_candidate_random_effects.csv'),SUBTYPE_TRUTH_CSV,
 '5-ClassifierValidation/publication/tables/step5_oncopredict_gdsc2_sample_centered_predictions.csv')
long<-fread(inputs[1],colClasses=c(Sample='character'));g<-fread(inputs[2]);candidates<-unique(g[,.(Drug,CandidateSubtype=Subtype)]);setorder(candidates,CandidateSubtype,Drug)
stopifnot(nrow(candidates)==9,!anyDuplicated(candidates$Drug))
main<-fread(SUBTYPE_TRUTH_CSV,colClasses=c(Sample='character'))[SubtypeClass=='Main',.(Sample,Study,FinalSubtype)]
setorder(main,FinalSubtype,Study,Sample);main[,PatientOrder:=seq_len(.N)]
z<-merge(long[Drug%chin%candidates$Drug],candidates,by='Drug')
stopifnot(nrow(z)==9*497,!anyDuplicated(z[,.(Drug,Sample)]),uniqueN(z$Sample)==497,all(is.finite(z$Z)))
ref<-main[match(z$Sample,Sample)];stopifnot(identical(z$FinalSubtype,ref$FinalSubtype),identical(z$Study,ref$Study),all(z$PatientOrder==ref$PatientOrder))
# Reconstruct the displayed score from the protected patient prediction interface.
owner<-fread(inputs[4],colClasses=c(Sample='character'));recon<-merge(melt(owner[,c('Sample',candidates$Drug),with=FALSE],id.vars='Sample',variable.name='Drug',value.name='OwnerScore'),main,by='Sample')
recon[,OwnerCentered:=OwnerScore-mean(OwnerScore),by=.(Drug,Study)];recon[,OwnerZ:=OwnerCentered/sd(OwnerCentered),by=Drug]
chk<-merge(z,recon[,.(Sample,Drug,OwnerScore,OwnerCentered,OwnerZ)],by=c('Sample','Drug'))
stopifnot(max(abs(chk$Score-chk$OwnerScore))<1e-10,max(abs(chk$Z-chk$OwnerZ))<1e-10)
z[,Drug:=factor(Drug,levels=rev(candidates$Drug))]
z[,CandidateGroup:=factor(paste(CandidateSubtype,'candidates'),levels=paste(c('C1','C2','C3'),'candidates'))]
z[,DisplayZ:=pmax(-2,pmin(2,Z))]
names_short<-setNames(sub('_[0-9]+$','',candidates$Drug),candidates$Drug)
cols<-c(C1='C1  |  n = 112',C2='C2  |  n = 223',C3='C3  |  n = 162')
base<-theme_classic(base_size=8.5,base_family='Arial')+theme(axis.text=element_text(size=8,color='#252525'),axis.title=element_text(size=8.5),plot.title=element_text(size=10,face='bold'),plot.subtitle=element_text(size=8),plot.caption=element_text(size=7.5,hjust=0),plot.title.position='plot',plot.caption.position='plot',strip.background=element_blank(),strip.text=element_text(size=8.5,face='bold'),legend.position='bottom',legend.text=element_text(size=8),legend.title=element_text(size=8),plot.margin=margin(3,3,3,3,'mm'))
theme_set(base)
p<-ggplot(z,aes(PatientOrder,Drug,fill=DisplayZ))+geom_tile(width=1,height=1)+
 facet_grid(CandidateGroup~FinalSubtype,scales='free',space='free',labeller=labeller(FinalSubtype=cols))+
 scale_x_continuous(expand=expansion(add=.5))+scale_y_discrete(labels=names_short,expand=expansion(add=.5))+
 scale_fill_gradient2(low='#397CA6',mid='white',high='#B66D7D',midpoint=0,limits=c(-2,2),breaks=c(-2,0,2),labels=c('≤ −2','0','≥ 2'),name='Relative predicted score (SD)',guide=guide_colorbar(barwidth=grid::unit(33,'mm'),barheight=grid::unit(2.5,'mm')))+
 labs(title='Patient-level patterns of drug candidates',subtitle='GDSC2-trained oncoPredict | one column per tumor',x=NULL,y=NULL,
 caption='Blue: lower relative score. Within subtypes: Study / sample order; no drug-score sorting.\nColor saturates at ±2 SD. Predictions, not measured patient responses or cross-drug potency.')+
 theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_blank(),panel.spacing=grid::unit(1.5,'mm'),strip.text.y=element_text(size=8))
manifest<-list();save<-function(id,p,w,h,d,track){
 pdf<-file.path(out,if(track=='Main')'figures' else 'supplement',paste0(id,'.pdf'));png<-file.path(out,'review',paste0(id,'.png'));tb<-file.path(out,'tables',paste0(id,'.csv'))
 cairo_pdf(pdf,width=w/25.4,height=h/25.4,family='Arial');print(p);dev.off()
 ragg::agg_png(png,width=w,height=h,units='mm',res=450,background='white');print(p);dev.off();fwrite(d,tb)
 manifest[[id]]<<-data.table(Panel=id,Track=track,PDF=pdf,PNG=png,SourceTable=tb,WidthMM=w,HeightMM=h,Reuse='New rendering of existing scores',Patients=uniqueN(d$Sample),Rows=nrow(d))
}
save('Figure6F_patient_drug_patterns',p,174,105,z,'Main')
# All nine candidates, not a visually selected representative subset.
a<-copy(z);a[,Drug:=factor(as.character(Drug),levels=candidates$Drug)]
p<-ggplot(a,aes(FinalSubtype,Z,color=FinalSubtype))+geom_hline(yintercept=0,color='#D6DADD',linewidth=.3)+
 geom_point(position=position_jitter(width=.18,height=0,seed=20261101),size=.45,alpha=.3)+
 geom_boxplot(width=.4,outlier.shape=NA,fill=NA,linewidth=.4)+facet_wrap(~Drug,ncol=3,labeller=labeller(Drug=names_short))+
 scale_color_manual(values=pal,guide='none')+scale_x_discrete(limits=c('C1','C2','C3'),labels=c(C1='C1\nn=112',C2='C2\nn=223',C3='C3\nn=162'))+
 labs(title='Patient distributions for all nine candidates',subtitle='Same scores and patients as Figure6F',x=NULL,y='Relative predicted score (SD)',caption='Every patient shown; no y-axis clipping. Boxes: median / IQR; whiskers: 1.5 × IQR.\nLower scores are relative predictions within a drug; no new significance tests.')+
 theme(panel.spacing=grid::unit(3,'mm'),strip.text=element_text(size=8))
save('Supplement_nine_candidate_distributions',p,174,168,a,'Supplement')
fwrite(main,file.path(out,'tables/patient_order.csv'));fwrite(candidates,file.path(out,'tables/nine_candidate_inventory.csv'))
# Evidence ledger documents available RNA evidence, without inventing protein/IHC validation.
expaths<-c('8-Visualization/refinements/fig6_duct_anchors_20261031_v11/tables/duct_anchor_expression_long.csv','8-Visualization/refinements/fig6_compact17_expression_20261101_v15/tables/compact17_expression_long.csv')
inputs<-c(inputs,expaths)
marker<-rbindlist(lapply(seq_along(expaths),function(j){d<-fread(expaths[j],colClasses=c(Sample='character'));r<-d[,.(Patients=uniqueN(Sample),Studies=uniqueN(Study),C1N=sum(FinalSubtype=='C1'),C2N=sum(FinalSubtype=='C2'),C3N=sum(FinalSubtype=='C3')),by=Gene];r[,`:=`(Group=if(j==1)'Duct-related reference' else 'Compact development candidate',RNAEvidence='Patient expression shown',ProteinEvidence='Not assessed in this bridge',IHCValidation='Not assessed in this bridge',SourceTable=expaths[j])];r}))
marker[,ClassificationUse:=fifelse(Group=='Compact development candidate','Separate fixed17 development candidate; D/E evaluate adaptive selection',fifelse(Gene%chin%c('S100P','CRP','CDH2'),'Three-marker RNA reference',fifelse(Gene=='MUC6','Four-marker sensitivity','Expression supplement only; missing in one source study')))]
stopifnot(nrow(marker)==24);fwrite(marker,file.path(out,'tables/marker_evidence_ledger.csv'))
verify<-data.table(Check=c('Nine candidates exact','All 497 patients each','Unique patient-drug pairs','Canonical labels and fixed patient order','Owner score maximum difference','Reconstructed Z maximum difference','Number of low saturated cells','Number of high saturated cells'),Value=c(nrow(candidates),uniqueN(z$Sample),nrow(z),1,max(abs(chk$Score-chk$OwnerScore)),max(abs(chk$Z-chk$OwnerZ)),sum(z$Z< -2),sum(z$Z>2)),Status='PASS')
fwrite(verify,file.path(out,'qa/numeric_verification.csv'));fwrite(rbindlist(manifest),file.path(out,'new_figure_manifest.csv'))
fwrite(data.table(Path=inputs,SHA256=vapply(inputs,function(p)digest(file=p,algo='sha256'),character(1))),file.path(out,'tables/input_hashes.csv'))
print(verify)
