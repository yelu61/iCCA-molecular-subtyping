# Figure6 evidence layers: existing source statistics; independent R PDF/PNG exports.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(grid);library(digest)})
source('8-Visualization/scripts/fig_style.R')
out<-'8-Visualization/refinements/fig6_evidence_20261103_v21';v20<-'8-Visualization/refinements/fig6_integrated_20261101_v20'
for(d in c('figures','supplement_pathology','supplement_drugs','figure5_handoff','review','tables','qa','provenance'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
inputs<-character();manifest<-list();checks<-list()
rd<-function(p){inputs<<-unique(c(inputs,p));cc<-names(fread(p,nrows=0));if('Sample'%in%cc)fread(p,colClasses=c(Sample='character'))else fread(p)}
r20<-function(n)rd(file.path(v20,'tables',paste0(n,'.csv')))
pal<-load_subtype_palette()$colors[c('C1','C2','C3')];invisible(load_study_palette());cls<-names(pal)
ink<-'#263E4B';muted<-'#7C8992';light<-'#E5EBEF'
base<-theme_classic(base_size=8,base_family='Arial')+theme(axis.text=element_text(size=7.5,color=ink),axis.title=element_text(size=8),plot.title=element_text(size=9,face='bold'),plot.subtitle=element_text(size=7.5),plot.caption=element_text(size=7,hjust=0),plot.title.position='plot',plot.caption.position='plot',legend.position='bottom',legend.title=element_blank(),legend.text=element_text(size=7),strip.background=element_blank(),strip.text=element_text(size=7.5,face='bold'),panel.spacing=unit(3,'mm'),plot.margin=margin(2,2,2,2,'mm'))
theme_set(base)
save_panel<-function(id,p,w,h,d,track='Main',more=NULL){folder<-switch(track,Main='figures',Pathology='supplement_pathology',Drugs='supplement_drugs',Figure5='figure5_handoff');pdf<-file.path(out,folder,paste0(id,'.pdf'));png<-file.path(out,'review',paste0(id,'.png'));tb<-file.path(out,'tables',paste0(id,'.csv'))
 draw<-function(){if(is.function(p))p()else print(p)};cairo_pdf(pdf,width=w/25.4,height=h/25.4,family='Arial');draw();dev.off();ragg::agg_png(png,width=w,height=h,units='mm',res=450,background='white');draw();dev.off();fwrite(d,tb)
 extra<-character();if(!is.null(more))for(n in names(more)){tt<-file.path(out,'tables',paste0(id,'_',n,'.csv'));fwrite(more[[n]],tt);extra<-c(extra,tt)}
 manifest[[id]]<<-data.table(Panel=id,Track=track,PDF=pdf,PNG=png,SourceTable=tb,AdditionalTables=paste(extra,collapse=';'),WidthMM=w,HeightMM=h,Rows=nrow(d),Reuse='New rendering');message(id)
}
ck<-function(n,b){stopifnot(all(b));checks[[n]]<<-data.table(Check=n,Pass=TRUE)}
models<-c('Anchor3','CompactWorkflow','JointWorkflow','FullWorkflow');mlabel<-c(Anchor3='RNA anchors',CompactWorkflow='Compact',JointWorkflow='Joint',FullWorkflow='Full workflow',Anchor4='Anchors + MUC6')
role<-c(Anchor3='#85919A',CompactWorkflow='#426C85',JointWorkflow='#243E50',FullWorkflow='#6F8C86',Anchor4='#AE8A58')
perf<-r20('joint_performance_ci');metrics<-c('BalancedAccuracy','C1Recall','C2Recall','C3Recall');metricnames<-c('Balanced accuracy','C1 recall','C2 recall','C3 recall')
performance<-function(d,title,subtitle){d[,Model:=factor(Model,levels=rev(names(mlabel)))];d[,Metric:=factor(Metric,levels=c(metrics,'C3Specificity','C3PPV'),labels=c(metricnames,'C3 specificity','C3 PPV'))]
 p<-ggplot(d,aes(100*Estimate,Model,color=Model))+geom_vline(xintercept=c(0,50,100),color=light,linewidth=.3)+geom_errorbar(aes(xmin=100*Lower95,xmax=100*Upper95),orientation='y',width=.12,linewidth=.45)+geom_point(size=1.5)+geom_text(aes(label=sprintf('%.1f',100*Estimate)),nudge_y=.25,size=7/.pt,show.legend=FALSE)+facet_wrap(~Metric,nrow=1)+scale_y_discrete(labels=mlabel)+scale_x_continuous(limits=c(-5,110),breaks=c(0,50,100))+scale_color_manual(values=role,guide='none')+labs(title=title,subtitle=subtitle,x='OOF performance (%)',y=NULL)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())
 list(p=p,data=d)
}
d<-perf[Model%chin%models&Metric%chin%metrics];a<-performance(copy(d),'Three-subtype identification across four workflows','Same 497 tumors · shared outer folds · conditional 95% intervals')
save_panel('Figure6D_three_subtype_identification',a$p,174,77,a$data);ck('D complete16 estimates and original intervals',nrow(d)==16)
s<-r20('Supplement_reference_sensitivity')[Model%chin%c('Anchor3','Anchor4')&Metric%chin%metrics];# source may carry display labels; use canonical source when necessary
if(nrow(s)!=8)s<-rd('5-ClassifierValidation/review/20261101_duct_compact_bridge_v1/tables/performance_ci.csv')[Model%chin%c('Anchor3','Anchor4')&Metric%chin%metrics]
a<-performance(copy(s),'Sensitivity to adding MUC6 to the RNA anchors','Baseline anchor trio versus four-marker reference · same 497 tumors')
save_panel('Path04_anchor_sensitivity',a$p,174,65,a$data,'Pathology');ck('Anchor sensitivity8 estimates',nrow(s)==8)
a<-performance(copy(perf[Model%chin%models&Metric%chin%c('C3Specificity','C3PPV')]),'C3 specificity and positive predictive value','Same reference-label prevalence · no new threshold selection')
save_panel('Path05_C3_secondary_metrics',a$p,120,79,a$data,'Pathology');ck('C3 secondary8 estimates',nrow(a$data)==8)

# Prediction error and calibration match current A owner, not historical publication models.
cl<-'6-ClinicalSurvival/review/20261030_full_subtype_clinical_value_v1'
cp<-rd(file.path(cl,'tables/model_performance_and_contrasts.csv'))[Metric%chin%c('Brier12','Brier24')];cp[,Horizon:=factor(sub('Brier','',Metric),levels=c('12','24'),labels=c('12 months','24 months'))]
cp[,Family:=factor(Family,levels=c('RidgePH','RidgeTimeVarying'),labels=c('Primary','Time-varying'))]
famcol<-c(Primary='#426C85','Time-varying'='#8C8176')
abs<-cp[Comparison%chin%c('M0','MC3','MAll')];abs[,Comparison:=factor(Comparison,levels=c('MAll','MC3','M0'))]
del<-cp[Comparison%chin%c('MC3-M0','MAll-M0','MAll-MC3')];del[,Comparison:=factor(Comparison,levels=c('MAll-MC3','MAll-M0','MC3-M0'))]
errplot<-function(a,title,labels){ggplot(a,aes(Estimate,Comparison,color=Family,shape=Family))+geom_vline(xintercept=0,color=light,linewidth=.3)+geom_errorbar(aes(xmin=Lower95,xmax=Upper95),orientation='y',position=position_dodge(width=.45),width=.1,linewidth=.45)+geom_point(position=position_dodge(width=.45),size=1.6)+facet_wrap(~Horizon,ncol=1)+scale_color_manual(values=famcol)+scale_shape_manual(values=c(Primary=16,'Time-varying'=2))+scale_y_discrete(labels=labels)+labs(title=title,x='Brier score / difference',y=NULL)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.title=element_text(size=8),panel.spacing=unit(5,'mm'))}
p<-errplot(abs,'Absolute prediction error',c(M0='Clinical',MC3='+ C3',MAll='+ C1/C2/C3'))+errplot(del,'Paired change in prediction error',c('MC3-M0'='Add C3','MAll-M0'='Add three subtypes','MAll-MC3'='Three subtypes vs C3'))+plot_layout(guides='collect')+plot_annotation(title='Prediction error for the current clinical models',subtitle='331 patients · same OOF risks as Figure6A · lower is better',theme=theme(plot.title=element_text(size=9,face='bold',family='Arial'),plot.subtitle=element_text(size=7.5,family='Arial')))
save_panel('Path10_clinical_prediction_error',p,174,112,cp,'Pathology');ck('Clinical error all24 rows',nrow(cp)==24)
cal<-rd(file.path(cl,'tables/calibration.csv'));cal[,Family:=factor(Family,levels=c('RidgePH','RidgeTimeVarying'),labels=c('Primary','Time-varying'))];cal[,Model:=factor(Model,levels=c('MAll','MC3','M0'))]
ca<-melt(cal,id.vars=c('Family','Model','Horizon','N'),measure.vars=c('Intercept','Slope'),variable.name='Parameter',value.name='Estimate');ca[,Facet:=factor(paste(Horizon,Parameter),levels=c('12 Intercept','24 Intercept','12 Slope','24 Slope'),labels=c('12 mo · joint intercept','24 mo · joint intercept','12 mo · slope','24 mo · slope'))]
ref<-unique(ca[,.(Facet,Ideal=ifelse(Parameter=='Slope',1,0))])
ca[,Y:=as.numeric(Model)+ifelse(Family=='Primary',-.16,.16)];ca[,LabelY:=Y+ifelse(Family=='Primary',-.16,.16)]
p<-ggplot(ca,aes(Estimate,Y,color=Family,shape=Family))+geom_vline(data=ref,aes(xintercept=Ideal),linetype=2,color='#B2BBC1',linewidth=.4)+geom_point(size=1.7)+geom_text(aes(y=LabelY,label=sprintf('%.2f',Estimate)),vjust=.5,size=7/.pt,show.legend=FALSE)+facet_wrap(~Facet,ncol=2,scales='free_x')+scale_y_continuous(breaks=3:1,labels=c('Clinical','+ C3','+ C1/C2/C3'),limits=c(.45,3.55))+scale_color_manual(values=famcol)+scale_shape_manual(values=c(Primary=16,'Time-varying'=2))+scale_x_continuous(expand=expansion(mult=.25))+labs(title='Calibration of the same out-of-fold risk predictions',subtitle='331 patients · jointly fitted IPCW logistic intercept / slope · point estimates only',x='Calibration parameter',y=NULL)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),panel.spacing=unit(7,'mm'))
save_panel('Path11_clinical_calibration',p,174,114,ca,'Pathology');ck('Calibration24 parameters',nrow(ca)==24)
# F: all patients and candidates; add verified annotation categories, no potency ranking.
aud<-file.path(out,'qa/candidate_audit');ev<-rd(file.path(aud,'candidate_evidence_summary.csv'))
z<-r20('Figure6F_drug_heterogeneity');ord<-ev$Drug
short<-setNames(sub('_[0-9]+$','',ord),ord);short['Sepantronium bromide_1941']<-'Sepantronium'
labels<-setNames(paste0(short[ord],'\n',ev$FunctionalClass),ord)
z[,Drug:=factor(Drug,levels=rev(ord))];z[,CandidateGroup:=factor(paste(CandidateSubtype,'nominees'),levels=paste(cls,'nominees'))]
summ<-z[,.(N=.N,Median=median(Z),Q1=quantile(Z,.25),Q3=quantile(Z,.75)),by=.(Drug,CandidateGroup,FinalSubtype)]
hm<-ggplot(z,aes(PatientOrder,Drug,fill=pmax(-2,pmin(2,Z))))+geom_tile(width=1,height=1)+facet_grid(CandidateGroup~FinalSubtype,scales='free',space='free',labeller=labeller(FinalSubtype=c(C1='C1 · 112',C2='C2 · 223',C3='C3 · 162')))+scale_y_discrete(labels=labels,expand=expansion(add=.5))+scale_x_continuous(expand=expansion(add=.5))+scale_fill_gradient2(low='#397CA6',mid='white',high='#B66D7D',limits=c(-2,2),breaks=c(-2,0,2),labels=c('≤−2','0','≥2'),name='Score (SD)',guide=guide_colorbar(title.position='top',barwidth=unit(20,'mm'),barheight=unit(2,'mm')))+labs(x=NULL,y=NULL,title='Individual tumors')+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=7),strip.text.y=element_blank(),plot.title=element_text(size=8),panel.spacing=unit(1,'mm'))
q<-ggplot(summ,aes(Median,Drug,color=FinalSubtype))+geom_vline(xintercept=0,linetype=2,color='#C0C8CD',linewidth=.3)+geom_errorbar(aes(xmin=Q1,xmax=Q3),orientation='y',position=position_dodge(width=.65),width=0,linewidth=.55)+geom_point(position=position_dodge(width=.65),size=1.3)+facet_grid(CandidateGroup~.,scales='free_y',space='free_y')+scale_y_discrete(expand=expansion(add=.5))+scale_x_continuous(limits=c(-2.3,2.3),breaks=c(-2,0,2))+scale_color_manual(values=pal)+labs(x='Score (SD)',y=NULL,title='Median / IQR')+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),strip.text.y=element_text(size=6.8),plot.title=element_text(size=8),legend.key.width=unit(3,'mm'),legend.spacing.x=unit(0,'mm'),panel.spacing=unit(1,'mm'))
p<-hm+q+plot_layout(widths=c(3.4,1))+plot_annotation(title='Predicted drug patterns across all three subtypes',subtitle='Nine nominees · 497 tumors · within-drug SD · lower = greater predicted sensitivity',theme=theme(plot.title=element_text(size=9,face='bold',family='Arial'),plot.subtitle=element_text(size=7.5,family='Arial')))
save_panel('Figure6F_annotated_drug_patterns',p,174,112,z,more=list(distribution_summary=summ,functional_annotations=ev[,.(Drug,FunctionalClass,SourceTarget,FunctionalCaveat)]));ck('F all4473 unchanged observations',nrow(z)==4473&uniqueN(z$Sample)==497)

# G: one row per drug; separate evidence types and retain every absent/discordant source.
progshort<-c('Unmapped','JAK / STAT','TGF-beta receptor','B-cell receptor','TP53 regulation','Mitosis (broad)','Checkpoints','MEK / MAPK','MAPK1 / MAPK3')
ev[,ProgramLabel:=progshort]
compactexpr<-function(s){s<-gsub('Up:','↑ ',s,fixed=TRUE);s<-gsub('Down:','↓ ',s,fixed=TRUE);s<-gsub('NA:','Missing: ',s,fixed=TRUE);s<-gsub('NS:','NS ',s,fixed=TRUE);s<-gsub('; ','\n',s,fixed=TRUE);s<-gsub('Missing: .*','Missing',s);s}
ev[,RNAPlot:=vapply(RNA,compactexpr,character(1))];ev[,ProteinPlot:=vapply(Protein,compactexpr,character(1))]
textbounds<-list()
drawG<-function(){grid.newpage();pushViewport(viewport(xscale=c(0,174),yscale=c(146,0)));on.exit(popViewport())
 tx<-function(s,x,y,size=7,bold=FALSE,col=ink,just='left'){gg<-textGrob(s,unit(x,'native'),unit(y,'native'),just=c(just,'center'),gp=gpar(fontfamily='Arial',fontsize=size,fontface=if(bold)'bold' else 'plain',col=col,lineheight=1.05));ww<-convertWidth(grobWidth(gg),'mm',TRUE);hh<-convertHeight(grobHeight(gg),'mm',TRUE);lo<-x-if(just=='center')ww/2 else if(just=='right')ww else 0;stopifnot(lo>=0,lo+ww<=174,y-hh/2>=0,y+hh/2<=146);grid.draw(gg)}
 tx('Candidate evidence: rationale, consistency and missing links',2,4,9,TRUE);tx('Original nine nominees · each evidence type remains separate',2,9,7.5)
 xs<-c(2,29,65,91,119,138,156)
 tx('Nominee',xs[1],17,7.5,TRUE);tx('Program context',xs[2],17,7.5,TRUE);tx('RNA',xs[3],17,7.5,TRUE);tx('Protein',xs[4],17,7.5,TRUE);tx('Meta CI',xs[5],17,7.1,TRUE);tx('Top 3',xs[6],17,7.1,TRUE);tx('PRISM',xs[7],17,7.1,TRUE)
 tx('Focal subtype',xs[2],21.5,6.8);tx('Focal subtype',xs[3],21.5,6.8);tx('Focal subtype',xs[4],21.5,6.8);tx('below 0',xs[5],21.5,6.8);tx('LOSO',xs[6],21.5,6.8);tx('C1/C2/C3',xs[7],21.5,6.8)
 for(i in seq_len(nrow(ev))){r<-ev[i];y<-29+(i-1)*11.9
  grid.rect(x=unit(87,'native'),y=unit(y,'native'),width=unit(172,'native'),height=unit(11.4,'native'),gp=gpar(fill=if(i%%2)'#F5F7F8' else 'white',col=NA))
  tx(short[r$Drug],xs[1],y-1.5,7.2,TRUE);tx(paste(r$CandidateSubtype,'nominee'),xs[1],y+2.6,6.8,FALSE,pal[r$CandidateSubtype])
  tx(r$ProgramLabel,xs[2],y-2.6,7)
  ps<-if(is.na(r$ProgramNES))'Not mapped' else sprintf('%s NES %.2f\nFDR %.3g',c(Up='↑',Down='↓',NS='NS')[r$ProgramState],r$ProgramNES,r$ProgramFDR)
  tx(ps,xs[2],y+2,6.8,FALSE,if(r$ProgramState=='Down')'#447D96' else ink)
  tx(r$RNAPlot,xs[3],y,6.8);tx(r$ProteinPlot,xs[4],y,6.8)
  tx(sprintf('%d / 2',r$PairwiseCIBelowZero),xs[5]+5,y-2,7.5,TRUE,just='center');tx(sprintf('I² %.0f–%.0f%%',r$MinI2,r$MaxI2),xs[5]+5,y+2.5,6.2,just='center')
  tx(sprintf('%d / %d',r$Top3Retained,r$Omissions),xs[6]+4,y,7.5,TRUE,just='center')
  tx(sprintf('%d / %d / %d',r$PRISMMeasuredC1,r$PRISMMeasuredC2,r$PRISMMeasuredC3),xs[7]+6,y-1.5,7.3,TRUE,just='center');tx(paste0(r$PRISMMeasured,' mapped / 3'),xs[7]+6,y+2.7,6.2,just='center')
 }
 tx('↑/↓: source FDR < 0.05; NS: not significant. Missing or unmapped is not a negative assay.',2,135,6.8)
 tx('LOSO: downstream study omissions. PRISM: measured lines with projected labels; coverage is not validation.',2,140,6.8)
 tx('Program mapping uses the full Figure2 GSEA background; no combined evidence score.',2,144,6.8)
}
save_panel('Figure6G_candidate_evidence_matrix',drawG,174,146,ev,more=list(program_all_subtypes=rd(file.path(aud,'candidate_program_all_subtypes.csv')),mechanism_map=rd(file.path(aud,'candidate_program_map.csv'))));ck('G nine separate evidence records',nrow(ev)==9&sum(ev$PRISMMeasured)==12&all(ev$PRISM24OriginalCandidateMeasurements==0))

# H: same candidate, now with study uncertainty and all available actual line assays.
hz<-z[as.character(Drug)=='Ulixertinib_1908'];hh<-r20('Figure6H_C3_hypothesis');hh[,Y:=ifelse(Gene=='MAPK1',2,1)+c(C1=.21,C2=0,C3=-.21)[Subtype]];hh[,Significant:=FDR<.05];hh[,Assay:=factor(Assay,levels=c('RNA','Protein'))]
x<-ggplot(hz,aes(FinalSubtype,Z,color=FinalSubtype))+geom_point(position=position_jitter(width=.16,height=0,seed=20261103),size=.4,alpha=.3)+geom_boxplot(width=.45,outlier.shape=NA,fill=NA,linewidth=.4)+scale_color_manual(values=pal,guide='none')+scale_x_discrete(labels=c(C1='C1 ·112',C2='C2 ·223',C3='C3 ·162'))+labs(title='Patient predictions',subtitle='497 tumors · lower relative score',x=NULL,y='Relative score (SD)')+theme(plot.title=element_text(size=8),plot.subtitle=element_text(size=7),axis.text.x=element_text(size=7))
hmta<-r20('SFig6F_drug_uncertainty')[Drug=='Ulixertinib_1908'];hmta[,Comp:=factor(paste('C3 vs',Other),levels=c('C3 vs C2','C3 vs C1'))]
m<-ggplot(hmta,aes(Effect,Comp))+geom_vline(xintercept=0,linetype=2,color=muted,linewidth=.3)+geom_errorbar(aes(xmin=Lower,xmax=Upper),orientation='y',width=.1,color=pal['C3'],linewidth=.55)+geom_point(color=pal['C3'],size=1.8)+geom_text(aes(label=sprintf('%.2f [%.2f, %.2f]  I² %.0f%%',Effect,Lower,Upper,I2_percent)),nudge_y=.25,size=7/.pt)+scale_x_continuous(limits=c(min(hmta$Lower)-.2,.15))+labs(title='Cross-study prediction effects',subtitle='Random effects · 95% intervals · Top 3 retained 7/7',x='Predicted score difference (SD)',y=NULL)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.title=element_text(size=8),plot.subtitle=element_text(size=7))
y<-ggplot(hh,aes(Effect,Y,color=Subtype,shape=Significant))+geom_vline(xintercept=0,linetype=2,linewidth=.3,color='#ACB8C0')+geom_point(size=1.7)+facet_wrap(~Assay,nrow=1,scales='free_x',labeller=labeller(Assay=c(RNA='RNA · log2 effect',Protein='Protein · SD effect')))+scale_color_manual(values=pal)+scale_shape_manual(values=c('TRUE'=16,'FALSE'=1),labels=c('FALSE'='FDR ≥ 0.05','TRUE'='FDR < 0.05'))+scale_y_continuous(breaks=c(2,1),labels=c('MAPK1','MAPK3'),limits=c(.5,2.5))+scale_x_continuous(expand=expansion(mult=.15))+labs(title='Both mapped genes',subtitle='C3 MAPK program: NES 1.10, FDR 0.497',x='Study-adjusted expression difference',y=NULL)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.title=element_text(size=8),plot.subtitle=element_text(size=7),legend.key.width=unit(2,'mm'),legend.text=element_text(size=6.5),strip.text=element_text(size=6.8),legend.box='vertical',legend.margin=margin(0,0,0,0))
pr<-rd(file.path(aud,'candidate_PRISM_all_cells.csv'))[Drug=='Ulixertinib_1908'];pr[,LineLabel:=factor(paste0(CellLineName,' · ',Prediction),levels=paste0(CellLineName,' · ',Prediction))]
pp<-ggplot(pr,aes(AUC,LineLabel,color=Prediction))+geom_point(size=2)+geom_text(aes(label=sprintf('%.3f',AUC)),nudge_y=.24,size=7/.pt,show.legend=FALSE)+scale_color_manual(values=pal,guide='none')+scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1))+labs(title='Observed PRISM responses',subtitle='Three lines · projected subtype labels',x='AUC (lower = greater sensitivity)',y=NULL)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.title=element_text(size=8),plot.subtitle=element_text(size=7),axis.text.y=element_text(size=7))
p<-(x+m+plot_layout(widths=c(1,1.45)))/(y+pp+plot_layout(widths=c(1.45,1)))+plot_layout(heights=c(1,1.2))+plot_annotation(title='Ulixertinib: a C3 hypothesis with unresolved evidence',subtitle='MAPK3 abundance supports follow-up; pathway enrichment and three-line assays do not confirm the prediction',theme=theme(plot.title=element_text(size=9,face='bold',family='Arial'),plot.subtitle=element_text(size=7,family='Arial')))
save_panel('Figure6H_C3_evidence_and_limits',p,174,130,hh,more=list(patient_predictions=hz,pairwise_meta=hmta,PRISM_observations=pr));ck('H all497 predictions12 molecular and3 line observations',nrow(hz)==497&nrow(hh)==12&nrow(pr)==3&all(is.finite(pr$AUC)))

# Drug supplement: original nine only, complete mapped measurement/absence grid.
prall<-rd(file.path(aud,'candidate_PRISM_all_cells.csv'));fwrite(prall,file.path(out,'tables/Drug06_PRISM_all_cells_original.csv'))
prall[,Drug:=factor(Drug,levels=rev(ord))];prall[,Line:=factor(paste0(CellLineName,'\nprojected ',Prediction),levels=c('SNU-1079\nprojected C2','HuH28\nprojected C3','HuCCT1\nprojected C3'))]
prall[,Label:=fifelse(is.finite(AUC),sprintf('%.3f',AUC),fifelse(ObservationStatus=='No accepted assay identity','No match','Missing'))]
p<-ggplot(prall,aes(Line,Drug,fill=AUC))+geom_tile(width=.96,height=.9,color='white')+geom_text(aes(label=Label),size=7.5/.pt,color=ink)+scale_y_discrete(labels=short)+scale_fill_gradient(low='#DCEEF5',high='#719FB5',limits=c(0,1),na.value='#EEF0F2',name='AUC (lower = greater sensitivity)',guide=guide_colorbar(title.position='top',barwidth=unit(25,'mm'),barheight=unit(2,'mm')))+labs(title='Original nine candidates: every mapped PRISM observation',subtitle='12 drug–line measurements / 27 combinations · no projected C1 line',x=NULL,y=NULL,caption='No match: no accepted identity mapping. Missing: matched assay without a value. Neither means resistance.')+theme(axis.line=element_blank(),axis.ticks=element_blank(),legend.title=element_text(size=7))
save_panel('Drug06_PRISM_observations',p,174,110,prall,'Drugs');ck('PRISM complete27 grid12 measured',nrow(prall)==27&sum(is.finite(prall$AUC))==12)
coverage<-melt(ev[,.(Drug,LegacyC1=PRISMMeasuredC1,LegacyC2=PRISMMeasuredC2,LegacyC3=PRISMMeasuredC3,New24Q2=PRISM24OriginalCandidateMeasurements)],id.vars='Drug',variable.name='Scope',value.name='Measured');coverage[,Drug:=factor(Drug,levels=rev(ord))]
p<-ggplot(coverage,aes(Scope,Drug))+geom_tile(fill='#F0F4F6',color='white',width=.94,height=.9)+geom_text(aes(label=Measured),size=8/.pt,color=ink)+scale_y_discrete(labels=short)+scale_x_discrete(position='top',labels=c(LegacyC1='Legacy C1\n0 lines',LegacyC2='Legacy C2\n1 line',LegacyC3='Legacy C3\n2 lines',New24Q2='24Q2\n19 lines'))+labs(title='Candidate assay coverage across the available cell-line resources',subtitle='Counts of mapped measurements for the original nine compounds',x=NULL,y=NULL,caption='Projected subtype labels. Zero indicates no mapped measurement, not resistance or a negative efficacy result.')+theme(axis.line=element_blank(),axis.ticks=element_blank())
save_panel('Drug07_PRISM_coverage',p,174,108,coverage,'Drugs');ck('New PRISM zero mappedoriginal9',sum(coverage[Scope=='New24Q2',Measured])==0)

fwrite(rbindlist(manifest,fill=TRUE),file.path(out,'new_figure_manifest.csv'));fwrite(rbindlist(checks),file.path(out,'qa/numeric_checks.csv'));fwrite(data.table(Path=inputs,SHA256=vapply(inputs,function(p)digest(file=p,algo='sha256'),character(1))),file.path(out,'tables/input_hashes.csv'))
