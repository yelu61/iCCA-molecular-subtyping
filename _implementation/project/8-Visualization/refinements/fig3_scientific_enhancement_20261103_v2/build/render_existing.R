suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr);library(readr);library(patchwork);library(scales)})
o<-normalizePath('8-Visualization/refinements/fig3_scientific_enhancement_20261103_v2');root<-normalizePath('.');tt<-file.path(o,'source_data');fig<-file.path(o,'figures')
studycols<-read_csv(file.path(root,'0-Data/qc/study_color_dictionary.csv'),show_col_types=FALSE)
ap<-read_csv(file.path(root,'7-SingleCell/runs/20261103_figure3_opening_umap_v3/source_data/annotation_palette.csv'),show_col_types=FALSE);cellcols<-setNames(ap$Color,ap$CellTypeDisplay);cellcols['Myeloid']<-'#B77B48';cellcols['T_NK']<-'#42669B'
reg<-read_csv(file.path(tt,'shared_program_registry.csv'),show_col_types=FALSE)
labels<-setNames(gsub('^ME','',reg$CompactLabel),reg$ProgramID)
labels['WGCNA_MEdarkturquoise']<-'darkturquoise | Antiviral / IFN'
labels['WGCNA_MEcyan']<-'cyan | Respiratory transport'
labels['WGCNA_MElightgreen']<-'lightgreen | Fatty-acid / peroxisome'
orderp<-reg$ProgramID
selected<-c('Pattern_P2','Pattern_P3','Pattern_P4','WGCNA_MEcyan','WGCNA_MEmagenta','WGCNA_MEdarkturquoise');fine<-c('Epithelial','Hepatocyte','CD8 T','Other T','NK-like','B cell','Plasma cell','Monocytes/macrophages','Dendritic-like','Fibroblast','Pericyte-like','Endothelial')
flabel<-setNames(fine,fine);flabel['Monocytes/macrophages']<-'Mono / macrophage';flabel['Other T']<-'Other T†';flabel['NK-like']<-'NK-like†'
blabel<-c(Epithelial='Epithelial',Myeloid='Myeloid',T_NK='T / NK',Fibroblast='Fibroblast')
theme_set(theme_classic(base_family='Arial',base_size=7)+theme(plot.title=element_text(size=10,face='bold',margin=margin(b=4)),plot.subtitle=element_text(size=6.5,colour='#555555',margin=margin(b=5)),plot.caption=element_text(size=5.5,hjust=0,colour='#555555'),axis.text=element_text(size=6,colour='#333333'),axis.title=element_text(size=7),strip.text=element_text(size=7,face='bold'),strip.background=element_rect(fill='#F1F2F3',colour=NA),legend.title=element_text(size=6.5),legend.text=element_text(size=6),legend.position='bottom',panel.spacing=grid::unit(3,'mm'),plot.margin=margin(5,5,5,5)))
exports<-list();emit<-function(p,name,w,h,source,role='main_candidate'){
 saveRDS(p,file.path(fig,paste0(name,'.rds')))

 cairo_pdf(file.path(fig,paste0(name,'.pdf')),width=w/25.4,height=h/25.4,family='Arial');print(p);dev.off()
 svglite::svglite(file.path(fig,paste0(name,'.svg')),width=w/25.4,height=h/25.4);print(p);dev.off()
 ragg::agg_png(file.path(fig,paste0(name,'.png')),width=w,height=h,units='mm',res=300,background='white');print(p);dev.off()
 ragg::agg_tiff(file.path(fig,paste0(name,'.tiff')),width=w,height=h,units='mm',res=600,compression='lzw',background='white');print(p);dev.off()
 exports[[length(exports)+1]]<<-data.frame(Stem=name,WidthMM=w,HeightMM=h,Source=source,Role=role);cat('rendered',name,'\n')
}
fill_score<-function(limits=NULL)scale_fill_gradientn(colours=c('#F1F1F1','#B7D7DA','#4B9CA8','#245366'),limits=limits,name='UCell score',na.value='#E6E6E6')
# Program dot plot by tissue. No expression-based positivity threshold.
s<-read_csv(file.path(tt,'dotplot_summary.csv'),show_col_types=FALSE)%>%filter(Level=='Fine',Method=='UCell')
d<-expand_grid(tissue=c('Adjacent','Tumor'),CellType=fine,ProgramID=orderp)%>%left_join(s,by=c('tissue','CellType','ProgramID'))
d$CellType<-factor(d$CellType,levels=fine);d$ProgramID<-factor(d$ProgramID,levels=rev(orderp));d$tissue<-factor(d$tissue,levels=c('Adjacent','Tumor'))
write_csv(d,file.path(tt,'plot_program_dotplot.csv'))
p<-ggplot(d,aes(CellType,ProgramID))+geom_tile(fill='#FAFAFA',colour='white',linewidth=.3)+geom_point(data=filter(d,!is.na(Score)),aes(fill=Score,size=NPatients),shape=21,colour='#555555',stroke=.2)+geom_text(data=filter(d,is.na(Score)),label='–',size=5/ggplot2::.pt,colour='#AAAAAA')+
 facet_grid(~tissue)+fill_score()+scale_size_area(max_size=3.4,breaks=c(1,3,6,9,12),name='Patients')+scale_y_discrete(labels=labels)+scale_x_discrete(labels=flabel,drop=FALSE)+
 labs(title='Cellular sources of Figure2 programs',subtitle='Adjacent and tumor tissues; patient-weighted program scores',x=NULL,y=NULL,caption='Colour: UCell score; size: patients. Dash: unsupported. † Descriptive phenotypes. Compare within each program.')+
 theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=55,hjust=1),panel.spacing=grid::unit(4,'mm'))
emit(p,'Figure3_program_dotplot',174,129,'plot_program_dotplot.csv')
# Cohort-separated sensitivity retains all 17 programs.
cs<-read_csv(file.path(tt,'cohort_summary.csv'),show_col_types=FALSE)%>%filter(Level=='Fine',Method=='UCell',tissue=='Tumor')
q<-expand_grid(cohort=c('GSE151530','GSE138709'),CellType=fine,ProgramID=orderp)%>%left_join(cs,by=c('cohort','CellType','ProgramID'));q$CellType<-factor(q$CellType,levels=fine);q$ProgramID<-factor(q$ProgramID,levels=rev(orderp));q$cohort<-factor(q$cohort,levels=c('GSE151530','GSE138709'))
p<-ggplot(q,aes(CellType,ProgramID))+geom_tile(fill='#FAFAFA',colour='white')+geom_point(data=filter(q,!is.na(Score)),aes(fill=Score,size=NPatients),shape=21,colour='#555555',stroke=.2)+facet_grid(~cohort)+fill_score()+scale_size_area(max_size=3.4,name='Patients')+scale_y_discrete(labels=labels)+scale_x_discrete(labels=flabel,drop=FALSE)+labs(title='Tumor program sources by cohort',subtitle='Cohort-specific check of the tissue-faceted main panel',x=NULL,y=NULL)+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=55,hjust=1))
write_csv(q,file.path(tt,'plot_cohort_dotplot.csv'));emit(p,'SFigure3_cohort_dotplot',174,126,'plot_cohort_dotplot.csv','supplement')
# Same coordinates, unsmoothed native scores; per-program scale shared across tissues.
f<-read_csv(file.path(tt,'featureplot_cells.csv.gz'),show_col_types=FALSE);f$tissue<-factor(f$tissue,levels=c('Adjacent','Tumor'))
feature<-function(pr){
 dd<-f%>%mutate(Value=.data[[pr]])%>%arrange(Value);lim<-range(dd$Value,na.rm=TRUE)
 ggplot(dd,aes(UMAP1,UMAP2))+geom_point(colour='#DEDEDE',size=.12,shape=16)+geom_point(data=filter(dd,!is.na(Value)),aes(colour=Value),size=.12,shape=16)+
 facet_grid(tissue~.)+scale_colour_gradientn(colours=c('#E7E7E7','#8FC4C8','#317F8C','#123E50'),limits=lim,name='UCell')+coord_equal()+
 labs(title=gsub(' | ','\n',labels[pr],fixed=TRUE))+theme_void(base_family='Arial',base_size=6)+theme(plot.title=element_text(size=7,face='bold',hjust=.5),strip.text=element_text(size=6.5,face='bold'),strip.background=element_blank(),legend.position='bottom',legend.title=element_text(size=5.5),legend.text=element_text(size=5),plot.margin=margin(2,2,2,2))+
 guides(colour=guide_colourbar(barwidth=grid::unit(25,'mm'),barheight=grid::unit(2,'mm'),title.position='top'))
}
for(block in 1:2){
 ps<-if(block==1)c('Pattern_P2','Pattern_P3','WGCNA_MEcyan')else c('Pattern_P4','WGCNA_MEmagenta','WGCNA_MEdarkturquoise')
 p<-wrap_plots(lapply(ps,feature),nrow=1)+plot_annotation(title=if(block==1)'Metabolic and epithelial-associated program localization' else 'Epithelial development, cell cycle and antiviral programs',subtitle='Shared integrated UMAP; identical colour scale across tissues for each program',caption='All atlas cells retained. Grey includes low-support libraries or insufficient gene coverage; scores are not smoothed.',theme=theme(plot.title=element_text(size=10,face='bold'),plot.subtitle=element_text(size=6.5),plot.caption=element_text(size=5.5,hjust=0)))
 emit(p,paste0('Figure3_program_featureplots_',if(block==1)'C1_C2'else'C3'),174,120,'featureplot_cells.csv.gz')
}
# Patient-matched source lines, with each C group on one row.
pairs<-read_csv(file.path(tt,'source_patient_pairs.csv'),show_col_types=FALSE)%>%filter(Method=='UCell')
sp<-tibble(ProgramID=c('Pattern_P2','WGCNA_MElightgreen','Pattern_P3','WGCNA_MEcyan','Pattern_P4','WGCNA_MEmagenta','WGCNA_MEdarkturquoise'),Other=c('Hepatocyte','Hepatocyte','Myeloid','Myeloid','Myeloid','T_NK','Myeloid'),tissue=c('Adjacent','Adjacent',rep('Tumor',5)))
pairplots<-list();pairdata<-list()
for(i in seq_len(nrow(sp))){
 rr<-sp[i,];d<-pairs%>%filter(ProgramID==rr$ProgramID,Other==rr$Other,tissue==rr$tissue)
 d<-bind_rows(transmute(d,PatientKey,patient,cohort,CellType='Epithelial',Score=EpithelialScore),transmute(d,PatientKey,patient,cohort,CellType=Other,Score=OtherScore));d$CellType<-factor(d$CellType,levels=c('Epithelial',rr$Other));d$ProgramID<-rr$ProgramID;d$tissue<-rr$tissue;pairdata[[i]]<-d
 title<-paste0(gsub(' | ','\n',labels[rr$ProgramID],fixed=TRUE),'\n',rr$tissue,'; n = ',n_distinct(d$PatientKey))
 pairplots[[i]]<-ggplot(d,aes(CellType,Score,group=PatientKey))+geom_line(colour='#B7B7B7',linewidth=.35)+geom_point(aes(colour=CellType,shape=cohort),size=1.8,stroke=.5,show.legend=TRUE)+scale_colour_manual(values=cellcols,guide='none')+scale_shape_manual(values=c(GSE151530=16,GSE138709=17),limits=c('GSE151530','GSE138709'),drop=FALSE,name='Cohort',guide=guide_legend(override.aes=list(colour='#333333')))+scale_x_discrete(labels=c(Epithelial='Epithelial',Hepatocyte='Hepatocyte',Myeloid='Myeloid',T_NK='T / NK'))+labs(title=title,x=NULL,y='UCell score')+theme(plot.title=element_text(size=7.2,face='bold'),axis.text.x=element_text(size=6.5))
}
write_csv(bind_rows(pairdata),file.path(tt,'plot_patient_source_pairs.csv'))
p<-((pairplots[[1]]+pairplots[[2]]) / (pairplots[[3]]+pairplots[[4]]) / (pairplots[[5]]+pairplots[[6]]+pairplots[[7]]))+plot_layout(guides='collect')+plot_annotation(title='Patient-matched sources of selected Figure2 programs',subtitle='Each line connects two cell compartments from the same patient and tissue',caption='P2/lightgreen use adjacent tissue; the other five contrasts use tumor tissue. Cohort-specific tests and all contrasts are in Source Data.')&theme(legend.position='bottom')
emit(p,'Figure3_patient_source_pairs',174,143,'plot_patient_source_pairs.csv')
# Sample program landscape; retain every tumor library, including unsupported columns.
lib<-read_csv(file.path(tt,'library_program_scores.csv'),show_col_types=FALSE);cells<-read_csv(file.path(tt,'cells.csv.gz'),show_col_types=FALSE)
samples<-cells%>%filter(tissue=='Tumor')%>%distinct(cohort,patient,PatientKey,sample)%>%arrange(desc(cohort),patient,sample)
samplelabels<-setNames(gsub('ICC_','ICC',gsub('_Tumor',' / T',samples$sample)),samples$sample)
q<-expand_grid(sample=samples$sample,CellType=c('Epithelial','Myeloid','T_NK'),ProgramID=selected)%>%left_join(samples,by='sample')%>%left_join(lib%>%filter(Level=='Clean',Method=='UCell',tissue=='Tumor')%>%select(sample,CellType,ProgramID,Score,Eligible,NCells),by=c('sample','CellType','ProgramID'))%>%mutate(Value=ifelse(is.na(Eligible)|!Eligible,NA,Score))
q$sample<-factor(q$sample,levels=samples$sample);q$ProgramID<-factor(q$ProgramID,levels=rev(selected));q$CellType<-factor(q$CellType,levels=c('Epithelial','Myeloid','T_NK'));q$cohort<-factor(q$cohort,levels=c('GSE151530','GSE138709'))
p<-ggplot(q,aes(sample,ProgramID))+geom_tile(aes(fill=Value),colour='white',linewidth=.25)+geom_text(data=filter(q,is.na(Value)),label='–',size=5/ggplot2::.pt,colour='#999999')+facet_grid(CellType~cohort,scales='free_x',space='free_x',labeller=labeller(CellType=blabel))+fill_score()+scale_y_discrete(labels=labels)+scale_x_discrete(labels=samplelabels,drop=TRUE)+labs(title='Program profiles across tumor samples',subtitle='19 tumor libraries from 16 patients; repeated regions remain separate columns',x=NULL,y=NULL,caption='Uncentered UCell scores; no positive/active threshold. Repeated regions remain visible; patient-level inference groups them.')+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=.5),panel.spacing.y=grid::unit(5,'mm'))
write_csv(q,file.path(tt,'plot_sample_program_profiles.csv'));emit(p,'Figure3_sample_program_profiles',174,133,'plot_sample_program_profiles.csv')
# Full epithelial programs: all17 retained, no selected-pattern-only conclusions.
q<-expand_grid(sample=samples$sample,ProgramID=orderp)%>%left_join(samples,by='sample')%>%left_join(lib%>%filter(Level=='Clean',Method=='UCell',tissue=='Tumor',CellType=='Epithelial')%>%select(sample,ProgramID,Score,Eligible,NCells),by=c('sample','ProgramID'))%>%mutate(Value=ifelse(is.na(Eligible)|!Eligible,NA,Score))
q$sample<-factor(q$sample,levels=samples$sample);q$ProgramID<-factor(q$ProgramID,levels=rev(orderp));q$cohort<-factor(q$cohort,levels=c('GSE151530','GSE138709'))
p<-ggplot(q,aes(sample,ProgramID))+geom_tile(aes(fill=Value),colour='white',linewidth=.25)+geom_text(data=filter(q,is.na(Value)),label='–',size=5/ggplot2::.pt,colour='#999999')+facet_grid(~cohort,scales='free_x',space='free_x')+fill_score()+scale_y_discrete(labels=labels)+scale_x_discrete(labels=samplelabels,drop=TRUE)+labs(title='All Figure2 programs in tumor-tissue epithelium',subtitle='Clean epithelial identity; every tumor library retained as a column',x=NULL,y=NULL,caption='Grey: insufficient cells, gene coverage or low-support library. Epithelial identity does not establish malignancy.')+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=.5))
write_csv(q,file.path(tt,'plot_all_epithelial_programs.csv'));emit(p,'SFigure3_all_epithelial_programs',174,119,'plot_all_epithelial_programs.csv','supplement')
# Six prespecified patient associations: raw scores, cohort shapes, adjusted rank coefficients.
a<-read_csv(file.path(tt,'patient_associations.csv'),show_col_types=FALSE)%>%filter(Method=='UCell');pts<-read_csv(file.path(tt,'patient_association_points.csv'),show_col_types=FALSE)%>%filter(Method=='UCell');spec<-read_csv(file.path(tt,'association_specification.csv'),show_col_types=FALSE)
plist<-list()
for(i in seq_len(nrow(spec))){
 n<-spec$Association[i];rr<-a%>%filter(Association==n);d<-pts%>%filter(Association==n);fraction<-rr$YProgram=='Fraction'
 yl<-if(fraction)paste0(blabel[rr$YCellType],' captured fraction')else paste0(blabel[rr$YCellType],' | ',sub('C[123] \\| ','',labels[rr$YProgram]))
 if(is.na(yl))yl<-paste0(rr$YCellType,' captured fraction')
 ttl<-paste0(gsub(' | ','\n',labels[rr$XProgram],fixed=TRUE),' / ',if(fraction)rr$YCellType else if(rr$YCellType=='Myeloid')'Myeloid IFN' else 'T/NK cycle')
 subtitle<-if(is.na(rr$AdjustedR))paste0('n=',rr$NPatients,'; association not estimable')else sprintf('n=%d; adjusted r=%.2f; FDR=%.2f',rr$NPatients,rr$AdjustedR,rr$FDR)
 plist[[i]]<-ggplot(d,aes(X,Y))+geom_point(aes(shape=cohort),size=2,colour='#315F73')+scale_shape_manual(values=c(GSE151530=16,GSE138709=17),limits=c('GSE151530','GSE138709'),drop=FALSE,name='Cohort',guide=guide_legend(override.aes=list(colour='#333333')))+labs(title=ttl,subtitle=subtitle,x='Epithelial UCell score',y=yl)+theme(plot.title=element_text(size=6.8,face='bold'),plot.subtitle=element_text(size=6),axis.title=element_text(size=6))
}
p<-wrap_plots(plist,ncol=3,guides='collect')+plot_annotation(title='Patient-level epithelial–TME associations',subtitle='Independent patients; rank correlation adjusted for cohort',caption='Exploratory. Permutations preserve cohort; BH across six planned associations. Captured fractions are not absolute infiltration.')&theme(legend.position='bottom')
emit(p,'Figure3_epithelial_TME_associations',174,132,'patient_association_points.csv')
# Same-cell coexpression: descriptive correlations, no cell-level significance tests.
cc<-read_csv(file.path(tt,'within_epithelial_program_correlations.csv'),show_col_types=FALSE);cc$Pair<-paste(sub('Pattern_','',cc$Program1),sub('Pattern_','',cc$Program2),sep=' / ');cc$cohort<-factor(cc$cohort,levels=c('GSE151530','GSE138709')); cp<-read_csv(file.path(tt,'patient_program_scores.csv'),show_col_types=FALSE)%>%filter(Level=='Clean',Method=='UCell',tissue=='Tumor',CellType=='Epithelial')%>%distinct(PatientKey,patient,cohort); cc<-expand_grid(PatientKey=cp$PatientKey,Pair=c('P2 / P3','P2 / P4','P3 / P4'))%>%left_join(cp,by='PatientKey')%>%left_join(select(cc,-patient,-cohort),by=c('PatientKey','Pair'));cc$cohort<-factor(cc$cohort,levels=c('GSE151530','GSE138709'));write_csv(cc,file.path(tt,'plot_within_epithelial_correlations.csv'))
p<-ggplot(cc,aes(patient,Pair))+geom_tile(aes(fill=AdjustedR),colour='white',linewidth=.3)+geom_text(aes(label=ifelse(is.na(AdjustedR),'–',sprintf('%.2f',AdjustedR))),size=5.5/ggplot2::.pt,family='Arial')+facet_grid(~cohort,scales='free_x',space='free_x')+scale_fill_gradient2(low='#4575A5',mid='white',high='#BB615A',limits=c(-1,1),na.value='#E6E6E6',name='Adjusted r')+labs(title='Within-epithelial program coexpression',subtitle='Per-patient cell correlations adjusted for library, detected genes and mitochondrial fraction',x=NULL,y=NULL,caption='Disjoint Pattern genes; ≥50 clean cells per patient/pair. Grey: unsupported. Descriptive correlations, not triple activation.')+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=45,hjust=1))
emit(p,'SFigure3_within_epithelial_coexpression',174,65,'plot_within_epithelial_correlations.csv','supplement')
write_csv(bind_rows(exports),file.path(tt,'export_dimensions.csv'));writeLines(capture.output(sessionInfo()),file.path(o,'review/R_session_info.txt'))
