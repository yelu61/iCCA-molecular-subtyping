# Final-size typography revision for the panels used on Figure1/SFigure1 slides.
# Display-only: reuse locked data and stored fits; do not overwrite the v2 assets.
if(!exists('STEP8_ROOT'))STEP8_ROOT<-normalizePath('.')
project_root<-STEP8_ROOT
source(file.path(STEP8_ROOT,'8-Visualization/scripts/fig_style.R'))
source(file.path(STEP8_ROOT,'8-Visualization/scripts/stat_annotate.R'))
source(file.path(STEP8_ROOT,'8-Visualization/scripts/render_bulk_figures.R'))
theme_set(theme_step8())
source(file.path(STEP8_ROOT,'8-Visualization/scripts/render_fig1_discovery_v2.R'))
suppressPackageStartupMessages({library(ggplot2);library(grid);library(data.table)})
TYPE_DIR<-file.path(STEP8_ROOT,'8-Visualization/refinements/fig1_sfig1_typography_20261030')
OLD_DIR<-file.path(STEP8_ROOT,'8-Visualization/refinements/fig1_sfig1_20261026_v2')
for(d in c('pdf','review','source_data'))dir.create(file.path(TYPE_DIR,d),recursive=TRUE,showWarnings=FALSE)
refine_data<-function(id,d,suffix='') {
 path<-file.path(TYPE_DIR,'source_data',paste0(V2_ID,suffix,'.csv'));fwrite(d,path,na='NA')
 V2_DATA<<-unique(c(V2_DATA,path));invisible(d)
}
sz<-function(x)x/.pt
wrap<-function(x,width)paste(strwrap(x,width=width),collapse='\n')
type_theme<-function(title=8,axis=7,tick=6.5,note=6,legend=6.5) {
 refine_theme()+theme(
  plot.title=element_text(size=title,face='bold',lineheight=1,margin=margin(b=1)),
  plot.subtitle=element_text(size=note,colour='grey40',lineheight=1,margin=margin(b=1)),
  axis.title=element_text(size=axis),axis.text=element_text(size=tick,colour='grey25'),
  strip.text=element_text(size=7,face='bold',margin=margin(.5,.5,.5,.5)),
  legend.title=element_text(size=legend),legend.text=element_text(size=legend),
  legend.key.size=grid::unit(2.5,'mm'),legend.key.height=grid::unit(2.5,'mm'),
  legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,0,0),
  legend.box.spacing=grid::unit(.5,'mm'),panel.spacing=grid::unit(1.3,'mm'),
  plot.margin=margin(1,1,1,1,'mm'))
}
text_size<-function(p,points=6.5) {
 for(i in seq_along(p$layers))if(inherits(p$layers[[i]]$geom,c('GeomText','GeomLabel'))) {
  p$layers[[i]]$aes_params$size<-sz(points);p$layers[[i]]$geom_params$size<-sz(points)
 }
 p
}
point_size<-function(p,points) {
 for(i in seq_along(p$layers))if(inherits(p$layers[[i]]$geom,'GeomPoint'))p$layers[[i]]$aes_params$size<-points
 p
}
short_umap_legend<-function(p,id) {
 d<-p$data;method<-c(SFig1K1='NMF',SFig1K2='KMeans')[id]
 if(id=='Figure1C') {
  key<-levels(droplevels(d$Group));p<-p+scale_colour_manual(values=load_subtype_palette()$colors,
   breaks=key,labels=key,drop=TRUE)+guides(colour=guide_legend(nrow=1,override.aes=list(size=1.3,alpha=1)))
 } else if(id=='SFig2D') {
  key<-STATE_LEGEND_ORDER[STATE_LEGEND_ORDER%in%levels(droplevels(d$Group))]
  p<-p+scale_colour_manual(values=load_subtype_palette()$colors,breaks=key,labels=key,drop=TRUE)+
   guides(colour=guide_legend(nrow=2,ncol=3,override.aes=list(size=1.1,alpha=1)))
 } else if(id=='SFig1K5') {
  key<-levels(droplevels(d$Group));p<-p+scale_colour_manual(values=load_study_palette()$colors,
   breaks=key,labels=key,drop=TRUE)+guides(colour=guide_legend(nrow=4,ncol=2,override.aes=list(size=1.1,alpha=1)))
 } else if(!is.na(method)) {
  key<-levels(droplevels(d$Group));counts<-table(factor(d$Group,levels=key))
  p<-p+scale_colour_manual(values=METHOD_COLORS[[method]],breaks=key,
   labels=paste0(key,' (',as.integer(counts),')'),drop=TRUE)+
   guides(colour=guide_legend(nrow=1,override.aes=list(size=1.1,alpha=1)))
 } else {
  key<-levels(droplevels(d$Group));counts<-table(factor(d$Group,levels=key))
  p<-p+scale_colour_manual(values=setNames(c('grey35','grey75'),key),breaks=key,
   labels=paste0(key,' (',as.integer(counts),')'),drop=TRUE)+
   guides(colour=guide_legend(nrow=1,override.aes=list(size=1.1,alpha=1)))
 }
 p
}
compact_taxonomy<-function(id) {
 original<-v2_taxonomy(id);e<-environment(original);d<-e$d;key<-e$key
 nr<-e$nr;row_levels<-e$row_levels;pal<-e$pal;colors<-e$colors
 function() {
  grid.newpage()
  txt<-function(s,x,y,size=6.5,bold=FALSE,color='black',just='centre')
   grid.text(s,x,y,just=just,gp=gpar(fontfamily='Arial',fontsize=size,lineheight=.95,
    fontface=if(bold)'bold' else 'plain',col=color))
  box<-function(x,y,w,h,fill)grid.rect(x,y,w,h,gp=gpar(fill=fill,col='white',lwd=.8))
  title<-wrap(key$Title,31);two_line_title<-grepl('\n',title,fixed=TRUE)
  txt(title,.5,.985,7.5,TRUE,just='top')
  header_y<-if(two_line_title).825 else .85
  x<-c(.35,.53,.71);w<-.17;top<-if(two_line_title).75 else .79;bottom<-.19;rh<-(top-bottom)/nr
  for(j in 1:3){box(x[j],header_y,w,.065,pal[j]);txt(paste0('C',j),x[j],header_y,7,TRUE,'grey15')}
  txt('Original',.015,header_y,6.1,TRUE,just='left');txt('Total',.90,header_y,6.1,TRUE)
  for(i in seq_len(nr)) {
   y<-top-(i-.5)*rh;lab<-row_levels[i];ss<-d[OriginalLabel==lab]
   txt(if(lab=='Not reported')'Not\nreported' else sub('class: ','',lab),.125,y,if(lab=='Not reported')5.8 else 6.5,TRUE)
   for(j in 1:3){r<-ss[OriginalLabel==lab&FinalSubtype==paste0('C',j)]
    fill<-if(lab=='Not reported')'grey94' else colors[1+round(100*r$RowFraction)]
    box(x[j],y,w,rh,fill);txt(sprintf('%d\n(%.0f%%)',r$N,100*r$RowFraction),x[j],y,6,FALSE,
     if(lab!='Not reported'&&r$RowFraction>.7)'white' else 'black')}
   txt(sprintf('%d\n(%.0f%%)',ss$RowN[1],100*ss$RowN[1]/unique(d$TotalN)),.90,y,6)
  }
  txt('Column\ntotal',.015,.125,6,TRUE,just='left')
  for(j in 1:3){n<-d[FinalSubtype==paste0('C',j),sum(N)];box(x[j],.125,w,.08,'grey96')
   txt(sprintf('%d\n(%.0f%%)',n,100*n/sum(d$N)),x[j],.125,6,TRUE)}
  txt(paste0('N = ',sum(d$N)),.90,.125,6.1,TRUE)
  txt(paste0('n (row %); labels ',unique(d$AvailableLabelsN),'/',unique(d$TotalN)),.5,.035,5.8)
 }
}
compact_clinical<-function() {
 id<-'Figure1F_compact';path<-file.path(OLD_DIR,'source_data',paste0(id,'.csv'));d<-fread(path)
 baseline<-fread(file.path(STEP8_ROOT,T6,'step6_clinical_baseline_by_main_state.csv'))
 for(i in seq_len(nrow(d))){r<-d[i];owner<-baseline[Variable==r$Variable&Level==r$Level&FinalSubtype==r$FinalSubtype]
  stopifnot(nrow(owner)==1L,owner$AvailableN==r$AvailableN,owner$Summary==r$Summary)}
 d[,FinalSubtype:=factor(FinalSubtype,levels=c('C1','C2','C3'))];setorder(d,TraitOrder,FinalSubtype)
 newpath<-file.path(TYPE_DIR,'source_data',paste0(id,'.csv'));fwrite(d,newpath,na='NA');V2_DATA<<-newpath
 v2_track(c(sub(paste0(STEP8_ROOT,'/'),'',path,fixed=TRUE),paste0(T6,'step6_clinical_baseline_by_main_state.csv')))
 function() {
  grid.newpage();txt<-function(label,x,y,size=6.5,bold=FALSE,col='grey20',just='centre')
   grid.text(label,x,y,just=just,gp=gpar(fontfamily='Arial',fontsize=size,
    fontface=if(bold)'bold' else 'plain',col=col))
  pal<-load_subtype_palette()$colors;percent_col<-circlize::colorRamp2(c(0,50,100),c('#F8FAFC','#D5E7F0','#246085'))
  x<-c(.405,.635,.865);ys<-c(.735,.555,.375,.195);w<-.22;h<-.16
  txt('Clinical characteristics',.025,.965,8.5,TRUE,just='left')
  for(j in 1:3){grid.rect(x[j],.87,w,.075,gp=gpar(fill=pal[paste0('C',j)],col=NA))
   txt(paste0('C',j),x[j],.87,8,TRUE,'grey15')}
  row_titles<-c('Tumor size','Stage III/IV','Vascular invasion','Large-duct type')
  row_notes<-c('Median [IQR], cm','Observed proportion','Observed proportion','Observed pathology')
  for(i in 1:4){ss<-d[TraitOrder==i];yy<-ys[i]
   txt(row_titles[i],.025,yy+.035,7.2,TRUE,just='left');txt(row_notes[i],.025,yy-.005,6,just='left')
   txt(paste0('Available n = ',sum(ss$AvailableN)),.025,yy-.045,6,col='grey40',just='left')
   for(j in 1:3){r<-ss[FinalSubtype==paste0('C',j)];fill<-if(i==1)'grey96' else percent_col(r$Value)
    grid.rect(x[j],yy,w,h,gp=gpar(fill=fill,col=NA));txt(r$Primary,x[j],yy+.027,8,TRUE)
    txt(r$Secondary,x[j],yy-.015,6.5);if(i==1)txt(r$Tertiary,x[j],yy-.054,6,col='grey40')}}
  txt('Proportion (%)',.025,.05,6,just='left')
  for(k in 1:101)grid.rect(.305+(k-.5)*.35/101,.055,.35/101,.017,gp=gpar(fill=percent_col(k-1),col=NA))
  txt('0',.305,.025,5.8);txt('50',.48,.025,5.8);txt('100',.655,.025,5.8)
  txt('n/N: positive / available',.98,.05,6,just='right')
 }
}
compact_km<-function() {
 stopifnot(requireNamespace('survminer',quietly=TRUE),requireNamespace('survival',quietly=TRUE))
 owner<-'6-ClinicalSurvival/review/20261026_figure1_h_options'
 fit<-readRDS(file.path(STEP8_ROOT,owner,'km_fit_loglog.rds'))
 d<-fread(file.path(STEP8_ROOT,owner,'tables/survival_population.csv'))
 d[,FinalSubtype:=factor(FinalSubtype,levels=c('C1','C2','C3'))]
 med<-fread(file.path(STEP8_ROOT,T6,'step6_km_median_survival.csv'));tab<-summary(fit)$table
 for(st in c('C1','C2','C3')){old<-med[FinalSubtype==st];rr<-tab[paste0('FinalSubtype=',st),]
  stopifnot(old$N==rr['records'],old$Events==rr['events'],isTRUE(all.equal(old$MedianMonths,unname(rr['median']),tolerance=1e-8)))}
 lr<-fread(file.path(STEP8_ROOT,T6,'step6_km_logrank_tests.csv'))
 pval<-lr[Test=='Study-stratified log-rank',PValue]
 plab<-if(pval<.001)'P < 0.001' else paste0('P = ',signif(pval,2))
 v2_track(c(paste0(owner,c('/km_fit_loglog.rds','/tables/survival_population.csv')),
  paste0(T6,c('step6_km_median_survival.csv','step6_km_logrank_tests.csv'))))
 file.copy(file.path(STEP8_ROOT,owner,'tables/survival_population.csv'),
  file.path(TYPE_DIR,'source_data','Figure1G_KM_population.csv'),overwrite=TRUE)
 V2_DATA<<-file.path(TYPE_DIR,'source_data','Figure1G_KM_population.csv')
 km<-survminer::ggsurvplot(fit,data=d,risk.table=TRUE,conf.int=FALSE,
  palette=unname(load_subtype_palette()$colors[c('C1','C2','C3')]),legend.title='',legend.labs=c('C1','C2','C3'),
  legend='top',xlab='Overall survival (months)',ylab='Survival probability',break.time.by=12,
  risk.table.height=.29,risk.table.fontsize=2.1,censor.shape=124,censor.size=2,
  ggtheme=theme_classic(base_size=6.5,base_family='Arial'),
  tables.theme=survminer::theme_cleantable(base_size=6.2,base_family='Arial'),
  title='Overall survival',
  subtitle=sprintf('Available-case KM | n = 331; events = 153\nStudy-stratified log-rank %s',plab))
 km$plot<-km$plot+type_theme(title=8,axis=6.5,tick=6.2,note=6,legend=6.5)+
  theme(legend.position='top',legend.direction='horizontal',plot.margin=margin(.8,.8,.5,.8,'mm'))
 km$table<-km$table+theme(text=element_text(family='Arial',size=6.2),
  axis.text=element_text(size=6.2),axis.title=element_text(size=6.2),
  plot.title=element_text(size=6.2),plot.margin=margin(.5,.5,.5,.5,'mm'))
 grid::grid.grabExpr(print(km))
}
compact_rmst<-function() {
 path<-'6-ClinicalSurvival/review/20261026_figure1_h_options/tables/rmst_24_differences.csv'
 d<-fread(file.path(STEP8_ROOT,path));stopifnot(nrow(d)==3L)
 newpath<-file.path(TYPE_DIR,'source_data','Figure1H_RMST24_difference.csv')
 file.copy(file.path(STEP8_ROOT,path),newpath,overwrite=TRUE);V2_DATA<<-newpath
 d<-copy(d);d[,Row:=c(3,2,1)];d[,Label:=sprintf('%+.2f  [%.2f, %.2f]',DifferenceMonths,Lower95,Upper95)]
 v2_track(c(path,paste0(T6,'step6_clinical_survival_release_gate_status.csv')))
 ggplot(d,aes(DifferenceMonths,Row))+geom_vline(xintercept=0,colour='grey65',linetype=2,linewidth=.3)+
  geom_segment(aes(x=Lower95,xend=Upper95,yend=Row),linewidth=.6,colour='grey45')+
  geom_point(shape=21,size=2,fill='#426782',colour='white',stroke=.35)+
  geom_text(aes(y=Row+.25,label=Label),size=sz(6),family='Arial')+
  scale_y_continuous(breaks=c(3,2,1),labels=c('C2 − C1','C3 − C1','C3 − C2'),limits=c(.65,3.65),expand=expansion(mult=0))+
  scale_x_continuous(limits=c(-12,9),breaks=c(-10,-5,0,5),expand=expansion(mult=.01))+
  labs(title='24-month RMST differences',subtitle='Difference [95% CI] | n = 331; events = 153',
   x='RMST difference (months)',y=NULL)+
  type_theme(title=8,axis=6.5,tick=6.2,note=6)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())
}
TYPE_SIZES<-list(
 Figure1B=c(123.3,88.6),Figure1C=c(58.1,49.3),Figure1D=c(53.7,39.9),
 Figure1E1=c(43.4,38),Figure1E2=c(43.4,38),Figure1E3=c(43.4,38),Figure1E4=c(43.4,38),Figure1E5=c(43.4,38),
 Figure1F_compact=c(82.7,47.9),Figure1G_KM=c(51.7,44.6),Figure1H_RMST24_difference=c(43.4,29.7),
 SFig1A=c(93,56),SFig1B=c(78,56),SFig1C=c(90,49.1),SFig1K5=c(44,44),SFig1E=c(97.7,62.8),
 SFig1F=c(48,42),SFig1K1=c(42,42),SFig1K2=c(42,42),SFig1K3=c(42,42),SFig1K4=c(42,42),
 SFig1G=c(50.2,44),SFig2D=c(45.7,44),SFig2E=c(73,44))
type_plot<-function(id) {
 if(id%in%paste0('Figure1E',1:5))return(compact_taxonomy(id))
 if(id=='Figure1F_compact')return(compact_clinical())
 if(id=='Figure1G_KM')return(compact_km())
 if(id=='Figure1H_RMST24_difference')return(compact_rmst())
 p<-v2_plot(id)
 if(id=='Figure1B') {
  body(p)<-parse(text=gsub('show_legend\\s*=\\s*FALSE\\)',
   'show_legend=FALSE,show_annotation_name=FALSE)',paste(deparse(body(p)),collapse='\n')))[[1]]
  return(p)
 }
 if(id=='Figure1C'){p<-short_umap_legend(p,id);p$labels$title<-wrap(p$labels$title,32)
  return(point_size(p,1)+type_theme(title=8,axis=7,tick=6.5,legend=6.5))}
 if(id=='Figure1D'){p$labels$title<-'Principal states across studies'
  return(text_size(p,6.3)+type_theme(title=8,axis=6.5,tick=6.2,legend=6.3)+
   scale_y_continuous(labels=scales::label_percent(),breaks=c(0,.25,.5,.75,1),expand=expansion(mult=c(0,.04)))+
   theme(axis.text.y=element_text(size=6),plot.margin=margin(.7,2,.7,.7,'mm')))}
 if(id%in%c('SFig1K1','SFig1K2','SFig1K3','SFig1K4','SFig1K5','SFig2D')) {
  p<-short_umap_legend(p,id)
  short_titles<-c(SFig1K1='NMF assignments',SFig1K2='k-means assignments',
   SFig1K3='ConsensusCluster (excluded)',SFig1K4='Hierarchical (excluded)',
   SFig1K5='Discovery UMAP by study',SFig2D='Principal and intermediate states')
  p$labels$title<-wrap(short_titles[[id]],28)
  return(point_size(p,.8)+type_theme(title=7.2,axis=6.5,tick=6,legend=6)+
   theme(plot.margin=margin(.7,.7,.7,.7,'mm')))
 }
 if(id=='SFig1A')return(text_size(p,6.3)+type_theme(title=8,axis=7,tick=6.3,legend=6.3))
 if(id=='SFig1B')return(text_size(p,5.9)+type_theme(title=8,axis=6.5,tick=5.9,legend=6.2)+
  scale_x_continuous(labels=function(x)ifelse(x>=1000,paste0(x/1000,'k'),as.character(x)),expand=expansion(mult=c(0,.42)))+
  theme(axis.text.y=element_text(size=5.9),panel.spacing=grid::unit(2.2,'mm')))
 if(id=='SFig1C')return(p+type_theme(title=8,axis=6.5,tick=6,legend=6)+
  guides(colour=guide_legend(nrow=2,override.aes=list(size=1.2,alpha=1)),shape=guide_legend(nrow=2)))
 if(id=='SFig1E') {
  for(i in 1:6)p[[i]]$labels$title<-wrap(p[[i]]$labels$title,23)
  p[[6]]<-p[[6]]+scale_y_continuous(expand=expansion(mult=c(0,.25)))+
   theme(plot.title=element_text(hjust=0))
  p[[6]]$labels$title<-'Hierarchical\nstored partition'
  p<-p&theme(plot.title=element_text(size=6.8,lineheight=1),plot.subtitle=element_text(size=6,lineheight=1),
   axis.title=element_text(size=6.2),axis.text=element_text(size=6),plot.margin=margin(.7,.7,.7,.7,'mm'),
   panel.spacing=grid::unit(1,'mm'))
  return(p+patchwork::plot_annotation(title='Cluster-number diagnostics and method selection',
   theme=theme(plot.title=element_text(size=8.5,face='bold',family='Arial'))))
 }
 if(id=='SFig1F')return(text_size(p,6)+type_theme(title=8,axis=6.5,tick=6,legend=6.2)+
  scale_x_discrete(labels=c(ConsensusCluster='Cons.',NMF='NMF',KMeans='k-means',HC='HC'))+
  scale_y_discrete(labels=c(ConsensusCluster='Cons.',NMF='NMF',KMeans='k-means',HC='HC'))+
  guides(fill=guide_colourbar(barwidth=grid::unit(25,'mm'),barheight=grid::unit(2.5,'mm'),title.position='top'))+
  labs(title='Method agreement (ARI)')+theme(axis.text.x=element_text(angle=30,hjust=1,size=5.8)))
 if(id=='SFig1G')return(text_size(p,6.3)+type_theme(title=8,axis=6.5,tick=6,legend=6.2)+
  scale_x_discrete(labels=c('1\n(C2)','2\n(C3)','3\n(C1)'))+
  scale_y_discrete(labels=c('1 (C1)','2 (C2)','3 (C3)'))+
  labs(title='NMF–k-means state crosswalk',x='k-means cluster (aligned state)',y='NMF cluster (aligned state)'))
 if(id=='SFig2E')return(text_size(p,5.9)+type_theme(title=8,axis=6.5,tick=6,legend=6)+
  scale_x_continuous(labels=scales::label_percent(),limits=c(0,1),expand=expansion(mult=c(0,.03)))+
  guides(fill=guide_legend(nrow=2,ncol=3))+theme(axis.text.y=element_text(size=5.9),plot.margin=margin(.7,2,.7,.7,'mm')))
 stop('No typography renderer for ',id)
}
render_type<-function(id) {
 check_owner_gate('Step0');check_owner_gate('Step1');if(id%in%c('Figure1F_compact','Figure1G_KM','Figure1H_RMST24_difference'))check_owner_gate('Step6')
 V2_ID<<-id;V2_DATA<<-character();refine_sources<<-c(
  '8-Visualization/scripts/render_fig1_sfig1_typography_20261030.R','8-Visualization/scripts/render_fig1_recovery_panels.R',
  '8-Visualization/scripts/render_fig1_discovery_v2.R','8-Visualization/scripts/render_fig1_refinement.R',
  '8-Visualization/scripts/fig_style.R','8-Visualization/scripts/stat_annotate.R','8-Visualization/scripts/render_bulk_figures.R',
  '8-Visualization/scripts/figure_config.json',SUBTYPE_PALETTE_CSV,STUDY_PALETTE_CSV,ASSAY_PALETTE_CSV,
  '0-Data/qc/step0_release_gate_checks.csv',paste0(T1,'step1_clustering_release_gate_status.csv'))
 p<-type_plot(id);size<-TYPE_SIZES[[id]];stopifnot(!is.null(size))
 pdf<-file.path(TYPE_DIR,'pdf',paste0(id,'.pdf'));save_panel(p,pdf,size[1],size[2])
 png<-file.path(TYPE_DIR,'review',paste0(id,'.png'))
 grDevices::png(png,width=size[1],height=size[2],units='mm',res=300,type='cairo',family='Arial')
 tryCatch({if(is.function(p))p() else if(inherits(p,'ggsurvplot'))print(p) else if(inherits(p,'grob'))grid.draw(p) else print(p)},finally=dev.off())
 data_names<-basename(V2_DATA);data_hash<-unname(tools::md5sum(V2_DATA));unchanged<-TRUE
 if(id!='Figure1G_KM')for(i in seq_along(V2_DATA)) {
  old<-if(id=='Figure1H_RMST24_difference')file.path(OLD_DIR,'h_options_20261026/source_data',paste0(id,'.csv')) else file.path(OLD_DIR,'source_data',data_names[i])
  if(file.exists(old))unchanged<-unchanged&&identical(unname(tools::md5sum(old)),data_hash[i])
 }
 stopifnot(unchanged)
 inputs<-sort(unique(refine_sources));hashes<-unname(tools::md5sum(file.path(STEP8_ROOT,inputs)));stopifnot(!anyNA(hashes))
 row<-data.table(Panel=id,Status='Rendered_typography_review',WidthMM=size[1],HeightMM=size[2],
  PDF=sub(paste0(STEP8_ROOT,'/'),'',pdf,fixed=TRUE),PDFMD5=unname(tools::md5sum(pdf)),
  SourcePaths=paste(inputs,collapse=';'),SourceMD5=paste(hashes,collapse=';'),
  SourceData=paste(data_names,collapse=';'),SourceDataMD5=paste(data_hash,collapse=';'),SourceDataUnchanged=unchanged)
 path<-file.path(TYPE_DIR,'figure_manifest.csv');old<-if(file.exists(path))fread(path) else data.table()
 if(nrow(old))old<-old[Panel!=id];m<-rbindlist(list(old,row),fill=TRUE);setorder(m,Panel);fwrite(m,path)
 message('Rendered final-size typography: ',id);invisible(pdf)
}
render_type_group<-function(ids=names(TYPE_SIZES))for(id in ids)render_type(id)
if(sys.nframe()==0)render_type_group()
