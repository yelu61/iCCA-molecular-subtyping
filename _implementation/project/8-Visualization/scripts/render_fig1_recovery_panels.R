# Focused completion of the interrupted Figure1 revision; no model refitting.
v2_method_selection <- function(id) {
  d<-refine_read(paste0(T1,'step1_method_selection_metrics_long.csv'))
  disposition<-refine_read(paste0(T1,'step1_method_disposition.csv'))
  parameters<-refine_read(paste0(T1,'step1_clustering_methods_parameters.csv'))
  summary<-refine_read('1-Clustering/data/clustering_methods_summary.csv')
  v2_track('archive/notebook_backups/1-clustering.ipynb')
  path<-'1-Clustering/dataclustering_results.RData';v2_track(path)
  e<-new.env();load(file.path(STEP8_ROOT,path),envir=e)
  old<-as.data.table(e$comparison_data)
  stopifnot(nrow(old)==586L,!anyDuplicated(old$Sample))
  counts<-old[,.(N=.N),by=HC][order(-N)]
  stopifnot(identical(counts$N,c(585L,1L)))
  refine_data(id,d);refine_data(id,disposition,'_disposition')
  refine_data(id,parameters,'_parameters');refine_data(id,summary,'_stored_k')
  refine_data(id,counts,'_hc_counts')
  metric<-function(method,metric,title,subtitle,chosen=NULL) {
    x<-d[Method==method & Metric==metric]
    stopifnot(nrow(x)==5L)
    p<-ggplot(x,aes(K,Value))+geom_line(colour='grey40',linewidth=.4)+
      geom_point(size=1.6,colour='grey35')+scale_x_continuous(breaks=2:6)+
      labs(title=title,subtitle=subtitle,x='Number of clusters (k)',y=NULL)+refine_theme()+
      theme(plot.title=element_text(size=8),plot.subtitle=element_text(size=7),
        axis.text=element_text(size=7),axis.title=element_text(size=7),plot.margin=margin(3,3,3,3,'mm'))
    if(!is.null(chosen))p<-p+geom_vline(xintercept=chosen,linetype=2,colour='grey70',linewidth=.3)+
      geom_point(data=x[K==chosen],shape=21,fill='white',size=2.8,stroke=.5)
    p
  }
  ps<-list(
    metric('NMF','Cophenetic (higher is better)','NMF · cophenetic','Maximum at k = 3',3),
    metric('NMF','PAC (lower is better)','NMF · PAC','Minimum at k = 3',3),
    metric('KMeans','Mean silhouette (higher is better)','k-means · silhouette','Maximum at k = 3',3),
    metric('Consensus','PAC (lower is better)','Consensus · standard PAC','Minimum at k = 2',2),
    metric('Consensus','Largest cluster fraction','Consensus · largest cluster','99–100% across tested k')+scale_y_continuous(limits=c(0,1.03),breaks=c(0,.5,1)),
    ggplot(counts,aes(factor(HC),N))+geom_col(width=.5,fill='grey65')+
      geom_text(aes(label=paste0('n = ',N)),vjust=-.4,size=2.5)+
      scale_y_continuous(expand=expansion(mult=c(0,.15)))+
      labs(title='Hierarchical · stored partition',subtitle='Complete linkage; k = 2',x='Stored cluster',y='Samples')+
      refine_theme()+theme(plot.title=element_text(size=8),plot.subtitle=element_text(size=7),
        axis.text=element_text(size=7),axis.title=element_text(size=7),plot.margin=margin(3,3,3,3,'mm')))
  patchwork::wrap_plots(ps,ncol=3)+patchwork::plot_annotation(
    title='Cluster-number diagnostics and method selection',
    theme=theme(plot.title=element_text(size=10,face='bold',family='Arial')))
}

v2_cox_timevarying <- function(id) {
  check_owner_gate('Step6')
  v2_track(c(paste0(T6,'step6_clinical_survival_release_gate_status.csv'),
    '6-ClinicalSurvival/scripts/run_step6_clinical_survival.R'))
  d<-refine_read(paste0(T6,'step6_cox_time_varying_subtype_effects.csv'))
  ph<-refine_read(paste0(T6,'step6_cox_proportional_hazards_diagnostics.csv'))
  stopifnot(nrow(d)==18L,all(d$N==331L),all(d$Events==153L),
    all(d$Lower95>0),all(d$Lower95<=d$HR),all(d$Upper95>=d$HR),
    !anyDuplicated(d[,.(Contrast,TimeMonths)]))
  refine_data(id,d);refine_data(id,ph,'_ph_diagnostics')
  d[,Contrast:=factor(Contrast,levels=c('C2 vs C1','C3 vs C1','C3 vs C2'))]
  d[,Time:=factor(TimeMonths,levels=rev(c(3,6,12,24,36,60)))]
  ggplot(d,aes(HR,Time))+geom_vline(xintercept=1,linetype=2,colour='grey65',linewidth=.35)+
    geom_segment(aes(x=Lower95,xend=Upper95,yend=Time),colour='grey45',linewidth=.5)+
    geom_point(shape=21,size=2.2,fill='#426782',colour='white',stroke=.3)+
    facet_wrap(~Contrast,nrow=1)+scale_x_log10(limits=c(.0625,64),breaks=c(.1,1,10),labels=c('0.1','1','10'))+
    labs(title='Adjusted overall-survival associations',subtitle='Time-varying Cox  |  n = 331; events = 153',
      x='Hazard ratio (95% CI; log scale)',y='Follow-up (months)')+refine_theme()+
    theme(plot.title=element_text(size=9),plot.subtitle=element_text(size=7),
      axis.text=element_text(size=7),axis.title=element_text(size=7.5),strip.text=element_text(size=7.5),
      panel.spacing=grid::unit(3,'mm'),plot.margin=margin(3,3,3,3,'mm'))
}

v2_state_counts <- function(id) {
  m<-refine_master(FALSE);d<-m[,.N,by=.(Study,FinalSubtype)]
  d<-merge(CJ(Study=sort(unique(m$Study)),FinalSubtype=STATE_ORDER),d,all.x=TRUE);d[is.na(N),N:=0L]
  d[,FinalSubtype:=factor(FinalSubtype,levels=STATE_ORDER)];setorder(d,Study,FinalSubtype)
  d[,Total:=sum(N),by=Study];d[,Fraction:=N/Total];d[,Midpoint:=cumsum(Fraction)-Fraction/2,by=Study]
  d[,Y:=match(Study,rev(sort(unique(Study))))]
  d[,LabelPlacement:=ifelse(N==0,'zero_no_segment',ifelse(Fraction<.04,'above_with_leader','inside'))]
  refine_data(id,d)
  totals<-unique(d[,.(Study,Total,Y)]);pal<-load_subtype_palette()$colors
  ggplot(d)+geom_rect(aes(xmin=Midpoint-Fraction/2,xmax=Midpoint+Fraction/2,
      ymin=Y-.28,ymax=Y+.28,fill=FinalSubtype),colour='white',linewidth=.2)+
    geom_text(data=d[LabelPlacement=='inside'],aes(Midpoint,Y,label=N,
      colour=FinalSubtype=='C1'),size=2.5)+
    geom_segment(data=d[LabelPlacement=='above_with_leader'],
      aes(x=Midpoint,xend=Midpoint,y=Y+.28,yend=Y+.39),colour='grey35',linewidth=.25)+
    geom_text(data=d[LabelPlacement=='above_with_leader'],aes(Midpoint,Y+.49,label=N),size=2.4)+
    scale_colour_manual(values=c('TRUE'='white','FALSE'='grey10'),guide='none')+
    scale_fill_manual(values=pal,breaks=STATE_LEGEND_ORDER)+
    scale_x_continuous(labels=scales::label_percent(),limits=c(0,1),expand=expansion(mult=c(0,.01)))+
    scale_y_continuous(breaks=totals$Y,labels=paste0(totals$Study,'  (n = ',totals$Total,')'),
      limits=c(.5,7.7),expand=expansion(mult=0))+
    labs(title='All discovery states across studies',x='Within-study proportion',y=NULL,fill=NULL)+
    guides(fill=guide_legend(nrow=1))+refine_theme()+
    theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(size=8))
}
