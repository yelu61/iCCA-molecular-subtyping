# Revised discovery narrative. Source tables and patient labels are read-only.
# Load the preceding renderer for shared styling, gates and validated primitives.
source(file.path(STEP8_ROOT,'8-Visualization/scripts/render_fig1_refinement.R'))
REFINE_DIR <- file.path(STEP8_ROOT,'8-Visualization/refinements/fig1_sfig1_20261026_v2')
for(d in c('pdf','source_data','review'))dir.create(file.path(REFINE_DIR,d),recursive=TRUE,showWarnings=FALSE)
V2_ID <- '';V2_DATA <- character();V2_CACHE<-new.env(parent=emptyenv())
STATE_ORDER<-c('C1','C1_C2','C1_C3','C2','C2_C3','C3')
STATE_LEGEND_ORDER<-c('C1','C2','C3','C1_C2','C1_C3','C2_C3')
# Algorithm IDs are categorical labels, not aliases of the governed state palette.
METHOD_COLORS<-list(NMF=c('1'='#D9D9D9','2'='#969696','3'='#404040'),
  KMeans=c('1'='#E0DBEB','2'='#A99AC4','3'='#66517F'))
AVAILABILITY_COLORS<-list(Genome=c(Yes='#466E65',No='#DAEAE6'),
  Protein=c(Yes='#856C4B',No='#EAE0D2'))
refine_data <- function(id,d,suffix='') {
  path<-file.path(REFINE_DIR,'source_data',paste0(V2_ID,suffix,'.csv'))
  fwrite(d,path,na='NA');V2_DATA<<-unique(c(V2_DATA,path));invisible(d)
}
v2_track <- function(paths)refine_sources<<-unique(c(refine_sources,paths))
v2_top_data <- function() {
  paths<-c(SUBTYPE_TRUTH_CSV,'1-Clustering/data/top_variable_genes.csv','1-Clustering/data/combined_expr_after_correction.rds')
  v2_track(paths)
  if(exists('top',V2_CACHE,inherits=FALSE))return(V2_CACHE$top)
  m<-refine_master(FALSE);top<-refine_read(paths[2]);stopifnot(nrow(top)==2000L,!anyDuplicated(top$Gene))
  e<-new.env();load(file.path(STEP8_ROOT,paths[3]),envir=e);expr<-e$corrected_expr_limma
  stopifnot(all(top$Gene %in% rownames(expr)),all(m$Sample %in% colnames(expr)))
  m[,FinalSubtype:=factor(FinalSubtype,levels=STATE_ORDER)];setorder(m,FinalSubtype,Study,Sample)
  raw<-expr[top$Gene,m$Sample,drop=FALSE];stopifnot(all(is.finite(raw)),all(apply(raw,1,sd)>0))
  mu<-rowMeans(raw);sig<-apply(raw,1,sd);z<-sweep(sweep(raw,1,mu,'-'),1,sig,'/')
  gene_order<-hclust(dist(z),method='ward.D2')$order
  sample_order<-unlist(lapply(STATE_ORDER,function(st) {
    ix<-which(m$FinalSubtype==st)
    if(length(ix)<2L)ix else ix[hclust(dist(t(z[,ix,drop=FALSE])),method='ward.D2')$order]
  }),use.names=FALSE)
  out<-list(z=z[gene_order,sample_order],m=m[sample_order],
    rows=data.table(Gene=rownames(z)[gene_order],Center=mu[gene_order],SD=sig[gene_order],RowOrder=seq_along(gene_order)))
  V2_CACHE$top<-out;out
}
v2_heatmap <- function(id,markers=FALSE,main=FALSE) {
  if(!markers) {
    x<-v2_top_data();z<-x$z;m<-copy(x$m);rows<-copy(x$rows);lim<-2
  } else {
    m<-refine_master(FALSE);m[,FinalSubtype:=factor(FinalSubtype,levels=STATE_ORDER)];setorder(m,FinalSubtype,Study,Sample)
    rows<-refine_read(paste0(T1,'step1_bulk_principal_state_markers_for_step8.csv'))[!duplicated(Gene)]
    path<-'1-Clustering/data/combined_expr_after_correction.rds';v2_track(path)
    e<-new.env();load(file.path(STEP8_ROOT,path),envir=e)
    raw<-e$corrected_expr_limma[rows$Gene,m$Sample,drop=FALSE];ref<-m$SubtypeClass=='Main'
    mu<-rowMeans(raw[,ref]);sig<-apply(raw[,ref],1,sd);stopifnot(all(sig>0),all(is.finite(raw)))
    z<-sweep(sweep(raw,1,mu,'-'),1,sig,'/');rows[,`:=`(Center=mu,SD=sig,RowOrder=seq_len(.N))];lim<-3
  }
  if(main){keep<-m$SubtypeClass=='Main';m<-m[keep];z<-z[,keep,drop=FALSE]}
  states<-STATE_ORDER[STATE_ORDER %in% as.character(m$FinalSubtype)]
  m[,FinalSubtype:=factor(FinalSubtype,levels=states)]
  refine_data(id,data.table(Gene=rownames(z),z),'_row_zscore')
  stopifnot(all(m[['Genome sequencing']] %in% c('YES','NO')),
    all(m[['Protein profiling']] %in% c('YES','NO')))
  m[,GenomeData:=ifelse(get('Genome sequencing')=='YES','Yes','No')]
  m[,ProteinData:=ifelse(get('Protein profiling')=='YES','Yes','No')]
  refine_data(id,rows,'_rows');refine_data(id,m[,.(Sample,Study,FinalSubtype,SubtypeClass,NMF,KMeans,GenomeData,ProteinData)],'_columns')
  refine_data(id,rbindlist(list(m[,.(N=.N),by=.(Availability=GenomeData)][,Track:='Genome data'],
    m[,.(N=.N),by=.(Availability=ProteinData)][,Track:='Protein data'])),'_annotation_counts')
  counts<-table(m$FinalSubtype);stopifnot(sum(counts)==if(main)497L else 586L)
  if(markers)rows[,LabelShown:=seq_len(.N) %in% unique(round(seq(1,.N,length.out=8))),by=Subtype]
  function() {
    pal<-load_subtype_palette()$colors
    state_keys<-STATE_LEGEND_ORDER[STATE_LEGEND_ORDER %in% states]
    study_colors<-load_study_palette()$colors
    study_keys<-names(study_colors)[names(study_colors) %in% m$Study]
    ann<-ComplexHeatmap::HeatmapAnnotation(State=as.character(m$FinalSubtype),
      NMF=paste0('N',m$NMF),KMeans=paste0('K',m$KMeans),Study=m$Study,
      `Genome data`=m$GenomeData,`Protein data`=m$ProteinData,
      col=list(State=pal,NMF=setNames(METHOD_COLORS$NMF,paste0('N',1:3)),
        KMeans=setNames(METHOD_COLORS$KMeans,paste0('K',1:3)),Study=study_colors,
        `Genome data`=AVAILABILITY_COLORS$Genome,`Protein data`=AVAILABILITY_COLORS$Protein),
      simple_anno_size=grid::unit(2,'mm'),annotation_name_gp=grid::gpar(fontsize=7),
      show_legend=FALSE)
    discrete_legend<-function(title,keys,colors,nrow=1,ncol=NULL,by_row=TRUE)ComplexHeatmap::Legend(title=title,at=keys,
      legend_gp=grid::gpar(fill=unname(colors[keys])),nrow=nrow,ncol=ncol,by_row=by_row,
      grid_height=grid::unit(2.6,'mm'),grid_width=grid::unit(2.6,'mm'),
      title_gp=grid::gpar(fontsize=7),labels_gp=grid::gpar(fontsize=7))
    effect_col<-circlize::colorRamp2(c(-lim,0,lim),c(SCALE_EFFECT_LOW,SCALE_EFFECT_MID,SCALE_EFFECT_HIGH))
    method_legends<-ComplexHeatmap::packLegend(
      discrete_legend('NMF',paste0('N',1:3),setNames(METHOD_COLORS$NMF,paste0('N',1:3))),
      discrete_legend('KMeans',paste0('K',1:3),setNames(METHOD_COLORS$KMeans,paste0('K',1:3))),
      ComplexHeatmap::Legend(title='Row z-score',col_fun=effect_col,at=c(-lim,0,lim),direction='horizontal',
        legend_width=grid::unit(18,'mm'),title_gp=grid::gpar(fontsize=7),labels_gp=grid::gpar(fontsize=7)),
      direction='horizontal',gap=grid::unit(6,'mm'))
    legends<-ComplexHeatmap::packLegend(discrete_legend('State',state_keys,pal),
      discrete_legend('Study',study_keys,study_colors),method_legends,
      ComplexHeatmap::packLegend(discrete_legend('Genome data',c('Yes','No'),AVAILABILITY_COLORS$Genome),
        discrete_legend('Protein data',c('Yes','No'),AVAILABILITY_COLORS$Protein),
        direction='horizontal',gap=grid::unit(6,'mm')),
      direction='vertical',gap=grid::unit(2,'mm'))
    right_legends<-!markers
    if(right_legends) {
      method_column<-ComplexHeatmap::packLegend(
        discrete_legend('NMF',paste0('N',1:3),setNames(METHOD_COLORS$NMF,paste0('N',1:3)),nrow=3),
        discrete_legend('KMeans',paste0('K',1:3),setNames(METHOD_COLORS$KMeans,paste0('K',1:3)),nrow=3),
        direction='horizontal',gap=grid::unit(4,'mm'))
      availability_column<-ComplexHeatmap::packLegend(
        discrete_legend('Genome data',c('Yes','No'),AVAILABILITY_COLORS$Genome,nrow=2),
        discrete_legend('Protein data',c('Yes','No'),AVAILABILITY_COLORS$Protein,nrow=2),
        direction='horizontal',gap=grid::unit(4,'mm'))
      legends<-ComplexHeatmap::packLegend(
        discrete_legend('State',state_keys,pal,nrow=3,ncol=if(main)1 else 2,by_row=FALSE),
        discrete_legend('Study',study_keys,study_colors,nrow=length(study_keys)),
        method_column,availability_column,
        ComplexHeatmap::Legend(title='Row z-score',col_fun=effect_col,at=c(-lim,0,lim),direction='horizontal',
          legend_width=grid::unit(20,'mm'),title_gp=grid::gpar(fontsize=7),labels_gp=grid::gpar(fontsize=7)),
        direction='vertical',gap=grid::unit(3,'mm'))
    }
    h<-ComplexHeatmap::Heatmap(z,name='Row z-score',col=circlize::colorRamp2(c(-lim,0,lim),c(SCALE_EFFECT_LOW,SCALE_EFFECT_MID,SCALE_EFFECT_HIGH)),
      cluster_rows=FALSE,cluster_columns=FALSE,cluster_row_slices=FALSE,cluster_column_slices=FALSE,
      row_split=if(markers)factor(rows$Subtype,levels=c('C1','C2','C3')) else NULL,
      column_split=m$FinalSubtype,column_title=NULL,
      column_title_gp=grid::gpar(fontsize=7,fontface='bold'),row_title_gp=grid::gpar(fontsize=8),
      column_gap=grid::unit(.45,'mm'),row_gap=grid::unit(1.5,'mm'),
      show_column_names=FALSE,show_row_names=markers,
      row_labels=if(markers)ifelse(rows$LabelShown,rows$Gene,'') else rownames(z),row_names_gp=grid::gpar(fontsize=7),
      top_annotation=ann,use_raster=!markers,raster_device='agg_png',raster_quality=6,raster_by_magick=FALSE,raster_resize_mat=FALSE,
      heatmap_legend_param=list(at=c(-lim,0,lim),direction='horizontal',title_position='topcenter',title_gp=grid::gpar(fontsize=7),labels_gp=grid::gpar(fontsize=7)))
    ComplexHeatmap::draw(h,show_heatmap_legend=FALSE,annotation_legend_side=if(right_legends)'right' else 'bottom',
      annotation_legend_list=list(legends),show_annotation_legend=TRUE,padding=grid::unit(c(4,4,4,4),'mm'))
  }
}
v2_umap_data <- function(main=FALSE) {
  u<-refine_read('governance/data/umap_coordinates_top2000.csv');u[,Sample:=as.character(Sample)]
  m<-refine_master(FALSE);d<-merge(m,u[,.(Sample,UMAP1,UMAP2)],by='Sample',all.x=TRUE)
  stopifnot(nrow(d)==586L,all(is.finite(d$UMAP1)),all(is.finite(d$UMAP2)))
  if(main)d<-d[SubtypeClass=='Main'];d
}
v2_umap <- function(id,main=FALSE,method=NULL) {
  all<-v2_umap_data(FALSE);d<-if(main)all[SubtypeClass=='Main'] else copy(all)
  pal<-load_subtype_palette()$colors;labels<-load_subtype_palette()$labels
  if(is.null(method)) {
    d[,Group:=factor(FinalSubtype,levels=STATE_ORDER[STATE_ORDER %in% as.character(d$FinalSubtype)])];colors<-pal
    title<-if(main)'Principal states in the discovery UMAP' else 'Principal and intermediate discovery states'
    display<-setNames(paste0(names(table(d$Group)),' (n = ',table(d$Group),')'),names(table(d$Group)))
    legend_shapes<-ifelse(names(table(d$Group)) %in% c('C1','C2','C3'),16,1)
    p<-ggplot(d,aes(UMAP1,UMAP2,color=Group,shape=SubtypeClass))+geom_point(size=1.5,alpha=.85,stroke=.65)+
      scale_shape_manual(values=c(Main=16,Intermediate=1),guide='none')+
      scale_colour_manual(values=colors,labels=display,drop=TRUE)
  } else {
    legend_shapes<-16
    if(method %in% c('ConsensusCluster','HC')) {
      path<-'1-Clustering/dataclustering_results.RData';v2_track(path);e<-new.env();load(file.path(STEP8_ROOT,path),envir=e)
      old<-as.data.table(e$comparison_data);old[,Sample:=as.character(Sample)]
      stopifnot(!anyDuplicated(old$Sample),setequal(old$Sample,d$Sample))
      d[,Group:=factor(old[[method]][match(Sample,old$Sample)])]
      colors<-setNames(c('grey35','grey75'),levels(d$Group));title<-paste0(method,': stored excluded partition')
    } else if(method=='Study') {
      d[,Group:=factor(Study)];colors<-load_study_palette()$colors;title<-'Discovery UMAP by study'
    } else {
      d[,Group:=factor(get(method))]
      colors<-METHOD_COLORS[[method]]
      title<-paste0(method,' assignments on the same UMAP')
    }
    stopifnot(!anyNA(d$Group));display<-setNames(paste0(levels(d$Group),' (n = ',table(d$Group),')'),levels(d$Group))
    p<-ggplot(d,aes(UMAP1,UMAP2,color=Group))+geom_point(size=1.5,alpha=.8)+scale_colour_manual(values=colors,labels=display)
  }
  refine_data(id,d[,.(Sample,Study,FinalSubtype,SubtypeClass,NMF,KMeans,UMAP1,UMAP2,Group)])
  # Equal spans make the panel square while preserving equal coordinate units.
  span<-max(diff(range(all$UMAP1)),diff(range(all$UMAP2)))
  xlim<-mean(range(all$UMAP1))+c(-.5,.5)*span
  ylim<-mean(range(all$UMAP2))+c(-.5,.5)*span
  p+coord_fixed(ratio=1,xlim=xlim,ylim=ylim)+labs(title=title,x='UMAP 1',y='UMAP 2',color=NULL)+
    refine_theme()+theme(panel.border=element_rect(colour='grey30',fill=NA,linewidth=.45),axis.line=element_blank())+
    guides(colour=guide_legend(nrow=if(is.null(method) || method=='Study')3 else 1,override.aes=list(alpha=1,size=2,shape=legend_shapes)))
}
TAXONOMIES<-data.table(ID=paste0('Figure1E',1:5),Study=c('2019Ahn','2022Dong','2022Dong','2023Deng','2023Jeon'),
  Field=c('Subclass','Subclass','Subclass-add','Subclass','Subclass'),
  Title=c('2019Ahn · original A/B classes','2022Dong · proteomic S1–S4','2022Dong · transcriptomic IG1–IG3',
    '2023Deng · proteomic S-I–S-III','2023Jeon · original LC classes'))
v2_taxonomy <- function(id) {
  key<-TAXONOMIES[ID==id];m<-refine_master(TRUE);m<-m[Study==key$Study]
  m[,OriginalLabel:=as.character(get(key$Field))];m[is.na(OriginalLabel)|OriginalLabel %in% c('NA','w/o',''),OriginalLabel:='Not reported']
  known<-sort(setdiff(unique(m$OriginalLabel),'Not reported'));row_levels<-c(known,if('Not reported' %in% m$OriginalLabel)'Not reported')
  d<-m[,.N,by=.(OriginalLabel,FinalSubtype)]
  d<-merge(CJ(OriginalLabel=row_levels,FinalSubtype=c('C1','C2','C3')),d,all.x=TRUE);d[is.na(N),N:=0L]
  d[,RowN:=sum(N),by=OriginalLabel];d[,RowFraction:=N/RowN];d[,ColumnN:=sum(N),by=FinalSubtype]
  d[,`:=`(Study=key$Study,SourceField=key$Field,TotalN=nrow(m),AvailableLabelsN=sum(m$OriginalLabel!='Not reported'))]
  d[,RowOrder:=match(OriginalLabel,row_levels)];setorder(d,RowOrder,FinalSubtype)
  stopifnot(sum(d$N)==nrow(m),all(abs(d[,sum(RowFraction),by=OriginalLabel]$V1-1)<1e-10))
  refine_data(id,d);refine_data(id,m[,.(Sample,Study,FinalSubtype,OriginalLabel)],'_samples')
  nr<-length(row_levels);pal<-load_subtype_palette()$colors;colors<-grDevices::colorRampPalette(c('white','lightyellow',unname(pal['C1'])))(101)
  function() {
    grid::grid.newpage()
    txt<-function(s,x,y,size=8,bold=FALSE,color='black',just='centre')grid::grid.text(s,x,y,just=just,gp=grid::gpar(fontfamily='Arial',fontsize=size,fontface=if(bold)'bold' else 'plain',col=color))
    box<-function(x,y,w,h,fill)grid::grid.rect(x,y,w,h,gp=grid::gpar(fill=fill,col='white',lwd=1.2))
    txt(key$Title,.5,.96,10,TRUE);txt('Correspondence with C1 / C2 / C3',.5,.91,8)
    x<-c(.37,.55,.73);w<-.18;top<-.79;bottom<-.27;rh<-(top-bottom)/nr
    for(j in 1:3){box(x[j],.835,w,.07,pal[j]);txt(paste0('C',j),x[j],.835,9,TRUE,'grey15')}
    txt('Original class',.035,.835,8,TRUE,just='left');txt('Row total',.91,.85,7,TRUE);txt('n (% total)',.91,.815,7)
    for(i in seq_len(nr)) {
      y<-top-(i-.5)*rh;lab<-row_levels[i];ss<-d[OriginalLabel==lab]
      txt(if(lab=='Not reported')'Not\nreported' else sub('class: ','',lab),.14,y,8,TRUE)
      for(j in 1:3){r<-ss[FinalSubtype==paste0('C',j)];fill<-if(lab=='Not reported')'grey94' else colors[1+round(100*r$RowFraction)]
        box(x[j],y,w,rh,fill);txt(sprintf('%d (%.0f%%)',r$N,100*r$RowFraction),x[j],y,8,FALSE,if(lab!='Not reported'&&r$RowFraction>.7)'white' else 'black')}
      txt(sprintf('%d\n(%.0f%%)',ss$RowN[1],100*ss$RowN[1]/nrow(m)),.91,y,8)
    }
    txt('Column total',.14,.215,7,TRUE)
    for(j in 1:3){n<-d[FinalSubtype==paste0('C',j),sum(N)];box(x[j],.215,w,.07,'grey96');txt(sprintf('%d (%.0f%%)',n,100*n/nrow(m)),x[j],.215,8,TRUE)}
    txt(paste0('N = ',nrow(m)),.91,.215,8,TRUE)
    txt('Cell: n (% within original class)',.5,.145,7)
    txt(paste0('Original label available: ',sum(m$OriginalLabel!='Not reported'),' / ',nrow(m)),.5,.105,7)
    for(k in 1:101)grid::grid.rect(.37+(k-.5)*.36/101,.055,.36/101,.015,gp=grid::gpar(fill=colors[k],col=NA))
    txt('Row %',.28,.055,7);txt('0',.37,.025,7);txt('50',.55,.025,7);txt('100',.73,.025,7)
  }
}
CLINICAL_PANELS<-data.table(ID=paste0('Figure1F',1:4),
  Variable=c('Tumor size (cm)','Stage group','Vascular invasion','Observed duct type'),
  Level=c('Median (IQR)','III/IV','1','large duct type'),
  Endpoint=c('Tumor size','Advanced stage','Vascular invasion','Observed large-duct morphology'),
  Title=c('Tumor size','Stage III/IV','Vascular invasion','Observed large-duct type'),
  ExpectedN=c(312L,341L,248L,93L))
v2_clinical <- function(id) {
  check_owner_gate('Step6');v2_track(paste0(T6,'step6_clinical_survival_release_gate_status.csv'))
  key<-CLINICAL_PANELS[ID==id]
  baseline<-refine_read(paste0(T6,'step6_clinical_baseline_by_main_state.csv'))
  availability<-refine_read(paste0(T6,'step6_clinical_phenotype_availability.csv'))
  d<-baseline[Variable==key$Variable & Level==key$Level]
  stopifnot(nrow(d)==3L,setequal(d$FinalSubtype,c('C1','C2','C3')),!anyDuplicated(d$FinalSubtype),
    sum(d$AvailableN)==key$ExpectedN,
    availability[Endpoint==key$Endpoint & Study=='All studies',AvailableN]==key$ExpectedN)
  d[,FinalSubtype:=factor(FinalSubtype,levels=c('C1','C2','C3'))];setorder(d,FinalSubtype)
  d[,Endpoint:=key$Endpoint]
  continuous<-id=='Figure1F1'
  if(continuous) {
    values<-lapply(strsplit(gsub('[^0-9.]',' ',d$Summary),' +'),function(x)as.numeric(x[nzchar(x)]))
    stopifnot(all(lengths(values)==3L))
    d[,`:=`(Median=vapply(values,`[`,numeric(1),1),Q1=vapply(values,`[`,numeric(1),2),Q3=vapply(values,`[`,numeric(1),3))]
    stopifnot(all(is.finite(d$Median)),all(d$Q1<=d$Median),all(d$Median<=d$Q3),all(d$Q1>=0))
    d[,`:=`(Value=Median,LabelY=Q3+.45,ValueLabel=sprintf('%.1f',Median),
      AxisLabel=paste0(FinalSubtype,'\nn = ',AvailableN))]
  } else {
    stopifnot(all(is.finite(d$N)),all(d$N>=0),all(d$N<=d$AvailableN))
    d[,`:=`(PositiveN=N,Value=100*N/AvailableN)]
    d[,`:=`(LabelY=Value+7,ValueLabel=sprintf('%.1f%%',Value),
      AxisLabel=paste0(FinalSubtype,'\n',PositiveN,'/',AvailableN))]
  }
  refine_data(id,d)
  refine_data(id,availability[Endpoint==key$Endpoint],'_availability')
  pal<-load_subtype_palette()$colors
  p<-ggplot(d,aes(FinalSubtype,Value))
  if(continuous)p<-p+geom_errorbar(aes(ymin=Q1,ymax=Q3,colour=FinalSubtype),width=.14,linewidth=.65)
  else p<-p+geom_segment(aes(xend=FinalSubtype,y=0,yend=Value),colour='grey83',linewidth=.45)
  p<-p+geom_point(aes(fill=FinalSubtype),shape=21,size=3,colour='white',stroke=.45)+
    geom_text(aes(y=LabelY,label=ValueLabel),size=2.65,colour='grey20')+
    scale_fill_manual(values=pal,guide='none')+scale_colour_manual(values=pal,guide='none')+
    scale_x_discrete(labels=setNames(d$AxisLabel,as.character(d$FinalSubtype)),expand=expansion(add=.55))+
    labs(title=key$Title,subtitle=paste0(if(continuous)'Median and IQR' else 'Observed proportion','  |  n = ',key$ExpectedN),
      x=NULL,y=if(continuous)'Tumor size (cm)' else 'Patients (%)')+
    refine_theme()+theme(plot.title=element_text(size=9,face='bold'),plot.subtitle=element_text(size=7,colour='grey40'),
      axis.text=element_text(size=8,colour='grey25'),axis.title=element_text(size=8),
      panel.grid.major.y=element_line(colour='grey93',linewidth=.25),panel.grid.major.x=element_blank(),
      panel.grid.minor=element_blank(),legend.position='none',plot.margin=margin(4,4,3,3,'mm'))
  if(continuous)p+scale_y_continuous(limits=c(0,10),breaks=seq(0,10,2),expand=expansion(mult=c(0,.02)))
  else p+scale_y_continuous(limits=c(0,100),breaks=seq(0,100,25),expand=expansion(mult=c(0,.02)))
}
source(file.path(STEP8_ROOT,'8-Visualization/scripts/render_fig1_recovery_panels.R'))
v2_plot <- function(id) {
  if(id=='Figure1G_cox_timevarying')return(v2_cox_timevarying(id))
  if(id=='Figure1B')return(v2_heatmap(id))
  if(id=='Figure1C')return(v2_umap(id,main=TRUE))
  if(id=='Figure1D')return(refine_plot('Figure1E'))
  if(id %in% TAXONOMIES$ID)return(v2_taxonomy(id))
  if(id %in% CLINICAL_PANELS$ID)return(v2_clinical(id))
  borrow<-c(SFig1A='SFig1B',SFig1C='SFig1E',SFig1D='SFig1F')
  if(id %in% names(borrow))return(refine_plot(borrow[[id]]))
  if(id=='SFig1B') {
    d<-refine_read('0-Data/qc/step0_expression_overview.csv');refine_data(id,d)
    long<-melt(d,id.vars=c('Dataset','Study','AnalysisRole'),measure.vars=c('FinalProcessedFeatures','Samples'),variable.name='Quantity',value.name='N')
    long[,Quantity:=factor(Quantity,levels=c('FinalProcessedFeatures','Samples'),labels=c('Processed features','Expression samples'))]
    long[,Dataset:=factor(Dataset,levels=rev(unique(d$Dataset)))]
    return(ggplot(long,aes(N,Dataset,fill=AnalysisRole))+geom_col(width=.65)+geom_text(aes(label=format(N,big.mark=',')),hjust=-.15,size=2.4)+
      facet_wrap(~Quantity,scales='free_x',nrow=1)+scale_fill_manual(values=load_assay_palette()$colors)+
      scale_x_continuous(expand=expansion(mult=c(0,.22)))+labs(title='Expression resources across studies',x=NULL,y=NULL,fill=NULL)+refine_theme())
  }
  if(id=='SFig1E')return(v2_method_selection(id))
  if(id=='SFig1F') {
    d<-refine_read(paste0(T1,'step1_method_agreement_matrix_long.csv'));refine_data(id,d)
    lv<-c('ConsensusCluster','NMF','KMeans','HC');d[,MethodA:=factor(MethodA,levels=lv)];d[,MethodB:=factor(MethodB,levels=rev(lv))]
    return(ggplot(d,aes(MethodA,MethodB,fill=AdjustedRandIndex))+geom_tile(color='white')+
      geom_text(aes(label=ifelse(AdjustedRandIndex>0 & AdjustedRandIndex<.01,'<0.01',sprintf('%.2f',AdjustedRandIndex)),color=AdjustedRandIndex>.7),size=3)+
      scale_colour_manual(values=c('FALSE'='black','TRUE'='white'),guide='none')+
      scale_fill_gradient(low='white',high=load_subtype_palette()$colors['C1'],limits=c(0,1),breaks=c(0,.5,1))+
      guides(fill=guide_colourbar(barwidth=grid::unit(35,'mm'),barheight=grid::unit(3,'mm')))+
      scale_x_discrete(labels=c(ConsensusCluster='Consensus',NMF='NMF',KMeans='k-means',HC='Hierarchical'))+
      scale_y_discrete(labels=c(ConsensusCluster='Consensus',NMF='NMF',KMeans='k-means',HC='Hierarchical'))+
      coord_equal()+labs(title='Agreement between stored method assignments',x=NULL,y=NULL,fill='ARI')+refine_theme()+theme(axis.line=element_blank(),axis.ticks=element_blank()))
  }
  if(id=='SFig1G') {
    d<-refine_read(paste0(T1,'step1_nmf_kmeans_state_definition_full.csv'));refine_data(id,d)
    pal<-c(load_subtype_palette()$colors,None='grey97')
    return(ggplot(d,aes(factor(KMeans),factor(NMF),fill=Ensemble))+geom_tile(color='white',linewidth=.7)+
      geom_text(aes(label=ifelse(N==0,'n = 0',paste0(Ensemble,'\nn = ',N)),color=Ensemble=='C1'),size=3)+
      scale_colour_manual(values=c('FALSE'='black','TRUE'='white'),guide='none')+scale_fill_manual(values=pal,guide='none')+
      scale_x_discrete(labels=c('1\n(C2-aligned)','2\n(C3-aligned)','3\n(C1-aligned)'))+
      scale_y_discrete(labels=c('1 (C1-aligned)','2 (C2-aligned)','3 (C3-aligned)'))+
      coord_equal()+labs(title='NMF–k-means crosswalk defines the six states',x='k-means cluster',y='NMF cluster')+refine_theme()+theme(axis.line=element_blank(),axis.ticks=element_blank()))
  }
  if(id=='SFig1H') {
    d<-refine_read(paste0(T1,'step1_kmeans_feature_number_sensitivity.csv'));refine_data(id,d)
    return(ggplot(d,aes(Features,ARIvsLockedKMeans))+geom_hline(yintercept=1,linetype=2,color='grey75')+geom_line()+geom_point(size=2)+
      scale_x_continuous(breaks=d$Features)+scale_y_continuous(limits=c(0,1.03))+labs(title='k-means: feature-number sensitivity',x='Variable genes',y='ARI versus fixed k-means')+refine_theme())
  }
  if(id=='SFig1I') {
    d<-refine_read(paste0(T1,'step1_kmeans_k3_seed_sensitivity.csv'));refine_data(id,d)
    return(ggplot(d,aes(Seed,ARIvsHistoricalKMeans))+geom_point(size=1.2)+scale_y_continuous(limits=c(0,1.03),breaks=c(0,.5,1))+
      labs(title='k-means: reproducibility across 100 seeds',x='Seed',y='ARI versus fixed k-means')+refine_theme())
  }
  if(id=='SFig1J') {
    d<-refine_read(paste0(T1,'step1_kmeans_leave_one_study_out.csv'));refine_data(id,d)
    long<-melt(d,id.vars=c('HeldOutStudy','TrainingN','HeldOutN'),measure.vars=c('ARIAllVsHistoricalKMeans','ARITrainingVsHistoricalKMeans','ARIHeldOutVsHistoricalKMeans'),variable.name='Scope',value.name='ARI')
    long[,Scope:=factor(Scope,levels=c('ARIAllVsHistoricalKMeans','ARITrainingVsHistoricalKMeans','ARIHeldOutVsHistoricalKMeans'),labels=c('All samples','Training samples','Held-out samples'))]
    long[,StudyLabel:=paste0(HeldOutStudy,' (n = ',HeldOutN,')')]
    return(ggplot(long,aes(ARI,StudyLabel,shape=Scope))+geom_vline(xintercept=1,color='grey80',linetype=2)+
      geom_point(position=position_dodge(width=.45),size=2)+scale_shape_manual(values=c(16,1,17))+scale_x_continuous(limits=c(0,1.03))+
      labs(title='k-means: leave-one-study-out sensitivity',x='ARI versus fixed k-means',y='Omitted study (held-out N)',shape=NULL)+refine_theme())
  }
  methods<-c(SFig1K1='NMF',SFig1K2='KMeans',SFig1K3='ConsensusCluster',SFig1K4='HC',SFig1K5='Study')
  if(id %in% names(methods))return(v2_umap(id,method=methods[[id]]))
  if(id=='SFig2A')return(v2_heatmap(id,main=TRUE))
  if(id=='SFig2B')return(v2_heatmap(id,markers=TRUE))
  if(id=='SFig2C')return(v2_heatmap(id,markers=TRUE,main=TRUE))
  if(id=='SFig2D')return(v2_umap(id))
  if(id=='SFig2E')return(v2_state_counts(id))
  stop('Unknown v2 panel: ',id)
}
v2_caption <- function(id) {
  if(id=='Figure1G_cox_timevarying')return('Study-stratified time-varying Cox associations with overall survival among 331 principal-state patients with complete core covariates and 153 events. All three pairwise contrasts at all six locked follow-up times (3, 6, 12, 24, 36 and 60 months) are shown. Points are HRs and horizontal intervals are source-reported pointwise 95% Wald CIs, not simultaneous intervals. The model adjusts for age, sex and stage, stratifies by Study, and includes subtype and age interactions with log(max(time, 0.1)). The constant-HR core model violated proportional hazards (subtype P=0.002945; global P=0.002191); one constant HR therefore does not describe the full time course. Source P values and BH FDR across the 18 contrasts remain in the CSV, without significance stars. No model is refitted or time point selected by significance. This is an adjusted association display, not an independent validation or incremental prediction result. Existing KM remains a separate panel.')
  if(id %in% CLINICAL_PANELS$ID) {
    key<-CLINICAL_PANELS[ID==id]
    return(paste0(key$Title,' among principal-state tumors with observed data (available n=',key$ExpectedN,
      ' of 497 eligible main-state tumors). ',
      if(id=='Figure1F1')'Points are source-reported medians and intervals are interquartile ranges, not confidence intervals. Numeric labels show medians in cm; subtype labels show available n. Medians/quartiles retain the precision of the locked summary table. '
      else 'Points are positive counts divided by endpoint-specific available cases within each subtype. Numeric labels show percentages and subtype labels show positive n / available N. All three binary panels share a 0–100% scale. No confidence intervals are calculated. ',
      'Missing observations are excluded only from the relevant endpoint denominator, not imputed or treated as negative. These are descriptive available-case summaries from the locked Step6 baseline table, not study-adjusted estimates. No new tests, significance stars or model fits are added. Study-adjusted effects remain in the clinical analysis. Observed duct labels are never supplemented with predicted labels.'))
  }
  if(id %in% c('Figure1B','SFig2A'))return(paste0('Top 2,000 locked variable genes in ',if(id=='Figure1B')'all 586 discovery samples' else 'the 497 principal-state subset',
    '. Row z-scores use means and SDs from all 586 samples; the color scale saturates at ±2 without capping source values. Ward.D2 with Euclidean distance orders genes and samples within each fixed state for display only. The main-only panel preserves exactly the all-state row order, scaling and relative sample order. No patient clustering, labels or state counts are refitted. Historical heatmap k6/k7 refers to gene-row grouping, not patient subtype count; no new six/seven-module partition is asserted. This displays the discovery structure, not independent validation.'))
  if(id %in% c('SFig2B','SFig2C'))return('All 120 frozen corrected markers, shown across all six states (586) or principal states only (497). The same main-state reference means and SDs define both z-score matrices; color saturation ±3 is display-only. Gene order is fixed by the source table. Samples are ordered by governed state, Study and Sample; all-state and main-only views preserve relative ordering. Eight evenly spaced row labels per originating subtype are shown; no observations are dropped for readability. Discovery-selected markers do not independently validate the same labels.')
  if(id %in% c('Figure1C','SFig2D'))return('Existing top2000 UMAP coordinates, joined by Sample to the current governed labels. Main Figure1C shows only the 497 principal-state samples (C1=112, C2=223, C3=162). Supplementary SFig2D retains all 586 samples and distinguishes principal (filled) and intermediate (open) points. The principal-only display is an exact coordinate subset with the same full-cohort limits and square, equal-unit frame. No UMAP fitting or hyperparameter search is performed. Geometry is descriptive and does not establish statistical separation, trajectories or independent subtype validation.')
  if(id=='Figure1D')return(refine_caption('Figure1E'))
  if(id %in% TAXONOMIES$ID)return('Study-specific correspondence between an original classification and C1/C2/C3 among principal-state samples. Each cell shows n and the row-normalized percentage; blue intensity follows this fraction. Colored column headers encode the governed subtype colors, independently of the cell color scale. Row/column margins show counts and percentages of all principal-state samples from that study. Not-reported original labels are retained as a separate row where present, so column totals use the full principal-state denominator. Known-class row percentages do not include the missing-label row. No new association test, accuracy metric or one-to-one label equivalence is claimed.')
  borrow<-c(SFig1A='SFig1B',SFig1C='SFig1E',SFig1D='SFig1F');if(id %in% names(borrow))return(refine_caption(borrow[[id]]))
  if(id=='SFig1B')return('Locked processed-expression feature and sample counts for every dataset. Analysis-role colors distinguish discovery and projection. Two Jusakual matrix resources remain separate at the feature-resource level; these bars must not be summed as independent patient cohorts. Platform-specific feature counts are descriptive and are not a common-quality ranking.')
  if(id=='SFig1E')return('Method-specific cluster-number selection across 586 discovery samples. Top row: the locked NMF 10-run audit across k=2–6 has maximal cophenetic and minimal standard PAC at k=3; k-means mean silhouette is maximal at k=3 (0.1014). These audit diagnostics support the retained operational NMF/k-means pair; they do not replace the historical selection procedure. Historical NMF rank estimation used 20 runs, encountered automatic-selection fallback and explicitly retained user_chosen_k=3; the final fit used 50 runs. Historical k-means used the WSS elbow with a silhouette fallback; the saved solution is k=3. Bottom row: matched consensus standard PAC is smallest at k=2, but the largest cluster remains 98.98–99.83% across k=2–6 (99.49% at k=3), so it is excluded as degenerate. Historical consensus selection used a non-standard PAC-labelled score and is not relabelled as the standard audit PAC. Hierarchical clustering historically compared Euclidean-distance ward.D2, complete and average linkage using maximum silhouette across candidate k; the saved complete-linkage k=2 partition has 585 and 1 samples and is excluded. The final bar chart displays that stored partition, not a reconstructed silhouette curve. No four-method vote defines C1/C2/C3. Rings mark k=3 in retained methods and k=2 for the consensus PAC minimum; axes have separate scales. No new clustering, intervals or tests are calculated. Use this single selection asset for screenshot panel D; omit the earlier duplicate retained-method-only chart.')
  if(id=='SFig1F')return('Stored pairwise adjusted Rand indices for ConsensusCluster, NMF, k-means and hierarchical assignments. Consensus/hierarchical are degenerate stored two-cluster solutions and are excluded from the operational pair. NMF–k-means ARI is approximately 0.619. Algorithm agreement is not a resampling stability estimate; H–J show executed k-means sensitivity tests. No historical four-method ensemble is used to define current states.')
  if(id=='SFig1G')return('Locked NMF-by-k-means crosswalk across 586 discovery samples. Algorithm IDs retain their original meaning; axis annotations identify their principal-state alignment. Concordant aligned cells define C1/C2/C3, and off-diagonal disagreements map to three explicit pairwise intermediate labels. Empty cells are shown as n=0. Counts sum to 586: principal 497, intermediate 89. Intermediate labels describe method discordance and are not proven temporal transitions.')
  if(id=='SFig1H')return('Locked k-means sensitivity to 1,000, 2,000 and 5,000 variable genes. The plotted ARI compares each assignment with fixed k-means labels across the discovery cohort. No feature number is selected anew for figure appearance; no interval or test is added.')
  if(id=='SFig1I')return('All 100 executed k-means seed results from the locked source table, each with ARI 1 against fixed k-means assignments. Seeds are computational replicates, not patients. This establishes the observed k-means seed reproducibility under this configuration, not multi-seed stability of NMF or the combined state definition.')
  if(id=='SFig1J')return('Locked leave-one-study-out k-means sensitivity, showing separate ARI values for all, training and omitted-study samples. Axis labels give omitted-study sample counts; full training and held-out counts are retained in source data. This is a discovery perturbation analysis, not independently labelled external accuracy validation. No intervals or significance tests are invented.')
  if(grepl('^SFig1K',id))return('Unchanged all-sample UMAP coordinates colored separately by an algorithm assignment or Study. NMF and k-means labels come from the canonical master table; consensus/hierarchical labels come from the saved historical comparison object and are explicitly excluded stored two-cluster partitions. The historical four-method Ensemble column is not plotted. Point counts are sample counts. These displays share coordinates and axis limits; geometric appearance is not an additional cluster-selection test.')
  if(id=='SFig2E')return('All 586 discovery samples, grouped by Study and the six governed states. Within-study proportions are descriptive. Every nonzero segment has an absolute count: fractions below 4% use a leader and an above-bar label; other counts are centered inside. Study labels give all-state totals. Zero combinations remain explicitly recorded in source data and have no segment. Figure1D uses principal-only denominators, so its bars must not be compared as if denominators were identical. Intermediate states are preserved rather than reassigned to a principal state.')
  stop('Missing caption')
}
V2_GROUPS<-list(main=c('Figure1B','Figure1C','Figure1D',paste0('Figure1E',1:5),paste0('Figure1F',1:4),'Figure1G_cox_timevarying'),
  sfig1=c(paste0('SFig1',LETTERS[1:10]),paste0('SFig1K',1:5)),paired=paste0('SFig2',LETTERS[1:5]))
render_discovery_v2 <- function(id) {
  check_owner_gate('Step0');check_owner_gate('Step1');V2_ID<<-id;V2_DATA<<-character()
  refine_sources<<-c('8-Visualization/scripts/render_fig1_recovery_panels.R','8-Visualization/scripts/render_fig1_discovery_v2.R','8-Visualization/scripts/render_fig1_refinement.R',
    '8-Visualization/scripts/fig_style.R','8-Visualization/scripts/stat_annotate.R','8-Visualization/scripts/render_bulk_figures.R',
    '8-Visualization/scripts/figure_config.json',SUBTYPE_PALETTE_CSV,STUDY_PALETTE_CSV,ASSAY_PALETTE_CSV,
    '0-Data/qc/step0_release_gate_checks.csv',paste0(T1,'step1_clustering_release_gate_status.csv'))
  size<-c(174,105)
  if(id %in% c('Figure1B','SFig2A'))size<-c(174,125)
  if(id %in% c('SFig2B','SFig2C'))size<-c(174,150)
  if(id %in% c('Figure1C','SFig2D')||grepl('^SFig1K',id))size<-c(114,114)
  if(id=='Figure1D')size<-c(114,85)
  if(id %in% TAXONOMIES$ID)size<-c(114,100)
  if(id %in% CLINICAL_PANELS$ID)size<-c(74,70)
  if(id=='SFig1B')size<-c(174,125)
  if(id %in% c('SFig1F','SFig1G','SFig1H','SFig1I'))size<-c(114,100)
  if(id=='SFig1E')size<-c(174,112)
  if(id=='Figure1G_cox_timevarying')size<-c(114,78)
  plot<-v2_plot(id);pdf<-file.path(REFINE_DIR,'pdf',paste0(id,'.pdf'));save_panel(plot,pdf,size[1],size[2])
  grDevices::png(file.path(REFINE_DIR,'review',paste0(id,'.png')),width=size[1],height=size[2],units='mm',res=300,type='cairo',family='Arial')
  tryCatch({if(is.function(plot))plot() else print(plot)},finally=grDevices::dev.off())
  caption<-v2_caption(id)
  if(id %in% c('Figure1B','SFig2A','SFig2B','SFig2C'))caption<-paste(caption,
    'Genome data and Protein data are the governed master-table YES/NO availability flags (Genome sequencing and Protein profiling), joined to the displayed samples; they do not encode assay results, quality or usable downstream matched-subset size. No is the recorded negative flag, not a missing-value imputation. NMF and k-means use categorical shades within separate color families; shade intensity does not encode confidence or biological magnitude. Top state/count labels are omitted; state identities remain in the annotation legend and counts remain in source data.')
  if(id %in% c('SFig1K1','SFig1K2'))caption<-paste(caption,
    'Algorithm IDs use the same categorical shade palettes as the heatmaps; shade intensity has no quantitative or ordinal interpretation.')
  writeLines(caption,file.path(REFINE_DIR,'source_data',paste0(id,'_caption.txt')))
  inputs<-sort(unique(refine_sources));hashes<-unname(tools::md5sum(file.path(STEP8_ROOT,inputs)));stopifnot(!anyNA(hashes))
  row<-data.table(Panel=id,Status='Rendered_review_draft',WidthMM=size[1],HeightMM=size[2],PDF=sub(paste0(STEP8_ROOT,'/'),'',pdf,fixed=TRUE),
    PDFMD5=unname(tools::md5sum(pdf)),SourcePaths=paste(inputs,collapse=';'),SourceMD5=paste(hashes,collapse=';'),
    SourceData=paste(basename(V2_DATA),collapse=';'),SourceDataMD5=paste(unname(tools::md5sum(V2_DATA)),collapse=';'))
  path<-file.path(REFINE_DIR,'figure_manifest.csv');old<-if(file.exists(path))fread(path) else data.table()
  if(nrow(old))old<-old[Panel!=id];m<-rbindlist(list(old,row),fill=TRUE);setorder(m,Panel);fwrite(m,path)
  message('Rendered v2: ',id);invisible(pdf)
}
render_discovery_group <- function(group) {
  ids<-V2_GROUPS[[group]];if(is.null(ids))stop('Unknown group')
  for(id in ids)render_discovery_v2(id)
  validate_fig1_refinement(ids)
  invisible(ids)
}
