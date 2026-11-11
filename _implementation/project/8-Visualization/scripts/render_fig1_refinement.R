# Figure 1 B–F and SFig1: rendering from locked sources, no upstream refitting.
# PDF-only manuscript exports; R-rendered PNGs are review derivatives.
if (!exists('STEP8_ROOT')) stop('Source fig_style.R, stat_annotate.R and render_bulk_figures.R first')
REFINE_DIR <- file.path(STEP8_ROOT, '8-Visualization/refinements/fig1_sfig1_20261026')
for (d in c('pdf', 'source_data', 'review')) dir.create(file.path(REFINE_DIR,d),recursive=TRUE,showWarnings=FALSE)
refine_sources <- character()
refine_read <- function(path) {
  refine_sources <<- unique(c(refine_sources,path))
  read_source(path)
}
refine_data <- function(id, d, suffix='') {
  data.table::fwrite(d,file.path(REFINE_DIR,'source_data',paste0(id,suffix,'.csv')),na='NA')
  invisible(d)
}
refine_master <- function(main=TRUE) {
  d <- refine_read(SUBTYPE_TRUTH_CSV)
  d[,Sample:=as.character(Sample)]
  stopifnot(nrow(d)==586L,!anyDuplicated(d$Sample))
  if(main)d<-d[SubtypeClass=='Main']
  if(main)stopifnot(nrow(d)==497L,identical(as.integer(table(factor(d$FinalSubtype,levels=c('C1','C2','C3')))),c(112L,223L,162L)))
  d
}
refine_theme <- function() {
  theme_step8()+theme(plot.title=element_text(size=10,face='bold'),
    strip.background=element_rect(fill='grey96',color=NA),strip.text=element_text(size=8,face='bold'),
    legend.position='bottom',legend.title=element_text(size=7),legend.text=element_text(size=7),
    panel.spacing=grid::unit(4,'mm'),plot.margin=margin(4,4,4,4,'mm'))
}
refine_caption <- function(id) {
  switch(id,
    Figure1B='Retained-method rank diagnostics in 586 discovery samples. NMF cophenetic and PAC and k-means mean silhouette are taken unchanged from the locked Step1 table. Each facet has its own vertical scale; k=3 is the retained operational choice. Points are metric estimates, not biological replicates; no confidence intervals or new tests are added. Consensus and hierarchical degeneration are exclusion evidence in S2. Modest silhouette values are retained without implying sharply separated biological classes.',
    Figure1C='Frozen corrected discovery markers across all 497 principal-state samples (C1 112, C2 223, C3 162). All unique genes from the 120-record Step1 marker set are shown. Rows are z-scored across these samples; the display color saturates at ±3 while source values are uncapped. Marker rows follow their frozen source order and are grouped by originating subtype; sample columns are ordered by subtype, Study and Sample, with no new clustering. Eight evenly spaced row labels per subtype are shown for readability; all genes and sample IDs are available in source data. This is a descriptive display of discovery-selected markers, not independent subtype validation.',
    Figure1D='Locked top-2,000-gene PCA coordinates for 497 principal-state samples; points show samples, colors denote C1/C2/C3. All principal samples are retained. No PCA is refitted, no ellipse is interpreted as uncertainty, and axis variance percentages are not invented. These coordinates are distinct from the all-gene preprocessing PCA in SFig1E. Study diagnostics are in SFig1E–F.',
    Figure1E='Within-study proportions of the 497 principal-state samples. Each bar sums to 100%; colored segments show C1/C2/C3, segment labels are sample counts and axis labels give the total main-state denominator. Explicit zero count rows are retained in source data. The 89 intermediate samples are excluded from this main-state denominator and retained for S3. No between-study significance test is added.',
    Figure1F='Available-case mapping between principal states and original classification fields, separated by Study and source field. Rows sum to 100%; cells show count and percentage, including explicit zero cells. Panel titles report available-label N / all principal-state N for that study; missing labels are excluded from row denominators and recorded in the coverage source table. Original systems are not pooled or renamed as Sia; no classification-accuracy interpretation or association test is added.',
    SFig1A='Curated expression samples by study, separated into RNA-seq discovery (586) and microarray projection (506). Bar endpoints show retained sample counts. Source data preserve input-matrix and excluded-column counts; excluded columns are not automatically interpreted as failed patient samples. No new filtering is performed.',
    SFig1B='RNA-matched samples among the recorded assay samples. Cell text is matched-to-curated-expression samples / recorded assay samples, and fill is this percentage. Expression rows use their own curated expression cohort as denominator. DNA denominators are samples represented in harmonized mutation records, not complete assayed-patient rosters or gene-level callability. Protein denominators are matrix columns. An × means absent from this source inventory, and ? means unknown coverage for an available entry; these are distinct from measured 0%. This is not the fraction of the expression cohort with an assay. Denominator definitions are retained in source data.',
    SFig1C='Processed feature retention by source dataset, separated by discovery and projection roles. The horizontal position is retained unique-source features (%); right-hand labels give final processed features / unique source feature IDs. Assay/platform-specific feature definitions differ, so percentages are preprocessing descriptions, not a cross-platform quality ranking.',
    SFig1D='Locked protein coordinates before and after ComBat for all 298 assay samples, including samples outside the matched expression subset. Study is encoded consistently by color and shape. Both stages use shared display limits. The locked coordinates are reused without fitting PCA or ComBat; visual overlap does not establish complete batch removal or independent subtype validation.',
    SFig1E='Locked RNA preprocessing PCA before and after Study correction, using all 10,531 genes and 586 discovery samples. Study colors and shapes are consistent between facets, with shared display limits. This all-gene PCA is distinct from the top-2,000-gene PCA used in Figure1D. Correction used limma::removeBatchEffect without a biological design term in the historical fit. Mixing alone does not establish absence of confounding.',
    SFig1F='Locked Study eta-squared for PCs 1–10 before and after RNA correction. Values and inferential metadata are copied from the owner table; no ANOVA or p-value is recalculated. These describe association with Study in the preprocessing representation and do not prove preservation of all biological effects. Both null and nonzero residual associations are retained.',
    stop('Unknown refinement panel'))
}
refine_plot <- function(id) {
  pal<-load_subtype_palette()$colors
  if(id=='Figure1B') {
    d<-refine_read(paste0(T1,'step1_method_selection_metrics_long.csv'))[Method %in% c('NMF','KMeans')]
    labs<-c('NMF | Cophenetic (higher is better)'='NMF · cophenetic\nHigher is better',
      'NMF | PAC (lower is better)'='NMF · PAC\nLower is better',
      'KMeans | Mean silhouette (higher is better)'='k-means · silhouette\nHigher is better')
    stopifnot(nrow(d)==15L,all(d$Facet %in% names(labs)))
    d[,DisplayFacet:=factor(labs[Facet],levels=unname(labs))];refine_data(id,d)
    return(ggplot(d,aes(K,Value))+geom_vline(xintercept=3,color='grey65',linetype=2,linewidth=.35)+
      geom_line(color='grey35',linewidth=.45)+geom_point(size=1.8,color='grey35')+
      geom_point(data=d[K==3],shape=21,size=3,stroke=.65,fill='white',color='black')+
      facet_wrap(~DisplayFacet,scales='free_y',nrow=1)+scale_x_continuous(breaks=2:6)+
      scale_y_continuous(expand=expansion(mult=c(.12,.16)))+
      labs(title='Rank selection for the retained methods',x='Number of clusters (k)',y=NULL)+refine_theme())
  }
  if(id=='Figure1C') {
    markers<-refine_read(paste0(T1,'step1_bulk_principal_state_markers_for_step8.csv'))
    markers<-markers[!duplicated(Gene)];m<-refine_master();m[,FinalSubtype:=factor(FinalSubtype,levels=c('C1','C2','C3'))];setorder(m,FinalSubtype,Study,Sample)
    rel<-'1-Clustering/data/combined_expr_after_correction.rds';refine_sources<<-unique(c(refine_sources,rel))
    env<-new.env();load(file.path(STEP8_ROOT,rel),envir=env);expr<-env$corrected_expr_limma
    stopifnot(all(markers$Gene %in% rownames(expr)),all(m$Sample %in% colnames(expr)))
    raw<-as.matrix(expr[markers$Gene,m$Sample,drop=FALSE]);stopifnot(all(is.finite(raw)),all(apply(raw,1,sd)>0))
    z<-t(scale(t(raw)));stopifnot(ncol(z)==497L)
    refine_data(id,data.table(Gene=rownames(z),z), '_row_zscore')
    refine_data(id,data.table(Gene=rownames(raw),raw),'_expression')
    refine_data(id,m[,.(Sample,Study,FinalSubtype)],'_columns')
    markers[,RowIndex:=seq_len(.N)];markers[,LabelShown:=seq_len(.N) %in% unique(round(seq(1,.N,length.out=8))),by=Subtype]
    refine_data(id,markers,'_rows')
    return(function() {
      top<-ComplexHeatmap::HeatmapAnnotation(Subtype=as.character(m$FinalSubtype),Study=m$Study,
        col=list(Subtype=pal,Study=load_study_palette()$colors),simple_anno_size=grid::unit(2.4,'mm'),
        annotation_name_gp=grid::gpar(fontsize=7),show_legend=c(Subtype=FALSE,Study=TRUE),
        annotation_legend_param=list(Study=list(nrow=2,title_gp=grid::gpar(fontsize=7),labels_gp=grid::gpar(fontsize=7))))
      h<-ComplexHeatmap::Heatmap(z,name='Row z-score',col=circlize::colorRamp2(c(-3,0,3),c(SCALE_EFFECT_LOW,SCALE_EFFECT_MID,SCALE_EFFECT_HIGH)),
        cluster_rows=FALSE,cluster_columns=FALSE,cluster_row_slices=FALSE,cluster_column_slices=FALSE,
        row_split=factor(markers$Subtype,levels=c('C1','C2','C3')),column_split=m$FinalSubtype,
        row_title=c('C1 markers','C2 markers','C3 markers'),row_title_gp=grid::gpar(fontsize=8),
        column_title=c('C1 · n = 112','C2 · n = 223','C3 · n = 162'),column_title_gp=grid::gpar(fontsize=8,fontface='bold'),
        row_gap=grid::unit(1.5,'mm'),column_gap=grid::unit(1.5,'mm'),show_row_names=TRUE,show_column_names=FALSE,
        row_labels=ifelse(markers$LabelShown,markers$Gene,''),row_names_gp=grid::gpar(fontsize=7,fontfamily='Arial'),
        top_annotation=top,use_raster=FALSE,
        heatmap_legend_param=list(at=c(-3,0,3),direction='horizontal',title_position='topcenter',
          title_gp=grid::gpar(fontsize=7),labels_gp=grid::gpar(fontsize=7)))
      ComplexHeatmap::draw(h,heatmap_legend_side='bottom',annotation_legend_side='bottom',merge_legends=TRUE,
        padding=grid::unit(c(4,4,4,4),'mm'))
    })
  }
  if(id=='Figure1D') {
    coords<-refine_read('governance/data/pca_coordinates_top2000.csv');coords[,Sample:=as.character(Sample)]
    m<-refine_master();d<-merge(m[,.(Sample,Study,FinalSubtype)],coords[,.(Sample,PC1,PC2)],by='Sample',all.x=TRUE)
    stopifnot(nrow(d)==497L,all(is.finite(d$PC1)),all(is.finite(d$PC2)))
    d<-d[order(Sample)];refine_data(id,d)
    counts<-table(d$FinalSubtype);lab<-setNames(paste0(names(counts),'  (n = ',counts,')'),names(counts))
    return(ggplot(d,aes(PC1,PC2,color=FinalSubtype))+geom_hline(yintercept=0,color='grey92',linewidth=.3)+
      geom_vline(xintercept=0,color='grey92',linewidth=.3)+geom_point(size=1.5,alpha=.8)+
      scale_colour_manual(values=pal,breaks=names(lab),labels=lab)+coord_fixed()+
      labs(title='Principal-state expression structure',x='PC1',y='PC2',color=NULL)+refine_theme()+
      guides(colour=guide_legend(override.aes=list(alpha=1,size=2.5))))
  }
  if(id=='Figure1E') {
    m<-refine_master();d<-m[,.N,by=.(Study,FinalSubtype)]
    d<-merge(CJ(Study=sort(unique(m$Study)),FinalSubtype=c('C1','C2','C3')),d,all.x=TRUE);d[is.na(N),N:=0L]
    d[,Total:=sum(N),by=Study];d[,Fraction:=N/Total];setorder(d,Study,FinalSubtype)
    d[,Midpoint:=cumsum(Fraction)-Fraction/2,by=Study];stopifnot(sum(d$N)==497L)
    refine_data(id,d);study_levels<-sort(unique(d$Study));d[,Study:=factor(Study,levels=rev(study_levels))]
    study_totals<-d$Total[match(study_levels,as.character(d$Study))]
    lab<-setNames(paste0(study_levels,'  (n = ',study_totals,')'),study_levels)
    return(ggplot(d,aes(Study,Fraction,fill=FinalSubtype))+geom_col(width=.64,color='white',linewidth=.35,position=position_stack(reverse=TRUE))+
      geom_text(aes(y=Midpoint,label=ifelse(N>0,N,'')),colour='grey15',
        size=2.6,show.legend=FALSE)+
      scale_fill_manual(values=pal,breaks=c('C1','C2','C3'))+scale_x_discrete(labels=lab)+
      scale_y_continuous(labels=scales::label_percent(),breaks=c(0,.25,.5,.75,1),expand=expansion(mult=0))+coord_flip()+
      labs(title='Principal-state composition across studies',x=NULL,y='Within-study proportion',fill=NULL)+refine_theme())
  }
  if(id=='Figure1F') {
    m<-refine_master();d<-melt(m,id.vars=c('Sample','Study','FinalSubtype'),measure.vars=c('Subclass','Subclass-add'),variable.name='SourceField',value.name='OriginalLabel')
    d[,Available:=!is.na(OriginalLabel)&!OriginalLabel %in% c('NA','w/o','')]
    coverage<-d[,.(AvailableN=sum(Available),MainN=.N,MissingN=sum(!Available)),by=.(Study,SourceField)]
    refine_data(id,coverage,'_coverage');d<-d[Available==TRUE]
    d<-d[,.N,by=.(Study,SourceField,OriginalLabel,FinalSubtype)]
    keys<-unique(d[,.(Study,SourceField,OriginalLabel)])
    complete<-keys[,.(FinalSubtype=c('C1','C2','C3')),by=.(Study,SourceField,OriginalLabel)]
    d<-merge(complete,d,all.x=TRUE);d[is.na(N),N:=0L];d[,RowN:=sum(N),by=.(Study,SourceField,OriginalLabel)];d[,Fraction:=N/RowN]
    d<-merge(d,coverage,by=c('Study','SourceField'))
    taxonomy_labels<-c('2019Ahn/Subclass'='2019Ahn · A/B classes','2022Dong/Subclass'='2022Dong · S1–S4',
      '2022Dong/Subclass-add'='2022Dong · IG1–IG3','2023Deng/Subclass'='2023Deng · S-I–S-III','2023Jeon/Subclass'='2023Jeon · LC classes')
    d[,TaxonomyLabel:=taxonomy_labels[paste(Study,SourceField,sep='/')]];stopifnot(!anyNA(d$TaxonomyLabel))
    d[,Taxonomy:=paste0(TaxonomyLabel,'\nAvailable n = ',AvailableN,' / ',MainN)]
    d[,DisplayLabel:=paste0(OriginalLabel,'  (n = ',RowN,')')]
    label_order<-rev(unique(d$DisplayLabel[order(d$Study,d$SourceField,d$OriginalLabel)]))
    d[,DisplayLabel:=factor(DisplayLabel,levels=label_order)];refine_data(id,d)
    stopifnot(all(abs(d[,sum(Fraction),by=.(Study,SourceField,OriginalLabel)]$V1-1)<1e-10))
    return(ggplot(d,aes(FinalSubtype,DisplayLabel,fill=Fraction))+geom_tile(color='white',linewidth=.6)+
      geom_text(aes(label=sprintf('%d\n%.0f%%',N,100*Fraction),color=Fraction>.58),size=2.5,lineheight=.95)+
      scale_colour_manual(values=c('FALSE'='black','TRUE'='white'),guide='none')+
      scale_fill_gradient(low='white',high='grey20',limits=c(0,1),labels=scales::label_percent(),breaks=c(0,.5,1))+
      facet_wrap(~Taxonomy,ncol=3,scales='free_y')+scale_x_discrete(drop=FALSE,position='top')+
      labs(title='Correspondence with original classifications',x=NULL,y=NULL,fill='Row proportion')+refine_theme()+
      theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_text(size=7),strip.text=element_text(size=7)))
  }
  if(id=='SFig1A') {
    d<-refine_read('0-Data/qc/step0_cohort_inclusion_manifest.csv');refine_data(id,d)
    d[,Study:=factor(Study,levels=rev(sort(unique(Study))))]
    return(ggplot(d,aes(Study,FinalSamples,fill=AnalysisRole))+geom_col(width=.65)+
      geom_text(aes(label=FinalSamples),hjust=-.2,size=2.8)+coord_flip()+
      scale_fill_manual(values=load_assay_palette()$colors,guide='none')+
      facet_wrap(~AnalysisRole,ncol=2,scales='free')+scale_y_continuous(expand=expansion(mult=c(0,.18)))+
      labs(title='Curated expression cohorts',x=NULL,y='Included samples')+refine_theme())
  }
  if(id=='SFig1B') {
    d<-refine_read('0-Data/qc/step0_assay_availability_matrix.csv');refine_data(id,d)
    d[,Study:=factor(Study,levels=rev(sort(unique(Study))))]
    d[,Display:=ifelse(!Available,'×',ifelse(is.na(CoveragePercent),'?',DisplayLabel))]
    d[,DisplayCoverage:=ifelse(Available,CoveragePercent,NA_real_)]
    d[,Ink:=ifelse(is.na(DisplayCoverage),'black',ifelse(DisplayCoverage>60,'white','black'))]
    return(ggplot(d,aes(Assay,Study,fill=DisplayCoverage))+geom_tile(color='white',linewidth=.6)+
      geom_text(aes(label=Display,color=Ink),size=2.7)+scale_colour_identity()+
      scale_fill_gradient(low='white',high='grey20',limits=c(0,100),na.value='grey88',breaks=c(0,50,100))+
      scale_x_discrete(position='top',labels=function(x)gsub('Protein \\(','Protein\n(',x))+
      labs(title='Recorded assay samples matched to expression',x=NULL,y=NULL,fill='RNA-matched (%)')+refine_theme()+
      theme(axis.line=element_blank(),axis.ticks=element_blank()))
  }
  if(id=='SFig1C') {
    d<-refine_read('0-Data/qc/step0_expression_overview.csv');refine_data(id,d)
    d[,Dataset:=factor(Dataset,levels=rev(unique(Dataset)))]
    return(ggplot(d,aes(RetainedPercentOfUniqueSource,Dataset))+geom_segment(aes(x=0,xend=RetainedPercentOfUniqueSource,yend=Dataset),color='grey80',linewidth=.6)+
      geom_point(aes(fill=AnalysisRole),shape=21,size=2.6,color='white',stroke=.3)+
      geom_text(aes(x=110,label=paste0(format(FinalProcessedFeatures,big.mark=','),' / ',format(UniqueSourceFeatureIDs,big.mark=','))),hjust=0,size=2.5)+
      scale_fill_manual(values=load_assay_palette()$colors,guide='none')+facet_wrap(~AnalysisRole,ncol=2,scales='free_y')+
      scale_x_continuous(breaks=c(0,50,100),limits=c(0,175),expand=expansion(mult=0))+
      labs(title='Processed expression feature retention',x='Retained source features (%)',y=NULL)+refine_theme())
  }
  if(id %in% c('SFig1D','SFig1E')) {
    protein<-id=='SFig1D'
    d<-refine_read(if(protein)'0-Data/qc/step0_protein_pca_scores.csv' else paste0(T1,'step1_batch_pca_scores.csv'))
    stage<-if(protein)'Stage' else 'CorrectionState';stages<-if(protein)c('Pre-ComBat','Post-ComBat') else c('Before correction','Locked limma correction')
    stopifnot(setequal(unique(d[[stage]]),stages))
    d[,DisplayStage:=factor(get(stage),levels=stages)];stopifnot(!anyNA(d$DisplayStage));refine_data(id,d)
    n<-if(protein)298L else 586L;stopifnot(all(d[,.N,by=DisplayStage]$N==n))
    studies<-sort(unique(d$Study));shapes<-setNames(c(16,17,15,18,3,7,8)[seq_along(studies)],studies)
    return(ggplot(d,aes(PC1,PC2,color=Study,shape=Study))+geom_point(size=1.2,alpha=.75,stroke=.45)+
      facet_wrap(~DisplayStage,ncol=2)+scale_colour_manual(values=load_study_palette()$colors)+scale_shape_manual(values=shapes)+
      labs(title=if(protein)'Protein preprocessing: study structure' else 'RNA preprocessing: study structure',x='PC1',y='PC2')+
      refine_theme()+guides(color=guide_legend(nrow=2,override.aes=list(alpha=1,size=2)),shape=guide_legend(nrow=2)))
  }
  if(id=='SFig1F') {
    d<-refine_read(paste0(T1,'step1_batch_pc_study_association.csv'));refine_data(id,d)
    stopifnot(setequal(unique(d$CorrectionState),c('Before correction','Locked limma correction')))
    d[,CorrectionState:=factor(CorrectionState,levels=c('Before correction','Locked limma correction'),labels=c('Before correction','After correction'))]
    stopifnot(!anyNA(d$CorrectionState),nrow(d)==20L,all(is.finite(d$StudyEtaSquared)))
    return(ggplot(d,aes(PCIndex,StudyEtaSquared,linetype=CorrectionState,shape=CorrectionState))+geom_line(color='grey35',linewidth=.45)+
      geom_point(size=2)+scale_shape_manual(values=c(16,1))+scale_linetype_manual(values=c('solid','dashed'))+
      scale_x_continuous(breaks=1:10,labels=paste0('PC',1:10))+scale_y_continuous(limits=c(0,1),breaks=c(0,.25,.5,.75,1))+
      labs(title='RNA principal components associated with Study',x=NULL,y='Study eta-squared',linetype=NULL,shape=NULL)+refine_theme())
  }
  stop('Unknown panel: ',id)
}
render_fig1_refinement <- function(id) {
  check_owner_gate('Step0');check_owner_gate('Step1')
  refine_sources <<- c('8-Visualization/scripts/render_fig1_refinement.R','8-Visualization/scripts/fig_style.R',
    '8-Visualization/scripts/stat_annotate.R','8-Visualization/scripts/render_bulk_figures.R','8-Visualization/scripts/figure_config.json',
    SUBTYPE_PALETTE_CSV,STUDY_PALETTE_CSV,'governance/data/assay_type_color_dictionary.csv',
    '0-Data/qc/step0_release_gate_checks.csv',paste0(T1,'step1_clustering_release_gate_status.csv'))
  sizes<-list(Figure1B=c(174,68),Figure1C=c(174,135),Figure1D=c(114,95),Figure1E=c(114,85),Figure1F=c(174,108),
    SFig1A=c(174,90),SFig1B=c(174,115),SFig1C=c(174,90),SFig1D=c(174,95),SFig1E=c(174,95),SFig1F=c(114,85))
  size<-sizes[[id]];if(is.null(size))stop('No output size contract')
  plot<-refine_plot(id);pdf<-file.path(REFINE_DIR,'pdf',paste0(id,'.pdf'))
  save_panel(plot,pdf,size[1],size[2])
  # Same R draw path and dimensions on the QA device; no final PNG export contract.
  png<-file.path(REFINE_DIR,'review',paste0(id,'.png'))
  grDevices::png(png,width=size[1],height=size[2],units='mm',res=180,type='cairo',family='Arial')
  tryCatch({if(is.function(plot))plot() else print(plot)},finally=grDevices::dev.off())
  writeLines(refine_caption(id),file.path(REFINE_DIR,'source_data',paste0(id,'_caption.txt')))
  inputs<-sort(unique(refine_sources));hashes<-unname(tools::md5sum(file.path(STEP8_ROOT,inputs)))
  stopifnot(!anyNA(hashes),file.info(pdf)$size>0)
  row<-data.table(Panel=id,Status='Rendered_review_draft',WidthMM=size[1],HeightMM=size[2],
    PDF=sub(paste0(STEP8_ROOT,'/'),'',pdf,fixed=TRUE),PDFMD5=unname(tools::md5sum(pdf)),
    SourcePaths=paste(inputs,collapse=';'),SourceMD5=paste(hashes,collapse=';'),
    SourceData=paste(list.files(file.path(REFINE_DIR,'source_data'),pattern=paste0('^',id,'[_.]'),full.names=FALSE),collapse=';'))
  path<-file.path(REFINE_DIR,'figure_manifest.csv')
  old<-if(file.exists(path))fread(path) else data.table();if(nrow(old))old<-old[Panel!=id]
  manifest<-rbindlist(list(old,row),fill=TRUE);setorder(manifest,Panel);fwrite(manifest,path)
  message('Rendered review draft: ',id)
  invisible(pdf)
}

validate_fig1_refinement <- function(ids) {
  m<-fread(file.path(REFINE_DIR,'figure_manifest.csv'))
  stopifnot(all(ids %in% m$Panel),!anyDuplicated(m$Panel))
  for(id in ids) {
    r<-m[Panel==id];pdf<-file.path(STEP8_ROOT,r$PDF)
    stopifnot(identical(unname(tools::md5sum(pdf)),r$PDFMD5))
    src<-strsplit(r$SourcePaths,';',fixed=TRUE)[[1]]
    stopifnot(identical(paste(unname(tools::md5sum(file.path(STEP8_ROOT,src))),collapse=';'),r$SourceMD5))
    info<-system2('pdfinfo',shQuote(pdf),stdout=TRUE)
    stopifnot(grepl('^Pages:[[:space:]]+1$',grep('^Pages:',info,value=TRUE)))
    size_line<-grep('^Page size:',info,value=TRUE)
    actual<-as.numeric(strsplit(trimws(sub('^Page size:[[:space:]]+','',size_line)),'[[:space:]]+')[[1]][c(1,3)])*25.4/72
    # Cairo records integer PDF points: one-point rounding is at most 0.353 mm.
    stopifnot(all(abs(actual-c(r$WidthMM,r$HeightMM))<.36))
  }
  message('Source hashes, output hashes, page counts and dimensions verified: ',paste(ids,collapse=', '))
  invisible(TRUE)
}
