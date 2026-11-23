# R-only Figure4 extension renderer. New statistical analyses are owner tables.
# Single-page PDFs, Arial 7-9 pt, 174 mm width; no in-plot panel letters.
if (!exists('STEP8_ROOT') || !exists('read_source')) stop('Source shared Step8 helpers first')
FIG4_GENOMICS <- '4-Multiomic/refinements/fig4_genomics_20261026/'
FIG4_REVIEW <- '8-Visualization/refinements/fig4_genomics_20261026/'
fig4_read <- function(name) read_source(paste0(FIG4_GENOMICS,name))
fig4_variant_colors <- c('Missense'='#4B9B6E','Nonsense'='#D97933',
  'Frameshift deletion'='#467AA8','Frameshift insertion'='#8C6BB1',
  'In-frame deletion'='#D5B44B','In-frame insertion'='#319E9A',
  'Splice'='#B34B5E','Other non-silent'='#9B8A7A','Multi-hit'='#333333',
  'No recorded variant'='#E9E9E9')

register_fig4_genomics <- function() {
  gate <- fig4_read('technical_checks.csv'); stopifnot(all(gate$Status=='PASS'))
  # Verify analysis code/input/output hashes before using any cached inference.
  hashes <- fig4_read('input_hashes.csv')
  stopifnot(identical(hashes$MD5,unname(tools::md5sum(file.path(STEP8_ROOT,hashes$Path)))))
  hashes <- fig4_read('output_hashes.csv')
  stopifnot(identical(hashes$MD5,unname(tools::md5sum(file.path(STEP8_ROOT,FIG4_GENOMICS,hashes$Path)))))
  m <- read_source('8-Visualization/tables/figure_manifest.csv')
  ids <- c('A','G','H','I')
  paths <- c(A='Figure4A_driver_oncoprint.pdf',G='Figure4G_mutation_frequency.pdf',
    H='Figure4H_cnv_frequency.pdf',I='Figure4I_mutation_cooccurrence.pdf')
  roles <- c(A='Driver mutation waterfall with assay-roster coverage',G='Mutation frequencies across principal subtypes',
    H='CNV gain/loss frequencies in 2022Dong',I='Exploratory within-subtype mutation co-occurrence / exclusivity')
  src <- list(A=c('mutation_display_long.csv','sample_order.csv','gene_order.csv'),
    G=c('mutation_frequency.csv','gene_order.csv'),H=c('cnv_frequency.csv','gene_order.csv'),
    I=c('pairwise_associations.csv','pairwise_burden_sensitivity.csv','gene_order.csv'))
  notes <- c(A='All 308 assayed principal samples retained; no recorded variant is not proof of callable wild type; top counts are not TMB.',
    G='Locked Study-stratified CMH subtype-vs-others q values; BH within each subtype in the original owner analysis.',
    H='2022Dong only, n=219. Gain >=0.3 and loss <=-0.3 continuous CNV; exploratory threshold calls, not GISTIC significance.',
    I='Exploratory exact conditional Study-stratified association; BH over 315 pair/subtype slots. Minimum 5 positive and 5 negative per gene. Gray cells untested. No test of between-subtype interaction; source callability unavailable.')
  scopes <- c(A='DNA-assayed main n=308 (63/153/92)',G='DNA-assayed main n=308 (63/153/92)',
    H='2022Dong CNV main n=219 (47/105/67)',I='DNA-assayed main n=308; broad-record burden sensitivity n=248')
  common <- c('4-Multiomic/scripts/prepare_fig4_genomics.R','8-Visualization/scripts/render_fig4_genomics.R',
    '8-Visualization/scripts/fig_style.R','8-Visualization/scripts/stat_annotate.R',
    '8-Visualization/scripts/figure_config.json','governance/data/subtype_label_color_dictionary.csv',
    '0-Data/qc/study_color_dictionary.csv',paste0(FIG4_GENOMICS,c('input_hashes.csv','output_hashes.csv','technical_checks.csv')))
  for (id in ids) {
    row <- copy(m[Figure=='Figure 4' & Panel=='A'][1])
    row[,`:=`(Panel=id,Role=roles[[id]],AnalysisOwner='Step4',SourceFigure='',
      SourceTable=paste(c(paste0(FIG4_GENOMICS,src[[id]]),common),collapse=';'),
      TargetOutput=paste0('8-Visualization/figures/main/',paths[[id]]),
      Status='AnalysisReady',SubtypeScope=scopes[[id]],Notes=paste(notes[[id]],'Standalone PDF; review draft.'),
      WidthMM='174',HeightMM=as.character(c(A=132,G=96,H=128,I=91)[[id]]),SourceMD5='',QAFile='',RenderedAt='')]
    m <- rbindlist(list(m[!(Figure=='Figure 4' & Panel==id)],row),use.names=TRUE,fill=TRUE)
  }
  # Write only the targeted rows into a freshly read manifest; preserve all other figures.
  current <- read_source('8-Visualization/tables/figure_manifest.csv')
  updated <- rbindlist(list(current[!(Figure=='Figure 4' & Panel %in% ids)],m[Figure=='Figure 4' & Panel %in% ids]),use.names=TRUE)
  setorder(updated,Figure,Panel)
  temp <- tempfile(tmpdir=file.path(STEP8_ROOT,'8-Visualization/tables'),fileext='.csv')
  fwrite(updated,temp,na='');stopifnot(file.rename(temp,file.path(STEP8_ROOT,'8-Visualization/tables/figure_manifest.csv')))
  invisible(ids)
}

build_fig4_genomics <- function(id) {
  genes <- fig4_read('gene_order.csv')[order(DisplayOrder),Gene]
  pal <- load_subtype_palette()$colors[c('C1','C2','C3')]
  if (id=='A') {
    d <- fig4_read('mutation_display_long.csv'); samples <- fig4_read('sample_order.csv')[order(DisplayOrder)]
    samples[,Sample:=as.character(Sample)]
    mat <- matrix('No recorded variant',length(genes),nrow(samples),dimnames=list(genes,samples$Sample))
    mat[cbind(match(d$Gene,genes),match(as.character(d$Sample),samples$Sample))] <- d$MutationClass
    types <- setdiff(names(fig4_variant_colors),'No recorded variant')
    bars <- sapply(types,function(ty)rowSums(mat==ty));rownames(bars)<-genes
    pct <- rowSums(mat!='No recorded variant')/ncol(mat)*100
    ns <- table(factor(samples$FinalSubtype,levels=names(pal)))
    return(function() {
      gp <- grid::gpar(fontfamily=STEP8_FONT$family,fontsize=7)
      top <- ComplexHeatmap::HeatmapAnnotation(
        `Recorded variants`=ComplexHeatmap::anno_barplot(samples$NonsynonymousVariantCount,
          border=FALSE,gp=grid::gpar(fill='grey40',col=NA),height=grid::unit(17,'mm'),
          axis_param=list(gp=gp)),annotation_name_gp=gp,annotation_name_side='left')
      bottom <- ComplexHeatmap::HeatmapAnnotation(Subtype=samples$FinalSubtype,Study=samples$Study,
        col=list(Subtype=pal,Study=load_study_palette()$colors),simple_anno_size=grid::unit(3,'mm'),
        annotation_name_gp=gp,annotation_name_side='left',annotation_legend_param=list(
          Subtype=list(title_gp=gp,labels_gp=gp,nrow=1),Study=list(title_gp=gp,labels_gp=gp,nrow=1)))
      right <- ComplexHeatmap::rowAnnotation(
        `Altered (%)`=ComplexHeatmap::anno_text(sprintf('%.1f%%',pct),gp=gp,just='left',width=grid::unit(10,'mm')),
        Samples=ComplexHeatmap::anno_barplot(bars,gp=grid::gpar(fill=fig4_variant_colors[types],col=NA),
          border=FALSE,width=grid::unit(19,'mm'),axis_param=list(gp=gp)),
        annotation_name_gp=gp,annotation_name_rot=0)
      h <- ComplexHeatmap::Heatmap(mat,name='Recorded consequence',col=fig4_variant_colors,
        cluster_rows=FALSE,cluster_columns=FALSE,cluster_column_slices=FALSE,
        column_split=factor(samples$FinalSubtype,levels=names(pal)),column_gap=grid::unit(1.3,'mm'),
        column_title=paste0(names(ns),'  n = ',as.integer(ns)),column_title_gp=grid::gpar(fontsize=8,fontface='bold'),
        show_column_names=FALSE,row_names_side='left',row_names_gp=grid::gpar(fontsize=8,fontface='italic'),
        rect_gp=grid::gpar(col=NA),top_annotation=top,bottom_annotation=bottom,right_annotation=right,
        heatmap_legend_param=list(title_gp=gp,labels_gp=gp,ncol=3,grid_height=grid::unit(2.6,'mm')))
      ComplexHeatmap::draw(h,heatmap_legend_side='bottom',annotation_legend_side='bottom',merge_legends=TRUE,
        padding=grid::unit(c(3,3,8,3),'mm'))
      grid::grid.text('Driver mutation landscape',x=grid::unit(3,'mm'),y=grid::unit(1,'npc')-grid::unit(3,'mm'),
        just=c('left','top'),gp=grid::gpar(fontsize=9,fontface='bold',fontfamily=STEP8_FONT$family))
    })
  }
  if (id=='G') {
    d <- fig4_read('mutation_frequency.csv');d[,Gene:=factor(Gene,levels=genes)]
    d[,FinalSubtype:=factor(FinalSubtype,levels=names(pal))]
    d[,Stars:=ifelse(CMHFDR<.05,fdr_label(CMHFDR),'')]
    n <- unique(d[,.(FinalSubtype,N)]);labels <- setNames(paste0(n$FinalSubtype,' (n = ',n$N,')'),n$FinalSubtype)
    pos <- position_dodge(width=.82)
    return(ggplot(d,aes(Gene,Fraction*100,fill=FinalSubtype,group=FinalSubtype))+
      geom_col(position=pos,width=.74)+geom_text(aes(label=Stars),position=pos,vjust=-.25,size=2.6)+
      scale_fill_manual(values=pal,labels=labels,name=NULL)+
      scale_y_continuous(expand=expansion(mult=c(0,.16)))+
      labs(title='Driver mutation frequencies',x=NULL,y='Samples with recorded mutation (%)')+
      theme_step8()+theme(axis.text.x=element_text(angle=45,hjust=1,face='italic',size=8),
        legend.position='top',legend.justification='left'))
  }
  if (id=='H') {
    d <- fig4_read('cnv_frequency.csv')
    d <- melt(d,id.vars=c('Gene','FinalSubtype','N'),measure.vars=c('GainFraction','LossFraction'),
      variable.name='Direction',value.name='Fraction')
    d[,Direction:=factor(Direction,levels=c('GainFraction','LossFraction'),labels=c('Copy-number gain','Copy-number loss'))]
    d[,Gene:=factor(Gene,levels=genes)];d[,FinalSubtype:=factor(FinalSubtype,levels=names(pal))]
    n <- unique(d[,.(FinalSubtype,N)]);labels <- setNames(paste0(n$FinalSubtype,' (n = ',n$N,')'),n$FinalSubtype)
    return(ggplot(d,aes(Gene,100*Fraction,fill=FinalSubtype))+geom_col(position=position_dodge(.82),width=.74)+
      facet_wrap(~Direction,ncol=1)+scale_fill_manual(values=pal,labels=labels,name=NULL)+
      scale_y_continuous(expand=expansion(mult=c(0,.08)))+
      labs(title='Driver-gene copy-number frequencies',subtitle='2022Dong cohort',x=NULL,y='Samples with CNV call (%)')+
      theme_step8()+theme(axis.text.x=element_text(angle=45,hjust=1,face='italic',size=8),
        legend.position='top',legend.justification='left',panel.spacing=grid::unit(5,'mm')))
  }
  if (id=='I') {
    d <- fig4_read('pairwise_associations.csv')
    d[,`:=`(x=match(Gene2,genes),y=match(Gene1,genes))]
    # Exact conditional OR boundary estimates are displayed at the color limit.
    d[,Log2OR:=ifelse(Status=='Tested',log2(OddsRatio),NA_real_)]
    d[,DisplayEffect:=pmax(-3,pmin(3,Log2OR))]
    d[,Stars:=ifelse(!is.na(FDR)&FDR<.05,ifelse(FDR<.01,'**','*'),'')]
    return(function() {
      grid::grid.newpage()
      grid::grid.text('Pairwise mutation associations within molecular subtypes',
        x=grid::unit(3,'mm'),y=grid::unit(87,'mm'),just='left',
        gp=grid::gpar(fontfamily=STEP8_FONT$family,fontsize=9,fontface='bold'))
      # Identical gene order and scale in all three triangular matrices.
      cell <- 2.50; side <- length(genes)*cell; left <- c(14,67,120); bottom <- 20
      colfun <- circlize::colorRamp2(c(-3,0,3),c('#2166AC','#FFFFFF','#B2182B'))
      for (k in seq_along(pal)) {
        st <- names(pal)[k];z <- d[Subtype==st];xx <- left[k]
        grid::grid.rect(x=grid::unit(xx+side/2,'mm'),y=grid::unit(bottom+side+7,'mm'),
          width=grid::unit(side,'mm'),height=grid::unit(5,'mm'),gp=grid::gpar(fill=pal[[st]],col=NA))
        grid::grid.text(paste0(st,'  (n = ',unique(z$N),')'),x=grid::unit(xx+side/2,'mm'),
          y=grid::unit(bottom+side+7,'mm'),gp=grid::gpar(fontsize=8,fontface='bold',col=if(st=='C1')'white' else 'black'))
        for (i in seq_len(nrow(z))) {
          fill <- if(is.na(z$DisplayEffect[i])) '#D6D6D6' else colfun(z$DisplayEffect[i])
          # Lower triangle: x = first gene, y = second gene (top-to-bottom order).
          cx <- xx+(z$y[i]-.5)*cell;cy<-bottom+side-(z$x[i]-.5)*cell
          grid::grid.rect(grid::unit(cx,'mm'),grid::unit(cy,'mm'),grid::unit(cell,'mm'),grid::unit(cell,'mm'),
            gp=grid::gpar(fill=fill,col='#C4C4C4',lwd=.25))
          if(nzchar(z$Stars[i])) grid::grid.text(z$Stars[i],grid::unit(cx,'mm'),grid::unit(cy,'mm'),
            gp=grid::gpar(fontsize=7,col=if(abs(z$DisplayEffect[i])>2)'white' else 'black'))
        }
        for (j in seq_along(genes)) {
          grid::grid.rect(grid::unit(xx+(j-.5)*cell,'mm'),grid::unit(bottom+side-(j-.5)*cell,'mm'),
            grid::unit(cell,'mm'),grid::unit(cell,'mm'),gp=grid::gpar(fill='#EFEFEF',col='#C4C4C4',lwd=.25))
          grid::grid.text(genes[j],grid::unit(xx-1,'mm'),grid::unit(bottom+side-(j-.5)*cell,'mm'),
            just='right',gp=grid::gpar(fontsize=7,fontface='italic',fontfamily=STEP8_FONT$family))
          grid::grid.text(genes[j],grid::unit(xx+(j-.5)*cell,'mm'),grid::unit(bottom-1,'mm'),
            just='right',rot=90,gp=grid::gpar(fontsize=7,fontface='italic',fontfamily=STEP8_FONT$family))
        }
      }
      # Compact shared color bar above matrices; no clipped side legends.
      vals <- seq(-3,3,length.out=60);cbx<-63;cby<-77;cbw<-43
      for(i in seq_along(vals))grid::grid.rect(grid::unit(cbx+(i-.5)*cbw/60,'mm'),grid::unit(cby,'mm'),
        grid::unit(cbw/60+.02,'mm'),grid::unit(2.5,'mm'),gp=grid::gpar(fill=colfun(vals[i]),col=NA))
      grid::grid.text('log2 odds ratio',grid::unit(59,'mm'),grid::unit(cby,'mm'),just='right',gp=grid::gpar(fontsize=7))
      for(i in 0:2)grid::grid.text(c('\u2264 -3','0','\u2265 3')[i+1],grid::unit(cbx+i*cbw/2,'mm'),grid::unit(73,'mm'),gp=grid::gpar(fontsize=7))
      grid::grid.text('Exclusivity',grid::unit(cbx,'mm'),grid::unit(69.5,'mm'),just='left',gp=grid::gpar(fontsize=7,col='#2166AC'))
      grid::grid.text('Co-occurrence',grid::unit(cbx+cbw,'mm'),grid::unit(69.5,'mm'),just='right',gp=grid::gpar(fontsize=7,col='#B2182B'))
      grid::grid.text('Gray: insufficient support   |   * q < 0.05; ** q < 0.01',
        grid::unit(87,'mm'),grid::unit(3.5,'mm'),gp=grid::gpar(fontsize=7,fontfamily=STEP8_FONT$family))
    })
  }
  stop('Unknown Figure4 genomics panel: ',id)
}

render_fig4_genomics <- function(id) {
  if (id=='A' && 'IRkernel' %in% loadedNamespaces()) {
    # The ComplexHeatmap 2.26.1 Jupyter grid.grabExpr replay loses nested annotation
    # viewports. Render on the SAME pinned R/Cairo runtime in an isolated process.
    # No private package state or global notebook graphics options are modified.
    script <- tempfile(fileext='.R');on.exit(unlink(script),add=TRUE)
    lines <- c(paste0('root <- ',paste(capture.output(dput(STEP8_ROOT)),collapse='')),
      sprintf('source(file.path(root, "%s"))',c('8-Visualization/scripts/fig_style.R',
        '8-Visualization/scripts/stat_annotate.R','8-Visualization/scripts/render_bulk_figures.R',
        '8-Visualization/scripts/render_fig4_genomics.R')),'render_fig4_genomics("A")')
    writeLines(lines,script)
    status <- system2(file.path(R.home('bin'),'exec/R'),c('--vanilla','--slave','-f',shQuote(script)))
    if(status!=0L)stop('Isolated R waterfall export failed')
    return(invisible(file.path(STEP8_ROOT,'8-Visualization/figures/main/Figure4A_driver_oncoprint.pdf')))
  }
  row <- preflight_subfigure('Figure 4',id); plot <- build_fig4_genomics(id)
  save_panel(plot,row$TargetOutput,as.numeric(row$WidthMM),as.numeric(row$HeightMM))
  sources <- strsplit(row$SourceTable,';',fixed=TRUE)[[1]]
  hashes <- paste(paste0(sources,'=',unname(tools::md5sum(file.path(STEP8_ROOT,sources)))),collapse=';')
  m <- read_source('8-Visualization/tables/figure_manifest.csv')
  m[Figure=='Figure 4' & Panel==id,`:=`(Status='Rendered',SourceMD5=hashes,QAFile='',RenderedAt=format(Sys.time(),tz='UTC',usetz=TRUE))]
  temp <- tempfile(tmpdir=file.path(STEP8_ROOT,'8-Visualization/tables'),fileext='.csv')
  fwrite(m,temp,na='');stopifnot(file.rename(temp,file.path(STEP8_ROOT,'8-Visualization/tables/figure_manifest.csv')))
  message('Rendered Figure4 genomics PDF: ',row$TargetOutput)
  invisible(file.path(STEP8_ROOT,row$TargetOutput))
}

preview_fig4_genomics <- function(ids=c('A','G','H','I'),dpi=300) {
  # R-native QA derivatives, never manuscript PNG deliverables.
  folder <- file.path(STEP8_ROOT,FIG4_REVIEW,'review');dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  for(id in ids) {
    row <- preflight_subfigure('Figure 4',id)
    grDevices::png(file.path(folder,paste0('Figure4',id,'.png')),width=as.numeric(row$WidthMM),
      height=as.numeric(row$HeightMM),units='mm',res=dpi,type='cairo',family=STEP8_FONT$family)
    tryCatch({p<-build_fig4_genomics(id);if(is.function(p))p() else print(p)},finally=grDevices::dev.off())
  }
}
