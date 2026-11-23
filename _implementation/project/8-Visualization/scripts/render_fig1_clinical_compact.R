# Optional compact F layout. Source the Figure1 v2 setup before calling.
render_fig1_clinical_compact <- function() {
  check_owner_gate('Step6')
  id<-'Figure1F_compact'
  csvs<-paste0('8-Visualization/refinements/fig1_sfig1_20261026_v2/source_data/Figure1F',1:4,'.csv')
  d<-rbindlist(lapply(seq_along(csvs),function(i){
    x<-fread(file.path(STEP8_ROOT,csvs[i]));x[,TraitOrder:=i];x
  }),fill=TRUE)
  baseline_path<-paste0(T6,'step6_clinical_baseline_by_main_state.csv')
  baseline<-read_source(baseline_path)
  stopifnot(nrow(d)==12L,all(d[, .N,by=TraitOrder]$N==3L))
  for(i in seq_len(nrow(d))) {
    r<-d[i];owner<-baseline[Variable==r$Variable & Level==r$Level & FinalSubtype==r$FinalSubtype]
    stopifnot(nrow(owner)==1L,owner$AvailableN==r$AvailableN,owner$Summary==r$Summary)
    if(r$TraitOrder>1L)stopifnot(owner$N==r$PositiveN,abs(r$Value-100*owner$N/owner$AvailableN)<1e-9)
  }
  d[,FinalSubtype:=factor(FinalSubtype,levels=c('C1','C2','C3'))];setorder(d,TraitOrder,FinalSubtype)
  d[,Primary:=ifelse(TraitOrder==1L,sprintf('%.1f',Median),sprintf('%.1f%%',Value))]
  d[,Secondary:=ifelse(TraitOrder==1L,sprintf('[%.1f–%.1f]',Q1,Q3),paste0(PositiveN,'/',AvailableN))]
  d[,Tertiary:=ifelse(TraitOrder==1L,paste0('n = ',AvailableN),'')]
  data_path<-file.path(REFINE_DIR,'source_data',paste0(id,'.csv'));fwrite(d,data_path,na='NA')
  draw<-function() {
    grid::grid.newpage()
    txt<-function(label,x,y,size=7,bold=FALSE,col='grey20',just='centre')
      grid::grid.text(label,x,y,just=just,gp=grid::gpar(fontfamily='Arial',fontsize=size,
        fontface=if(bold)'bold' else 'plain',col=col))
    pal<-load_subtype_palette()$colors
    percent_col<-circlize::colorRamp2(c(0,50,100),c('#F8FAFC','#D5E7F0','#246085'))
    x<-c(.405,.635,.865);ys<-c(.735,.555,.375,.195);w<-.22;h<-.16
    txt('Clinical characteristics',.025,.965,9,TRUE,just='left')
    for(j in 1:3) {
      grid::grid.rect(x[j],.87,w,.075,gp=grid::gpar(fill=pal[paste0('C',j)],col=NA))
      txt(paste0('C',j),x[j],.87,8.5,TRUE,'grey15')
    }
    row_titles<-c('Tumor size','Stage III/IV','Vascular invasion','Large-duct type')
    row_notes<-c('Median [IQR], cm','Observed proportion','Observed proportion','Observed pathology')
    for(i in 1:4) {
      ss<-d[TraitOrder==i];yy<-ys[i]
      txt(row_titles[i],.025,yy+.035,7.5,TRUE,just='left')
      txt(row_notes[i],.025,yy-.005,6.2,just='left')
      txt(paste0('Available n = ',sum(ss$AvailableN)),.025,yy-.045,6.2,col='grey40',just='left')
      for(j in 1:3) {
        r<-ss[FinalSubtype==paste0('C',j)]
        fill<-if(i==1)'grey96' else percent_col(r$Value)
        grid::grid.rect(x[j],yy,w,h,gp=grid::gpar(fill=fill,col=NA))
        txt(r$Primary,x[j],yy+.027,8.5,TRUE)
        txt(r$Secondary,x[j],yy-.015,7)
        if(i==1)txt(r$Tertiary,x[j],yy-.054,6.3,col='grey40')
      }
    }
    txt('Proportion (%)',.025,.05,6.5,just='left')
    for(k in 1:101)grid::grid.rect(.305+(k-.5)*.35/101,.055,.35/101,.017,
      gp=grid::gpar(fill=percent_col(k-1),col=NA))
    txt('0',.305,.021,6);txt('50',.48,.021,6);txt('100',.655,.021,6)
    txt('n/N: positive / available',.98,.05,6.5,just='right')
  }
  pdf<-file.path(REFINE_DIR,'pdf',paste0(id,'.pdf'));save_panel(draw,pdf,114,66)
  grDevices::png(file.path(REFINE_DIR,'review',paste0(id,'.png')),width=114,height=66,units='mm',res=300,type='cairo',family='Arial')
  tryCatch(draw(),finally=grDevices::dev.off())
  caption<-paste('Compact descriptive clinical overview of principal states. The first row reports the locked median [IQR] in cm and available n; its neutral background has no quantitative color encoding.',
    'The other three rows report positive percentages and positive n / available N. One fixed 0–100% color scale is shared across all three binary traits, without row normalization.',
    'Endpoint-specific available totals are 312, 341, 248 and 93. Missing observations are not imputed or treated as negative. Only observed duct labels are used.',
    'The values are identical to the four standalone F1–F4 panels; no new statistics, confidence intervals or significance tests are introduced. This is an alternative compact layout, not an additional independent result.')
  writeLines(caption,file.path(REFINE_DIR,'source_data',paste0(id,'_caption.txt')))
  inputs<-sort(c(csvs,baseline_path,paste0(T6,'step6_clinical_survival_release_gate_status.csv'),
    '8-Visualization/scripts/render_fig1_clinical_compact.R','8-Visualization/scripts/render_fig1_discovery_v2.R',
    '8-Visualization/scripts/render_fig1_refinement.R','8-Visualization/scripts/render_bulk_figures.R',
    '8-Visualization/scripts/fig_style.R','8-Visualization/scripts/stat_annotate.R','8-Visualization/scripts/figure_config.json',SUBTYPE_PALETTE_CSV))
  hashes<-unname(tools::md5sum(file.path(STEP8_ROOT,inputs)));stopifnot(!anyNA(hashes))
  row<-data.table(Panel=id,Status='Rendered_review_draft',WidthMM=114,HeightMM=66,
    PDF=sub(paste0(STEP8_ROOT,'/'),'',pdf,fixed=TRUE),PDFMD5=unname(tools::md5sum(pdf)),
    SourcePaths=paste(inputs,collapse=';'),SourceMD5=paste(hashes,collapse=';'),
    SourceData=basename(data_path),SourceDataMD5=unname(tools::md5sum(data_path)))
  path<-file.path(REFINE_DIR,'figure_manifest.csv');m<-rbindlist(list(fread(path)[Panel!=id],row));setorder(m,Panel);fwrite(m,path)
  validate_fig1_refinement(id)
  message('Rendered compact alternative: ',id)
  invisible(pdf)
}
