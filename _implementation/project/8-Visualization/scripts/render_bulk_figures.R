# Step8 rendering only: no hypothesis tests, model fits, bootstrap CIs or new cutoffs.
suppressPackageStartupMessages(library(data.table))
if (!exists('STEP8_ROOT')) stop('Source fig_style.R and stat_annotate.R first')
# ComplexHeatmap uses a null PostScript device for text measurement. Arial is
# metrically compatible with Helvetica there; actual Cairo exports embed Arial.
if (!'Arial' %in% names(grDevices::pdfFonts())) grDevices::pdfFonts(Arial=grDevices::pdfFonts('Helvetica')[[1]])
T1 <- '1-Clustering/publication/tables/'
T2 <- '2-TranscriptomicCharacterization/publication/tables/'
T3 <- '3-DuctType/publication/tables/'
T4 <- '4-Multiomic/publication/tables/'
T5 <- '5-ClassifierValidation/publication/tables/'
T6 <- '6-ClinicalSurvival/publication/tables/'
read_source <- function(path, columns=character()) {
  p <- file.path(STEP8_ROOT,path)
  if (!file.exists(p)) stop('Missing source: ',path)
  d <- if (grepl('figure_manifest[.]csv$',path)) fread(p,na.strings=c('NA',''),check.names=FALSE,colClasses='character') else fread(p,na.strings=c('NA',''),check.names=FALSE)
  if (!nrow(d) || !all(columns %in% names(d))) stop('Empty source or missing columns: ',path,' [',paste(setdiff(columns,names(d)),collapse=', '),']')
  d
}
master_table <- function(main=TRUE) {
  d <- read_source(SUBTYPE_TRUTH_CSV,c('Sample','Study','FinalSubtype','SubtypeClass'))
  if (anyDuplicated(d$Sample) || nrow(d)!=586) stop('Canonical sample interface changed')
  if (main) d <- d[SubtypeClass=='Main']
  d
}
subtype_scale <- function(aesthetic='colour',main=TRUE) {
  p <- load_subtype_palette()
  df <- if(main)p$df[p$df$SubtypeClass=='Main',] else p$df
  if (aesthetic=='fill') scale_fill_manual(values=setNames(df$ColorHex,df$FinalSubtype),drop=FALSE)
  else scale_colour_manual(values=setNames(df$ColorHex,df$FinalSubtype),drop=FALSE)
}
effect_scale <- function(values,name='Effect') {
  lim <- max(abs(values[is.finite(values)]))
  if (!length(lim) || !is.finite(lim) || lim==0) lim <- 1
  scale_fill_gradient2(low=SCALE_EFFECT_LOW,mid=SCALE_EFFECT_MID,high=SCALE_EFFECT_HIGH,midpoint=0,limits=c(-lim,lim),name=name)
}
# Display aggregation (counts/proportions/mean scores) is explicitly descriptive.
# Inferential effects, adjusted p/q values and CIs always come from owner tables.
forest_plot <- function(d, label, estimate, lower=NULL, upper=NULL, facet=NULL, log_axis=FALSE,ref=0,title='') {
  d <- copy(d); d[, .label := as.character(get(label))]; d[, .estimate := get(estimate)]
  d[, .valid := is.finite(.estimate) & (!log_axis | .estimate > 0)]
  # Preserve unplottable rows visibly rather than silently removing them.
  valid <- d[.valid==TRUE]; missing <- d[.valid==FALSE]
  p <- ggplot(d,aes(x=.estimate,y=.label))
  if(is.finite(ref))p<-p+geom_vline(xintercept=ref,color='grey65',linewidth=.3)
  if (!is.null(lower)) {
    d[, .lo := get(lower)];d[, .hi := get(upper)]
    intervals <- d[.valid & is.finite(.lo) & is.finite(.hi) & (!log_axis | (.lo > 0 & .hi > 0))]
    p <- p+geom_segment(data=intervals,aes(x=.lo,xend=.hi,yend=.label),linewidth=.4)
  }
  p <- p+geom_point(data=valid,size=1.7)
  if(nrow(missing))p<-p+geom_text(data=missing,aes(x=ref,label='NE/boundary'),size=2.1)
  if (!is.null(facet)) p <- p+facet_wrap(as.formula(paste('~',facet)),scales='free')
  if (log_axis) p <- p+scale_x_log10()
  p+labs(title=title,x=estimate,y=NULL)+theme_step8()
}

box_plot <- function(d,value,facet,title,main=TRUE) {
  ggplot(d,aes(x=FinalSubtype,y=.data[[value]],fill=FinalSubtype))+geom_boxplot(width=.6,outlier.shape=NA,linewidth=.3)+
    geom_point(position=position_jitter(width=.15,seed=123),size=.35,alpha=.35)+subtype_scale('fill',main)+
    facet_wrap(as.formula(paste('~',facet)),scales='free_y')+labs(title=title,x=NULL,y=value)+theme_step8()+theme(legend.position='none')
}
plot_tile <- function(d,x,y,z,title,diverging=TRUE) {
  p <- ggplot(d,aes(x=.data[[x]],y=.data[[y]],fill=.data[[z]]))+geom_tile()+labs(title=title,x=NULL,y=NULL)+theme_step8()
  if (diverging) p <- p+effect_scale(d[[z]],z)
  else p <- p+scale_fill_gradient(low='white',high='grey25',name=z)
  p
}
check_owner_gate <- function(owner) {
  gates <- list(Step0='0-Data/qc/step0_release_gate_checks.csv',Step1=paste0(T1,'step1_clustering_release_gate_status.csv'),
    Step2='2-TranscriptomicCharacterization/reports/transcriptomic_characterization_release_gate_status.csv',
    Step3=paste0(T3,'step3_release_gate_status.csv'),Step4=paste0(T4,'step4_multiomic_release_gate_status.csv'),
    Step5=paste0(T5,'step5_release_gate_status.csv'),Step6=paste0(T6,'step6_clinical_survival_release_gate_status.csv'))
  if (owner %in% names(gates)) {
    d <- read_source(gates[[owner]])
    status_col <- intersect(c('Status','Pass','Passed'),names(d))[1]
    if (is.na(status_col)) stop('Gate status missing: ',owner)
    blocking <- if ('Blocking' %in% names(d))tolower(as.character(d$Blocking)) %in% c('true','1') else rep(TRUE,nrow(d))
    if (any(blocking & !grepl('^(PASS|QUALIFIED|DISCLOSED|TRUE)',toupper(as.character(d[[status_col]]))))) stop('Owner blocking gate failed: ',owner)
  }
}
preflight_subfigure <- function(figure,subfigure) {
  m <- read_source('8-Visualization/tables/figure_manifest.csv')
  row <- m[Figure==figure & Panel==subfigure]
  if (nrow(row)!=1) stop('Missing/duplicate subfigure contract: ',figure,' ',subfigure)
  if (row$Status=='Blocked' || row$AnalysisOwner=='Step7') stop('BLOCKED: ',row$Notes)
  check_owner_gate(row$AnalysisOwner)
  specialized <- c(wgcna=paste0(T2,'step2_wgcna_release_gate_status.csv'),
                    tide=paste0(T2,'step2_tide_release_gate_status.csv'),
                    oncopredict=paste0(T5,'step5_oncopredict_release_gate_status.csv'),
                    c3_clinical_panel=paste0(T5,'step5_c3_clinical_panel_release_gate_status.csv'),
                    incremental_prognostic=paste0(T6,'step6_targeted_upgrade_release_gate_status.csv'))
  for (name in names(specialized)) if (grepl(name,row$SourceTable,fixed=TRUE)) {
    gate <- read_source(specialized[[name]])
    col <- intersect(c('Status','Pass','Passed'),names(gate))[1]
    block <- if('Blocking' %in% names(gate))tolower(as.character(gate$Blocking)) %in% c('true','1') else rep(TRUE,nrow(gate))
    if(any(block & !grepl('^(PASS|QUALIFIED|TRUE)',toupper(as.character(gate[[col]])))))stop('Specialized gate failed: ',name)
  }
  for (s in strsplit(row$SourceTable,';',fixed=TRUE)[[1]]) {
    p <- file.path(STEP8_ROOT,s)
    if (!file.exists(p) || file.info(p)$size==0) stop('Missing/empty source: ',s)
  }
  row
}
# Individual subfigures only. Facets compare like quantities within one graph;
# no patchwork/cowplot/composite main-figure assembly is performed.
build_subfigure <- function(figure,id) {
  key <- paste(figure,id)
  if (figure=='Figure 1') {
    if (id=='A') {
      counts<-read_source('governance/data/subtype_counts.csv');cohort<-read_source('0-Data/qc/step0_cohort_inclusion_manifest.csv')
      discovery<-cohort[AnalysisRole=='Discovery RNA-seq',sum(FinalSamples)];projection<-cohort[AnalysisRole=='Microarray projection',sum(FinalSamples)]
      states<-counts[SubtypeClass=='Main'];pal<-load_subtype_palette()$colors;apal<-load_assay_palette()$colors
      stopifnot(sum(counts$N)==discovery)
      return(function() {
        grid::grid.newpage()
        box<-function(x,y,w,h,text='',fill='white',fontsize=7,bold=FALSE,color='black') {
          grid::grid.roundrect(x,y,w,h,r=grid::unit(1.4,'mm'),gp=grid::gpar(fill=fill,col='grey65',lwd=.6))
          if(nzchar(text))grid::grid.text(text,x,y,gp=grid::gpar(fontfamily=STEP8_FONT$family,fontsize=fontsize,fontface=if(bold)'bold' else 'plain',col=color))
        }
        text<-function(label,x,y,size=7,bold=FALSE)grid::grid.text(label,x,y,gp=grid::gpar(fontfamily=STEP8_FONT$family,fontsize=size,fontface=if(bold)'bold' else 'plain'))
        arrow<-function(x,y,xend,yend,dashed=FALSE)grid::grid.lines(c(x,xend),c(y,yend),arrow=grid::arrow(length=grid::unit(1.3,'mm')),gp=grid::gpar(col='grey45',lwd=.7,lty=if(dashed)2 else 1))
        text('Study design',.5,.965,9,TRUE)
        box(.085,.54,.15,.72,paste0('Curated iCCA\nexpression\nresource\n\n',format(discovery+projection,big.mark=','),' samples\n',nrow(cohort),' studies'),fill='grey97',fontsize=7.5,bold=TRUE)
        text('RNA-seq discovery',.4,.895,8,TRUE)
        text('External microarray projection',.43,.445,8,TRUE)
        grid::grid.lines(c(.16,.175,.175),c(.54,.54,.74),gp=grid::gpar(col='grey50',lwd=.7))
        grid::grid.lines(c(.175,.175),c(.54,.28),gp=grid::gpar(col='grey50',lwd=.7))
        arrow(.175,.74,.185,.74);arrow(.175,.28,.185,.28)
        box(.25,.74,.13,.22,paste0('RNA-seq\nn = ',discovery,'\n',cohort[AnalysisRole=='Discovery RNA-seq',.N],' studies'),fill=apal[['Discovery RNA-seq']],bold=TRUE,color='white')
        box(.42,.74,.15,.22,'Expression\nharmonization\nGene filtering\nStudy adjustment')
        box(.60,.74,.17,.22,'NMF + k-means\nStability\nCross-study sensitivity')
        arrow(.315,.74,.345,.74);arrow(.495,.74,.515,.74);arrow(.685,.74,.7,.74)
        box(.84,.74,.28,.30,'',fill='grey98')
        text(paste0('Principal states: n = ',sum(states$N)),.84,.85,7.5,TRUE)
        for(i in seq_len(nrow(states)))box(.75+(i-1)*.09,.75,.08,.12,paste0(states$FinalSubtype[i],'\nn = ',states$N[i]),fill=pal[states$FinalSubtype[i]],bold=TRUE,color=if(i==1)'white' else 'black')
        grid::grid.lines(c(.71,.97),c(.675,.675),gp=grid::gpar(col='grey75',lwd=.5))
        text(paste0('Intermediate states: n = ',counts[SubtypeClass=='Intermediate',sum(N)],'\nSupplementary continuity'),.84,.63,6.5)
        # Label source is explicitly the principal-state cards, not intermediates.
        grid::grid.lines(c(.715,.715,.65),c(.75,.51,.51),gp=grid::gpar(col='grey45',lwd=.7,lty=2))
        arrow(.65,.51,.65,.385,TRUE)
        text('Trained on principal states',.83,.51,6.5)
        box(.25,.28,.13,.22,paste0('Microarray\nn = ',projection,'\n',cohort[AnalysisRole=='Microarray projection',.N],' studies'),fill=apal[['Microarray projection']],bold=TRUE,color='black')
        box(.44,.28,.19,.22,'Gene matching\nCohort-wise\nstandardization')
        box(.65,.28,.17,.22,'Fixed 123-gene\nlinear SVM')
        box(.86,.28,.19,.22,'External projection\nProgram preservation')
        arrow(.315,.28,.345,.28);arrow(.535,.28,.565,.28);arrow(.735,.28,.765,.28)
        text('Complementary evidence',.13,.10,7,TRUE)
        text('Genomics / proteomics / MOFA',.43,.10,6.5)
        text('Clinical / survival',.67,.10,6.5)
        text('Separate scRNA cohorts',.87,.10,6.5)
      })
    }
    if (id=='B') {
      d <- read_source(paste0(T1,'step1_method_selection_metrics_long.csv'),c('Method','K','Value','Facet'))
      d <- d[Method %in% c('NMF','KMeans','K-means')]
      return(ggplot(d,aes(K,Value))+geom_line(color='grey40')+geom_point(size=1.4)+facet_wrap(~Facet,scales='free_y')+
        labs(title='Retained-method rank diagnostics',caption='Metrics use different scales; consensus/hierarchical degeneration is reported in S2.',x='Rank / cluster number',y=NULL)+theme_step8())
    }
    if (id=='C') {
      marker <- read_source(paste0(T1,'step1_bulk_principal_state_markers_for_step8.csv'),c('Gene','Subtype'))
      master <- master_table();master[,FinalSubtype:=factor(FinalSubtype,levels=names(load_subtype_palette()$colors)[1:3])];setorder(master,FinalSubtype,Study,Sample)
      path <- file.path(STEP8_ROOT,'1-Clustering/data/combined_expr_after_correction.rds')
      expr <- tryCatch(readRDS(path),error=function(e){ env<-new.env();load(path,envir=env); if(!exists('corrected_expr_limma',env,inherits=FALSE))stop('Expression object name not found');env$corrected_expr_limma })
      genes <- unique(marker$Gene);stopifnot(all(genes %in% rownames(expr)),all(master$Sample %in% colnames(expr)))
      mat <- as.matrix(expr[genes,as.character(master$Sample),drop=FALSE]);s <- apply(mat,1,sd)
      if(any(!is.finite(s)|s==0))stop('Constant/nonfinite marker row; cannot row-z-score')
      mat <- t(scale(t(mat)))
      return(function() {
      study_palette <- load_study_palette()$colors
      ann <- ComplexHeatmap::HeatmapAnnotation(Subtype=as.character(master$FinalSubtype),Study=master$Study,
        col=list(Subtype=load_subtype_palette()$colors,Study=study_palette),annotation_name_gp=grid::gpar(fontsize=7))
      h <- ComplexHeatmap::Heatmap(mat,name='Row z-score',cluster_rows=FALSE,cluster_columns=FALSE,show_column_names=FALSE,
        row_names_gp=grid::gpar(fontsize=7,fontfamily=STEP8_FONT$family),top_annotation=ann,
        column_split=master$FinalSubtype,column_title_gp=grid::gpar(fontsize=8),
        col=circlize::colorRamp2(c(-3,0,3),c(SCALE_EFFECT_LOW,SCALE_EFFECT_MID,SCALE_EFFECT_HIGH)),
        column_title='Corrected discovery markers; main-state samples; row z-score')
      ComplexHeatmap::draw(h)
      })
    }
    if (id=='D') {
      d <- read_source('governance/data/pca_coordinates_top2000.csv');d <- merge(d[,.(Sample,PC1,PC2)],master_table()[,.(Sample,FinalSubtype,Study)],by='Sample')
      return(ggplot(d,aes(PC1,PC2,color=FinalSubtype))+geom_point(size=1,alpha=.7)+subtype_scale()+
        labs(title='Main-state expression structure',subtitle='Fixed top-2,000-gene PCA; axis variance is not supplied in the coordinate table',x='PC1',y='PC2')+theme_step8())
    }
    if (id=='E') {
      d <- read_source(paste0(T1,'step1_main_subtype_by_study.csv'));d[,FinalSubtype:=Ensemble]
      return(ggplot(d,aes(Study,Proportion,fill=FinalSubtype))+geom_col()+subtype_scale('fill')+coord_flip()+
        geom_text(aes(label=paste0('n=',N)),position=position_stack(vjust=.5),size=2.3)+
        scale_y_continuous(labels=scales::label_percent())+labs(title='Principal states across studies',x=NULL,y='Within-study fraction')+theme_step8())
    }
    if (id=='F') {
      m <- master_table();d <- melt(m,id.vars=c('Sample','Study','FinalSubtype'),measure.vars=c('Subclass','Subclass-add'),variable.name='SourceField',value.name='OriginalLabel')
      d <- d[!is.na(OriginalLabel)& !OriginalLabel %in% c('NA','w/o','')];d[,Taxonomy:=paste(Study,SourceField,sep=' / ')]
      d <- d[,.N,by=.(Taxonomy,OriginalLabel,FinalSubtype)];d[,Fraction:=N/sum(N),by=.(Taxonomy,OriginalLabel)]
      return(plot_tile(d,'FinalSubtype','OriginalLabel','Fraction','Study-specific original classifications',FALSE)+
        geom_text(aes(label=paste0(N,'\n',round(Fraction*100),'%')),size=2.3)+facet_wrap(~Taxonomy,scales='free_y')+
        labs(caption='Available-case row fractions. Original systems remain separate; no pooled “Sia” taxonomy.'))
    }
  }
  if (figure=='Figure 2') {
    if(id=='A') {
      d <- read_source(paste0(T1,'step1_main_only_study_adjusted_deg.csv'));markers <- read_source(paste0(T1,'step1_bulk_principal_state_markers_for_step8.csv'))
      selected <- unique(markers[,head(.SD,8),by=Subtype]$Gene);d<-d[Gene %in% selected]
      return(plot_tile(d,'Subtype','Gene','logFC_subtype_minus_others','Corrected marker effects')+
        geom_text(aes(label=fdr_label(FDR)),size=2)+labs(caption='First eight genes per subtype from the frozen marker table; Study-adjusted subtype-minus-others effects.'))
    }
    if(id=='B') {
      d<-read_source(paste0(T1,'step1_ranked_gsea_main_only_study_adjusted.csv'));d<-d[Collection=='Hallmark']
      d[,DisplayPathway:=tools::toTitleCase(tolower(gsub('_',' ',sub('^HALLMARK_','',pathway))))]
      return(plot_tile(d,'Subtype','DisplayPathway','NES','Hallmark pathway enrichment')+geom_text(aes(label=ifelse(padj<.05,fdr_label(padj),'')),size=2)+labs(caption='All available Hallmark pathways; no per-row rescaling; q from the locked GSEA family.'))
    }
    if(id=='C') {
      d<-read_source('2-TranscriptomicCharacterization/data/metabolic_main_subtypes_scores_long.csv')
      return(box_plot(d,'Score','Signature','Metabolic program distributions')+labs(caption='Patient/sample-level descriptive distributions; no new tests. Other program scores are supplementary.'))
    }
    if(id=='D') {
      d<-read_source(paste0(T1,'step1_dorothea_top_tf_by_subtype.csv'));cov<-read_source(paste0(T1,'step1_dorothea_regulon_coverage.csv'));d<-merge(d,cov,by.x='TF',by.y='source',all.x=TRUE)
      return(ggplot(d,aes(Subtype,TF,color=ActivityDifferenceSubtypeMinusOthers,size=TargetsAvailable))+geom_point()+
        scale_colour_gradient2(low=SCALE_EFFECT_LOW,mid=SCALE_EFFECT_MID,high=SCALE_EFFECT_HIGH)+labs(title='Inferred TF activity',color='Adjusted effect',size='Targets available',x=NULL,y=NULL)+theme_step8())
    }
    if(id=='E') {
      d<-read_source(paste0(T2,'step2_wgcna_module_subtype_associations.csv'));d<-d[grepl('Validation',Dataset)]
      if(!nrow(d))stop('No held-out WGCNA associations')
      return(plot_tile(d,'Subtype','Module','EffectSD','Held-out WGCNA module associations')+geom_text(aes(label=fdr_label(FDR_within_dataset_subtype)),size=2)+labs(caption='Supportive co-expression evidence; no C2 biological suffix. Preservation is shown separately in S5.'))
    }
    if(id=='F') {
      d<-read_source(paste0(T1,'step1_tme_concordance_scores_long.csv'));d[,Standardized:=as.numeric(scale(Value)),by=.(Method,Score)]
      a<-d[,.(MeanStandardized=mean(Standardized,na.rm=TRUE),N=.N),by=.(Method,Score,Subtype)]
      return(plot_tile(a,'Subtype','Score','MeanStandardized','Microenvironment score profiles')+facet_wrap(~Method,scales='free_y')+
        labs(caption='Descriptive within-score standardization; methods have different estimands. TIDE is separate in S6.'))
    }
  }
  if(figure=='Figure 4') {
    if(id=='A') {
      d<-read_source(paste0(T4,'step4_driver_snv_indel_sample_long.csv'));d<-d[SubtypeClass=='Main'];setorder(d,FinalSubtype,Study,Sample)
      samples<-unique(d$Sample);d[,Sample:=factor(Sample,levels=samples)];d[,Call:=ifelse(is.na(Mutated),'Unknown',ifelse(Mutated,'Variant recorded','No variant recorded'))]
      return(ggplot(d,aes(Sample,Gene,fill=Call))+geom_tile()+facet_grid(~FinalSubtype,scales='free_x',space='free_x')+
        scale_fill_manual(values=c('Variant recorded'='grey20','No variant recorded'='grey95','Unknown'=MISSING_UNKNOWN_FILL))+
        labs(title='Prespecified driver variants',caption=paste('Assay main-state samples:',length(samples),'; no variant recorded is not proof of gene-level callability.'),x=NULL,y=NULL)+theme_step8()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()))
    }
    if(id=='B') {
      d<-read_source(paste0(T4,'step4_driver_snv_indel_study_stratified_tests.csv'))
      return(forest_plot(d,'Gene','CMHOddsRatio','CMHLower95','CMHUpper95','Subtype',TRUE,1,'Study-stratified alteration associations')+labs(x='CMH odds ratio (95% CI)',caption='Prespecified comparisons; within-subtype FDR remains in source data.'))
    }
    if(id=='C') {
      d<-read_source(paste0(T4,'step4_matched_rna_protein_subtype_contrasts_long.csv'),c('RNAModeratedT','ProteinModeratedT'))
      summary<-read_source(paste0(T4,'step4_matched_rna_protein_subtype_concordance.csv'))
      labels<-setNames(paste0(summary$Subtype,'; genes=',summary$Genes,'; rho=',round(summary$SpearmanRho,3)),summary$Subtype)
      return(ggplot(d,aes(RNAModeratedT,ProteinModeratedT))+geom_bin_2d(bins=45)+scale_fill_gradient(low='grey95',high='grey20')+
        facet_wrap(~Subtype,labeller=labeller(Subtype=labels))+labs(title='RNA–protein contrast concordance',x='RNA moderated t',y='Protein moderated t',fill='Genes',caption='Matched-gene contrast statistics, not logFC; no new correlation test.')+theme_step8())
    }
    if(id=='D') {
      d<-read_source(paste0(T4,'step4_matched_rna_protein_pathway_concordance.csv'));d<-d[Collection=='Hallmark'];d[,Significance:=paste0('RNA q<.05: ',padj_RNA<.05,'; protein q<.05: ',padj_Protein<.05)]
      return(ggplot(d,aes(NES_RNA,NES_Protein,shape=Significance))+geom_hline(yintercept=0,color='grey70')+geom_vline(xintercept=0,color='grey70')+geom_point(size=1.5)+facet_wrap(~Subtype)+
        labs(title='RNA–protein pathway direction',x='RNA NES',y='Protein NES',caption='Quadrants encode direction concordance; symbols encode the two locked q-value criteria.')+theme_step8())
    }
    if(id=='E') {
      d<-read_source(paste0(T4,'step4_mofa_variance_explained_per_factor.csv'))
      return(plot_tile(d,'View','Factor','VarianceExplainedPercent','MOFA explained variance by factor and view',FALSE)+labs(caption='Fit population: 204 triple-matched tumors. One seed; no multi-seed stability claim.'))
    }
    if(id=='F') {
      d<-read_source(paste0(T4,'step4_multiomic_cluster_comparison.csv'))
      return(forest_plot(d,'Comparison','ARI','BootstrapLower95','BootstrapUpper95',title='Agreement with locked principal states')+labs(x='ARI (locked bootstrap 95% interval)',caption='182 principal-state matched tumors; MOFA and integrated PCA provide different supporting evidence.'))
    }
  }
  if(figure=='Figure 5') {
    if(id=='A') {
      d<-data.table(x=c(.1,.35,.65,.9),y=.55,text=c('Outer training data','Inner folds\nFeature selection\nScaling / tuning','Refit on outer training\nFrozen preprocessing','Outer-test predictions\nOOF aggregation'))
      folds<-read_source(paste0(T5,'step5_nested_cv_fold_definitions.csv'));nfolds<-uniqueN(folds[,.(Repeat,Fold)])
      return(ggplot(d,aes(x,y))+geom_label(aes(label=text),size=2.5,label.padding=grid::unit(3,'mm'))+
        geom_segment(data=data.table(x=c(.19,.46,.76),xend=c(.24,.54,.8),y=.55),aes(x=x,xend=xend,y=y,yend=y),arrow=grid::arrow(length=grid::unit(2,'mm')))+
        coord_cartesian(xlim=c(-.02,1.04),ylim=c(.2,.85),clip='off')+theme_void()+labs(title='Leakage-controlled label reproducibility',subtitle=paste(nfolds,'outer repeat-fold pairs; all selection and tuning stay inside training folds')))
    }
    if(id=='B') {
      d<-read_source(paste0(T5,'step5_nested_cv_confusion_matrix.csv'));d<-d[Panel=='Operational123']
      return(plot_tile(d,'Prediction','Truth','RowFraction','Operational123 out-of-fold assignments',FALSE)+
        geom_text(aes(label=paste0(Freq,'\n',round(RowFraction*100,1),'%')),size=2.5)+labs(caption='497 principal-state samples; row fractions; mean of five outer-test probabilities.'))
    }
    if(id=='C') {
      d<-read_source(paste0(T5,'step5_classifier_performance_with_bootstrap_ci.csv'));d<-d[Metric %in% c('BalancedAccuracy','MacroF1','MulticlassBrier')]
      return(forest_plot(d,'Panel','Estimate','Lower95','Upper95','Metric',title='Same-fold molecular panel comparison')+labs(caption='Locked patient-bootstrap intervals; repeated folds are not independent patients.'))
    }
    if(id=='D') {
      d<-read_source(paste0(T5,'step5_nested_cv_calibration_bins.csv'));d<-d[Panel=='Operational123']
      return(ggplot(d,aes(Predicted,Observed,color=Class))+geom_abline(slope=1,intercept=0,color='grey60',linetype=2)+geom_line()+geom_point(aes(size=N))+subtype_scale()+
        facet_wrap(~Class)+coord_cartesian(xlim=c(0,1),ylim=c(0,1))+labs(title='OOF probability calibration',x='Mean predicted probability',y='Observed fraction',caption='Locked bins and counts; probabilities are not deployment-calibrated. No new ECE or threshold.')+theme_step8())
    }
    if(id=='E') {
      d<-read_source(paste0(T5,'step5_loso_per_study_performance.csv'));d[,Label:=paste0(Study,' (n=',N,')')]
      return(forest_plot(d,'Label','BalancedAccuracy',title='Leave-one-study-out label reproducibility')+labs(caption='Point estimates only: source contains no CIs. Held-out cohort centering is transductive.',x='Balanced accuracy'))
    }
    if(id=='F') {
      d<-read_source(paste0(T5,'step5_external_projection_preservation.csv'));d<-melt(d,id.vars=c('Platform','N'),measure.vars=c('PredictedC1','PredictedC2','PredictedC3'),variable.name='Subtype',value.name='Count');d[,FinalSubtype:=sub('Predicted','',Subtype)];d[,Fraction:=Count/N];d[,Track:=ifelse(Platform=='DepMap/CCLE','Cell lines','Microarray samples')]
      return(ggplot(d,aes(Platform,Fraction,fill=FinalSubtype))+geom_col()+subtype_scale('fill')+coord_flip()+facet_wrap(~Track,scales='free_y')+scale_y_continuous(labels=scales::label_percent())+
        labs(title='External subtype projection',x=NULL,y='Projected fraction',caption='No independent subtype truth. Samples and cell lines are separate; no low-margin cutoff.')+theme_step8())
    }
  }
  if(figure=='Figure 6') {
    if(id=='A')return(structure(list(source=paste0('6-ClinicalSurvival/publication/figures/','step6_km_main_states.pdf')),class='locked_pdf'))
    if(id=='B') {
      d<-read_source(paste0(T6,'step6_cox_time_varying_subtype_effects.csv'));d<-d[Contrast=='C3 vs C1']
      return(ggplot(d,aes(TimeMonths,HR))+geom_ribbon(aes(ymin=Lower95,ymax=Upper95),fill='grey85')+geom_line(color=load_subtype_palette()$colors['C3'])+geom_point()+geom_hline(yintercept=1,linetype=2)+scale_y_log10()+labs(title='Time-varying C3 survival association',x='Months',y='Adjusted HR (95% CI)',caption='331 patients / 153 events. Adverse early association; no constant lifetime hazard.')+theme_step8())
    }
    if(id=='C') {
      d<-read_source(paste0(T6,'step6_study_adjusted_clinical_phenotype_contrasts.csv'));d[,Contrast:=paste(Target,'vs',Reference)]
      nulls<-unique(d[,.(Endpoint,EffectScale)]);nulls[,Null:=ifelse(EffectScale=='Odds ratio',1,0)]
      return(forest_plot(d,'Contrast','Estimate','Lower95','Upper95','Endpoint',ref=NA_real_,title='Study-adjusted clinical phenotype contrasts')+geom_vline(data=nulls,aes(xintercept=Null),color='grey45',linetype=2)+labs(caption='Endpoint-specific available-case N and effect scales are supplied by the locked table; observed duct counts are separate in S12.'))
    }
    if(id=='D') {
      d<-read_source(paste0(T6,'step6_incremental_prognostic_performance.csv'));d<-d[grepl('Delta|delta|Increment',Comparison)]
      if(!nrow(d))stop('Locked delta rows not found')
      return(forest_plot(d,'Comparison','Estimate','Lower95','Upper95','Metric',title='Incremental prognostic performance')+labs(caption='331 patients / 153 events. Positive delta improves discrimination; negative delta improves Brier/IBS.'))
    }
    if(id=='E') {
      d<-read_source(paste0(T5,'step5_c3_clinical_panel_performance.csv'));d<-d[Reference=='Absolute' & Metric %in% c('AUROC','AUPRC')]
      return(forest_plot(d,'Comparison','Estimate','Lower95','Upper95','Metric',title='Compact C3 detection feasibility')+labs(caption='Same 341 patients and repeated 5x5 outer folds; discovery-label reproducibility, not independent diagnostic accuracy. Duct incremental null is in S12.'))
    }
    if(id=='F') {
      d<-read_source(paste0(T5,'step5_oncopredict_candidate_evidence_matrix.csv'));d<-d[ScoreType=='SampleCentered'];d[,Label:=paste(Drug,Contrast)]
      return(forest_plot(d,'Label','Estimate','CI_low','CI_high',title='Prespecified therapeutic hypotheses')+labs(x='Study-adjusted predicted-score difference',caption='497 main-state tumors. Lower scores mean greater predicted sensitivity; no clinical efficacy claim. Quality and method sensitivity are in S14.'))
    }
  }
  if(figure=='Supplement')return(build_supplement(id))
  stop('No renderer for ',key)
}
# Exports carry their own source hashes and remain Rendered until real visual QA.
render_subfigure <- function(figure,id) {
  row<-preflight_subfigure(figure,id);plot<-build_subfigure(figure,id)
  if (inherits(plot,'locked_pdf')) {
    src<-file.path(STEP8_ROOT,plot$source);out<-file.path(STEP8_ROOT,row$TargetOutput)
    dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
    # The locked KM PDF has a blank device-opening page followed by its graph.
    # Remove only verified blank pages; do not refit or alter the graph.
    info <- system2('pdfinfo',shQuote(src),stdout=TRUE)
    pages <- as.integer(sub('^Pages:[[:space:]]+','',grep('^Pages:',info,value=TRUE)))
    if(pages==1) {
      if(!file.copy(src,out,overwrite=TRUE))stop('Cannot copy locked source PDF')
    } else {
      for(page in seq_len(pages-1)) {
        txt <- system2('pdftotext',c('-f',page,'-l',page,shQuote(src),'-'),stdout=TRUE)
        if(nzchar(gsub('[[:space:]]','',paste(txt,collapse=''))))stop('Locked PDF contains multiple nonblank pages; do not silently discard')
      }
      tmpdir <- tempfile();dir.create(tmpdir)
      on.exit(unlink(tmpdir,recursive=TRUE),add=TRUE)
      pattern <- file.path(tmpdir,'page-%d.pdf')
      status <- system2('pdfseparate',c('-f',pages,'-l',pages,shQuote(src),shQuote(pattern)))
      if(status!=0 || !file.copy(file.path(tmpdir,paste0('page-',pages,'.pdf')),out,overwrite=TRUE))stop('KM PDF extraction failed')
    }
    # Preserve existing physical size; validator can inspect it, not rescale blindly.
    width<-height<-NA_real_
  } else {
    width<-as.numeric(row$WidthMM);height<-as.numeric(row$HeightMM)
    save_panel(plot,row$TargetOutput,width,height)
  }
  m<-read_source('8-Visualization/tables/figure_manifest.csv');i<-which(m$Figure==figure&m$Panel==id)
  sources<-strsplit(row$SourceTable,';',fixed=TRUE)[[1]]
  hashes<-paste(paste0(sources,'=',unname(tools::md5sum(file.path(STEP8_ROOT,sources)))),collapse=';')
  m[i,`:=`(Status='Rendered',SourceMD5=hashes,QAFile='',RenderedAt=format(Sys.time(),tz='UTC',usetz=TRUE))]
  if(inherits(plot,'locked_pdf'))m[i,`:=`(WidthMM=NA_real_,HeightMM=NA_real_)]
  tmp<-tempfile(tmpdir=file.path(STEP8_ROOT,'8-Visualization/tables'),fileext='.csv')
  fwrite(m,tmp,na='');if(!file.rename(tmp,file.path(STEP8_ROOT,'8-Visualization/tables/figure_manifest.csv')))stop('Manifest write failed')
  message('Rendered standalone PDF: ',row$TargetOutput)
  invisible(file.path(STEP8_ROOT,row$TargetOutput))
}
render_figure <- function(figure) {
  m<-read_source('8-Visualization/tables/figure_manifest.csv');ids<-m[Figure==figure&Status!='Blocked',Panel]
  if(!length(ids))stop('No unblocked subfigures: ',figure)
  lapply(ids,function(id)render_subfigure(figure,id))
}

build_supplement <- function(id) {
  if(id=='S1A') {
    d<-read_source('0-Data/qc/step0_cohort_inclusion_manifest.csv')
    return(ggplot(d,aes(Study,FinalSamples))+geom_col(fill='grey45')+facet_wrap(~AnalysisRole,scales='free')+coord_flip()+labs(title='Curated cohort inclusion',x=NULL,y='Included expression samples')+theme_step8())
  }
  if(id=='S1B') {
    d<-read_source('0-Data/qc/step0_assay_availability_matrix.csv')
    return(plot_tile(d,'Assay','Study','CoveragePercent','Assay availability relative to expression',FALSE)+geom_text(aes(label=DisplayLabel),size=2.3)+labs(caption='Cell labels give matched/curated-expression counts; assays have different denominators.'))
  }
  if(id=='S1C') {
    d<-read_source('0-Data/qc/step0_expression_overview.csv')
    return(ggplot(d,aes(RetainedPercentOfUniqueSource,Dataset))+geom_point(size=2)+labs(title='Processed expression feature retention',x='Retained unique-source features (%)',y=NULL)+theme_step8())
  }
  if(id=='S1D') {
    d<-read_source('0-Data/qc/step0_protein_pca_scores.csv')
    return(ggplot(d,aes(PC1,PC2,color=Study))+geom_point(size=.7)+facet_wrap(~Stage)+scale_colour_manual(values=load_study_palette()$colors)+labs(title='Protein batch-adjustment diagnostics',caption='All assay samples; matched-expression status is recorded in the source. No visual claim of complete batch removal.')+theme_step8())
  }
  if(id=='S2A') {
    d<-read_source(paste0(T1,'step1_kmeans_rank_sensitivity.csv'))
    return(ggplot(d,aes(K,MeanSilhouette))+geom_line()+geom_point()+labs(title='K-means rank sensitivity',x='k',y='Mean silhouette')+theme_step8())
  }
  if(id=='S2B') {
    d<-read_source(paste0(T1,'step1_kmeans_feature_number_sensitivity.csv'))
    return(ggplot(d,aes(Features,ARIvsLockedKMeans))+geom_line()+geom_point()+labs(title='Feature-number sensitivity',x='Features',y='ARI vs locked k-means')+theme_step8())
  }
  if(id=='S2C') {
    d<-read_source(paste0(T1,'step1_kmeans_k3_seed_sensitivity.csv'))
    return(ggplot(d,aes(Seed,ARIvsHistoricalKMeans))+geom_point(size=1)+labs(title='K-means seed sensitivity',x='Seed',y='ARI vs locked assignments')+theme_step8())
  }
  if(id=='S2D') {
    d<-read_source(paste0(T1,'step1_kmeans_leave_one_study_out.csv'))
    return(forest_plot(d,'HeldOutStudy','ARIAllVsHistoricalKMeans',title='Discovery leave-one-study-out sensitivity')+labs(x='All-sample ARI vs locked k-means'))
  }
  if(id=='S2E') {
    d<-read_source(paste0(T1,'step1_standard_consensus_cluster_sizes.csv'))
    return(ggplot(d,aes(factor(K),Fraction,fill=factor(Cluster)))+geom_col()+scale_fill_grey()+labs(title='Consensus degeneration: exclusion evidence',x='k',y='Cluster fraction',fill='Cluster',caption='Low PAC alone does not imply a useful partition; consensus was excluded.')+theme_step8())
  }
  if(id=='S2F') {
    d<-read_source(paste0(T1,'step1_nmf_kmeans_state_definition_full.csv'))
    return(ggplot(d,aes(factor(NMF),factor(KMeans),fill=Ensemble))+geom_tile()+geom_text(aes(label=N),size=3)+
      scale_fill_manual(values=load_subtype_palette()$colors)+labs(title='NMF / k-means state definition',x='NMF cluster',y='K-means cluster',fill='Governed state')+theme_step8())
  }
  if(id=='S3A') {
    d<-read_source('governance/data/subtype_marker_scores.csv');d<-melt(d,id.vars=c('Sample','FinalSubtype','SubtypeClass'),measure.vars=c('C1_marker_score','C2_marker_score','C3_marker_score'),variable.name='Program',value.name='Score')
    return(box_plot(d,'Score','Program','Main and intermediate program scores',FALSE)+labs(caption='Descriptive continuity; not a proven transition trajectory.'))
  }
  if(id=='S3B') {
    d<-read_source('governance/data/subtype_marker_scores.csv')
    return(ggplot(d,aes(C1_marker_score,C3_marker_score,color=FinalSubtype))+geom_point(size=.7,alpha=.7)+subtype_scale(main=FALSE)+labs(title='Continuous bulk program space',x='C1 marker score',y='C3 marker score',caption='All six governed states; two displayed score axes, not trajectory inference.')+theme_step8())
  }
  if(id %in% c('S4A','S4B','S4C')) {
    subtype<-c(S4A='C1',S4B='C2',S4C='C3')[[id]]
    d<-read_source(paste0(T1,'step1_deg_',subtype,'_minus_others_main_only_study_adjusted.csv'))
    # Preserve every row; only the display axis is capped, with disclosure.
    d[,LogQ:=pmin(-log10(pmax(FDR,.Machine$double.xmin)),50)]
    return(ggplot(d,aes(logFC_subtype_minus_others,LogQ))+geom_point(size=.3,alpha=.3,color='grey35')+
      labs(title=paste(subtype,'minus others: Study-adjusted DEG'),x='log2 fold change',y='-log10(q), display capped at 50',caption='All genes retained; numeric underflow zeros are plotted at the disclosed display cap.')+theme_step8())
  }
  if(id=='S4D') {
    d<-read_source(paste0(T1,'step1_cross_study_deg_meta_summary.csv'))
    return(ggplot(d,aes(MetaLogFC,DirectionConsistency))+geom_bin_2d(bins=35)+scale_fill_gradient(low='grey95',high='grey20')+facet_wrap(~Subtype)+labs(title='Cross-study DEG direction concordance',x='Meta logFC',y='Direction consistency',fill='Genes')+theme_step8())
  }
  if(id=='S4E') {
    d<-read_source(paste0(T1,'step1_dorothea_regulon_coverage.csv'))
    return(ggplot(d,aes(TargetsAvailable))+geom_histogram(bins=30,fill='grey50')+labs(title='Available TF regulon coverage',x='Targets available per regulon',y='Regulons')+theme_step8())
  }
  if(id=='S5A') {
    d<-read_source(paste0(T2,'step2_wgcna_soft_threshold_diagnostics.csv'))
    return(ggplot(d,aes(Power,SignedR2))+geom_line()+geom_point(aes(shape=Selected),size=2)+geom_hline(yintercept=.8,linetype=2)+labs(title='WGCNA soft-threshold diagnostics',x='Power',y='Signed scale-free R²',caption='2022Dong source n=223; connectivity-QC network n=220; held-out validation n=274.')+theme_step8())
  }
  if(id=='S5B') {
    d<-read_source(paste0(T2,'step2_wgcna_module_preservation.csv'))
    return(ggplot(d,aes(ZsummaryPreservation,Module))+geom_point(aes(size=ModuleSize))+labs(title='Held-out module preservation',x='Zsummary preservation',y=NULL,size='Module genes',caption='Locked 200-permutation full-module evaluation.')+theme_step8())
  }
  if(id=='S5C') {
    d<-read_source(paste0(T2,'step2_wgcna_subtype_random_effects_meta.csv'))
    return(forest_plot(d,'Module','RandomEffect','CI_low','CI_high','TargetSubtype',title='Symmetric WGCNA validation effects')+labs(caption='Supportive co-expression evidence; no biological suffix is inferred from a module alone.'))
  }
  if(id=='S5D') {
    d<-read_source(paste0(T2,'step2_wgcna_module_gene_set_overlap.csv'))
    return(plot_tile(d,'GeneSet','Module','Overlap','Module / gene-set overlap',FALSE)+theme(axis.text.x=element_text(angle=45,hjust=1)))
  }
  if(id=='S6A') {
    d<-read_source(paste0(T1,'step1_tme_concordance_scores_long.csv'));d[,FinalSubtype:=Subtype]
    return(box_plot(d,'Value','Score','Microenvironment score distributions')+labs(caption='Native score scales; different methods estimate different quantities. No added tests.'))
  }
  if(id=='S6B') {
    d<-read_source(paste0(T1,'step1_tme_score_spearman_concordance.csv'))
    return(plot_tile(d,'ScoreA','ScoreB','SpearmanRho','Locked inter-method score correlations')+theme(axis.text.x=element_text(angle=90,hjust=1)))
  }
  if(id=='S6C') {
    d<-read_source(paste0(T2,'step2_tide_sample_results.csv'));d[,Metric:='TIDE']
    return(box_plot(d,'TIDE','Metric','TIDE computational prediction')+labs(caption='497 principal-state tumors; algorithmic predictions, not observed immunotherapy outcomes.'))
  }
  if(id=='S6D') {
    d<-read_source(paste0(T2,'step2_tide_meta_analysis.csv'))
    return(forest_plot(d,'Contrast','RandomEffect','CI_low','CI_high','Metric',title='TIDE cross-study sensitivity'))
  }
  if(id=='S10A') {
    d<-read_source(paste0(T4,'step4_cnv_burden_2022dong.csv'));d<-d[SubtypeClass=='Main']
    return(box_plot(d,'MeanAbsoluteCNV','Study','CNV burden proxy')+labs(caption='2022Dong assay-specific main-state subset; not harmonized genomic burden.'))
  }
  if(id=='S10B') {
    d<-read_source(paste0(T4,'step4_nonsynonymous_variant_count_by_sample.csv'));d<-d[SubtypeClass=='Main']
    return(box_plot(d,'NonsynonymousVariantCount','Study','Nonsynonymous variant counts')+labs(caption='Assay/study-specific counts; not harmonized TMB.'))
  }
  if(id=='S10C') {
    d<-read_source(paste0(T4,'step4_protein_feature_coverage_by_study.csv'));d[,c('Coverage_2022Dong','Coverage_2023Deng'):=lapply(.SD,as.numeric),.SDcols=c('Coverage_2022Dong','Coverage_2023Deng')];d<-melt(d,id.vars='Gene',measure.vars=c('Coverage_2022Dong','Coverage_2023Deng'),variable.name='Study',value.name='Coverage')
    return(ggplot(d,aes(Coverage))+geom_histogram(bins=20,fill='grey50')+facet_wrap(~Study)+labs(title='Protein feature coverage',x='Within-study nonmissing fraction',y='Features')+theme_step8())
  }
  if(id=='S10D') {
    d<-read_source(paste0(T4,'step4_tcga_rppa_subtype_associations.csv'))
    return(ggplot(d,aes(Difference))+geom_histogram(bins=30,fill='grey50')+facet_wrap(~Subtype)+labs(title='TCGA RPPA effect distributions',x='Subtype-minus-others RPPA effect',y='Targets',caption='RPPA is separate from mass spectrometry; small assay-specific subtype groups, exploratory.')+theme_step8())
  }
  if(id=='S10E') {
    d<-read_source(paste0(T4,'step4_mofa_training_statistics.csv'))
    return(ggplot(d,aes(IterationIndex,elbo))+geom_line()+labs(title='MOFA optimization trace',x='Iteration',y='ELBO',caption='One executed seed. This trace does not establish multi-seed stability.')+theme_step8())
  }
  if(id=='S10F') {
    d<-read_source(paste0(T4,'step4_mofa_top_feature_weights.csv'))
    return(ggplot(d,aes(Weight,WithinViewFactorRank))+geom_point(size=.7)+facet_grid(View~Factor,scales='free_x')+labs(title='MOFA top feature-weight profiles',x='Weight (factor signs arbitrary)',y='Within-view factor rank')+theme_step8())
  }
  if(id=='S10G') {
    d<-read_source(paste0(T4,'step4_mofa_factor_scores.csv'));d<-d[SubtypeClass=='Main']
    return(ggplot(d,aes(Factor1,Factor2,color=FinalSubtype))+geom_point(size=1)+subtype_scale()+labs(title='MOFA principal-state factor map',caption='182 principal-state tumors from the 204-sample fitted model; factors 1 and 2 shown without post-hoc axis optimization.')+theme_step8())
  }
  if(id=='S11A') {
    d<-read_source(paste0(T5,'step5_nested_cv_calibration_metrics.csv'))
    return(ggplot(d,aes(Class,CalibrationSlope,shape=Panel))+geom_hline(yintercept=1,linetype=2)+geom_point(position=position_dodge(width=.4),size=2)+labs(title='OOF calibration slopes',x=NULL,y='Calibration slope',caption='Locked estimates only; no deployment calibration or new confidence intervals.')+theme_step8())
  }
  if(id=='S11B') {
    d<-read_source(paste0(T5,'step5_nested_cv_feature_stability.csv'))
    return(ggplot(d,aes(SelectionFraction))+geom_histogram(breaks=seq(0,1,.1),fill='grey50')+facet_wrap(~Panel)+labs(title='Fold-selected feature stability',x='Selection fraction across 25 outer folds',y='Genes',caption='Repeated 5x5 CV; zero-selection features absent from the source are not added.')+theme_step8())
  }
  if(id=='S11C') {
    d<-read_source(paste0(T5,'step5_external_projection_preservation.csv'))
    return(ggplot(d,aes(MatchedFeatures/TotalFeatures,Platform))+geom_point(size=2)+labs(title='External classifier gene coverage',x='Matched / 123 model genes',y=NULL)+theme_step8())
  }
  if(id=='S11D') {
    d<-read_source(paste0(T5,'step5_external_projection_sample_predictions.csv'))
    return(ggplot(d,aes(ProbabilityMargin))+geom_histogram(breaks=seq(0,1,.05),fill='grey50')+facet_wrap(~Platform)+labs(title='External probability-difference distributions',x='Largest minus second-largest uncalibrated probability',y='Samples / cell lines',caption='Continuous output only; no reliability threshold or independent accuracy interpretation.')+theme_step8())
  }
  if(id=='S11E') {
    d<-read_source(paste0(T5,'step5_nested_cv_predictions_by_repeat.csv'));d<-d[Panel=='Operational123']
    curves<-rbindlist(lapply(c('C1','C2','C3'),function(cl) {
      d[,.(Probability=get(paste0('Probability',cl)),Positive=Truth==cl,Repeat)][,{
        totals<-c(sum(Positive),sum(!Positive)); t<-data.table(p=Probability,pos=as.integer(Positive),neg=as.integer(!Positive))[,.(pos=sum(pos),neg=sum(neg)),by=p];setorder(t,-p)
        data.table(FPR=c(0,cumsum(t$neg)/totals[2]),TPR=c(0,cumsum(t$pos)/totals[1]),Class=cl)
      },by=Repeat]
    }))
    return(ggplot(curves,aes(FPR,TPR,group=Repeat))+geom_line(color='grey50',alpha=.6)+geom_abline(slope=1,intercept=0,linetype=2)+facet_wrap(~Class)+labs(title='Per-repeat empirical OOF ROC curves',caption='Curve coordinates are descriptive aggregations of locked predictions; no new AUC, CI, or test.')+theme_step8())
  }
  if(id=='S12A') {
    d<-read_source(paste0(T3,'step3_label_provenance_by_study.csv'))
    return(ggplot(d,aes(Study,N,fill=DuctObserved))+geom_col(position='dodge')+scale_fill_grey()+labs(title='Observed pathology label provenance',x=NULL,y='Observed tumors',caption='115 observed labels; terminology harmonization caveat retained.')+theme_step8())
  }
  if(id=='S12B') {
    d<-read_source(paste0(T3,'step3_classifier_performance_summary.csv'))
    return(forest_plot(d,'Analysis','Estimate','Lower95','Upper95','Metric',title='Observed-duct model OOF performance')+labs(caption='Locked CIs; observed labels are not pooled with full-model projections.'))
  }
  if(id=='S12C') {
    d<-read_source(paste0(T3,'step3_nested_cv_calibration_bins.csv'))
    return(ggplot(d,aes(Predicted,Observed))+geom_abline(slope=1,intercept=0,linetype=2)+geom_line()+geom_point(aes(size=N))+labs(title='Observed-duct OOF calibration',x='Predicted large-duct probability',y='Observed large-duct fraction')+theme_step8())
  }
  if(id=='S12D') {
    d<-read_source(paste0(T3,'step3_leave_one_study_out_metrics.csv'))
    return(forest_plot(d,'HeldOutStudy','Estimate',facet='Metric',title='Duct-model LOSO performance')+labs(caption='Point estimates only; locked source contains no intervals.'))
  }
  if(id %in% c('S12E','S12F')) {
    observed<-id=='S12E';file<-if(observed)'step3_observed_duct_distribution_by_subtype.csv' else 'step3_predicted_duct_distribution_by_subtype.csv'
    d<-read_source(paste0(T3,file));type<-if(observed)'DuctObserved' else 'DuctPredicted';frac<-if(observed)'WithinSubtypeFraction' else 'WithinSubtypeFraction'
    return(ggplot(d,aes(FinalSubtype,.data[[frac]],fill=.data[[type]]))+geom_col()+geom_text(aes(label=paste0('n=',N)),position=position_stack(vjust=.5),size=2.5)+scale_fill_grey()+labs(title=if(observed)'Observed duct morphology' else 'Predicted duct morphology',x=NULL,y='Within-subtype fraction',caption=unique(d$EvidenceType)[1])+theme_step8())
  }
  if(id=='S12G') {
    d<-read_source(paste0(T5,'step5_c3_clinical_panel_feature_stability.csv'))
    return(ggplot(d,aes(SelectionFrequency))+geom_histogram(breaks=seq(0,1,.1),fill='grey50')+labs(title='Compact C3 feature stability',x='Selection frequency',y='Genes')+theme_step8())
  }
  if(id=='S12H') {
    d<-read_source(paste0(T5,'step5_c3_incremental_value_summary.csv'));d<-d[grepl('Delta|delta',Comparison)]
    return(forest_plot(d,'Comparison','Estimate','Lower95','Upper95','Metric',title='Cross-fitted duct incremental value')+labs(caption='93 matched observed-label patients; negative result retained.'))
  }
  if(id=='S13A') {
    d<-read_source(paste0(T6,'step6_cox_proportional_hazards_diagnostics.csv'));d[,LogP:=-log10(pmax(PValue,.Machine$double.xmin))]
    return(ggplot(d,aes(LogP,Term))+geom_point()+facet_wrap(~Model,scales='free_y')+labs(title='Proportional-hazards diagnostic tests',x='-log10(unadjusted diagnostic p)',y=NULL,caption='Locked test summaries; no residual curve is reconstructed from this table.')+theme_step8())
  }
  if(id %in% c('S13B','S13C','S13D')) {
    f<-c(S13B='step6_cox_main_states_core.csv',S13C='step6_cox_main_states_extended_sensitivity.csv',S13D='step6_cox_all_states_exploratory.csv')[[id]]
    d<-read_source(paste0(T6,f));return(forest_plot(d,'Term','HR','Lower95','Upper95',log_axis=TRUE,ref=1,title=unique(d$Model)[1])+labs(x='HR (95% CI)',caption=paste('N=',unique(d$N)[1],'; events=',unique(d$Events)[1],'; average association; non-PH caveat retained.')))
  }
  if(id=='S13E') {
    d<-read_source(paste0(T6,'step6_c3_vs_c1_leave_one_study_out.csv'));return(forest_plot(d,'HeldOutStudy','HR','Lower95','Upper95',log_axis=TRUE,ref=1,title='C3 survival leave-one-study-out sensitivity'))
  }
  if(id=='S13F') {
    d<-read_source(paste0(T6,'step6_subtype_cox_random_effects_meta.csv'));d[,Contrast:=paste(Target,'vs',Reference)]
    return(forest_plot(d,'Contrast','HR','Lower95','Upper95','Scope',TRUE,1,'Study-specific Cox meta-analysis sensitivity'))
  }
  if(id=='S13G') {
    d<-read_source(paste0(T6,'step6_incremental_prognostic_calibration.csv'))
    return(ggplot(d,aes(HorizonMonths,CalibrationSlope,shape=Model))+geom_hline(yintercept=1,linetype=2)+geom_point(size=2)+geom_line(aes(group=Model))+labs(title='OOF prognostic calibration slopes',x='Horizon (months)',y='IPCW calibration slope',caption='Locked slopes only; six-month calibration is sparse. Not a reliability curve.')+theme_step8())
  }
  if(id=='S14A') {
    d<-read_source(paste0(T5,'step5_oncopredict_gdsc2_cv_quality.csv'));return(ggplot(d,aes(Spearman_rho))+geom_histogram(bins=25,fill='grey50')+labs(title='GDSC2 training-CV quality',x='Training OOF Spearman rho',y='Drug models')+theme_step8())
  }
  if(id=='S14B') {
    d<-read_source(paste0(T5,'step5_oncopredict_batch_method_sensitivity.csv'));d<-melt(d,id.vars='Drug',measure.vars=c('Primary_vs_pooled_standardize_Spearman','Primary_vs_pooled_eb_Spearman'),variable.name='Method',value.name='Correlation')
    return(ggplot(d,aes(Correlation))+geom_histogram(bins=25,fill='grey50')+facet_wrap(~Method)+labs(title='oncoPredict batch-method sensitivity',x='Prediction Spearman correlation',y='Drug models')+theme_step8())
  }
  if(id=='S14C') {
    d<-read_source(paste0(T5,'step5_oncopredict_candidate_evidence_matrix.csv'));d<-d[ScoreType=='SampleCentered']
    return(ggplot(d,aes(Fraction_C3_more_sensitive,Drug))+geom_point(aes(size=Studies))+facet_wrap(~Contrast)+labs(title='Prespecified candidate cross-study direction',x='Fraction of studies predicting greater C3 sensitivity',y=NULL,caption='Hypothesis-generating; repeated contrast rows do not add independent studies.')+theme_step8())
  }
  if(id=='S14E') {
    d<-fread(file.path(STEP8_ROOT,'5-ClassifierValidation/CCLE_drug_data_with_predicted_subtypes.csv'),sep=';',dec=',')
    candidates<-read_source(paste0(T5,'step5_oncopredict_candidate_evidence_matrix.csv'))
    keep<-unique(toupper(sub('_[0-9]+$','',candidates$Drug)))
    d<-d[Drug %in% keep]
    if(!nrow(d))stop('No prespecified candidate-level PRISM measurements')
    return(ggplot(d,aes(AUC,Drug))+geom_point(size=2)+facet_wrap(~DepMap_ID)+labs(title='Individual iCCA cell-line PRISM observations',x='Observed PRISM AUC',y=NULL,caption='Individual cell lines only; no subtype means, significance tests, or efficacy validation.')+theme_step8())
  }
  if(id=='S14D') {
    d<-read_source(paste0(T5,'step5_oncopredict_candidate_evidence_matrix.csv'));d<-unique(d[,.(Drug,N_training,Spearman_rho,ReliabilityFlag)])
    return(ggplot(d,aes(Spearman_rho,Drug))+geom_point(aes(size=N_training))+labs(title='Prespecified candidate model quality',x='GDSC2 training-CV Spearman rho',y=NULL,size='Training lines')+theme_step8())
  }
  stop('No supplementary renderer: ',id)
}
