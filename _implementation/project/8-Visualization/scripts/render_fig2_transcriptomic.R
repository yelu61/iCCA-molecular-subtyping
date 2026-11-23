# Figure 2 A/B: display-only transforms of locked Study-adjusted statistics.
# Source after render_bulk_figures.R. Other figures and Figure 2 C-F are unchanged.
if (!exists('build_subfigure')) stop('Source render_bulk_figures.R first')
if (!exists('build_subfigure_before_fig2')) build_subfigure_before_fig2 <- build_subfigure
FIG2_DIR <- '8-Visualization/refinements/fig2_transcriptomic_20261026_ora_split_v6'

fig2_pathway_plan <- function() {
  # Display selection, not a new gene-set test. All three contrasts remain in B.
  data.table::data.table(
    Group=c(rep('Metabolism',5),rep('Immune and extracellular programs',7),
            rep('Cell cycle and growth',5),rep('Cell state and microenvironment',5)),
    ShortName=c('XENOBIOTIC_METABOLISM','BILE_ACID_METABOLISM','FATTY_ACID_METABOLISM',
      'OXIDATIVE_PHOSPHORYLATION','PEROXISOME','COAGULATION','COMPLEMENT',
      'INTERFERON_ALPHA_RESPONSE','INTERFERON_GAMMA_RESPONSE','INFLAMMATORY_RESPONSE',
      'TNFA_SIGNALING_VIA_NFKB','IL6_JAK_STAT3_SIGNALING','G2M_CHECKPOINT','E2F_TARGETS',
      'MYC_TARGETS_V1','MTORC1_SIGNALING','DNA_REPAIR','APICAL_JUNCTION',
      'EPITHELIAL_MESENCHYMAL_TRANSITION','TGF_BETA_SIGNALING','ANGIOGENESIS','HYPOXIA'))[
        , `:=`(pathway=paste0('HALLMARK_',ShortName), DisplayOrder=.I)]
}

prepare_fig2_data <- function(write_tables=TRUE) {
  deg <- read_source(paste0(T1,'step1_main_only_study_adjusted_deg.csv'))
  all_gsea <- read_source(paste0(T1,'step1_ranked_gsea_main_only_study_adjusted.csv'))
  gsea <- all_gsea[Collection=='Hallmark']
  ora <- read_source(paste0(T1,'step1_thresholded_ora_main_only_study_adjusted.csv'))[Collection=='GO_BP']
  ora[, `:=`(pathway=ID,Side=ifelse(Direction=='Up','Higher','Lower'),padj=p.adjust,HitGenes=geneID,Method='ORA')]
  stopifnot(!anyDuplicated(ora[,.(Subtype,Side,pathway)]))
  markers <- read_source(paste0(T1,'step1_bulk_principal_state_markers_for_step8.csv'))
  stopifnot(!anyDuplicated(deg[,.(Subtype,Gene)]),!anyDuplicated(gsea[,.(Subtype,pathway)]))
  stopifnot(all(is.finite(deg$logFC_subtype_minus_others)),all(is.finite(deg$FDR)),
    all(deg$FDR>=0 & deg$FDR<=1),all(is.finite(gsea$NES)),all(gsea$padj>=0 & gsea$padj<=1))
  main <- master_table()
  ns <- main[,.N,by=FinalSubtype]
  for (cl in c('C1','C2','C3')) {
    stopifnot(all(deg[Subtype==cl,NSubtype]==ns[FinalSubtype==cl,N]),
      all(deg[Subtype==cl,NTotal]==nrow(main)),all(deg[Subtype==cl,StudyAdjusted]))
  }
  # Match the owning Step1 volcano and thresholded ORA contract: FDR<.05 and |log2FC|>=1.
  # B GSEA uses full moderated-t ranks; A ORA-hit labels obey the DEG cutoff.
  deg[, Direction:=fifelse(FDR<.05 & logFC_subtype_minus_others>=1,'Higher',
                   fifelse(FDR<.05 & logFC_subtype_minus_others<=-1,'Lower','Not selected'))]
  deg[, NegLog10FDR:=-log10(pmax(FDR,.Machine$double.xmin))]
  plan <- fig2_pathway_plan()
  b <- merge(plan,gsea,by='pathway',all.x=TRUE,sort=FALSE)
  stopifnot(nrow(b)==nrow(plan)*3L, !anyNA(b$NES),
    all(b[,.N,by=pathway]$N==3L))
  b[, LeadingEdgeN:=lengths(strsplit(leadingEdge,';',fixed=TRUE))]
  b[, DisplayPathway:=tools::toTitleCase(tolower(gsub('_',' ',ShortName))) ]
  display_names <- c(G2M_CHECKPOINT='G2M checkpoint',E2F_TARGETS='E2F targets',
    MYC_TARGETS_V1='MYC targets V1',MTORC1_SIGNALING='mTORC1 signaling',DNA_REPAIR='DNA repair',
    TNFA_SIGNALING_VIA_NFKB='TNF-alpha / NF-kB signaling',
    IL6_JAK_STAT3_SIGNALING='IL6 / JAK / STAT3 signaling',TGF_BETA_SIGNALING='TGF-beta signaling')
  for (key in names(display_names)) b[ShortName==key,DisplayPathway:=display_names[[key]]]
  b[, Significant:=padj<.05]
  setorder(b,DisplayOrder,Subtype)

  # Display-only term list; any number per side is supported by the layout.
  aplan <- read_source('8-Visualization/scripts/fig2_gobp_display_plan.csv',
    c('Subtype','AOrder','Side','ShortName'))
  stopifnot(!anyDuplicated(aplan[,.(Subtype,ShortName)]),
    all(aplan$Side %in% c('Lower','Higher')),
    all(c('C1','C2','C3') %in% aplan$Subtype))
  aplan[,pathway:=paste0('GOBP_',ShortName)]
  a <- merge(aplan,ora,by=c('Subtype','Side','pathway'),sort=FALSE)
  a[, `:=`(Significant=padj<.05,
    DisplayPathway=paste0(toupper(substr(gsub('_',' ',tolower(ShortName)),1,1)),
                         substr(gsub('_',' ',tolower(ShortName)),2,1000)))]
  a[ShortName=='PROTON_MOTIVE_FORCE_DRIVEN_ATP_SYNTHESIS',
    DisplayPathway:='Proton motive force-driven ATP synthesis']
  a[ShortName=='POSITIVE_REGULATION_OF_WNT_SIGNALING_PATHWAY',
    DisplayPathway:='Positive regulation of WNT signaling']
  setorder(a,Subtype,AOrder)
  stopifnot(nrow(a)==nrow(aplan),all(a$Significant),
    all(a$Method=='ORA'))
  labels <- rbindlist(lapply(seq_len(nrow(a)),function(i) {
    row <- a[i]; cl <- row$Subtype
    genes <- strsplit(row$HitGenes,'/',fixed=TRUE)[[1]]
    candidates <- copy(deg[Subtype==cl & Gene %in% genes & FDR<.05 &
                       abs(logFC_subtype_minus_others)>=1 &
                       sign(logFC_subtype_minus_others)==ifelse(row$Side=='Higher',1,-1)])
    stable <- markers[Subtype==cl & DirectionConsistency>=.8 & MetaFDR<.05,Gene]
    candidates[, `:=`(StableMarker=Gene %in% stable,AbsT=abs(t))]
    setorder(candidates,-StableMarker,-AbsT,Gene)
    if(nrow(candidates)<2L) stop('Fewer than two eligible ORA-hit labels: ',row$pathway)
    chosen <- candidates[1:2]
    chosen[, `:=`(pathway=row$pathway,AOrder=row$AOrder,Side=row$Side,WithinTermRank=1:2,
      SelectionRule='ORA-hit membership + DEG FDR<0.05 + abs(log2FC)>=1 + matching Up/Down direction; prefer frozen marker with direction consistency>=0.8 and meta-FDR<0.05, then abs(t), gene')]
    chosen
  }))
  # Compact first-row volcano: three unique ORA-hit labels in each direction.
  # Two hit genes per term remain in the source table; both appear above the plot.
  labels[,VolcanoLabel:=FALSE]
  for (cl in c('C1','C2','C3')) for (side in c('Lower','Higher')) {
    candidates <- unique(labels[Subtype==cl & Side==side][order(-AbsT,Gene)],by='Gene')
    selected <- head(candidates$Gene,3L)
    labels[Subtype==cl & Side==side & Gene %in% selected,VolcanoLabel:=TRUE]
  }
  gene_text <- labels[,.(RepresentativeGenes=paste(Gene,collapse=', ')),by=.(Subtype,pathway)]
  a <- merge(a,gene_text,by=c('Subtype','pathway'),sort=FALSE);setorder(a,Subtype,AOrder)
  audit <- deg[,.(InputGenes=.N,PlottedGenes=.N,ExcludedGenes=0L,
    ZeroFDR=sum(FDR==0),Higher=sum(Direction=='Higher'),Lower=sum(Direction=='Lower'),
    NotSelected=sum(Direction=='Not selected')),by=Subtype]
  overlap <- rbindlist(lapply(c('C1','C2','C3'),function(cl) {
    sel<-a[Subtype==cl];pairs<-combn(seq_len(nrow(sel)),2)
    rbindlist(lapply(seq_len(ncol(pairs)),function(j) {
      x<-sel[pairs[1,j]];y<-sel[pairs[2,j]]
      gx<-strsplit(x$HitGenes,'/',fixed=TRUE)[[1]]
      gy<-strsplit(y$HitGenes,'/',fixed=TRUE)[[1]]
      data.table(Subtype=cl,Term1=x$pathway,Term2=y$pathway,
        ORAHitJaccard=length(intersect(gx,gy))/length(union(gx,gy)))
    }))
  }))
  if(write_tables) {
    out <- file.path(STEP8_ROOT,FIG2_DIR,'tables');dir.create(out,recursive=TRUE,showWarnings=FALSE)
    fwrite(a,file.path(out,'A_selected_pathways.csv'))
    fwrite(labels,file.path(out,'A_representative_ORA_hit_genes.csv'))
    fwrite(deg,file.path(out,'A_volcano_source.csv'))
    fwrite(b,file.path(out,'B_pathway_bubble_source.csv'))
    fwrite(gsea,file.path(out,'all_Hallmark_results.csv'))
    fwrite(audit,file.path(out,'display_audit.csv'))
    fwrite(overlap,file.path(out,'A_GOBP_ORA_hit_overlap.csv'))
  }
  list(deg=deg,a=a,b=b,labels=labels,audit=audit)
}

fig2_volcano_programs <- function(cl,d) {
  color <- load_subtype_palette()$colors[[cl]]
  lower_color <- '#526779'
  a <- copy(d$a[Subtype==cl]);v <- copy(d$deg[Subtype==cl])
  labels <- unique(d$labels[Subtype==cl & VolcanoLabel],by='Gene')
  aliases <- c(EPIDERMIS_DEVELOPMENT='Epidermis development',
    PLACENTA_DEVELOPMENT='Placenta development',EPITHELIAL_CELL_DIFFERENTIATION='Epithelial differentiation',
    INTERMEDIATE_FILAMENT_ORGANIZATION='Intermediate filament organization',FORMATION_OF_PRIMARY_GERM_LAYER='Primary germ layer formation',
    AMINO_ACID_CATABOLIC_PROCESS='Amino acid catabolism',FATTY_ACID_METABOLIC_PROCESS='Fatty acid metabolism',
    XENOBIOTIC_METABOLIC_PROCESS='Xenobiotic metabolism',STEROL_METABOLIC_PROCESS='Sterol metabolism',
    COMPLEMENT_ACTIVATION='Complement activation',HUMORAL_IMMUNE_RESPONSE='Humoral immune response',
    ACUTE_PHASE_RESPONSE='Acute-phase response',MONOCARBOXYLIC_ACID_TRANSPORT='Monocarboxylic acid transport',
    RETINOL_METABOLIC_PROCESS='Retinol metabolism',SODIUM_ION_TRANSMEMBRANE_TRANSPORT='Sodium ion transport',
    GLUCOCORTICOID_BIOSYNTHETIC_PROCESS='Glucocorticoid synthesis',REGULATION_OF_HORMONE_LEVELS='Hormone-level regulation',
    FOREBRAIN_GENERATION_OF_NEURONS='Forebrain neurogenesis',SENSORY_PERCEPTION='Sensory perception',
    STEROID_METABOLIC_PROCESS='Steroid metabolism',MONOCARBOXYLIC_ACID_METABOLIC_PROCESS='Monocarboxylic acid metabolism',
    TAXIS='Taxis',REGULATION_OF_WOUND_HEALING='Wound-healing regulation',
    PROTEIN_O_LINKED_GLYCOSYLATION_VIA_N_ACETYL_GALACTOSAMINE='O-linked GalNAc glycosylation')
  a[, DisplayPathway:=unname(aliases[ShortName])]
  stopifnot(!anyNA(a$DisplayPathway))
  a[, DisplayGene:=RepresentativeGenes]
  qlim <- ceiling(max(-log10(pmax(d$a$padj,.Machine$double.xmin))))
  # Fixed physical rows: lower enrichment left, higher enrichment right.
  linebreaks <- c('Epidermis development'='Epidermis\ndevelopment',
    'Placenta development'='Placenta\ndevelopment',
    'Epithelial differentiation'='Epithelial\ndifferentiation',
    'Intermediate filament organization'='Intermediate\nfilament organization',
    'Primary germ layer formation'='Primary germ layer\nformation',
    'Amino acid catabolism'='Amino acid\ncatabolism',
    'Fatty acid metabolism'='Fatty acid\nmetabolism',
    'Xenobiotic metabolism'='Xenobiotic\nmetabolism',
    'Complement activation'='Complement\nactivation',
    'Humoral immune response'='Humoral immune\nresponse',
    'Acute-phase response'='Acute-phase\nresponse',
    'Monocarboxylic acid transport'='Monocarboxylic acid\ntransport',
    'Sodium ion transport'='Sodium ion\ntransport',
    'Glucocorticoid synthesis'='Glucocorticoid\nsynthesis',
    'Hormone-level regulation'='Hormone-level\nregulation',
    'Forebrain neurogenesis'='Forebrain\nneurogenesis',
    'Monocarboxylic acid metabolism'='Monocarboxylic acid\nmetabolism',
    'Wound-healing regulation'='Wound-healing\nregulation',
    'O-linked GalNAc glycosylation'='O-linked GalNAc\nglycosylation')
  term_grobs <- list()
  for (side in c('Lower','Higher')) {
    z<-copy(a[Side==side]);setorder(z,AOrder)
    accent<-if(side=='Higher')color else lower_color
    left<-if(side=='Lower')2 else 29
    term_grobs<-c(term_grobs,list(grid::textGrob(if(side=='Higher')'Up DEGs' else 'Down DEGs',
      x=grid::unit(left,'mm'),y=grid::unit(98,'mm'),just='left',
      gp=grid::gpar(fontfamily=STEP8_FONT$family,fontsize=7,fontface='bold',col=accent))))
    for(i in seq_len(nrow(z))) {
      y<-96-(i-1)*10.6
      label<-if(z$DisplayPathway[i] %in% names(linebreaks))linebreaks[[z$DisplayPathway[i]]] else z$DisplayPathway[i]
      term_grobs<-c(term_grobs,list(
        grid::textGrob(label,x=grid::unit(left,'mm'),y=grid::unit(y,'mm'),just=c('left','top'),
          gp=grid::gpar(fontfamily=STEP8_FONT$family,fontsize=6.3,lineheight=1,fontface='bold',col='grey15')),
        grid::textGrob(z$DisplayGene[i],x=grid::unit(left,'mm'),y=grid::unit(y-6,'mm'),just='left',
          gp=grid::gpar(fontfamily=STEP8_FONT$family,fontsize=6,fontface='italic',col='grey35')),
        grid::textGrob(sprintf('FDR %.1e',z$padj[i]),x=grid::unit(left,'mm'),y=grid::unit(y-8.1,'mm'),just='left',
          gp=grid::gpar(fontfamily=STEP8_FONT$family,fontsize=6,col='grey35')),
        grid::segmentsGrob(x0=grid::unit(left,'mm'),x1=grid::unit(left+24,'mm'),
          y0=grid::unit(y-9.4,'mm'),y1=grid::unit(y-9.4,'mm'),gp=grid::gpar(col='grey93',lwd=1.5)),
        grid::segmentsGrob(x0=grid::unit(left,'mm'),x1=grid::unit(left+24*(-log10(z$padj[i]))/qlim,'mm'),
          y0=grid::unit(y-9.4,'mm'),y1=grid::unit(y-9.4,'mm'),gp=grid::gpar(col=accent,lwd=1.5))
      ))
    }
  }
  xlim<-ceiling(max(abs(d$deg$logFC_subtype_minus_others)));ylim<-ceiling(max(d$deg$NegLog10FDR)/10)*10
  bottom<-ggplot(v,aes(logFC_subtype_minus_others,NegLog10FDR))+
    geom_hline(yintercept=-log10(.05),linetype=2,linewidth=.2,color='grey60')+
    geom_vline(xintercept=c(-1,1),linetype=2,linewidth=.2,color='grey60')+
    geom_point(data=v[Direction=='Not selected'],color='grey80',size=.25,alpha=.45)+
    geom_point(data=v[Direction!='Not selected'],aes(color=Direction),size=.4,alpha=.7)+
    geom_point(data=labels,aes(fill=Direction),shape=21,color='grey15',size=1.0,stroke=.2)+
    ggrepel::geom_text_repel(data=labels,aes(label=Gene),seed=20261026,size=2.3,
      family=STEP8_FONT$family,fontface='italic',box.padding=.25,point.padding=.15,
      bg.colour='white',bg.r=.08,min.segment.length=0,max.overlaps=Inf,
      segment.color='grey50',segment.size=.18,max.time=2,max.iter=20000)+
    scale_color_manual(values=c(Higher=color,Lower=lower_color),guide='none')+
    scale_fill_manual(values=c(Higher=color,Lower=lower_color),guide='none')+
    scale_x_continuous(limits=c(-xlim,xlim),breaks=c(-4,0,4),expand=expansion(mult=.08))+
    scale_y_continuous(limits=c(0,ylim+14),breaks=c(0,50,100),expand=expansion(mult=c(0,0)))+
    annotate('text',x=-xlim,y=ylim+12,hjust=0,vjust=1,label=sum(v$Direction=='Lower'),
      size=2.3,color=lower_color,family=STEP8_FONT$family)+
    annotate('text',x=xlim,y=ylim+12,hjust=1,vjust=1,label=sum(v$Direction=='Higher'),
      size=2.3,color=color,family=STEP8_FONT$family)+
    labs(x=expression(log[2]*'FC'),y=expression(-log[10]*'FDR'))+
    theme_step8()+theme(axis.title=element_text(size=7),axis.text=element_text(size=6.5),
      axis.title.x=element_text(margin=margin(t=1)),axis.title.y=element_text(margin=margin(r=1)),
      axis.line=element_line(linewidth=.25),axis.ticks=element_line(linewidth=.2),
      legend.position='none',plot.margin=margin(.5,.6,.3,.3,'mm'))
  rgb<-as.numeric(grDevices::col2rgb(color))/255
  lin<-ifelse(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)^2.4)
  title_color<-if(sum(lin*c(.2126,.7152,.0722))>.179)'black' else 'white'
  band<-grid::grobTree(grid::rectGrob(gp=grid::gpar(fill=color,col=NA)),
    grid::textGrob(cl,gp=grid::gpar(fontfamily=STEP8_FONT$family,fontsize=9,fontface='bold',col=title_color)),
    vp=grid::viewport(x=grid::unit(1,'mm'),y=grid::unit(101,'mm'),width=grid::unit(54,'mm'),
      height=grid::unit(4,'mm'),just=c('left','bottom')))
  volcano<-ggplotGrob(bottom)
  volcano$vp<-grid::viewport(x=grid::unit(1,'mm'),y=grid::unit(2,'mm'),width=grid::unit(54,'mm'),
    height=grid::unit(41.5,'mm'),just=c('left','bottom'))
  grid::grobTree(children=do.call(grid::gList,c(list(band),term_grobs,list(volcano))))
}

render_fig2_first_row <- function() {
  d<-prepare_fig2_data()
  plots<-lapply(seq_along(c('C1','C2','C3')),function(i) {
    p<-fig2_volcano_programs(c('C1','C2','C3')[i],d)
    p$vp<-grid::viewport(x=grid::unit((i-1)*59,'mm'),y=grid::unit(0,'mm'),
      width=grid::unit(56,'mm'),height=grid::unit(106,'mm'),just=c('left','bottom'))
    p
  })
  row<-grid::grobTree(children=do.call(grid::gList,plots))
  save_panel(row,file.path(FIG2_DIR,'Figure2A_first_row.pdf'),174,106)
}

fig2_pathway_bubbles <- function(d) {
  b <- copy(d$b)
  b[, DisplayPathway:=factor(DisplayPathway,levels=rev(unique(DisplayPathway)))]
  b[, Group:=factor(Group,levels=unique(Group))]
  b[, Subtype:=factor(Subtype,levels=c('C1','C2','C3'))]
  b[, Significance:=factor(ifelse(Significant,'FDR < 0.05','FDR >= 0.05'),
    levels=c('FDR < 0.05','FDR >= 0.05'))]
  lim <- ceiling(max(abs(b$NES))*2)/2
  ggplot(b,aes(Subtype,DisplayPathway))+
    geom_point(aes(size=LeadingEdgeN,fill=NES,color=Significance),shape=21,stroke=.65)+
    scale_fill_gradient2(low=SCALE_EFFECT_LOW,mid=SCALE_EFFECT_MID,high=SCALE_EFFECT_HIGH,
      midpoint=0,limits=c(-lim,lim),name='NES')+
    scale_size_area(max_size=5.2,name='Leading-edge\ngenes')+
    scale_color_manual(values=c('FDR < 0.05'='black','FDR >= 0.05'='grey80'),name='Outline')+
    scale_x_discrete(position='top',expand=expansion(add=.6))+
    facet_grid(Group~.,scales='free_y',space='free_y',
      labeller=labeller(Group=function(x)vapply(x,function(s)paste(strwrap(s,width=19),collapse='\n'),character(1))))+
    labs(title='Selected Hallmark pathway overview',
      subtitle='Identical pathways across C1, C2 and C3 | Study-adjusted ranked GSEA',x=NULL,y=NULL,
      caption='22 selected pathways; all three contrasts retained, including non-significant results.\nPositive NES: enriched toward the indicated subtype.\nNegative NES: enriched toward the other principal states.')+
    guides(fill=guide_colorbar(order=1),size=guide_legend(order=2),
      color=guide_legend(order=3,override.aes=list(fill='white',size=3)))+
    theme_step8()+theme(axis.line=element_blank(),axis.ticks=element_blank(),
      panel.grid.major=element_line(color='grey93',linewidth=.25),
      strip.background=element_rect(fill='grey97',color=NA),
      strip.text.y=element_text(angle=0,size=7),axis.text.x=element_text(face='bold'),
      legend.position='bottom',legend.box='vertical',legend.margin=margin(0,0,0,0),
      plot.margin=margin(4,4,3,3,'mm'))
}

build_subfigure <- function(figure,id) {
  if(figure=='Figure 2' && id %in% c('A_C1','A_C2','A_C3','B')) {
    d <- prepare_fig2_data()
    if(id=='B') return(fig2_pathway_bubbles(d))
    return(fig2_volcano_programs(sub('^A_','',id),d))
  }
  build_subfigure_before_fig2(figure,id)
}
