# Compact A and Reactome/GO BP C: display-only changes of frozen source statistics.
if(!exists('build_subfigure')) {source('8-Visualization/scripts/fig_style.R');source('8-Visualization/scripts/stat_annotate.R');source('8-Visualization/scripts/render_bulk_figures.R')}
source('8-Visualization/scripts/render_fig2_transcriptomic.R')
suppressPackageStartupMessages({library(ComplexHeatmap);library(circlize);library(grid)})
V8_ROOT <- '8-Visualization/refinements/fig2_transcriptomic_20261026_v8'
prepare_A_v8<-function(){
  f<-prepare_fig2_data
  body(f)<-parse(text=gsub('8-Visualization/scripts/fig2_gobp_display_plan.csv',paste0(V8_ROOT,'/tables/A_display_plan.csv'),paste(deparse(body(f)),collapse='\n'),fixed=TRUE))[[1]]
  d<-f(FALSE)
  stopifnot(nrow(d$a)==18L,all(d$a$Count>=5),all(d$a$padj<.05))
  fwrite(d$a,file.path(V8_ROOT,'tables/A_selected_ORA.csv'));fwrite(d$labels,file.path(V8_ROOT,'tables/A_ORA_hit_genes.csv'));fwrite(d$audit,file.path(V8_ROOT,'tables/A_DEG_counts.csv'))
  ov<-rbindlist(lapply(c('C1','C2','C3'),function(cl)rbindlist(lapply(c('Lower','Higher'),function(side){
    z<-d$a[Subtype==cl & Side==side];p<-combn(1:3,2)
    rbindlist(lapply(1:ncol(p),function(i){a<-strsplit(z$HitGenes[p[1,i]],'/',fixed=TRUE)[[1]];b<-strsplit(z$HitGenes[p[2,i]],'/',fixed=TRUE)[[1]]
      data.table(Subtype=cl,Side=side,Term1=z$pathway[p[1,i]],Term2=z$pathway[p[2,i]],HitJaccard=length(intersect(a,b))/length(union(a,b)))}))
  }))))
  stopifnot(all(ov$HitJaccard<.5));fwrite(ov,file.path(V8_ROOT,'tables/A_hit_overlap.csv'));d
}
fig2_compact_A_v8 <- function(cl,d) {
  color <- load_subtype_palette()$colors[[cl]]
  lower_color <- '#526779'
  a <- copy(d$a[Subtype==cl]);v <- copy(d$deg[Subtype==cl])
  labels <- unique(d$labels[Subtype==cl & VolcanoLabel],by='Gene')
  # Separate the three close C1 negative labels in the empty left margin.
  manual<-labels[0];manual[,`:=`(LabelX=numeric(),LabelY=numeric())]
  if(cl=='C1') {manual<-copy(labels[Direction=='Lower']);manual[,LabelX:=-4.0];manual[,LabelY:=unname(c(DDR1=48,SPINT1=31,ESRP1=14)[Gene])];stopifnot(!anyNA(manual$LabelY))}
  repel_labels<-labels[!Gene %in% manual$Gene]
  a[, DisplayPathway:=DisplayLabel]
  a[, DisplayGene:=RepresentativeGenes]
  qlim <- ceiling(max(-log10(pmax(d$a$padj,.Machine$double.xmin))))
  # Fixed physical rows: lower enrichment left, higher enrichment right.
  term_grobs <- list()
  for (side in c('Lower','Higher')) {
    z<-copy(a[Side==side]);setorder(z,AOrder)
    accent<-if(side=='Higher')color else lower_color
    left<-if(side=='Lower')1 else 29
    term_grobs<-c(term_grobs,list(grid::textGrob(if(side=='Higher')'Up DEGs' else 'Down DEGs',
      x=grid::unit(left,'mm'),y=grid::unit(70,'mm'),just='left',
      gp=grid::gpar(fontfamily=STEP8_FONT$family,fontsize=7,fontface='bold',col=accent))))
    for(i in seq_len(nrow(z))) {
      y<-67.5-(i-1)*9
      label<-z$DisplayPathway[i]
      term_grobs<-c(term_grobs,list(
        grid::textGrob(label,x=grid::unit(left,'mm'),y=grid::unit(y,'mm'),just=c('left','top'),
          gp=grid::gpar(fontfamily=STEP8_FONT$family,fontsize=6.5,lineheight=1,fontface='plain',col='grey15')),
        grid::textGrob(z$DisplayGene[i],x=grid::unit(left,'mm'),y=grid::unit(y-3.5,'mm'),just='left',
          gp=grid::gpar(fontfamily=STEP8_FONT$family,fontsize=6.5,fontface='italic',col='grey35')),
        grid::textGrob(sprintf('FDR %.1e',z$padj[i]),x=grid::unit(left,'mm'),y=grid::unit(y-5.6,'mm'),just='left',
          gp=grid::gpar(fontfamily=STEP8_FONT$family,fontsize=6.5,col='grey35')),
        grid::segmentsGrob(x0=grid::unit(left,'mm'),x1=grid::unit(left+24,'mm'),
          y0=grid::unit(y-6.9,'mm'),y1=grid::unit(y-6.9,'mm'),gp=grid::gpar(col='grey93',lwd=1.5)),
        grid::segmentsGrob(x0=grid::unit(left,'mm'),x1=grid::unit(left+24*(-log10(z$padj[i]))/qlim,'mm'),
          y0=grid::unit(y-6.9,'mm'),y1=grid::unit(y-6.9,'mm'),gp=grid::gpar(col=accent,lwd=1.5))
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
    geom_segment(data=manual,aes(x=LabelX+1.3,y=LabelY,xend=logFC_subtype_minus_others,yend=NegLog10FDR),inherit.aes=FALSE,color='grey50',linewidth=.18)+
    geom_text(data=manual,aes(x=LabelX,y=LabelY,label=Gene),inherit.aes=FALSE,family=STEP8_FONT$family,fontface='italic',size=2.3)+
    ggrepel::geom_text_repel(data=repel_labels,aes(label=Gene),seed=20261026,size=2.3,
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
  title_color<-'black'
  band<-grid::grobTree(grid::rectGrob(gp=grid::gpar(fill=color,col=NA)),
    grid::textGrob(cl,gp=grid::gpar(fontfamily=STEP8_FONT$family,fontsize=9,fontface='bold',col=title_color)),
    vp=grid::viewport(x=grid::unit(1,'mm'),y=grid::unit(73,'mm'),width=grid::unit(54,'mm'),
      height=grid::unit(4,'mm'),just=c('left','bottom')))
  volcano<-ggplotGrob(bottom)
  volcano$vp<-grid::viewport(x=grid::unit(1,'mm'),y=grid::unit(2,'mm'),width=grid::unit(54,'mm'),
    height=grid::unit(36.5,'mm'),just=c('left','bottom'))
  grid::grobTree(children=do.call(grid::gList,c(list(band),term_grobs,list(volcano))))
}

render_A_C_v8<-function(){
  dir.create(file.path(V8_ROOT,'figures'),recursive=TRUE,showWarnings=FALSE);dir.create(file.path(V8_ROOT,'tables'),recursive=TRUE,showWarnings=FALSE)
  common<-c('8-Visualization/scripts/fig_style.R','8-Visualization/scripts/stat_annotate.R','8-Visualization/scripts/render_bulk_figures.R','8-Visualization/scripts/figure_config.json','8-Visualization/scripts/render_fig2_A_C_v8.R','governance/data/subtype_master_table.csv','governance/data/subtype_label_color_dictionary.csv','0-Data/qc/study_color_dictionary.csv')
  records<-list();emit<-function(p,id,w,h,sources,notes){
    path<-file.path(V8_ROOT,'figures',paste0(id,'.pdf'));save_panel(p,path,w,h)
    src<-unique(c(common,sources));records[[length(records)+1]]<<-data.table(Panel=id,TargetOutput=path,WidthMM=w,HeightMM=h,SourceTable=paste(src,collapse=';'),SourceMD5=paste(paste0(src,'=',unname(tools::md5sum(src))),collapse=';'),Notes=notes)
  }
  a<-prepare_A_v8()
  for(cl in c('C1','C2','C3'))emit(fig2_compact_A_v8(cl,a),paste0('Figure2A_',cl,'_volcano_programs'),56,78,
    c('8-Visualization/scripts/render_fig2_transcriptomic.R',paste0(T1,c('step1_main_only_study_adjusted_deg.csv','step1_thresholded_ora_main_only_study_adjusted.csv','step1_bulk_principal_state_markers_for_step8.csv')),file.path(V8_ROOT,'tables',c('A_display_plan.csv','A_selected_ORA.csv','A_ORA_hit_genes.csv','A_DEG_counts.csv','A_hit_overlap.csv'))),
    'Three curated GO BP ORA terms/direction; FDR<.05,>=5 hits, within-direction hit Jaccard<.5. All10531 volcano genes; unchanged FDR<.05 & |log2FC|>=1 counts. 497 main samples. No tissue-origin inference.')
  plan<-data.table(
    pathway=c(paste0('REACTOME_',c('BIOLOGICAL_OXIDATIONS','FATTY_ACID_METABOLISM','METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES','COMPLEMENT_CASCADE','FORMATION_OF_FIBRIN_CLOT_CLOTTING_CASCADE','EUKARYOTIC_TRANSLATION_ELONGATION','RESPIRATORY_ELECTRON_TRANSPORT','CILIUM_ASSEMBLY','CELL_CYCLE_MITOTIC','RHO_GTPASE_EFFECTORS','INTERFERON_SIGNALING','NEUTROPHIL_DEGRANULATION','COLLAGEN_FORMATION')),
      paste0('GOBP_',c('POSITIVE_REGULATION_OF_WNT_SIGNALING_PATHWAY','EPIDERMAL_CELL_DIFFERENTIATION','GLUCOSE_CATABOLIC_PROCESS'))),
    Label=c('Biological oxidations','Fatty acid metabolism','Amino acid metabolism','Complement cascade','Fibrin clot formation','Translation elongation','Respiratory electron transport','Cilium assembly','Mitotic cell cycle','Rho GTPase effectors','Interferon signaling','Neutrophil degranulation','Collagen formation','Positive regulation of WNT signaling','Epidermal cell differentiation','Glucose catabolism'),
    Collection=c(rep('Reactome',13),rep('GO_BP',3)),
    Theme=c('Xenobiotic metabolism','Lipid metabolism','Amino acid metabolism','Complement','Coagulation','Translation','Mitochondrial respiration','Ciliary organization','Proliferation','Cytoskeletal signaling','Immune signaling','Myeloid-associated expression','Extracellular matrix','WNT regulation','Epithelial differentiation','Glucose metabolism'))
  plan[,DisplayOrder:=.I];plan[,SelectionRule:='Curated distinct functional themes from the discussion; require original FDR<.05 in at least one contrast; all three contrasts retained. Not independent or subtype-specific evidence.']
  g<-read_source(paste0(T1,'step1_ranked_gsea_main_only_study_adjusted.csv'))
  x<-merge(plan,g,by=c('pathway','Collection'),all.x=TRUE,sort=FALSE);setorder(x,DisplayOrder,Subtype)
  stopifnot(nrow(x)==48,!anyNA(x$NES),all(x[,.(ok=any(padj<.05)),by=pathway]$ok),all(x[,.N,by=pathway]$N==3))
  x[,LeadingEdgeN:=lengths(strsplit(leadingEdge,';',fixed=TRUE))]
  fwrite(plan,file.path(V8_ROOT,'tables/C_display_plan.csv'));fwrite(x,file.path(V8_ROOT,'tables/C_GSEA_selected.csv'))
  # Lead-edge overlap audit is descriptive and cannot prove full-reference-set independence.
  pairs<-combn(plan$pathway,2)
  ov<-rbindlist(lapply(c('C1','C2','C3'),function(cl)rbindlist(lapply(seq_len(ncol(pairs)),function(i){
    a<-strsplit(x[Subtype==cl & pathway==pairs[1,i],leadingEdge],';',fixed=TRUE)[[1]];b<-strsplit(x[Subtype==cl & pathway==pairs[2,i],leadingEdge],';',fixed=TRUE)[[1]]
    data.table(Subtype=cl,Term1=pairs[1,i],Term2=pairs[2,i],LeadingEdgeJaccard=length(intersect(a,b))/length(union(a,b)))
  }))))
  fwrite(ov,file.path(V8_ROOT,'tables/C_leading_edge_overlap.csv'))
  m<-q<-matrix(NA_real_,nrow(plan),3,dimnames=list(plan$pathway,c('C1','C2','C3')))
  for(i in seq_len(nrow(x))){m[x$pathway[i],x$Subtype[i]]<-x$NES[i];q[x$pathway[i],x$Subtype[i]]<-x$padj[i]}
  lim<-ceiling(max(abs(m))*2)/2;pal<-load_subtype_palette()$colors[c('C1','C2','C3')]
  labels<-vapply(plan$Label,function(z)paste(strwrap(z,32),collapse='\n'),character(1))
  hm<-Heatmap(m,name='NES',col=colorRamp2(c(-lim,0,lim),c(SCALE_EFFECT_LOW,'white',SCALE_EFFECT_HIGH)),
    cluster_rows=FALSE,cluster_columns=FALSE,cluster_row_slices=FALSE,row_split=factor(plan$Collection,levels=c('Reactome','GO_BP'),labels=c('Reactome','GO BP')),
    row_labels=labels,row_names_side='left',row_names_max_width=unit(62,'mm'),row_names_gp=gpar(fontsize=7,fontfamily='Arial'),
    row_title_gp=gpar(fontsize=7,fontfamily='Arial'),row_title_rot=90,row_gap=unit(2,'mm'),column_names_side='top',column_names_rot=0,column_names_gp=gpar(fontsize=8,fontface='bold',fontfamily='Arial'),
    width=unit(25,'mm'),rect_gp=gpar(col='white',lwd=.5),
    top_annotation=HeatmapAnnotation(Subtype=colnames(m),col=list(Subtype=pal),show_annotation_name=FALSE,show_legend=FALSE,simple_anno_size=unit(1.5,'mm')),
    cell_fun=function(j,i,x,y,w,h,fill){if(q[i,j]<.05)grid.points(x,y,pch=16,size=unit(.8,'mm'),gp=gpar(col=if(abs(m[i,j])>3)'white' else 'black'))},
    heatmap_legend_param=list(at=c(-lim,0,lim),direction='horizontal',legend_width=unit(30,'mm'),title_gp=gpar(fontsize=7),labels_gp=gpar(fontsize=7)))
  p<-function()draw(hm,heatmap_legend_side='bottom',padding=unit(c(3,3,3,3),'mm'),heatmap_legend_list=list(Legend(labels='FDR < 0.05',type='points',pch=16,size=unit(1,'mm'),labels_gp=gpar(fontsize=7))))
  emit(p,'Figure2C_Reactome_GO_GSEA',112,112,c(paste0(T1,'step1_ranked_gsea_main_only_study_adjusted.csv'),file.path(V8_ROOT,'tables',c('C_display_plan.csv','C_GSEA_selected.csv','C_leading_edge_overlap.csv'))),
    '13 Reactome +3 GO BP representative pathways; all3 original NES/FDR results retained including negatives/nonsignificant. Dot FDR<.05. Collection blocks distinct. Curated explanatory selection, not independent validation.')
  fwrite(rbindlist(records,fill=TRUE),file.path(V8_ROOT,'A_C_manifest.csv'));invisible(rbindlist(records,fill=TRUE))
}
if(sys.nframe()==0)render_A_C_v8()
