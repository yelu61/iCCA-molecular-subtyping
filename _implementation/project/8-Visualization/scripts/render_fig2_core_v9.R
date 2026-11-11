# Display-only refinement of frozen DEG, GSEA and WGCNA results.
if(!exists('STEP8_ROOT'))source('8-Visualization/scripts/fig_style.R')
suppressPackageStartupMessages({library(data.table);library(ComplexHeatmap);library(circlize);library(grid);library(ggplot2)})
V9<-'8-Visualization/refinements/fig2_transcriptomic_20261026_v9'
render_core_v9<-function(){
  dir.create(file.path(V9,'figures'),recursive=TRUE,showWarnings=FALSE);dir.create(file.path(V9,'tables'),showWarnings=FALSE)
  records<-list();common<-c('8-Visualization/scripts/render_fig2_core_v9.R','8-Visualization/scripts/fig_style.R','8-Visualization/scripts/figure_config.json','governance/data/subtype_master_table.csv','governance/data/subtype_label_color_dictionary.csv','0-Data/qc/study_color_dictionary.csv')
  emit<-function(p,id,w,h,src,notes){
    target<-file.path(V9,'figures',paste0(id,'.pdf'));save_panel(p,target,w,h);src<-unique(c(src,common))
    records[[length(records)+1L]]<<-data.table(Panel=id,TargetOutput=target,WidthMM=w,HeightMM=h,SourceTable=paste(src,collapse=';'),SourceMD5=paste(paste0(src,'=',unname(tools::md5sum(src))),collapse=';'),Notes=notes,Status='Rendered',OutputMD5=unname(tools::md5sum(target)))
  }
  pal<-load_subtype_palette()$colors[c('C1','C2','C3')]
  # Down is left/slate; Up is right/subtype accent, matching the accepted volcanoes.
  degsrc<-'1-Clustering/publication/tables/step1_main_only_study_adjusted_deg.csv';d<-fread(degsrc)
  stopifnot(all(c('Subtype','FDR','logFC_subtype_minus_others')%in%names(d)))
  cnt<-rbindlist(lapply(c(0,.5,1),function(cut)d[,.(Down=sum(FDR<.05 & logFC_subtype_minus_others<0 & abs(logFC_subtype_minus_others)>=cut),Up=sum(FDR<.05 & logFC_subtype_minus_others>0 & abs(logFC_subtype_minus_others)>=cut)),by=Subtype][,Cutoff:=cut]))
  cnt<-melt(cnt,id.vars=c('Subtype','Cutoff'),variable.name='Direction',value.name='Count');cnt[,SignedCount:=ifelse(Direction=='Down',-Count,Count)];cnt[,Fill:=ifelse(Direction=='Down','#526779',unname(pal[Subtype]))]
  cnt[,Threshold:=factor(Cutoff,levels=c(1,.5,0),labels=c('1.0','0.5','0'))];cnt[,Subtype:=factor(Subtype,levels=c('C1','C2','C3'))]
  fwrite(cnt,file.path(V9,'tables/A_DEG_threshold_counts.csv'))
  p<-ggplot(cnt,aes(SignedCount,Threshold,fill=Fill))+geom_vline(xintercept=0,color='grey50',linewidth=.3)+geom_col(width=.56)+
    geom_text(aes(label=format(Count,big.mark=','),hjust=ifelse(Direction=='Down',1.1,-.1)),family=STEP8_FONT$family,size=2.6)+
    scale_fill_identity()+scale_x_continuous(limits=c(-8200,8200),breaks=c(-6000,0,6000),labels=c('6,000','0','6,000'))+
    facet_wrap(~Subtype,nrow=1)+labs(x='DEG count (left: Down; right: Up)',y=expression('Minimum |log'[2]*'FC|'))+theme_step8(base_size=8)+
    theme(strip.text=element_text(size=9,face='bold'),strip.background=element_blank(),axis.text=element_text(size=7.5),axis.title=element_text(size=8),panel.spacing=unit(5,'mm'),axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=margin(3,4,3,3,'mm'))
  emit(p,'SFigure2A_DEG_threshold_counts',174,65,c(degsrc,file.path(V9,'tables/A_DEG_threshold_counts.csv')),'497 main samples; source subtype-minus-pooled-others Study-adjusted DEG; FDR<.05, FC cutoff0/0.5/1.0. Down left/slate and Up right/subtype color, counts labelled as positive counts. Counts unchanged.')
  gsrc<-'1-Clustering/publication/tables/step1_ranked_gsea_main_only_study_adjusted.csv';g<-fread(gsrc)
  plan<-data.table(pathway=c(paste0('GOBP_',c('AMINO_ACID_CATABOLIC_PROCESS','FATTY_ACID_METABOLIC_PROCESS','XENOBIOTIC_METABOLIC_PROCESS','STEROL_METABOLIC_PROCESS','CYTOPLASMIC_TRANSLATION','CILIUM_ORGANIZATION','POSITIVE_REGULATION_OF_WNT_SIGNALING_PATHWAY','GLANDULAR_EPITHELIAL_CELL_DIFFERENTIATION','MITOTIC_NUCLEAR_DIVISION','EPIDERMAL_CELL_DIFFERENTIATION','GLUCOSE_CATABOLIC_PROCESS')),paste0('REACTOME_',c('COMPLEMENT_CASCADE','RESPIRATORY_ELECTRON_TRANSPORT','RHO_GTPASE_EFFECTORS','COLLAGEN_FORMATION'))),
    Label=c('Amino acid catabolism','Fatty acid metabolism','Xenobiotic metabolism','Sterol metabolism','Cytoplasmic translation','Cilium organization','Positive regulation of WNT signaling','Glandular epithelial cell differentiation','Mitotic nuclear division','Epidermal cell differentiation','Glucose catabolism','Complement cascade','Respiratory electron transport','Rho GTPase effectors','Collagen formation'),Collection=c(rep('GO_BP',11),rep('Reactome',4)),Context=c(rep('C1',4),rep('C2',4),rep('C3',3),'C1','C2','C3','C3'),DisplayOrder=1:15)
  plan[,SelectionRule:='Question-driven non-exhaustive functional themes, GO BP primary and Reactome complementary; original FDR<.05 in context; all three contrasts retained. No independent validation or exclusive subtype identity.']
  x<-merge(plan,g,by=c('pathway','Collection'),all.x=TRUE);setorder(x,DisplayOrder,Subtype);stopifnot(nrow(x)==45,!anyNA(x$NES),all(x[Subtype==Context,padj]<.05))
  fwrite(plan,file.path(V9,'tables/B_display_plan.csv'));fwrite(x,file.path(V9,'tables/B_GSEA_selected.csv'));fwrite(g,file.path(V9,'tables/B_GSEA_complete_results.csv'))
  pairs<-combn(plan$pathway,2);ov<-rbindlist(lapply(c('C1','C2','C3'),function(cl)rbindlist(lapply(seq_len(ncol(pairs)),function(i){a<-strsplit(x[Subtype==cl & pathway==pairs[1,i],leadingEdge],';',fixed=TRUE)[[1]];b<-strsplit(x[Subtype==cl & pathway==pairs[2,i],leadingEdge],';',fixed=TRUE)[[1]];data.table(Subtype=cl,Term1=pairs[1,i],Term2=pairs[2,i],LeadingEdgeJaccard=length(intersect(a,b))/length(union(a,b)))}))))
  fwrite(ov,file.path(V9,'tables/B_leading_edge_overlap.csv'))
  draw_gsea<-function(z,split=FALSE){
    ids<-unique(z$pathway);label<-z[match(ids,pathway),Label];m<-q<-matrix(NA_real_,length(ids),3,dimnames=list(ids,c('C1','C2','C3')))
    for(i in seq_len(nrow(z))){m[z$pathway[i],z$Subtype[i]]<-z$NES[i];q[z$pathway[i],z$Subtype[i]]<-z$padj[i]}
    stopifnot(!anyNA(m),!anyNA(q));lim<-4.5
    hm<-Heatmap(m,name='NES',col=colorRamp2(c(-lim,0,lim),c(SCALE_EFFECT_LOW,'white',SCALE_EFFECT_HIGH)),cluster_rows=FALSE,cluster_columns=FALSE,cluster_row_slices=FALSE,
      row_split=if(split)factor(z[match(ids,pathway),Collection],levels=c('GO_BP','Reactome'),labels=c('GO BP','Reactome')) else NULL,
      row_labels=vapply(label,function(s)paste(strwrap(s,42),collapse='\n'),character(1)),row_names_side='left',row_names_max_width=unit(73,'mm'),row_names_gp=gpar(fontsize=7.5,fontfamily='Arial'),row_title_gp=gpar(fontsize=8,fontfamily='Arial'),row_gap=unit(3,'mm'),column_names_side='top',column_names_rot=0,column_names_gp=gpar(fontsize=9,fontface='bold',fontfamily='Arial'),width=unit(30,'mm'),rect_gp=gpar(col='white',lwd=.6),
      top_annotation=HeatmapAnnotation(Subtype=colnames(m),col=list(Subtype=pal),show_annotation_name=FALSE,show_legend=FALSE,simple_anno_size=unit(2,'mm')),
      cell_fun=function(j,i,x,y,w,h,fill){if(q[i,j]<.05)grid.points(x,y,pch=16,size=unit(.85,'mm'),gp=gpar(col=if(abs(m[i,j])>3)'white' else 'black'))},heatmap_legend_param=list(at=c(-lim,0,lim),direction='horizontal',legend_width=unit(32,'mm'),title_gp=gpar(fontsize=7.5),labels_gp=gpar(fontsize=7.5)))
    function()draw(hm,heatmap_legend_side='bottom',padding=unit(rep(4,4),'mm'),heatmap_legend_list=list(Legend(labels='FDR < 0.05',type='points',pch=16,size=unit(1,'mm'),labels_gp=gpar(fontsize=7.5))))
  }
  emit(draw_gsea(x,TRUE),'Figure2B_GO_Reactome_GSEA',122,117,c(gsrc,file.path(V9,'tables',c('B_display_plan.csv','B_GSEA_selected.csv','B_leading_edge_overlap.csv'))),'11 GO BP +4 Reactome functional themes; all45 frozen NES/FDR; fixed +/-4.5 NES scale; dot FDR<.05; includes negative and nonsignificant contrasts; no row rescaling.')
  # More comprehensive supplementary selection: all Hallmark; top8 FDR per direction/contrast for other collections.
  supplans<-list()
  for(co in c('GO_BP','Reactome','Hallmark','KEGG')) {
    z<-copy(g[Collection==co]);z[,Side:=ifelse(NES>0,'Positive','Negative')];setorder(z,padj,-AbsNES,pathway)
    ids<-if(co=='Hallmark')sort(unique(z$pathway)) else unique(z[padj<.05,head(.SD,8),by=.(Subtype,Side)]$pathway)
    sel<-g[Collection==co & pathway%in%ids];sel[,PeakSubtype:=c('C1','C2','C3')[which.max(NES)],by=pathway];sel[,BestFDR:=min(padj),by=pathway];setorder(sel,PeakSubtype,BestFDR,pathway,Subtype)
    ordered<-unique(sel$pathway);sel[,Label:=vapply(pathway,function(s)tools::toTitleCase(tolower(gsub('_',' ',sub('^(GOBP|REACTOME|HALLMARK|KEGG)_','',s)))),character(1))]
    spp<-unique(sel[,.(Collection,pathway,PeakSubtype,BestFDR)]);spp[,SelectionRule:=if(co=='Hallmark')'All50 Hallmark sets, including nonsignificant' else 'Union of top8 FDR-significant positive and negative terms per subtype; ranked by frozen FDR then AbsNES; no redundancy filtering in comprehensive supplement'];supplans[[co]]<-spp
    for(k in seq_len(ceiling(length(ordered)/22))){idx<-ordered[((k-1)*22+1):min(k*22,length(ordered))];part<-sel[pathway%in%idx];part[,Ord:=match(pathway,idx)];setorder(part,Ord,Subtype)
      id<-paste0('SFigure2B_',co,'_',k);src<-file.path(V9,'tables',paste0(id,'_source.csv'));fwrite(part,src)
      emit(draw_gsea(part),id,132,max(72,28+length(idx)*6.2),c(gsrc,src),paste0(co,': complete Hallmark or top8 per direction/subtype union; all3 source contrasts incl negative/nonsignificant; same NES scale as main. KEGG is KEGG MEDICUS, distinct from GSVA metabolic reference.'))
    }
  }
  fwrite(rbindlist(supplans),file.path(V9,'tables/B_supplement_display_plan.csv'))
  # WGCNA main: eight functional modules, one cross-study estimate per subtype.
  old<-'8-Visualization/refinements/fig2_transcriptomic_20261026_v8/tables/'
  sel<-fread(paste0(old,'wgcna_v8_display_selection.csv'));setorder(sel,DisplayOrder)
  ms<-'2-TranscriptomicCharacterization/publication/tables/step2_wgcna_subtype_random_effects_meta.csv';allmeta<-fread(ms);meta<-merge(sel,allmeta,by='Module');setorder(meta,DisplayOrder,TargetSubtype)
  stopifnot(nrow(meta)==24,all(is.finite(meta$RandomEffect)));fwrite(meta,file.path(V9,'tables/C_wgcna_main_effects.csv'));fwrite(allmeta,file.path(V9,'tables/C_wgcna_all_cross_study.csv'))
  ramp<-colorRamp(c(SCALE_EFFECT_LOW,'white',SCALE_EFFECT_HIGH));lim<-2.5;colfun<-function(v)rgb(ramp(pmin(1,pmax(0,(v/lim+1)/2))),maxColorValue=255)
  txt<-function(s,x,y,size=7.5,just='left',face='plain',col='grey15')textGrob(s,x=unit(x,'mm'),y=unit(y,'mm'),just=just,gp=gpar(fontfamily=STEP8_FONT$family,fontsize=size,fontface=face,col=col))
  grobs<-list(txt('Module',3,87,face='bold'),txt('Program / hub candidates',29,87,face='bold'),txt('Cross-study effect',102,88,7.5,'center','bold'),txt('Positive',125,88,7,'center','bold'),txt('studies',125,84,7,'center'),txt('Preservation',143,88,7,'center','bold'),txt('Z summary',143,84,6.7,'center'))
  for(j in 1:3)grobs<-c(grobs,list(txt(paste0('C',j),102+(j-2)*8,80,8,'center','bold',pal[[paste0('C',j)]])))
  ys<-c(72,64,54,46,38,28,20,12)
  for(i in 1:8){r<-sel[i];y<-ys[i];context<-meta[Module==r$Module & TargetSubtype==r$DisplayClass]
    grobs<-c(grobs,list(rectGrob(x=unit(3.5,'mm'),y=unit(y,'mm'),width=unit(1.4,'mm'),height=unit(6,'mm'),gp=gpar(fill=pal[[r$DisplayClass]],col=NA)),txt(sub('^ME','',r$Module),6,y,7.2),txt(r$Program,29,y+1.8,7.5,face='bold'),txt(r$HubCandidates,29,y-1.9,7,face='italic',col='grey40'),txt(paste0(context$PositiveEligibleStudies,'/',context$DirectionEligibleStudies),125,y,7.5,'center'),txt(sprintf('%.1f',r$ZsummaryPreservation),143,y,7.5,'center')))
    for(j in 1:3){z<-meta[Module==r$Module & TargetSubtype==paste0('C',j)];xx<-102+(j-2)*8;grobs<-c(grobs,list(rectGrob(x=unit(xx,'mm'),y=unit(y,'mm'),width=unit(7.6,'mm'),height=unit(6.2,'mm'),gp=gpar(fill=colfun(z$RandomEffect),col=NA))))
      if(z$FDR<.05)grobs<-c(grobs,list(pointsGrob(x=unit(xx,'mm'),y=unit(y,'mm'),pch=16,size=unit(.8,'mm'),gp=gpar(col=if(abs(z$RandomEffect)>1.7)'white' else 'grey15'))))}
  }
  lg<-Legend(col_fun=colorRamp2(c(-lim,0,lim),c(SCALE_EFFECT_LOW,'white',SCALE_EFFECT_HIGH)),title='Effect (SD)',direction='horizontal',at=c(-lim,0,lim),legend_width=unit(28,'mm'),title_position='lefttop',title_gp=gpar(fontsize=7),labels_gp=gpar(fontsize=7))
  pl<-function(){grid.draw(grobTree(children=do.call(gList,grobs)));draw(lg,x=unit(4,'mm'),y=unit(4,'mm'),just=c('left','center'));grid.points(x=unit(74,'mm'),y=unit(4,'mm'),pch=16,size=unit(.8,'mm'));grid.text('FDR < 0.05',x=unit(77,'mm'),y=unit(4,'mm'),just='left',gp=gpar(fontfamily='Arial',fontsize=7))}
  emit(pl,'Figure2C_WGCNA_programs',154,94,c(ms,paste0(old,c('wgcna_v8_display_selection.csv','wgcna_v8_function_anchors.csv','wgcna_v8_hub_candidates.csv')),file.path(V9,'tables/C_wgcna_main_effects.csv')),'Eight representative modules, frozen cross-study random-effects estimates; all3 contrasts. Positive studies column refers to color-coded context subtype, not all3 or significance. Preservation Z summary separate from association. Discovery/heldout moved toS2C; qualified network sensitivity remains.')
  # Complete cross-study module results with heterogeneity explicitly retained.
  cls<-fread('2-TranscriptomicCharacterization/publication/tables/step2_wgcna_module_subtype_classification.csv');ord<-c(sel$Module,setdiff(cls$Module,sel$Module));mat<-iq<-q<-matrix(NA_real_,length(ord),3,dimnames=list(ord,c('C1','C2','C3')))
  for(i in 1:nrow(allmeta)){r<-allmeta[i];mat[r$Module,r$TargetSubtype]<-r$RandomEffect;q[r$Module,r$TargetSubtype]<-r$FDR;iq[r$Module,r$TargetSubtype]<-r$I2Percent}
  for(kind in c('effects','heterogeneity')){
    zz<-copy(allmeta);zz[,ModuleLabel:=factor(sub('^ME','',Module),levels=rev(sub('^ME','',ord)))];zz[,Subtype:=factor(TargetSubtype,levels=c('C1','C2','C3'))]
    zz[,Value:=if(kind=='effects')RandomEffect else I2Percent]
    p<-ggplot(zz,aes(Subtype,ModuleLabel,fill=Value))+geom_tile(color='white',linewidth=.35)+
      scale_y_discrete(expand=c(0,0))+scale_x_discrete(position='top',expand=c(0,0))+labs(x=NULL,y=NULL)+theme_step8(base_size=8)+
      theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(size=9,face='bold'),axis.text.y=element_text(size=8),legend.position='bottom',plot.margin=margin(4,4,4,4,'mm'))
    if(kind=='effects')p<-p+geom_point(data=zz[FDR<.05],size=.7)+scale_fill_gradient2(low=SCALE_EFFECT_LOW,mid='white',high=SCALE_EFFECT_HIGH,limits=c(-2.5,2.5),oob=scales::squish,name='Effect (SD)')
    else p<-p+scale_fill_gradient(low='white',high='#274A61',limits=c(0,100),name='I2 (%)')
    emit(p,paste0('SFigure2C_all_modules_',kind),87,160,c(ms,'2-TranscriptomicCharacterization/publication/tables/step2_wgcna_module_subtype_classification.csv'),paste0('All24 modules x3 target contrasts; frozen random-effects ',kind,'; source eligibility, CIs and direction retained in complete source table.'))
  }
  fwrite(rbindlist(records,fill=TRUE),file.path(V9,'core_manifest.csv'));invisible(records)
}
if(sys.nframe()==0)render_core_v9()
