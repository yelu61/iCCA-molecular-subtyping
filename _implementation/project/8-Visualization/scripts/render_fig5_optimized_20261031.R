# Scientific display revision. The frozen 114-gene primary model is unchanged.
source('8-Visualization/scripts/render_fig5_completion_20261030.R')
suppressPackageStartupMessages({library(ggrepel);library(survival)})
FIG5_OPT <- '8-Visualization/refinements/fig5_optimized_20261031_v5'
opt_path <- function(...)file.path(STEP8_ROOT,FIG5_OPT,...)
completion_path <- opt_path
width_mm <- 174
opt_read <- function(n)fread(opt_path('source_data',paste0(n,'.csv')))
size_run <- '5-ClassifierValidation/runs/bulk-rnaseq__20261030_step5_marker_size_v1'
opt_names <- c(B='Figure5B_feature_evidence_RF',C='Figure5C_strategy_and_marker_count',E='Figure5E_external_expression',F='Figure5F_external_survival',S11F='Supplement_S11F_feature_strategies',S11G='Supplement_S11G_marker_count',S11I='Supplement_S11I_development_expression')
opt_dims <- list(B=c(174,110),C=c(174,105),E=c(174,240),F=c(174,112),S11F=c(174,105),S11G=c(174,140),S11I=c(174,240))

opt_evidence <- function(){
 d<-opt_read('feature_evidence_114');all<-opt_read('feature_evidence_all_genes')
 stopifnot(all(all$RFRank>0),nrow(all)==10531,nrow(d)==114)
 left<-function(){
  grid.newpage();pushViewport(viewport(xscale=c(0,84),yscale=c(0,84)))
  txt<-function(s,x,y,size=6.5,bold=FALSE,just='left',col=FIG5_INK)grid.text(s,unit(x,'native'),unit(y,'native'),just=just,gp=gpar(fontfamily='Arial',fontsize=size,fontface=if(bold)'bold' else 'plain',col=col))
  seg<-function(x,y,col='grey85')grid.lines(unit(x,'native'),unit(y,'native'),gp=gpar(col=col,lwd=.45))
  rect<-function(x,y,w,h,fill)grid.rect(unit(x,'native'),unit(y,'native'),unit(w,'native'),unit(h,'native'),gp=gpar(fill=fill,col=NA))
  x0<-11;x1<-81;xp<-function(i)x0+(i-.5)*(x1-x0)/114
  txt('Selection stability',x0,79,7.5,TRUE)
  for(v in c(0,.5,1)){yy<-54+20*v;seg(c(x0,x1),c(yy,yy),'grey92');txt(sprintf('%d%%',100*v),x0-1,yy,6,just='right')}
  for(i in seq_len(nrow(d))){seg(rep(xp(i),2),c(54,54+20*d$ConsensusFraction[i]),'#B7C6D0');grid.points(unit(xp(i),'native'),unit(54+20*d$ConsensusFraction[i],'native'),pch=16,size=unit(.7,'mm'),gp=gpar(col=FIG5_INK))}
  txt('Support in the final development fit',x0,49,7,TRUE)
  for(j in 1:3){yy<-44-(j-1)*3;txt(c('limma','RF','LASSO')[j],x0-1,yy,6,just='right');for(i in seq_len(nrow(d)))rect(xp(i),yy,(x1-x0)/114,2.6,if(d[[c('LimmaFull','RFFull','LASSOFull')[j]]][i])FIG5_INK else '#E7EDF1')}
  txt('Mean expression by reference class',x0,32,7,TRUE)
  ramp<-circlize::colorRamp2(c(-1.5,0,1.5),STEP8_CONFIG$effect_colors)
  for(j in 1:3){cl<-c('C1','C2','C3')[j];yy<-27-(j-1)*3;txt(cl,x0-1,yy,6,TRUE,'right',FIG5_SUBTYPE[cl]);for(i in seq_len(nrow(d)))rect(xp(i),yy,(x1-x0)/114,2.6,ramp(d[[paste0(cl,'MeanTrainingZ')]][i]))}
  for(k in c(1,20,40,60,80,100,114))txt(k,xp(k),17,5.8,just='centre')
  txt('114 panel genes ordered by RF importance',x0,12,6.2)
  rect(12,6,1.6,1.6,FIG5_INK);txt('Selected',14,6,5.8);rect(30,6,1.6,1.6,'#E7EDF1');txt('Not selected',32,6,5.8)
  for(i in 1:20)rect(63+(i-.5)*18/20,6,18/20,1.8,ramp(-1.5+(i-.5)*3/20))
  txt('Mean z',59,6,5.8,just='right');txt('-1.5',63,2.5,5.5,just='centre');txt('1.5',81,2.5,5.5,just='centre');popViewport()
 }
 d[,Support:=factor(Methods,levels=c(2,3),labels=c('2 methods','3 methods'))]
 labels<-rbind(d[order(RFRank)][1:6],tail(d[order(RFRank)],2));labels<-unique(labels,by='Gene')
 p<-ggplot(all,aes(RFRank,RFImportance))+geom_line(colour='#C9D1D7',linewidth=.35)+
  geom_point(data=d,aes(fill=ConsensusFraction,shape=Support),size=2.1,colour=FIG5_INK,stroke=.3)+
  geom_text_repel(data=labels,aes(label=Gene),size=2.25,seed=20261031,min.segment.length=0,segment.size=.2,box.padding=.3,max.overlaps=Inf)+
  scale_x_log10(breaks=c(1,10,100,1000,10000),labels=c('1','10','100','1,000','10,000'))+
  scale_shape_manual(values=c(21,24))+scale_fill_gradient(low='#E7EEF3',high='#3F718E',limits=c(0,1),breaks=c(0,.5,1),labels=scales::label_percent())+
  labs(title='RF ranking and consensus support',x='Rank among all screened genes (log scale)',y='Scaled permutation importance',fill='Selection frequency',shape=NULL)+fig5_theme()+
  theme(plot.title=element_text(size=7.8),axis.title=element_text(size=6.8),axis.text=element_text(size=6.2),legend.position='bottom',legend.box='vertical',legend.text=element_text(size=6),legend.title=element_text(size=6),legend.key.size=unit(2.5,'mm'))+
  guides(fill=guide_colourbar(display='rectangles',barwidth=unit(22,'mm'),barheight=unit(2,'mm')),shape=guide_legend(override.aes=list(fill='white')))
 left_grob<-grid.grabExpr(left(),width=84/25.4,height=84/25.4)
 (wrap_elements(full=left_grob)+p+plot_layout(widths=c(1,1)))+plot_annotation(title='Evidence supporting the 114-gene classifier panel',subtitle='20 genes supported by all three methods; 71/114 selected in at least 80% of outer fits',caption='Left: 25 outer training fits; equal-width panel-gene columns in RF order. Right: all screened genes; highlighted points are the final panel.\nRF values describe the full development fit. Neither importance nor selection frequency establishes an optimal panel size.',theme=fig5_theme()+theme(plot.caption=element_text(size=6.2)))
}

opt_size_data <- function(){
 d<-fread(file.path(size_run,'tables/marker_size_performance.csv'));fwrite(d,opt_path('source_data/marker_size_performance.csv'));d
}
opt_size_curve <- function(metric='BalancedAccuracy'){
 d<-opt_size_data()[Metric==metric];ref<-fig5_read('nested_performance_ci')[Model==FIG5_PRIMARY & Metric==metric]
 ggplot(d,aes(K,Estimate))+geom_hline(yintercept=ref$Estimate,linetype=2,colour='#9B6176',linewidth=.5)+geom_ribbon(aes(ymin=Lower95,ymax=Upper95),fill='#D6E4EC',alpha=.6)+geom_line(colour=FIG5_INK,linewidth=.65)+geom_point(colour=FIG5_INK,size=2)+scale_x_continuous(breaks=c(5,10,20,40,60,80))+scale_y_continuous(limits=c(.91,1.003),breaks=c(.92,.94,.96,.98,1),labels=scales::label_percent(accuracy=1))+
  annotate('text',x=79,y=ref$Estimate+.006,label=sprintf('Complete pipeline: %.2f%%',100*ref$Estimate),hjust=1,size=2.25,colour='#8B5068')+
  labs(title='Cost of reducing the feature set',x='Training-selected genes (K)',y='OOF balanced accuracy')+fig5_theme()+theme(plot.title=element_text(size=8),axis.title=element_text(size=7),axis.text=element_text(size=6.5))
}
opt_comparison <- function(){
 p<-fread(file.path(FIG5_STRATEGY_RUN,'tables/performance_ci.csv'))[Metric=='BalancedAccuracy' & Budget=='Native']
 n<-fread(file.path(FIG5_STRATEGY_RUN,'tables/feature_count_summary.csv'));p<-merge(p,n,by=c('Spec','Budget','Strategy'))
 labs<-c(Limma='limma',RF='Random forest',LASSO='LASSO',Consensus2='Consensus (primary)',Core3='All three')
 p[,Strategy:=factor(Strategy,levels=rev(names(labs)))];p[,Primary:=Strategy=='Consensus2']
 fwrite(p,opt_path('source_data/5C_native_strategy.csv'))
 q<-ggplot(p,aes(Estimate,Strategy,colour=Primary))+geom_errorbar(aes(xmin=Lower95,xmax=Upper95),orientation='y',width=.18,linewidth=.5)+geom_point(size=2)+geom_text(aes(x=1.014,label=sprintf('%g',MedianFeatures)),colour=FIG5_INK,size=2.3)+scale_colour_manual(values=c('FALSE'='#8DA0AC','TRUE'=FIG5_INK),guide='none')+scale_y_discrete(labels=labs)+scale_x_continuous(limits=c(.95,1.026),breaks=c(.95,.97,.99),labels=scales::label_percent(accuracy=1))+
  annotate('text',x=1.014,y=5.4,label='Genes',size=2.3)+labs(title='Feature-selection strategies',x='OOF balanced accuracy',y=NULL)+fig5_theme()+theme(axis.title=element_text(size=7),axis.text=element_text(size=6.5),axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.title=element_text(size=8))
 (q+opt_size_curve()+plot_layout(widths=c(1,1.3)))+plot_annotation(title='Feature strategy and the cost of a smaller panel',subtitle='The complete consensus classifier remains primary; no optimal K is selected',caption='Identical repeated nested CV; points and 95% conditional patient-bootstrap intervals. Left: median gene count across outer fits.\nRight: training-only consensus genes ranked by method support, then limma F rank. Complete pipeline is a reference, not fixed K = 114.',theme=fig5_theme()+theme(plot.caption=element_text(size=6.2)))
}

opt_marker_supplement <- function(){
 d<-opt_size_data();d[,MetricLabel:=factor(Metric,levels=c('BalancedAccuracy','MacroF1','MacroAUROC'),labels=c('Balanced accuracy','Macro-F1','Macro-AUROC'))]
 ref<-fig5_read('nested_performance_ci')[Model==FIG5_PRIMARY & Metric%in%d$Metric];ref<-rbind(ref,data.table(Metric='MacroAUROC',Estimate=mean(fig5_read('nested_performance_ci')[Model==FIG5_PRIMARY & Metric%in%paste0(c('C1','C2','C3'),'_AUROC'),Estimate])),fill=TRUE);ref[,MetricLabel:=factor(Metric,levels=c('BalancedAccuracy','MacroF1','MacroAUROC'),labels=c('Balanced accuracy','Macro-F1','Macro-AUROC'))]
 p<-ggplot(d,aes(K,Estimate))+geom_hline(data=ref,aes(yintercept=Estimate),linetype=2,colour='#9B6176',linewidth=.45)+geom_ribbon(aes(ymin=Lower95,ymax=Upper95),fill='#D6E4EC',alpha=.65)+geom_line(colour=FIG5_INK,linewidth=.55)+geom_point(colour=FIG5_INK,size=1.7)+facet_wrap(~MetricLabel,nrow=1,scales='free_y')+scale_x_continuous(breaks=c(5,20,40,60,80))+labs(x='Training-selected consensus genes (K)',y='OOF performance')+fig5_theme()+theme(strip.text=element_text(size=7),axis.text=element_text(size=6.5))
 # Paired differences use the same patients for each K and complete pipeline.
 pr<-fread(file.path(size_run,'tables/oof_predictions.csv'),colClasses=c(Sample='character'));full<-fig5_read('nested_oof_predictions')[Model==FIG5_PRIMARY];full[,Sample:=as.character(Sample)]
 setorder(full,Sample);ks<-sort(unique(pr$K));correct<-sapply(ks,function(k){x<-pr[K==k];x<-x[match(full$Sample,x$Sample)];stopifnot(identical(x$Truth,full$Truth));as.numeric(x$Prediction==x$Truth)})
 base<-as.numeric(full$Truth==full$Prediction);gr<-split(seq_len(nrow(full)),interaction(full$Study,full$Truth,drop=TRUE));cls<-split(seq_len(nrow(full)),full$Truth)
 # R recycling is explicit across rows for the patient-matched difference matrix.
 delta<-sweep(correct,1,base,'-');getdiff<-function(ii)colMeans(do.call(rbind,lapply(c('C1','C2','C3'),function(cl)colMeans(delta[ii[full$Truth[ii]==cl],,drop=FALSE]))))
 set.seed(20261031);bt<-replicate(1000,getdiff(unlist(lapply(gr,function(ii)sample(ii,length(ii),replace=TRUE)),use.names=FALSE)))
 loss<-data.table(K=ks,Estimate=getdiff(seq_len(nrow(full))),Lower95=apply(bt,1,quantile,.025),Upper95=apply(bt,1,quantile,.975));fwrite(loss,opt_path('source_data/marker_size_paired_BA_difference.csv'))
 q<-ggplot(loss,aes(K,100*Estimate))+geom_hline(yintercept=0,linetype=2,colour='#9B6176')+geom_errorbar(aes(ymin=100*Lower95,ymax=100*Upper95),width=1.5,linewidth=.4)+geom_line(colour=FIG5_INK,linewidth=.5)+geom_point(colour=FIG5_INK,size=1.8)+scale_x_continuous(breaks=ks)+labs(title='Paired difference from the complete pipeline',x='Training-selected genes (K)',y='Balanced-accuracy difference\n(percentage points)')+fig5_theme()+theme(plot.title=element_text(size=8),axis.title=element_text(size=7))
 (p/q+plot_layout(heights=c(1,1)))+plot_annotation(title='Exploratory marker-count sensitivity',subtitle='Complete 114-gene development panel retained; no compact classifier selected',caption='Top: saved OOF estimates and conditional 95% intervals; dashed lines show the complete training-fold selection pipeline.\nBottom: 1,000 paired patient bootstraps stratified by study and reference class. Overlap with zero does not demonstrate equivalence.\nRanking is training-only method support, then limma F rank. Not random deletion from a fixed full-data list; no compact LOSO evaluation.',theme=fig5_theme()+theme(plot.caption=element_text(size=6.2)))
}

opt_heatmap <- function(development=FALSE){function(){
 d<-readRDS('8-Visualization/refinements/fig5_complete_20261030_v4/source_data/5F_display_data.rds');ht<-NULL
 display_names<-c(log2dat_2012Andersen='2012Andersen',log2dat_2013Sia='2013Sia',log2dat_2017Jusakual1='Jusakual\nP1',log2dat_2017Jusakual2='Jusakual\nP2',log2dat_2020Job='2020Job',log2dat_TIGER='TIGER')
 if(development){old<-readRDS('8-Visualization/refinements/fig5_marker_story_20261030_v2/source_data/display_data.rds');meta<-copy(old$meta);ix<-order(meta$Subtype,meta$Study,meta$Sample);meta<-meta[ix];mat<-t(old$z[ix,d$ref$Gene,drop=FALSE]);d$cohorts<-'Development';d$data<-list(Development=list(mat=mat,meta=meta));total<-nrow(meta);fwrite(data.table(Gene=rownames(mat),as.data.table(mat)),opt_path('source_data/development_expression_114.csv'));fwrite(meta,opt_path('source_data/development_sample_order.csv'))}else{
 total<-506
 widths<-rbindlist(lapply(d$cohorts,function(st)data.table(Cohort=st,N=ncol(d$data[[st]]$mat),WidthMM=133*ncol(d$data[[st]]$mat)/total,SampleWidthMM=133/total)))
 fwrite(widths,opt_path('source_data/external_column_widths.csv'))
 }
 for(cohort in d$cohorts){x<-d$data[[cohort]];first<-is.null(ht);last<-cohort==tail(d$cohorts,1)
  if(development){anno<-HeatmapAnnotation(Reference=x$meta$Subtype,Study=x$meta$Study,col=list(Reference=FIG5_SUBTYPE,Study=FIG5_STUDY),show_annotation_name=TRUE,simple_anno_size=unit(2.4,'mm'),annotation_name_gp=gpar(fontfamily='Arial',fontsize=6),annotation_legend_param=list(Reference=list(title_gp=gpar(fontsize=6),labels_gp=gpar(fontsize=6),nrow=1),Study=list(title_gp=gpar(fontsize=6),labels_gp=gpar(fontsize=6),nrow=2)));title<-NULL;nm<-'Training z'}else{anno<-HeatmapAnnotation(Projected=x$meta$Prediction,col=list(Projected=FIG5_SUBTYPE),show_annotation_name=FALSE,simple_anno_size=unit(2.4,'mm'),show_legend=first,annotation_legend_param=list(title_gp=gpar(fontfamily='Arial',fontsize=7,fontface='bold'),labels_gp=gpar(fontfamily='Arial',fontsize=6.5),nrow=1));title<-sprintf('%s\nn=%d\n%d/114',display_names[cohort],nrow(x$meta),unique(x$meta$MatchedFeatures));nm<-'Within-study z'}
  one<-Heatmap(x$mat,name=if(first)nm else cohort,col=circlize::colorRamp2(c(-2.5,0,2.5),STEP8_CONFIG$effect_colors),na_col='#B7B7B7',cluster_rows=FALSE,cluster_columns=FALSE,row_split=d$ref$ReferenceClass,row_title_rot=0,row_title_gp=gpar(fontfamily='Arial',fontsize=7,fontface='bold'),show_row_names=last,row_names_side='right',row_names_gp=gpar(fontfamily='Arial',fontsize=5.5),row_names_max_width=unit(20,'mm'),show_column_names=FALSE,column_title=title,column_title_gp=gpar(fontfamily='Arial',fontsize=5.5,fontface='bold'),top_annotation=anno,width=unit(133*ncol(x$mat)/total,'mm'),height=unit(178,'mm'),row_gap=unit(1,'mm'),show_heatmap_legend=first,use_raster=FALSE,heatmap_legend_param=list(direction='horizontal',legend_width=unit(28,'mm'),title_gp=gpar(fontfamily='Arial',fontsize=6.5,fontface='bold'),labels_gp=gpar(fontfamily='Arial',fontsize=6),at=c(-2.5,0,2.5)))
  ht<-if(first)one else ht+one
 }
 draw(ht,heatmap_legend_side='bottom',annotation_legend_side='bottom',merge_legends=TRUE,ht_gap=unit(1.2,'mm'),padding=unit(c(24,3,28,3),'mm'))
 grid.text(if(development)'Expression of the complete classifier panel in development samples' else 'External microarray projection of the 114-gene classifier',x=unit(3,'mm'),y=unit(1,'npc')-unit(3,'mm'),just=c('left','top'),gp=gpar(fontfamily='Arial',fontsize=8.5,fontface='bold',col=FIG5_INK))
 grid.text(if(development)'Reference classes and study; supervised feature space, descriptive only' else '506 samples; constant sample-column width across all study/platform blocks',x=unit(3,'mm'),y=unit(1,'npc')-unit(8,'mm'),just=c('left','top'),gp=gpar(fontfamily='Arial',fontsize=6.8,col=FIG5_INK))
 grid.text(if(development)'All 114 genes; same row order as the external heatmap. Training z colors saturate at ±2.5.' else 'All 114 genes labelled. Grey: unmeasured. Within-study z colors saturate at ±2.5.',x=unit(3,'mm'),y=unit(7,'mm'),just='left',gp=gpar(fontfamily='Arial',fontsize=6.2,col=FIG5_INK))
 grid.text(if(development)'Rows grouped by highest development reference-class mean; not independent validation.' else 'Projected labels and expression share classifier features; this display does not measure external accuracy.',x=unit(3,'mm'),y=unit(3,'mm'),just='left',gp=gpar(fontfamily='Arial',fontsize=6.2,col=FIG5_INK))
}}

opt_survival <- function(){
 curves<-opt_read('external_survival_curves');risk<-opt_read('external_survival_at_risk');tests<-opt_read('external_survival_tests');counts<-opt_read('external_survival_class_counts')
 plots<-lapply(c('2017Jusakual','2020Job'),function(st){z<-curves[Study==st];rr<-risk[Study==st];tt<-tests[Study==st];cc<-counts[Study==st];xmax<-max(rr$Time);ticks<-sort(unique(rr$Time));pval<-function(p)if(p<.001)'<0.001' else sprintf('= %.3f',p)
  labels<-setNames(sprintf('%s (n=%d; events=%d)',cc$Prediction,cc$N,cc$Events),cc$Prediction)
  band<-rbindlist(lapply(split(z,z$Class),function(v){setorder(v,Time);b<-v[rep(seq_len(nrow(v)),each=2)];b[,Time:=as.vector(rbind(v$Time,c(v$Time[-1],tail(v$Time,1))))];b}))
  p<-ggplot(z,aes(Time,Survival,colour=Class))+geom_ribbon(data=band,aes(ymin=Lower95,ymax=Upper95,fill=Class),alpha=.08,colour=NA,show.legend=FALSE,na.rm=TRUE)+geom_step(linewidth=.65)+geom_point(data=z[Censored>0],shape=3,size=1.2,stroke=.35,show.legend=FALSE)+scale_colour_manual(values=FIG5_SUBTYPE,labels=labels)+scale_fill_manual(values=FIG5_SUBTYPE)+scale_x_continuous(limits=c(0,xmax),breaks=ticks,expand=expansion(mult=0))+scale_y_continuous(limits=c(0,1),breaks=c(0,.25,.5,.75,1),expand=expansion(mult=c(.01,.02)))+labs(title=sprintf('%s | n=%d; %d deaths',st,tt$N,tt$Events),subtitle=sprintf('Log-rank P %s; Holm P %s',pval(tt$P),pval(tt$P_Holm)),x='Overall survival (months)',y='Survival probability',colour=NULL)+fig5_theme()+theme(plot.title=element_text(size=8),plot.subtitle=element_text(size=6.5),axis.title=element_text(size=7),axis.text=element_text(size=6.5),legend.position='top',legend.direction='vertical',legend.text=element_text(size=6),legend.key.height=unit(2.2,'mm'),legend.margin=margin(0,0,0,0))
  rr[,Class:=factor(Class,levels=c('C3','C2','C1'))]
  q<-ggplot(rr,aes(Time,Class,label=AtRisk,colour=Class))+geom_text(size=2.15)+scale_colour_manual(values=FIG5_SUBTYPE,guide='none')+scale_x_continuous(limits=c(0,xmax),breaks=ticks,expand=expansion(mult=0))+labs(title='Number at risk',x=NULL,y=NULL)+fig5_theme()+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=6.5,margin=margin(r=8)),plot.title=element_text(size=6.5),plot.margin=margin(0,3,1,3,'mm'))+coord_cartesian(clip='off')
  p/q+plot_layout(heights=c(4,1))
 })
 wrap_plots(plots,ncol=2)+plot_annotation(title='Survival associations of externally projected subtypes',subtitle='Frozen 114-gene consensus SVM; each study evaluated separately',caption='Available-case overall survival; 95% pointwise log-log intervals; crosses mark censoring. Holm correction across two global tests.\n2017Jusakual: 81/83 records eligible; 2020Job: 75/78. Unadjusted prognostic associations, not subtype-accuracy validation.',theme=fig5_theme()+theme(plot.caption=element_text(size=6.2)))
}

render_fig5_optimized_panel <- function(id){
 stopifnot(id%in%names(opt_names));fig5_check()
 p<-switch(id,B=opt_evidence(),C=opt_comparison(),E=opt_heatmap(),F=opt_survival(),S11F=completion_strategy(),S11G=opt_marker_supplement(),S11I=opt_heatmap(TRUE))
 drawplot<-function(){if(is.function(p))p() else print(p)};sz<-opt_dims[[id]];stem<-opt_path('figures',opt_names[[id]])
 cairo_pdf(paste0(stem,'.pdf'),width=sz[1]/25.4,height=sz[2]/25.4,family='Arial');drawplot();dev.off()
 svglite(paste0(stem,'.svg'),width=sz[1]/25.4,height=sz[2]/25.4);drawplot();dev.off()
 ragg::agg_png(paste0(stem,'.png'),width=sz[1],height=sz[2],units='mm',res=300,background='white');drawplot();dev.off()
 stopifnot(system2('pdftoppm',c('-png','-singlefile','-r','150',shQuote(paste0(stem,'.pdf')),shQuote(opt_path('review/pdf_rendered',opt_names[[id]]))))==0)
 cat('Rendered',id,'\n');invisible(paste0(stem,'.pdf'))
}
render_fig5_optimized <- function(){for(id in names(opt_names))render_fig5_optimized_panel(id)}

# Current notebook entry: fixed A/D are returned without modifying their bytes or records.
render_fig5_current_panel <- function(id){
 mpath<-file.path(STEP8_ROOT,'8-Visualization/tables/figure_manifest.csv')
 m<-fread(mpath,colClasses='character',na.strings=NULL)
 ix<-which(m$Figure==if(startsWith(id,'S11'))'Supplement' else 'Figure 5')
 ix<-ix[m$Panel[ix]==id];stopifnot(length(ix)==1)
 target<-file.path(STEP8_ROOT,m$TargetOutput[ix])
 if(id%in%c('A','D')){stopifnot(file.exists(target));return(invisible(target))}
 if(id=='S11H'){
  src<-file.path(STEP8_ROOT,'8-Visualization/refinements/fig5_complete_20261030_v4/figures/Figure5E_study_transfer.pdf')
 }else src<-render_fig5_optimized_panel(id)
 for(ext in c('pdf','svg'))stopifnot(file.copy(sub('[.]pdf$',paste0('.',ext),src),sub('[.]pdf$',paste0('.',ext),target),overwrite=TRUE))
 sources<-strsplit(m$SourceTable[ix],';',fixed=TRUE)[[1]];stopifnot(all(file.exists(file.path(STEP8_ROOT,sources))))
 m[ix,`:=`(Status='Rendered',SourceMD5=paste(paste0(sources,'=',unname(tools::md5sum(file.path(STEP8_ROOT,sources)))),collapse=';'),QAFile='',RenderedAt=format(Sys.time(),tz='UTC',usetz=TRUE))]
 tmp<-tempfile(tmpdir=dirname(mpath),fileext='.csv');fwrite(m,tmp,na='');stopifnot(file.rename(tmp,mpath));invisible(target)
}
