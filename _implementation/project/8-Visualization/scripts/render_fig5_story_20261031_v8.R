# Figure5 v8: display-only adaptation of verified source tables.
if(!exists('fig5_v7_current_panel',inherits=FALSE)){
 source('8-Visualization/scripts/render_fig5_BCD_20261031_v7.R')
 fig5_v7_current_panel<-render_fig5_current_panel
}
V8<-'8-Visualization/refinements/fig5_story_20261031_v8'
v8_path<-function(...)file.path(STEP8_ROOT,V8,...)
V8_NAMES<-c(A='Figure5A_training_schematic',B='Figure5B_subtype_feature_evidence',C='Figure5C_subtype_gene_ranking',D='Figure5D_SVM_reduction',E='Figure5E_oof_performance',F='Figure5F_microarray_projection',G='Figure5G_external_survival',S11I='Supplement_S11I_aligned_expression',S11K='Supplement_S11K_SVM_RFE_detail')
V8_DIMS<-list(A=c(174,68),B=c(174,86),C=c(174,122),D=c(174,90),E=c(174,94),F=c(174,96),G=c(174,112),S11I=c(210,220),S11K=c(174,125))
v8_prepare<-function(){
 d<-fread(v7_path('source_data/B_subtype_gene_order.csv'))
 stopifnot(nrow(d)==114,!anyDuplicated(d$Gene),identical(d$DisplayOrder,1:114))
 fwrite(d,v8_path('source_data/B_subtype_gene_order.csv'))
 for(n in c('full_model_weight_ranking','marker_size_performance','paired_BA_difference'))
  file.copy(file.path(STEP8_ROOT,RFE_RUN,'tables',paste0(n,'.csv')),v8_path('source_data',paste0(n,'.csv')),overwrite=TRUE)
 x<-fig5_read('external_sample_predictions')[Model==FIG5_PRIMARY & Track=='Microarray']
 cohorts<-c('log2dat_2012Andersen','log2dat_2013Sia','log2dat_2017Jusakual1','log2dat_2017Jusakual2','log2dat_2020Job','log2dat_TIGER')
 labels<-c('2012Andersen','2013Sia','2017Jusakual P1','2017Jusakual P2','2020Job','TIGER')
 stopifnot(nrow(x)==506,all(x$IndependentTruthAvailable==FALSE),all(x$ProjectionOnly),!anyDuplicated(x[,.(Cohort,Sample)]))
 x[,CohortOrder:=match(Cohort,cohorts)];x[,CohortLabel:=labels[CohortOrder]]
 sm<-x[,.(N=.N,MatchedFeatures=unique(MatchedFeatures),MedianGap=median(Margin),Q1=quantile(Margin,.25),Q3=quantile(Margin,.75)),by=.(Cohort,CohortOrder,CohortLabel)]
 setorder(sm,CohortOrder);co<-x[,.N,by=.(Cohort,CohortOrder,CohortLabel,Prediction)]
 co[,Fraction:=N/sum(N),by=Cohort]
 fwrite(x,v8_path('source_data/microarray_samples.csv'));fwrite(sm,v8_path('source_data/microarray_summary.csv'));fwrite(co,v8_path('source_data/microarray_composition.csv'))
}

v8_B<-function(){function(){
 d<-fread(v8_path('source_data/B_subtype_gene_order.csv'))
 grid.newpage();pushViewport(viewport(xscale=c(0,174),yscale=c(0,86)))
 bc6_text('Consensus gene evidence',4,82,10,TRUE)
 x0<-23;cw<-(170-x0-4)/114;d[,X:=x0+(DisplayOrder-.5)*cw+(SubtypeGroup-1)*2]
 groups<-d[,.(N=.N,Left=min(X)-cw/2,Right=max(X)+cw/2),by=SubtypeGroup]
 for(j in 1:3){g<-groups[SubtypeGroup==j];mid<-(g$Left+g$Right)/2
  bc6_rect(mid,74,g$Right-g$Left,6,adjustcolor(FIG5_SUBTYPE[j],alpha.f=.15))
  bc6_text(sprintf('C%d  |  %d genes',j,g$N),mid,74,7.2,TRUE,'centre')
 }
 bc6_text('Method support',x0,66.5,7.5,TRUE)
 for(j in 1:3){yy<-62-(j-1)*3.6;bc6_text(c('limma','RF','LASSO')[j],x0-2,yy,7,FALSE,'right')
  for(i in 1:114)bc6_rect(d$X[i],yy,cw-.06,3,if(d[[c('LimmaFull','RFFull','LASSOFull')[j]]][i])FIG5_INK else '#E7EDF1')}
 bc6_text('Selection frequency (25 fits)',x0,48.5,7.5,TRUE)
 for(v in c(0,.5,1)){yy<-29+14*v;for(j in 1:3){g<-groups[SubtypeGroup==j];bc6_line(c(g$Left,g$Right),rep(yy,2),'grey90')};bc6_text(sprintf('%d%%',100*v),x0-2,yy,6.5,FALSE,'right')}
 for(i in 1:114){yy<-29+14*d$ConsensusFraction[i];bc6_line(rep(d$X[i],2),c(29,yy),'#A7BCCB',.6);grid.points(unit(d$X[i],'native'),unit(yy,'native'),pch=16,size=unit(.65,'mm'),gp=gpar(col=FIG5_INK))}
 bc6_text('Mean expression',x0,24,7.5,TRUE)
 ramp<-circlize::colorRamp2(c(-1.5,0,1.5),STEP8_CONFIG$effect_colors)
 for(j in 1:3){cl<-c('C1','C2','C3')[j];yy<-19.7-(j-1)*3.5;bc6_text(cl,x0-2,yy,7,TRUE,'right',FIG5_SUBTYPE[cl]);for(i in 1:114)bc6_rect(d$X[i],yy,cw-.06,3,ramp(d[[paste0(cl,'MeanTrainingZ')]][i]))}
 bc6_rect(24,5,1.8,1.8,FIG5_INK);bc6_text('Selected',27,5,6.5);bc6_rect(49,5,1.8,1.8,'#E7EDF1');bc6_text('Not selected',52,5,6.5)
 bc6_text('Mean z',111,5,6.5);for(i in 1:30)bc6_rect(132+(i-.5)*.7,5,.7,2,ramp(-1.5+(i-.5)*3/30));bc6_text('-1.5',130,5,6.5,FALSE,'right');bc6_text('1.5',155,5,6.5)
 popViewport()
}}

v8_D<-function(){
 w<-fread(v8_path('source_data/full_model_weight_ranking.csv'))[Rank<=15];w[,Gene:=factor(Gene,levels=rev(Gene))]
 p<-ggplot(w,aes(RelativeWeight*100,Gene))+geom_col(fill=FIG5_INK,width=.65)+scale_y_discrete(labels=function(x)ifelse(x=='C2','C2 (gene)',x))+
  labs(title='Gene weights',subtitle='Top 15 of the 114-gene model',x='Relative squared weight (%)',y=NULL)+fig5_theme()+theme(axis.text.y=element_text(face='italic',size=7),plot.title=element_text(size=8.5),plot.subtitle=element_text(size=7),axis.title=element_text(size=7))
 d<-fread(v8_path('source_data/marker_size_performance.csv'))[Metric=='BalancedAccuracy'];ref<-fig5_read('nested_performance_ci')[Model==FIG5_PRIMARY & Metric=='BalancedAccuracy',Estimate]
 q<-ggplot(d,aes(K,Estimate*100))+geom_ribbon(aes(ymin=Lower95*100,ymax=Upper95*100),fill='#DDE8EF',alpha=.85)+
  geom_hline(yintercept=100*ref,colour='#98647D',linetype=2,linewidth=.45)+geom_line(colour=FIG5_INK,linewidth=.6)+geom_point(colour=FIG5_INK,size=1.9)+
  annotate('text',x=78,y=100*ref+.6,label=sprintf('Full workflow: %.2f%%',100*ref),hjust=1,size=2.45,colour='#865269')+
  scale_x_continuous(breaks=c(2,10,20,40,60,80))+scale_y_continuous(breaks=seq(88,100,2))+
  labs(title='Performance as genes are removed',subtitle='Training-fold SVM recursive elimination',x='Retained genes',y='OOF balanced accuracy (%)')+fig5_theme()+theme(plot.title=element_text(size=8.5),plot.subtitle=element_text(size=7),axis.title=element_text(size=7),axis.text=element_text(size=7))
 (p|q)+plot_layout(widths=c(.85,1.4))+plot_annotation(title='Gene contribution and feature reduction',caption='Shading: 95% conditional intervals. Nine prespecified sizes; no optimal size selected.',theme=fig5_theme()+theme(plot.caption=element_text(size=6.5)))
}

v8_F<-function(){
 x<-fread(v8_path('source_data/microarray_samples.csv'));s<-fread(v8_path('source_data/microarray_summary.csv'));d<-fread(v8_path('source_data/microarray_composition.csv'))
 s[,Y:=7-CohortOrder];d[,Y:=7-CohortOrder];x[,Y:=7-CohortOrder]
 yscale<-function(labels=NULL)scale_y_continuous(breaks=s$Y,labels=if(is.null(labels))rep('',6) else labels,limits=c(.4,6.6),expand=c(0,0))
 base<-fig5_theme()+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.title=element_text(size=7),axis.text=element_text(size=7),plot.title=element_text(size=8),plot.subtitle=element_text(size=6.5))
 cover<-ggplot(s,aes(x=MatchedFeatures/114,y=Y))+geom_segment(aes(x=0,xend=1,yend=Y),linewidth=5,colour='#E8EDF1')+geom_segment(aes(x=0,xend=MatchedFeatures/114,yend=Y),linewidth=5,colour='#90A6B5')+
  geom_text(aes(x=.49,label=paste0(MatchedFeatures,' / 114')),size=2.35,colour=FIG5_INK)+yscale(sprintf('%s  (n=%d)',s$CohortLabel,s$N))+scale_x_continuous(breaks=c(0,.5,1),labels=c('0','50','100'),limits=c(0,1))+
  labs(title='Gene coverage',x='Available genes (%)',y=NULL)+base
 d[,Prediction:=factor(Prediction,levels=c('C1','C2','C3'))];setorder(d,CohortOrder,Prediction);d[,Center:=cumsum(Fraction)-Fraction/2,by=Cohort]
 comp<-ggplot(d,aes(Fraction,Y,fill=Prediction))+geom_col(orientation='y',width=.57,position=position_stack(reverse=TRUE))+geom_text(aes(x=Center,label=N),size=2.2,colour='white')+yscale()+scale_fill_manual(values=FIG5_SUBTYPE)+scale_x_continuous(breaks=c(0,.5,1),labels=c('0','50','100'),limits=c(0,1))+
  labs(title='Projected subtypes',x='Samples (%)',y=NULL,fill=NULL)+base+theme(legend.position='bottom',legend.text=element_text(size=6.5),legend.key.size=unit(2.5,'mm'))
 # Deterministic vertical offsets preserve every sample without random jitter.
 setorder(x,CohortOrder,Sample);x[,Offset:=((seq_len(.N)-1)%%9-4)*.035,by=Cohort]
 gap<-ggplot(x,aes(Margin,Y+Offset))+geom_point(size=.7,colour='#9DAFBC',alpha=.55)+geom_segment(data=s,aes(x=Q1,xend=Q3,y=Y,yend=Y),inherit.aes=FALSE,linewidth=1,colour=FIG5_INK)+geom_point(data=s,aes(MedianGap,Y),inherit.aes=FALSE,size=1.7,colour=FIG5_INK)+yscale()+scale_x_continuous(breaks=c(0,.5,1),limits=c(0,1.02))+
  labs(title='Probability gap',x='Largest − second largest',y=NULL)+base
 (cover|comp|gap)+plot_layout(widths=c(.95,1.2,1.1))+plot_annotation(title='Microarray projection',subtitle='506 samples across six study/platform blocks',caption='Gap: all samples, median and interquartile range. No independent subtype labels; projection only.',theme=fig5_theme()+theme(plot.subtitle=element_text(size=7),plot.caption=element_text(size=6.5)))
}

v8_heatmap<-function(){function(){
 ext<-readRDS(file.path(STEP8_ROOT,'8-Visualization/refinements/fig5_complete_20261030_v4/source_data/5F_display_data.rds'))
 dev<-readRDS(file.path(STEP8_ROOT,'8-Visualization/refinements/fig5_marker_story_20261030_v2/source_data/display_data.rds'))
 meta<-copy(dev$meta);ix<-order(meta$Subtype,meta$Study,meta$Sample);meta<-meta[ix];mat<-t(dev$z[ix,ext$ref$Gene,drop=FALSE])
 stopifnot(nrow(mat)==114,ncol(mat)==497,identical(rownames(mat),ext$ref$Gene))
 fwrite(data.table(Gene=rownames(mat),as.data.table(mat)),v8_path('source_data/development_expression_114.csv'));fwrite(meta,v8_path('source_data/development_sample_order.csv'));fwrite(ext$ref,v8_path('source_data/heatmap_gene_order.csv'))
 all<-c(list(Development=list(mat=mat,meta=meta)),ext$data[ext$cohorts]);total<-sum(vapply(all,function(z)ncol(z$mat),integer(1)));stopifnot(total==1003)
 widths<-data.table(Block=names(all),N=vapply(all,function(z)ncol(z$mat),integer(1)));widths[,WidthMM:=138*N/total];fwrite(widths,v8_path('source_data/heatmap_column_widths.csv'))
 ht<-NULL;ramp<-circlize::colorRamp2(c(-2.5,0,2.5),STEP8_CONFIG$effect_colors)
 gp<-gpar(fontfamily='Arial',fontsize=6.5)
 study_keys<-c('2012Andersen','2013Sia','2017Jusakual','2017Jusakual','2020Job','TIGER-LC')
 study_labels<-c('1  Andersen','2  Sia','3  Jusakual P1','4  Jusakual P2','5  Job','6  TIGER')
 stopifnot(all(study_keys%in%names(FIG5_STUDY)))
 tracks<-rbindlist(lapply(seq_along(all),function(j){z<-all[[j]];data.table(Block=names(all)[j],Sample=as.character(z$meta$Sample),Subtype=as.character(if(j==1)z$meta$Subtype else z$meta$Prediction),Study=if(j==1)as.character(z$meta$Study) else rep(study_keys[j-1],ncol(z$mat)),SubtypeSource=if(j==1)'Reference' else 'Projected',TrackHeightMM=2)}))
 fwrite(tracks,v8_path('source_data/heatmap_annotation_tracks.csv'))
 for(j in seq_along(all)){z<-all[[j]];stopifnot(identical(rownames(z$mat),ext$ref$Gene));isdev<-j==1
  if(isdev){anno<-HeatmapAnnotation(Reference=z$meta$Subtype,Study=z$meta$Study,col=list(Reference=FIG5_SUBTYPE,Study=FIG5_STUDY),show_annotation_name=FALSE,simple_anno_size=unit(2,'mm'),gap=unit(.5,'mm'),show_legend=FALSE);title<-'Development\nn = 497';nm<-'Training z'}else{
   anno<-HeatmapAnnotation(Projected=z$meta$Prediction,Study=rep(study_keys[j-1],ncol(z$mat)),col=list(Projected=FIG5_SUBTYPE,Study=FIG5_STUDY),show_annotation_name=FALSE,simple_anno_size=unit(2,'mm'),gap=unit(.5,'mm'),show_legend=FALSE)
   title<-sprintf('%d\nn=%d',j-1,ncol(z$mat));nm<-if(j==2)'Within-study z' else names(all)[j]
  }
  one<-Heatmap(z$mat,name=nm,col=ramp,na_col='#B7B7B7',cluster_rows=FALSE,cluster_columns=FALSE,row_split=ext$ref$ReferenceClass,row_title=NULL,show_row_names=isdev,row_names_side='left',row_names_gp=gpar(fontfamily='Arial',fontsize=5.5),row_names_max_width=unit(18,'mm'),show_column_names=FALSE,column_title=title,column_title_gp=gpar(fontfamily='Arial',fontsize=6,fontface='bold'),top_annotation=anno,width=unit(widths$WidthMM[j],'mm'),height=unit(186,'mm'),row_gap=unit(1,'mm'),show_heatmap_legend=FALSE,use_raster=FALSE,heatmap_legend_param=list(direction='horizontal',legend_width=unit(25,'mm'),title_gp=gp,labels_gp=gp,at=c(-2.5,0,2.5)))
  ht<-if(is.null(ht))one else ht+one
 }
 # One compact right-side legend column; both panels have equal-height subtype/Study tracks.
 dev_studies<-names(FIG5_STUDY)[names(FIG5_STUDY)%in%unique(meta$Study)]
 lg<-list(
  Legend(title='Subtype',at=names(FIG5_SUBTYPE),legend_gp=gpar(fill=FIG5_SUBTYPE),title_gp=gp,labels_gp=gp,grid_height=unit(3,'mm')),
  Legend(title='Development study',at=dev_studies,legend_gp=gpar(fill=FIG5_STUDY[dev_studies]),title_gp=gp,labels_gp=gp,grid_height=unit(3,'mm')),
  Legend(title='Microarray study',at=study_labels,legend_gp=gpar(fill=FIG5_STUDY[study_keys]),title_gp=gp,labels_gp=gp,grid_height=unit(3,'mm')),
  Legend(title='Training z',col_fun=ramp,at=c(-2.5,0,2.5),direction='vertical',legend_height=unit(20,'mm'),title_gp=gp,labels_gp=gp),
  Legend(title='Within-study z',col_fun=ramp,at=c(-2.5,0,2.5),direction='vertical',legend_height=unit(20,'mm'),title_gp=gp,labels_gp=gp),
  Legend(title='Missing',at='Unmeasured',legend_gp=gpar(fill='#B7B7B7'),title_gp=gp,labels_gp=gp,grid_height=unit(3,'mm')))
 draw(ht,heatmap_legend_side='right',annotation_legend_side='right',annotation_legend_list=lg,merge_legends=TRUE,ht_gap=unit(c(3,rep(.8,5)),'mm'),padding=unit(c(3,2,20,2),'mm'))
 grid.text('Classifier-gene expression across development and microarray cohorts',x=unit(3,'mm'),y=unit(217,'mm'),just=c('left','top'),gp=gpar(fontfamily='Arial',fontsize=9,fontface='bold',col=FIG5_INK))
 grid.text('Microarray',x=unit(135.5,'mm'),y=unit(204,'mm'),just=c('centre','top'),gp=gpar(fontfamily='Arial',fontsize=7,fontface='bold',col=FIG5_INK))
}}

render_fig5_v8_panel<-function(id){
 stopifnot(id%in%names(V8_NAMES));fig5_check();v8_prepare();stem<-v8_path('figures',V8_NAMES[[id]])
 old<-c(A='8-Visualization/figures/main/Figure5A_training_schematic',C='8-Visualization/figures/main/Figure5C_subtype_gene_ranking',E='8-Visualization/refinements/fig5_BCD_20261031_v7/figures/Supplement_S11K_oof_performance',G='8-Visualization/refinements/fig5_optimized_20261031_v5/figures/Figure5F_external_survival',S11K='8-Visualization/refinements/fig5_BCD_20261031_v7/figures/Figure5D_SVM_RFE')
 if(id%in%names(old)){for(ext in c('pdf','svg'))stopifnot(file.copy(file.path(STEP8_ROOT,paste0(old[[id]],'.',ext)),paste0(stem,'.',ext),overwrite=TRUE))}else{
  p<-switch(id,B=v8_B(),D=v8_D(),F=v8_F(),S11I=v8_heatmap());sz<-V8_DIMS[[id]];drawplot<-function(){if(is.function(p))p() else print(p)}
  cairo_pdf(paste0(stem,'.pdf'),width=sz[1]/25.4,height=sz[2]/25.4,family='Arial');drawplot();dev.off()
  svglite::svglite(paste0(stem,'.svg'),width=sz[1]/25.4,height=sz[2]/25.4);drawplot();dev.off()
  ragg::agg_png(paste0(stem,'.png'),width=sz[1],height=sz[2],units='mm',res=300,background='white');drawplot();dev.off()
 }
 stopifnot(system2('pdftoppm',c('-png','-singlefile','-r','150',shQuote(paste0(stem,'.pdf')),shQuote(v8_path('review/pdf_rendered',V8_NAMES[[id]]))))==0)
 cat('Rendered',id,'\n');invisible(paste0(stem,'.pdf'))
}

render_fig5_current_panel<-function(id){
 if(!id%in%names(V8_NAMES))return(fig5_v7_current_panel(id))
 mpath<-file.path(STEP8_ROOT,'8-Visualization/tables/figure_manifest.csv');m<-fread(mpath,colClasses='character',na.strings=NULL)
 ix<-which(m$Figure==if(startsWith(id,'S11'))'Supplement' else 'Figure 5');ix<-ix[m$Panel[ix]==id];stopifnot(length(ix)==1)
 target<-file.path(STEP8_ROOT,m$TargetOutput[ix]);if(id%in%c('A','C'))return(invisible(target))
 src<-render_fig5_v8_panel(id);for(ext in c('pdf','svg'))stopifnot(file.copy(sub('[.]pdf$',paste0('.',ext),src),sub('[.]pdf$',paste0('.',ext),target),overwrite=TRUE))
 sources<-strsplit(m$SourceTable[ix],';',fixed=TRUE)[[1]];stopifnot(all(file.exists(file.path(STEP8_ROOT,sources))))
 m[ix,`:=`(Status='Rendered',SourceMD5=paste(paste0(sources,'=',unname(tools::md5sum(file.path(STEP8_ROOT,sources)))),collapse=';'),QAFile='',RenderedAt=format(Sys.time(),tz='UTC',usetz=TRUE))]
 tmp<-tempfile(tmpdir=dirname(mpath),fileext='.csv');fwrite(m,tmp,na='');stopifnot(file.rename(tmp,mpath));invisible(target)
}
