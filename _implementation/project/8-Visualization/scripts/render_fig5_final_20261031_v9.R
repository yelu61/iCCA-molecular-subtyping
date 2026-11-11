# Final-size Figure5 revisions; scientific definitions are in the external legends.
if(!exists('fig5_v8_current_panel',inherits=FALSE)){
 source('8-Visualization/scripts/render_fig5_story_20261031_v8.R')
 fig5_v8_current_panel<-render_fig5_current_panel
}
V9<-'8-Visualization/refinements/fig5_final_layout_20261031_v9';v9_path<-function(...)file.path(STEP8_ROOT,V9,...)
MINRUN<-'5-ClassifierValidation/runs/bulk-rnaseq__20261031_step5_minpanel_v1'
mn<-function(n)fread(file.path(STEP8_ROOT,MINRUN,'tables',paste0(n,'.csv')))
completion_path<-function(...)v9_path(...)
v9_theme<-function(){fig5_theme()+theme(plot.title=element_text(size=8,face='bold'),plot.subtitle=element_text(size=6.5),axis.title=element_text(size=7),axis.text=element_text(size=6.5),legend.text=element_text(size=6.5),legend.title=element_text(size=6.5),plot.caption=element_text(size=6.5),plot.margin=margin(2,2,2,2,'mm'))}
strip_annotation<-function(p,title){p$patches$annotation$title<-title;p$patches$annotation$subtitle<-NULL;p$patches$annotation$caption<-NULL;p}
v9_B<-function(){f<-v8_B();s<-paste(deparse(body(f)),collapse='\n');s<-gsub('"Mean z", 111','"Mean z", 101',s,fixed=TRUE);body(f)<-parse(text=s)[[1]];f}
v9_D<-function(){
 stopifnot(isTRUE(fromJSON(file.path(STEP8_ROOT,MINRUN,'review/independent_verification.json'))$passed))
 w<-mn('full_candidate_weights');n<-nrow(w);w[,Gene:=factor(Gene,levels=rev(Gene))]
 left<-ggplot(w,aes(100*RelativeWeight,Gene))+geom_col(fill=FIG5_INK,width=.65)+labs(title=sprintf('%d-gene candidate',n),subtitle='Full-development weights',x='Relative weight (%)',y=NULL)+v9_theme()+theme(axis.text.y=element_text(face='italic',size=6.5))
 d<-mn('performance_ci')[Metric%in%c('BalancedAccuracy','C1Recall','C2Recall','C3Recall')];d[,Y:=match(Metric,c('C3Recall','C2Recall','C1Recall','BalancedAccuracy'))+ifelse(Model=='Full',.12,-.12)]
 right<-ggplot(d,aes(100*Estimate,Y,shape=Model,colour=Model))+geom_errorbar(aes(xmin=100*Lower95,xmax=100*Upper95),orientation='y',width=.1,linewidth=.35)+geom_point(size=1.6)+
  scale_y_continuous(breaks=1:4,labels=c('C3','C2','C1','BA'),limits=c(.5,4.5))+scale_x_continuous(breaks=c(94,97,100),limits=c(93,100.4))+scale_colour_manual(values=c(Full='#8197A6',Compact=FIG5_INK),labels=c(Full='Full',Compact='Inner-selected'))+scale_shape_manual(values=c(Full=1,Compact=16),labels=c(Full='Full',Compact='Inner-selected'))+
  labs(title='Outer OOF',subtitle='BA and class recall',x='Performance (%)',y=NULL,shape=NULL,colour=NULL)+v9_theme()+theme(legend.position='bottom',legend.direction='vertical',legend.key.height=unit(2,'mm'),axis.line.y=element_blank(),axis.ticks.y=element_blank())
 sizes<-mn('choice')$K;delta<-mn('paired_differences')[Metric=='BalancedAccuracy',Estimate]
 (left|right)+plot_layout(widths=c(1,1.05))+plot_annotation(title='Criterion-guided feature reduction',subtitle=sprintf('Fold-selected genes: %g [%g–%g]',median(sizes),min(sizes),max(sizes)),caption=sprintf('OOF BA loss %.2f pp; exceeds the 1 pp target.',-delta*100),theme=v9_theme())
}
v9_E<-function(){
 r<-fig5_read('roc_coordinates')[Model==FIG5_PRIMARY];r[,FPR:=1-Specificity];m<-fig5_read('nested_performance_ci')[Model==FIG5_PRIMARY]
 labels<-setNames(vapply(c('C1','C2','C3'),function(cl)sprintf('%s %.4f',cl,m[Metric==paste0(cl,'_AUROC'),Estimate]),character(1)),c('C1','C2','C3'))
 p<-ggplot(r,aes(FPR,Sensitivity,colour=Class))+geom_abline(slope=1,intercept=0,colour='grey80',linetype=2,linewidth=.25)+geom_step(linewidth=.5)+scale_colour_manual(values=FIG5_SUBTYPE,labels=labels)+scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1))+scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1))+coord_equal()+labs(title='ROC',x='False-positive rate',y='True-positive rate',colour='AUC')+v9_theme()+theme(legend.position='bottom',legend.direction='vertical',legend.key.height=unit(2,'mm'),legend.key.width=unit(3,'mm'))
 cm<-fig5_read('confusion_matrix')[Model==FIG5_PRIMARY];cm[,Truth:=factor(Truth,levels=c('C3','C2','C1'))];cm[,Prediction:=factor(Prediction,levels=c('C1','C2','C3'))]
 q<-ggplot(cm,aes(Prediction,Truth,fill=RowFraction))+geom_tile(colour='white',linewidth=.5)+geom_text(aes(label=sprintf('%d\n%.1f%%',Freq,100*RowFraction),colour=ifelse(RowFraction>.5,'white',FIG5_INK)),size=2.3,lineheight=.95)+scale_colour_identity()+scale_fill_gradient(low='#F0F5F8',high='#4E8DB1',limits=c(0,1),guide='none')+coord_equal()+labs(title='Confusion matrix',x='Predicted',y='Reference')+v9_theme()+theme(axis.line=element_blank(),axis.ticks=element_blank())
 (p|q)+plot_layout(widths=c(.85,1.2))+plot_annotation(title='114-gene classifier: internal OOF',subtitle=sprintf('n = 497  |  Accuracy %.2f%%',100*m[Metric=='Accuracy',Estimate]),theme=v9_theme())
}
v9_F<-function(){
 x<-fread(v8_path('source_data/microarray_samples.csv'));s<-fread(v8_path('source_data/microarray_summary.csv'));d<-fread(v8_path('source_data/microarray_composition.csv'))
 short<-c('Andersen','Sia','Jusakual P1','Jusakual P2','Job','TIGER');s[,Y:=7-CohortOrder];d[,Y:=7-CohortOrder];x[,Y:=7-CohortOrder]
 ys<-function(labels=NULL)scale_y_continuous(breaks=s$Y,labels=if(is.null(labels))rep('',6) else labels,limits=c(.4,6.6),expand=c(0,0))
 base<-v9_theme()+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.title=element_text(size=7.5))
 p<-ggplot(s,aes(MatchedFeatures/114,Y))+geom_segment(aes(x=0,xend=1,yend=Y),linewidth=4,colour='#E8EDF1')+geom_segment(aes(x=0,xend=MatchedFeatures/114,yend=Y),linewidth=4,colour='#90A6B5')+geom_text(aes(x=.49,label=MatchedFeatures),size=2.2)+ys(sprintf('%s (%d)',short,s$N))+scale_x_continuous(breaks=c(0,1),labels=c('0','100'),limits=c(0,1))+labs(title='Coverage',x='Genes (%)',y=NULL)+base
 d[,Prediction:=factor(Prediction,levels=c('C1','C2','C3'))];setorder(d,CohortOrder,Prediction);d[,Center:=cumsum(Fraction)-Fraction/2,by=Cohort]
 q<-ggplot(d,aes(Fraction,Y,fill=Prediction))+geom_col(orientation='y',width=.5,position=position_stack(reverse=TRUE))+geom_text(aes(x=Center,label=N),size=2.05,colour='white')+ys()+scale_fill_manual(values=FIG5_SUBTYPE)+scale_x_continuous(limits=c(0,1),breaks=c(0,1),labels=c('0','100'))+labs(title='Subtype',x='Samples (%)',y=NULL,fill=NULL)+base+theme(legend.position='bottom',legend.key.size=unit(2,'mm'),legend.spacing.x=unit(1,'mm'))
 setorder(x,CohortOrder,Sample);x[,Offset:=((seq_len(.N)-1)%%9-4)*.035,by=Cohort]
 g<-ggplot(x,aes(Margin,Y+Offset))+geom_point(size=.6,colour='#9DAFBC',alpha=.55)+geom_segment(data=s,aes(x=Q1,xend=Q3,y=Y,yend=Y),inherit.aes=FALSE,linewidth=.8,colour=FIG5_INK)+geom_point(data=s,aes(MedianGap,Y),inherit.aes=FALSE,size=1.4,colour=FIG5_INK)+ys()+scale_x_continuous(limits=c(0,1.02),breaks=c(0,1))+labs(title='Gap',x='Probability',y=NULL)+base
 (p|q|g)+plot_layout(widths=c(.8,1.1,.9))+plot_annotation(title='Microarray projection (n = 506)',caption='Projection only; no independent subtype labels.',theme=v9_theme())
}
v9_G<-function(){
 curves<-opt_read('external_survival_curves');risk<-opt_read('external_survival_at_risk');tests<-opt_read('external_survival_tests')
 ps<-lapply(c('2017Jusakual','2020Job'),function(st){z<-curves[Study==st];rr<-risk[Study==st];tt<-tests[Study==st];xmax<-max(rr$Time);ticks<-seq(0,xmax,by=48);rr<-rr[Time%in%ticks];band<-rbindlist(lapply(split(z,z$Class),function(v){setorder(v,Time);b<-v[rep(seq_len(nrow(v)),each=2)];b[,Time:=as.vector(rbind(v$Time,c(v$Time[-1],tail(v$Time,1))))];b}))
  p<-ggplot(z,aes(Time,Survival,colour=Class))+geom_ribbon(data=band,aes(ymin=Lower95,ymax=Upper95,fill=Class),alpha=.08,colour=NA,show.legend=FALSE,na.rm=TRUE)+geom_step(linewidth=.5)+geom_point(data=z[Censored>0],shape=3,size=.9,stroke=.3,show.legend=FALSE)+scale_colour_manual(values=FIG5_SUBTYPE)+scale_fill_manual(values=FIG5_SUBTYPE)+scale_x_continuous(limits=c(0,xmax),breaks=ticks,expand=expansion(mult=c(.035,.045)))+scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1),expand=expansion(mult=c(.01,.02)))+labs(title=sprintf('%s (n=%d)',if(st=='2017Jusakual')'Jusakual' else 'Job',tt$N),subtitle=sprintf('Holm P = %.3f',tt$P_Holm),x='OS (months)',y='Survival',colour=NULL)+v9_theme()+theme(legend.position='top',legend.key.width=unit(2,'mm'),legend.spacing.x=unit(1,'mm'),legend.key.height=unit(2,'mm'),plot.margin=margin(2,1,1,1,'mm'))
  rr[,Class:=factor(Class,levels=c('C3','C2','C1'))]
  q<-ggplot(rr,aes(Time,Class,label=AtRisk,colour=Class))+geom_text(size=2.15)+scale_colour_manual(values=FIG5_SUBTYPE,guide='none')+scale_x_continuous(limits=c(0,xmax),breaks=ticks,expand=expansion(mult=c(.035,.045)))+labs(title='At risk',x=NULL,y=NULL)+v9_theme()+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_blank(),plot.title=element_text(size=6.5),plot.margin=margin(0,1,1,1,'mm'))+coord_cartesian(clip='off')
  p/q+plot_layout(heights=c(4,1))
 })
 wrap_plots(ps,ncol=2)+plot_annotation(title='External survival associations',theme=v9_theme())
}
v9_SA<-function(){
 d<-fread('5-ClassifierValidation/runs/bulk-rnaseq__20261031_step5_threshold_grid_v1/tables/threshold_summary.csv');d[,Y:=rev(seq_len(.N))];d[,Reference:=Spec=='L30_RF200']
 p<-ggplot(d,aes(BADifference*100,Y))+geom_vline(xintercept=0,linetype=2,colour='grey65',linewidth=.3)+geom_errorbar(aes(xmin=BADifferenceLower95*100,xmax=BADifferenceUpper95*100),orientation='y',width=.15,linewidth=.4,colour='#667E8F')+geom_point(size=1.6,colour=FIG5_INK)+scale_y_continuous(breaks=d$Y,labels=paste(d$LimmaTop,d$RFTop,sep=' / '),limits=c(.5,10),expand=c(0,0))+labs(title='Paired difference',x='Δ BA (pp)',y='limma / RF')+v9_theme()
 q<-ggplot(d,aes(y=Y))+geom_text(aes(x=1,label=FullGenes),size=2.3)+geom_text(aes(x=2.5,label=sprintf('%.2f',FullJaccard)),size=2.3)+annotate('text',x=c(1,2.5),y=9.8,label=c('Genes','Jaccard'),size=2.2,fontface='bold')+scale_y_continuous(limits=c(.5,10),expand=c(0,0))+scale_x_continuous(limits=c(0,3.3),expand=c(0,0))+theme_void()+theme(plot.margin=margin(2,2,2,2,'mm'))
 (p|q)+plot_layout(widths=c(1.5,1))+plot_annotation(title='Screening-threshold sensitivity',theme=v9_theme())
}
v9_SB<-function(){function(){
 d<-readRDS('8-Visualization/refinements/fig5_marker_story_20261030_v2/source_data/display_data.rds');grid.newpage();pushViewport(viewport(xscale=c(0,1),yscale=c(0,1)))
 tx<-function(t,x,y,size=8,bold=FALSE)grid.text(t,x,y,gp=gpar(fontfamily='Arial',fontsize=size,fontface=if(bold)'bold' else 'plain',col=FIG5_INK))
 tx('Three-method overlap',.5,.96,9,TRUE);fill<-c('#A8C5D6','#BDD5C4','#D9C3CF');centres<-list(c(.37,.64),c(.63,.64),c(.5,.43))
 for(i in 1:3)grid.circle(centres[[i]][1],centres[[i]][2],r=.268,gp=gpar(fill=adjustcolor(fill[i],alpha.f=.55),col=FIG5_INK,lwd=.6))
 get<-function(l,r,a)d$regions[Limma==l & RF==r & LASSO==a,N]
 tx('limma\n129',.18,.88,7,TRUE);tx('RF\n200',.82,.88,7,TRUE);tx('LASSO: 95',.5,.13,7,TRUE)
 tx(get(TRUE,FALSE,FALSE),.25,.66);tx(get(FALSE,TRUE,FALSE),.75,.66);tx(get(FALSE,FALSE,TRUE),.5,.28);tx(get(TRUE,TRUE,FALSE),.5,.76);tx(get(TRUE,FALSE,TRUE),.35,.46);tx(get(FALSE,TRUE,TRUE),.65,.46);tx(get(TRUE,TRUE,TRUE),.5,.56,9,TRUE)
 tx('≥2 methods: 114 genes',.5,.065,7,TRUE);popViewport()
}}
v9_SC<-function(){p<-completion_transfer();strip_annotation(p,'Study-held-out performance')}
v9_SD<-function(){completion_strategy()+labs(title='Feature-selection strategies',subtitle=NULL,caption=NULL)+theme(plot.title=element_text(size=9),strip.text=element_text(size=7))}
v9_SE<-function(){
 d<-fread(file.path(RFE_RUN,'tables/marker_size_performance.csv'))[Metric%in%c('BalancedAccuracy','MacroF1','MacroAUROC')];labels<-c(BalancedAccuracy='Balanced accuracy',MacroF1='Macro-F1',MacroAUROC='Macro-AUROC');d[,Label:=factor(labels[Metric],levels=unname(labels))]
 p<-ggplot(d,aes(K,Estimate))+geom_ribbon(aes(ymin=Lower95,ymax=Upper95),fill='#DDE8EF')+geom_line(colour=FIG5_INK,linewidth=.5)+geom_point(size=1.5,colour=FIG5_INK)+facet_wrap(~Label,nrow=1,scales='free_y')+scale_x_continuous(breaks=c(2,20,40,60,80))+labs(x='Retained genes',y='OOF performance')+v9_theme()
 z<-fread(file.path(RFE_RUN,'tables/paired_BA_difference.csv'));q<-ggplot(z,aes(K,100*Estimate))+geom_hline(yintercept=0,linetype=2,colour='grey65')+geom_errorbar(aes(ymin=100*Lower95,ymax=100*Upper95),width=1.5,linewidth=.35)+geom_line(colour=FIG5_INK,linewidth=.45)+geom_point(size=1.5,colour=FIG5_INK)+scale_x_continuous(breaks=c(2,10,20,40,60,80))+labs(x='Retained genes',y='Δ BA (pp)')+v9_theme()
 (p/q)+plot_layout(heights=c(1,1))+plot_annotation(title='RFE sensitivity and performance loss',theme=v9_theme())
}
V9_NAMES<-c(A='Figure5A_training_schematic',B='Figure5B_subtype_feature_evidence',C='Figure5C_subtype_gene_ranking',D='Figure5D_criterion_guided_reduction',E='Figure5E_oof_performance',F='Figure5F_microarray_projection',G='Figure5G_external_survival',SA='SupplementaryFigure5A_threshold_sensitivity',SB='SupplementaryFigure5B_method_overlap',SC='SupplementaryFigure5C_study_transfer',SD='SupplementaryFigure5D_feature_strategies',SE='SupplementaryFigure5E_RFE_sensitivity',SF='SupplementaryFigure5F_aligned_expression')
V9_DIMS<-list(A=c(174,68),B=c(90,80),C=c(174,122),D=c(90,100),E=c(90,82),F=c(90,84),G=c(90,96),SA=c(90,82),SB=c(90,76),SC=c(174,88),SD=c(174,82),SE=c(174,104),SF=c(210,220))
render_fig5_v9_panel<-function(id){
 stopifnot(id%in%names(V9_NAMES));fig5_check();stem<-v9_path('figures',V9_NAMES[[id]]);sz<-V9_DIMS[[id]]
 if(id%in%c('A','C','SF')){old<-if(id=='SF')v8_path('figures','Supplement_S11I_aligned_expression') else v8_path('figures',V9_NAMES[[id]]);for(ext in c('pdf','svg'))stopifnot(file.copy(paste0(old,'.',ext),paste0(stem,'.',ext),overwrite=TRUE))}else{
  p<-switch(id,B=v9_B(),D=v9_D(),E=v9_E(),F=v9_F(),G=v9_G(),SA=v9_SA(),SB=v9_SB(),SC=v9_SC(),SD=v9_SD(),SE=v9_SE());drawplot<-function(){if(is.function(p))p() else print(p)}
  cairo_pdf(paste0(stem,'.pdf'),width=sz[1]/25.4,height=sz[2]/25.4,family='Arial');drawplot();dev.off();svglite::svglite(paste0(stem,'.svg'),width=sz[1]/25.4,height=sz[2]/25.4);drawplot();dev.off()
  ragg::agg_png(paste0(stem,'.png'),width=sz[1],height=sz[2],units='mm',res=300,background='white');drawplot();dev.off()
 }
 stopifnot(system2('pdftoppm',c('-png','-singlefile','-r','150',shQuote(paste0(stem,'.pdf')),shQuote(v9_path('review/pdf_rendered',V9_NAMES[[id]]))))==0);cat('Rendered',id,'\n');invisible(paste0(stem,'.pdf'))
}

# Stable analysis IDs remain separate from manuscript display labels.
V9_IDS<-c(A='A',B='B',C='C',D='D',E='E',F='F',G='G',S11L='SA',S11M='SB',S11H='SC',S11F='SD',S11K='SE',S11I='SF')
render_fig5_current_panel<-function(id){
 if(!id%in%names(V9_IDS))return(fig5_v8_current_panel(id))
 mp<-file.path(STEP8_ROOT,'8-Visualization/tables/figure_manifest.csv');m<-fread(mp,colClasses='character',na.strings=NULL)
 ix<-which(m$Figure==if(startsWith(id,'S11'))'Supplement' else 'Figure 5');ix<-ix[m$Panel[ix]==id];stopifnot(length(ix)==1)
 target<-file.path(STEP8_ROOT,m$TargetOutput[ix]);if(id%in%c('A','C','S11I'))return(invisible(target))
 src<-render_fig5_v9_panel(V9_IDS[[id]]);for(ext in c('pdf','svg','png'))stopifnot(file.copy(sub('[.]pdf$',paste0('.',ext),src),sub('[.]pdf$',paste0('.',ext),target),overwrite=TRUE))
 sources<-strsplit(m$SourceTable[ix],';',fixed=TRUE)[[1]];stopifnot(all(file.exists(file.path(STEP8_ROOT,sources))))
 m[ix,`:=`(Status='Rendered',SourceMD5=paste(paste0(sources,'=',unname(tools::md5sum(file.path(STEP8_ROOT,sources)))),collapse=';'),QAFile='',RenderedAt=format(Sys.time(),tz='UTC',usetz=TRUE))]
 tmp<-tempfile(tmpdir=dirname(mp),fileext='.csv');fwrite(m,tmp,na='');stopifnot(file.rename(tmp,mp));invisible(target)
}
