# Figure5 visual revision v7; original B retained, new class-grouped alternative and compact C.
if(!exists('fig5_v6_current_panel',inherits=FALSE)){
 source('8-Visualization/scripts/render_fig5_BC_20261031_v6.R')
 fig5_v6_current_panel<-render_fig5_current_panel
}
V7<-'8-Visualization/refinements/fig5_BCD_20261031_v7'
v7_path<-function(...)file.path(STEP8_ROOT,V7,...)
RFE_RUN<-'5-ClassifierValidation/runs/bulk-rnaseq__20261031_step5_svm_rfe_v1'
v7_prepare<-function(){
 d<-fread(file.path(STEP8_ROOT,BC6,'source_data/B_gene_order.csv'))
 means<-as.matrix(d[,c('C1MeanTrainingZ','C2MeanTrainingZ','C3MeanTrainingZ'),with=FALSE])
 d[,SubtypeGroup:=max.col(means,ties.method='first')]
 d[,ExpressionMargin:=apply(means,1,function(v){v<-sort(v,decreasing=TRUE);v[1]-v[2]})]
 setorder(d,SubtypeGroup,-ExpressionMargin,-ConsensusFraction,Gene);d[,DisplayOrder:=seq_len(.N)]
 fwrite(d,v7_path('source_data/B_subtype_gene_order.csv'))
}

v7_D <- function(){
 rd<-function(n)fread(file.path(STEP8_ROOT,RFE_RUN,'tables',paste0(n,'.csv')))
 weights<-rd('full_model_weight_ranking');z<-weights[Rank<=15]
 z[,Gene:=factor(Gene,levels=rev(Gene))]
 left<-ggplot(z,aes(RelativeWeight*100,Gene))+geom_col(fill=FIG5_INK,width=.65)+
  scale_y_discrete(labels=function(x)ifelse(x=='C2','C2 (gene)',x))+
  labs(title='Full-model gene weights',subtitle='15 highest of 114; three SVM boundaries',x='Relative squared weight (%)',y=NULL)+fig5_theme()+
  theme(plot.title=element_text(size=8),plot.subtitle=element_text(size=6.4),axis.text=element_text(size=6.8),axis.text.y=element_text(face='italic'),axis.title=element_text(size=7))
 d<-rd('marker_size_performance')[Metric%in%c('BalancedAccuracy','MacroF1','MacroAUROC')]
 labs<-c(BalancedAccuracy='Balanced accuracy',MacroF1='Macro-F1',MacroAUROC='Macro-AUROC')
 d[,MetricLabel:=factor(labs[Metric],levels=unname(labs))]
 ref<-fig5_read('nested_performance_ci')[Model==FIG5_PRIMARY & Metric%in%c('BalancedAccuracy','MacroF1','C1_AUROC','C2_AUROC','C3_AUROC')]
 ref<-rbind(ref[Metric%in%c('BalancedAccuracy','MacroF1'),.(Metric,Estimate)],data.table(Metric='MacroAUROC',Estimate=mean(ref[grepl('AUROC',Metric),Estimate])))
 ref[,MetricLabel:=factor(labs[Metric],levels=unname(labs))]
 colours<-setNames(c('#315B78','#AE6157','#6A9178'),unname(labs))
 curve<-ggplot(d,aes(K,Estimate,colour=MetricLabel,fill=MetricLabel))+
  geom_hline(data=ref,aes(yintercept=Estimate,colour=MetricLabel),linetype=2,linewidth=.35,alpha=.65)+
  geom_ribbon(aes(ymin=Lower95,ymax=Upper95),alpha=.1,colour=NA)+geom_line(linewidth=.55)+
  geom_point(aes(shape=MetricLabel),size=1.7)+scale_colour_manual(values=colours)+scale_fill_manual(values=colours)+
  scale_shape_manual(values=c(16,17,15))+scale_x_continuous(breaks=c(2,10,20,40,60,80))+
  scale_y_continuous(breaks=seq(0,1,.02),labels=scales::label_percent(accuracy=1))+
  labs(title='Recursive feature elimination',subtitle='Training-only ranking; repeated nested CV',x='Retained genes (K)',y='OOF performance',colour=NULL,fill=NULL,shape=NULL)+fig5_theme()+
  theme(plot.title=element_text(size=8),plot.subtitle=element_text(size=6.4),axis.title=element_text(size=7),axis.text=element_text(size=6.7),legend.position='bottom',legend.text=element_text(size=6.3),legend.key.width=unit(4,'mm'))+
  guides(colour=guide_legend(nrow=1),fill='none',shape=guide_legend(nrow=1))
 delta<-rd('paired_BA_difference')
 loss<-ggplot(delta,aes(K,Estimate*100))+geom_hline(yintercept=0,linetype=2,colour='grey60',linewidth=.3)+
  geom_errorbar(aes(ymin=Lower95*100,ymax=Upper95*100),width=1.3,linewidth=.35)+geom_line(colour=FIG5_INK,linewidth=.4)+geom_point(colour=FIG5_INK,size=1.6)+
  scale_x_continuous(breaks=c(2,10,20,40,60,80))+labs(x='Retained genes (K)',y='Δ BA (percentage points)',title='Difference from the full workflow')+
  fig5_theme()+theme(plot.title=element_text(size=7.5),axis.title=element_text(size=6.7),axis.text=element_text(size=6.4))
 (left | (curve/loss+plot_layout(heights=c(1.7,1))))+plot_layout(widths=c(.82,1.5))+
  plot_annotation(title='SVM weights and the cost of feature reduction',
    caption='Left: frozen 114-gene model. Right: fold-specific consensus → SVM-RFE → held-out predictions.\nDashed references: full selection workflow. Intervals: conditional 95% patient bootstrap. No optimal K selected.',
    theme=fig5_theme()+theme(plot.caption=element_text(size=6.3)))
}

V7_NAMES<-c(B_subtype='Figure5B_subtype_order_alternative',C='Figure5C_subtype_gene_ranking',D='Figure5D_SVM_RFE')
V7_DIMS<-list(B_subtype=c(174,114),C=c(174,122),D=c(174,125))
render_fig5_v7_panel<-function(id){
 stopifnot(id%in%names(V7_NAMES));fig5_check();v7_prepare()
 p<-switch(id,B_subtype=v7_Bsubtype(),C=v7_C(),D=v7_D());sz<-V7_DIMS[[id]];stem<-v7_path('figures',V7_NAMES[[id]])
 draw<-function(){if(is.function(p))p() else print(p)}
 cairo_pdf(paste0(stem,'.pdf'),width=sz[1]/25.4,height=sz[2]/25.4,family='Arial');draw();dev.off()
 svglite::svglite(paste0(stem,'.svg'),width=sz[1]/25.4,height=sz[2]/25.4);draw();dev.off()
 ragg::agg_png(paste0(stem,'.png'),width=sz[1],height=sz[2],units='mm',res=300,background='white');draw();dev.off()
 stopifnot(system2('pdftoppm',c('-png','-singlefile','-r','150',shQuote(paste0(stem,'.pdf')),shQuote(v7_path('review/pdf_rendered',V7_NAMES[[id]]))))==0)
 cat('Rendered',id,'\n');invisible(paste0(stem,'.pdf'))
}

render_fig5_current_panel<-function(id){
 mpath<-file.path(STEP8_ROOT,'8-Visualization/tables/figure_manifest.csv');m<-fread(mpath,colClasses='character',na.strings=NULL)
 ix<-which(m$Figure==if(startsWith(id,'S11'))'Supplement' else 'Figure 5');ix<-ix[m$Panel[ix]==id];stopifnot(length(ix)==1)
 target<-file.path(STEP8_ROOT,m$TargetOutput[ix])
 # Selected B and former D/OOF are frozen; never rerender with obsolete panel mappings.
 if(id%in%c('A','B','S11K')){stopifnot(file.exists(target));return(invisible(target))}
 if(!id%in%c('C','D'))return(fig5_v6_current_panel(id))
 src<-render_fig5_v7_panel(id)
 for(ext in c('pdf','svg'))stopifnot(file.copy(sub('[.]pdf$',paste0('.',ext),src),sub('[.]pdf$',paste0('.',ext),target),overwrite=TRUE))
 sources<-strsplit(m$SourceTable[ix],';',fixed=TRUE)[[1]]
 m[ix,`:=`(Status='Rendered',SourceMD5=paste(paste0(sources,'=',unname(tools::md5sum(file.path(STEP8_ROOT,sources)))),collapse=';'),QAFile='',RenderedAt=format(Sys.time(),tz='UTC',usetz=TRUE))]
 tmp<-tempfile(tmpdir=dirname(mpath),fileext='.csv');fwrite(m,tmp,na='');stopifnot(file.rename(tmp,mpath));invisible(target)
}
v7_Bsubtype <- function(){function(){
  d<-fread(v7_path('source_data/B_subtype_gene_order.csv'))
  grid.newpage();pushViewport(viewport(xscale=c(0,174),yscale=c(0,114)))
  bc6_text('Consensus features grouped by reference subtype',4,109,10,TRUE)
  bc6_text('All 114 genes; subtype blocks ordered by expression margin',4,103,7.5)
  x0<-23;x1<-170;gap<-2;cw<-(x1-x0-2*gap)/114
  d[,X:=x0+(DisplayOrder-.5)*cw+(SubtypeGroup-1)*gap]
  groups<-d[,.(N=.N,Left=min(X)-cw/2,Right=max(X)+cw/2),by=SubtypeGroup]
  labels<-c('C1','C2','C3')
  for(j in 1:3){g<-groups[SubtypeGroup==j];mid<-(g$Left+g$Right)/2
    bc6_rect(mid,92,g$Right-g$Left,10,grDevices::adjustcolor(FIG5_SUBTYPE[j],alpha.f=.15))
    bc6_text(labels[j],mid,if(j==4)94 else 94,if(j==4)6 else 7.2,TRUE,'centre')
    bc6_text(paste0('n = ',g$N),mid,88.8,6.4,FALSE,'centre')
  }
  bc6_text('Method support in the final development fit',x0,83,7.6,TRUE)
  for(j in 1:3){yy<-77-(j-1)*4
    bc6_text(c('limma','RF','LASSO')[j],x0-2,yy,7,FALSE,'right')
    for(i in seq_len(nrow(d)))bc6_rect(d$X[i],yy,cw-.07,3.3,if(d[[c('LimmaFull','RFFull','LASSOFull')[j]]][i])FIG5_INK else '#E7EDF1')
  }
  bc6_text('Consensus selection frequency across 25 outer training fits',x0,62,7.6,TRUE)
  for(v in c(0,.5,1)){yy<-39+17*v
    for(j in 1:3){g<-groups[SubtypeGroup==j];bc6_line(c(g$Left,g$Right),c(yy,yy),'grey90')}
    bc6_text(sprintf('%d%%',100*v),x0-2,yy,6.6,FALSE,'right')
  }
  for(i in seq_len(nrow(d))){yy<-39+17*d$ConsensusFraction[i];bc6_line(rep(d$X[i],2),c(39,yy),'#A7BCCB',.6)
    grid.points(unit(d$X[i],'native'),unit(yy,'native'),pch=16,size=unit(.65,'mm'),gp=gpar(col=FIG5_INK))}
  bc6_text('Mean training-standardized expression',x0,33.5,7.6,TRUE)
  ramp<-circlize::colorRamp2(c(-1.5,0,1.5),STEP8_CONFIG$effect_colors)
  for(j in 1:3){cl<-c('C1','C2','C3')[j];yy<-28-(j-1)*3.8
    bc6_text(cl,x0-2,yy,7,TRUE,'right',FIG5_SUBTYPE[cl])
    for(i in seq_len(nrow(d)))bc6_rect(d$X[i],yy,cw-.07,3.3,ramp(d[[paste0(cl,'MeanTrainingZ')]][i]))}
  bc6_rect(24,14,1.8,1.8,FIG5_INK);bc6_text('Selected',27,14,6.5)
  bc6_rect(47,14,1.8,1.8,'#E7EDF1');bc6_text('Not selected',50,14,6.5)
  bc6_text('Mean z',105,14,6.5)
  for(i in 1:30)bc6_rect(123+(i-.5)*.7,14,.7,2,ramp(-1.5+(i-.5)*3/30))
  bc6_text('-1.5',121,14,6.5,FALSE,'right');bc6_text('+1.5',146,14,6.5)
  bc6_text('One column per gene, aligned across tracks. Within each subtype: expression margin high to low; selection frequency breaks ties.',4,7,6.5)
  bc6_text('Expression and selection describe the development resource; they do not establish independent marker validity.',4,3,6.5)
  popViewport()
}}

v7_C <- function(){
  d<-fread(file.path(STEP8_ROOT,BC6,'source_data/C_subtype_ranking_all114.csv'))
  d[,Support:=factor(Methods,levels=c(2,3),labels=c('2 methods','3 methods'))]
  counts<-d[,.(Positive=sum(HigherThanBoth)),by=Class]
  strip<-setNames(sprintf('%s  |  %d / 114 higher than both other classes',counts$Class,counts$Positive),counts$Class)
  ymax<-max(d$ExpressionMargin)+.5;ymin<-min(d$ExpressionMargin)-.15
  plots<-lapply(c('C1','C2','C3'),function(cl){
    z<-d[Class==cl];top<-z[Labelled==TRUE]
    top[,LabelY:=seq(ymax-.13,ymax-1.95,length.out=.N)]
    ggplot(z,aes(WithinClassRank,ExpressionMargin))+
      geom_hline(yintercept=0,colour='#A8B3BD',linetype=2,linewidth=.3)+
      geom_line(colour='#CED5DA',linewidth=.4)+
      geom_point(colour='#CDD4DA',size=.9)+
      geom_point(data=z[HigherThanBoth==TRUE],colour=FIG5_SUBTYPE[cl],size=1.4)+
      geom_segment(data=top,aes(xend=57,yend=LabelY),colour='#677D8E',linewidth=.25)+
      geom_text(data=top,aes(x=59,y=LabelY,label=Gene),size=2.45,fontface='italic',colour=FIG5_INK,hjust=0)+
      scale_x_continuous(limits=c(1,120),breaks=c(1,40,80,114),expand=expansion(mult=c(.02,.02)))+
      coord_cartesian(ylim=c(ymin,ymax))+
      labs(title=cl,subtitle=paste0(counts[Class==cl,Positive],' genes with Δz > 0'),
        x='Gene rank',y=if(cl=='C1')'Expression margin (Δz)' else NULL)+
      fig5_theme()+theme(plot.title=element_text(colour=FIG5_SUBTYPE[cl],size=9),plot.subtitle=element_text(size=6.3),
        axis.title=element_text(size=7),axis.text=element_text(size=6.5),legend.position='none',legend.text=element_text(size=6.5),legend.title=element_text(size=6.5),
        plot.margin=margin(2,2,2,2,'mm'))
  })
  tab<-function(){
    grid.newpage();pushViewport(viewport(xscale=c(0,174),yscale=c(0,51)))
    for(j in 1:3){cl<-c('C1','C2','C3')[j];z<-d[Class==cl & Labelled==TRUE];x0<-3+(j-1)*57.5
      # Six annotations retained; repeated section headings removed.
      xs<-x0+c(0,17,25,32,36,40,49)
      heads<-c('Gene','Margin','Sel. %','L','R','A','RF rank')
      for(k in seq_along(xs))bc6_text(heads[k],xs[k],41.5,6.2,TRUE,if(k==1)'left' else 'centre')
      bc6_line(c(x0,x0+54),c(38.5,38.5),'#BBC8D1')
      for(i in 1:6){yy<-35-(i-1)*5.3
        if(i%%2)bc6_rect(x0+27,yy,54,5.1,'#F2F5F7')
        bc6_text(z$Gene[i],xs[1],yy,6.4)
        bc6_text(sprintf('%.2f',z$ExpressionMargin[i]),xs[2],yy,6.4,FALSE,'centre')
        bc6_text(sprintf('%g',100*z$ConsensusFraction[i]),xs[3],yy,6.4,FALSE,'centre')
        for(k in 1:3)bc6_rect(xs[3+k],yy,1.7,1.7,if(z[[c('LimmaFull','RFFull','LASSOFull')[k]]][i])FIG5_INK else '#DCE4EA')
        bc6_text(z$RFRank[i],xs[7],yy,6.4,FALSE,'centre')
      }
    }
    popViewport()
  }
  top<-wrap_plots(plots,nrow=1)
  (top / wrap_elements(full=grid.grabExpr(tab(),width=174/25.4,height=51/25.4)) + plot_layout(heights=c(1.6,.85)))+
    plot_annotation(title='Subtype-representative genes',
      caption='All 114 genes. Δz: class mean − max(other class means). Sel.: selection frequency. L / R / A: limma / RF / LASSO.',
      theme=fig5_theme()+theme(plot.title=element_text(size=9.5),plot.caption=element_text(size=6.2)))
}
