# Display-only B/C revision: fixed 114-gene classifier; no model fitting.
if(!exists('fig5_v5_current_panel',inherits=FALSE)){
  # Repeated notebook source calls must not re-register Arial after it is used.
  source('8-Visualization/scripts/render_fig5_optimized_20261031.R')
  fig5_v5_current_panel <- render_fig5_current_panel
}
BC6 <- '8-Visualization/refinements/fig5_BC_20261031_v6'
bc6_path <- function(...) file.path(STEP8_ROOT, BC6, ...)
BC6_NAMES <- c(B='Figure5B_consensus_feature_evidence', C='Figure5C_subtype_gene_ranking', S11J='Supplement_S11J_RF_ranking')
BC6_DIMS <- list(B=c(174,114), C=c(174,158), S11J=c(130,105))
bc6_prepare <- function(){
  d <- opt_read('feature_evidence_114')
  stopifnot(nrow(d)==114, !anyDuplicated(d$Gene), all(d$Methods>=2))
  d[,SupportGroup:=fifelse(Methods==3,1L,fifelse(LimmaFull & RFFull,2L,fifelse(RFFull & LASSOFull,3L,4L)))]
  setorder(d,SupportGroup,-ConsensusFraction,RFRank,Gene)
  d[,DisplayOrder:=seq_len(.N)]
  fwrite(d,bc6_path('source_data/B_gene_order.csv'))
  ranks <- rbindlist(lapply(c('C1','C2','C3'), function(cl){
    z <- copy(d); others <- setdiff(c('C1','C2','C3'),cl)
    z[,Class:=cl]
    z[,ClassMean:=get(paste0(cl,'MeanTrainingZ'))]
    z[,OtherMaximum:=pmax(get(paste0(others[1],'MeanTrainingZ')),get(paste0(others[2],'MeanTrainingZ')))]
    z[,ExpressionMargin:=ClassMean-OtherMaximum]
    # Symbol resolves exact ties; stability and RF do not alter expression ranking.
    setorder(z,-ExpressionMargin,Gene);z[,WithinClassRank:=seq_len(.N)]
    z[,HigherThanBoth:=ExpressionMargin>0];z[,Labelled:=WithinClassRank<=6]
    z
  }))
  stopifnot(nrow(ranks)==342,all(is.finite(ranks$ExpressionMargin)),all(ranks[Labelled==TRUE]$HigherThanBoth))
  fwrite(ranks,bc6_path('source_data/C_subtype_ranking_all114.csv'))
  fwrite(ranks[Labelled==TRUE],bc6_path('source_data/C_top6_evidence.csv'))
  invisible(list(panel=d,ranks=ranks))
}

bc6_text <- function(s,x,y,size=7,bold=FALSE,just='left',col=FIG5_INK){
  grid.text(s,unit(x,'native'),unit(y,'native'),just=just,
    gp=gpar(fontfamily='Arial',fontsize=size,fontface=if(bold)'bold' else 'plain',col=col,lineheight=1.12))
}
bc6_rect <- function(x,y,w,h,fill,col=NA){
  grid.rect(unit(x,'native'),unit(y,'native'),unit(w,'native'),unit(h,'native'),gp=gpar(fill=fill,col=col,lwd=.4))
}
bc6_line <- function(x,y,col='grey85',lwd=.5){grid.lines(unit(x,'native'),unit(y,'native'),gp=gpar(col=col,lwd=lwd))}

bc6_B <- function(){function(){
  d<-fread(bc6_path('source_data/B_gene_order.csv'))
  grid.newpage();pushViewport(viewport(xscale=c(0,174),yscale=c(0,114)))
  bc6_text('How the 114 consensus genes were selected',4,109,10,TRUE)
  bc6_text('Final-fit method support, repeated selection and reference-class expression',4,103,7.5)
  x0<-23;x1<-170;gap<-2;cw<-(x1-x0-3*gap)/114
  d[,X:=x0+(DisplayOrder-.5)*cw+(SupportGroup-1)*gap]
  groups<-d[,.(N=.N,Left=min(X)-cw/2,Right=max(X)+cw/2),by=SupportGroup]
  labels<-c('All three','limma + RF','RF + LASSO','limma\n+ LASSO')
  for(j in 1:4){g<-groups[SupportGroup==j];mid<-(g$Left+g$Right)/2
    bc6_rect(mid,92,g$Right-g$Left,10,if(j==1)'#DDE8EF' else '#EEF2F5')
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
    for(j in 1:4){g<-groups[SupportGroup==j];bc6_line(c(g$Left,g$Right),c(yy,yy),'grey90')}
    bc6_text(sprintf('%d%%',100*v),x0-2,yy,6.6,FALSE,'right')
  }
  for(i in seq_len(nrow(d))){yy<-39+17*d$ConsensusFraction[i];bc6_line(rep(d$X[i],2),c(39,yy),'#A7BCCB',.6)
    grid.points(unit(d$X[i],'native'),unit(yy,'native'),pch=16,size=unit(.65,'mm'),gp=gpar(col=FIG5_INK))}
  bc6_text('Mean training-standardized expression by reference class',x0,33.5,7.6,TRUE)
  ramp<-circlize::colorRamp2(c(-1.5,0,1.5),STEP8_CONFIG$effect_colors)
  for(j in 1:3){cl<-c('C1','C2','C3')[j];yy<-28-(j-1)*3.8
    bc6_text(cl,x0-2,yy,7,TRUE,'right',FIG5_SUBTYPE[cl])
    for(i in seq_len(nrow(d)))bc6_rect(d$X[i],yy,cw-.07,3.3,ramp(d[[paste0(cl,'MeanTrainingZ')]][i]))}
  bc6_rect(24,14,1.8,1.8,FIG5_INK);bc6_text('Selected',27,14,6.5)
  bc6_rect(47,14,1.8,1.8,'#E7EDF1');bc6_text('Not selected',50,14,6.5)
  bc6_text('Mean z',105,14,6.5)
  for(i in 1:30)bc6_rect(123+(i-.5)*.7,14,.7,2,ramp(-1.5+(i-.5)*3/30))
  bc6_text('-1.5',121,14,6.5,FALSE,'right');bc6_text('+1.5',146,14,6.5)
  bc6_text('One column per gene, aligned across tracks. Within each support group: selection frequency high to low; RF rank breaks ties.',4,7,6.5)
  bc6_text('Expression and selection describe the development resource; they do not establish independent marker validity.',4,3,6.5)
  popViewport()
}}

bc6_C <- function(){
  d<-fread(bc6_path('source_data/C_subtype_ranking_all114.csv'))
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
      geom_point(data=z[HigherThanBoth==TRUE],aes(size=ConsensusFraction,shape=Support),fill=FIG5_SUBTYPE[cl],colour=FIG5_INK,stroke=.3)+
      geom_segment(data=top,aes(xend=57,yend=LabelY),colour='#677D8E',linewidth=.25)+
      geom_text(data=top,aes(x=59,y=LabelY,label=Gene),size=2.45,fontface='italic',colour=FIG5_INK,hjust=0)+
      scale_shape_manual(values=c(21,24))+scale_size_continuous(range=c(1.2,2.7),limits=c(0,1),breaks=c(.25,.5,1),labels=c('25%','50%','100%'))+
      scale_x_continuous(limits=c(1,120),breaks=c(1,40,80,114),expand=expansion(mult=c(.02,.02)))+
      coord_cartesian(ylim=c(ymin,ymax))+
      labs(title=cl,subtitle=paste0(counts[Class==cl,Positive],' / 114 higher than both other classes'),
        x='Within-class gene rank',y=if(cl=='C1')'Expression margin (standardized units)' else NULL,size='Consensus selection',shape='Final-fit support')+
      fig5_theme()+theme(plot.title=element_text(colour=FIG5_SUBTYPE[cl],size=9),plot.subtitle=element_text(size=6.3),
        axis.title=element_text(size=7),axis.text=element_text(size=6.5),legend.position='none',legend.text=element_text(size=6.5),legend.title=element_text(size=6.5),
        plot.margin=margin(2,2,2,2,'mm'))
  })
  tab<-function(){
    grid.newpage();pushViewport(viewport(xscale=c(0,174),yscale=c(0,51)))
    for(j in 1:3){cl<-c('C1','C2','C3')[j];z<-d[Class==cl & Labelled==TRUE];x0<-3+(j-1)*57.5
      bc6_text(paste0(cl,': evidence for the six highest-ranked genes'),x0,47,6.7,TRUE)
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
  leg<-function(){
    grid.newpage();pushViewport(viewport(xscale=c(0,174),yscale=c(0,9)))
    bc6_text('Consensus selection',10,4.5,6.7)
    for(k in 1:3){xx<-48+(k-1)*14
      grid.points(unit(xx,'native'),unit(4.5,'native'),pch=21,size=unit(1.2+sqrt(c(.25,.5,1)[k])*1.5,'mm'),gp=gpar(fill='#B5C5D1',col=FIG5_INK))
      bc6_text(c('25%','50%','100%')[k],xx+3,4.5,6.7)
    }
    bc6_text('Final-fit support',96,4.5,6.7)
    for(k in 1:2){xx<-125+(k-1)*23
      grid.points(unit(xx,'native'),unit(4.5,'native'),pch=c(21,24)[k],size=unit(2.3,'mm'),gp=gpar(fill='white',col=FIG5_INK))
      bc6_text(c('2 methods','3 methods')[k],xx+3,4.5,6.7)
    };popViewport()
  }
  top<-wrap_plots(plots,nrow=1)
  (top / wrap_elements(full=grid.grabExpr(leg(),width=174/25.4,height=9/25.4)) /
     wrap_elements(full=grid.grabExpr(tab(),width=174/25.4,height=51/25.4)) + plot_layout(heights=c(1.6,.17,1)))+
    plot_annotation(title='Subtype-representative expression within the 114-gene classifier',
      subtitle='All 114 genes ranked separately for C1, C2 and C3; labels show the six highest expression margins',
      caption='Margin = class mean z minus the larger mean z of the other two classes; positive values indicate higher expression than both.\nSel. %: consensus selection in 25 outer fits. L / R / A: limma / RF / LASSO support (dark: selected; pale: not selected).\nRF rank is global across 10,531 genes, not subtype-specific. Display ranking only; no new panel selection or independent validation.',
      theme=fig5_theme()+theme(plot.title=element_text(size=9.5),plot.subtitle=element_text(size=7),plot.caption=element_text(size=6.4)))
}

bc6_RF <- function(){
  d<-opt_read('feature_evidence_114');all<-opt_read('feature_evidence_all_genes')
  stopifnot(all(all$RFRank>0));d[,Support:=factor(Methods,levels=c(2,3),labels=c('2 methods','3 methods'))]
  labels<-rbind(head(d[order(RFRank)],6),tail(d[order(RFRank)],2))
  ggplot(all,aes(RFRank,RFImportance))+geom_line(colour='#C9D1D7',linewidth=.4)+
    geom_point(data=d,aes(fill=ConsensusFraction,shape=Support),size=2.2,colour=FIG5_INK,stroke=.3)+
    geom_text_repel(data=labels,aes(label=Gene),size=2.55,seed=20261031,nudge_y=.25,min.segment.length=0,segment.size=.25,box.padding=.5,max.overlaps=Inf)+
    scale_x_log10(breaks=c(1,10,100,1000,10000),labels=c('1','10','100','1,000','10,000'))+
    scale_shape_manual(values=c(21,24))+scale_fill_gradient(low='#E7EEF3',high='#3F718E',limits=c(0,1),labels=scales::label_percent())+
    labs(title='RF importance in the full screening universe',subtitle='114 consensus-panel genes highlighted among 10,531 genes',
      x='Global RF importance rank (log scale)',y='Scaled permutation importance',fill='Consensus selection',shape='Final-fit support',
      caption='Ranking uses RF only; color and shape overlay consensus evidence.\nFull-development scaled MeanDecreaseAccuracy; not SVM contribution or weighted consensus ranking.')+
    fig5_theme()+theme(legend.position='bottom',legend.box='vertical',legend.text=element_text(size=6.5),legend.title=element_text(size=6.5),plot.caption=element_text(size=6.2),plot.title=element_text(size=9),plot.subtitle=element_text(size=7))+
    guides(fill=guide_colourbar(display='rectangles',barwidth=unit(26,'mm'),barheight=unit(2,'mm')),shape=guide_legend(override.aes=list(fill='white')))
}

render_fig5_bc6_panel <- function(id){
  stopifnot(id %in% names(BC6_NAMES));fig5_check();bc6_prepare()
  p<-switch(id,B=bc6_B(),C=bc6_C(),S11J=bc6_RF());sz<-BC6_DIMS[[id]]
  stem<-bc6_path('figures',BC6_NAMES[[id]])
  draw<-function(){if(is.function(p))p() else print(p)}
  cairo_pdf(paste0(stem,'.pdf'),width=sz[1]/25.4,height=sz[2]/25.4,family='Arial');draw();dev.off()
  svglite::svglite(paste0(stem,'.svg'),width=sz[1]/25.4,height=sz[2]/25.4);draw();dev.off()
  ragg::agg_png(paste0(stem,'.png'),width=sz[1],height=sz[2],units='mm',res=300,background='white');draw();dev.off()
  stopifnot(system2('pdftoppm',c('-png','-singlefile','-r','150',shQuote(paste0(stem,'.pdf')),shQuote(bc6_path('review/pdf_rendered',BC6_NAMES[[id]]))))==0)
  cat('Rendered',id,'at',paste(sz,collapse=' x '),'mm\n')
  invisible(paste0(stem,'.pdf'))
}

# Updated notebook entry keeps A/D/E/F and existing supplements on their prior routes.
render_fig5_current_panel <- function(id){
  if(!id %in% names(BC6_NAMES))return(fig5_v5_current_panel(id))
  src<-render_fig5_bc6_panel(id)
  mpath<-file.path(STEP8_ROOT,'8-Visualization/tables/figure_manifest.csv')
  m<-fread(mpath,colClasses='character',na.strings=NULL)
  ix<-which(m$Figure==if(startsWith(id,'S11'))'Supplement' else 'Figure 5');ix<-ix[m$Panel[ix]==id]
  stopifnot(length(ix)==1);target<-file.path(STEP8_ROOT,m$TargetOutput[ix])
  for(ext in c('pdf','svg'))stopifnot(file.copy(sub('[.]pdf$',paste0('.',ext),src),sub('[.]pdf$',paste0('.',ext),target),overwrite=TRUE))
  sources<-strsplit(m$SourceTable[ix],';',fixed=TRUE)[[1]]
  m[ix,`:=`(Status='Rendered',SourceMD5=paste(paste0(sources,'=',unname(tools::md5sum(file.path(STEP8_ROOT,sources)))),collapse=';'),QAFile='',RenderedAt=format(Sys.time(),tz='UTC',usetz=TRUE))]
  tmp<-tempfile(tmpdir=dirname(mpath),fileext='.csv');fwrite(m,tmp,na='');stopifnot(file.rename(tmp,mpath));invisible(target)
}
