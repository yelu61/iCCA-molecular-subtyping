# Cohort-specific CNV profiles and burden: native R/Cairo, standalone PDFs.
F4CNV <- '4-Multiomic/refinements/fig4_cnv_20261030/'
F4CNV_REVIEW <- '8-Visualization/refinements/fig4_cnv_20261030_v4/'
F4CNV_TESTS <- '4-Multiomic/refinements/fig4_cnv_tests_20261030/'
f4cnv_read <- function(name)read_source(paste0(F4CNV,name,'.csv'))
f4cnv_test_read <- function(name)read_source(paste0(F4CNV_TESTS,name,'.csv'))
F4CNV_PAIRS <- '4-Multiomic/refinements/fig4_cnv_pairwise_20261030/'
f4cnv_pair_read <- function(name)read_source(paste0(F4CNV_PAIRS,name,'.csv'))
F4CNV_GENES <- '4-Multiomic/refinements/fig4_cnv_gene_annotations_20261030/'
f4cnv_gene_read <- function(name)read_source(paste0(F4CNV_GENES,name,'.csv'))
build_fig4_cnv <- function(id) {
  pal<-load_subtype_palette()$colors[c('C1','C2','C3')]
  if(id=='S10Q') {
    d<-f4cnv_read('cytoband_frequency');summary<-f4cnv_read('burden_summary')
    region<-f4cnv_test_read('cytoband_tests')
    line<-rbind(d[,.(FinalSubtype,Chromosome,X,Direction='Gain',Frequency=100*GainFrequency)],
      d[,.(FinalSubtype,Chromosome,X,Direction='Loss',Frequency=-100*LossFrequency)])
    ymax<-max(25,ceiling(max(abs(line$Frequency))/10)*10)
    cols<-c(Gain='#4C9C78',Loss='#D36A62')
    return(function() {
      grid::grid.newpage()
      fig4_v2_text('Chromosome-ordered CNV profiles',3,128,9,TRUE,just='left')
      fig4_v2_text('2022Dong cohort  |  219 tumors  |  CNV thresholds ±0.3',3,122,7.3,col='grey35',just='left')
      for(j in 1:2) {
        xx<-127+(j-1)*23
        fig4_v2_line(c(xx,xx+4),c(127.5,127.5),col=cols[j],lwd=2)
        fig4_v2_text(names(cols)[j],xx+5.5,127.5,7.2,just='left')
      }
      grobs<-list()
      for(k in seq_along(pal)) {
        st<-names(pal)[k];ss<-summary[FinalSubtype==st]
        p<-ggplot(line[FinalSubtype==st],aes(X,Frequency,group=interaction(Chromosome,Direction)))+
          annotate('rect',xmin=seq(0,20,2),xmax=seq(1,21,2),ymin=-ymax,ymax=ymax,fill='#F3F3F3',color=NA)+
          geom_vline(xintercept=0:22,color='grey85',linewidth=.2)+
          geom_hline(yintercept=0,color='grey45',linewidth=.25)+
          geom_ribbon(aes(ymin=pmin(Frequency,0),ymax=pmax(Frequency,0),fill=Direction),alpha=.28,color=NA)+
          geom_line(aes(color=Direction),linewidth=.35)+
          scale_color_manual(values=cols,guide='none')+scale_fill_manual(values=cols,guide='none')+
          scale_x_continuous(limits=c(0,22),breaks=(1:22)-.5,labels=1:22,expand=c(0,0))+
          scale_y_continuous(limits=c(-ymax,ymax),breaks=c(-ymax,-ymax/2,0,ymax/2,ymax),labels=abs,expand=c(0,0))+
          labs(x=NULL,y=NULL)+theme_step8()+
          theme(axis.text=element_text(size=6.8),axis.line=element_blank(),axis.ticks=element_line(linewidth=.2),
            plot.margin=margin(.7,.4,.5,.4,'mm'))
        if(k<3)p<-p+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
        grobs[[k]]<-ggplotGrob(p)
        fig4_v2_rect(11.5,114.5-(k-1)*28.5,18,4.4,pal[[st]])
        fig4_v2_text(sprintf('%s  n = %d',st,ss$N),11.5,114.5-(k-1)*28.5,7.2,TRUE)
      }
      widths<-do.call(grid::unit.pmax,lapply(grobs,function(g)g$widths))
      for(k in seq_along(grobs)) {
        grobs[[k]]$widths<-widths
        grobs[[k]]$vp<-grid::viewport(x=grid::unit(21,'mm'),y=grid::unit(89-(k-1)*28.5,'mm'),
          width=grid::unit(151,'mm'),height=grid::unit(26,'mm'),just=c('left','bottom'))
        grid::grid.draw(grobs[[k]])
      }
      fig4_v2_text('Mean gene-level call frequency (%)',24,117,7.2,col='grey35',just='left')
      fig4_v2_text('Omnibus differences: S10Q; pairwise contrasts: S10H',172,117,6.8,col='grey35',just='right')
      fig4_v2_text('Chromosome (cytoband order; non-physical spacing)',97,25,7.4)
      grid::grid.force()
      namesvp<-grid::grid.ls(viewports=TRUE,grobs=FALSE,print=FALSE)$name
      pv<-namesvp[grepl('^panel\\.',namesvp)][1];grid::seekViewport(pv)
      left<-grid::deviceLoc(grid::unit(0,'npc'),grid::unit(0,'npc'),valueOnly=TRUE)$x*25.4
      right<-grid::deviceLoc(grid::unit(1,'npc'),grid::unit(0,'npc'),valueOnly=TRUE)$x*25.4
      grid::upViewport(0)
      maxlog<-max(2,ceiling(max(-log10(pmax(region$FDR,1e-12)))))
      gray<-grDevices::colorRampPalette(c('white','#383838'))(101)
      for(j in 1:2) {
        direction<-c('Gain','Loss')[j];zz<-region[Direction==direction]
        value<-pmin(maxlog,-log10(pmax(zz$FDR,1e-12)))
        fill<-gray[1+pmin(100,floor(value/maxlog*100))]
        yy<-20-(j-1)*5
        grid::grid.rect(x=grid::unit(left+zz$X/22*(right-left),'mm'),y=grid::unit(yy,'mm'),
          width=grid::unit((right-left)/22/zz$BandsInChromosome+.01,'mm'),height=grid::unit(3.3,'mm'),
          gp=grid::gpar(fill=fill,col=NA))
        sig<-zz[FDR<.05]
        if(nrow(sig))grid::grid.points(grid::unit(left+sig$X/22*(right-left),'mm'),grid::unit(rep(yy, nrow(sig)),'mm'),
          pch=15,size=grid::unit(.7,'mm'),gp=grid::gpar(col=cols[[direction]]))
        fig4_v2_text(paste0(direction,' q'),left-2,yy,6.8,just='right')
      }
      fig4_v2_text(sprintf('BH q < 0.05: gain %d / 767; loss %d / 767',
        sum(region$Direction=='Gain'&region$FDR<.05),sum(region$Direction=='Loss'&region$FDR<.05)),3,5.5,6.7,col='grey30',just='left')
      for(i in 1:101)fig4_v2_rect(109+(i-.5)*.35,5.5,.36,2.3,gray[i])
      fig4_v2_text('−log10 q',105,5.5,6.7,just='right')
      fig4_v2_text('0',109,2,6.5);fig4_v2_text(as.character(maxlog),144.35,2,6.5)
      fig4_v2_text('Colored ticks: q < 0.05',172,5.5,6.5,just='right')
    })
  }
  if(id=='S10H') {
    d<-f4cnv_pair_read('pairwise_display')
    stopifnot(nrow(d)==54L)
    levelsrow<-unique(d$Cytoband)
    d[,`:=`(Cytoband=factor(Cytoband,levels=rev(levelsrow)),
      Contrast=factor(Contrast,levels=c('C2 vs C1','C3 vs C1','C3 vs C2')),
      Direction=factor(Direction,levels=c('Gain','Loss')),
      Label=sprintf('%+.1f%s',DifferencePP,ifelse(FDR<.05,'*','')))]
    d[,ContrastIndex:=as.integer(Contrast)]
    ann<-f4cnv_gene_read('region_gene_summary')
    ann[,`:=`(Cytoband=factor(Cytoband,levels=rev(levelsrow)),Direction=factor(Direction,levels=c('Gain','Loss')))]
    stopifnot(nrow(ann)==18L,sum(ann$MeasuredGeneFeatures)==633L)
    lim<-max(10,ceiling(max(abs(d$DifferencePP))/5)*5)
    return(ggplot(d,aes(ContrastIndex,Cytoband,fill=DifferencePP))+
      geom_tile(width=.94,height=.94,color='white',linewidth=.2)+
      geom_text(aes(label=Label),size=2.65,family=STEP8_FONT$family,color='grey10')+
      geom_text(data=ann,aes(x=3.7,y=Cytoband,label=GeneExamples),inherit.aes=FALSE,hjust=0,
        size=2.55,fontface='italic',family=STEP8_FONT$family,color='grey20')+
      geom_text(data=ann,aes(x=7.0,y=Cytoband,label=MeasuredGeneFeatures),inherit.aes=FALSE,
        size=2.55,family=STEP8_FONT$family,color='grey35')+
      facet_grid(Direction~.,scales='free_y',space='free_y')+
      scale_x_continuous(position='top',breaks=c(1,2,3,5.1,7),
        labels=c('C2 − C1','C3 − C1','C3 − C2','Gene examples','n genes'),limits=c(.45,7.35),expand=c(0,0))+
      scale_y_discrete(expand=expansion(add=.25))+
      scale_fill_gradient2(low='#527C9D',mid='#FCFCFA',high='#C87D53',midpoint=0,
        limits=c(-lim,lim),breaks=c(-lim,0,lim),name='Frequency difference (percentage points)')+
      labs(title='Pairwise regional CNV differences',subtitle='2022Dong | 18 omnibus-significant band/direction entries',
        x=NULL,y=NULL,caption='* BH q < 0.05 across all 4,602 pairwise tests')+
      theme_step8()+theme(plot.title=element_text(size=9,face='bold'),
        plot.subtitle=element_text(size=7,color='grey35'),axis.text=element_text(size=7.4),
        axis.ticks=element_blank(),axis.line=element_blank(),strip.text=element_text(size=7.5),
        strip.background=element_rect(fill='#F5F5F5',color=NA),strip.text.y=element_text(angle=0,size=7.3),
        panel.spacing=grid::unit(2,'mm'),legend.position='bottom',
        legend.title=element_text(size=7),legend.text=element_text(size=6.8),
        plot.caption=element_text(size=6.8,color='grey30',hjust=.5),
        plot.margin=margin(3,3,2,3,'mm'))+
      guides(fill=guide_colorbar(title.position='top',title.hjust=.5,barwidth=grid::unit(44,'mm'),barheight=grid::unit(2.8,'mm'))))
  }
  if(id=='S10A') {
    d<-f4cnv_read('burden_samples');s<-f4cnv_read('burden_summary')
    tests<-f4cnv_read('burden_tests');global<-f4cnv_read('burden_global_test')
    stopifnot(length(unique(round(tests$FDR,12)))==1L)
    top<-max(d$MeanAbsoluteCNV)
    d[,FinalSubtype:=factor(FinalSubtype,levels=names(pal))]
    ggplot(d,aes(FinalSubtype,MeanAbsoluteCNV,color=FinalSubtype))+
      geom_boxplot(aes(fill=FinalSubtype),width=.5,outlier.shape=NA,alpha=.16,linewidth=.35)+
      geom_point(position=position_jitter(width=.16,height=0,seed=20261030),size=1.05,alpha=.62,stroke=0)+
      scale_color_manual(values=pal,guide='none')+scale_fill_manual(values=pal,guide='none')+
      scale_x_discrete(labels=setNames(sprintf('%s\nn = %d',s$FinalSubtype,s$N),s$FinalSubtype))+
      scale_y_continuous(limits=c(0,top*1.08),expand=expansion(mult=c(0,.015)))+
      labs(title='CNV burden across subtypes',subtitle=sprintf('2022Dong  |  Kruskal–Wallis P = %.2f',global$PValue),
        x=NULL,y='Mean absolute gene-level CNV',
        caption=sprintf('Pairwise BH q = %.2f for all three contrasts',tests$FDR[1]))+theme_step8()+
      theme(plot.title=element_text(size=9,face='bold'),plot.subtitle=element_text(size=7.2,color='grey35'),
        axis.title=element_text(size=8),axis.text=element_text(size=7.5),
        plot.caption=element_text(size=6.8,color='grey35',hjust=.5),
        plot.margin=margin(3,3,2,3,'mm'))
  } else stop('Unknown CNV panel: ',id)
}
register_fig4_cnv <- function() {
  m<-read_source('8-Visualization/tables/figure_manifest.csv');base<-copy(m[Figure=='Supplement'&Panel=='S10A'][1])
  sources<-unique(c(f4cnv_read('input_hashes')$Path,
    paste0(F4CNV,c('burden_samples','burden_summary','burden_tests','burden_global_test','gene_frequency',
      'cytoband_frequency','cytoband_layout','gene_annotation_audit','release_gate_checks','input_hashes','output_hashes'),'.csv'),
    '8-Visualization/scripts/render_fig4_cnv.R','8-Visualization/scripts/render_fig4_genomics_v2.R',
    '4-Multiomic/scripts/prepare_fig4_cnv_statistics.R',f4cnv_test_read('input_hashes')$Path,
    paste0(F4CNV_TESTS,c('input_hashes','cytoband_tests','cytoband_sample_fractions'),'.csv'),
    '4-Multiomic/scripts/prepare_fig4_cnv_pairwise.R',f4cnv_pair_read('input_hashes')$Path,
    f4cnv_gene_read('input_hashes')$Path,
    paste0(F4CNV_GENES,c('input_hashes','all_measured_region_genes','figure_gene_examples','region_gene_summary'),'.csv'),
    paste0(F4CNV_PAIRS,c('pairwise_all','pairwise_display','region_patterns','summary','input_hashes'),'.csv'),
    '8-Visualization/scripts/fig_style.R','8-Visualization/scripts/stat_annotate.R',
    '8-Visualization/scripts/figure_config.json',SUBTYPE_PALETTE_CSV,STUDY_PALETTE_CSV))
  for(id in c('S10A','S10H','S10Q')) {
    row<-copy(base);row[,`:=`(Panel=id,Role=if(id=='S10A')'CNV burden, frozen mean absolute gene-level CNV' else if(id=='S10H')'Omnibus-significant CNV regions: pairwise contrasts and measured gene annotations' else 'Cytoband-ordered CNV profiles with sample-level omnibus FDR tracks',
      SourceFigure='',SourceTable=paste(sources,collapse=';'),
      TargetOutput=if(id=='S10A')'8-Visualization/figures/supplement/Supplement_S10A.pdf' else if(id=='S10H')'8-Visualization/figures/supplement/Supplement_S10H_cnv_region_contrasts.pdf' else '8-Visualization/figures/supplement/Supplement_S10Q_cnv_chromosome_profiles.pdf',
      Status='AnalysisReady',SubtypeScope='2022Dong Main n=219 (C1=47,C2=105,C3=67)',
      Notes='Single-cohort exploratory CNV. Frozen mean absolute burden; two-sided rank tests with BH over three pairwise comparisons. Band profiles average gene-level call frequencies; ordinal equal-band chromosome axis, not physical lengths. Thresholds +/-0.3; not FGA or GISTIC peaks. Regional KW BH1534; all4602 pairwise Wilcoxon tests before omnibus display selection. H shows all54 contrasts across18 global-significant band/direction entries, including nonsignificance. Previous driver display retired to review history. H adds measured gene examples and full633-feature NCBI reference table; annotations do not alter selection or statistics.',
      WidthMM=if(id=='S10A')'85' else '174',HeightMM=if(id=='S10A')'90' else if(id=='S10H')'130' else '132',SourceMD5='',QAFile='',RenderedAt='')]
    idx<-which(m$Figure=='Supplement'&m$Panel==id)
    if(length(idx))m[idx]<-row else m<-rbindlist(list(m,row),use.names=TRUE)
  }
  tmp<-tempfile(tmpdir=file.path(STEP8_ROOT,'8-Visualization/tables'),fileext='.csv')
  fwrite(m,tmp,na='');stopifnot(file.rename(tmp,file.path(STEP8_ROOT,'8-Visualization/tables/figure_manifest.csv')))
}
render_fig4_cnv <- function(id) {
  row<-preflight_subfigure('Supplement',id)
  stopifnot(all(f4cnv_read('release_gate_checks')$Status=='PASS'))
  save_panel(build_fig4_cnv(id),row$TargetOutput,as.numeric(row$WidthMM),as.numeric(row$HeightMM))
  src<-strsplit(row$SourceTable,';',fixed=TRUE)[[1]]
  hashes<-paste(paste0(src,'=',unname(tools::md5sum(file.path(STEP8_ROOT,src)))),collapse=';')
  m<-read_source('8-Visualization/tables/figure_manifest.csv')
  m[Figure=='Supplement'&Panel==id,`:=`(Status='Rendered',SourceMD5=hashes,QAFile=paste0(F4CNV_REVIEW,'QA.md'),RenderedAt=format(Sys.time(),tz='UTC',usetz=TRUE))]
  tmp<-tempfile(tmpdir=file.path(STEP8_ROOT,'8-Visualization/tables'),fileext='.csv')
  fwrite(m,tmp,na='');stopifnot(file.rename(tmp,file.path(STEP8_ROOT,'8-Visualization/tables/figure_manifest.csv')))
}
