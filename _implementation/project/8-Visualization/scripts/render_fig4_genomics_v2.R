# Main genomics: full driver landscape + aligned frequencies/CMH estimates.
# R grid/Cairo PDF via save_panel; Arial 7-9 pt, 174 mm wide. No panel letters.
FIG4_V2 <- '8-Visualization/refinements/fig4_genomics_20261030_v2/'
fig4_v2_text <- function(s,x,y,size=7,bold=FALSE,col='black',just='centre',italic=FALSE) {
  grid::grid.text(s,grid::unit(x,'mm'),grid::unit(y,'mm'),just=just,
    gp=grid::gpar(fontfamily=STEP8_FONT$family,fontsize=size,
      fontface=if(italic)'italic' else if(bold)'bold' else 'plain',col=col))
}
fig4_v2_rect <- function(x,y,w,h,fill,col=NA,lwd=.3) grid::grid.rect(
  grid::unit(x,'mm'),grid::unit(y,'mm'),grid::unit(w,'mm'),grid::unit(h,'mm'),
  gp=grid::gpar(fill=fill,col=col,lwd=lwd))
fig4_v2_line <- function(x,y,col='grey60',lwd=.5,lty=1) grid::grid.lines(
  grid::unit(x,'mm'),grid::unit(y,'mm'),gp=grid::gpar(col=col,lwd=lwd,lty=lty))

build_fig4_v2 <- function(id) {
  genes <- fig4_read('gene_order.csv')[order(DisplayOrder),Gene]
  pal <- load_subtype_palette()$colors[c('C1','C2','C3')]
  if(id=='A') {
    d <- fig4_read('mutation_display_long.csv');s <- fig4_read('sample_order.csv')[order(DisplayOrder)]
    s[,Sample:=as.character(Sample)]
    return(function() {
      grid::grid.newpage();fig4_v2_text('Driver mutation landscape',3,98,9,TRUE,just='left')
      left<-20;right<-111;gap<-1.5;usable<-right-left-2*gap;top<-83;bottom<-30
      rh<-(top-bottom)/length(genes)
      for(st in names(pal)) {
        ss<-s[FinalSubtype==st];width<-nrow(ss)/nrow(s)*usable;cw<-width/nrow(ss)
        fig4_v2_rect(left+width/2,(top+bottom)/2,width,top-bottom,fig4_variant_colors[['No recorded variant']])
        z<-d[FinalSubtype==st&Mutated==TRUE]
        grid::grid.rect(grid::unit(left+(match(as.character(z$Sample),ss$Sample)-.5)*cw,'mm'),
          grid::unit(top-(match(z$Gene,genes)-.5)*rh,'mm'),grid::unit(cw+.008,'mm'),grid::unit(rh,'mm'),
          gp=grid::gpar(fill=fig4_variant_colors[z$MutationClass],col=NA))
        fig4_v2_rect(left+width/2,88,width,5,pal[[st]])
        fig4_v2_text(paste0(st,'  n = ',nrow(ss)),left+width/2,88,8,TRUE,col='black')
        grid::grid.rect(grid::unit(left+(seq_len(nrow(ss))-.5)*cw,'mm'),grid::unit(26,'mm'),
          grid::unit(cw+.008,'mm'),grid::unit(2.8,'mm'),
          gp=grid::gpar(fill=load_study_palette()$colors[ss$Study],col=NA))
        left<-left+width+gap
      }
      # Overall and subtype frequencies share the frozen roster and a 0-40% scale.
      freq<-fig4_read('mutation_frequency.csv');centers<-c(121,135,149,163)
      allfreq<-freq[,.(N=sum(N),Altered=sum(Altered)),by=Gene]
      allfreq[,Fraction:=Altered/N]
      fig4_v2_text('Mutation frequency (%)',142,94,7,TRUE)
      for(k in seq_along(centers))fig4_v2_text(c('All',names(pal))[k],centers[k],88,7,TRUE)
      for(j in seq_along(genes)) {
        yy<-top-(j-.5)*rh
        fig4_v2_text(genes[j],18,yy,8,just='right',italic=TRUE)
        if(j%%2==0)fig4_v2_rect(142,yy,56,rh,'#F5F5F5')
        for(k in seq_along(centers)) {
          if(k==1L) {
            row<-allfreq[Gene==genes[j]];ss<-d[Gene==genes[j]];color<-'grey50'
          } else {
            st<-names(pal)[k-1L];row<-freq[Gene==genes[j]&FinalSubtype==st]
            ss<-d[Gene==genes[j]&FinalSubtype==st];color<-pal[[st]]
          }
          stopifnot(nrow(row)==1L,row$N==nrow(ss),row$Altered==sum(ss$Mutated))
          bw<-14*row$Fraction/.4
          if(bw>0)fig4_v2_rect(centers[k]-7+bw/2,yy,bw,3.4,grDevices::adjustcolor(color,alpha.f=.24))
          fig4_v2_text(sprintf('%.1f',100*row$Fraction),centers[k],yy,7)
        }
      }
      fig4_v2_text('Study',18,26,7,just='right')
      # Separate mutation and Study legend blocks, with aligned swatches/text.
      legend_order<-c('Missense','Nonsense','Splice','Frameshift deletion','Frameshift insertion',
        'Multi-hit','In-frame deletion','In-frame insertion','No recorded variant','Other non-silent')
      types<-legend_order[legend_order%in%unique(d$MutationClass)]
      label<-c('Frameshift deletion'='Frameshift del.','Frameshift insertion'='Frameshift ins.',
        'In-frame deletion'='In-frame del.','In-frame insertion'='In-frame ins.')
      fig4_v2_text('Mutation class',4,22,7,TRUE,just='left')
      fig4_v2_text('Study',130,22,7,TRUE,just='left')
      fig4_v2_line(c(121,121),c(4,23),col='grey80',lwd=.4)
      for(i in seq_along(types)) {
        xx<-5+((i-1)%%3)*38; yy<-17-floor((i-1)/3)*5
        fig4_v2_rect(xx,yy,2.4,2.4,fig4_variant_colors[[types[i]]])
        lab<-if(types[i]%in%names(label))label[[types[i]]] else types[i]
        fig4_v2_text(lab,xx+2,yy,7,just='left')
      }
      studies<-unique(s$Study);cols<-load_study_palette()$colors
      for(i in seq_along(studies)) {
        yy<-17-(i-1)*5;fig4_v2_rect(131,yy,2.4,2.4,cols[[studies[i]]]);fig4_v2_text(studies[i],133,yy,7,just='left')
      }
    })
  }
  if(id %in% c('B','S10M')) {
    d<-fig4_read('mutation_frequency.csv')
    if(id=='S10M') {
      f<-read_source('4-Multiomic/refinements/fig4_frequency_20261030/dong_tcga_frequency.csv')
      t<-read_source('4-Multiomic/refinements/fig4_frequency_20261030/dong_tcga_frequency_associations.csv')
      d<-merge(f,t,by.x=c('Gene','FinalSubtype'),by.y=c('Gene','Subtype'))
    }
    return(function() {
      grid::grid.newpage()
      if(id=='S10M')fig4_v2_text('Driver associations: Dong + TCGA sensitivity',3,124,9,TRUE,just='left')
      if(id=='S10M')fig4_v2_text('Mutation frequency (%)',45,114,8,TRUE)
      fig4_v2_text('Study-stratified odds ratio (95% CI)',if(id=='B')98 else 124,114,8,TRUE)
      top<-96;bottom<-27;rh<-(top-bottom)/length(genes)
      centers<-c(28,44,60);starts<-if(id=='B')c(30,78,126) else c(80,111,142)
      fw<-if(id=='B')40 else 24;lims<-c(.02,200)
      trans<-function(v,k)starts[k]+(log10(v)-log10(lims[1]))/diff(log10(lims))*fw
      for(j in seq_along(genes)) {
        yy<-top-(j-.5)*rh
        if(j%%2==0)fig4_v2_rect(87,yy,166,rh,'#F5F5F5')
        fig4_v2_text(genes[j],18,yy,8,just='right',italic=TRUE)
      }
      for(k in seq_along(pal)) {
        st<-names(pal)[k];z<-d[FinalSubtype==st];n<-unique(z$N)
        if(id=='S10M')fig4_v2_text(paste0(st,'\nn = ',n),centers[k],104,7,TRUE,col=pal[[st]])
        fig4_v2_text(if(id=='B')paste0(st,' (n = ',n,')') else st,starts[k]+fw/2,104,8,TRUE,col=pal[[st]])
        fig4_v2_line(c(trans(1,k),trans(1,k)),c(bottom,top),col='grey65',lty=2)
        for(j in seq_along(genes)) {
          row<-z[Gene==genes[j]];stopifnot(nrow(row)==1)
          yy<-top-(j-.5)*rh;frac<-row$Fraction
          if(id=='S10M') {
            bw<-14*frac/.4
            if(bw>0)fig4_v2_rect(centers[k]-7+bw/2,yy,bw,3.4,grDevices::adjustcolor(pal[[st]],alpha.f=.24))
            fig4_v2_text(sprintf('%.1f',frac*100),centers[k],yy,7)
          }
          est<-row$CMHOddsRatio;lo<-row$CMHLower95;hi<-row$CMHUpper95;q<-row$CMHFDR
          fill<-if(!is.na(q)&&q<.05)pal[[st]] else 'white'
          if(is.finite(est)&&est>0) {
            xx<-trans(pmax(lims[1],pmin(lims[2],est)),k)
            if(is.finite(lo)&&is.finite(hi)&&lo>0&&hi>0)
              fig4_v2_line(trans(pmax(lims[1],pmin(lims[2],c(lo,hi))),k),c(yy,yy),pal[[st]],.8)
            grid::grid.points(grid::unit(xx,'mm'),grid::unit(yy,'mm'),pch=21,size=grid::unit(1.8,'mm'),
              gp=grid::gpar(col=pal[[st]],fill=fill,lwd=.8))
          } else if(!is.na(est)&&(est==0||is.infinite(est))) {
            xx<-if(est==0)starts[k] else starts[k]+fw
            offsets<-if(est==0)c(0,1.5,1.5) else c(0,-1.5,-1.5)
            grid::grid.polygon(grid::unit(xx+offsets,'mm'),grid::unit(yy+c(0,.85,-.85),'mm'),
              gp=grid::gpar(col=pal[[st]],fill=fill,lwd=.8))
          } else fig4_v2_text('NA',starts[k]+fw/2,yy,7,col='grey50')
        }
        fig4_v2_line(c(starts[k],starts[k]+fw),c(bottom-1,bottom-1),'grey40')
        for(v in c(.1,1,10,100)) {
          xx<-trans(v,k);fig4_v2_line(c(xx,xx),c(bottom-1,bottom-2),'grey40')
          fig4_v2_text(as.character(v),xx,bottom-5,7)
        }
      }
      # Statistical key close to graph; inference definitions remain in the legend.
      grid::grid.points(grid::unit(c(7,40),'mm'),grid::unit(c(11,11),'mm'),pch=21,
        size=grid::unit(1.8,'mm'),gp=grid::gpar(col='grey30',fill=c('grey30','white')))
      fig4_v2_text('q < 0.05',10,11,7,just='left');fig4_v2_text('q \u2265 0.05',43,11,7,just='left')
      grid::grid.polygon(grid::unit(c(80,81.5,81.5),'mm'),grid::unit(c(11,11.85,10.15),'mm'),
        gp=grid::gpar(col='grey30',fill='white'))
      fig4_v2_text('OR = 0; CI unavailable',84,11,7,just='left')
      if(any(is.infinite(d$CMHOddsRatio))) {
        grid::grid.polygon(grid::unit(c(132,130.5,130.5),'mm'),grid::unit(c(11,11.85,10.15),'mm'),
          gp=grid::gpar(col='grey30',fill='white'))
        fig4_v2_text('OR = Inf',135,11,7,just='left')
      }
      if(id=='S10M')fig4_v2_text('Each subtype versus the other two; BH correction over 15 genes per subtype',87,4,7)
    })
  }
  if(id=='S10H')return(build_fig4_genomics('H'))
  if(id=='S10I')return(build_fig4_genomics('I'))
  if(id=='S10J')return(build_fig4_genomics('G'))
  if(id=='S10L')return(build_fig4_genomics('A'))
  stop('Unknown Figure4 v2 panel: ',id)
}

register_fig4_v2 <- function() {
  ids<-c('A','B','S10H','S10I','S10J','S10L','S10M')
  paths<-c(A='figures/main/Figure4A_driver_oncoprint.pdf',B='figures/main/Figure4B_mutation_frequency_effects.pdf',
    S10H='figures/supplement/Supplement_S10H_cnv_frequency.pdf',
    S10I='figures/supplement/Supplement_S10I_mutation_cooccurrence.pdf',
    S10J='figures/supplement/Supplement_S10J_mutation_frequency.pdf',
    S10L='figures/supplement/Supplement_S10L_mutation_record_counts.pdf',
    S10M='figures/supplement/Supplement_S10M_dong_tcga_associations.pdf')
  roles<-c(A='Driver landscape with overall and subtype frequencies',B='Subtype-associated mutation effects (Study-stratified OR/CI)',
    S10H='CNV frequency (2022Dong)',S10I='Exploratory pairwise mutation associations',
    S10J='Detailed grouped mutation frequencies',S10L='Full landscape with source variant counts',S10M='Dong + TCGA frequency association sensitivity')
  src<-c(paste0(FIG4_GENOMICS,c('mutation_display_long.csv','sample_order.csv','gene_order.csv',
    'mutation_frequency.csv','cnv_frequency.csv','pairwise_associations.csv','input_hashes.csv','output_hashes.csv')),
    '4-Multiomic/refinements/fig4_frequency_20261030/dong_tcga_frequency_associations.csv',
    '4-Multiomic/refinements/fig4_frequency_20261030/dong_tcga_frequency.csv',
    '4-Multiomic/refinements/fig4_frequency_20261030/input_hashes.csv',
    '4-Multiomic/scripts/prepare_fig4_frequency_sensitivity.R','8-Visualization/scripts/render_fig4_genomics_v2.R',
    '8-Visualization/scripts/render_fig4_genomics.R','8-Visualization/scripts/fig_style.R',
    '8-Visualization/scripts/stat_annotate.R','8-Visualization/scripts/figure_config.json',
    'governance/data/subtype_label_color_dictionary.csv','0-Data/qc/study_color_dictionary.csv')
  old<-read_source('8-Visualization/tables/figure_manifest.csv');base<-copy(old[Figure=='Figure 4'&Panel=='A'][1])
  for(id in ids) {
    row<-copy(base);main<-id%in%c('A','B')
    row[,`:=`(Figure=if(main)'Figure 4' else 'Supplement',Panel=id,Role=roles[[id]],SourceFigure='',
      SourceTable=paste(src,collapse=';'),TargetOutput=paste0('8-Visualization/',paths[[id]]),
      Status='AnalysisReady',SubtypeScope=if(id=='S10H')'2022Dong CNV main n=219' else if(id=='S10M')'Dong + TCGA DNA main n=248 (54/122/72)' else 'DNA main n=308 (63/153/92)',
      WidthMM='174',HeightMM=as.character(c(A=102,B=128,S10H=128,S10I=91,S10J=96,S10L=132,S10M=128)[[id]]),
      Notes='20261030 v2 review draft. All 15 drivers retained. Frozen primary CMH/FDR unchanged; source coverage limits remain in legends.',
      SourceMD5='',QAFile='',RenderedAt='')]
    idx<-which(old$Figure==row$Figure & old$Panel==id)
    if(length(idx))old[idx]<-row else old<-rbindlist(list(old,row),use.names=TRUE)
  }
  old<-old[!(Figure=='Figure 4'&Panel%in%c('G','H','I')&
    grepl('Figure4[GHI]_(mutation_frequency|cnv_frequency|mutation_cooccurrence)[.]pdf$',TargetOutput))]
  tmp<-tempfile(tmpdir=file.path(STEP8_ROOT,'8-Visualization/tables'),fileext='.csv')
  fwrite(old,tmp,na='');stopifnot(file.rename(tmp,file.path(STEP8_ROOT,'8-Visualization/tables/figure_manifest.csv')))
}

render_fig4_v2 <- function(id) {
  if(id=='S10L'&&'IRkernel'%in%loadedNamespaces()) {
    script<-tempfile(fileext='.R');on.exit(unlink(script),add=TRUE)
    writeLines(c(paste0('root <- ',paste(capture.output(dput(STEP8_ROOT)),collapse='')),
      sprintf('source(file.path(root,"%s"))',paste0('8-Visualization/scripts/',c('fig_style.R','stat_annotate.R',
        'render_bulk_figures.R','render_fig4_genomics.R','render_fig4_genomics_v2.R'))),'render_fig4_v2("S10L")'),script)
    stopifnot(system2(file.path(R.home('bin'),'exec/R'),c('--vanilla','--slave','-f',shQuote(script)))==0)
    return(invisible(NULL))
  }
  figure<-if(id%in%c('A','B'))'Figure 4' else 'Supplement'
  row<-preflight_subfigure(figure,id)
  save_panel(build_fig4_v2(id),row$TargetOutput,as.numeric(row$WidthMM),as.numeric(row$HeightMM))
  src<-strsplit(row$SourceTable,';',fixed=TRUE)[[1]]
  hashes<-paste(paste0(src,'=',unname(tools::md5sum(file.path(STEP8_ROOT,src)))),collapse=';')
  m<-read_source('8-Visualization/tables/figure_manifest.csv')
  m[Figure==figure&Panel==id,`:=`(Status='Rendered',SourceMD5=hashes,QAFile='',RenderedAt=format(Sys.time(),tz='UTC',usetz=TRUE))]
  tmp<-tempfile(tmpdir=file.path(STEP8_ROOT,'8-Visualization/tables'),fileext='.csv')
  fwrite(m,tmp,na='');stopifnot(file.rename(tmp,file.path(STEP8_ROOT,'8-Visualization/tables/figure_manifest.csv')))
  message('Rendered: ',row$TargetOutput)
}
