# Style-only revision: 3 rows by 2 columns; frozen sample scores and summary unchanged.
suppressPackageStartupMessages({library(data.table);library(grid)})
source('8-Visualization/scripts/fig_style.R')
out <- '8-Visualization/refinements/sfig2_patterns_3x2_20261031_v36'
for(d in c('figures','tables','review'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
pal <- load_subtype_palette()$colors
states <- c('C1','C2','C3','C1_C2','C1_C3','C2_C3')
src <- '8-Visualization/refinements/sfig2_compact_20261031_v35/tables/pattern_scores.csv'
hashes <- tools::md5sum(src)
scores <- fread(src)
stopifnot(nrow(scores)==3516,uniqueN(scores$Sample)==586)
mm <- function(x)unit(x,'mm'); texts<-list();marks<-list(); manifest<-list()
txt <- function(s,x,y,size=6,face='plain',just='left',col='#24282B',rot=0){grid.text(s,mm(x),mm(y),just=just,rot=rot,gp=gpar(fontfamily='Arial',fontsize=size,fontface=face,col=col));texts[[length(texts)+1L]]<<-data.table(Text=as.character(s),X=x,Y=y,Size=size)}
rect <- function(x,y,w,h,fill,border=NA)grid.rect(mm(x),mm(y),mm(w),mm(h),gp=gpar(fill=fill,col=border,lwd=.4))
line <- function(x,y,col='#D3D8DC',lwd=.5,lty=1)grid.lines(mm(x),mm(y),gp=gpar(col=col,lwd=lwd,lty=lty))
point <- function(x,y,pch=16,size=1,col='#303840',fill=col)grid.points(mm(x),mm(y),pch=pch,size=mm(size),gp=gpar(col=col,fill=fill,lwd=.55))
ramp<-colorRamp(c(SCALE_EFFECT_LOW,'#FAFAFA',SCALE_EFFECT_HIGH))
colval<-function(v,lim=2)rgb(ramp(pmax(0,pmin(1,(v/lim+1)/2))),maxColorValue=255)
bar <- function(x,y,label,w=27,lim=2){txt(label,x,y+4,5.8);for(k in 1:100)rect(x+(k-.5)*w/100,y,w/100+.01,1.8,colval(-lim+(k-1)/99*2*lim,lim));for(v in c(-lim,0,lim))txt(v,x+(v+lim)/(2*lim)*w,y-3,5.3,just='center')}
emit <- function(name,draw,w,h,figure,panel){
 pdf<-file.path(out,'figures',paste0(name,'.pdf'))
 grDevices::cairo_pdf(pdf,width=w/25.4,height=h/25.4,family='Arial');draw();dev.off()
 fwrite(rbindlist(texts),file.path(out,'tables',paste0(name,'_text.csv')));if(length(marks))fwrite(rbindlist(marks,fill=TRUE),file.path(out,'tables',paste0(name,'_marks.csv')))
 ragg::agg_png(sub('pdf$','png',pdf),width=w,height=h,units='mm',res=600,background='white');draw();dev.off()
 ragg::agg_png(file.path(out,'review',paste0(name,'.png')),width=w,height=h,units='mm',res=180,background='white');draw();dev.off()
 manifest[[length(manifest)+1L]]<<-data.table(Figure=figure,Panel=panel,TargetOutput=pdf,WidthMM=w,HeightMM=h,Status='Rendered review candidate',SourceTable=paste(src,collapse=';'),SourceMD5=paste(hashes,collapse=';'))
}
draw_pattern <- function(){
 grid.newpage();texts<<-list();marks<<-list();rg<-range(scores$Score);rg<-rg+c(-.15,.15)
 for(i in 1:6){p<-paste0('P',i);x0<-10+((i-1)%%2)*60;y0<-106-floor((i-1)/2)*46;ww<-49;hh<-31;fy<-function(v)y0+(v-rg[1])/diff(rg)*hh
  txt(p,x0,y0+hh+3.4,7,'bold');txt(paste0(unique(scores[Pattern==p,Genes]),' genes'),x0+6,y0+hh+3.4,5.6,col='#68727A')
  line(c(x0,x0,x0+ww),c(y0+hh,y0,y0),'#69757D');line(c(x0,x0+ww),rep(fy(0),2),'#CED5DB',.5,2)
  for(v in pretty(rg,3))if(v>=rg[1]&v<=rg[2])txt(v,x0-1.1,fy(v),5.2,just='right')
  xcs<-x0+c(3.6,11.3,19,30,37.7,45.4);line(rep(x0+24.5,2),c(y0,y0+hh),'#E6E9EB',.4)
  for(j in 1:6){z<-scores[Pattern==p&FinalSubtype==states[j]];setorder(z,Sample);xx<-xcs[j]+((seq_len(nrow(z))*.61803398875)%%1-.5)*4.8;yy<-fy(z$Score)
   point(xx,yy,size=.55,col=adjustcolor(pal[states[j]],alpha.f=.43));q<-quantile(z$Score,c(.25,.5,.75));line(rep(xcs[j],2),fy(q[c(1,3)]),'#20272D',1);line(xcs[j]+c(-1.1,1.1),rep(fy(q[2]),2),'#20272D',1.3)
   txt(gsub('_','/',states[j]),xcs[j],y0-3,5.2,just='center');marks[[length(marks)+1L]]<<-data.table(Pattern=p,Sample=z$Sample,State=states[j],Score=z$Score,X=xx,Y=yy)
  }
 }
 txt('Pattern score (mean reference z-score)',2,77,5.7,rot=90,just='center')
 txt('Samples: C1 112   C2 223   C3 162   C1/C2 41   C1/C3 39   C2/C3 9',10,5.2,5.6)
 txt('Points: samples; black marks: median and IQR',10,1.5,5.2)
}
emit('SFig_states_B_pattern_distributions',draw_pattern,123,146,'Intermediate states','B')

fwrite(rbindlist(manifest),file.path(out,'figure_manifest.csv'))
fwrite(scores,file.path(out,'tables/pattern_scores.csv'))
fwrite(scores[,.(N=.N,Q25=quantile(Score,.25),Median=median(Score),Q75=quantile(Score,.75)),by=.(Pattern,FinalSubtype)],file.path(out,'tables/pattern_summary.csv'))
stopifnot(identical(hashes,tools::md5sum(src)))
fwrite(data.table(Source=src,BeforeMD5=unname(hashes),AfterMD5=unname(tools::md5sum(src))),file.path(out,'review/source_hashes.csv'))
