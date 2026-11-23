# Style-only narrowing of the principal all-module supplementary panel.
# Frozen statistics, labels, fonts and vertical geometry are unchanged.
suppressPackageStartupMessages({library(data.table);library(grid)})
source('8-Visualization/scripts/fig_style.R')
out <- '8-Visualization/refinements/sfig2_wgcna_narrow_20261031_v37'
for(d in c('figures','tables','review'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
pal <- load_subtype_palette()$colors
states <- c('C1','C2','C3')
src <- c('8-Visualization/refinements/sfig2_compact_20261031_v35/tables/module_order.csv','8-Visualization/refinements/sfig2_compact_20261031_v35/tables/network_principal_effects.csv')
hashes <- tools::md5sum(src)
cl <- fread(src[1]);me <- fread(src[2])
stopifnot(nrow(cl)==24,nrow(me)==72,sum(cl$MainFigure)==11)
groups <- c('C1_specific','C2_specific','C3_specific','SharedHighCompatible_C1_C2','SharedHighCompatible_C2_C3','Replicated_gradient_or_ambiguous','Non_specific')
group_labels <- c('C1-specific','C2-specific','C3-specific','Shared-high C1 / C2','Shared-high C2 / C3','Gradient / ambiguous','Non-specific')
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
module_draw <- function(is_state=FALSE){function(){
 grid.newpage();texts<<-list();marks<<-list();H<-153;W<-if(is_state)115 else 70
 xs<-if(is_state)c(43,54.5,66,81,92.5,104) else c(43,53,63);cw<-if(is_state)10.8 else 9.2
 for(j in seq_along(xs)){txt(gsub('_','/',states[j]),xs[j],H-5,6.4,'bold','center');rect(xs[j],H-8,cw,1,pal[states[j]]);if(is_state)txt(paste0('n=',ns[j]),xs[j],H-11.5,5.2,just='center')}
 txt('Module',4,H-5,6.3,'bold');txt('Genes',33,H-5,5.8,just='right');y<-H-18
 for(g in seq_along(groups)){
  txt(group_labels[g],4,y,5.8,'bold',col='#56636D');line(c(4,W-4),rep(y-2,2),'#E4E8EA');y<-y-4
  for(i in cl[GroupOrder==g,RowOrder]){r<-cl[RowOrder==i];rect(4.8,y,1.5,3.2,r$ModuleColor)
   txt(r$ModuleColor,7,y,5.9,face=if(r$MainFigure)'bold'else'plain');txt(r$ModuleSize,33,y,5.3,just='right',col='#66717A')
   if(r$MainFigure)point(36,y,size=.65)
   for(j in seq_along(xs)){
    q<-if(is_state)state[Module==r$Module & FinalSubtype==states[j]] else me[Module==r$Module & TargetSubtype==states[j]]
    v<-if(is_state)q$Median else q$RandomEffect;stopifnot(length(v)==1,is.finite(v));rect(xs[j],y,cw,3.1,colval(v))
    if(is_state)txt(sprintf('%.2f',if(abs(v)<.005)0 else v),xs[j],y,5.2,just='center',col=if(abs(v)>1.3)'white'else'#24282B') else if(q$FDR<.05)point(xs[j],y,size=.65,col=if(abs(v)>1.3)'white'else'#24282B')
    marks[[length(marks)+1L]]<<-data.table(Module=r$Module,State=states[j],Value=v,X=xs[j],Y=y)
   };y<-y-3.4
  };y<-y-1
 }
 stopifnot(y>15)
 bar(5,8,if(is_state)'Median frozen score (SD)'else'Meta-analysis effect (SD)',w=27)
 point(43,12,size=.65);txt('Main Fig. 2C module',45,12,5.5)
 if(!is_state){point(43,7,size=.65);txt('In tile: FDR < 0.05',45,7,5.5)}
}}
emit('SFig_principal_B_all_WGCNA',module_draw(),70,153,'Principal support','B')

fwrite(cl,file.path(out,'tables/module_order.csv'));fwrite(me,file.path(out,'tables/network_principal_effects.csv'))
fwrite(rbindlist(manifest),file.path(out,'figure_manifest.csv'))
stopifnot(identical(hashes,tools::md5sum(src)))
fwrite(data.table(Source=src,BeforeMD5=unname(hashes),AfterMD5=unname(tools::md5sum(src))),file.path(out,'review/source_hashes.csv'))
