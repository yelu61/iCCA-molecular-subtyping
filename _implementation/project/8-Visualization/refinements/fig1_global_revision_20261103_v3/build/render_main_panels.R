# Display-only redraw of correspondence and survival panels from audited tables.
# No expression transformation, clustering, survival model or contrast is refitted.
options(stringsAsFactors=FALSE, warn=2)
library(grid)
library(svglite)
library(ragg)
library(jsonlite)
out <- normalizePath('8-Visualization/refinements/fig1_global_revision_20261103_v3')
old <- '8-Visualization/refinements/fig1_GH_enlarged_20261102_v2/source_data'
e <- read.csv('8-Visualization/refinements/fig1_slide_review_20261102_v1/source_data/Figure1E_combined.csv')
curve <- read.csv(file.path(old,'G_KM_coordinates.csv'))
risk <- read.csv(file.path(old,'G_risk_counts.csv'))
rmst <- read.csv(file.path(old,'H_RMST_differences.csv'))
pal <- read.csv('governance/data/subtype_label_color_dictionary.csv')
pal <- setNames(pal$Color,pal$FinalSubtype)
study_dictionary <- read.csv('0-Data/qc/study_color_dictionary.csv')
stopifnot(nrow(e)==39, nrow(risk)==27, nrow(rmst)==3,
          all(e$N/e$RowN==e$RowFraction | abs(e$N/e$RowN-e$RowFraction)<1e-14),
          sum(risk$AtRisk[risk$Month==0])==331, all(rmst$TauMonths==24))
write.csv(e,file.path(out,'source_data/Figure1E_all_cells.csv'),row.names=FALSE)
write.csv(curve,file.path(out,'source_data/Figure1G_KM_coordinates.csv'),row.names=FALSE)
write.csv(risk,file.path(out,'source_data/Figure1G_risk_counts.csv'),row.names=FALSE)
write.csv(rmst,file.path(out,'source_data/SFigure1L_RMST_contrasts.csv'),row.names=FALSE)
mm <- function(x) unit(x,'mm')
txt <- function(t,x,y,size=6.5,bold=FALSE,col='#222222',just='left',rot=0)
  grid.text(t,mm(x),mm(y),just=just,rot=rot,gp=gpar(fontfamily='Arial',fontsize=size,fontface=if(bold)'bold' else 'plain',col=col,lineheight=1.06))
ln <- function(x,y,col='#555555',lwd=.5,lty=1) grid.lines(mm(x),mm(y),gp=gpar(col=col,lwd=lwd,lty=lty))
rect <- function(x,y,w,h,fill) grid.rect(mm(x),mm(y),mm(w),mm(h),just=c('left','bottom'),gp=gpar(fill=fill,col=NA))
page <- function(w,h){grid.newpage();pushViewport(viewport(x=0,y=0,width=mm(w),height=mm(h),just=c('left','bottom')))}
export <- function(draw,name,w,h){
  cairo_pdf(file.path(out,'figures',paste0(name,'.pdf')),width=w/25.4,height=h/25.4,family='Arial');draw();dev.off()
  svglite(file.path(out,'figures',paste0(name,'.svg')),width=w/25.4,height=h/25.4,bg='white');draw();dev.off()
  agg_png(file.path(out,'figures',paste0(name,'.png')),width=w,height=h,units='mm',res=600,background='white');draw();dev.off()
}
drawE <- function(){
  page(174,44)
  ids <- c('Figure1E1','Figure1E3','Figure1E2','Figure1E4')
  heads <- c('Ahn','Dong: RNA','Dong: protein','Deng: protein')
  subs <- c('Original A/B classes','IG1–IG3','S1–S4','S-I–S-III')
  shade <- colorRamp(c('#FAFCFD','#0072A8'))
  for(i in seq_along(ids)){
    d <- e[e$PanelSource==ids[i],]; x <- (i-1)*44
    txt(heads[i],x,42.3,6.8,TRUE)
    txt(paste0('n = ',unique(d$TotalN)),x+41,42.3,6.2,just='right')
    txt(subs[i],x,38.8,6.4)
    txt('Class',x+.5,34.4,6.1)
    for(j in 1:3){rect(x+9+(j-1)*8,32.7,7.9,3.4,pal[paste0('C',j)]);txt(paste0('C',j),x+12.95+(j-1)*8,34.4,6.2,TRUE,just='centre')}
    txt('Row n',x+38.1,34.4,5.8,just='centre')
    rows <- unique(d[order(d$RowOrder),c('OriginalLabel','RowOrder')])
    cellh<-23/nrow(rows)
    for(k in seq_len(nrow(rows))){
      r <- d[d$OriginalLabel==rows$OriginalLabel[k],]; y <- 32-cellh*(k-.5)
      label <- sub('class: ','',rows$OriginalLabel[k]);label <- gsub('Not reported','Not\nreported',label)
      txt(label,x+4.3,y,6.0,just='centre')
      for(j in 1:3){z<-r[r$FinalSubtype==paste0('C',j),]; rgbv<-shade(z$RowFraction)/255
        fill<-if(grepl('Not',label))'#F0F1F2' else rgb(rgbv[1],rgbv[2],rgbv[3])
        rect(x+9+(j-1)*8,y-cellh/2+.08,7.9,cellh-.16,fill)
        txt(sprintf('%d(%d%%)',z$N,round(z$RowFraction*100)),x+12.95+(j-1)*8,y,6.2,
            col=if(z$RowFraction>.65 && !grepl('Not',label))'white' else '#222222',just='centre')
      }
      txt(unique(r$RowN),x+38.1,y,6.2,just='centre')
    }
    txt('Total',x+.3,6.5,6.0)
    for(j in 1:3)txt(unique(d$ColumnN[d$FinalSubtype==paste0('C',j)]),x+12.95+(j-1)*8,6.5,6.2,TRUE,just='centre')
    txt(unique(d$TotalN),x+38.1,6.5,6.2,TRUE,just='centre')
  }
  txt('Cells: n (% within original class). Dong RNA/protein panels share the same 223 samples.',0,2.5,6.2,col='#555555')
  popViewport()
}
drawG <- function(){
  page(76,53)
  txt('Overall survival',0,51.2,7.5,TRUE)
  txt('n = 331; events = 153',0,47.5,6.4)
  txt('Study-stratified log-rank P < 0.001',0,44.1,6.4)
  for(j in 1:3){x<-42+(j-1)*10;ln(c(x,x+3),c(48.1,48.1),pal[paste0('C',j)],1);txt(paste0('C',j),x+3.7,48.1,6.3)}
  xx<-function(v)10+v/96*63; yy<-function(v)21.5+v*19
  ln(c(10,73),c(21.5,21.5));ln(c(10,10),c(21.5,40.5))
  for(v in c(0,.25,.5,.75,1)){ln(c(9.2,10),rep(yy(v),2));txt(sprintf('%.2f',v),8.5,yy(v),6,just='right')}
  times <- seq(0,96,12)
  for(t in times){ln(rep(xx(t),2),c(21.5,20.8));txt(t,xx(t),19.1,6.1,just='centre')}
  txt('Survival probability',1,30,6.3,rot=90,just='centre')
  txt('Overall survival (months)',41.5,15.6,6.3,just='centre')
  for(c in paste0('C',1:3)){
    q<-curve[curve$Subtype==c,];prev<-c(1,head(q$Survival,-1))
    ln(xx(c(0,rep(q$Month,each=2))),yy(c(1,as.vector(rbind(prev,q$Survival)))),pal[c],.85)
    cen<-q[q$Censored>0,];for(k in seq_len(nrow(cen)))ln(rep(xx(cen$Month[k]),2),yy(cen$Survival[k])+c(-.3,.3),pal[c],.5)
  }
  txt('Number at risk',0,12.5,6.3)
  for(j in 1:3){y<-9.5-(j-1)*3;c<-paste0('C',j);txt(c,0,y,6.3,col=pal[c]);r<-risk[risk$Subtype==c,];for(k in seq_len(nrow(r)))txt(r$AtRisk[k],xx(r$Month[k]),y,6.2,just='centre')}
  popViewport()
}
drawL <- function(){
  page(174,47)
  txt('Restricted mean survival over 0–24 months',1,44.3,8.1,TRUE)
  txt('Available-case population: n = 331; events = 153. Unadjusted pairwise differences.',1,39.7,7)
  xx<-function(v)39+(v+10)/16*80
  ln(rep(xx(0),2),c(12.5,35),'#AAAAAA',.6,2)
  ys<-c(31,23,15)
  txt('Difference [95% CI], months',132,35.6,7)
  for(i in 1:3){r<-rmst[i,];y<-ys[i]
    txt(gsub(' - ',' − ',r$Contrast),2,y,7.5)
    ln(xx(c(r$Lower95,r$Upper95)),c(y,y),'#555555',1)
    grid.points(mm(xx(r$DifferenceMonths)),mm(y),pch=16,size=mm(1.5),gp=gpar(col='#426782'))
    txt(sprintf('%+.2f [%.2f, %.2f]',r$DifferenceMonths,r$Lower95,r$Upper95),132,y,7.2)
  }
  ln(c(xx(-10),xx(6)),c(11,11))
  for(v in c(-10,-5,0,5)){ln(rep(xx(v),2),c(11,10.2));txt(v,xx(v),8.2,7,just='centre')}
  txt('RMST difference in months (first − second)',79,3.5,7,just='centre')
  popViewport()
}
export(drawE,'Figure1E_correspondence',174,44)
export(drawG,'Figure1G_survival',76,53)
export(drawL,'SFigure1L_RMST',174,47)
ari <- read.csv('8-Visualization/refinements/fig1_sfig1_typography_20261030/source_data/SFig1F.csv')
cross <- read.csv('8-Visualization/refinements/fig1_sfig1_typography_20261030/source_data/SFig1G.csv')
stopifnot(nrow(ari)==16,nrow(cross)==9,sum(cross$N)==586)
write.csv(ari,file.path(out,'source_data/SFigure1E_ARI.csv'),row.names=FALSE)
write.csv(cross,file.path(out,'source_data/SFigure1F_crosswalk.csv'),row.names=FALSE)
drawSE <- function(){
 page(59,52)
 txt('Agreement between stored',0,50,7.2,TRUE);txt('method assignments',0,46.5,7.2,TRUE)
 methods<-c('ConsensusCluster','NMF','KMeans','HC');labels<-c('Consensus','NMF','k-means','Hierarchical')
 short<-c('Cons.','NMF','k-means','HC');shade<-colorRamp(c('#FFFFFF','#0072B2'))
 for(i in 1:4){y<-39.5-(i-1)*8;txt(labels[i],15,y,6.5,just='right')
  for(j in 1:4){a<-ari$AdjustedRandIndex[ari$MethodA==methods[i]&ari$MethodB==methods[j]];col<-shade(max(a,0))/255
   rect(17+(j-1)*10,y-3.85,9.8,7.7,rgb(col[1],col[2],col[3]));txt(if(a<.01)'<0.01' else sprintf('%.2f',a),21.9+(j-1)*10,y,6.5,col=if(a>.6)'white' else '#222222',just='centre')
  }
 }
 for(j in 1:4)txt(short[j],21.9+(j-1)*10,9.4,6.3,just='centre')
 txt('ARI: 0 (white) to 1 (blue)',0,5.3,6.5)
 txt('Cons. = consensus; HC = hierarchical',0,1.7,6.1)
 popViewport()
}
drawSF <- function(){
 page(52,45.6)
 txt('NMF–k-means state definition',0,43.8,7.2,TRUE)
 txt('n = 586; labels aligned to C1/C2/C3',0,40.2,6.5)
 kval<-c(3,1,2)
 for(i in 1:3){y<-32-(i-1)*9;txt(paste0('N',i),1,y,6.5,TRUE)
  for(j in 1:3){d<-cross[cross$NMF==i&cross$KMeans==kval[j],];c<-if(d$Ensemble=='None')'#F4F4F4' else pal[d$Ensemble]
   rect(9+(j-1)*13.6,y-4.3,13.4,8.6,c)
   lab<-if(d$Ensemble=='None')'—' else gsub('_','/',d$Ensemble)
   txt(lab,15.7+(j-1)*13.6,y+1.8,6.2,just='centre')
   txt(paste0('n = ',d$N),15.7+(j-1)*13.6,y-1.8,6.5,just='centre')
  }
 }
 for(j in 1:3)txt(paste0('K',kval[j]),15.7+(j-1)*13.6,6.9,6.5,TRUE,just='centre')
 txt('N: NMF cluster; K: k-means cluster',0,2.2,6.3)
 popViewport()
}
export(drawSE,'SFigure1E_method_agreement',59,52)
export(drawSF,'SFigure1F_state_definition',52,45.6)
write_json(list(no_refit=TRUE,E_cells=39,E_total_populations=c(25,223,223,64),
                G_patients=331,G_events=153,G_risk_rows=27,L_contrasts=3),
           file.path(out,'review/main_panel_checks.json'),pretty=TRUE,auto_unbox=TRUE)
