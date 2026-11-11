# R-only, editable-text PDF rendering of the exploratory duct-by-C3 extension.
suppressPackageStartupMessages({library(data.table);library(grid);library(jsonlite)})
source('8-Visualization/scripts/fig_style.R')
out <- '6-ClinicalSurvival/review/20261101_duct_c3_stratified_v1'
z <- fread(file.path(out,'tables/analysis_population.csv'))
curves <- fread(file.path(out,'tables/km_curves.csv'))
risk <- fread(file.path(out,'tables/km_risk_table.csv'))
effects <- fread(file.path(out,'tables/small_duct_c3_effects.csv'))
counts <- fread(file.path(out,'tables/group_counts.csv'))
pal <- load_subtype_palette();invisible(load_study_palette())
cols <- c('Non-C3'='#596873','C3'=pal$colors[['C3']])
active <- '';height <- 0;checking <- FALSE;bounds <- list()
start <- function(h){height<<-h;grid.newpage();pushViewport(viewport(xscale=c(0,174),yscale=c(h,0)))}
txt <- function(s,x,y,size=7.2,bold=FALSE,just='left',color='#263641',rot=0){
 g<-textGrob(s,x=unit(x,'native'),y=unit(y,'native'),just=c(just,'center'),rot=rot,
 gp=gpar(fontfamily='Arial',fontsize=size,fontface=if(bold)'bold' else 'plain',col=color))
 if(checking && rot==0){w<-convertWidth(grobWidth(g),'mm',valueOnly=TRUE);h<-convertHeight(grobHeight(g),'mm',valueOnly=TRUE)
 lo<-x-switch(just,left=0,center=w/2,right=w)
 if(lo<0 || lo+w>174 || y-h/2<0 || y+h/2>height)stop('Text outside bounds: ',s)
 stopifnot(size>=7)
 bounds[[length(bounds)+1L]]<<-data.table(Panel=active,Text=s,XMin=lo,XMax=lo+w,YMin=y-h/2,YMax=y+h/2,FontPt=size)}
 grid.draw(g)
}
ln <- function(x,y,col='#A8B1B8',lwd=.65,lty=1)grid.lines(unit(x,'native'),unit(y,'native'),gp=gpar(col=col,lwd=lwd,lty=lty))
pt <- function(x,y,col,pch=16,sz=1.8)grid.points(unit(x,'native'),unit(y,'native'),pch=pch,size=unit(sz,'mm'),gp=gpar(col=col,fill=col))
step <- function(d,field){v<-d[[field]];data.table(t=c(0,rep(d$Time,each=2)),v=c(1,as.vector(rbind(c(1,head(v,-1)),v))))}
draw_km <- function(duct,x0,x1,y0,y1,risk_y,show_y=TRUE){
 xmax<-max(risk$Month);xx<-function(v)x0+(x1-x0)*v/xmax;yy<-function(v)y1-(y1-y0)*v
 for(v in c(0,.5,1)){ln(c(x0,x1),rep(yy(v),2),'#E7EBEE',.5)
  if(show_y)txt(sprintf('%.1f',v),x0-2,yy(v),7,just='right')}
 for(m in c('Non-C3','C3')){
  d<-curves[Duct==duct & Molecular==m];setorder(d,Time)
  lo<-step(d,'Lower95');hi<-step(d,'Upper95');ok<-is.finite(lo$v)&is.finite(hi$v)
  grid.polygon(unit(c(xx(lo$t[ok]),rev(xx(hi$t[ok]))),'native'),unit(c(yy(lo$v[ok]),rev(yy(hi$v[ok]))),'native'),
   gp=gpar(fill=adjustcolor(cols[m],alpha.f=.13),col=NA))
  a<-step(d,'Survival');ln(xx(a$t),yy(a$v),cols[m],1.25)
  c<-d[Censors>0];if(nrow(c))pt(xx(c$Time),yy(c$Survival),cols[m],pch=3,sz=1.3)
 }
 ln(c(x0,x1),c(y1,y1),'#65747E');ln(c(x0,x0),c(y0,y1),'#65747E')
 for(t in risk[!duplicated(Month),Month]){ln(rep(xx(t),2),c(y1,y1+1),'#65747E');txt(as.character(t),xx(t),y1+4,7,just='center')}
 txt('OS (months)',mean(c(x0,x1)),y1+9,7.2,just='center')
 txt('At risk',x0-18,risk_y,7.1,bold=TRUE)
 for(j in seq_along(cols)){m<-names(cols)[j];y<-risk_y+5+j*4.5
  txt(m,x0-18,y,7,color=cols[m]);a<-risk[Duct==duct & Molecular==m]
  for(i in seq_len(nrow(a)))txt(a$AtRisk[i],xx(a$Month[i]),y,7,just='center',color=cols[m])
 }
}
draw_four <- function(){start(113);on.exit(popViewport())
 txt('Overall survival by observed duct morphology and C3',4,6,9,TRUE)
 txt('79 patients / 45 deaths | descriptive Kaplan–Meier curves with pointwise 95% CIs',4,12,7.2)
 for(j in 1:2){duct<-c('Small duct','Large duct')[j];x0<-c(23,110)[j];x1<-c(81,168)[j]
  a<-counts[Duct==duct];txt(sprintf('%s  |  n = %d',duct,sum(a$N)),x0,20,8,TRUE)
  txt(sprintf('C3: %d (%d deaths)   Non-C3: %d (%d deaths)',a[Molecular=='C3',N],a[Molecular=='C3',Events],a[Molecular=='Non-C3',N],a[Molecular=='Non-C3',Events]),x0,26,7)
  draw_km(duct,x0,x1,32,73,87,TRUE)
 }
 txt('Survival',5,53,7.2,rot=90,just='center')
 txt('C3',80,107,7.2,color=cols['C3']);ln(c(73,78),c(107,107),cols['C3'],1.2)
 txt('Non-C3',102,107,7.2,color=cols['Non-C3']);ln(c(95,100),c(107,107),cols['Non-C3'],1.2)
}
draw_small <- function(){start(112);on.exit(popViewport())
 txt('C3 within observed small-duct iCCA',4,6,9,TRUE)
 txt('61 patients / 33 deaths | C3 n = 15; non-C3 n = 46 | exploratory subgroup',4,12,7.2)
 txt('Observed survival',23,21,8,TRUE)
 draw_km('Small duct',23,83,30,72,85,TRUE)
 txt('Survival',5,51,7.2,rot=90,just='center')
 txt('Adjusted association',99,21,8,TRUE)
 txt('C3 versus non-C3: HR [95% CI]',99,27,7.2)
 xx<-function(v){stopifnot(all(is.finite(v)),all(v>0));102+(165-102)*(log(v)-log(.5))/(log(12)-log(.5))}
 labs<-c('Proportional hazards','Time-varying: 12 months','Time-varying: 24 months')
 for(i in 1:3){e<-effects[i];y<-43+(i-1)*17
  ln(rep(xx(1),2),c(y-4,y+4),'#ABB3BA',.65,2)
  txt(labs[i],99,y-7,7.2,bold=i==1)
  ln(xx(c(e$Lower95,e$Upper95)),c(y,y),cols['C3'],1.1)
  pt(xx(e$HR),y,cols['C3'],if(i==1)16 else 1)
  txt(sprintf('%.2f [%.2f, %.2f]',e$HR,e$Lower95,e$Upper95),168,y+5,7.1,just='right')
 }
 ln(xx(c(.5,12)),c(86,86),'#65747E')
 for(v in c(.5,1,2,4,8)){ln(rep(xx(v),2),c(86,87),'#65747E');txt(format(v),xx(v),91,7,just='center')}
 txt('Hazard ratio (log scale)',mean(xx(c(.5,12))),97,7.2,just='center')
 txt('Adjusted for age, sex and stage; study-specific baseline hazards. Time-varying C3 and age sensitivity.',4,108,7)
}
draws<-list(Duct_C3_four_group_KM=draw_four,Small_duct_C3_association=draw_small)
heights<-c(113,112);manifest<-list()
for(i in seq_along(draws)){
 active<-names(draws)[i];pdf<-file.path(out,'figures',paste0(active,'.pdf'));checking<-TRUE
 save_panel(draws[[i]],pdf,174,heights[i]);checking<-FALSE
 png<-file.path(out,'review',paste0(active,'.png'))
 ragg::agg_png(png,width=174,height=heights[i],units='mm',res=300,background='white');draws[[i]]();dev.off()
 manifest[[i]]<-data.table(Panel=active,PDF=pdf,Preview=png,WidthMM=174,HeightMM=heights[i],
 Status='RENDERED_EXPLORATORY_REVIEW',PDFMD5=unname(tools::md5sum(pdf)))
}
fwrite(rbindlist(manifest),file.path(out,'figure_manifest.csv'))
fwrite(rbindlist(bounds),file.path(out,'review/figure_text_bounds.csv'))
