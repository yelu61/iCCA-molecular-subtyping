suppressPackageStartupMessages({library(data.table);library(grid);library(ragg)})
source('8-Visualization/scripts/fig_style.R')
pal<-load_subtype_palette()$colors;invisible(load_study_palette())
out<-'8-Visualization/refinements/fig6_manuscript_support_20261103_v23'
source_path<-'8-Visualization/refinements/fig6_evidence_20261103_v21/tables/Path09_small_duct_adjusted_association.csv'
d<-fread(source_path);stopifnot(nrow(d)==3,all(d$N==61),all(d$Events==33))
fwrite(d,file.path(out,'tables/SClinID_B.csv'))
ink<-'#263E4B'
tx<-function(label,x,y,pt=6.5,bold=FALSE,col=ink){grid.text(label,x=unit(x,'mm'),y=unit(52-y,'mm'),just=c('left','center'),gp=gpar(fontfamily='Arial',fontsize=pt,fontface=if(bold)'bold'else'plain',col=col))}
xx<-function(hr)22+34*log(hr/.5)/log(24)
draw<-function(){grid.newpage();tx('Adjusted C3 association within small-duct disease',2,4,7,TRUE);tx('61 patients / 33 deaths; C3 15, non-C3 46',2,9,6.5)
  grid.lines(unit(c(xx(1),xx(1)),'mm'),unit(52-c(14,36),'mm'),gp=gpar(col='#B6C2C9',lty=2,lwd=.6))
  labs<-c('PH reference','12 months','24 months');ys<-c(17,25,33)
  for(i in 1:3){tx(labs[i],2,ys[i]);grid.lines(unit(xx(c(d$Lower95[i],d$Upper95[i])),'mm'),unit(rep(52-ys[i],2),'mm'),gp=gpar(col=pal['C3'],lwd=1));grid.points(unit(xx(d$HR[i]),'mm'),unit(52-ys[i],'mm'),pch=if(i==1)1 else 16,size=unit(1.7,'mm'),gp=gpar(col=pal['C3']));tx(sprintf('%.2f [%.2f, %.2f]',d$HR[i],d$Lower95[i],d$Upper95[i]),59,ys[i],6.2)}
  ticks<-c(.5,1,2,4,8);grid.lines(unit(c(xx(.5),xx(12)),'mm'),unit(c(14,14),'mm'),gp=gpar(col=ink,lwd=.6))
  for(t in ticks)tx(as.character(t),xx(t)-1,41,6.2)
  tx('Hazard ratio (95% CI; log scale)',22,46,6.5);tx('Age/sex/stage-adjusted; study-stratified. Exploratory.',2,50,6)
}
cairo_pdf(file.path(out,'figures/SClinID_B.pdf'),width=85/25.4,height=52/25.4,family='Arial');draw();dev.off()
agg_png(file.path(out,'review/SClinID_B.png'),width=85,height=52,units='mm',res=450,background='white');draw();dev.off()
