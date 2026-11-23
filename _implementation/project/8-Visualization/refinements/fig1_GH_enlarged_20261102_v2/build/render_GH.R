# Display-only enlargement within the original Figure1 bottom row.
# Reuse the locked KM fit and RMST table; no model fitting or new comparisons.
library(grid)
library(survival)
library(ragg)
library(svglite)
library(jsonlite)
out<-normalizePath('8-Visualization/refinements/fig1_GH_enlarged_20261102_v2')
owner<-'6-ClinicalSurvival/review/20261026_figure1_h_options'
fit<-readRDS(file.path(owner,'km_fit_loglog.rds'))
d<-read.csv(file.path(owner,'tables/survival_population.csv'))
rmst<-read.csv(file.path(owner,'tables/rmst_24_differences.csv'))
lr<-read.csv('6-ClinicalSurvival/publication/tables/step6_km_logrank_tests.csv')
pal<-read.csv('governance/data/subtype_label_color_dictionary.csv')
pal<-setNames(pal$Color,pal$FinalSubtype)
study_colors<-read.csv('0-Data/qc/study_color_dictionary.csv')
stopifnot(nrow(study_colors)==12,nrow(d)==331,sum(fit$n.event)==153,identical(as.integer(fit$n),c(67L,165L,99L)))
stopifnot(nrow(rmst)==3,all(rmst$TauMonths==24))
times<-seq(0,96,12)
risk<-summary(fit,times=times,extend=TRUE)
rt<-data.frame(Subtype=sub('FinalSubtype=','',risk$strata),Month=risk$time,AtRisk=risk$n.risk)
expected<-rbind(c(67,57,35,17,8,2,0,0,0),c(165,148,99,60,18,2,1,1,0),c(99,59,23,15,5,1,1,0,0))
for(j in 1:3)stopifnot(identical(as.numeric(rt$AtRisk[rt$Subtype==paste0('C',j)]),expected[j,]))
write.csv(rt,file.path(out,'source_data/G_risk_counts.csv'),row.names=FALSE)
file.copy(file.path(owner,'tables/rmst_24_differences.csv'),file.path(out,'source_data/H_RMST_differences.csv'),overwrite=TRUE)
curve<-data.frame(Subtype=rep(sub('FinalSubtype=','',names(fit$strata)),fit$strata),Month=fit$time,Survival=fit$surv,Censored=fit$n.censor)
write.csv(curve,file.path(out,'source_data/G_KM_coordinates.csv'),row.names=FALSE)
tx<-function(t,x,y,size=5.8,bold=FALSE,just='centre',col='#222222',rot=0)
 grid.text(t,x,y,just=just,rot=rot,gp=gpar(fontfamily='Arial',fontsize=size,fontface=if(bold)'bold' else 'plain',col=col,lineheight=1.05))
ln<-function(x,y,col='#222222',lwd=.6,lty=1)grid.lines(x,y,gp=gpar(col=col,lwd=lwd,lty=lty))
drawG<-function(){
 grid.newpage()
 tx('Overall survival across\nprincipal molecular states',.02,.99,6.8,TRUE,just=c('left','top'))
 tx('Available-case KM; n = 331; events = 153',.02,.882,5.1,just='left')
 tx('Study-stratified log-rank P < 0.001',.02,.851,5.1,just='left')
 for(j in 1:3){x<-.3+(j-1)*.23;ln(c(x-.04,x+.005),c(.807,.807),pal[paste0('C',j)],1);tx(paste0('C',j),x+.025,.807,5.6,just='left')}
 xx<-function(t).155+t/96*.785; yy<-function(s).34+s*.435
 ln(c(.155,.94),c(.34,.34));ln(c(.155,.155),c(.34,.775))
 for(v in c(0,.25,.5,.75,1)){ln(c(.145,.155),rep(yy(v),2));tx(sprintf('%.2f',v),.128,yy(v),5.2,just='right')}
 for(t in times){ln(rep(xx(t),2),c(.34,.329));tx(as.character(t),xx(t),.307,5.2)}
 tx('Survival probability',.025,.56,5.8,rot=90)
 tx('Overall survival (months)',.57,.263,5.8)
 for(c in paste0('C',1:3)){
  q<-curve[curve$Subtype==c,];prev<-c(1,head(q$Survival,-1))
  ln(xx(c(0,rep(q$Month,each=2))),yy(c(1,as.vector(rbind(prev,q$Survival)))),pal[c],.8)
  cen<-q[q$Censored>0,];for(k in seq_len(nrow(cen)))ln(rep(xx(cen$Month[k]),2),yy(cen$Survival[k])+c(-.006,.006),pal[c],.6)
 }
 tx('Number at risk',.02,.213,5.8,just='left')
 for(j in 1:3){y<-.153-(j-1)*.061;c<-paste0('C',j);tx(c,.035,y,5.5,just='left',col=pal[c]);for(k in seq_along(times))tx(expected[j,k],xx(times[k]),y,5.5)}
}
drawH<-function(){
 grid.newpage()
 tx('Restricted mean survival:\n0–24 months',.02,.99,6.8,TRUE,just=c('left','top'))
 tx('Differences and 95% CI',.02,.88,5.3,just='left')
 tx('n = 331; events = 153',.02,.847,5.3,just='left')
 xx<-function(v).24+(v+12)/21*.735
 ln(rep(xx(0),2),c(.145,.785),'#AAAAAA',.6,2)
 ys<-c(.665,.445,.225)
 for(i in 1:3){r<-rmst[i,];y<-ys[i]
  tx(gsub(' - ',' − ',r$Contrast),.025,y,5.6,just='left')
  grid.rect(.61,y+.077,.67,.09,gp=gpar(fill='white',col=NA))
  tx(sprintf('%+.2f',r$DifferenceMonths),.61,y+.095,6.2,TRUE)
  tx(sprintf('[%.2f, %.2f]',r$Lower95,r$Upper95),.61,y+.054,5.8)
  ln(xx(c(r$Lower95,r$Upper95)),c(y,y),'#777777',1)
  grid.circle(xx(r$DifferenceMonths),y,r=unit(.65,'mm'),gp=gpar(fill='#426782',col='white',lwd=.4))
 }
 ln(c(.24,.975),c(.145,.145))
 for(v in c(-10,-5,0,5)){ln(rep(xx(v),2),c(.145,.135));tx(v,xx(v),.109,5.3)}
 tx('RMST difference (months)\n(first − second)',.60,.045,5.4)
}
sizes<-list(G=c(51.45,59.25),H=c(39,59.25))
for(panel in c('G','H')){
 w<-sizes[[panel]][1]/25.4;h<-sizes[[panel]][2]/25.4;draw<-if(panel=='G')drawG else drawH
 cairo_pdf(file.path(out,'pdf',paste0('Figure1',panel,'_enlarged.pdf')),width=w,height=h,family='Arial');draw();dev.off()
 svglite(file.path(out,'svg',paste0('Figure1',panel,'_enlarged.svg')),width=w,height=h);draw();dev.off()
 agg_png(file.path(out,'png',paste0('Figure1',panel,'_enlarged.png')),width=w,height=h,units='in',res=600,background='white');draw();dev.off()
}
write_json(list(locked_KM_fit_reused=TRUE,patients=331,events=153,risk_numbers_exact=TRUE,RMST_table_unchanged=TRUE,RMST_contrasts=3,original_KM_x_range=c(0,96),original_RMST_x_range=c(-12,9),no_refitting=TRUE),file.path(out,'review/data_checks.json'),pretty=TRUE,auto_unbox=TRUE)
