# Scoped display-only revision from frozen evidence; no upstream fits or writes.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(grid);library(ragg);library(png);library(jsonlite);library(digest)})
source('8-Visualization/scripts/fig_style.R')
out <- '8-Visualization/refinements/fig5_integrated_20261101_v11'
old <- '8-Visualization/refinements/fig5_three_rows_20261101_v10'
for(d in c('figures','supplement','source_data','review','qa','provenance'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
ink <- '#293F50'; mid <- '#758B9B'; pal <- load_subtype_palette()$colors[c('C1','C2','C3')];invisible(load_study_palette())
src <- function(n)fread(file.path(old,'source_data',paste0(n,'.csv')))
mins <- '5-ClassifierValidation/runs/bulk-rnaseq__20261031_step5_minpanel_v1/tables'
mi <- function(n)fread(file.path(mins,paste0(n,'.csv')))
manifest <- list(); source_index <- list()
tx <- function(s,x,y,size=6,bold=FALSE,just='left',col=ink,italic=FALSE)grid.text(s,unit(x,'native'),unit(y,'native'),just=just,gp=gpar(fontfamily='Arial',fontsize=size,fontface=if(italic)'italic' else if(bold)'bold' else 'plain',col=col,lineheight=1.05))
rc <- function(x,y,w,h,col,edge=NA)grid.rect(unit(x,'native'),unit(y,'native'),unit(w,'native'),unit(h,'native'),gp=gpar(fill=col,col=edge,lwd=.35))
ln <- function(x,y,col=ink,lwd=.5,lty=1,arrow=FALSE)grid.lines(unit(x,'native'),unit(y,'native'),gp=gpar(col=col,lwd=lwd,lty=lty),arrow=if(arrow)grid::arrow(length=unit(1,'mm'),type='closed')else NULL)
pt <- function(x,y,col=ink,size=.8,pch=16)grid.points(unit(x,'native'),unit(y,'native'),pch=pch,size=unit(size,'mm'),gp=gpar(col=col))
canvas <- function(w,h){grid.newpage();pushViewport(viewport(xscale=c(0,w),yscale=c(0,h)))}
save <- function(id,stem,draw,w,h,track='Main'){
 directory<-if(track=='Main')'figures' else 'supplement';p<-file.path(out,directory,paste0(stem,'.pdf'));im<-file.path(out,'review',paste0(stem,'.png'))
 grDevices::cairo_pdf(p,width=w/25.4,height=h/25.4,family='Arial');draw();dev.off()
 ragg::agg_png(im,width=w,height=h,units='mm',res=450,background='white');draw();dev.off()
 manifest[[id]]<<-data.table(Panel=id,Track=track,PDF=p,PNG=im,WidthMM=w,HeightMM=h,Reuse='Redrawn from fixed source data')
}
reused<-fread(file.path(old,'presentation_manifest.csv'))
for(i in seq_len(nrow(reused))){
 z<-reused[i];id<-sub('Figure 5 ','',z$DisplayLabel,fixed=TRUE);id<-sub('Supplementary ','S',id,fixed=TRUE)
 if(!id%in%c('B','C','F',paste0('S',LETTERS[1:8])))next
 stem<-sub('[.]pdf$','',basename(z$PDF));track<-if(startsWith(id,'S'))'Supplement' else 'Main';directory<-if(track=='Main')'figures' else 'supplement'
 target<-file.path(out,directory,paste0(stem,'.pdf'));im<-file.path(out,'review',paste0(stem,'.png'))
 stopifnot(file.copy(z$PDF,target,overwrite=TRUE),file.copy(sub('[.]pdf$','.png',z$PDF),im,overwrite=TRUE))
 manifest[[id]]<-data.table(Panel=id,Track=track,PDF=target,PNG=im,WidthMM=z$WidthMM,HeightMM=z$HeightMM,Reuse='Exact PDF and R-PNG reuse')
}
drawA <- function(){canvas(174,32)
 starts<-c(1,31,76,107,139);ww<-c(27,42,28,29,34);heads<-c('Data','Selection within training data','Consensus','Classifier','Evaluation')
 for(j in 1:5){x<-starts[j];w<-ww[j];rc(x+w/2,17.8,w,24.5,'#F6F8FA','#DCE4EA');rc(x+w/2,27.5,w,4.5,'#DCE5EC');tx(heads[j],x+w/2,27.5,6.6,TRUE,'centre');if(j<5)ln(c(x+w+.3,starts[j+1]-.5),c(17.8,17.8),arrow=TRUE)}
 tx('497 principal tumors',14.5,22.5,6.2,FALSE,'centre');tx('Expression + reference',14.5,18.3,5.8,FALSE,'centre')
 for(j in 1:3){rc(6.5+(j-1)*8,13.3,7.5,2.8,pal[j]);tx(names(pal)[j],6.5+(j-1)*8,9.1,6,TRUE,'centre')}
 tx('limma',33,22,6.2,TRUE);tx('Study-adjusted DE',48,22,5.9);tx('RF',33,16,6.2,TRUE);tx('Permutation importance',48,16,5.9);tx('LASSO',33,10,6.2,TRUE);tx('Multinomial selection',48,10,5.9)
 tx('Supported by ≥2 methods',90,22,5.9,FALSE,'centre');tx('114 genes',90,15.9,7.5,TRUE,'centre');tx('Full-development fit',90,9.5,5.7,FALSE,'centre')
 tx('Linear SVM',121.5,21,7,TRUE,'centre');tx('Primary model',121.5,15.6,6,FALSE,'centre');tx('RF / kNN comparators',121.5,9.5,5.7,FALSE,'centre')
 tx('Repeated nested CV',156,22,6,FALSE,'centre');tx('Study-held-out transfer',156,16,6,FALSE,'centre');tx('Microarray projection',156,10,6,FALSE,'centre')
 tx('Validation reselects genes within training folds · Outer 5-fold × 5 repeats; inner 4-fold',87,2,6.1,FALSE,'centre');popViewport()
}
cm<-src('D_confusion');st<-fread('8-Visualization/refinements/fig5_final_layout_20261031_v9/source_data/5E_study_balanced_accuracy.csv');setorder(st,Study)
drawD <- function(){canvas(60,78);ramp<-circlize::colorRamp2(c(0,1),c('#F0F5F8','#4E8DB1'))
 tx('Internal and study-held-out performance',1,75.2,6.8,TRUE);tx('Internal OOF · n = 497 · BA 98.79%',1,71.1,5.8)
 for(i in 1:3)for(j in 1:3){z<-cm[Truth==paste0('C',i)&Prediction==paste0('C',j)];x<-17+(j-.5)*10;y<-38+(3-i+.5)*10;rc(x,y,9.7,9.7,ramp(z$RowFraction));tx(sprintf('%d\n%.1f%%',z$Freq,100*z$RowFraction),x,y,6.4,FALSE,'centre',if(z$RowFraction>.5)'white'else ink)}
 for(j in 1:3){tx(paste0('C',j),17+(j-.5)*10,35.7,5.8,FALSE,'centre');tx(paste0('C',j),15,38+(3-j+.5)*10,5.8,FALSE,'right')}
 grid.text('Reference',unit(7,'native'),unit(53,'native'),rot=90,gp=gpar(fontfamily='Arial',fontsize=5.8,col=ink));tx('Predicted',49,35.7,5.6)
 tx('Held-out studies · transductive centering',1,30.7,6.2,TRUE)
 xp<-function(v)31+(v-.84)/.16*26
 for(v in c(.85,.9,.95,1)){ln(rep(xp(v),2),c(6,28),'#E5EAEE',.4);tx(round(100*v),xp(v),4.2,5.3,FALSE,'centre')}
 for(i in 1:7){z<-st[i];y<-27-(i-1)*3.3;tx(sprintf('%s (%d)',z$Study,z$N),1,y,5.4);ln(xp(c(z$OOF,z$LOSO)),rep(y,2),'#A5B3BD',.7);pt(xp(z$OOF),y,'#879CA9',1.1,1);pt(xp(z$LOSO),y,ink,1.1,16)}
 pt(3,1.5,'#879CA9',1.1,1);tx('Internal',5,1.5,5.3);pt(19,1.5,ink,1.1);tx('Held-out',21,1.5,5.3);tx('BA (%)',57,1.5,5.3,FALSE,'right');popViewport()
}
drawE<-function(){canvas(60,60);s<-src('E_summary');d<-src('E_composition');short<-c('Andersen','Sia','Jusakual P1','Jusakual P2','Job','TIGER')
 tx('Microarray projection (n = 506)',1,57.7,7,TRUE);tx('Study (n) · genes / 114',1,53.4,5.7);tx('Prob. gap',51,53.4,5.5,FALSE,'centre',mid)
 for(i in 1:6){yy<-49.5-(i-1)*7.3;z<-s[CohortOrder==i];q<-d[CohortOrder==i];setorder(q,Prediction);tx(sprintf('%s (%d) · %d/114',short[i],z$N,z$MatchedFeatures),1,yy,5.7)
  cum<-0;for(j in 1:3){a<-q[Prediction==paste0('C',j)];wid<-36*a$Fraction;rc(1+cum+wid/2,yy-3,wid,2.7,pal[j]);tx(a$N,1+cum+wid/2,yy-3,5.5,FALSE,'centre',if(j==2)ink else 'white');cum<-cum+wid}
  ln(c(44,58),rep(yy-3,2),'#E4E9ED',.4);ln(44+14*c(z$Q1,z$Q3),rep(yy-3,2),mid,1.1);pt(44+14*z$MedianGap,yy-3,mid,.85)
 }
 tx('0',1,7.6,5.5);tx('100%',37,7.6,5.5,FALSE,'right');tx('0',44,7.6,5.5,FALSE,'centre',mid);tx('1',58,7.6,5.5,FALSE,'centre',mid)
 for(j in 1:3){xx<-3+(j-1)*10;rc(xx,3.3,1.8,1.8,pal[j]);tx(paste0('C',j),xx+1.7,3.3,5.5)};tx('Median [IQR]',51,3.3,5.3,FALSE,'centre',mid);popViewport()
}
w<-mi('full_candidate_weights');delta<-mi('paired_differences')[Metric=='BalancedAccuracy'];perf<-mi('performance_ci');choice<-mi('choice')
drawG<-function(){canvas(53,106);tx('17-gene development candidate',1,103,7,TRUE);tx('Model weights; fixed panel not validated',1,98.7,5.7)
 xp<-function(v)16+v/12*35
 for(i in 1:17){y<-94-(i-1)*3.1;tx(w$Gene[i],14.5,y,5.8,FALSE,'right',italic=TRUE);ln(c(xp(0),xp(100*w$RelativeWeight[i])),rep(y,2),'#C6D1D9',.75);pt(xp(100*w$RelativeWeight[i]),y,ink,1.25)}
 ln(c(16,51),c(41,41),'#7D8A92',.4);for(v in c(0,5,10)){ln(rep(xp(v),2),c(40.3,41),'#7D8A92',.4);tx(v,xp(v),38.6,5.5,FALSE,'centre')};tx('Relative squared weight (%)',33.5,35,5.8,FALSE,'centre')
 tx('Adaptive compact − full workflow',1,29.5,6.5,TRUE);tx(sprintf('Δ BA %.2f pp [%.2f, %.2f]',100*delta$Estimate,100*delta$Lower95,100*delta$Upper95),1,25.3,5.9)
 tx(sprintf('Fold-selected K: %g [%g–%g]',median(choice$K),min(choice$K),max(choice$K)),1,21.1,5.7)
 xx<-function(v)5+(v+2.6)/2.9*45
 ln(rep(xx(0),2),c(9,18),'#9AA8B2',.5);ln(rep(xx(-1),2),c(9,18),'#9AA8B2',.5,2);tx('−1 pp target',xx(-1),18.7,5.3,FALSE,'centre',mid)
 ln(xx(100*c(delta$Lower95,delta$Upper95)),rep(13.2,2),ink,1.2);for(v in 100*c(delta$Lower95,delta$Upper95))ln(rep(xx(v),2),c(12.2,14.2),ink,.65);pt(xx(100*delta$Estimate),13.2,ink,1.8)
 ln(xx(c(-2.6,.3)),c(8.5,8.5),'#7D8A92',.4);for(v in c(-2,-1,0)){tx(v,xx(v),6.5,5.4,FALSE,'centre')};tx('BA-loss target not met',26.5,2,5.9,TRUE,'centre');popViewport()
}
save('A','Figure5A_workflow',drawA,174,32);save('D','Figure5D_internal_and_study_transfer',drawD,60,78);save('E','Figure5E_microarray_projection',drawE,60,60);save('G','Figure5G_compact_tradeoff',drawG,53,106)
drawSI<-function(){d<-copy(perf[Metric%in%c('BalancedAccuracy','C1Recall','C2Recall','C3Recall')]);d[,Y:=match(Metric,c('C3Recall','C2Recall','C1Recall','BalancedAccuracy'))+ifelse(Model=='Full',.12,-.12)]
 p<-ggplot(d,aes(100*Estimate,Y,shape=Model,colour=Model))+geom_errorbar(aes(xmin=100*Lower95,xmax=100*Upper95),orientation='y',width=.12,linewidth=.4)+geom_point(size=1.8)+scale_y_continuous(breaks=1:4,labels=c('C3 recall','C2 recall','C1 recall','Balanced accuracy'),limits=c(.5,4.5))+scale_x_continuous(breaks=c(94,96,98,100),limits=c(93,100.4))+scale_colour_manual(values=c(Full='#8197A6',Compact=ink),labels=c(Full='Full workflow',Compact='Adaptive compact'))+scale_shape_manual(values=c(Full=1,Compact=16),labels=c(Full='Full workflow',Compact='Adaptive compact'))+labs(title='Adaptive compact and full workflow',subtitle='Same 497 patients; repeated nested CV',x='OOF performance (%)',y=NULL,shape=NULL,colour=NULL,caption='95% intervals conditional on saved OOF predictions.\nThese predictions do not validate a fixed 17-gene panel.')+theme_classic(base_size=7,base_family='Arial')+theme(plot.title=element_text(size=8,face='bold'),plot.subtitle=element_text(size=6.5),plot.caption=element_text(size=6,hjust=0),axis.text=element_text(size=6.5),legend.position='bottom',legend.text=element_text(size=6.5),legend.key.width=unit(3,'mm'),axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=margin(2,2,2,2,'mm'));print(p)}
save('SI','SupplementaryFigure5I_compact_absolute_performance',drawSI,90,68,'Supplement')
for(p in list.files(file.path(old,'source_data'),pattern='[.]csv$',full.names=TRUE)){if(basename(p)!='layout_positions_mm.csv')file.copy(p,file.path(out,'source_data',basename(p)),overwrite=TRUE)}
fwrite(st,file.path(out,'source_data/D_study_transfer.csv'));fwrite(choice,file.path(out,'source_data/G_training_fold_choices.csv'));fwrite(perf,file.path(out,'source_data/SI_absolute_performance.csv'))
mf<-rbindlist(manifest);mf[,Order:=match(Panel,c(LETTERS[1:7],paste0('S',LETTERS[1:9])))];setorder(mf,Order);mf[,Order:=NULL];fwrite(mf,file.path(out,'figure_manifest.csv'))
fwrite(reused,file.path(out,'provenance/previous_panel_sources.csv'))
# Hash every input indexed by the preceding source manifest, including non-CSV owners.
input_paths<-unique(c(unlist(strsplit(reused$SourceTable,';',fixed=TRUE)),list.files(file.path(old,'source_data'),full.names=TRUE),reused$PDF,sub('[.]pdf$','.png',reused$PDF),file.path(mins,'choice.csv'),'8-Visualization/refinements/fig5_final_layout_20261031_v9/source_data/5E_study_balanced_accuracy.csv'))
stopifnot(all(file.exists(input_paths)));hashes<-data.table(Path=input_paths,SHA256=vapply(input_paths,function(p)digest(p,algo='sha256',file=TRUE),character(1)));fwrite(hashes,file.path(out,'provenance/input_sha256.csv'))
for(i in seq_len(nrow(reused))){p<-unique(strsplit(reused$SourceTable[i],';',fixed=TRUE)[[1]]);for(f in p[grepl('[.]csv$',p)]){name<-paste0(substr(digest(f,algo='sha256'),1,10),'_',basename(f));dest<-file.path(out,'source_data',name);stopifnot(file.copy(f,dest,overwrite=TRUE));source_index[[length(source_index)+1]]<-data.table(Panel=reused$DisplayLabel[i],Source=f,Packaged=dest)}}
fwrite(rbindlist(source_index),file.path(out,'source_data/source_index.csv'))
# R-only image placement; layout preview is not a publication master.
ragg::agg_png(file.path(out,'review/Figure5_three_row_preview.png'),width=182,height=270,units='mm',res=300,background='white');canvas(182,270)
place<-data.table(Panel=c('A','B','C','D','E','F','G'),X=c(4,4,80,4,4,68,125),Top=c(3,41,41,123,205,123,123))
for(i in 1:nrow(place)){z<-mf[Panel==place$Panel[i]];im<-png::readPNG(z$PNG);grid.raster(im,x=unit(place$X[i],'native'),y=unit(270-place$Top[i],'native'),width=unit(z$WidthMM,'native'),height=unit(z$HeightMM,'native'),just=c('left','top'));tx(place$Panel[i],place$X[i]-3,270-place$Top[i]-2,8,TRUE)};popViewport();dev.off();fwrite(place,file.path(out,'source_data/layout_positions_mm.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'provenance/sessionInfo.txt'))
file.copy('8-Visualization/scripts/render_fig5_integrated_20261101_v11.R',file.path(out,'provenance/render_fig5_integrated_20261101_v11.R'),overwrite=TRUE)
cat('Figure5 integrated suite: ',nrow(mf),' panels\n')
