suppressPackageStartupMessages({library(ComplexHeatmap);library(circlize);library(grid)})
o<-normalizePath('8-Visualization/refinements/sfig4_cnv_D_20261104_v1');tt<-file.path(o,'source_data');ff<-file.path(o,'figures')
studycols<-read.csv('0-Data/qc/study_color_dictionary.csv')
d<-read.csv(file.path(tt,'ICC20_rows.csv'));m<-as.matrix(read.csv(gzfile(file.path(tt,'ICC20_cnv.csv.gz')),header=FALSE));w<-read.csv(file.path(tt,'ICC20_windows.csv'));support<-read.csv(file.path(tt,'patient_support.csv'))
d$StrictCandidate<-tolower(as.character(d$StrictCandidate))=='true'
stopifnot(nrow(d)==2265,nrow(m)==2265,ncol(m)==nrow(w),sum(d$StrictCandidate)==27)
metrics_device<-file.path(o,'review/font_metrics.pdf');cairo_pdf(metrics_device,width=174/25.4,height=61/25.4,family='Arial')
chr<-factor(sub('chr','',w$Chromosome),levels=as.character(1:22));g<-unique(d$Group);lab<-setNames(c('Ref. controls (220)','Tumor epi. (2,017)','Mixed / review (28)'),g);split<-factor(lab[d$Group],levels=lab[g]);pal<-colorRamp2(c(-.3,0,.3),c('#3D6E9B','#FAFAFA','#B54D4B'));gp<-gpar(fontfamily='Arial',fontsize=5.2)
flag<-rowAnnotation(Flag=factor(ifelse(d$StrictCandidate,'Strict','Other / ref.'),levels=c('Other / ref.','Strict')),col=list(Flag=c('Other / ref.'='#F2F2F2',Strict='#292929')),simple_anno_size=unit(1.6,'mm'),show_annotation_name=FALSE,annotation_legend_param=list(Flag=list(title='Strict flag',title_gp=gp,labels_gp=gp,direction='horizontal')))
ht<-Heatmap(m,name='RNA-CNV',col=pal,cluster_rows=FALSE,cluster_columns=FALSE,row_split=split,cluster_row_slices=FALSE,column_split=chr,cluster_column_slices=FALSE,show_row_names=FALSE,show_column_names=FALSE,row_title_rot=0,row_title_gp=gp,column_title=ifelse(1:22 %in% c(1,3,5,7,9,11,13,15,17,19,22),as.character(1:22),''),column_title_rot=90,column_title_gp=gpar(fontfamily='Arial',fontsize=5),row_gap=unit(1,'mm'),column_gap=unit(.15,'mm'),use_raster=TRUE,raster_device='png',raster_quality=3,raster_resize_mat=FALSE,raster_by_magick=FALSE,width=unit(77,'mm'),height=unit(32,'mm'),heatmap_legend_param=list(title_gp=gp,labels_gp=gp,at=c(-.3,0,.3),direction='horizontal',legend_width=unit(25,'mm')))+flag
W<-174;H<-61
textmm<-function(t,x,y,size=5.5,bold=FALSE,just='left',col='#202020'){grid.text(t,x=unit(x,'mm'),y=unit(H-y,'mm'),just=just,gp=gpar(fontfamily='Arial',fontsize=size,fontface=if(bold)'bold'else'plain',col=col))}
drawpanel<-function(){grid.newpage();grid.rect(gp=gpar(fill='white',col=NA));textmm('Representative RNA-CNV evidence and patient support',5,3.5,7,TRUE);
 pushViewport(viewport(x=unit(0,'mm'),y=unit(H-11,'mm'),width=unit(110,'mm'),height=unit(48,'mm'),just=c('left','top')));draw(ht,newpage=FALSE,heatmap_legend_side='bottom',annotation_legend_side='bottom',merge_legends=TRUE,padding=unit(c(.5,1,.5,.5),'mm'),legend_grouping='original');popViewport()
 textmm('ICC20: complete saved profiles and strict-candidate flag',4,8,5.5,TRUE)
 textmm('Current strict-screen support',114,9,5.5,TRUE)
 xs<-c(114,151,171);textmm('Patient / region',xs[1],15,5.2,TRUE);textmm('Clean epi.',xs[2],15,5.2,TRUE,just='right');textmm('Strict',xs[3],15,5.2,TRUE,just='right');grid.lines(x=unit(c(114,172),'mm'),y=unit(rep(H-17,2),'mm'),gp=gpar(col='#B8B8B8',lwd=.5))
 for(i in 1:nrow(support)){yy<-20+(i-1)*4.2;textmm(support$Region[i],xs[1],yy,5.4);textmm(format(support$CleanEpithelial[i],big.mark=','),xs[2],yy,5.4,just='right');textmm(support$StrictCandidates[i],xs[3],yy,5.4,just='right')}
 textmm('42 candidates / 3 patients',114,39,5.5,TRUE)
 textmm('Only ICC20 has at least 20 strict candidates.',114,44,5)
 textmm('ICC24 regions belong to one patient.',114,48,5)
 textmm('Saved controls are a display subset;\nfitting used all eligible same-patient controls.',114,53,5)
 textmm('RNA-derived candidates; no DNA confirmation.',114,59,5)
}
for(ext in c('pdf','svg','png','tiff')){
 p<-file.path(ff,paste0('SFigure4D_RNA_CNV_patient_support.',ext));
 if(ext=='pdf'){tmp<-tempfile(tmpdir=ff,fileext='.pdf');cairo_pdf(tmp,width=W/25.4,height=H/25.4,family='Arial')}
 if(ext=='svg')svglite::svglite(p,width=W/25.4,height=H/25.4)
 if(ext=='png')ragg::agg_png(p,width=W,height=H,units='mm',res=600,background='white')
 if(ext=='tiff')ragg::agg_tiff(p,width=W,height=H,units='mm',res=600,compression='lzw',background='white')
 drawpanel();dev.off();if(ext=='pdf')stopifnot(file.rename(tmp,p))
}
dev.off();unlink(metrics_device)
writeLines(capture.output(sessionInfo()),file.path(o,'review/R_session_info.txt'))

writeLines(capture.output(warnings()),file.path(o,'review/R_warnings.txt'))
