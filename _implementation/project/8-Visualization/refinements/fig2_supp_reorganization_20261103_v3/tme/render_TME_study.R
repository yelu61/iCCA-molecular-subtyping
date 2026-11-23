# Display-only compact principal-state support. Frozen values and populations.
# R grid redraw: physical coordinates are recomposed; no uniform panel scaling.
suppressPackageStartupMessages({library(data.table);library(grid)})
source('8-Visualization/scripts/fig_style.R')
out <- '8-Visualization/refinements/fig2_supp_reorganization_20261103_v3/tme'
for(d in c('figures','tables','review'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
src <- c(
 '8-Visualization/refinements/sfig2_wgcna_narrow_20261031_v37/tables/module_order.csv',
 '8-Visualization/refinements/sfig2_wgcna_narrow_20261031_v37/tables/network_principal_effects.csv',
 '8-Visualization/refinements/sfig2_compact_20261031_v35/tables/TME_primary_LOSO_sensitivity.csv',
 '8-Visualization/refinements/sfig2_compact_20261031_v35/tables/TME_per_study.csv',
 '8-Visualization/refinements/sfig2_compact_20261031_v35/tables/TME_LOSO_all.csv',
 'governance/data/subtype_label_color_dictionary.csv','0-Data/qc/study_color_dictionary.csv')
hashes <- tools::md5sum(src)
for(p in src[1:5])stopifnot(file.copy(p,file.path(out,'tables',basename(p)),overwrite=TRUE))
cl<-fread(src[1]);me<-fread(src[2]);tme<-fread(src[3]);per<-fread(src[4]);lo<-fread(src[5])
stopifnot(nrow(cl)==24,nrow(me)==72,sum(cl$MainFigure)==11,nrow(tme)==99,nrow(per)==594)
stopifnot(uniqueN(tme$RowID)==33,all(tme$EstimableOmissions>0))
pal<-load_subtype_palette()$colors;study_pal<-load_study_palette();states<-c('C1','C2','C3')
groups<-c('C1_specific','C2_specific','C3_specific','SharedHighCompatible_C1_C2','SharedHighCompatible_C2_C3','Replicated_gradient_or_ambiguous','Non_specific')
group_labels<-c('C1-specific','C2-specific','C3-specific','Shared-high C1 / C2','Shared-high C2 / C3','Gradient / ambiguous','Non-specific')
mm<-function(x)unit(x,'mm');texts<-list();marks<-list();manifest<-list()
txt<-function(s,x,y,size=6,face='plain',just='left',col='#24282B'){
 grid.text(s,mm(x),mm(y),just=just,gp=gpar(fontfamily='Arial',fontsize=size,fontface=face,col=col))
 texts[[length(texts)+1L]]<<-data.table(Text=as.character(s),X=x,Y=y,Size=size,Face=face,Just=just)
}
rect<-function(x,y,w,h,fill,border=NA)grid.rect(mm(x),mm(y),mm(w),mm(h),gp=gpar(fill=fill,col=border,lwd=.4))
line<-function(x,y,col='#D3D8DC',lwd=.5)grid.lines(mm(x),mm(y),gp=gpar(col=col,lwd=lwd))
point<-function(x,y,pch=16,size=1,col='#303840',fill=col)grid.points(mm(x),mm(y),pch=pch,size=mm(size),gp=gpar(col=col,fill=fill,lwd=.55))
ramp<-colorRamp(c(SCALE_EFFECT_LOW,'#FAFAFA',SCALE_EFFECT_HIGH))
colval<-function(v,lim=2)rgb(ramp(pmax(0,pmin(1,(v/lim+1)/2))),maxColorValue=255)
emit<-function(name,draw,w,h,panel){
 f<-file.path(out,'figures',name)
 grDevices::cairo_pdf(paste0(f,'.pdf'),width=w/25.4,height=h/25.4,family='Arial');draw();dev.off()
 fwrite(rbindlist(texts,fill=TRUE),file.path(out,'tables',paste0(name,'_text.csv')))
 fwrite(rbindlist(marks,fill=TRUE),file.path(out,'tables',paste0(name,'_marks.csv')))
 svglite::svglite(paste0(f,'.svg'),width=w/25.4,height=h/25.4);draw();dev.off()
 ragg::agg_png(paste0(f,'.png'),width=w,height=h,units='mm',res=600,background='white');draw();dev.off()
 ragg::agg_png(file.path(out,'review',paste0(name,'_preview.png')),width=w,height=h,units='mm',res=300,background='white');draw();dev.off()
 manifest[[length(manifest)+1L]]<<-data.table(Figure='Supplementary Figure 2',Panel=panel,TargetOutput=paste0(f,'.pdf'),WidthMM=w,HeightMM=h,Status='Scientific-content draft review',SourceTable=paste(src[1:5],collapse=';'))
}
draw_d<-function(){
 grid.newpage();texts<<-list();marks<<-list();W<-174;H<-121
 txt('Study-specific effects and study-omission consistency',2,118.5,7,'bold')
 studies<-c('2019Ahn','2022Dong','2023Jeon','2023Lee','2024Fan','TCGA-CHOL');abbr<-c('A','D','J','L','F','T')
 starts<-c(54,93,132);step<-3.9
 sel<-unique(tme[,.(RowID,Method,Score,Label,Block,N,Studies)]);sel[,PlotGroup:=ifelse(Block=='Context','Context',Method)];setorder(sel,RowID)
 txt('Readout',2,114,6.3,'bold');txt('n / studies',48,114,6,just='right')
 for(j in 1:3){s<-starts[j];txt(states[j],s+14.5,114,6.6,'bold','center');rect(s+14.5,111.7,35.5,.8,pal[states[j]])
  for(k in 1:6)txt(abbr[k],s+(k-1)*step,108.8,6,just='center')
  txt('LOSO',s+29,108.8,6,just='center')
 }
 y<-105.5
 for(g in c('Context','MCP-counter','EPIC','xCell','CIBERSORT')){
  rect(W/2,y,W-4,2.5,'#EDF1F3')
  txt(if(g=='Context')'ESTIMATE context'else if(g=='CIBERSORT')'CIBERSORT  |  relative immune composition'else g,2,y,6,'bold');y<-y-2.8
  for(id in sel[PlotGroup==g,RowID]){
   r<-sel[RowID==id];if(id%%2==0)rect(W/2,y,W-4,2.45,'#FAFBFC')
   lab<-r$Label
   txt(lab,2,y,6);txt(paste0(r$N,' / ',r$Studies),48,y,6,just='right',col='#65727C')
   for(j in 1:3){s<-starts[j];z<-tme[RowID==id&Subtype==states[j]];stopifnot(nrow(z)==1)
    for(k in 1:6){q<-per[RowID==id&Subtype==states[j]&Study==studies[k]];stopifnot(nrow(q)==1);xx<-s+(k-1)*step
     if(!isTRUE(q$Estimable))txt('\u2013',xx,y,6,just='center',col='#A6AFB7')else{
      cc<-if(q$Effect>=0)'#B65058'else'#477DA6';solid<-q$Lower>0||q$Upper<0
      point(xx,y,pch=if(q$Effect>=0)24 else 25,size=1.15,col=cc,fill=if(solid)cc else 'white')
     }
     marks[[length(marks)+1L]]<<-data.table(MarkType='Study',RowID=id,Subtype=states[j],Study=studies[k],Effect=q$Effect,Lower=q$Lower,Upper=q$Upper,StudyEstimable=q$Estimable,X=xx,Y=y)
    }
    ok<-z$Retained==z$EstimableOmissions
    txt(paste0(z$Retained,'/',z$EstimableOmissions),s+29,y,6,face=if(ok)'plain'else'bold',col=if(ok)'#65727C'else'#24282B',just='center')
    marks[[length(marks)+1L]]<<-data.table(MarkType='LOSO',RowID=id,Subtype=states[j],Retained=z$Retained,EstimableOmissions=z$EstimableOmissions,TotalOmissions=z$TotalOmissions,N=z$N,Studies=z$Studies,X=s+29,Y=y)
   };y<-y-2.5
  };y<-y-.2
 }
 stopifnot(y>7)
 point(2.8,6.8,24,1.15,'#B65058');txt('Higher',4.6,6.8,5.5)
 point(17,6.8,25,1.15,'#477DA6');txt('Lower',18.8,6.8,5.5)
 txt('Filled: 95% CI excludes 0; open: includes 0; dash: not estimable',31,6.8,5.5)
 txt('LOSO: retained / informative omissions',119,6.8,5.5)
 txt('A Ahn   D Dong   J Jeon   L Lee   F Fan   T TCGA',2,3.2,5.5)
 txt('Bold LOSO fraction: at least one direction change',99,3.2,5.5)
}
emit('SFigure2J_TME_study_robustness',draw_d,174,121,'J')
fwrite(rbindlist(manifest),file.path(out,'figure_manifest.csv'))
stopifnot(identical(unname(hashes),unname(tools::md5sum(src))))
fwrite(data.table(Source=src,BeforeMD5=unname(hashes),AfterMD5=unname(tools::md5sum(src))),file.path(out,'review/source_hashes.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'review/sessionInfo.txt'))
cat('TME J:174x121mm; EPIC before xCell to match K; all original study/LOSO values retained.\n')
