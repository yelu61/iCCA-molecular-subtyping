# Frozen-membership composition, not enrichment or a new program definition.
# R-only quantitative grid with all zero/ambiguous categories preserved.
suppressPackageStartupMessages({library(data.table);library(grid)})
source('8-Visualization/scripts/fig_style.R')
out<-'8-Visualization/refinements/fig2_global_revision_20261103_v2/pattern_bridge'
m<-fread(file.path(out,'tables/pattern_module_counts_rowpercent.csv'))
cols<-fread(file.path(out,'tables/module_column_registry.csv'))
reg<-fread(file.path(out,'tables/shared_program_registry.csv'))[Family=='Pattern']
stopifnot(nrow(m)==156,nrow(cols)==26,nrow(reg)==6,sum(m$IntersectionGenes)==842)
study_pal<-load_study_palette()
mm<-function(x)unit(x,'mm');texts<-list();marks<-list()
txt<-function(s,x,y,size=6,face='plain',just='left',col='#26333D',rot=0){
 grid.text(s,mm(x),mm(y),just=just,rot=rot,gp=gpar(fontfamily='Arial',fontsize=size,fontface=face,col=col))
 texts[[length(texts)+1L]]<<-data.table(Text=as.character(s),X=x,Y=y,Size=size,Face=face,Just=just,Rotation=rot)
}
rect<-function(x,y,w,h,fill,border=NA)grid.rect(mm(x),mm(y),mm(w),mm(h),gp=gpar(fill=fill,col=border,lwd=.35))
line<-function(x,y,col='#CDD6DC',lwd=.5)grid.lines(mm(x),mm(y),gp=gpar(col=col,lwd=lwd))
ramp<-colorRamp(c('#FAFBFC','#BFD6E2','#5794B1','#16506D'))
colval<-function(p)rgb(ramp(pmax(0,pmin(1,p/80))),maxColorValue=255)
draw<-function(){
 grid.newpage();texts<<-list();marks<<-list();W<-174;H<-123;left<-49;cw<-4.65;rh<-9
 txt('Pattern–WGCNA membership',2.5,119,8,'bold')
 txt('All 842 pattern genes; each cell shows gene count (top) and within-pattern % (bottom)',2.5,114.5,6.1)
 for(i in 1:26){c<-cols[ColumnOrder==i];x<-left+(i-.5)*cw
  lab<-if(c$Module=='grey')'grey'else if(c$Module=='Not_in_network')'Outside'else c$ModuleLabel
  txt(lab,x,97.1,5.8,face=if(c$InCurrentFigure3)'bold'else'plain',just='left',rot=90)
  rect(x,94.6,cw-.2,1.35,c$ModuleColor)
 }
 for(j in 1:6){p<-paste0('P',j);r<-reg[StableID==p];y<-89.2-(j-1)*rh
  txt(p,2.5,y+1.8,6.8,'bold');txt(r$CompactDescription,8.5,y+1.8,6.2)
  txt(sprintf('n = %s; low affinity = %s',r$NGenes,r$LowAffinityGenes),2.5,y-1.8,5.8,col='#66737D')
  for(i in 1:26){q<-m[Pattern==p&ColumnOrder==i];x<-left+(i-.5)*cw;v<-q$RowPercent
   rect(x,y,cw-.17,rh-.35,colval(v),border='#E4EAEE')
   cc<-if(v>=40)'white'else if(q$IntersectionGenes==0)'#A9B0B5'else'#253744'
   txt(q$IntersectionGenes,x,y+1.5,6.1,just='center',col=cc)
   txt(sprintf('%.1f',v),x,y-1.5,5.5,just='center',col=cc)
   marks[[length(marks)+1L]]<<-data.table(Pattern=p,Module=q$Module,IntersectionGenes=q$IntersectionGenes,RowPercent=v,PatternN=q$PatternN,X=x,Y=y)
  }
 }
 line(rep(left+24*cw-.03,2),c(39.5,95.5),'#8E9BA4',.7)
 line(c(2.5,170),rep(38.3,2))
 txt('Genes in the 842-gene display',2.5,34.4,6)
 txt('Complete module genes',2.5,29.5,6)
 for(i in 1:26){c<-cols[ColumnOrder==i];x<-left+(i-.5)*cw
  txt(c$GenesIn842,x,34.4,5.5,just='center')
  txt(if(is.na(c$FullModuleN))'–'else c$FullModuleN,x,29.5,5.5,just='center')
 }
 txt('Within-pattern %',2.5,21.3,6.1)
 for(i in 1:100)rect(36+(i-.5)*41/100,21.3,41/100+.01,1.9,colval((i-1)/99*80))
 for(v in seq(0,80,20))txt(v,36+v/80*41,17.8,5.5,just='center')
 txt('Bold module labels: current Fig. 2C / Fig. 3 programs',85,22.1,5.8)
 txt('Descriptions are working bulk annotations; P IDs are unchanged.',85,18.2,5.5)
 txt('Grey is unassigned; outside-network genes are retained. Complete-module totals use a different denominator.',2.5,10.8,6)
 txt('Pattern rows partition 842 genes; the 24 named modules contain 6,753 genes. No enrichment test was performed.',2.5,6.8,6)
}
base<-file.path(out,'figures','Pattern_WGCNA_membership_bridge')
grDevices::cairo_pdf(paste0(base,'.pdf'),width=174/25.4,height=123/25.4,family='Arial');draw();dev.off()
fwrite(rbindlist(texts),file.path(out,'tables/bridge_text.csv'));fwrite(rbindlist(marks),file.path(out,'tables/bridge_marks.csv'))
svglite::svglite(paste0(base,'.svg'),width=174/25.4,height=123/25.4);draw();dev.off()
ragg::agg_png(paste0(base,'.png'),width=174,height=123,units='mm',res=600,background='white');draw();dev.off()
ragg::agg_png(file.path(out,'review','Pattern_WGCNA_membership_bridge_preview.png'),width=174,height=123,units='mm',res=300,background='white');draw();dev.off()
fwrite(data.table(Panel='Pattern–WGCNA descriptive membership bridge',PDF=paste0(base,'.pdf'),SVG=paste0(base,'.svg'),PNG=paste0(base,'.png'),WidthMM=174,HeightMM=123,SourceTable=file.path(out,'tables/pattern_module_counts_rowpercent.csv'),Status='Scientific-content draft review'),file.path(out,'figure_manifest.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'review/sessionInfo.txt'))
cat('All 156 cells exported with counts and within-pattern percentages.\n')
