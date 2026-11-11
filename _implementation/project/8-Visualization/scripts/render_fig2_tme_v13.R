# Anatomically grouped method comparison; native nonadditive scores stay separate.
if(!exists('STEP8_ROOT'))source('8-Visualization/scripts/fig_style.R')
suppressPackageStartupMessages({library(data.table);library(grid)})
render_tme_v13<-function(){
 out<-'8-Visualization/refinements/fig2_transcriptomic_20261027_v13_tme';run<-'2-TranscriptomicCharacterization/runs/20261027_tme_display_v3'
 for(z in c('figures','tables','review','before'))dir.create(file.path(out,z),recursive=TRUE,showWarnings=FALSE)
 d<-fread(file.path(run,'tables/primary_contrasts.csv'));stopifnot(nrow(d)==447,all(d$Unit=='within-study SD'))
 m<-data.table(Group='Overall context',Label=c('Purity proxy','Immune score','Stromal score'),Method='ESTIMATE',Score=c('PurityProxy','ImmuneScore','StromalScore'),Y=c(22,28,34),Block='Context')
 pop<-data.table(Group=c(rep('Lymphoid',3),rep('Myeloid',2),'Dendritic','Granulocytic','Stromal','Endothelial'),Label=c('CD8 T cells','NK cells','B cells','Monocytes / lineage','Macrophages','Dendritic cells','Neutrophils','Fibroblasts / CAFs','Endothelial cells'),Y=c(46,53,60,70,77,87,97,107,117),M=c('CD8 T cells','NK cells','B lineage','Monocytic lineage',NA,'Myeloid dendritic cells','Neutrophils','Fibroblasts','Endothelial cells'),E=c('CD8_Tcells','NKcells','Bcells',NA,'Macrophages',NA,NA,'CAFs','Endothelial'),X=c('CD8+ T-cells','NK cells','B-cells','Monocytes','Macrophages','cDC','Neutrophils','Fibroblasts','Endothelial cells'),C=c('T cells CD8','Total NK cells','B cells total','Monocytes','Macrophages total','Dendritic total','Neutrophils',NA,NA))
 p<-melt(pop,id.vars=c('Group','Label','Y'),variable.name='Key',value.name='Score',na.rm=TRUE)
 p[,Method:=c(M='MCP-counter',E='EPIC',X='xCell',C='CIBERSORT')[as.character(Key)]];p[,Block:=ifelse(Method=='CIBERSORT','Composition','Native cells')]
 m<-rbindlist(list(m,p),fill=TRUE);stopifnot(nrow(m)==33)
 fwrite(m,file.path(out,'tables/E_selection_map.csv'));z<-merge(m,d,by=c('Method','Score'),all.x=TRUE,allow.cartesian=TRUE);stopifnot(nrow(z)==99,!anyNA(z$N))
 fwrite(z,file.path(out,'tables/E_selected_contrasts.csv'))
 pal<-load_subtype_palette()$colors[c('C1','C2','C3')];invisible(load_study_palette())
 shape<-c('MCP-counter'=21,'EPIC'=22,'xCell'=24,'CIBERSORT'=23,'ESTIMATE'=25)
 offset<-c('MCP-counter'=-1.8,'EPIC'=-.6,'xCell'=.6,'CIBERSORT'=1.8,'ESTIMATE'=0)
 draw<-function(){
  grid.newpage();pushViewport(viewport(xscale=c(0,174),yscale=c(149,0)));on.exit(popViewport())
  txt<-function(s,x,y,size=7.2,face='plain',just='left',col='#252B31')grid.text(s,unit(x,'native'),unit(y,'native'),just=just,gp=gpar(fontfamily='Arial',fontsize=size,fontface=face,col=col))
  rect<-function(x,y,w,h,fill)grid.rect(unit(x,'native'),unit(y,'native'),unit(w,'mm'),unit(h,'mm'),gp=gpar(fill=fill,col=NA))
  line<-function(x,y,col='#D7DCE1',lwd=.6,lty=1)grid.lines(unit(x,'native'),unit(y,'native'),gp=gpar(col=col,lwd=lwd,lty=lty))
  dot<-function(x,y,method,fill,col,size=1.5)grid.points(unit(x,'native'),unit(y,'native'),pch=shape[method],size=unit(size,'mm'),gp=gpar(col=col,fill=fill,lwd=.8))
  txt('Compartment',3,7,7.5,'bold');txt('Readout',36,7,7.5,'bold')
  centers<-c(94,124,154)
  for(j in 1:3){txt(paste0('C',j),centers[j],6.5,9,'bold','center');rect(centers[j],10.4,26,1.1,pal[j]);line(rep(centers[j],2),c(15,122),lty=2)}
  groups<-unique(m$Group)
  for(i in seq_along(groups)){
   g<-groups[i];ys<-m[Group==g,Y];lo<-min(ys)-3.4;hi<-max(ys)+3.4
   if(i%%2==1)rect(87,(lo+hi)/2,170,hi-lo,'#F5F7F9')
   rect(2,(lo+hi)/2,.7,hi-lo,'#B8C3CC');txt(if(g=='Overall context')'Overall\ncontext' else g,4,(lo+hi)/2,7.5,'bold')
  }
  for(j in 1:3)line(rep(centers[j],2),c(15,122),lty=2)
  labels<-unique(m[,.(Label,Y)])
  for(i in seq_len(nrow(labels)))txt(labels$Label[i],36,labels$Y[i],7)
  for(i in seq_len(nrow(z))){r<-z[i];j<-match(r$Subtype,c('C1','C2','C3'));yy<-r$Y+offset[r$Method];xx<-function(v)centers[j]+6*v
   if(!isTRUE(r$Estimable)){txt('–',centers[j],yy,7,just='center');next}
   stopifnot(r$Lower>=-2,r$Upper<=2)
   line(xx(c(r$Lower,r$Upper)),rep(yy,2),adjustcolor(pal[j],alpha.f=.45),.65)
   dot(xx(r$Effect),yy,r$Method,if(r$FDR<.05)pal[j] else 'white',pal[j])
  }
  for(j in 1:3){line(centers[j]+c(-12,12),c(123,123),'#8A9299',.6);for(v in c(-2,0,2)){xx<-centers[j]+6*v;line(c(xx,xx),c(123,124),'#8A9299');txt(v,xx,127,7,just='center')}}
  txt('Study-adjusted difference (within-study SD)',124,132,7.5,just='center')
  methods<-names(shape);lx<-c(5,42,65,92,130)
  for(i in seq_along(methods)){dot(lx[i],139,methods[i],'#444C53','#444C53',1.6);txt(methods[i],lx[i]+3,139,7)}
  dot(5,145,'MCP-counter','#444C53','#444C53');txt('FDR < 0.05',8,145,6.8)
  dot(40,145,'MCP-counter','white','#444C53');txt('FDR ≥ 0.05',43,145,6.8)
  line(c(79,86),c(145,145),'#808991');txt('95% CI',89,145,6.8)
 }
 target<-file.path(out,'figures/Figure2E_TME_compartments.pdf');save_panel(draw,target,174,149)
 png<-sub('[.]pdf$','.png',target);ragg::agg_png(png,width=174,height=149,units='mm',res=300);draw();dev.off();file.copy(png,file.path(out,'review/Figure2E_preview.png'),overwrite=TRUE)
 sources<-c(file.path(run,'tables/primary_contrasts.csv'),file.path(out,'tables/E_selection_map.csv'),'8-Visualization/scripts/render_fig2_tme_v13.R','8-Visualization/scripts/fig_style.R','governance/data/subtype_label_color_dictionary.csv')
 mf<-data.table(Panel='Figure2E_TME_compartments',TargetOutput=target,WidthMM=174,HeightMM=149,SourceTable=paste(sources,collapse=';'),SourceMD5=paste(paste0(sources,'=',unname(tools::md5sum(sources))),collapse=';'),Notes='GSVA removed from main E; ESTIMATE relative purity/immune/stromal context and six anatomical population families. All methods within-study SD with recomputed CIBERSORT contrasts; method shape and FDR fill. Grouping does not average nonadditive scores. Purity is negative ESTIMATE score, not percent.',Status='Rendered',RasterOutput=png)
 fwrite(mf,file.path(out,'tme_manifest.csv'));invisible(mf)
}
if(sys.nframe()==0)render_tme_v13()
