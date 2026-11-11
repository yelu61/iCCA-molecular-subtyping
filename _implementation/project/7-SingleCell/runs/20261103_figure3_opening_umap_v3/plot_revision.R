suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(readr);library(cowplot);library(ggrastr)})
cowplot::set_null_device('cairo')
run<-normalizePath(commandArgs(trailingOnly=TRUE)[1]);root<-normalizePath(file.path(run,'../../..'));figdir<-file.path(run,'figures')
invisible(read_csv(file.path(root,'0-Data/qc/study_color_dictionary.csv'),show_col_types=FALSE))
cells<-read_csv(file.path(run,'source_data/plotted_cells.csv.gz'),show_col_types=FALSE)
paldf<-read_csv(file.path(run,'source_data/annotation_palette.csv'),show_col_types=FALSE);pal<-setNames(paldf$Color,paldf$CellTypeDisplay)
cells$PlotCellType<-factor(cells$PlotCellType,levels=names(pal));cells<-cells%>%add_count(PlotCellType,name='GroupN')%>%arrange(desc(GroupN),CellKey)
xlim<-c(-10.8,22.1);ylim<-c(-13.2,20.0)
pos<-data.frame(PlotCellType=c('Epithelial','Hepatocyte','CD8 T','Other T','NK-like','B cell','Plasma cell','Mono/macrophage','Dendritic-like','Fibroblast','Pericyte-like','Endothelial'),
 x=c(-6.5,-7.0,18.1,15.2,11.4,18.6,14.8,10.0,17.3,2.0,3.7,1.1),
 y=c(7.6,-5.7,.5,-5.7,5.6,8.2,-11.4,18.0,11.7,-5.7,1.1,-10.4),
 target_x=c(-3.7,-6.4,15.5,13.5,11.4,15.6,12,10.2,12.5,4.6,4,4),
 target_y=c(9.5,-4.2,1.4,-3.5,3.8,8.5,-9.2,15.2,10,-3.8,-1,-9))
for(i in seq_len(nrow(pos))){
 z<-cells[cells$PlotCellType==pos$PlotCellType[i],];j<-which.min((z$UMAP1-pos$target_x[i])^2+(z$UMAP2-pos$target_y[i])^2)
 pos$AnchorX[i]<-z$UMAP1[j];pos$AnchorY[i]<-z$UMAP2[j];pos$AnchorCellKey[i]<-z$CellKey[j]
}
base_theme<-theme_classic(base_size=6,base_family='Arial')+theme(
 plot.title=element_text(size=11,face='bold',margin=margin(b=3),colour='#222222'),
 plot.subtitle=element_text(size=5.8,colour='#666666',margin=margin(b=4)),
 axis.title=element_text(size=5.5,colour='#555555'),axis.text=element_text(size=5,colour='#666666'),
 axis.ticks=element_line(linewidth=.20,colour='#888888'),axis.ticks.length=grid::unit(.6,'mm'),
 legend.position='none',panel.grid=element_blank(),plot.margin=margin(3,3,3,3),
 plot.background=element_rect(fill='white',colour=NA))
make_map<-function(s,title,subtitle,main=FALSE,small=FALSE){
 lp<-pos
 if(title=='Adjacent'){lp$x[lp$PlotCellType=='Epithelial']<--6;lp$y[lp$PlotCellType=='Epithelial']<-5.0}
 # Keep the rare tumour hepatocytes explicitly labelled; only three are plotted.
 if(title=='Tumor'){lp$y[lp$PlotCellType=='Hepatocyte']<--6.0}
 lp$Label<-as.character(lp$PlotCellType)
 lp$Label[lp$Label=='Mono/macrophage']<-'Monocytes/macrophages'
 if(main)write_csv(lp,file.path(run,'source_data','main_direct_labels.csv'))
 p<-ggplot(s,aes(UMAP1,UMAP2,colour=PlotCellType))+
  geom_point_rast(data=s[!s$PlotMixedOpenSymbol,],shape=16,size=if(small).19 else .23,alpha=.88,raster.dpi=600)+
  geom_point_rast(data=s[s$PlotMixedOpenSymbol,],shape=1,size=if(small).25 else .32,stroke=.10,alpha=.88,raster.dpi=600)+
  scale_colour_manual(values=pal,drop=FALSE)+scale_x_continuous(breaks=c(-5,5,15))+scale_y_continuous(breaks=c(-10,0,10))+
  coord_equal(xlim=xlim,ylim=ylim,expand=FALSE,clip='off')+
  labs(title=title,subtitle=subtitle,x='UMAP 1',y='UMAP 2')+base_theme
 if(main){
  p<-p+geom_segment(data=lp,aes(x=AnchorX,y=AnchorY,xend=x,yend=y),inherit.aes=FALSE,colour='#777777',linewidth=.22)+
   geom_label(data=lp,aes(x=x,y=y,label=Label),inherit.aes=FALSE,family='Arial',fontface='bold',
    size=6.8/ggplot2::.pt,colour='#333333',fill='white',linewidth=.3,
    lineheight=.95,label.padding=grid::unit(.12,'lines'),label.r=grid::unit(.11,'lines'))+
   theme(axis.line=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank())+
   annotate('segment',x=-9.3,xend=-5.6,y=-11.7,yend=-11.7,linewidth=.35,colour='#444444',arrow=grid::arrow(length=grid::unit(1.4,'mm'),type='closed'))+
   annotate('segment',x=-9.3,xend=-9.3,y=-11.7,yend=-8.0,linewidth=.35,colour='#444444',arrow=grid::arrow(length=grid::unit(1.4,'mm'),type='closed'))+
   annotate('text',x=-7.2,y=-12.5,label='UMAP 1',family='Arial',size=6/ggplot2::.pt,colour='#444444')+
   annotate('text',x=-10.1,y=-9.6,label='UMAP 2',angle=90,family='Arial',size=6/ggplot2::.pt,colour='#444444')
 }else{p<-p+theme(axis.line=element_blank(),panel.border=element_rect(fill=NA,colour='#777777',linewidth=.3))}
 p
}
subtitle<-function(s)paste0(format(nrow(s),big.mark=','),' cells shown | ',n_distinct(s$patient),' patients')
pa<-make_map(cells,'Integrated cell atlas',paste0(subtitle(cells),' | 22 libraries'),TRUE)
adj<-filter(cells,tissue=='Adjacent');tum<-filter(cells,tissue=='Tumor')
pb<-make_map(adj,'Adjacent',subtitle(adj),small=TRUE);pc<-make_map(tum,'Tumor',subtitle(tum),small=TRUE)
geometry<-lapply(list(Overall=pa,Adjacent=pb,Tumor=pc),function(p){b<-ggplot_build(p);d<-bind_rows(b$data[1:2]);pp<-b$layout$panel_params[[1]];data.frame(NPoints=nrow(d),Finite=all(is.finite(d$x)&is.finite(d$y)),XMin=pp$x.range[1],XMax=pp$x.range[2],YMin=pp$y.range[1],YMax=pp$y.range[2])})%>%bind_rows(.id='Panel')
write_csv(geometry,file.path(run,'review/rendered_geometry_checks.csv'));stopifnot(sum(geometry$NPoints[-1])==nrow(cells))
composite<-ggdraw()+draw_plot(pa,0,.035,.625,.965)+draw_plot(pb,.635,.51,.365,.49)+draw_plot(pc,.635,.01,.365,.49)
pair<-plot_grid(make_map(adj,'Adjacent',subtitle(adj)),make_map(tum,'Tumor',subtitle(tum)),nrow=1)
save_pub<-function(p,stem,w,h){
 width<-w/25.4;height<-h/25.4
 grDevices::cairo_pdf(file.path(figdir,paste0(stem,'.pdf')),width=width,height=height,family='Arial');print(p);dev.off()
 svglite::svglite(file.path(figdir,paste0(stem,'.svg')),width=width,height=height);print(p);dev.off()
 ragg::agg_png(file.path(figdir,paste0(stem,'.png')),width=width,height=height,units='in',res=300,background='white');print(p);dev.off()
 ragg::agg_tiff(file.path(figdir,paste0(stem,'.tiff')),width=width,height=height,units='in',res=600,background='white',compression='lzw');print(p);dev.off()
}
save_pub(pa,'Figure3_annotation_UMAP',120,124)
save_pub(pair,'Figure3_tissue_UMAP',183,98)
save_pub(composite,'Figure3_opening_atlas',183,130)
write_csv(data.frame(Stem=c('Figure3_annotation_UMAP','Figure3_tissue_UMAP','Figure3_opening_atlas'),WidthMM=c(120,183,183),HeightMM=c(124,98,130)),file.path(run,'source_data/export_dimensions.csv'))
writeLines(capture.output(sessionInfo()),file.path(run,'review/R_session_info.txt'))
cat('Rendered revised display subset:',nrow(cells),'cells\n')
