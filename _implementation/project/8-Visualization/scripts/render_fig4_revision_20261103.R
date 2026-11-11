source('8-Visualization/scripts/fig_style.R')
source('8-Visualization/scripts/stat_annotate.R')
source('8-Visualization/scripts/render_bulk_figures.R')
source('8-Visualization/scripts/render_fig4_genomics_v2.R')
source('8-Visualization/scripts/render_fig4_alluvial.R')
suppressPackageStartupMessages(library(svglite))
OUT <- '8-Visualization/refinements/fig4_revision_20261103_v1'
BR <- '4-Multiomic/refinements/fig4_bridge_20261103_v1'
for(d in c('main','supplement','source_data','review'))dir.create(file.path(OUT,d),recursive=TRUE,showWarnings=FALSE)
pal<-load_subtype_palette()$colors[c('C1','C2','C3')]
tx<-fig4_v2_text;rr<-fig4_v2_rect
save_draw<-function(fun,stem,w,h){
 for(ext in c('pdf','svg')){
  path<-file.path(OUT,paste0(stem,'.',ext))
  if(ext=='pdf')cairo_pdf(path,width=w/25.4,height=h/25.4,family=STEP8_FONT$family) else
   svglite(path,width=w/25.4,height=h/25.4,system_fonts=list(sans=STEP8_FONT$family))
  fun();dev.off()
 }
 system2('pdftoppm',c('-singlefile','-r','600','-png',shQuote(file.path(OUT,paste0(stem,'.pdf'))),shQuote(file.path(OUT,paste0(stem,'')))),stdout=FALSE)
 system2('pdftoppm',c('-singlefile','-r','220','-png',shQuote(file.path(OUT,paste0(stem,'.pdf'))),shQuote(file.path(OUT,'review',paste0(basename(stem),'_preview')))),stdout=FALSE)
}
# Two stacked two-axis alluvials preserve actual adjacent margins of the same182 cases.
mofa<-function(){
 paths<-f4all_read('path_layout');nodes<-f4all_read('node_layout');stats<-f4all_read('agreement')
 for(j in 1:2){
  config<-c('RP','RPC')[j];top<-c(58,25.5)[j];height<-20
  tx(c('RNA + protein','RNA + protein + CNV')[j],23,top+8.5,6.7,TRUE)
  tx(sprintf('ARI = %.3f',stats$ARI[j]),23,top+4.8,6.3)
  tx('RNA',8,top+1.7,5.7,col='grey35');tx('MOFA',38,top+1.7,5.7,col='grey35')
  mapy<-function(y)top-(85-y)/64*height
  z<-copy(paths);setorder(z,-N,PathID)
  for(i in 1:nrow(z))f4all_ribbon(11,34,mapy(z$PrimaryTop[i]),mapy(z[[paste0(config,'Top')]][i]),
    mapy(z$PrimaryBottom[i]),mapy(z[[paste0(config,'Bottom')]][i]),adjustcolor(pal[[z$FinalSubtype[i]]],alpha.f=.42))
  for(axis in c('Primary',config))for(i in which(nodes$Axis==axis)){
   n<-nodes[i];yc<-mapy((n$Top+n$Bottom)/2);hc<-(n$Top-n$Bottom)/64*height
   x<-if(axis=='Primary')8 else 37
   rr(x,yc,6,hc,if(axis=='Primary')adjustcolor(pal[[n$ClusterID]],alpha.f=.24)else '#EEEEEE',col='grey60',lwd=.3)
   tx(if(axis=='Primary')n$ClusterID else paste0(config,n$ClusterID),x,yc,5.4,TRUE)
  }
 }
 fwrite(f4all_read('correspondence_tables'),file.path(OUT,'source_data/D_correspondence_counts.csv'))
 fwrite(paths,file.path(OUT,'source_data/D_original_triple_paths.csv'))
 fwrite(stats,file.path(OUT,'source_data/D_agreement.csv'))
}
save_draw(mofa,'main/Figure4D_MOFA_vertical',46,69)
# Same17 selected comparisons; only physical layout and label coordinates change.
sel<-read_source('8-Visualization/refinements/fig4_concordance_20261030_v6/source_data/selected_pathways.csv')
setorder(sel,ID);stopifnot(nrow(sel)==17)
quadrant<-function(){
 tx('RNA–protein pathway concordance',3,57.2,7.5,TRUE,just='left')
 top<-51.8;pitch<-2.68;gap<-1.1;row<-0L;layout<-list()
 for(st in names(pal)){
  z<-sel[Subtype==st];ys<-top-(row+seq_len(nrow(z))-1)*pitch-(match(st,names(pal))-1)*gap
  rr(2,(max(ys)+min(ys))/2,1,max(ys)-min(ys)+2.4,pal[[st]])
  for(i in seq_len(nrow(z))){
   label<-paste0(z$Label[i],' [',z$SourceCode[i],']')
   tx(sprintf('%02d',z$ID[i]),4,ys[i],6.6,TRUE,col=pal[[st]],just='left')
   tx(label,10,ys[i],6.6,just='left')
   width<-grid::convertWidth(grid::grobWidth(grid::textGrob(label,gp=grid::gpar(fontfamily=STEP8_FONT$family,fontsize=6.6))),'mm',TRUE)
   stopifnot(width<100)
   layout[[length(layout)+1]]<-data.table(ID=z$ID[i],Label=label,Y=ys[i],Width=width)
  };row<-row+nrow(z)
 }
 L<-122;R<-170;B<-7;T<-55;xm<-146;ym<-31
 rr(134,43,24,24,'#FCF7F2');rr(158,43,24,24,'#F0F5EF')
 rr(134,19,24,24,'#F0F2F4');rr(158,19,24,24,'#FCF7F2')
 line<-function(x0,y0,x1,y1,col='grey65',lwd=.4)grid::grid.lines(grid::unit(c(x0,x1),'mm'),grid::unit(c(y0,y1),'mm'),gp=grid::gpar(col=col,lwd=lwd))
 line(L,ym,R,ym);line(xm,B,xm,T)
 for(v in seq(-4,4,2)){
  x<-xm+v*6;y<-ym+v*6
  tx(as.character(v),x,5.1,5.8);tx(as.character(v),119.6,y,5.8,just='right')
 }
 tx('RNA NES',xm,1.6,6.4);grid::grid.text('Protein NES',grid::unit(114.8,'mm'),grid::unit(ym,'mm'),rot=90,gp=grid::gpar(fontfamily=STEP8_FONT$family,fontsize=6.4))
 tx('RNA− / protein+',L+1,T-1.2,5.4,col='grey40',just=c('left','top'))
 tx('Both higher',R-1,T-1.2,5.4,col='grey40',just=c('right','top'))
 tx('Both lower',L+1,B+1.2,5.4,col='grey40',just=c('left','bottom'))
 tx('RNA+ / protein−',R-1,B+1.2,5.4,col='grey40',just=c('right','bottom'))
 # Physical text boxes, tested against every point and previously placed label.
 px<-xm+sel$NES_RNA*6;py<-ym+sel$NES_Protein*6;boxes<-list();labels<-list()
 directions<-rbind(c(1,1),c(-1,-1),c(1,-1),c(-1,1),c(0,1),c(0,-1),c(1,0),c(-1,0))
 for(i in seq_len(nrow(sel))){
  candidates<-rbindlist(lapply(c(2.1,3,4.2,5.2),function(r)data.table(X=px[i]+directions[,1]*r,Y=py[i]+directions[,2]*r)))
  candidates[,Score:=mapply(function(x,y){
   outside<-sum(c(max(0,L+1.3-x),max(0,x-R+1.3),max(0,B+2.8-y),max(0,y-T+2.8)))*1000
   other<-sum(abs(px-x)<2.2&abs(py-y)<1.8)*100
   overlap<-sum(vapply(boxes,function(b)abs(b[1]-x)<2.6&abs(b[2]-y)<2.2,logical(1)))*1000
   outside+other+overlap+sqrt((x-px[i])^2+(y-py[i])^2)
  },X,Y)]
  best<-candidates[which.min(Score)];boxes[[i]]<-c(best$X,best$Y);labels[[i]]<-data.table(ID=sel$ID[i],X=best$X,Y=best$Y,Score=best$Score)
  line(px[i],py[i],best$X,best$Y,'grey65',.3)
  tx(sprintf('%02d',sel$ID[i]),best$X,best$Y,6.3)
 }
 for(i in seq_len(nrow(sel)))grid::grid.points(grid::unit(px[i],'mm'),grid::unit(py[i],'mm'),pch=21,size=grid::unit(1.65,'mm'),
  gp=grid::gpar(col=pal[[sel$Subtype[i]]],fill=if(sel$BothSignificant[i])pal[[sel$Subtype[i]]]else 'white',lwd=.65))
 for(k in 1:2){
  xx<-c(119,145.5)[k];grid::grid.points(grid::unit(xx,'mm'),grid::unit(58,'mm'),pch=21,size=grid::unit(1.4,'mm'),gp=grid::gpar(col='grey30',fill=if(k==1)'grey30'else 'white'))
  tx(c('Both q < .05','One layer q < .05')[k],xx+1.8,58,5.7,just='left')
 }
 fwrite(rbindlist(layout),file.path(OUT,'review/E_label_metrics.csv'));fwrite(rbindlist(labels),file.path(OUT,'review/E_point_label_positions.csv'))
}
save_draw(quadrant,'main/Figure4E_pathway_quadrant',174,60)
fwrite(sel,file.path(OUT,'source_data/E_same17_values.csv'))
# Supplement H: unmodified gene-level concordance, native horizontal alternative.
source('8-Visualization/scripts/render_fig4_concordance.R')
save_draw(function()build_fig4_concordance('D1')(),'supplement/S5H_gene_concordance',174,90)
if(file.exists(file.path(BR,'displayed_programs.csv'))){
 d<-fread(file.path(BR,'displayed_programs.csv'));d[,Stage:=factor(Stage,levels=c('R0','R1','R2','R3','P3','P4'))]
 co<-circlize::colorRamp2(c(-3.5,0,3.5),c('#527FA9','#F7F7F7','#C56760'))
 bridge<-function(){
  tx('Cohort, gene coverage and RNA–protein comparison',3,103,8,TRUE,just='left')
  labels<-c('R0','R1','R2','R3','P3','P4')
  x0<-53;step<-6.5;panel<-40
  for(k in 1:3){
   left<-x0+(k-1)*panel;rr(left+16.25,95,38,4,adjustcolor(pal[k],alpha.f=.16));tx(names(pal)[k],left+16.25,95,7,TRUE,col=pal[k])
   for(j in 1:6)tx(labels[j],left+(j-1)*step,87.2,6.1,TRUE)
  }
  plan<-unique(d[,.(pathway,Program,ProgramOrder)])[order(ProgramOrder)]
  for(i in 1:nrow(plan)){
   yy<-77-(i-1)*4.2;tx(plan$Program[i],2,yy,6.1,just='left')
   for(k in 1:3)for(j in 1:6){z<-d[pathway==plan$pathway[i]&Subtype==names(pal)[k]&as.integer(Stage)==j];xx<-x0+(k-1)*panel+(j-1)*step
    rr(xx,yy,6.1,3.8,co(max(-3.5,min(3.5,z$NES))))
    if(z$QCommon<.05)tx('\u2022',xx,yy,6.5,TRUE,col=if(abs(z$NES)>2)'white'else 'grey20')
   }
  }
  tx('R0: full RNA, pooled contrast   R1: full RNA, equal contrast   R2: matched RNA, all genes',3,10,5.2,just='left')
  tx('R3: matched RNA, shared genes   P3: protein, shared genes   P4: protein, all genes',3,7.3,5.2,just='left')
  tx('NES',48,4.0,6.2,just='right');vals<-seq(-3.5,3.5,length.out=70)
  for(i in 1:70)rr(51+(i-.5)*.6,4.0,.61,2.1,co(vals[i]))
  tx('−3.5',51,1.0,5.5);tx('0',72,1.0,5.5);tx('3.5',93,1.0,5.5)
  tx('\u2022 q < .05 across the common pathway family',102,4.0,5.7,just='left')
 }
 save_draw(bridge,'supplement/S5I_program_bridge',174,107)
 fwrite(d,file.path(OUT,'source_data/I_displayed_programs.csv'))
 # Composition is descriptive: the same patients recur across full/matched bars.
 roster<-fread(file.path(BR,'sample_rosters.csv'));count<-roster[,.N,by=.(Population,Study,FinalSubtype)]
 studycols<-load_study_palette()$colors
 composition<-ggplot(count,aes(FinalSubtype,N,fill=Study))+geom_col(width=.68)+facet_wrap(~Population)+
  scale_fill_manual(values=studycols)+labs(x=NULL,y='Samples',title='Full and matched cohort composition')+
  theme_step8()+theme(plot.title=element_text(size=8),axis.text=element_text(size=6.5),legend.text=element_text(size=5.8),legend.title=element_blank(),legend.position='bottom')
 save_draw(function()print(composition,newpage=FALSE),'supplement/S5J_cohort_composition',85,59)
 fwrite(count,file.path(OUT,'source_data/J_population_counts.csv'))
 # Counts cover all2980 paired pathways and remain descriptive of overlapping sets.
 all<-read_source('4-Multiomic/refinements/fig4_matched_four_gsea_20261030/paired_gsea.csv')
 all[,Pattern:=ifelse(padj_RNA<.05&padj_Protein<.05,ifelse(sign(NES_RNA)==sign(NES_Protein),'Both significant\nsame direction','Both significant\nopposite direction'),ifelse(padj_RNA<.05,'RNA only',ifelse(padj_Protein<.05,'Protein only','Neither')))]
 cnt<-all[,.N,by=.(Subtype,Pattern)];cnt[,Pattern:=factor(Pattern,levels=c('Both significant\nsame direction','Both significant\nopposite direction','RNA only','Protein only','Neither'))]
 p<-ggplot(cnt,aes(Subtype,N,fill=Pattern))+geom_col(width=.65)+scale_fill_manual(values=c('#6F998B','#C18B6B','#719CC1','#B79AB5','#DDDDDD'),labels=c('Both: same sign','Both: opposite signs','RNA only','Protein only','Neither'))+
  labs(x=NULL,y='Pathway comparisons',title='All tested matched pathway patterns')+theme_step8()+theme(plot.title=element_text(size=8),axis.text=element_text(size=6.5),legend.text=element_text(size=5.8),legend.title=element_blank(),legend.position='right')
 save_draw(function()print(p,newpage=FALSE),'supplement/S5K_pathway_pattern_coverage',85,59)
 fwrite(cnt,file.path(OUT,'source_data/K_all_pathway_patterns.csv'))
}
