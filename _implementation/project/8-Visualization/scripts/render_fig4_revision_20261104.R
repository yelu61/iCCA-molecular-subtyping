# R/grid/Cairo: focused main bridge, quantitative concordance overview, continuous triple alluvial.
source('8-Visualization/scripts/fig_style.R')
source('8-Visualization/scripts/stat_annotate.R')
source('8-Visualization/scripts/render_bulk_figures.R')
source('8-Visualization/scripts/render_fig4_genomics_v2.R')
source('8-Visualization/scripts/render_fig4_alluvial.R')
suppressPackageStartupMessages(library(svglite))
OUT<-'8-Visualization/refinements/fig4_revision_20261104_v2'
PREV<-'8-Visualization/refinements/fig4_revision_20261103_v1'
BR<-'4-Multiomic/refinements/fig4_bridge_20261103_v1'
pal<-load_subtype_palette()$colors[c('C1','C2','C3')]
tx<-fig4_v2_text;rr<-fig4_v2_rect
ln<-function(x0,y0,x1,y1,col='grey65',lwd=.5,lty=1)grid::grid.lines(grid::unit(c(x0,x1),'mm'),grid::unit(c(y0,y1),'mm'),gp=grid::gpar(col=col,lwd=lwd,lty=lty))
pt<-function(x,y,col,filled=TRUE,size=1.6)grid::grid.points(grid::unit(x,'mm'),grid::unit(y,'mm'),pch=21,size=grid::unit(size,'mm'),gp=grid::gpar(col=col,fill=if(filled)col else 'white',lwd=.65))
save_draw<-function(fun,stem,w,h){
 for(ext in c('pdf','svg')){
  if(ext=='pdf')cairo_pdf(file.path(OUT,paste0(stem,'.pdf')),width=w/25.4,height=h/25.4,family=STEP8_FONT$family)else
   svglite(file.path(OUT,paste0(stem,'.svg')),width=w/25.4,height=h/25.4,system_fonts=list(sans=STEP8_FONT$family))
  fun();dev.off()
 }
 for(dpi in c(600,220))system2('pdftoppm',c('-singlefile','-r',dpi,'-png',
  shQuote(file.path(OUT,paste0(stem,'.pdf'))),shQuote(if(dpi==600)file.path(OUT,stem)else file.path(OUT,'review',paste0(basename(stem),'_preview')))),stdout=FALSE)
}
# Fixed question-driven main examples; the complete16-program six-stage matrix remains supplementary.
d<-fread(file.path(BR,'displayed_programs.csv'))
plan<-data.table(Subtype=c('C1','C2','C3','C3','C3'),Program=c('Fatty-acid metabolism','Translation','Inflammatory response','E2F targets','G2M checkpoint'),
 Label=c('C1 Fatty-acid\nmetabolism','C2 Translation','C3 Inflammatory\nresponse','C3 E2F targets','C3 G2M checkpoint'),Order=1:5)
main<-merge(d,plan,by=c('Subtype','Program'))[Stage%in%c('R1','R2','R3','P3')]
fwrite(main,file.path(OUT,'source_data/D_main_programs.csv'))
fwrite(plan,file.path(OUT,'source_data/D_fixed_program_plan.csv'))
bridge<-function(){
 tx('Key-program comparison',2,59,7.2,TRUE,just='left')
 xs<-c(32,48,64,80);stages<-c('R1','R2','R3','P3')
 labels<-c('Full RNA\n497 / all','Matched RNA\n248 / all','Matched RNA\n248 / shared','Protein\n248 / shared')
 for(j in 1:4)tx(labels[j],xs[j],52.5,5.3)
 tx('NES',21,52.5,5.3,col='grey40')
 for(i in 1:5){
  y0<-43.5-(i-1)*8.6;z<-main[Order==i][match(stages,Stage)]
  yy<-y0+z$NES*.91
  if(i%%2==1)rr(44,y0,86,7.9,'#F7F8F9')
  tx(plan$Label[i],2,y0,5.8,TRUE,col=pal[[plan$Subtype[i]]],just='left')
  ln(27,y0,84,y0,'grey75',.35,2)
  for(j in 1:3)ln(xs[j],yy[j],xs[j+1],yy[j+1],pal[[plan$Subtype[i]]],.7,if(j==3)2 else 1)
  for(j in 1:4)pt(xs[j],yy[j],pal[[plan$Subtype[i]]],z$QCommon[j]<.05,1.5)
  for(v in c(-3,3))ln(26,y0+v*.91,27,y0+v*.91,'grey55',.4)
  if(i==1){tx('3',25.3,y0+2.73,5.3,just='right');tx('0',25.3,y0,5.3,just='right');tx('−3',25.3,y0-2.73,5.3,just='right')}
 }
 pt(3,3,'grey35',TRUE,1.3);tx('q < .05',5,3,5.5,just='left')
 pt(24,3,'grey35',FALSE,1.3);tx('q ≥ .05',26,3,5.5,just='left')
 tx('RNA → protein',64,3,5.5,col='grey40')
}
save_draw(bridge,'main/Figure4D_key_program_bridge',88,61)
# Original17 rule-selected entries plus population-level summary of all2,980 sets.
all<-read_source('4-Multiomic/refinements/fig4_matched_four_gsea_20261030/paired_gsea.csv')
all[,Pattern:=ifelse(padj_RNA<.05&padj_Protein<.05,ifelse(sign(NES_RNA)==sign(NES_Protein),'Both: same sign','Both: opposite signs'),ifelse(padj_RNA<.05|padj_Protein<.05,'One layer','Neither'))]
counts<-all[,.N,by=.(Subtype,Pattern)];counts[,Fraction:=N/2980]
levels<-c('Both: same sign','Both: opposite signs','One layer','Neither')
cols<-c('#658F83','#BD8968','#AB94AF','#E1E1E1');names(cols)<-levels
counts[,Pattern:=factor(Pattern,levels=levels)];setorder(counts,Subtype,Pattern)
fwrite(counts,file.path(OUT,'source_data/E_all2980_pattern_counts.csv'))
sel<-read_source('8-Visualization/refinements/fig4_concordance_20261030_v6/source_data/selected_pathways.csv');setorder(sel,ID)
fwrite(sel,file.path(OUT,'source_data/E_representative_pathway_key.csv'))
overview<-function(){
 tx('Matched pathway concordance',2,59.2,7.0,TRUE,just='left')
 tx('All tested pathways',19,54.7,5.8,TRUE);tx('Representative comparisons',62,54.7,5.8,TRUE)
 for(k in 1:3){
  z<-counts[Subtype==names(pal)[k]];x<-c(9,19,29)[k];yy<-18
  for(i in 1:nrow(z)){
   h<-z$Fraction[i]*32;rr(x,yy+h/2,6,h,cols[as.character(z$Pattern[i])]);yy<-yy+h
  }
  tx(names(pal)[k],x,15.6,5.8,TRUE,col=pal[k])
 }
 for(v in c(0,.5,1))tx(paste0(v*100,'%'),5,18+v*32,5.2,just='right')
 L<-43;R<-81;B<-15.2;T<-53.2;xm<-62;ym<-34.2;s<-4.75
 rr(52.5,43.7,19,19,'#FCF7F2');rr(71.5,43.7,19,19,'#F0F5EF');rr(52.5,24.7,19,19,'#F0F2F4');rr(71.5,24.7,19,19,'#FCF7F2')
 ln(L,ym,R,ym);ln(xm,B,xm,T)
 for(v in seq(-4,4,2)){tx(v,xm+v*s,13.5,5.3);tx(v,41.8,ym+v*s,5.3,just='right')}
 tx('RNA NES',xm,10.2,5.8);grid::grid.text('Protein NES',grid::unit(37.8,'mm'),grid::unit(ym,'mm'),rot=90,gp=grid::gpar(fontfamily=STEP8_FONT$family,fontsize=5.8))
 # Same representative IDs; deterministic physical offsets avoid crowded points.
 p0<-fread(file.path(PREV,'review/E_point_label_positions.csv'))
 for(i in 1:nrow(sel)){
  x<-xm+sel$NES_RNA[i]*s;y<-ym+sel$NES_Protein[i]*s
  label<-p0[ID==sel$ID[i]];lx<-xm+(label$X-146)/6*s;ly<-ym+(label$Y-31)/6*s
  ln(x,y,lx,ly,'grey65',.3);tx(sprintf('%02d',sel$ID[i]),lx,ly,5.6)
  pt(x,y,pal[[sel$Subtype[i]]],sel$BothSignificant[i],1.45)
 }
 for(i in 1:4){x<-c(2,29,2,29)[i];y<-c(8,8,4.3,4.3)[i];rr(x+1,y,2,2,cols[i]);tx(levels[i],x+3,y,5.3,just='left')}
 pt(57,6,'grey35',TRUE,1.25);tx('Both q < .05',59,6,5.3,just='left')
 pt(57,2.8,'grey35',FALSE,1.25);tx('One layer q < .05',59,2.8,5.3,just='left')
}
save_draw(overview,'main/Figure4E_concordance_overview',88,61)
# One continuous three-column ribbon system; all nonzero triple paths retained.
paths<-f4all_read('path_layout');nodes<-f4all_read('node_layout');stats<-f4all_read('agreement')
triple<-function(){
 tx('MOFA correspondence',23,67,6.8,TRUE)
 tx('RNA +\nprotein',5,61.8,5.6,TRUE);tx('Primary\nRNA',23,61.8,5.6,TRUE);tx('RNA +\nprotein +\nCNV',41,61.8,5.6,TRUE)
 mapy<-function(y)12+(y-21)/64*43
 z<-copy(paths);setorder(z,-N,PathID)
 for(i in 1:nrow(z)){
  c<-adjustcolor(pal[[z$FinalSubtype[i]]],alpha.f=.44)
  f4all_ribbon(7.5,19.5,mapy(z$RPTop[i]),mapy(z$PrimaryTop[i]),mapy(z$RPBottom[i]),mapy(z$PrimaryBottom[i]),c)
  f4all_ribbon(26.5,38.5,mapy(z$PrimaryTop[i]),mapy(z$RPCTop[i]),mapy(z$PrimaryBottom[i]),mapy(z$RPCBottom[i]),c)
 }
 for(i in 1:nrow(nodes)){
  n<-nodes[i];x<-c(RP=5,Primary=23,RPC=41)[[n$Axis]];y<-mapy((n$Top+n$Bottom)/2);h<-mapy(n$Top)-mapy(n$Bottom)
  rr(x,y,if(n$Axis=='Primary')7 else 5,h,if(n$Axis=='Primary')adjustcolor(pal[[n$ClusterID]],alpha.f=.22)else '#EEEEEE',col='grey60',lwd=.35)
  tx(if(n$Axis=='Primary')n$ClusterID else paste0(n$Axis,n$ClusterID),x,y,5.0,TRUE)
 }
 tx(sprintf('ARI %.3f',stats$ARI[1]),9,7,5.6,TRUE);tx(sprintf('ARI %.3f',stats$ARI[2]),37,7,5.6,TRUE)
 tx('182 shared tumors',23,2.8,5.4,col='grey40')
}
save_draw(triple,'main/Figure4F_MOFA_continuous',46,69)
fwrite(paths,file.path(OUT,'source_data/F_original_triple_paths.csv'));fwrite(nodes,file.path(OUT,'source_data/F_original_node_layout.csv'));fwrite(stats,file.path(OUT,'source_data/F_agreement.csv'))
# Full readable17-term key retained as an independent supplement and Source Data.
for(ext in c('pdf','svg','png'))file.copy(file.path(PREV,paste0('main/Figure4E_pathway_quadrant.',ext)),file.path(OUT,paste0('supplement/Representative_pathway_key.',ext)),overwrite=TRUE)
for(p in list.files(file.path(PREV,'supplement'),full.names=TRUE))if(grepl('\\.(pdf|svg|png)$',p))file.copy(p,file.path(OUT,'supplement',basename(p)),overwrite=TRUE)
