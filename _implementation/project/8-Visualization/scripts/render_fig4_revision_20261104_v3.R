# R/grid/Cairo: focused main bridge, quantitative concordance overview, continuous triple alluvial.
source('8-Visualization/scripts/fig_style.R')
source('8-Visualization/scripts/stat_annotate.R')
source('8-Visualization/scripts/render_bulk_figures.R')
source('8-Visualization/scripts/render_fig4_genomics_v2.R')
source('8-Visualization/scripts/render_fig4_alluvial.R')
suppressPackageStartupMessages(library(svglite))
OUT<-'8-Visualization/refinements/fig4_revision_20261104_v3'
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
 xs<-c(31,44,57,70);stages<-c('R1','R2','R3','P3')
 labels<-c('Full RNA\n497 / all','Matched RNA\n248 / all','Matched RNA\n248 / shared','Protein\n248 / shared')
 for(j in 1:4)tx(labels[j],xs[j],52.5,5.3)
 tx('NES',21,52.5,5.3,col='grey40')
 for(i in 1:5){
  y0<-43.5-(i-1)*8.6;z<-main[Order==i][match(stages,Stage)]
  yy<-y0+z$NES*.91
  if(i%%2==1)rr(39,y0,76,7.9,'#F7F8F9')
  tx(plan$Label[i],2,y0,5.8,TRUE,col=pal[[plan$Subtype[i]]],just='left')
  ln(27,y0,74,y0,'grey75',.35,2)
  for(j in 1:3)ln(xs[j],yy[j],xs[j+1],yy[j+1],pal[[plan$Subtype[i]]],.7,if(j==3)2 else 1)
  for(j in 1:4)pt(xs[j],yy[j],pal[[plan$Subtype[i]]],z$QCommon[j]<.05,1.5)
  for(v in c(-3,3))ln(26,y0+v*.91,27,y0+v*.91,'grey55',.4)
  if(i==1){tx('3',25.3,y0+2.73,5.3,just='right');tx('0',25.3,y0,5.3,just='right');tx('−3',25.3,y0-2.73,5.3,just='right')}
 }
 pt(3,3,'grey35',TRUE,1.3);tx('q < .05',5,3,5.5,just='left')
 pt(24,3,'grey35',FALSE,1.3);tx('q ≥ .05',26,3,5.5,just='left')
 tx('RNA → protein',59,3,5.5,col='grey40')
}
save_draw(bridge,'main/Figure4D_key_program_bridge',78,61)

# Original17 comparisons: detached bars do not change the statistical selection.
sel<-read_source('8-Visualization/refinements/fig4_concordance_20261030_v6/source_data/selected_pathways.csv');setorder(sel,ID)
fwrite(sel,file.path(OUT,'source_data/E_representative_pathway_key.csv'))
quad<-function(){
 tx('Matched pathway concordance',2,59.2,7,TRUE,just='left')
 L<-22;R<-66;B<-11;T<-55;xm<-44;ym<-33;s<-5.5
 rr(33,44,22,22,'#FCF7F2');rr(55,44,22,22,'#F0F5EF');rr(33,22,22,22,'#F0F2F4');rr(55,22,22,22,'#FCF7F2')
 ln(L,ym,R,ym);ln(xm,B,xm,T)
 for(v in seq(-4,4,2)){tx(v,xm+v*s,9.2,5.8);tx(v,20.7,ym+v*s,5.8,just='right')}
 tx('RNA NES',xm,5.7,6.2);grid::grid.text('Protein NES',grid::unit(16.4,'mm'),grid::unit(ym,'mm'),rot=90,gp=grid::gpar(fontfamily=STEP8_FONT$family,fontsize=6.2))
 tx('RNA− / protein+',L+1,T-1.2,5.0,col='grey40',just=c('left','top'))
 tx('Both higher',xm+1,T-1.2,5.0,col='grey40',just=c('left','top'))
 tx('Both lower',L+1,B+1.2,5.0,col='grey40',just=c('left','bottom'))
 tx('RNA+ / protein−',R-1,B+1.2,5.0,col='grey40',just=c('right','bottom'))
 p0<-fread(file.path(PREV,'review/E_point_label_positions.csv'))
 for(i in 1:nrow(sel)){
  x<-xm+sel$NES_RNA[i]*s;y<-ym+sel$NES_Protein[i]*s
  label<-p0[ID==sel$ID[i]];lx<-xm+(label$X-146)/6*s;ly<-ym+(label$Y-31)/6*s
  ln(x,y,lx,ly,'grey65',.3);tx(sprintf('%02d',sel$ID[i]),lx,ly,6.0)
  pt(x,y,pal[[sel$Subtype[i]]],sel$BothSignificant[i],1.6)
 }
 pt(17,1.8,'grey35',TRUE,1.3);tx('Both q < .05',19,1.8,5.4,just='left')
 pt(43,1.8,'grey35',FALSE,1.3);tx('One layer q < .05',45,1.8,5.4,just='left')
}
save_draw(quad,'main/Figure4E_pathway_quadrant',76,61)
counts<-read_source('8-Visualization/refinements/fig4_revision_20261104_v2/source_data/E_all2980_pattern_counts.csv')
lev<-c('Both: same sign','Both: opposite signs','One layer','Neither')
cols<-setNames(c('#658F83','#BD8968','#AB94AF','#E1E1E1'),lev)
counts[,Pattern:=factor(Pattern,levels=lev)];setorder(counts,Subtype,Pattern)
fwrite(counts,file.path(OUT,'source_data/E_all2980_pattern_counts.csv'))
bars<-function(){
 tx('Pathway concordance patterns',2,59,7,TRUE,just='left')
 for(k in 1:3){
  z<-counts[Subtype==names(pal)[k]];x<-c(21,35,49)[k];yy<-20
  for(i in 1:nrow(z)){h<-z$Fraction[i]*32;rr(x,yy+h/2,9,h,cols[as.character(z$Pattern[i])]);yy<-yy+h}
  tx(names(pal)[k],x,16.8,6.4,TRUE,col=pal[k])
 }
 for(v in c(0,.5,1))tx(paste0(v*100,'%'),13,20+v*32,5.8,just='right')
 for(i in 1:4){x<-c(3,30,3,30)[i];y<-c(10,10,6,6)[i];rr(x+1,y,2,2,cols[i]);tx(lev[i],x+3,y,5.6,just='left')}
}
save_draw(bars,'independent/Figure4E_pathway_pattern_bars',62,61)
