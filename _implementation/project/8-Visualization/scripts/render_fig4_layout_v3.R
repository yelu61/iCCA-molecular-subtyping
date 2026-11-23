# Versioned Figure4 layout: current source tables, native R/grid/Cairo only.
source('8-Visualization/scripts/fig_style.R')
source('8-Visualization/scripts/stat_annotate.R')
source('8-Visualization/scripts/render_bulk_figures.R')
for(f in c('genomics','genomics_v2','proteomics','concordance','cnv','alluvial','summary_v2'))
  source(paste0('8-Visualization/scripts/render_fig4_',f,'.R'))
OUT4<-'8-Visualization/refinements/fig4_layout_20261031/'
for(d in c('main','supplement','source_data','review'))dir.create(paste0(OUT4,d),recursive=TRUE,showWarnings=FALSE)
P4<-load_subtype_palette()$colors[c('C1','C2','C3')]
text4<-fig4_v2_text;rect4<-fig4_v2_rect
# Preserve exact previous source matrices; no statistical reanalysis.
plan4<-unique(f4summary_read('display_plan')[,.(Pathway,Label)],by='Pathway')
plan4[,Row:=seq_len(.N)]
cells4<-merge(plan4,f4summary_read('all_contrasts'),by='Pathway')
cells4[,Glyph:=ifelse(FDR>=.05,'\u00b7',ifelse(NES>0,ifelse(abs(NES)>=2,'\u2191\u2191','\u2191'),ifelse(abs(NES)>=2,'\u2193\u2193','\u2193')))]
fwrite(cells4,paste0(OUT4,'source_data/G_all78_cells.csv'))
summary4<-f4summary_read('subtype_cards');fwrite(summary4,paste0(OUT4,'source_data/G_summary_cards.csv'))
summary_plot4<-function(){
 grid::grid.newpage()
 text4('Biological program',2,76,6.4,TRUE,just='left')
 for(v in 1:2){
  xx<-c(49,74)[v];rect4(xx,76,23,7,c('#EDF3FA','#EDF5EF')[v])
  text4(c('RNA','Protein')[v],xx,76,6.7,TRUE)
  for(k in 1:3)text4(names(P4)[k],xx-7+(k-1)*7,69.5,6,TRUE,col=P4[k])
 }
 text4('Genomic associations',106,76,6.4,TRUE)
 text4('Subtype synthesis',154.5,76,6.4,TRUE)
 for(i in 1:13){
  z<-plan4[Row==i];y<-65-(i-1)*4.2
  if(i%%2==1)rect4(43.5,y,85,3.9,'#F6F7F8')
  text4(z$Label,2,y,5.5,just='left')
  for(v in 1:2)for(k in 1:3){
   d<-cells4[Pathway==z$Pathway&Subtype==names(P4)[k]&View==c('RNA','Protein')[v]]
   fill<-if(d$FDR>=.05)'#EEEEEE' else if(d$NES>0)'#F7DCCB' else '#D4E4F0'
   rect4(c(42,67)[v]+(k-1)*7,y,6.2,3.7,fill)
   text4(d$Glyph,c(42,67)[v]+(k-1)*7,y,7.5,TRUE,col=if(d$FDR>=.05)'grey55' else if(d$NES>0)'#A95334' else '#2C719E')
  }
 }
 syn<-list(
  c('Metabolic specialization','RNA + protein: metabolic programs',
    'Inferred TF: HNF4A / HNF1A / CEBPA / PPARA','Immune / myeloid signals (method-dependent)'),
  c('Morphogenesis / translation','WNT / cilium concordant; translation discordant',
    'Inferred TF: HNF1B / ONECUT1 (not exclusive)','Lower immune / stromal scores'),
  c('Hypoxia / epithelial remodeling','Hypoxia / epithelial support; cell cycle: full-RNA',
    'Inferred TF: KLF5 / GRHL2 / HIF1A','CAF: MCP-counter / EPIC; xCell unclear'))
 fwrite(rbindlist(lapply(1:3,function(k)data.table(Subtype=names(P4)[k],Line=1:4,Text=syn[[k]]))),paste0(OUT4,'source_data/G_displayed_synthesis.csv'))
 for(k in 1:3){
  top<-67-(k-1)*19.2
  rect4(106,top-8,35,17.8,'#FAFAFA',col='grey85')
  rect4(154.5,top-8,53,17.8,grDevices::adjustcolor(P4[k],alpha.f=.05),col=grDevices::adjustcolor(P4[k],alpha.f=.4))
  rect4(106,top,35,3,grDevices::adjustcolor(P4[k],alpha.f=.2))
  text4(names(P4)[k],90,top,6,TRUE,col=P4[k],just='left')
  dna<-switch(k,'No FDR-significant\ndriver enrichment',
    'BAP1, IDH1/2 enriched\n1p36 loss > C3', 'KRAS, SMAD4, TP53\nenriched; 1p36 loss < C2')
  text4(dna,90,top-8,5.3,just='left')
  text4(paste(names(P4)[k],syn[[k]][1]),129.5,top,5.9,TRUE,col=P4[k],just='left')
  for(j in 2:4)text4(syn[[k]][j],129.5,top-4-(j-2)*3.6,5.05,just='left')
 }
}
# Rescale native grob physical geometry while preserving a final-size text floor.
# This is not a PDF raster/composite: text and shapes remain native Cairo objects.
u4<-function(x,s){
 if(!grid::is.unit(x))return(x)
 for(i in seq_along(x))if(grid::unitType(x[i])%in%c('mm','cm','inches','points','bigpts'))x[i]<-x[i]*s
 x
}
resize4<-function(g,sx,sy,fs=.8){
 if(!is.list(g)||grid::is.unit(g))return(g)
 if(identical(class(g),'viewport'))g$mask<-TRUE
 for(n in names(g)){
  if(is.null(g[[n]]))next
  if(n%in%c('x','x0','x1','width','widths'))g[[n]]<-u4(g[[n]],sx)
  else if(n%in%c('y','y0','y1','height','heights'))g[[n]]<-u4(g[[n]],sy)
  else if(n%in%c('r','size'))g[[n]]<-u4(g[[n]],min(sx,sy))
  else if(n=='gp'&&!is.null(g[[n]]$fontsize))g[[n]]$fontsize<-if(inherits(g,'text'))pmax(5,g[[n]]$fontsize*fs) else g[[n]]$fontsize*min(sx,sy)
  else if(is.list(g[[n]])&&!grid::is.unit(g[[n]]))g[[n]]<-resize4(g[[n]],sx,sy,fs)
 }
 # gList has names but vpList may be unnamed.
 if(is.null(names(g)))for(i in seq_along(g))if(is.list(g[[i]]))g[[i]]<-resize4(g[[i]],sx,sy,fs)
 g
}
capture4<-function(p,w,h){
 if(is.function(p)){b<-paste(deparse(body(p)),collapse='\n');b<-gsub('grid::grid.draw(grobs[[k]])',"grobs[[k]]$name<-paste0('plot_',k);grid::grid.draw(grobs[[k]])",b,fixed=TRUE);body(p)<-parse(text=b)[[1]]}
 grid::grid.grabExpr({if(is.function(p))p() else print(p)},width=w/25.4,height=h/25.4,wrap=FALSE,wrap.grobs=FALSE)}
# Move explanatory text to the caption while retaining axes and decoding keys.
trim_text4<-function(g,panel){
 if(inherits(g,'text')){
  lab<-paste(as.character(g$label),collapse=' ')
  remove<-switch(panel,
   'main/Figure4A'=lab=='Driver mutation landscape',
   'main/Figure4C'=grepl('^NES bars:',lab),
   'main/Figure4D'=grepl('248 matched tumors',lab,fixed=TRUE),
   'main/Figure4E'=grepl('^Matched cohort:',lab),
   'main/Figure4F'=grepl('^(MOFA clusters across|Same 182 Main tumors|Both fits:)',lab),
   FALSE)
  if(remove)g$label<-''
  if(panel=='main/Figure4B'&&lab=='Study-stratified odds ratio (95% CI)')g$label<-'Odds ratio (95% CI)'
  if(panel=='main/Figure4C'&&lab=='Study-adjusted standardized difference')g$label<-'Standardized difference'
 }
 if(!is.null(g$children))for(i in seq_along(g$children))g$children[[i]]<-trim_text4(g$children[[i]],panel)
 if(!is.null(g$grobs))for(i in seq_along(g$grobs))g$grobs[[i]]<-trim_text4(g$grobs[[i]],panel)
 g
}
export4<-function(p,name,ow,oh,w,h,fs=.8){
 cat(name,'capture\n');g<-trim_text4(capture4(p,ow,oh),name);cat(name,'resize\n');padleft<-if(name%in%c('supplement/Supplement_A','supplement/Supplement_C'))2.5 else 0
 g<-resize4(g,(w-padleft)/ow,h/oh,fs)
 if(padleft>0)g$vp<-grid::viewport(x=grid::unit(padleft,'mm'),width=grid::unit(w-padleft,'mm'),just='left')
 cat(name,'draw\n');path<-paste0(OUT4,name,'.pdf');grDevices::cairo_pdf(path,width=w/25.4,height=h/25.4,family=STEP8_FONT$family)
 grid::grid.newpage();grid::grid.draw(g);grDevices::dev.off()
 saveRDS(g,paste0(OUT4,name,'.rds'));g
}
# Supplement pathway heatmap: every one of the73 existing terms retained.
supp_path4<-function(){
 grid::grid.newpage();d<-f4p_read('supplement_gsea_selection')
 p<-unique(d[,.(Pathway,Label,RowOrder,Collection)])[order(RowOrder)]
 collections<-c('Hallmark','GO_BP','Reactome','KEGG')
 p[,Column:=match(Collection,collections)];p[,Within:=seq_len(.N),by=Column]
 p[,PlotLabel:=vapply(Label,f4p_wrap_term,character(1),width=39.8,size=5.1)]
 p[,Lines:=lengths(strsplit(PlotLabel,'\n',fixed=TRUE))]
 stopifnot(nrow(p)==73L,all(p$Lines<=2L))
 fwrite(p,paste0(OUT4,'source_data/S_E_73_pathway_layout.csv'))
 fwrite(d,paste0(OUT4,'source_data/S_E_219_displayed_values.csv'))
 cf<-circlize::colorRamp2(c(-3.5,0,3.5),c(SCALE_EFFECT_LOW,SCALE_EFFECT_MID,SCALE_EFFECT_HIGH))
 sc<-c(Hallmark='#807694',GO_BP='#6E91AB',Reactome='#6F9B86',KEGG='#B19A6A')
 text4('Proteomic pathway programs',2,106,7.5,TRUE,just='left')
 for(b in 1:4){
  left<-(b-1)*55.5
  rect4(left+27,101,53.5,4.2,grDevices::adjustcolor(sc[b],alpha.f=.13))
  text4(c('Hallmark','GO BP','Reactome','KEGG MEDICUS')[b],left+27,101,6.3,TRUE)
  for(k in 1:3){rect4(left+43+(k-1)*4,96,3.6,3.5,P4[k]);text4(names(P4)[k],left+43+(k-1)*4,96,5.2,TRUE)}
  for(i in which(p$Column==b)){
   y<-90-(p$Within[i]-1)*3.75
   text4(p$PlotLabel[i],left+39.8,y,5.1,just='right')
   rect4(left+40.7,y,.65,3.5,sc[[p$Collection[i]]])
   for(k in 1:3){z<-d[Pathway==p$Pathway[i]&Subtype==names(P4)[k]]
    rect4(left+43+(k-1)*4,y,3.6,3.5,cf(pmax(-3.5,pmin(3.5,z$NES))))
    if(z$padj<.05)text4(fdr_label(z$padj),left+43+(k-1)*4,y,5.2,TRUE,col=if(abs(z$NES)>2)'white' else 'black')
   }
  }
 }
 text4('Protein NES',68,7,5.5,just='right')
 vals<-seq(-3.5,3.5,length.out=60);for(i in 1:60)rect4(71+(i-.5)*.65,7,.66,2.3,cf(vals[i]))
 for(i in 1:3)text4(c('-3.5','0','3.5')[i],c(71,90.5,110)[i],3.5,5.2)
 text4('* q < .05; ** q < .01; *** q < .001',118,7,5.2,just='left')
}

# Local layout variants only; canonical renderers and analytical tables are unchanged.
variant4<-function(p,replacements){
 b<-paste(deparse(body(p)),collapse='\n')
 for(k in seq_along(replacements))b<-gsub(names(replacements)[k],replacements[[k]],b,fixed=TRUE)
 body(p)<-parse(text=b)[[1]];p
}
b4<-variant4(build_fig4_v2('B'),c('top <- 96'='top <- 83','bottom <- 27'='bottom <- 30',
 ', 114, 8'=', 98, 8',', 104, 8'=', 88, 8'))
c4<-build_fig4_proteomics('C')
# Compact conventional aliases preserve single-line terms and exact pathway IDs.
alias4<-c('Interferon-gamma response'='IFN-gamma response','Cell junction organization'='Cell junctions',
 'Cholesterol homeostasis'='Cholesterol balance','RAB7 minus-end transport'='RAB7 minus-end traffic',
 'Organic acid metabolism'='Organic acid metabolism','Oxidative phosphorylation'='OXPHOS')
c4<-variant4(c4,c('lab <- z$Label[j]'="lab <- z$Label[j];if(lab %in% names(alias4))lab <- alias4[[lab]]",
 'file.path(STEP8_ROOT, F4P, "term_text_fit.csv")'="paste0(OUT4,'review/C_term_text_fit.csv')"))
# The shorter quadrant sidebar uses the same exact pathway identities and all17 points.
wrapE4<-function(text) f4p_wrap_term(text,width=44,size=7.8)
e4<-variant4(build_fig4_concordance('D2'),c('f4d_wrap_rule_label('='wrapE4(',
 'rowgap > 1.5'='rowgap > .2',
 'file.path(STEP8_ROOT, F4D, "D2_text_metrics.csv")'="paste0(OUT4,'review/E_text_metrics.csv')",
 'file.path(STEP8_ROOT, F4D, "D2_layout_metrics.csv")'="paste0(OUT4,'review/E_layout_metrics.csv')"))
# Avoid freezing ggrepel native coordinates on the old, larger drawing device.
bE4<-paste(deparse(body(e4)),collapse='\n')
bE4<-sub('grid::grid.force\\(\\)[\\s\\S]*?grid::upViewport\\(0\\)',
 'left<-74;right<-166;bottompanel<-23;toppanel<-105',bE4,perl=TRUE)
body(e4)<-parse(text=gsub('F4D','OUT4',bE4,fixed=TRUE))[[1]]
c4<-variant4(c4,c('F4P'='OUT4'))
sa4<-variant4(build_fig4_v2('S10I'),c('Exclusivity'='Exclusive',
 'grid::unit(cbx, "mm")'='grid::unit(cbx-5, "mm")',
 'grid::unit(cbx + cbw, "mm")'='grid::unit(cbx+cbw+5, "mm")',
 'if (st == "C1") "white" else "black"'='"black"'))
sb4<-variant4(build_fig4_cnv('S10Q'),c(
 'Omnibus differences: S10Q; pairwise contrasts: S10H'='',
 'Colored ticks: q < 0.05'='q < 0.05',
 'Mean gene-level call frequency (%)'='Gene-level call frequency (%)'))
h4<-build_subfigure('Figure 4','E')
h4$data$Factor<-factor(h4$data$Factor,levels=paste0('Factor',10:1))
fillscale4<-h4$scales$get_scales('fill');fillscale4$name<-'Variance (%)'
h4<-h4+guides(fill=guide_colorbar(title.position='top',title.hjust=.5))+labs(title='MOFA variance explained',fill='Variance (%)',caption='Three-view fit: 204 tumors; one seed.')+
 theme(legend.position='bottom',legend.direction='horizontal',legend.key.width=grid::unit(8,'mm'))

fwrite(data.table(Original=names(alias4),Display=unname(alias4)),paste0(OUT4,'source_data/C_label_aliases.csv'))
# Main: explicit source-to-new-letter mapping, separate from legacy panel IDs.
only_suppE4<-isTRUE(getOption('fig4.layout.only_supplement_e',FALSE))
only_suppB4<-isTRUE(getOption('fig4.layout.only_supplement_b',FALSE))
targeted4<-only_suppE4||only_suppB4
stopifnot(!(only_suppE4&&only_suppB4))
if(targeted4){
 supp4<-setNames(lapply(c('A','B','C','D','E','F','G_CNV_burden','H_MOFA_variance'),function(x)
  readRDS(paste0(OUT4,'supplement/Supplement_',x,'.rds'))),LETTERS[1:8])
 if(only_suppE4)supp4$E<-export4(supp_path4,'supplement/Supplement_E',220,110,220,110,1)
 if(only_suppB4)supp4$B<-export4(sb4,'supplement/Supplement_B',174,132,110,60,.78)
}else{
main4<-list()
main4$A<-export4(build_fig4_v2('A'),'main/Figure4A',174,102,96,59,.75)
main4$B<-export4(b4,'main/Figure4B',174,102,81,59,.72)
main4$C<-export4(c4,'main/Figure4C',212,120,130,84,.78)
main4$D<-export4(build_fig4_concordance('D1_vertical'),'main/Figure4D',85,216,47,84,.72)
main4$E<-export4(e4,'main/Figure4E',168,120,96,65,.64)
main4$F<-export4(build_fig4_alluvial('G'),'main/Figure4F',174,116,81,65,.72)
main4$G<-export4(summary_plot4,'main/Figure4G',183,80,183,80,1)
supp4<-list()
supp4$A<-export4(sa4,'supplement/Supplement_A',174,91,88,65,.78)
supp4$B<-export4(sb4,'supplement/Supplement_B',174,132,110,60,.78)
supp4$C<-export4(build_fig4_cnv('S10H'),'supplement/Supplement_C',174,144,88,87,.78)
supp4$D<-export4(build_fig4_proteomics('S10N'),'supplement/Supplement_D',174,160,88,126,.78)
supp4$E<-export4(supp_path4,'supplement/Supplement_E',220,110,220,110,1)
supp4$F<-export4(build_fig4_alluvial('S10S'),'supplement/Supplement_F',174,90,88,57,.8)
supp4$G<-export4(build_fig4_cnv('S10A'),'supplement/Supplement_G_CNV_burden',100,85,88,62,.8)
supp4$H<-export4(h4,'supplement/Supplement_H_MOFA_variance',110,85,88,62,.8)
}
place4<-function(g,x,y,w,h){grid::pushViewport(grid::viewport(x=grid::unit(x,'mm'),y=grid::unit(y,'mm'),width=grid::unit(w,'mm'),height=grid::unit(h,'mm'),just=c('left','bottom')));grid::grid.draw(g);grid::popViewport()}
page4<-function(path,w,h,items,grobs,title=NULL){
 grDevices::cairo_pdf(path,width=w/25.4,height=h/25.4,family=STEP8_FONT$family);grid::grid.newpage()
 for(i in seq_len(nrow(items))){z<-items[i];place4(grobs[[z$Panel]],z$X,z$Y,z$W,z$H);text4(z$Panel,z$X,z$Y+z$H+2,9,TRUE,just='left')}
 grDevices::dev.off()
}
lay<-data.table(Panel=LETTERS[1:7],X=c(3,105,3,139,3,105,3),Y=c(258,258,167,167,95,95,7),W=c(96,81,130,47,96,81,183),H=c(59,59,84,84,65,65,80))
fwrite(lay,paste0(OUT4,'source_data/main_layout.csv'));if(!targeted4)page4(paste0(OUT4,'Figure4_updated.pdf'),189,324,lay,main4)
# Two supplementary pages preserve readable73-term coverage and add the existing null/basis results.
s1<-data.table(Panel=c('A','B','C','G'),X=c(3,97,3,97),Y=c(114,119,10,35),W=c(88,110,88,110),H=c(65,60,87,62))
s2<-data.table(Panel=c('D','F','H','E'),X=c(3,118,118,3),Y=c(131,200,131,15),W=c(105,105,105,220),H=c(126,57,62,110))
genomic_page4<-supp4;genomic_page4$G<-resize4(supp4$G,110/88,1,1)
if(!only_suppE4)page4(paste0(OUT4,'Supplement_genomic.pdf'),210,186,s1,genomic_page4)
supp_page4<-supp4
for(k in c('D','F','H'))supp_page4[[k]]<-resize4(supp4[[k]],105/88,1,1)
if(!only_suppB4)page4(paste0(OUT4,'Supplement_proteomic_integration.pdf'),226,264,s2,supp_page4)
fwrite(rbind(cbind(Page='genomic',s1),cbind(Page='proteomic_integration',s2)),paste0(OUT4,'source_data/supplement_layout.csv'))
for(f in if(only_suppB4)'Supplement_genomic' else if(only_suppE4)'Supplement_proteomic_integration' else c('Figure4_updated','Supplement_genomic','Supplement_proteomic_integration'))
 system2('pdftoppm',c('-scale-to','2200','-png','-singlefile',shQuote(paste0(OUT4,f,'.pdf')),shQuote(paste0(OUT4,'review/',f))))

for(f in if(only_suppB4)'Supplement_genomic' else if(only_suppE4)'Supplement_proteomic_integration' else c('Figure4_updated','Supplement_genomic','Supplement_proteomic_integration'))
 system2('pdftoppm',c('-r','300','-png','-singlefile',shQuote(paste0(OUT4,f,'.pdf')),shQuote(paste0(OUT4,f))))
if(!targeted4)system2('pdftoppm',c('-r','180','-png','-singlefile',shQuote(paste0(OUT4,'main/Figure4G.pdf')),shQuote(paste0(OUT4,'review/Figure4G'))))

if(!only_suppB4)system2('pdftoppm',c('-r','220','-png','-singlefile',shQuote(paste0(OUT4,'supplement/Supplement_E.pdf')),shQuote(paste0(OUT4,'review/Supplement_E'))))

if(!only_suppE4)system2('pdftoppm',c('-r','220','-png','-singlefile',shQuote(paste0(OUT4,'supplement/Supplement_B.pdf')),shQuote(paste0(OUT4,'review/Supplement_B'))))
