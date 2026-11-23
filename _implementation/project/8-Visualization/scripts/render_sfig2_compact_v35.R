# Four standalone supplementary panels; frozen source statistics only.
# Quantitative grid: complete network context, study robustness and six-state distributions.
# R-only drawing; PDF editable text plus 600-dpi PNG; old figures remain unchanged.
suppressPackageStartupMessages({library(data.table);library(grid)})
source('8-Visualization/scripts/fig_style.R')
out <- '8-Visualization/refinements/sfig2_compact_20261031_v35'
for(d in c('figures','tables','review'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
pal <- load_subtype_palette()$colors; study_pal <- load_study_palette()
states <- c('C1','C2','C3','C1_C2','C1_C3','C2_C3'); ns <- c(112,223,162,41,39,9)
src <- c(classification='2-TranscriptomicCharacterization/publication/tables/step2_wgcna_module_subtype_classification.csv',
 meta='2-TranscriptomicCharacterization/publication/tables/step2_wgcna_subtype_random_effects_meta.csv',
 main='8-Visualization/refinements/fig2_selection_20261031_v32/tables/C_all11_modules.csv',
 state='2-TranscriptomicCharacterization/runs/fig2_intermediate_20261026_v1/wgcna_state_summary.csv',
 state_scores='2-TranscriptomicCharacterization/runs/fig2_intermediate_20261026_v1/wgcna_all_states_long.csv',
 patterns='8-Visualization/refinements/sfig2_states_20261030_v26/tables/pattern_scores_all_states.csv',
 tme='8-Visualization/refinements/fig2_transcriptomic_20261027_v13_tme/tables/robustness_selected_effects.csv',
 study='8-Visualization/refinements/fig2_transcriptomic_20261027_v13_tme/tables/robustness_selected_perstudy.csv',
 loso='2-TranscriptomicCharacterization/runs/20261027_tme_display_v3/tables/loso_contrasts.csv',
 palette='governance/data/subtype_label_color_dictionary.csv',study_palette='0-Data/qc/study_color_dictionary.csv')
hashes <- tools::md5sum(src)
cl <- fread(src['classification']); main <- fread(src['main']); meta <- fread(src['meta']); st <- fread(src['state']); ss <- fread(src['state_scores']); scores <- fread(src['patterns'])
stopifnot(nrow(cl)==24,nrow(main)==11,nrow(meta)==72,nrow(st)==144,nrow(ss)==14064,nrow(scores)==3516)
groups <- c('C1_specific','C2_specific','C3_specific','SharedHighCompatible_C1_C2','SharedHighCompatible_C2_C3','Replicated_gradient_or_ambiguous','Non_specific')
group_labels <- c('C1-specific','C2-specific','C3-specific','Shared-high C1 / C2','Shared-high C2 / C3','Gradient / ambiguous','Non-specific')
cl[,MainFigure:=Module %in% main$Module];cl[,GroupOrder:=match(SpecificityClass,groups)]
cl[,WithinOrder:=match(Module,main$Module)];cl[is.na(WithinOrder),WithinOrder:=100+match(ModuleColor,sort(ModuleColor))]
setorder(cl,GroupOrder,WithinOrder);cl[,RowOrder:=.I]
stopifnot(sum(cl$MainFigure)==11,all(!is.na(cl$GroupOrder)))
fwrite(cl,file.path(out,'tables/module_order.csv'))
me <- merge(cl,meta,by='Module'); setorder(me,RowOrder,TargetSubtype);fwrite(me,file.path(out,'tables/network_principal_effects.csv'))
state <- merge(cl,st,by='Module');setorder(state,RowOrder,FinalSubtype);fwrite(state,file.path(out,'tables/network_state_medians.csv'))
chk <- ss[,.(N=.N,Median=median(Score),Q25=quantile(Score,.25),Q75=quantile(Score,.75)),by=.(Module,FinalSubtype)]
ck <- merge(chk,st,by=c('Module','FinalSubtype'));stopifnot(all(ck$N.x==ck$N.y),max(abs(ck$Median.x-ck$Median.y))<1e-12)
fwrite(scores,file.path(out,'tables/pattern_scores.csv'))
ps <- scores[,.(N=.N,Q25=quantile(Score,.25),Median=median(Score),Q75=quantile(Score,.75)),by=.(Pattern,FinalSubtype)]
fwrite(ps,file.path(out,'tables/pattern_summary.csv'))
tme <- fread(src['tme']); per <- fread(src['study']); lo <- fread(src['loso'])
stopifnot(nrow(tme)==99,nrow(per)==594)
lo <- merge(lo,tme[,.(Method,Score,Subtype,PrimaryEffect=Effect,PrimaryN=N)],by=c('Method','Score','Subtype'))
lo[,InformativeOmission:=N<PrimaryN]
ls <- lo[InformativeOmission==TRUE,.(Retained=sum(Estimable & is.finite(Effect) & sign(Effect)==sign(PrimaryEffect)),EstimableOmissions=sum(Estimable & is.finite(Effect)),TotalOmissions=.N),by=.(Method,Score,Subtype)]
tme <- merge(tme,ls,by=c('Method','Score','Subtype'));stopifnot(nrow(tme)==99,all(tme$EstimableOmissions>0))
fwrite(tme,file.path(out,'tables/TME_primary_LOSO_sensitivity.csv'));fwrite(per,file.path(out,'tables/TME_per_study.csv'));fwrite(lo,file.path(out,'tables/TME_LOSO_all.csv'))
mm <- function(x)unit(x,'mm'); texts<-list();marks<-list(); manifest<-list()
txt <- function(s,x,y,size=6,face='plain',just='left',col='#24282B',rot=0){grid.text(s,mm(x),mm(y),just=just,rot=rot,gp=gpar(fontfamily='Arial',fontsize=size,fontface=face,col=col));texts[[length(texts)+1L]]<<-data.table(Text=as.character(s),X=x,Y=y,Size=size)}
rect <- function(x,y,w,h,fill,border=NA)grid.rect(mm(x),mm(y),mm(w),mm(h),gp=gpar(fill=fill,col=border,lwd=.4))
line <- function(x,y,col='#D3D8DC',lwd=.5,lty=1)grid.lines(mm(x),mm(y),gp=gpar(col=col,lwd=lwd,lty=lty))
point <- function(x,y,pch=16,size=1,col='#303840',fill=col)grid.points(mm(x),mm(y),pch=pch,size=mm(size),gp=gpar(col=col,fill=fill,lwd=.55))
ramp<-colorRamp(c(SCALE_EFFECT_LOW,'#FAFAFA',SCALE_EFFECT_HIGH))
colval<-function(v,lim=2)rgb(ramp(pmax(0,pmin(1,(v/lim+1)/2))),maxColorValue=255)
bar <- function(x,y,label,w=27,lim=2){txt(label,x,y+4,5.8);for(k in 1:100)rect(x+(k-.5)*w/100,y,w/100+.01,1.8,colval(-lim+(k-1)/99*2*lim,lim));for(v in c(-lim,0,lim))txt(v,x+(v+lim)/(2*lim)*w,y-3,5.3,just='center')}
emit <- function(name,draw,w,h,figure,panel){
 pdf<-file.path(out,'figures',paste0(name,'.pdf'))
 grDevices::cairo_pdf(pdf,width=w/25.4,height=h/25.4,family='Arial');draw();dev.off()
 fwrite(rbindlist(texts),file.path(out,'tables',paste0(name,'_text.csv')));if(length(marks))fwrite(rbindlist(marks,fill=TRUE),file.path(out,'tables',paste0(name,'_marks.csv')))
 ragg::agg_png(sub('pdf$','png',pdf),width=w,height=h,units='mm',res=600,background='white');draw();dev.off()
 ragg::agg_png(file.path(out,'review',paste0(name,'.png')),width=w,height=h,units='mm',res=180,background='white');draw();dev.off()
 manifest[[length(manifest)+1L]]<<-data.table(Figure=figure,Panel=panel,TargetOutput=pdf,WidthMM=w,HeightMM=h,Status='Rendered review candidate',SourceTable=paste(src,collapse=';'),SourceMD5=paste(hashes,collapse=';'))
}
module_draw <- function(is_state=FALSE){function(){
 grid.newpage();texts<<-list();marks<<-list();H<-153;W<-if(is_state)115 else 90
 xs<-if(is_state)c(43,54.5,66,81,92.5,104) else c(49,64,79);cw<-if(is_state)10.8 else 14.2
 for(j in seq_along(xs)){txt(gsub('_','/',states[j]),xs[j],H-5,6.4,'bold','center');rect(xs[j],H-8,cw,1,pal[states[j]]);if(is_state)txt(paste0('n=',ns[j]),xs[j],H-11.5,5.2,just='center')}
 txt('Module',4,H-5,6.3,'bold');txt('Genes',33,H-5,5.8,just='right');y<-H-18
 for(g in seq_along(groups)){
  txt(group_labels[g],4,y,5.8,'bold',col='#56636D');line(c(4,W-4),rep(y-2,2),'#E4E8EA');y<-y-4
  for(i in cl[GroupOrder==g,RowOrder]){r<-cl[RowOrder==i];rect(4.8,y,1.5,3.2,r$ModuleColor)
   txt(r$ModuleColor,7,y,5.9,face=if(r$MainFigure)'bold'else'plain');txt(r$ModuleSize,33,y,5.3,just='right',col='#66717A')
   if(r$MainFigure)point(36,y,size=.65)
   for(j in seq_along(xs)){
    q<-if(is_state)state[Module==r$Module & FinalSubtype==states[j]] else me[Module==r$Module & TargetSubtype==states[j]]
    v<-if(is_state)q$Median else q$RandomEffect;stopifnot(length(v)==1,is.finite(v));rect(xs[j],y,cw,3.1,colval(v))
    if(is_state)txt(sprintf('%.2f',if(abs(v)<.005)0 else v),xs[j],y,5.2,just='center',col=if(abs(v)>1.3)'white'else'#24282B') else if(q$FDR<.05)point(xs[j],y,size=.65,col=if(abs(v)>1.3)'white'else'#24282B')
    marks[[length(marks)+1L]]<<-data.table(Module=r$Module,State=states[j],Value=v,X=xs[j],Y=y)
   };y<-y-3.4
  };y<-y-1
 }
 stopifnot(y>15)
 bar(5,8,if(is_state)'Median frozen score (SD)'else'Meta-analysis effect (SD)',w=27)
 point(43,12,size=.65);txt('Main Fig. 2C module',45,12,5.5)
 if(!is_state){point(43,7,size=.65);txt('In tile: FDR < 0.05',45,7,5.5)}
}}
emit('SFig_principal_B_all_WGCNA',module_draw(),90,153,'Principal support','B')
emit('SFig_states_C_all_WGCNA',module_draw(TRUE),115,153,'Intermediate states','C')
draw_pattern <- function(){
 grid.newpage();texts<<-list();marks<<-list();rg<-range(scores$Score);rg<-rg+c(-.15,.15)
 for(i in 1:6){p<-paste0('P',i);x0<-10+((i-1)%%3)*60;y0<-60-floor((i-1)/3)*46;ww<-49;hh<-31;fy<-function(v)y0+(v-rg[1])/diff(rg)*hh
  txt(p,x0,y0+hh+3.4,7,'bold');txt(paste0(unique(scores[Pattern==p,Genes]),' genes'),x0+6,y0+hh+3.4,5.6,col='#68727A')
  line(c(x0,x0,x0+ww),c(y0+hh,y0,y0),'#69757D');line(c(x0,x0+ww),rep(fy(0),2),'#CED5DB',.5,2)
  for(v in pretty(rg,3))if(v>=rg[1]&v<=rg[2])txt(v,x0-1.1,fy(v),5.2,just='right')
  xcs<-x0+c(3.6,11.3,19,30,37.7,45.4);line(rep(x0+24.5,2),c(y0,y0+hh),'#E6E9EB',.4)
  for(j in 1:6){z<-scores[Pattern==p&FinalSubtype==states[j]];setorder(z,Sample);xx<-xcs[j]+((seq_len(nrow(z))*.61803398875)%%1-.5)*4.8;yy<-fy(z$Score)
   point(xx,yy,size=.55,col=adjustcolor(pal[states[j]],alpha.f=.43));q<-quantile(z$Score,c(.25,.5,.75));line(rep(xcs[j],2),fy(q[c(1,3)]),'#20272D',1);line(xcs[j]+c(-1.1,1.1),rep(fy(q[2]),2),'#20272D',1.3)
   txt(gsub('_','/',states[j]),xcs[j],y0-3,5.2,just='center');marks[[length(marks)+1L]]<<-data.table(Pattern=p,Sample=z$Sample,State=states[j],Score=z$Score,X=xx,Y=yy)
  }
 }
 txt('Pattern score (mean reference z-score)',2,54,5.7,rot=90,just='center')
 txt('Samples: C1 112   C2 223   C3 162   C1/C2 41   C1/C3 39   C2/C3 9',10,5.2,5.6)
 txt('Points: samples; black marks: median and IQR',110,5.2,5.6)
}
emit('SFig_states_B_pattern_distributions',draw_pattern,183,100,'Intermediate states','B')
draw_tme <- function(){
 grid.newpage();texts<<-list();marks<<-list();W<-183;H<-158
 studies<-c('2019Ahn','2022Dong','2023Jeon','2023Lee','2024Fan','TCGA-CHOL');abbr<-c('A','D','J','L','F','T');starts<-c(56,97,138)
 sel<-unique(tme[,.(RowID,Method,Score,Label,Block,N,Studies)]);sel[,PlotGroup:=ifelse(Block=='Context','Context',Method)];setorder(sel,RowID)
 txt('Readout',3,H-5,6.5,'bold');txt('n / studies',49,H-5,5.8,just='right')
 for(j in 1:3){s<-starts[j];txt(states[j],s+16,H-5,6.8,'bold','center');rect(s+16,H-8,37,1,pal[states[j]])
  for(k in 1:6)txt(abbr[k],s+(k-1)*4.3,H-12,5.7,just='center')
  txt('LOSO',s+31,H-12,5.5,just='center')
 }
 y<-H-18
 for(g in c('Context','MCP-counter','xCell','EPIC','CIBERSORT')){
  rect(W/2,y,W-6,3.6,'#EDF1F3');txt(if(g=='Context')'ESTIMATE context'else if(g=='CIBERSORT')'CIBERSORT  |  relative immune composition'else g,3,y,6,'bold');y<-y-4
  for(id in sel[PlotGroup==g,RowID]){r<-sel[RowID==id];if(id%%2==0)rect(W/2,y,W-6,3.1,'#FAFBFC')
   lab<-if(r$Method=='MCP-counter'&&r$Score=='Monocytic lineage')'Monocytic lineage'else r$Label
   txt(lab,3,y,5.8);txt(paste0(r$N,' / ',r$Studies),49,y,5.2,just='right',col='#65727C')
   for(j in 1:3){s<-starts[j];z<-tme[RowID==id&Subtype==states[j]]
    for(k in 1:6){q<-per[RowID==id&Subtype==states[j]&Study==studies[k]];stopifnot(nrow(q)==1);xx<-s+(k-1)*4.3
     if(!isTRUE(q$Estimable))txt('\u2013',xx,y,5.8,just='center',col='#A6AFB7')else{cc<-if(q$Effect>=0)'#B65058'else'#477DA6';solid<-q$Lower>0||q$Upper<0;point(xx,y,pch=if(q$Effect>=0)24 else 25,size=1.3,col=cc,fill=if(solid)cc else 'white')}
    }
    ok<-z$Retained==z$EstimableOmissions;txt(paste0(z$Retained,'/',z$EstimableOmissions),s+31,y,5.7,face=if(ok)'plain'else'bold',col=if(ok)'#65727C'else'#24282B',just='center')
    marks[[length(marks)+1L]]<<-data.table(RowID=id,Subtype=states[j],Retained=z$Retained,Estimable=z$EstimableOmissions,Total=z$TotalOmissions,X=s+31,Y=y)
   };y<-y-3.1
  };y<-y-1
 }
 stopifnot(y>12)
 point(4,10,24,1.3,'#B65058');txt('Higher',6,10,5.5);point(22,10,25,1.3,'#477DA6');txt('Lower',24,10,5.5)
 txt('Filled: 95% CI excludes 0; open: includes 0; dash: not estimable',43,10,5.4)
 txt('LOSO: retained / informative omissions',126,10,5.2)
 txt('A Ahn   D Dong   J Jeon   L Lee   F Fan   T TCGA',3,5,5.5)
 txt('Bold LOSO fraction: at least one direction change',99,5,5.4)
}
emit('SFig_principal_C_TME_study_robustness',draw_tme,183,158,'Principal support','C')
fwrite(rbindlist(manifest),file.path(out,'figure_manifest.csv'))
stopifnot(identical(unname(hashes),unname(tools::md5sum(src))))
fwrite(data.table(Source=src,BeforeMD5=unname(hashes),AfterMD5=unname(tools::md5sum(src))),file.path(out,'review/source_hashes.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'review/sessionInfo.txt'))
cat('Four panels rendered. All source hashes unchanged.\n')
