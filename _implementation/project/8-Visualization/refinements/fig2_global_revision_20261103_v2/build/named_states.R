suppressPackageStartupMessages({library(data.table);library(grid);library(svglite);library(ragg);library(jsonlite)})
source('8-Visualization/scripts/fig_style.R')
out <- '8-Visualization/refinements/fig2_global_revision_20261103_v2/named_states'
inputs <- c(model='2-TranscriptomicCharacterization/runs/fig2_soft_patterns_20261030_v20/K6/frozen_display_model.rds',
  annotations='8-Visualization/refinements/sfig2_states_20261030_v26/tables/annotation_evidence.csv',
  rows='8-Visualization/refinements/sfig2_states_20261030_v26/tables/row_order.csv',
  columns='8-Visualization/refinements/sfig2_states_20261030_v26/tables/sample_order.csv',
  classification='2-TranscriptomicCharacterization/publication/tables/step2_wgcna_module_subtype_classification.csv',
  scores='8-Visualization/refinements/sfig2_patterns_3x2_20261031_v36/tables/pattern_scores.csv',
  summaries='8-Visualization/refinements/sfig2_patterns_3x2_20261031_v36/tables/pattern_summary.csv',
  modules='8-Visualization/refinements/sfig2_compact_20261031_v35/tables/module_order.csv',
  medians='8-Visualization/refinements/sfig2_compact_20261031_v35/tables/network_state_medians.csv',
  canonical_medians='2-TranscriptomicCharacterization/runs/fig2_intermediate_20261026_v1/wgcna_state_summary.csv',
  subtype_palette=SUBTYPE_PALETTE_CSV,study_palette=STUDY_PALETTE_CSV,
  style='8-Visualization/scripts/fig_style.R',style_config='8-Visualization/scripts/figure_config.json')
hashes <- unname(tools::md5sum(inputs));stopifnot(!anyNA(hashes))
model <- readRDS(inputs['model']);aa <- fread(inputs['annotations']);oldrows<-fread(inputs['rows']);oldcols<-fread(inputs['columns'])
rr <- copy(model$membership);columns <- copy(model$all);counts <- copy(model$counts);porder<-model$pattern_order
cls<-fread(inputs['classification']);scores<-fread(inputs['scores']);ps<-fread(inputs['summaries'])
cl<-fread(inputs['modules']);state<-fread(inputs['medians']);canonical_state<-fread(inputs['canonical_medians'])
pal<-load_subtype_palette()$colors;stpal<-load_study_palette()$colors
states<-c('C1','C2','C3','C1_C2','C1_C3','C2_C3');ns<-c(112,223,162,41,39,9)
stopifnot(nrow(rr)==842,nrow(columns)==586,nrow(aa)==29,nrow(scores)==3516,nrow(ps)==36,nrow(cl)==24,nrow(state)==144,
          identical(rr$Gene,oldrows$Gene),identical(rr$Pattern,oldrows$Pattern),identical(columns$Sample,oldcols$Sample),
          identical(columns$FinalSubtype,oldcols$FinalSubtype),all(aa$FDR<.05),all(aa$Count>=5),
          all(aa$UniverseGenes==10531),all(counts$Genes==c(242,104,155,74,175,92)))
for(i in seq_along(states))stopifnot(sum(columns$FinalSubtype==states[i])==ns[i])
z <- model$z[rr$Gene,columns$Sample,drop=FALSE]
stopifnot(all(is.finite(z)),identical(dim(z),c(842L,586L)))
for(i in seq_len(nrow(ps))){s<-ps[i];v<-scores[Pattern==s$Pattern & FinalSubtype==s$FinalSubtype,Score];q<-quantile(v,c(.25,.5,.75));stopifnot(length(v)==s$N,max(abs(q-c(s$Q25,s$Median,s$Q75)))<1e-12)}
check_state<-merge(state[,.(Module,FinalSubtype,Median,N)],canonical_state[,.(Module,FinalSubtype,Median,N)],by=c('Module','FinalSubtype'))
stopifnot(nrow(check_state)==144,all(check_state$N.x==check_state$N.y),max(abs(check_state$Median.x-check_state$Median.y))<1e-12)
for(n in setdiff(names(inputs),'model'))file.copy(inputs[n],file.path(out,'source_data',paste0(n,'_',basename(inputs[n]))),overwrite=TRUE)
con<-gzfile(file.path(out,'source_data/A_expression_reference_zscore.csv.gz'),'wt');write.csv(data.frame(Gene=rownames(z),z,check.names=FALSE),con,row.names=FALSE);close(con)
fwrite(data.table(Source=unname(inputs),MD5=hashes),file.path(out,'source_data/input_hashes.csv'))
modcols<-c(setNames(cls$ModuleColor,cls$Module),grey='grey',Not_in_network='#FFFFFF')
rr[,ModuleColor:=unname(modcols[WGCNAModule])];stopifnot(!anyNA(rr$ModuleColor))
pcols<-setNames(pal[counts$Context],counts$Pattern)
zcols<-colorRampPalette(c(SCALE_EFFECT_LOW,SCALE_EFFECT_MID,SCALE_EFFECT_HIGH))(256)
zcolor<-function(a)zcols[pmax(1L,pmin(256L,1L+floor((pmax(-2,pmin(2,a))+2)/4*255)))]
mm<-function(x)unit(x,'mm');texts<-list();marks<-list();tiles<-list();bars<-list();summaries<-list()
txt<-function(s,x,y,size=6,face='plain',just='left',col='#25282B',rot=0,lineheight=1){grid.text(s,mm(x),mm(y),just=just,rot=rot,gp=gpar(fontfamily='Arial',fontsize=size,fontface=face,col=col,lineheight=lineheight));texts[[length(texts)+1L]]<<-data.table(Text=as.character(s),X=x,Y=y,Size=size)}
rect<-function(x,y,w,h,fill,border=NA)grid.rect(mm(x),mm(y),mm(w),mm(h),just=c('left','bottom'),gp=gpar(fill=fill,col=border,lwd=.3))
line<-function(x,y,col='#BAC3CB',lwd=.45,lty=1)grid.lines(mm(x),mm(y),gp=gpar(col=col,lwd=lwd,lty=lty))
point<-function(x,y,size=.55,col='#303840',pch=16)grid.points(mm(x),mm(y),pch=pch,size=mm(size),gp=gpar(col=col,lwd=.6))
clear_geometry<-function(){texts<<-list();marks<<-list();tiles<<-list();bars<<-list();summaries<<-list()}

# A: display geometry changes only; all 842 x 586 source values are retained.
plan<-copy(counts);plan[,PatternOrder:=match(Pattern,porder)];setorder(plan,PatternOrder)
top<-118;bottom<-31;gap<-.7;left<-10.4;heatwidth<-87.2
usable<-top-bottom-gap*5;cursor<-top
for(i in seq_len(nrow(plan))){h<-usable*plan$Genes[i]/842;plan[i,`:=`(HeatTop=cursor,HeatBottom=cursor-h,HeatCenter=cursor-h/2)];cursor<-cursor-h-gap}
gaps<-c(.65,.65,1.5,.65,.65);sw<-(heatwidth-sum(gaps))/586;xx<-numeric(586);sl<-list();cursor<-left
for(i in seq_along(states)){ids<-which(columns$FinalSubtype==states[i]);xx[ids]<-cursor+(seq_along(ids)-1)*sw;sl[[i]]<-data.table(State=states[i],Start=cursor,End=cursor+length(ids)*sw,N=length(ids));cursor<-cursor+length(ids)*sw+if(i<6)gaps[i]else 0}
slices<-rbindlist(sl);columns[,`:=`(XStart=xx,XEnd=xx+sw)]
for(i in seq_len(nrow(plan))){ids<-which(rr$Pattern==plan$Pattern[i]);rh<-(plan$HeatTop[i]-plan$HeatBottom[i])/length(ids);rr[ids,`:=`(YTop=plan$HeatTop[i]-(seq_along(ids)-1)*rh,YBottom=plan$HeatTop[i]-seq_along(ids)*rh)]}
# Side labels avoid six additional header rows while retaining every term.
annotation_cursor<-118.8
for(i in seq_len(nrow(plan))){ids<-which(aa$Pattern==plan$Pattern[i]);blockheight<-length(ids)*2.75;plan[i,`:=`(AnnotationTop=annotation_cursor,AnnotationBottom=annotation_cursor-blockheight,AnnotationCenter=annotation_cursor-blockheight/2)]
  for(j in seq_along(ids)){yy<-annotation_cursor-(j-1)*2.75;aa[ids[j],`:=`(TextY=yy-.85,BarY=yy-2.45)]};annotation_cursor<-annotation_cursor-blockheight-.8}
stopifnot(min(plan$AnnotationBottom)>31)
draw_a<-function(newpage=TRUE){if(newpage)grid.newpage();clear_geometry()
  for(i in 1:6){s<-slices[i];rect(s$Start,123.7,s$End-s$Start,1.5,pal[s$State])}
  for(j in seq_len(nrow(columns)))rect(xx[j],120.7,sw,1.5,stpal[columns$Study[j]])
  txt('State',9.1,124.45,5.8,just='right');txt('Study',9.1,121.45,5.8,just='right')
  for(i in seq_len(nrow(plan))){p<-plan$Pattern[i];genes<-rr[Pattern==p,Gene]
    for(k in 1:6){s<-slices[k];ids<-which(columns$FinalSubtype==s$State);values<-z[genes,columns$Sample[ids],drop=FALSE];colors<-matrix(zcolor(values),nrow=nrow(values),ncol=ncol(values));expanded<-colors[rep(seq_len(nrow(values)),each=4),rep(seq_len(ncol(values)),each=5),drop=FALSE]
      grid.raster(as.raster(expanded),x=mm(s$Start),y=mm(plan$HeatBottom[i]),width=mm(s$End-s$Start),height=mm(plan$HeatTop[i]-plan$HeatBottom[i]),just=c('left','bottom'),interpolate=FALSE)
      tiles[[length(tiles)+1L]]<<-data.table(Pattern=p,State=s$State,Genes=nrow(values),Samples=ncol(values),XRepeat=5,YRepeat=4)}
    r<-rr[Pattern==p];for(j in seq_len(nrow(r)))rect(8.2,r$YBottom[j],1.3,r$YTop[j]-r$YBottom[j],r$ModuleColor[j],if(r$WGCNAModule[j]=='Not_in_network')'grey60'else NA)
    line(c(98.2,100.0,102.1),c(plan$HeatCenter[i],plan$HeatCenter[i],plan$AnnotationCenter[i]))
    txt(p,104.1,plan$AnnotationCenter[i],6.5,'bold',just='center')
    a<-aa[Pattern==p]
    for(j in seq_len(nrow(a))){bw<-61*a$FoldEnrichment[j]/30;rect(109.5,a$BarY[j],bw,.45,pcols[p]);txt(a$DisplayLabel[j],109.5,a$TextY[j],6.1)
      bars[[length(bars)+1L]]<<-data.table(Pattern=p,Pathway=a$Pathway[j],FoldEnrichment=a$FoldEnrichment[j],FDR=a$FDR[j],BarWidth=bw,Label=a$DisplayLabel[j])}
  }
  line(c(109.5,170.5),c(28.7,28.7),'grey55');for(v in c(0,15,30)){x<-109.5+61*v/30;line(c(x,x),c(28.7,28.15),'grey55');txt(v,x,26.3,5.4,just='center')};txt('Fold enrichment',140,23.4,6.1,just='center')
  txt('State',2,20,6.1,'bold');txt('Study',45,20,6.1,'bold');txt('WGCNA',86,20,6.1,'bold')
  for(j in 1:6){x<-2+((j-1)%%3)*13.5;y<-16.6-floor((j-1)/3)*3.5;rect(x,y-.65,1.3,1.3,pal[states[j]]);txt(gsub('_','/',states[j]),x+2,y,5.6)}
  txt('Z-score',2,8.8,5.6);for(j in 1:256)rect(2+(j-1)*29/256,5.4,29/256,1.4,zcols[j]);for(v in c(-2,0,2))txt(v,2+(v+2)/4*29,2.9,5.3,just='center')
  studies<-names(stpal)[names(stpal)%in%columns$Study]
  for(j in seq_along(studies)){x<-45+((j-1)%%2)*19.7;y<-16.6-floor((j-1)/2)*3.5;rect(x,y-.65,1.3,1.3,stpal[studies[j]]);txt(studies[j],x+2,y,5.5)}
  mods<-names(modcols)[names(modcols)%in%rr$WGCNAModule]
  for(j in seq_along(mods)){x<-86+((j-1)%%5)*18;y<-16.6-floor((j-1)/5)*3.5;rect(x,y-.65,1.3,1.3,modcols[mods[j]],if(mods[j]=='Not_in_network')'grey60'else NA);txt(if(mods[j]=='Not_in_network')'Outside network'else sub('^ME','',mods[j]),x+2,y,5.4)}
}

# B: explicit 3 rows x 2 columns, all points and all original quartile summaries.
reg<-fread('8-Visualization/refinements/fig2_global_revision_20261103_v2/pattern_bridge/tables/shared_program_registry.csv')
draw_b<-function(newpage=TRUE){if(newpage)grid.newpage();clear_geometry();rg<-range(scores$Score)+c(-.15,.15)
  for(i in 1:6){p<-paste0('P',i);x0<-9.5+((i-1)%%2)*49.0;y0<-92.2-floor((i-1)/2)*38.4;ww<-38;hh<-24.8;fy<-function(v)y0+(v-rg[1])/diff(rg)*hh
    txt(p,x0,y0+hh+5.8,6.5,'bold');txt(paste0(unique(scores[Pattern==p,Genes]),' genes'),x0+6,y0+hh+5.8,5.4,col='#68727A');txt(reg[StableID==p,CompactDescription],x0,y0+hh+2.7,5.4)
    line(c(x0,x0,x0+ww),c(y0+hh,y0,y0),'#69757D');line(c(x0,x0+ww),rep(fy(0),2),'#CED5DB',.5,2)
    for(v in pretty(rg,3))if(v>=rg[1]&v<=rg[2])txt(v,x0-1.1,fy(v),5.7,just='right')
    xs<-x0+c(2.8,8.8,14.8,23.2,29.2,35.2);line(rep(x0+19.0,2),c(y0,y0+hh),'#E6E9EB',.4)
    for(j in 1:6){s<-scores[Pattern==p&FinalSubtype==states[j]];setorder(s,Sample);x<-xs[j]+((seq_len(nrow(s))*.61803398875)%%1-.5)*3.6;y<-fy(s$Score)
      point(x,y,size=.55,col=adjustcolor(pal[states[j]],alpha.f=.43));q<-ps[Pattern==p&FinalSubtype==states[j]];line(rep(xs[j],2),fy(c(q$Q25,q$Q75)),'#20272D',1);line(xs[j]+c(-.85,.85),rep(fy(q$Median),2),'#20272D',1.3)
      label<-if(j<=3)states[j]else sub('_','/\n',states[j]);txt(label,xs[j],y0-3.5,5.8,just='center',lineheight=.9)
      marks[[length(marks)+1L]]<<-data.table(Pattern=p,Sample=s$Sample,FinalSubtype=states[j],Score=s$Score,X=x,Y=y)
      summaries[[length(summaries)+1L]]<<-copy(q)
    }
  }
  txt('Pattern score (mean reference z-score)',2.1,67,5.8,rot=90,just='center')
  txt('n: C1 112  C2 223  C3 162  C1/C2 41  C1/C3 39  C2/C3 9',5,6.3,5.6)
  txt('Points: samples; black marks: median and IQR',5,2.5,5.5)
}

# C: full 24 x 6 median matrix, with the prior order, group headings and main-module markers.
groups<-c('C1_specific','C2_specific','C3_specific','SharedHighCompatible_C1_C2','SharedHighCompatible_C2_C3','Replicated_gradient_or_ambiguous','Non_specific')
group_labels<-c('C1-specific','C2-specific','C3-specific','Shared-high C1 / C2','Shared-high C2 / C3','Gradient / ambiguous','Non-specific')
ramp<-colorRamp(c(SCALE_EFFECT_LOW,'#FAFAFA',SCALE_EFFECT_HIGH));colval<-function(v)rgb(ramp(pmax(0,pmin(1,(v/2+1)/2))),maxColorValue=255)
draw_c<-function(newpage=TRUE){if(newpage)grid.newpage();clear_geometry();xs<-c(32.8,39.45,46.1,54.25,60.9,67.55);cw<-6.3
  txt('Module',3,120.2,6.3,'bold');txt('Genes',27,120.2,5.8,just='right')
  for(j in 1:6){txt(if(j<=3)states[j]else sub('_','/\n',states[j]),xs[j],120.1,6.0,'bold','center',lineheight=.9);rect(xs[j]-cw/2,115.95,cw,1,pal[states[j]]);txt(paste0('n=',ns[j]),xs[j],113.5,5.2,just='center')}
  y<-109.2
  for(g in seq_along(groups)){
    txt(group_labels[g],3,y,5.8,'bold',col='#56636D');line(c(3,70.8),rep(y-1.6,2),'#E4E8EA');y<-y-2.9
    for(i in cl[GroupOrder==g,RowOrder]){r<-cl[RowOrder==i];rect(3.5,y-1.35,1.2,2.7,r$ModuleColor);txt(r$ModuleColor,5.3,y,5.9,face=if(r$MainFigure)'bold'else'plain');txt(r$ModuleSize,27,y,5.3,just='right',col='#66717A');if(r$MainFigure)point(28.6,y,size=.65)
      for(j in 1:6){q<-state[Module==r$Module&FinalSubtype==states[j]];stopifnot(nrow(q)==1);v<-q$Median;rect(xs[j]-cw/2,y-1.35,cw,2.7,colval(v));txt(sprintf('%.2f',if(abs(v)<.005)0 else v),xs[j],y,5.2,just='center',col=if(abs(v)>1.3)'white'else'#24282B');marks[[length(marks)+1L]]<<-data.table(Module=r$Module,FinalSubtype=states[j],Median=v,N=q$N,X=xs[j],Y=y)}
      y<-y-2.9
    };y<-y-.3
  }
  stopifnot(y>15)
  txt('Median frozen score (SD)',3,11.9,5.8);for(k in 1:100)rect(3+(k-1)*24/100,7.2,24/100+.01,1.7,colval(-2+(k-1)/99*4));for(v in c(-2,0,2))txt(v,3+(v+2)/4*24,4.6,5.3,just='center')
  point(34,8.1,size=.65);txt('Main Fig. 2C module',36,8.1,5.5)
}

export<-function(name,draw,w,h,save_geometry=TRUE){base<-file.path(out,'figures',name)
  svglite(paste0(base,'.svg'),width=w/25.4,height=h/25.4,bg='white');draw();dev.off()
  cairo_pdf(paste0(base,'.pdf'),width=w/25.4,height=h/25.4,family='Arial',bg='white');draw();dev.off()
  agg_png(paste0(base,'.png'),width=w,height=h,units='mm',res=600,background='white');draw();dev.off()
  agg_png(file.path(out,'review',paste0(name,'_preview.png')),width=w,height=h,units='mm',res=300,background='white');draw();dev.off()
  if(!save_geometry)return(invisible(NULL))
  fwrite(rbindlist(texts),file.path(out,'review',paste0(name,'_texts.csv')))
  if(length(marks))fwrite(rbindlist(marks),file.path(out,'source_data',paste0(name,'_marks.csv')))
  if(length(tiles))fwrite(rbindlist(tiles),file.path(out,'source_data',paste0(name,'_tiles.csv')))
  if(length(bars))fwrite(rbindlist(bars),file.path(out,'source_data',paste0(name,'_bars.csv')))
  if(length(summaries))fwrite(rbindlist(summaries),file.path(out,'source_data',paste0(name,'_summaries.csv')))
}
export('SFigure3A_expression_patterns',draw_a,176,128)
export('SFigure3B_pattern_distributions',draw_b,100,125)
export('SFigure3C_module_medians',draw_c,73,125)
fwrite(rr,file.path(out,'source_data/A_row_geometry.csv'));fwrite(columns,file.path(out,'source_data/A_sample_geometry.csv'));fwrite(plan,file.path(out,'source_data/A_pattern_geometry.csv'))
fwrite(aa,file.path(out,'source_data/A_annotation_geometry.csv'))
slots<-data.table(Panel=c('A','B','C'),X=c(2,2,105),Y=c(139,6,6),Width=c(176,100,73),Height=c(128,125,125))
fwrite(slots,file.path(out,'source_data/layout_slots.csv'))
assembled<-function(){grid.newpage();draws<-list(A=draw_a,B=draw_b,C=draw_c)
  for(i in 1:3){s<-slots[i];pushViewport(viewport(x=mm(s$X),y=mm(s$Y),width=mm(s$Width),height=mm(s$Height),just=c('left','bottom')));draws[[s$Panel]](FALSE);popViewport()}
}
export('SFigure3_full_layout',assembled,180,270,save_geometry=FALSE)
stopifnot(identical(hashes,unname(tools::md5sum(inputs))))
write_json(list(Status='RENDERED_PENDING_VISUAL_QA',A=list(Genes=842,Samples=586,ORA=29,Tiles=36,WidthMM=176,HeightMM=128),B=list(Points=3516,Summaries=36,Rows=3,Columns=2,WidthMM=100,HeightMM=125),C=list(Modules=24,States=6,Medians=144,WidthMM=73,HeightMM=125),Page=list(WidthMM=180,HeightMM=270),Refit=FALSE,SourceHashesUnchanged=TRUE),file.path(out,'review/render_summary.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'review/sessionInfo.txt'))
