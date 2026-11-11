# Fixed v20 expression geometry; actual-member full-tested-background ORA.
source('8-Visualization/scripts/fig_style.R')
source('8-Visualization/scripts/render_bulk_figures.R')
suppressPackageStartupMessages({library(data.table);library(grid)})

render_sfig2_states_v26 <- function(){
  check_owner_gate('Step2')
  out <- '8-Visualization/refinements/sfig2_states_20261030_v26'
  for(d in c('figures','tables','review'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
  paths <- c('2-TranscriptomicCharacterization/runs/fig2_soft_patterns_20261030_v20/K6/frozen_display_model.rds',
    '2-TranscriptomicCharacterization/runs/fig2_ORA_background_20261030_v21/full10531_significant_ORA.csv',
    '2-TranscriptomicCharacterization/publication/tables/step2_wgcna_module_subtype_classification.csv',
    SUBTYPE_TRUTH_CSV,SUBTYPE_PALETTE_CSV,STUDY_PALETTE_CSV,
    '8-Visualization/scripts/fig_style.R','8-Visualization/scripts/render_bulk_figures.R',
    '8-Visualization/scripts/render_sfig2_states_v26.R',
    '2-TranscriptomicCharacterization/runs/fig2_soft_patterns_20261030_v20/K6/pattern_scores_all_states.csv',
    '2-TranscriptomicCharacterization/runs/fig2_soft_patterns_20261030_v20/K6/pattern_state_summary.csv',
    '8-Visualization/refinements/fig2A_axis_bottom_20261030_v25/tables/annotation_evidence.csv')
  hashes <- unname(tools::md5sum(paths));stopifnot(!anyNA(hashes))
  model<-readRDS(paths[1]);ora<-fread(paths[2]);cls<-fread(paths[3])
  counts<-model$counts;rr<-copy(model$membership);columns<-copy(model$all);z<-model$z
  porder<-model$pattern_order;allpal<-load_subtype_palette()$colors;stpal<-load_study_palette()$colors
  pcols<-setNames(allpal[counts$Context],counts$Pattern)
  modcols<-c(setNames(cls$ModuleColor,cls$Module),grey='grey',Not_in_network='#FFFFFF')
  rr[,ModuleColor:=unname(modcols[WGCNAModule])];stopifnot(!anyNA(rr$ModuleColor))
  aa<-fread(paths[length(paths)])
  stopifnot(nrow(aa)==29,all(aa$UniverseGenes==10531),all(aa$FDR<.05 & aa$Count>=5))
  aa[,Score:=FoldEnrichment];stopifnot(all(is.finite(aa$Score)),all(aa$Score>0))
  aa[,SelectionReason:='The exact v23 same29 members, labels, statistics and FE order; only print styling changed']
  mm<-function(x)unit(x,'mm');font<-'Arial';termfont<-6.1;lh<-2.2;width<-183
  tmp<-tempfile(tmpdir=file.path(out,'review'),fileext='.pdf');cairo_pdf(tmp,width=183/25.4,height=20/25.4,family=font)
  wrap<-function(s,maxwidth=61){words<-strsplit(s,' +')[[1]];ans<-character();line<-''
    for(w in words){trial<-if(nzchar(line))paste(line,w)else w
      if(convertWidth(grobWidth(textGrob(trial,gp=gpar(fontfamily=font,fontsize=termfont))),'mm',valueOnly=TRUE)>maxwidth && nzchar(line)){ans<-c(ans,line);line<-w}else line<-trial}
    paste(c(ans,line),collapse='\n')}
  aa[,WrappedText:=vapply(DisplayLabel,wrap,character(1))];aa[,Lines:=lengths(strsplit(WrappedText,'\n',fixed=TRUE))]
  aa[,RowHeight:=Lines*lh+1.05];dev.off();unlink(tmp)
  # Keep heatmap row area proportional to gene count; reserve readable annotation space.
  bottom<-33;gap<-.7;left<-12;heatwidth<-96
  for(height in c(170)){
    top<-height-14;plan<-copy(counts);plan[,PatternOrder:=match(Pattern,porder)];setorder(plan,PatternOrder)
    usable<-top-bottom-gap*(nrow(plan)-1);pos<-top
    for(i in seq_len(nrow(plan))){h<-usable*plan$Genes[i]/842;plan[i,`:=`(HeatTop=pos,HeatBottom=pos-h,HeatCenter=pos-h/2)];pos<-pos-h-gap}
    heights<-vapply(plan$Pattern,function(p)3.3+sum(aa[Pattern==p,RowHeight]),numeric(1))
    offsets<-c(0,cumsum(head(heights,-1)+1.0))
    lower<-bottom+tail(heights,1)+tail(offsets,1)
    if(lower>top)next
    # Isotonic placement preserves order and fits all blocks inside the heat area.
    desired<-plan$HeatCenter+heights/2+offsets
    fitted<-pmax(lower,pmin(top,-isoreg(seq_along(desired),-desired)$yf))
    for(i in seq_len(nrow(plan))){p<-plan$Pattern[i];ix<-which(aa$Pattern==p);ht<-heights[i]
      at<-fitted[i]-offsets[i]
      plan[i,`:=`(HeaderY=at-.8,AnnotationTop=at,AnnotationBottom=at-ht)]
      cursor<-at-3.3
      for(j in ix){aa[j,`:=`(TextY=cursor-Lines*lh/2+.1,BarY=cursor-Lines*lh-.65,YTop=cursor,YBottom=cursor-RowHeight)];cursor<-cursor-aa$RowHeight[j]}
    }
    if(min(plan$AnnotationBottom)>=bottom)break
  }
  stopifnot(min(plan$AnnotationBottom)>=bottom,height==170)
  # Reuse the space freed by the old top axis; reserve seven mm below P6.
  plan[,c('HeaderY','AnnotationTop','AnnotationBottom'):=lapply(.SD,function(x)x+7),.SDcols=c('HeaderY','AnnotationTop','AnnotationBottom')]
  aa[,c('TextY','BarY','YTop','YBottom'):=lapply(.SD,function(x)x+7),.SDcols=c('TextY','BarY','YTop','YBottom')]
  states<-c('C1','C2','C3','C1_C2','C1_C3','C2_C3');gaps<-rep(gap,5);gaps[3]<-1.5;sw<-(heatwidth-sum(gaps))/nrow(columns);xx<-numeric(nrow(columns));slices<-list();pos<-left
  for(st in states){ix<-which(columns$FinalSubtype==st);xx[ix]<-pos+(seq_along(ix)-1)*sw
    slices[[st]]<-data.table(State=st,Start=pos,End=pos+length(ix)*sw,N=length(ix));pos<-pos+length(ix)*sw+if(match(st,states)<length(states))gaps[match(st,states)]else 0}
  slices<-rbindlist(slices);columns[,`:=`(XStart=xx,XEnd=xx+sw)]
  for(i in seq_len(nrow(plan))){ix<-which(rr$Pattern==plan$Pattern[i]);h<-(plan$HeatTop[i]-plan$HeatBottom[i])/length(ix)
    rr[ix,`:=`(YTop=plan$HeatTop[i]-(seq_along(ix)-1)*h,YBottom=plan$HeatTop[i]-seq_along(ix)*h)]}
  zcols<-colorRampPalette(c(SCALE_EFFECT_LOW,SCALE_EFFECT_MID,SCALE_EFFECT_HIGH))(256)
  zcolor<-function(a)zcols[pmax(1L,pmin(256L,1L+floor((pmax(-2,pmin(2,a))+2)/4*255)))]
  maxscore<-ceiling(max(aa$Score)/10)*10;stopifnot(maxscore==30,max(aa$Score)<=maxscore)
  texts<-list();tiles<-list();bars<-list()
  draw<-function(){grid.newpage();texts<<-list();tiles<<-list();bars<<-list()
    txt<-function(s,x,y,size=7,face='plain',col='grey15',just='left',lineheight=1){grid.text(s,mm(x),mm(y),just=just,gp=gpar(fontfamily=font,fontsize=size,fontface=face,col=col,lineheight=lineheight));texts[[length(texts)+1L]]<<-data.table(Text=s,X=x,Y=y,Size=size,Just=just)}
    rect<-function(x,y,w,h,c,b=NA)grid.rect(mm(x),mm(y),mm(w),mm(h),just=c('left','bottom'),gp=gpar(fill=c,col=b,lwd=.3))
    ln<-function(x,y,col='#BAC3CB')grid.lines(mm(x),mm(y),gp=gpar(col=col,lwd=.45))
    for(k in seq_len(nrow(slices))){sl<-slices[k]
      rect(sl$Start,height-5.3,sl$End-sl$Start,1.5,allpal[[sl$State]])}
    for(j in seq_len(nrow(columns)))rect(xx[j],height-8.3,sw,1.5,stpal[[columns$Study[j]]])
    txt('Subtype',10.5,height-4.55,5.8,just='right')
    txt('Study',10.5,height-7.55,5.8,just='right')
    for(i in seq_len(nrow(plan))){p<-plan$Pattern[i];gs<-rr[Pattern==p,Gene]
      yr<-4L
      for(k in seq_len(nrow(slices))){sl<-slices[k];ix<-which(columns$FinalSubtype==sl$State);a<-z[gs,columns$Sample[ix],drop=FALSE]
        colors<-matrix(zcolor(a),nrow=nrow(a),ncol=ncol(a));expanded<-colors[rep(seq_len(nrow(a)),each=yr),rep(seq_len(ncol(a)),each=5),drop=FALSE]
        grid.raster(as.raster(expanded),x=mm(sl$Start),y=mm(plan$HeatBottom[i]),width=mm(sl$End-sl$Start),height=mm(plan$HeatTop[i]-plan$HeatBottom[i]),just=c('left','bottom'),interpolate=FALSE)
        tiles[[length(tiles)+1L]]<<-data.table(Pattern=p,State=sl$State,Genes=nrow(a),Samples=ncol(a),XRepeat=5,YRepeat=yr)}
      a<-rr[Pattern==p]
      for(j in seq_len(nrow(a)))rect(8.5,a$YBottom[j],1.4,a$YTop[j]-a$YBottom[j],a$ModuleColor[j],if(a$WGCNAModule[j]=='Not_in_network')'grey60'else NA)
      ln(c(108.7,112,115.5),c(plan$HeatCenter[i],plan$HeatCenter[i],plan$HeaderY[i]-.3))
      txt(p,117,plan$HeaderY[i],6.5,'bold')
      a<-aa[Pattern==p]
      for(j in seq_len(nrow(a))){bw<-61*a$Score[j]/maxscore
        rect(117,a$BarY[j],bw,.55,pcols[[p]])
        txt(a$WrappedText[j],117,a$TextY[j],termfont,lineheight=lh/(termfont*25.4/72))
        bars[[length(bars)+1L]]<<-data.table(Pattern=p,Pathway=a$Pathway[j],Score=a$Score[j],WidthMM=bw,Color=pcols[[p]])}
    }
    ln(c(117,178),c(37.5,37.5),'grey55')
    for(v in c(0,maxscore/2,maxscore)){x<-117+61*v/maxscore;ln(c(x,x),c(37.5,36.9),'grey55');txt(as.character(v),x,35,5.4,just='center')}
    txt('Fold enrichment',147.5,31.8,6.1,just='center')
    # Compact data keys; counts, selection and test details are in the caption.
    txt('Subtype',8,27.4,6.1,'bold');txt('Study',49,27.4,6.1,'bold');txt('WGCNA',98,27.4,6.1,'bold')
    for(j in seq_along(states)){st<-states[j];y<-23.3-((j-1)%%3)*3.25;x<-8+floor((j-1)/3)*18
      rect(x,y-.65,1.3,1.3,allpal[[st]]);txt(gsub('_','/',st),x+2.3,y,5.6)}
    txt('Z-score',8,8.7,5.6)
    for(j in 1:256)rect(8+(j-1)*29/256,5.6,29/256,1.4,zcols[j])
    for(v in c(-2,0,2))txt(as.character(v),8+(v+2)/4*29,3.2,5.3,just='center')
    studies<-names(stpal)[names(stpal) %in% columns$Study]
    for(j in seq_along(studies)){y<-23.3-(j-1)*3.25;rect(49,y-.65,1.3,1.3,stpal[[studies[j]]]);txt(studies[j],51.3,y,5.5)}
    mods<-names(modcols)[names(modcols) %in% rr$WGCNAModule]
    for(j in seq_along(mods)){x<-98+((j-1)%%4)*20;y<-23.3-floor((j-1)/4)*3.25
      rect(x,y-.65,1.3,1.3,modcols[[mods[j]]],if(mods[j]=='Not_in_network')'grey60'else NA)
      label<-if(mods[j]=='Not_in_network')'Outside network'else sub('^ME','',mods[j]);txt(label,x+2,y,5.4)}

  }
  name<-'SFigure2_states_A_expression_patterns_586';pdf<-file.path(out,'figures',paste0(name,'.pdf'))
  save_panel(draw,pdf,width,height)
  ragg::agg_png(file.path(out,'figures',paste0(name,'.png')),width=width,height=height,units='mm',res=600);draw();dev.off()
  system2('/opt/homebrew/bin/pdftoppm',c('-scale-to','1900','-singlefile','-png',shQuote(pdf),shQuote(file.path(out,'review',paste0(name,'_preview')))))
  fwrite(aa,file.path(out,'tables/annotation_evidence.csv'));fwrite(plan,file.path(out,'tables/pattern_layout.csv'))
  fwrite(rr,file.path(out,'tables/row_order.csv'));fwrite(columns,file.path(out,'tables/sample_order.csv'))
  fwrite(rbindlist(list(
    slices[,.(Kind='Subtype',Label=State,X=Start,Width=End-Start,Y=height-5.3,Height=1.5,Color=unname(allpal[State]))],
    columns[,.(Kind='Study',Label=Study,X=XStart,Width=XEnd-XStart,Y=height-8.3,Height=1.5,Color=unname(stpal[Study]))])),
    file.path(out,'tables/top_annotation_layout.csv'))
  fwrite(rbindlist(texts),file.path(out,'tables/text_layout.csv'));fwrite(rbindlist(tiles),file.path(out,'tables/matrix_tiles.csv'))
  fwrite(rbindlist(bars),file.path(out,'tables/bar_geometry.csv'))
  fwrite(data.table(Source=paths,MD5=hashes),file.path(out,'review/input_provenance.csv'))
  fwrite(data.table(Figure='SFigure 2 (Intermediate states)',Panel='A',Role='Frozen principal-derived expression patterns in all586 samples; matched Figure2A v25 annotations',
    SourceTable=paste(paths[1:3],collapse=';'),TargetOutput=pdf,WidthMM=width,HeightMM=height,Status='Rendered exploratory candidate',
    QAFile=file.path(out,'review/verification.json'),ArtifactMD5=unname(tools::md5sum(pdf))),file.path(out,'figure_manifest.csv'))
  fwrite(columns[,.(N=.N),by=.(FinalSubtype,SubtypeClass)],file.path(out,'tables/state_counts.csv'))
  stopifnot(identical(hashes,unname(tools::md5sum(paths))))
  jsonlite::write_json(list(Genes=nrow(rr),Samples=nrow(columns),PatternGenes=as.list(setNames(plan$Genes,plan$Pattern)),
    Terms=nrow(aa),BarMetric='FoldEnrichment',BarMaximum=maxscore,AnnotationStripHeightMM=1.5,WidthMM=width,HeightMM=height,HeatmapValuesUnchanged=TRUE,MembershipUnchanged=TRUE,SourceHashesUnchanged=TRUE),
    file.path(out,'review/render_diagnostics.json'),auto_unbox=TRUE,pretty=TRUE)
}

render_state_summaries_v26<-function(){
  out<-'8-Visualization/refinements/sfig2_states_20261030_v26'
  run<-'2-TranscriptomicCharacterization/runs/fig2_soft_patterns_20261030_v20/K6'
  sourcefiles<-file.path(run,c('pattern_scores_all_states.csv','pattern_state_summary.csv'))
  oldpages<-fread('8-Visualization/refinements/fig2_soft_patterns_20261030_v20/K6/supplement_bundle_index.csv')[Page>=4]
  allsources<-c(sourcefiles,oldpages$SourcePDF);hashes<-unname(tools::md5sum(allsources))
  scores<-fread(sourcefiles[1]);summary<-fread(sourcefiles[2]);pal<-load_subtype_palette()$colors
  states<-c('C1','C2','C3','C1_C2','C1_C3','C2_C3');patterns<-paste0('P',1:6)
  stopifnot(nrow(scores)==3516,nrow(summary)==36,all(abs(summary$Mean)<1.5))
  mm<-function(x)unit(x,'mm');font<-'Arial'
  txt<-function(s,x,y,size=6,face='plain',just='left',col='grey15',rot=0)grid.text(s,mm(x),mm(y),just=just,rot=rot,gp=gpar(fontfamily=font,fontsize=size,fontface=face,col=col))
  rect<-function(x,y,w,h,color)grid.rect(mm(x),mm(y),mm(w),mm(h),just=c('left','bottom'),gp=gpar(fill=color,col=NA))
  ramp<-colorRampPalette(c(SCALE_EFFECT_LOW,SCALE_EFFECT_MID,SCALE_EFFECT_HIGH))(256)
  zcolor<-function(v)ramp[pmax(1L,pmin(256L,1L+floor((pmax(-1.5,pmin(1.5,v))+1.5)/3*255)))]
  exports<-list()
  export<-function(draw,name,h,panel){pdf<-file.path(out,'figures',paste0(name,'.pdf'));save_panel(draw,pdf,183,h)
    ragg::agg_png(file.path(out,'figures',paste0(name,'.png')),width=183,height=h,units='mm',res=600);draw();dev.off()
    system2('/opt/homebrew/bin/pdftoppm',c('-scale-to','1600','-singlefile','-png',shQuote(pdf),shQuote(file.path(out,'review',paste0(name,'_preview')))))
    exports[[panel]]<<-data.table(Figure='SFigure 2 (Intermediate states)',Panel=panel,Role=if(panel=='B')'Fixed-pattern state means'else 'Fixed-pattern per-sample distributions',
      SourceTable=paste(sourcefiles,collapse=';'),TargetOutput=pdf,WidthMM=183,HeightMM=h,Status='Rendered descriptive frozen-reference scores',QAFile=file.path(out,'review/verification.json'),ArtifactMD5=unname(tools::md5sum(pdf)))
    pdf}
  drawmeans<-function(){grid.newpage();x0<-25;cw<-25.5;top<-66;rh<-7.4
    for(j in seq_along(states)){st<-states[j];x<-x0+(j-.5)*cw
      txt(gsub('_','/',st),x,72,6.2,'bold',just='center');rect(x0+(j-1)*cw,68.2,cw-.35,.9,pal[[st]])}
    for(i in seq_along(patterns)){p<-patterns[i];y<-top-(i-1)*rh;txt(p,22,y-rh/2,6.3,'bold',just='right')
      for(j in seq_along(states)){v<-summary[Pattern==p & FinalSubtype==states[j],Mean];color<-zcolor(v)
        rect(x0+(j-1)*cw,y-rh,cw-.35,rh-.35,color)
        rgb<-grDevices::col2rgb(color)/255;fg<-if(sum(rgb*c(.2126,.7152,.0722))<.52)'white'else 'grey10'
        txt(sprintf('%.2f',v),x0+(j-.5)*cw,y-rh/2,6,just='center',col=fg)}}
    txt('Mean z-score',25,13,6)
    for(k in 1:256)rect(25+(k-1)*40/256,8.5,40/256,1.6,ramp[k])
    for(v in c(-1.5,0,1.5))txt(as.character(v),25+(v+1.5)/3*40,5.5,5.5,just='center')}
  meanpdf<-export(drawmeans,'SFigure2_states_B_pattern_means',82,'B')
  points<-list();iqr<-list();ranges<-list()
  drawdist<-function(){grid.newpage();points<<-list();iqr<<-list();ranges<<-list()
    txt('Mean z-score',2,55,6,just='center',rot=90)
    for(i in seq_along(patterns)){p<-patterns[i];a<-scores[Pattern==p];x0<-11+((i-1)%%3)*60;y0<-68-floor((i-1)/3)*50;ww<-47.7;hh<-32
      rg<-range(a$Score);padding<-max(.1,diff(rg)*.08);rg<-rg+c(-padding,padding);fy<-function(v)y0+(v-rg[1])/diff(rg)*hh
      ranges[[i]]<<-data.table(Pattern=p,YMin=rg[1],YMax=rg[2],X=x0,Y=y0,Width=ww,Height=hh)
      txt(p,x0,y0+hh+3,6.5,'bold')
      grid.lines(mm(c(x0,x0,x0+ww)),mm(c(y0+hh,y0,y0)),gp=gpar(col='grey45',lwd=.4))
      if(rg[1]<0 && rg[2]>0)grid.lines(mm(c(x0,x0+ww)),mm(rep(fy(0),2)),gp=gpar(col='grey88',lwd=.35))
      for(v in pretty(rg,3))if(v>=rg[1] && v<=rg[2]){grid.lines(mm(c(x0-.5,x0)),mm(rep(fy(v),2)),gp=gpar(col='grey45',lwd=.4));txt(format(v,digits=2),x0-1,fy(v),5.2,just='right')}
      for(j in seq_along(states)){st<-states[j];d<-a[FinalSubtype==st];x<-x0+(j-.5)*ww/6;v<-d$Score
        jitter<-((seq_along(v)*.61803398875)%%1-.5)*(ww/6*.64)
        grid.points(mm(x+jitter),mm(fy(v)),pch=16,size=unit(.32,'mm'),gp=gpar(col=grDevices::adjustcolor(pal[[st]],alpha.f=.35)))
        q<-quantile(v,c(.25,.5,.75));grid.lines(mm(c(x,x)),mm(fy(q[c(1,3)])),gp=gpar(col='grey20',lwd=.65));grid.lines(mm(x+c(-.8,.8)),mm(rep(fy(q[2]),2)),gp=gpar(col='grey10',lwd=1))
        txt(gsub('_','/',st),x,y0-3.5,5.3,just='center')
        points[[length(points)+1L]]<<-data.table(Pattern=p,Sample=d$Sample,FinalSubtype=st,Score=v,X=x+jitter,Y=fy(v))
        iqr[[length(iqr)+1L]]<<-data.table(Pattern=p,FinalSubtype=st,N=length(v),Q25=q[1],Median=q[2],Q75=q[3])}}
  }
  distpdf<-export(drawdist,'SFigure2_states_C_pattern_distributions',110,'C')
  fwrite(rbindlist(points),file.path(out,'tables/distribution_points.csv'));fwrite(rbindlist(iqr),file.path(out,'tables/distribution_quartiles.csv'));fwrite(rbindlist(ranges),file.path(out,'tables/distribution_ranges.csv'))
  file.copy(sourcefiles[1],file.path(out,'tables/pattern_scores_all_states.csv'),overwrite=TRUE)
  file.copy(sourcefiles[2],file.path(out,'tables/pattern_state_summary.csv'),overwrite=TRUE)
  manifest<-rbindlist(c(list(fread(file.path(out,'figure_manifest.csv'))),exports),use.names=TRUE,fill=TRUE)
  fwrite(manifest,file.path(out,'figure_manifest.csv'))
  bundlefiles<-c(file.path(out,'figures/SFigure2_states_A_expression_patterns_586.pdf'),meanpdf,distpdf,oldpages$SourcePDF)
  system2('/opt/homebrew/bin/pdfunite',c(vapply(bundlefiles,shQuote,character(1)),shQuote(file.path(out,'figures/SFigure2_intermediate_states_review_bundle.pdf'))))
  fwrite(data.table(Page=seq_along(bundlefiles),SourcePDF=bundlefiles,SourceMD5=unname(tools::md5sum(bundlefiles)),Role=c('Updated fixed842 heatmap across586','Frozen pattern state means','Frozen pattern per-sample distributions',oldpages$Role)),file.path(out,'supplement_bundle_index.csv'))
  fwrite(data.table(Source=allsources,MD5=hashes),file.path(out,'review/summary_input_provenance.csv'))
  stopifnot(identical(hashes,unname(tools::md5sum(allsources))))
}
if(sys.nframe()==0){render_sfig2_states_v26();render_state_summaries_v26()}

