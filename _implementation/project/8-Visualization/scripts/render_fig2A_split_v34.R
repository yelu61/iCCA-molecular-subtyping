# Fixed v20 expression geometry; actual-member full-tested-background ORA.
source('8-Visualization/scripts/fig_style.R')
source('8-Visualization/scripts/render_bulk_figures.R')
suppressPackageStartupMessages({library(data.table);library(grid)})

render_fig2A_split_v34 <- function(){
  check_owner_gate('Step2')
  out <- '8-Visualization/refinements/fig2_assembly_20261031_v34/A'
  for(d in c('figures','tables','review'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
  paths <- c('2-TranscriptomicCharacterization/runs/fig2_soft_patterns_20261030_v20/K6/frozen_display_model.rds',
    '2-TranscriptomicCharacterization/runs/fig2_ORA_background_20261030_v21/full10531_significant_ORA.csv',
    '2-TranscriptomicCharacterization/publication/tables/step2_wgcna_module_subtype_classification.csv',
    SUBTYPE_TRUTH_CSV,SUBTYPE_PALETTE_CSV,STUDY_PALETTE_CSV,
    '8-Visualization/scripts/fig_style.R','8-Visualization/scripts/render_bulk_figures.R',
    '8-Visualization/scripts/render_fig2A_split_v34.R',
    '8-Visualization/refinements/fig2A_fold_enrichment_20261030_v23/same29/tables/annotation_evidence.csv')
  hashes <- unname(tools::md5sum(paths));stopifnot(!anyNA(hashes))
  model<-readRDS(paths[1]);ora<-fread(paths[2]);cls<-fread(paths[3])
  counts<-model$counts;rr<-copy(model$membership);columns<-copy(model$main);z<-model$z
  porder<-model$pattern_order;allpal<-load_subtype_palette()$colors;stpal<-load_study_palette()$colors
  pcols<-setNames(allpal[counts$Context],counts$Pattern)
  modcols<-c(setNames(cls$ModuleColor,cls$Module),grey='grey',Not_in_network='#FFFFFF')
  rr[,ModuleColor:=unname(modcols[WGCNAModule])];stopifnot(!anyNA(rr$ModuleColor))
  aa<-fread(paths[length(paths)])
  stopifnot(nrow(aa)==29,all(aa$UniverseGenes==10531),all(aa$FDR<.05 & aa$Count>=5))
  aa[,Score:=FoldEnrichment];stopifnot(all(is.finite(aa$Score)),all(aa$Score>0))
  aa[,SelectionReason:='The exact v23 same29 members, labels, statistics and FE order; only print styling changed']
  mm<-function(x)unit(x,'mm');font<-'Arial';termfont<-5.8;lh<-1.95;width<-120
  tmp<-tempfile(tmpdir=file.path(out,'review'),fileext='.pdf');cairo_pdf(tmp,width=120/25.4,height=20/25.4,family=font)
  wrap<-function(s,maxwidth=45){words<-strsplit(s,' +')[[1]];ans<-character();line<-''
    for(w in words){trial<-if(nzchar(line))paste(line,w)else w
      if(convertWidth(grobWidth(textGrob(trial,gp=gpar(fontfamily=font,fontsize=termfont))),'mm',valueOnly=TRUE)>maxwidth && nzchar(line)){ans<-c(ans,line);line<-w}else line<-trial}
    paste(c(ans,line),collapse='\n')}
  aa[,WrappedText:=vapply(DisplayLabel,wrap,character(1))];aa[,Lines:=lengths(strsplit(WrappedText,'\n',fixed=TRUE))]
  aa[,RowHeight:=Lines*lh+.65];dev.off();unlink(tmp)
  # Keep heatmap row area proportional to gene count; reserve readable annotation space.
  bottom<-13;gap<-.6;left<-5;heatwidth<-61
  for(height in c(123)){
    top<-height-12;plan<-copy(counts);plan[,PatternOrder:=match(Pattern,porder)];setorder(plan,PatternOrder)
    usable<-top-bottom-gap*(nrow(plan)-1);pos<-top
    for(i in seq_len(nrow(plan))){h<-usable*plan$Genes[i]/842;plan[i,`:=`(HeatTop=pos,HeatBottom=pos-h,HeatCenter=pos-h/2)];pos<-pos-h-gap}
    heights<-vapply(plan$Pattern,function(p)2.9+sum(aa[Pattern==p,RowHeight]),numeric(1))
    offsets<-c(0,cumsum(head(heights,-1)+.75))
    lower<-bottom+tail(heights,1)+tail(offsets,1)
    if(lower>top)next
    # Isotonic placement preserves order and fits all blocks inside the heat area.
    desired<-plan$HeatCenter+heights/2+offsets
    fitted<-pmax(lower,pmin(top,-isoreg(seq_along(desired),-desired)$yf))
    for(i in seq_len(nrow(plan))){p<-plan$Pattern[i];ix<-which(aa$Pattern==p);ht<-heights[i]
      at<-fitted[i]-offsets[i]
      plan[i,`:=`(HeaderY=at-.8,AnnotationTop=at,AnnotationBottom=at-ht)]
      cursor<-at-2.9
      for(j in ix){aa[j,`:=`(TextY=cursor-Lines*lh/2+.1,BarY=cursor-Lines*lh-.6,YTop=cursor,YBottom=cursor-RowHeight)];cursor<-cursor-aa$RowHeight[j]}
    }
    if(min(plan$AnnotationBottom)>=bottom-1e-8)break
  }
  stopifnot(min(plan$AnnotationBottom)>=bottom,height==123)
  states<-c('C1','C2','C3');sw<-(heatwidth-2*gap)/nrow(columns);xx<-numeric(nrow(columns));slices<-list();pos<-left
  for(st in states){ix<-which(columns$FinalSubtype==st);xx[ix]<-pos+(seq_along(ix)-1)*sw
    slices[[st]]<-data.table(State=st,Start=pos,End=pos+length(ix)*sw,N=length(ix));pos<-pos+length(ix)*sw+gap}
  slices<-rbindlist(slices);columns[,`:=`(XStart=xx,XEnd=xx+sw)]
  for(i in seq_len(nrow(plan))){ix<-which(rr$Pattern==plan$Pattern[i]);h<-(plan$HeatTop[i]-plan$HeatBottom[i])/length(ix)
    rr[ix,`:=`(YTop=plan$HeatTop[i]-(seq_along(ix)-1)*h,YBottom=plan$HeatTop[i]-seq_along(ix)*h)]}
  zcols<-colorRampPalette(c(SCALE_EFFECT_LOW,SCALE_EFFECT_MID,SCALE_EFFECT_HIGH))(256)
  zcolor<-function(a)zcols[pmax(1L,pmin(256L,1L+floor((pmax(-2,pmin(2,a))+2)/4*255)))]
  maxscore<-ceiling(max(aa$Score)/10)*10;stopifnot(maxscore==30,max(aa$Score)<=maxscore)
  texts<-list();tiles<-list();bars<-list()
  draw<-function(newpage=TRUE){if(newpage)grid.newpage();texts<<-list();tiles<<-list();bars<<-list()
    txt<-function(s,x,y,size=7,face='plain',col='grey15',just='left',lineheight=1){grid.text(s,mm(x),mm(y),just=just,gp=gpar(fontfamily=font,fontsize=size,fontface=face,col=col,lineheight=lineheight));texts[[length(texts)+1L]]<<-data.table(Text=s,X=x,Y=y,Size=size,Just=just)}
    rect<-function(x,y,w,h,c,b=NA)grid.rect(mm(x),mm(y),mm(w),mm(h),just=c('left','bottom'),gp=gpar(fill=c,col=b,lwd=.3))
    ln<-function(x,y,col='#BAC3CB')grid.lines(mm(x),mm(y),gp=gpar(col=col,lwd=.45))
    for(k in seq_len(nrow(slices))){sl<-slices[k]
      rect(sl$Start,height-5.3,sl$End-sl$Start,1.5,allpal[[sl$State]])}
    for(j in seq_len(nrow(columns)))rect(xx[j],height-8.3,sw,1.5,stpal[[columns$Study[j]]])
    txt('Subtype',68.3,height-4.55,5.5)
    txt('Study',68.3,height-7.55,5.5)
    for(i in seq_len(nrow(plan))){p<-plan$Pattern[i];gs<-rr[Pattern==p,Gene]
      yr<-4L
      for(k in seq_len(nrow(slices))){sl<-slices[k];ix<-which(columns$FinalSubtype==sl$State);a<-z[gs,columns$Sample[ix],drop=FALSE]
        colors<-matrix(zcolor(a),nrow=nrow(a),ncol=ncol(a));expanded<-colors[rep(seq_len(nrow(a)),each=yr),rep(seq_len(ncol(a)),each=5),drop=FALSE]
        grid.raster(as.raster(expanded),x=mm(sl$Start),y=mm(plan$HeatBottom[i]),width=mm(sl$End-sl$Start),height=mm(plan$HeatTop[i]-plan$HeatBottom[i]),just=c('left','bottom'),interpolate=FALSE)
        tiles[[length(tiles)+1L]]<<-data.table(Pattern=p,State=sl$State,Genes=nrow(a),Samples=ncol(a),XRepeat=5,YRepeat=yr)}
      a<-rr[Pattern==p]
      for(j in seq_len(nrow(a)))rect(2.3,a$YBottom[j],1.1,a$YTop[j]-a$YBottom[j],a$ModuleColor[j],if(a$WGCNAModule[j]=='Not_in_network')'grey60'else NA)
      ln(c(66.7,68.3,70.5),c(plan$HeatCenter[i],plan$HeatCenter[i],plan$HeaderY[i]-.3))
      txt(p,72,plan$HeaderY[i],6.1,'bold')
      a<-aa[Pattern==p]
      for(j in seq_len(nrow(a))){bw<-45*a$Score[j]/maxscore
        rect(72,a$BarY[j],bw,.45,pcols[[p]])
        txt(a$WrappedText[j],72,a$TextY[j],termfont,lineheight=lh/(termfont*25.4/72))
        bars[[length(bars)+1L]]<<-data.table(Pattern=p,Pathway=a$Pathway[j],Score=a$Score[j],WidthMM=bw,Color=pcols[[p]])}
    }
    ln(c(72,117),c(9.5,9.5),'grey55')
    for(v in c(0,maxscore/2,maxscore)){x<-72+45*v/maxscore;ln(c(x,x),c(9.5,8.9),'grey55');txt(as.character(v),x,6.9,5.1,just='center')}
    txt('Fold enrichment',94.5,3,5.8,just='center')

  }
  name<-'Figure2A_without_legend';pdf<-file.path(out,'figures',paste0(name,'.pdf'))
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
  fwrite(data.table(Figure='Figure 2',Panel='A',Role='Full29 main A body; separate legend; geometry translated by 43mm; all rows and values preserved',
    SourceTable=paste(paths[1:3],collapse=';'),TargetOutput=pdf,WidthMM=width,HeightMM=height,Status='Rendered exploratory candidate',
    QAFile=file.path(out,'review/verification.json'),ArtifactMD5=unname(tools::md5sum(pdf))),file.path(out,'figure_manifest.csv'))
  stopifnot(identical(hashes,unname(tools::md5sum(paths))))
  jsonlite::write_json(list(Genes=nrow(rr),Samples=nrow(columns),PatternGenes=as.list(setNames(plan$Genes,plan$Pattern)),
    Terms=nrow(aa),BarMetric='FoldEnrichment',BarMaximum=maxscore,AnnotationStripHeightMM=1.5,WidthMM=width,HeightMM=height,HeatmapValuesUnchanged=TRUE,MembershipUnchanged=TRUE,SourceHashesUnchanged=TRUE),
    file.path(out,'review/render_diagnostics.json'),auto_unbox=TRUE,pretty=TRUE)
  invisible(list(draw=draw,width=width,height=height,allpal=allpal,pcols=pcols,stpal=stpal,modcols=modcols,studies=names(stpal)[names(stpal)%in%columns$Study],modules=names(modcols)[names(modcols)%in%rr$WGCNAModule],zcols=zcols))
}
if(sys.nframe()==0)render_fig2A_split_v34()
