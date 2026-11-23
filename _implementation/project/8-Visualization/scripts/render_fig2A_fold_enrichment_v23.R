# Fixed v20 expression geometry; actual-member full-tested-background ORA.
source('8-Visualization/scripts/fig_style.R')
source('8-Visualization/scripts/render_bulk_figures.R')
suppressPackageStartupMessages({library(data.table);library(grid)})

render_fig2A_fold_enrichment_v23 <- function(variant='same29'){
  check_owner_gate('Step2')
  out <- file.path('8-Visualization/refinements/fig2A_fold_enrichment_20261030_v23',variant)
  for(d in c('figures','tables','review'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
  paths <- c('2-TranscriptomicCharacterization/runs/fig2_soft_patterns_20261030_v20/K6/frozen_display_model.rds',
    '2-TranscriptomicCharacterization/runs/fig2_ORA_background_20261030_v21/full10531_significant_ORA.csv',
    '2-TranscriptomicCharacterization/publication/tables/step2_wgcna_module_subtype_classification.csv',
    SUBTYPE_TRUTH_CSV,SUBTYPE_PALETTE_CSV,STUDY_PALETTE_CSV,
    '8-Visualization/scripts/fig_style.R','8-Visualization/scripts/render_bulk_figures.R',
    '8-Visualization/scripts/render_fig2A_fold_enrichment_v23.R',
    '8-Visualization/refinements/fig2A_full_background_20261030_v22/tables/annotation_evidence.csv')
  hashes <- unname(tools::md5sum(paths));stopifnot(!anyNA(hashes))
  model<-readRDS(paths[1]);ora<-fread(paths[2]);cls<-fread(paths[3])
  counts<-model$counts;rr<-copy(model$membership);columns<-copy(model$main);z<-model$z
  porder<-model$pattern_order;allpal<-load_subtype_palette()$colors;stpal<-load_study_palette()$colors
  pcols<-setNames(allpal[counts$Context],counts$Pattern)
  modcols<-c(setNames(cls$ModuleColor,cls$Module),grey='grey',Not_in_network='#FFFFFF')
  rr[,ModuleColor:=unname(modcols[WGCNAModule])];stopifnot(!anyNA(rr$ModuleColor))
  # All candidates already passed the same FDR/hit/size rules. No ORA rerun.
  stopifnot(all(ora$SelectedGenes>0),all(ora$SetSize>0),all(ora$UniverseGenes==10531))
  ora[,RecomputedFE:=(Count/SelectedGenes)/(SetSize/UniverseGenes)]
  stopifnot(max(abs(ora$RecomputedFE-ora$FoldEnrichment))<1e-10)
  setorder(ora,Pattern,-FoldEnrichment,FDR,-Count,Pathway)
  ora[,FERank:=seq_len(.N),by=Pattern]
  fwrite(ora,file.path(out,'tables/all_significant_FE_ranked.csv'))
  fwrite(ora[,head(.SD,5),by=Pattern],file.path(out,'tables/strict_top5_FE.csv'))
  oldann<-fread(paths[length(paths)])
  if(variant=='same29'){
    aa<-merge(oldann[,.(Pattern,Pathway,DisplayLabel)],ora,by=c('Pattern','Pathway'),sort=FALSE)
    stopifnot(nrow(aa)==29)
    reason<-'Exactly the existing29 representatives, reordered by descending fold enrichment'
  } else {
    selected<-list();skipped<-list()
    for(p in porder){candidate<-ora[Pattern==p];chosen<-list();indices<-integer()
      for(i in seq_len(nrow(candidate))){g<-strsplit(candidate$HitGenes[i],'/',fixed=TRUE)[[1]]
        jac<-if(length(chosen))max(vapply(chosen,function(x)length(intersect(g,x))/length(union(g,x)),numeric(1)))else 0
        if(jac<.5){indices<-c(indices,i);chosen[[length(chosen)+1L]]<-g}else{
          skipped[[length(skipped)+1L]]<-candidate[i,.(Pattern,Pathway,FoldEnrichment,FERank,MaxHitJaccard=jac)]}
        if(length(indices)==5L)break}
      selected[[p]]<-candidate[indices]
    }
    aa<-rbindlist(selected)
    labels<-unique(oldann[,.(Pathway,DisplayLabel)])
    aa<-merge(aa,labels,by='Pathway',all.x=TRUE,sort=FALSE)
    pretty_term<-function(s){s<-sub('^(GOBP_|HALLMARK_|REACTOME_)','',s);s<-sub(' - Homo sapiens \\(human\\)$','',s)
      s<-gsub('_',' ',s);if(s==toupper(s))s<-tolower(s);paste0(toupper(substr(s,1,1)),substr(s,2,nchar(s)))}
    aa[is.na(DisplayLabel),DisplayLabel:=vapply(Pathway,pretty_term,character(1))]
    aa[Pathway=='REACTOME_DIFFERENTIATION_OF_KERATINOCYTES_IN_INTERFOLLICULAR_EPIDERMIS_IN_MAMMALIAN_SKIN',DisplayLabel:='Interfollicular keratinocyte differentiation']
    fwrite(rbindlist(skipped),file.path(out,'tables/redundant_terms_skipped.csv'))
    reason<-'All significant candidates ranked by FE; greedy hit-set Jaccard<0.5, up to5 per pattern'
  }
  aa[,PatternOrder:=match(Pattern,porder)];setorder(aa,PatternOrder,-FoldEnrichment,FDR,-Count,Pathway)
  aa[,DisplayOrder:=seq_len(.N),by=Pattern];aa[,SelectionReason:=reason]
  aa[,Score:=FoldEnrichment];stopifnot(all(is.finite(aa$Score)),all(aa$Score>0))
  mm<-function(x)unit(x,'mm');font<-'Arial';termfont<-6.6;lh<-2.45;width<-174
  tmp<-tempfile(tmpdir=file.path(out,'review'),fileext='.pdf');cairo_pdf(tmp,width=174/25.4,height=20/25.4,family=font)
  wrap<-function(s,maxwidth=57){words<-strsplit(s,' +')[[1]];ans<-character();line<-''
    for(w in words){trial<-if(nzchar(line))paste(line,w)else w
      if(convertWidth(grobWidth(textGrob(trial,gp=gpar(fontfamily=font,fontsize=termfont))),'mm',valueOnly=TRUE)>maxwidth && nzchar(line)){ans<-c(ans,line);line<-w}else line<-trial}
    paste(c(ans,line),collapse='\n')}
  aa[,WrappedText:=vapply(DisplayLabel,wrap,character(1))];aa[,Lines:=lengths(strsplit(WrappedText,'\n',fixed=TRUE))]
  aa[,RowHeight:=Lines*lh+1.25];dev.off();unlink(tmp)
  # Keep heatmap row area proportional to gene count; reserve readable annotation space.
  bottom<-64;gap<-.7;left<-12;heatwidth<-91
  for(height in seq(212,280,2)){
    top<-height-29;plan<-copy(counts);plan[,PatternOrder:=match(Pattern,porder)];setorder(plan,PatternOrder)
    usable<-top-bottom-gap*(nrow(plan)-1);pos<-top
    for(i in seq_len(nrow(plan))){h<-usable*plan$Genes[i]/842;plan[i,`:=`(HeatTop=pos,HeatBottom=pos-h,HeatCenter=pos-h/2)];pos<-pos-h-gap}
    heights<-vapply(plan$Pattern,function(p)4.2+sum(aa[Pattern==p,RowHeight]),numeric(1))
    offsets<-c(0,cumsum(head(heights,-1)+1.4))
    lower<-bottom+tail(heights,1)+tail(offsets,1)
    if(lower>top)next
    # Isotonic placement preserves order and fits all blocks inside the heat area.
    desired<-plan$HeatCenter+heights/2+offsets
    fitted<-pmax(lower,pmin(top,-isoreg(seq_along(desired),-desired)$yf))
    for(i in seq_len(nrow(plan))){p<-plan$Pattern[i];ix<-which(aa$Pattern==p);ht<-heights[i]
      at<-fitted[i]-offsets[i]
      plan[i,`:=`(HeaderY=at-1,AnnotationTop=at,AnnotationBottom=at-ht)]
      cursor<-at-4.2
      for(j in ix){aa[j,`:=`(TextY=cursor-Lines*lh/2+.1,BarY=cursor-Lines*lh-.65,YTop=cursor,YBottom=cursor-RowHeight)];cursor<-cursor-aa$RowHeight[j]}
    }
    if(min(plan$AnnotationBottom)>=bottom)break
  }
  stopifnot(min(plan$AnnotationBottom)>=bottom,height<=280)
  states<-c('C1','C2','C3');sw<-(heatwidth-2*gap)/nrow(columns);xx<-numeric(nrow(columns));slices<-list();pos<-left
  for(st in states){ix<-which(columns$FinalSubtype==st);xx[ix]<-pos+(seq_along(ix)-1)*sw
    slices[[st]]<-data.table(State=st,Start=pos,End=pos+length(ix)*sw,N=length(ix));pos<-pos+length(ix)*sw+gap}
  slices<-rbindlist(slices);columns[,`:=`(XStart=xx,XEnd=xx+sw)]
  for(i in seq_len(nrow(plan))){ix<-which(rr$Pattern==plan$Pattern[i]);h<-(plan$HeatTop[i]-plan$HeatBottom[i])/length(ix)
    rr[ix,`:=`(YTop=plan$HeatTop[i]-(seq_along(ix)-1)*h,YBottom=plan$HeatTop[i]-seq_along(ix)*h)]}
  zcols<-colorRampPalette(c(SCALE_EFFECT_LOW,SCALE_EFFECT_MID,SCALE_EFFECT_HIGH))(256)
  zcolor<-function(a)zcols[pmax(1L,pmin(256L,1L+floor((pmax(-2,pmin(2,a))+2)/4*255)))]
  maxscore<-50;stopifnot(max(aa$Score)<=maxscore)
  texts<-list();tiles<-list();bars<-list()
  draw<-function(){grid.newpage();texts<<-list();tiles<<-list();bars<<-list()
    txt<-function(s,x,y,size=7,face='plain',col='grey15',just='left',lineheight=1){grid.text(s,mm(x),mm(y),just=just,gp=gpar(fontfamily=font,fontsize=size,fontface=face,col=col,lineheight=lineheight));texts[[length(texts)+1L]]<<-data.table(Text=s,X=x,Y=y,Size=size,Just=just)}
    rect<-function(x,y,w,h,c,b=NA)grid.rect(mm(x),mm(y),mm(w),mm(h),just=c('left','bottom'),gp=gpar(fill=c,col=b,lwd=.3))
    ln<-function(x,y,col='#BAC3CB')grid.lines(mm(x),mm(y),gp=gpar(col=col,lwd=.45))
    txt('Principal-state expression patterns',4,height-5,9,'bold')
    txt('842 genes   |   497 samples   |   fixed principal-reference scaling',12,height-11,7,col='grey35')
    for(k in seq_len(nrow(slices))){sl<-slices[k]
      rect(sl$Start,height-22.5,sl$End-sl$Start,1.5,allpal[[sl$State]])}
    for(j in seq_len(nrow(columns)))rect(xx[j],height-25.5,sw,1.5,stpal[[columns$Study[j]]])
    txt('Subtype',10.5,height-21.75,6.1,just='right')
    txt('Study',10.5,height-24.75,6.1,just='right')
    txt('Pattern-member enrichment',112,height-15,7.3,'bold')
    txt('Fold enrichment',112,height-19.5,6.2,col='grey40')
    ln(c(112,169),rep(height-22.5,2),'grey60')
    for(v in c(0,10,20,30,40,50)){x<-112+57*v/maxscore;ln(c(x,x),c(height-22.5,height-23.1),'grey60');txt(as.character(v),x,height-25.6,5.6,just='center')}
    for(i in seq_len(nrow(plan))){p<-plan$Pattern[i];gs<-rr[Pattern==p,Gene]
      yr<-max(4L,as.integer(ceiling(usable/842/25.4*600)))
      for(k in seq_len(nrow(slices))){sl<-slices[k];ix<-which(columns$FinalSubtype==sl$State);a<-z[gs,columns$Sample[ix],drop=FALSE]
        colors<-matrix(zcolor(a),nrow=nrow(a),ncol=ncol(a));expanded<-colors[rep(seq_len(nrow(a)),each=yr),rep(seq_len(ncol(a)),each=5),drop=FALSE]
        grid.raster(as.raster(expanded),x=mm(sl$Start),y=mm(plan$HeatBottom[i]),width=mm(sl$End-sl$Start),height=mm(plan$HeatTop[i]-plan$HeatBottom[i]),just=c('left','bottom'),interpolate=FALSE)
        tiles[[length(tiles)+1L]]<<-data.table(Pattern=p,State=sl$State,Genes=nrow(a),Samples=ncol(a),XRepeat=5,YRepeat=yr)}
      a<-rr[Pattern==p]
      for(j in seq_len(nrow(a)))rect(8.5,a$YBottom[j],1.4,a$YTop[j]-a$YBottom[j],a$ModuleColor[j],if(a$WGCNAModule[j]=='Not_in_network')'grey60'else NA)
      ln(c(103.7,106.5,110.5),c(plan$HeatCenter[i],plan$HeatCenter[i],plan$HeaderY[i]-.3))
      txt(paste0(p,'  |  ',length(gs),' genes'),112,plan$HeaderY[i],7.1,'bold',col=pcols[[p]])
      a<-aa[Pattern==p]
      for(j in seq_len(nrow(a))){bw<-57*a$Score[j]/maxscore
        rect(112,a$BarY[j],bw,.8,pcols[[p]])
        txt(a$WrappedText[j],112,a$TextY[j],termfont,lineheight=lh/(termfont*25.4/72))
        bars[[length(bars)+1L]]<<-data.table(Pattern=p,Pathway=a$Pathway[j],Score=a$Score[j],WidthMM=bw,Color=pcols[[p]])}
    }
    # Three aligned legend columns with a shared header and row baselines.
    ln(c(8,169),c(60,60),'grey80')
    txt('Subtype',8,56,6.8,'bold');txt('Study',54,56,6.8,'bold');txt('WGCNA membership',104,56,6.8,'bold')
    for(j in seq_along(states)){st<-states[j];y<-50.5-(j-1)*4.3
      rect(8,y-.8,1.7,1.7,allpal[[st]]);txt(paste0(st,' (n = ',sum(columns$FinalSubtype==st),')'),11,y,6.1)}
    txt('Pattern accents',8,35.5,6.2,col='grey40')
    for(j in 1:2){p<-c('P2','P4')[j];x<-8+(j-1)*18;rect(x,30.7,1.7,1.7,pcols[[p]]);txt(p,x+3,31.5,6.1)}
    txt('Row z-score',8,24.5,6.5)
    for(j in 1:256)rect(8+(j-1)*29/256,20,29/256,1.6,zcols[j])
    for(v in c(-2,0,2))txt(as.character(v),8+(v+2)/4*29,17,5.9,just='center')
    studies<-names(stpal)[names(stpal) %in% columns$Study]
    for(j in seq_along(studies)){y<-50.5-(j-1)*4.3;rect(54,y-.8,1.7,1.7,stpal[[studies[j]]]);txt(studies[j],57,y,6.1)}
    mods<-names(modcols)[names(modcols) %in% rr$WGCNAModule]
    for(j in seq_along(mods)){x<-104+((j-1)%%3)*22;y<-50.5-floor((j-1)/3)*4.3
      rect(x,y-.75,1.55,1.55,modcols[[mods[j]]],if(mods[j]=='Not_in_network')'grey60'else NA)
      label<-if(mods[j]=='Not_in_network')'Outside network'else sub('^ME','',mods[j]);txt(label,x+2.2,y,5.8)}
    txt('ORA: all 10,531 tested genes | FDR < 0.05 | hits >= 5',54,17.2,6,col='grey35')
    txt('GO BP / Hallmark / Reactome / KEGG metabolic',54,13.3,6,col='grey35')

  }
  name<-paste0('Figure2A_fold_enrichment_',variant);pdf<-file.path(out,'figures',paste0(name,'.pdf'))
  save_panel(draw,pdf,width,height)
  ragg::agg_png(file.path(out,'figures',paste0(name,'.png')),width=width,height=height,units='mm',res=600);draw();dev.off()
  system2('/opt/homebrew/bin/pdftoppm',c('-scale-to','1900','-singlefile','-png',shQuote(pdf),shQuote(file.path(out,'review',paste0(name,'_preview')))))
  fwrite(aa,file.path(out,'tables/annotation_evidence.csv'));fwrite(plan,file.path(out,'tables/pattern_layout.csv'))
  fwrite(rr,file.path(out,'tables/row_order.csv'));fwrite(columns,file.path(out,'tables/sample_order.csv'))
  fwrite(rbindlist(list(
    slices[,.(Kind='Subtype',Label=State,X=Start,Width=End-Start,Y=height-22.5,Height=1.5,Color=unname(allpal[State]))],
    columns[,.(Kind='Study',Label=Study,X=XStart,Width=XEnd-XStart,Y=height-25.5,Height=1.5,Color=unname(stpal[Study]))])),
    file.path(out,'tables/top_annotation_layout.csv'))
  fwrite(rbindlist(texts),file.path(out,'tables/text_layout.csv'));fwrite(rbindlist(tiles),file.path(out,'tables/matrix_tiles.csv'))
  fwrite(rbindlist(bars),file.path(out,'tables/bar_geometry.csv'))
  fwrite(data.table(Source=paths,MD5=hashes),file.path(out,'review/input_provenance.csv'))
  fwrite(data.table(Figure='Figure 2',Panel='A',Role=paste0('Principal patterns with fold-enrichment-ranked ORA: ',variant),
    SourceTable=paste(paths[1:3],collapse=';'),TargetOutput=pdf,WidthMM=width,HeightMM=height,Status='Rendered exploratory candidate',
    QAFile=file.path(out,'review/verification.json'),ArtifactMD5=unname(tools::md5sum(pdf))),file.path(out,'figure_manifest.csv'))
  stopifnot(identical(hashes,unname(tools::md5sum(paths))))
  jsonlite::write_json(list(Genes=nrow(rr),Samples=nrow(columns),PatternGenes=as.list(setNames(plan$Genes,plan$Pattern)),
    Terms=nrow(aa),Variant=variant,BarMetric='FoldEnrichment',BarMaximum=maxscore,AnnotationStripHeightMM=1.5,WidthMM=width,HeightMM=height,HeatmapValuesUnchanged=TRUE,MembershipUnchanged=TRUE,SourceHashesUnchanged=TRUE),
    file.path(out,'review/render_diagnostics.json'),auto_unbox=TRUE,pretty=TRUE)
}
if(sys.nframe()==0){render_fig2A_fold_enrichment_v23('same29');render_fig2A_fold_enrichment_v23('FE_selected')}
