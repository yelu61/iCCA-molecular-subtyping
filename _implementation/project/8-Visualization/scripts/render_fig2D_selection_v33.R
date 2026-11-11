# Inferred TF activity: frozen effects, explicit display selection and sensitivity.
source('8-Visualization/scripts/fig_style.R')
suppressPackageStartupMessages({library(data.table);library(grid)})
out<-'8-Visualization/refinements/fig2D_selection_20261031_v33'
pal<-load_subtype_palette()$colors[c('C1','C2','C3')];invisible(load_study_palette())
mm<-function(x)unit(x,'mm');texts<-list();marks<-list()
txt<-function(s,x,y,size=5.8,face='plain',just='left',col='#25282B'){grid.text(s,mm(x),mm(y),just=just,gp=gpar(fontfamily='Arial',fontsize=size,fontface=face,col=col));texts[[length(texts)+1L]]<<-data.table(Text=s,X=x,Y=y,Size=size)}
rect<-function(x,y,w,h,fill)grid.rect(mm(x),mm(y),mm(w),mm(h),gp=gpar(fill=fill,col=NA))
ramp<-colorRamp(c(SCALE_EFFECT_LOW,'#FAFAFA',SCALE_EFFECT_HIGH));color<-function(v)rgb(ramp(pmax(0,pmin(1,(v/8+1)/2))),maxColorValue=255)
emit<-function(name,draw,w,h,source){pdf<-file.path(out,'figures',paste0(name,'.pdf'));save_panel(draw,pdf,w,h);ragg::agg_png(sub('pdf$','png',pdf),width=w,height=h,units='mm',res=600);draw();dev.off();fwrite(rbindlist(texts),file.path(out,'tables',paste0(name,'_texts.csv')));if(length(marks))fwrite(rbindlist(marks),file.path(out,'tables',paste0(name,'_marks.csv')));system2('/opt/homebrew/bin/pdftoppm',c('-scale-to','1700','-singlefile','-png',shQuote(pdf),shQuote(file.path(out,'review',paste0(name,'_preview')))));data.table(Panel=name,TargetOutput=pdf,WidthMM=w,HeightMM=h,SourceTable=file.path(out,source),SourceMD5=unname(tools::md5sum(file.path(out,source))),MD5=unname(tools::md5sum(pdf)),Status='Rendered review candidate')}
b<-fread(file.path(out,'tables/D_selected36.csv'));e<-fread(file.path(out,'tables/SFigure2D_extended54.csv'));stopifnot(nrow(b)==36,nrow(e)==54,all(abs(e$ActivityDifferenceSubtypeMinusOthers)<8))
make_d<-function(z,extended=FALSE){
 plan<-unique(z[,.(TF,AnchorSubtype,AnchorDirection,SelectionOrder,NonzeroSignedTargets)]);plan[,DirOrder:=match(AnchorDirection,c('Positive','Negative'))];setorder(plan,AnchorSubtype,DirOrder,SelectionOrder)
 H<-if(extended)108 else 73;xs<-if(extended)c(49,65,81)else c(27,39,51);cw<-if(extended)14 else 11;step<-if(extended)4 else 3.2
 function(){grid.newpage();texts<<-list();marks<<-list();txt('TF',2,H-4,6.3,'bold');txt('+/-',16,H-4,5.4,just='center',col='grey40');if(extended)txt('Targets',29,H-4,5.7,just='center')
 for(j in 1:3){txt(names(pal)[j],xs[j],H-4,6.3,'bold','center');rect(xs[j],H-6.8,cw,1.1,pal[j])};cursor<-H-12
 for(s in names(pal)){txt(paste0('Selected via ',s,' contrast'),2,cursor,5.5,'bold',col='grey40');cursor<-cursor-3
 for(tf in plan[AnchorSubtype==s,TF]){rr<-plan[TF==tf];txt(tf,2,cursor,6.1,'italic');txt(if(rr$AnchorDirection=='Positive')'+'else'-',16,cursor,5.8,just='center',col=pal[s]);if(extended)txt(rr$NonzeroSignedTargets,29,cursor,5.8,just='center')
 for(j in 1:3){r<-z[TF==tf&Subtype==names(pal)[j]];v<-r$ActivityDifferenceSubtypeMinusOthers;stopifnot(nrow(r)==1);rect(xs[j],cursor,cw,step-.25,color(v));txt(sprintf('%.1f',v),xs[j]-.5,cursor,5.6,just='center',col=if(abs(v)>5)'white'else'grey15');if(r$FDR<.05)grid.points(mm(xs[j]+cw/2-1.3),mm(cursor+.8),pch=16,size=mm(.48),gp=gpar(col=if(abs(v)>5)'white'else'grey15'));marks[[length(marks)+1L]]<<-data.table(TF=tf,Subtype=r$Subtype,Effect=v,FDR=r$FDR,X=xs[j],Y=cursor)}
 cursor<-cursor-step
 };cursor<-cursor-.8};stopifnot(cursor>10)
 txt('ULM difference',2,8,5.3);for(i in 1:100)rect(2+(i-.5)*.18,4.7,.19,1.4,color(-8+(i-1)/99*16));for(i in 1:3)txt(c('-8','0','8')[i],c(2,11,20)[i],2,5,just='center');grid.points(mm(33),mm(4.7),pch=16,size=mm(.5));txt('FDR < 0.05',35,4.7,5.2)
 }
}
m1<-emit('Figure2D_rule_selected_TF_activity',make_d(b),60,73,'tables/D_selected36.csv')
m2<-emit('SFigure2D_rule_selected_extended',make_d(e,TRUE),92,108,'tables/SFigure2D_extended54.csv')
s<-fread(file.path(out,'tables/D_sensitivity_membership.csv'));sc<-unique(s[,.(Mode,Threshold,Budget)]);sc[,Mo:=match(Mode,c('signed','unsigned'))];setorder(sc,Mo,Threshold,Budget);main<-unique(b[order(SelectionOrder),TF]);ids<-c(main,sort(setdiff(unique(s$TF),main)));H<-30+length(ids)*4.4
sdraw<-function(){grid.newpage();texts<<-list();marks<<-list();txt('TF selection sensitivity',2,H-4,8,'bold');xs<-seq(28,134,length.out=nrow(sc))
 for(k in 1:nrow(sc)){txt(if(sc$Mode[k]=='signed')'Signed'else'Unsigned',xs[k],H-11,5.4,just='center');txt(sprintf('%.2f',sc$Threshold[k]),xs[k],H-15,5.7,just='center');txt(sc$Budget[k],xs[k],H-19,5.7,just='center')};txt('Overlap',2,H-11,5.5);txt('Cutoff',2,H-15,5.5);txt('Count',2,H-19,5.5)
 for(i in seq_along(ids)){y<-H-26-(i-1)*4.4;txt(ids[i],2,y,6,face=if(ids[i]%in%main)'bold.italic'else'italic');for(k in 1:nrow(sc)){yes<-nrow(s[TF==ids[i]&Mode==sc$Mode[k]&Threshold==sc$Threshold[k]&Budget==sc$Budget[k]])>0;rect(xs[k],y,7,3.2,if(yes)'#506D82'else'#EEF1F3')}}
 txt('Bold: main 12 TFs. Filled: selected. Display-count nesting is not independent validation.',2,2.8,5.5,col='grey40')
}
m3<-emit('SFigure2D_selection_sensitivity',sdraw,141,H,'tables/D_sensitivity_membership.csv')
fwrite(rbindlist(list(m1,m2,m3)),file.path(out,'figure_manifest.csv'))
system2('/opt/homebrew/bin/pdfunite',shQuote(c(m1$TargetOutput,m2$TargetOutput,m3$TargetOutput,file.path(out,'Figure2D_selection_review.pdf'))))
