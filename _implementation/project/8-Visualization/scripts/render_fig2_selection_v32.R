# Frozen-statistic display refresh; selection provenance is stored separately.
source('8-Visualization/scripts/fig_style.R')
suppressPackageStartupMessages({library(data.table);library(grid)})
out<-'8-Visualization/refinements/fig2_selection_20261031_v32'
for(d in c('figures','tables','review'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
pal<-load_subtype_palette()$colors[c('C1','C2','C3')];invisible(load_study_palette())
mm<-function(x)unit(x,'mm');texts<-list();marks<-list()
txt<-function(s,x,y,size=5.8,face='plain',just='left',col='#25282B'){grid.text(s,mm(x),mm(y),just=just,gp=gpar(fontfamily='Arial',fontsize=size,fontface=face,col=col));texts[[length(texts)+1L]]<<-data.table(Text=s,X=x,Y=y,Size=size)}
rect<-function(x,y,w,h,fill,border=NA)grid.rect(mm(x),mm(y),mm(w),mm(h),gp=gpar(fill=fill,col=border,lwd=.5))
line<-function(x,y,col='#CDD3D8',lwd=.5,lty=1)grid.lines(mm(x),mm(y),gp=gpar(col=col,lwd=lwd,lty=lty))
point<-function(x,y,fill,col,size=1.45)grid.points(mm(x),mm(y),pch=21,size=mm(size),gp=gpar(fill=fill,col=col,lwd=.65))
ramp<-colorRamp(c(SCALE_EFFECT_LOW,'#FAFAFA',SCALE_EFFECT_HIGH));color<-function(v,lim=4.5)rgb(ramp(pmax(0,pmin(1,(v/lim+1)/2))),maxColorValue=255)
label<-function(s){s<-gsub('^(GOBP|HALLMARK|REACTOME|KEGG_MEDICUS)_','',s);s<-tolower(gsub('_',' ',s));s<-paste0(toupper(substr(s,1,1)),substring(s,2));for(w in c('wnt','nmd','g2m','rna','dna','nfkb','tnfa'))s<-gsub(paste0('\\b',w,'\\b'),toupper(w),s,ignore.case=TRUE);s}
emit<-function(name,draw,w,h,source){pdf<-file.path(out,'figures',paste0(name,'.pdf'));save_panel(draw,pdf,w,h);ragg::agg_png(sub('pdf$','png',pdf),width=w,height=h,units='mm',res=600);draw();dev.off();fwrite(rbindlist(texts),file.path(out,'tables',paste0(name,'_text_geometry.csv')));if(length(marks))fwrite(rbindlist(marks),file.path(out,'tables',paste0(name,'_mark_geometry.csv')));system2('/opt/homebrew/bin/pdftoppm',c('-scale-to','1800','-singlefile','-png',shQuote(pdf),shQuote(file.path(out,'review',paste0(name,'_preview')))));data.table(Panel=name,TargetOutput=pdf,WidthMM=w,HeightMM=h,SourceTable=file.path(out,source),SourceMD5=unname(tools::md5sum(file.path(out,source))),Status='Rendered review candidate',MD5=unname(tools::md5sum(pdf)))}
b<-fread(file.path(out,'tables/B_selected54.csv'));stopifnot(nrow(b)==54)
make_b<-function(z,w=60,h=90){
 plan<-unique(z[,.(pathway,AnchorSubtype,AnchorDirection,SelectionOrder,Collection)]);plan[,DirOrder:=match(AnchorDirection,c('Positive','Negative'))];setorder(plan,AnchorSubtype,DirOrder,SelectionOrder);plan[,Label:=label(pathway)]
 # Measure exact labels at the chosen point size; no pathway-specific abbreviation.
 metric<-tempfile(fileext='.pdf');grDevices::cairo_pdf(metric,family='Arial');pushViewport(viewport(gp=gpar(fontfamily='Arial',fontsize=5.2)))
 plan[,LabelWidthMM:=vapply(Label,function(s)convertWidth(stringWidth(s),'mm',valueOnly=TRUE),numeric(1))];popViewport();dev.off();unlink(metric)
 stopifnot(max(plan$LabelWidthMM)<w-24)
 step<-if(nrow(plan)==18)3.05 else 3.2
 function(){grid.newpage();texts<<-list();marks<<-list();xs<-c(w-14,w-8.5,w-3)
 for(j in 1:3){txt(names(pal)[j],xs[j],h-4,6.3,'bold','center');rect(xs[j],h-7,5.1,1.1,pal[j])}
 txt('DB',w-19,h-4,4.9,just='center',col='grey45');txt('+/-',w-23,h-4,4.9,just='center',col='grey45');cursor<-h-12
 for(s in names(pal)){txt(paste0('Selected via ',s,' contrast'),1,cursor,5.4,'bold',col='grey35');cursor<-cursor-3
  for(id in plan[AnchorSubtype==s,pathway]){r<-plan[pathway==id];txt(r$Label,1,cursor,5.2);txt(if(r$AnchorDirection=='Positive')'+'else'-',w-23,cursor,5.2,just='center',col=pal[s]);txt(c(GO_BP='G',Reactome='R',Hallmark='H',KEGG='K')[[r$Collection]],w-19,cursor,5.1,just='center',col='grey45')
  for(j in 1:3){a<-z[pathway==id&Subtype==names(pal)[j]];stopifnot(nrow(a)==1);rect(xs[j],cursor,5.1,step-.3,color(a$NES));if(a$padj<.05)grid.points(mm(xs[j]),mm(cursor),pch=16,size=mm(.55),gp=gpar(col=if(abs(a$NES)>3)'white'else'grey20'));marks[[length(marks)+1L]]<<-data.table(Pathway=id,Subtype=a$Subtype,NES=a$NES,FDR=a$padj,X=xs[j],Y=cursor)}
  cursor<-cursor-step
 };cursor<-cursor-1.2}
 stopifnot(cursor>10)
 txt('NES',1,8,5.5);for(i in 1:100)rect(8+(i-.5)*.15,8,.16,1.4,color(-4.5+(i-1)/99*9));for(i in 1:3)txt(c('-4.5','0','4.5')[i],c(8,15.5,23)[i],5.2,5,just='center');grid.points(mm(30),mm(8),pch=16,size=mm(.55));txt('FDR < 0.05',32,8,5.2)
 txt('G GO BP   R Reactome   H Hallmark',1,1.8,5,col='grey40')
 }
}
mb<-emit('Figure2B_rule_selected_GSEA',make_b(b),60,90,'tables/B_selected54.csv')
full<-fread('1-Clustering/publication/tables/step1_ranked_gsea_main_only_study_adjusted.csv');s<-fread(file.path(out,'tables/B_sensitivity_membership.csv'));p24<-s[Threshold==.5&Budget==24,.(pathway=Pathway,AnchorSubtype,AnchorDirection,SelectionOrder)];b24<-merge(p24,full,by='pathway');fwrite(b24,file.path(out,'tables/SFigure2B_extended72.csv'))
ms<-emit('SFigure2B_rule_selected_extended',make_b(b24,100,114),100,114,'tables/SFigure2B_extended72.csv')
cc<-fread(file.path(out,'tables/C_all33_effects.csv'));cp<-fread(file.path(out,'tables/C_all11_modules.csv'));stopifnot(nrow(cc)==33,nrow(cp)==11,all(cc$CI_low> -2.5),all(cc$CI_high<2.5),max(cp$ZsummaryPreservation)<25)
cdraw<-function(){grid.newpage();texts<<-list();marks<<-list();xs<-c(61,94,127);ys<-c(87,80,69,62,55,48,37,30,23,16,9);scale<-28/5
 txt('Co-expression program',2,102,7,'bold');for(j in 1:3){txt(names(pal)[j],xs[j],102,7,'bold','center');rect(xs[j],99,29,1.1,pal[j]);line(rep(xs[j],2),c(6,94),'#CCD2D8',.5,2)};txt('Preservation',160,102,6.3,'bold','center');txt('Z summary',160,98,6.1,just='center',col='grey40');for(v in c(2,10))line(rep(150+v/25*20,2),c(6,94),'#D4D9DD',.5,2)
 for(i in seq_len(nrow(cp))){y<-ys[i];r<-cp[i];rect(2,y,1.25,4.7,r$ModuleColor);txt(r$Program,4,y+1.1,6.5,'bold');txt(paste0(r$ModuleColor,'  |  ',r$ModuleSize,' genes'),4,y-1.8,5.5,col='grey45')
  for(j in 1:3){z<-cc[Module==r$Module&TargetSubtype==names(pal)[j]];line(xs[j]+c(z$CI_low,z$CI_high)*scale,c(y,y),pal[j],.85);point(xs[j]+z$RandomEffect*scale,y,if(z$FDR<.05)pal[j]else'white',pal[j],1.7);marks[[length(marks)+1L]]<<-data.table(Module=r$Module,Subtype=z$TargetSubtype,Effect=z$RandomEffect,CI_low=z$CI_low,CI_high=z$CI_high,FDR=z$FDR,X=xs[j]+z$RandomEffect*scale,Y=y)}
 line(c(150,150+r$ZsummaryPreservation/25*20),c(y,y),'#AEBAC2',.65);point(150+r$ZsummaryPreservation/25*20,y,'#455B69','#455B69',1.45)
 }
 for(j in 1:3){line(xs[j]+c(-2.5,2.5)*scale,c(4.5,4.5),'grey50');for(v in c(-2,0,2))txt(v,xs[j]+v*scale,2.5,5.5,just='center')};for(v in c(0,10,20))txt(v,150+v/25*20,2.5,5.5,just='center')
 # Keep explanatory legend separate from numeric axes, within the upper margin.
 txt('Difference (SD), 95% CI',3,94,5.7,col='grey40');point(82,94,'#455B69','#455B69',1.1);txt('FDR < 0.05',84,94,5.4);point(111,94,'white','#455B69',1.1);txt('FDR >= 0.05',113,94,5.4)
}
mc<-emit('Figure2C_all11_network_modules',cdraw,174,106,'tables/C_all33_effects.csv')
# Selection sensitivity is a diagnostic display, not an independent enrichment result.
sc<-unique(s[,.(Threshold,Budget)]);setorder(sc,Threshold,Budget);freq<-fread(file.path(out,'tables/B_selection_frequency.csv'));mainids<-unique(b[order(SelectionOrder),pathway]);rest<-freq[!Pathway%in%mainids][order(-SelectedInScenarios,Pathway),Pathway];ids<-c(mainids,rest);H<-20+length(ids)*4.2+13
sdraw<-function(){grid.newpage();texts<<-list();marks<<-list();txt('Gene-set selection sensitivity',2,H-4,8,'bold');xs<-seq(100,164,length.out=9)
 for(k in 1:9){txt(sprintf('%.2f',sc$Threshold[k]),xs[k],H-10,5.9,just='center');txt(sc$Budget[k],xs[k],H-14,5.9,just='center')};txt('Jaccard cutoff',73,H-10,5.8,just='right');txt('Display count',73,H-14,5.8,just='right')
 for(i in seq_along(ids)){y<-H-21-(i-1)*4.2;txt(label(ids[i]),3,y,6,face=if(ids[i]%in%mainids)'bold'else'plain');for(k in 1:9){yes<-nrow(s[Pathway==ids[i]&Threshold==sc$Threshold[k]&Budget==sc$Budget[k]])>0;rect(xs[k],y,5.5,3.2,if(yes)'#506D82'else'#EEF1F3')}}
 txt('Bold: main 18 sets; filled cells: selected. Nested display budgets are not independent validation.',3,5,5.8,col='grey40')
}
mt<-emit('SFigure2B_selection_sensitivity',sdraw,174,H,'tables/B_sensitivity_membership.csv')
fwrite(rbindlist(list(mb,mc,ms,mt)),file.path(out,'figure_manifest.csv'))
# Review pages preserve vector text; do not resize the prior A/D/E panels.
system2('/opt/homebrew/bin/pdfunite',shQuote(c(mb$TargetOutput,mc$TargetOutput,ms$TargetOutput,mt$TargetOutput,file.path(out,'Figure2_BC_selection_review.pdf'))))
