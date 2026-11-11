# Source-preserving R layout assembly; all five current main panels, separate A legend.
source('8-Visualization/scripts/render_fig2A_split_v34.R')
suppressPackageStartupMessages({library(data.table);library(grid)})
out<-'8-Visualization/refinements/fig2_assembly_20261031_v34'
mm<-function(x)unit(x,'mm');texts<-list();marks<-list()
pal<-load_subtype_palette()$colors[c('C1','C2','C3')];invisible(load_study_palette())
txt<-function(s,x,y,size=5.8,face='plain',just='left',col='#25282B'){grid.text(s,mm(x),mm(y),just=just,gp=gpar(fontfamily='Arial',fontsize=size,fontface=face,col=col));texts[[length(texts)+1L]]<<-data.table(Text=s,X=x,Y=y,Size=size)}
rect<-function(x,y,w,h,fill,border=NA)grid.rect(mm(x),mm(y),mm(w),mm(h),gp=gpar(fill=fill,col=border,lwd=.5))
line<-function(x,y,col='#CDD3D8',lwd=.5,lty=1)grid.lines(mm(x),mm(y),gp=gpar(col=col,lwd=lwd,lty=lty))
point<-function(x,y,fill,col,size=1.2)grid.points(mm(x),mm(y),pch=21,size=mm(size),gp=gpar(fill=fill,col=col,lwd=.65))
ramp<-colorRamp(c(SCALE_EFFECT_LOW,'#FAFAFA',SCALE_EFFECT_HIGH));basecolor<-function(v,lim)rgb(ramp(pmax(0,pmin(1,(v/lim+1)/2))),maxColorValue=255)
label<-function(s){s<-gsub('^(GOBP|HALLMARK|REACTOME|KEGG_MEDICUS)_','',s);s<-tolower(gsub('_',' ',s));s<-paste0(toupper(substr(s,1,1)),substring(s,2));for(w in c('wnt','nmd','g2m','rna','dna','nfkb','tnfa'))s<-gsub(paste0('\\b',w,'\\b'),toupper(w),s,ignore.case=TRUE);s}
emit<-function(name,draw,w,h,source){pdf<-file.path(out,'figures',paste0(name,'.pdf'));save_panel(draw,pdf,w,h);ragg::agg_png(sub('pdf$','png',pdf),width=w,height=h,units='mm',res=600);draw();dev.off();fwrite(rbindlist(texts),file.path(out,'tables',paste0(name,'_texts.csv')));if(length(marks))fwrite(rbindlist(marks),file.path(out,'tables',paste0(name,'_marks.csv')));system2('/opt/homebrew/bin/pdftoppm',c('-scale-to','2400','-singlefile','-png',shQuote(pdf),shQuote(file.path(out,'review',paste0(name,'_preview')))));data.table(Panel=name,TargetOutput=pdf,WidthMM=w,HeightMM=h,SourceTable=paste(vapply(strsplit(source,';',fixed=TRUE)[[1]],function(p)if(file.exists(p))p else file.path(out,p),character(1)),collapse=';'),Status='Rendered layout candidate',MD5=unname(tools::md5sum(pdf)))}
make_b<-function(z,w=60,h=82){
 color<-function(v)basecolor(v,4.5)
 plan<-unique(z[,.(pathway,AnchorSubtype,AnchorDirection,SelectionOrder,Collection)]);plan[,DirOrder:=match(AnchorDirection,c('Positive','Negative'))];setorder(plan,AnchorSubtype,DirOrder,SelectionOrder);plan[,Label:=label(pathway)]
 # Measure exact labels at the chosen point size; no pathway-specific abbreviation.
 metric<-tempfile(fileext='.pdf');grDevices::cairo_pdf(metric,family='Arial');pushViewport(viewport(gp=gpar(fontfamily='Arial',fontsize=5.2)))
 plan[,LabelWidthMM:=vapply(Label,function(s)convertWidth(stringWidth(s),'mm',valueOnly=TRUE),numeric(1))];popViewport();dev.off();unlink(metric)
 stopifnot(max(plan$LabelWidthMM)<w-24)
 step<-if(nrow(plan)==18)2.65 else 3.2
 function(newpage=TRUE){if(newpage)grid.newpage();texts<<-list();marks<<-list();xs<-c(w-14,w-8.5,w-3)
 for(j in 1:3){txt(names(pal)[j],xs[j],h-4,6.3,'bold','center');rect(xs[j],h-7,5.1,1.1,pal[j])}
 txt('DB',w-19,h-4,5,just='center',col='grey45');txt('+/-',w-23,h-4,5,just='center',col='grey45');cursor<-h-12
 for(s in names(pal)){txt(paste0('Selected via ',s,' contrast'),1,cursor,5.4,'bold',col='grey35');cursor<-cursor-3
  for(id in plan[AnchorSubtype==s,pathway]){r<-plan[pathway==id];txt(r$Label,1,cursor,5.2);txt(if(r$AnchorDirection=='Positive')'+'else'-',w-23,cursor,5.2,just='center',col=pal[s]);txt(c(GO_BP='G',Reactome='R',Hallmark='H',KEGG='K')[[r$Collection]],w-19,cursor,5.1,just='center',col='grey45')
  for(j in 1:3){a<-z[pathway==id&Subtype==names(pal)[j]];stopifnot(nrow(a)==1);rect(xs[j],cursor,5.1,step-.3,color(a$NES));if(a$padj<.05)grid.points(mm(xs[j]),mm(cursor),pch=16,size=mm(.55),gp=gpar(col=if(abs(a$NES)>3)'white'else'grey20'));marks[[length(marks)+1L]]<<-data.table(Pathway=id,Subtype=a$Subtype,NES=a$NES,FDR=a$padj,X=xs[j],Y=cursor)}
  cursor<-cursor-step
 };cursor<-cursor-.8}
 stopifnot(cursor>10)
 txt('NES',1,8,5.5);for(i in 1:100)rect(8+(i-.5)*.15,8,.16,1.4,color(-4.5+(i-1)/99*9));for(i in 1:3)txt(c('-4.5','0','4.5')[i],c(8,15.5,23)[i],5.2,5,just='center');grid.points(mm(30),mm(8),pch=16,size=mm(.55));txt('FDR < 0.05',32,8,5.2)
 txt('G GO BP   R Reactome   H Hallmark',1,1.8,5,col='grey40')
 }
}
make_d<-function(z,extended=FALSE){
 color<-function(v)basecolor(v,8)
 plan<-unique(z[,.(TF,AnchorSubtype,AnchorDirection,SelectionOrder,NonzeroSignedTargets)]);plan[,DirOrder:=match(AnchorDirection,c('Positive','Negative'))];setorder(plan,AnchorSubtype,DirOrder,SelectionOrder)
 H<-if(extended)108 else 76;xs<-if(extended)c(49,65,81)else c(27,39,51);cw<-if(extended)14 else 11;step<-if(extended)4 else 3.35
 function(newpage=TRUE){if(newpage)grid.newpage();texts<<-list();marks<<-list();txt('TF',5,H-4,6.3,'bold');txt('+/-',16,H-4,5.4,just='center',col='grey40');if(extended)txt('Targets',29,H-4,5.7,just='center')
 for(j in 1:3){txt(names(pal)[j],xs[j],H-4,6.3,'bold','center');rect(xs[j],H-6.8,cw,1.1,pal[j])};cursor<-H-12
 for(s in names(pal)){txt(paste0('Selected via ',s,' contrast'),2,cursor,5.5,'bold',col='grey40');cursor<-cursor-3
 for(tf in plan[AnchorSubtype==s,TF]){rr<-plan[TF==tf];txt(tf,2,cursor,6.1,'italic');txt(if(rr$AnchorDirection=='Positive')'+'else'-',16,cursor,5.8,just='center',col=pal[s]);if(extended)txt(rr$NonzeroSignedTargets,29,cursor,5.8,just='center')
 for(j in 1:3){r<-z[TF==tf&Subtype==names(pal)[j]];v<-r$ActivityDifferenceSubtypeMinusOthers;stopifnot(nrow(r)==1);rect(xs[j],cursor,cw,step-.25,color(v));txt(sprintf('%.1f',v),xs[j]-.5,cursor,5.6,just='center',col=if(abs(v)>5)'white'else'grey15');if(r$FDR<.05)grid.points(mm(xs[j]+cw/2-1.3),mm(cursor+.8),pch=16,size=mm(.48),gp=gpar(col=if(abs(v)>5)'white'else'grey15'));marks[[length(marks)+1L]]<<-data.table(TF=tf,Subtype=r$Subtype,Effect=v,FDR=r$FDR,X=xs[j],Y=cursor)}
 cursor<-cursor-step
 };cursor<-cursor-.8};stopifnot(cursor>10)
 txt('ULM difference',2,8,5.3);for(i in 1:100)rect(2+(i-.5)*.18,4.7,.19,1.4,color(-8+(i-1)/99*16));for(i in 1:3)txt(c('-8','0','8')[i],c(2,11,20)[i],2,5,just='center');grid.points(mm(33),mm(4.7),pch=16,size=mm(.5));txt('FDR < 0.05',35,4.7,5.2)
 }
}
aa<-render_fig2A_split_v34()
legendkeys<-list()
ldraw<-function(newpage=TRUE){if(newpage)grid.newpage();texts<<-list();marks<<-list();legendkeys<<-list()
 key<-function(group,label,color,x,y,border=NA){rect(x+.55,y,1.1,1.1,color,border);txt(label,x+1.7,y,5);legendkeys[[length(legendkeys)+1L]]<<-data.table(Group=group,Label=label,Color=color,X=x,Y=y)}
 txt('Subtype / pattern',2,36,5,'bold');labs<-c('C1','C2','C3','P2','P4');cols<-c(aa$allpal[c('C1','C2','C3')],aa$pcols[c('P2','P4')])
 for(i in 1:5)key('SubtypePattern',labs[i],cols[i],18+(i-1)*8,36)
 txt('Study',2,32.5,5.2,'bold');for(i in seq_along(aa$studies))key('Study',aa$studies[i],aa$stpal[aa$studies[i]],2+((i-1)%%3)*20,29.8-((i-1)%/%3)*2.5)
 txt('WGCNA',2,21.8,5.2,'bold');for(i in seq_along(aa$modules)){m<-aa$modules[i];key('WGCNA',if(m=='Not_in_network')'Outside network'else sub('^ME','',m),aa$modcols[m],2+((i-1)%%4)*14,19-((i-1)%/%4)*2.5,if(m=='Not_in_network')'grey60'else NA)}
 txt('Z-score',2,3.2,5);for(i in 1:256)rect(17+(i-.5)*25/256,3.2,25/256+.01,1.2,aa$zcols[i]);for(v in c(-2,0,2))txt(v,17+(v+2)/4*25,.9,5,just='center')
}
b<-fread('8-Visualization/refinements/fig2_selection_20261031_v32/tables/B_selected54.csv');dd<-fread('8-Visualization/refinements/fig2D_selection_20261031_v33/tables/D_selected36.csv');bdraw<-make_b(b);ddraw<-make_d(dd)
cc<-fread('8-Visualization/refinements/fig2_selection_20261031_v32/tables/C_all33_effects.csv');cp<-fread('8-Visualization/refinements/fig2_selection_20261031_v32/tables/C_all11_modules.csv')
cdraw<-function(newpage=TRUE){if(newpage)grid.newpage();texts<<-list();marks<<-list();xs<-c(77,108,139);ys<-40.5-(0:10)*3-c(0,0,rep(1.5,4),rep(3,5));scale<-28/5
 txt('Co-expression program',6,49,6.2,'bold');txt('Module | genes',43,49,5.2,'bold');for(j in 1:3){txt(names(pal)[j],xs[j],49,6.3,'bold','center');rect(xs[j],46.3,27,1.1,pal[j]);line(rep(xs[j],2),c(6,43),'#CCD2D8',.5,2)};txt('Preservation',169,49,6,'bold','center');txt('Z summary',169,46,5.2,just='center',col='grey40')
 for(v in c(2,10))line(rep(158+v/25*20,2),c(6,43),'#D4D9DD',.5,2)
 txt('Difference (SD), 95% CI',6,44,5.1,col='grey40');point(86,43.4,'#455B69','#455B69',.9);txt('FDR < 0.05',88,43.4,5);point(114,43.4,'white','#455B69',.9);txt('FDR >= 0.05',116,43.4,5)
 for(i in seq_len(nrow(cp))){y<-ys[i];r<-cp[i];rect(2,y,1.1,2.4,r$ModuleColor);txt(r$Program,4,y,5.8,'bold');txt(paste0(r$ModuleColor,' | ',r$ModuleSize),43,y,5,col='grey45')
 for(j in 1:3){z<-cc[Module==r$Module&TargetSubtype==names(pal)[j]];line(xs[j]+c(z$CI_low,z$CI_high)*scale,c(y,y),pal[j],.75);point(xs[j]+z$RandomEffect*scale,y,if(z$FDR<.05)pal[j]else'white',pal[j],1.25);marks[[length(marks)+1L]]<<-data.table(Module=r$Module,Subtype=z$TargetSubtype,Effect=z$RandomEffect,CI_low=z$CI_low,CI_high=z$CI_high,FDR=z$FDR,X=xs[j]+z$RandomEffect*scale,Y=y)}
 line(c(158,158+r$ZsummaryPreservation/25*20),c(y,y),'#AEBAC2',.65);point(158+r$ZsummaryPreservation/25*20,y,'#455B69','#455B69',1.15)
 }
 for(j in 1:3){line(xs[j]+c(-2.5,2.5)*scale,c(4.2,4.2),'grey50');for(v in c(-2,0,2))txt(v,xs[j]+v*scale,2,5,just='center')};for(v in c(0,10,20))txt(v,158+v/25*20,2,5,just='center')
}
em<-fread('8-Visualization/refinements/fig2_layout_audit_20261031_v31/E/tables/E_selection_map.csv');ez<-fread('8-Visualization/refinements/fig2_layout_audit_20261031_v31/E/tables/E_selected_contrasts.csv')
em[,Y:=63-(Y-14)/.72*.39];ez[,Y:=63-(Y-14)/.72*.39]
fwrite(em,file.path(out,'tables/E_selection_map.csv'));fwrite(ez,file.path(out,'tables/E_selected_contrasts.csv'))
eshape<-c('MCP-counter'=21,'EPIC'=22,'xCell'=24,'CIBERSORT'=23,'ESTIMATE'=25);eoffset<-c('MCP-counter'=.8,'EPIC'=.2667,'xCell'=-.2667,'CIBERSORT'=-.8,'ESTIMATE'=0)
edraw<-function(newpage=TRUE){if(newpage)grid.newpage();texts<<-list();marks<<-list();centers<-c(65,86,107)
 dot<-function(x,y,method,fill,col,size=1.05)grid.points(mm(x),mm(y),pch=eshape[method],size=mm(size),gp=gpar(col=col,fill=fill,lwd=.6))
 txt('Compartment',5,71.5,5.7,'bold');txt('Readout',25,71.5,5.7,'bold');for(j in 1:3){txt(names(pal)[j],centers[j],71.5,6.3,'bold','center');rect(centers[j],68.8,18,1.1,pal[j])}
 groups<-unique(em$Group)
 for(i in seq_along(groups)){g<-groups[i];yy<-em[Group==g,Y];lo<-min(yy)-1.5;hi<-max(yy)+1.5;if(i%%2==1)rect(60,(lo+hi)/2,118,hi-lo,'#F5F7F9');rect(1.8,(lo+hi)/2,.6,hi-lo,'#B8C3CC');txt(if(g=='Overall context')'Overall\ncontext'else g,3,(lo+hi)/2,5.3,'bold')}
 for(j in 1:3)line(rep(centers[j],2),c(24,65),'#CCD2D8',.5,2)
 labs<-unique(em[,.(Label,Y)]);for(i in seq_len(nrow(labs)))txt(labs$Label[i],25,labs$Y[i],5.5)
 for(i in seq_len(nrow(ez))){r<-ez[i];j<-match(r$Subtype,names(pal));y<-r$Y+eoffset[r$Method];xx<-function(v)centers[j]+4.4*v;stopifnot(isTRUE(r$Estimable),r$Lower>=-2,r$Upper<=2);line(xx(c(r$Lower,r$Upper)),c(y,y),adjustcolor(pal[j],alpha.f=.5),.6);dot(xx(r$Effect),y,r$Method,if(r$FDR<.05)pal[j]else'white',pal[j]);marks[[length(marks)+1L]]<<-data.table(Method=r$Method,Score=r$Score,Subtype=r$Subtype,Effect=r$Effect,Lower=r$Lower,Upper=r$Upper,FDR=r$FDR,N=r$N,X=xx(r$Effect),Y=y)}
 for(j in 1:3){line(centers[j]+c(-8.8,8.8),c(23.7,23.7),'grey50');for(v in c(-2,0,2))txt(v,centers[j]+4.4*v,21.3,5,just='center')};txt('Study-adjusted difference (within-study SD)',86,17.6,5.3,just='center')
 methods<-names(eshape);lx<-c(3,30,49,70,98);for(i in seq_along(methods)){dot(lx[i],11.5,methods[i],'#444C53','#444C53',1.2);txt(methods[i],lx[i]+2.2,11.5,5.3)}
 dot(3,6,'MCP-counter','#444C53','#444C53');txt('FDR < 0.05',5.2,6,5.2);dot(32,6,'MCP-counter','white','#444C53');txt('FDR >= 0.05',34.2,6,5.2);line(c(65,71),c(6,6),'grey50');txt('95% CI',74,6,5.2)
}
ml<-emit('Figure2A_legend',ldraw,60,39,'A/tables/row_order.csv;A/tables/sample_order.csv');fwrite(rbindlist(legendkeys),file.path(out,'tables/A_legend_keys.csv'))
mb<-emit('Figure2B_layout',bdraw,60,82,'8-Visualization/refinements/fig2_selection_20261031_v32/tables/B_selected54.csv')
mc<-emit('Figure2C_layout',cdraw,183,53,'8-Visualization/refinements/fig2_selection_20261031_v32/tables/C_all33_effects.csv')
md<-emit('Figure2D_layout',ddraw,60,76,'8-Visualization/refinements/fig2D_selection_20261031_v33/tables/D_selected36.csv')
me<-emit('Figure2E_layout',edraw,120,76,'tables/E_selected_contrasts.csv')
slots<-data.table(Panel=c('A','B','A_legend','C','D','E'),X=c(0,123,123,0,0,63),Y=c(137,178,137,80,0,0),Width=c(120,60,60,183,60,120),Height=c(123,82,39,53,76,76))
fwrite(slots,file.path(out,'tables/layout_slots.csv'))
assembled<-function(){grid.newpage();draws<-list(A=aa$draw,B=bdraw,A_legend=ldraw,C=cdraw,D=ddraw,E=edraw)
 for(i in seq_len(nrow(slots))){s<-slots[i];pushViewport(viewport(x=mm(s$X),y=mm(s$Y),width=mm(s$Width),height=mm(s$Height),just=c('left','bottom')));draws[[s$Panel]](FALSE);popViewport()}
 for(i in which(slots$Panel!='A_legend')){s<-slots[i];txt(s$Panel,s$X+.2,s$Y+s$Height-2.3,8,'bold')}
}
mf<-emit('Figure2_full_layout',assembled,183,260,'tables/layout_slots.csv')
ma<-fread(file.path(out,'A/figure_manifest.csv'));fwrite(rbindlist(list(ma,ml,mb,mc,md,me,mf),fill=TRUE),file.path(out,'figure_manifest.csv'))
