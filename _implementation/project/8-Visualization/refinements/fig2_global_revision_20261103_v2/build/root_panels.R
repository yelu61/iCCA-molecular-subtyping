suppressPackageStartupMessages({library(data.table);library(grid)})
source('8-Visualization/scripts/fig_style.R')
out<-'8-Visualization/refinements/fig2_global_revision_20261103_v2';for(d in c('figures','source_data','review'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
pal<-load_subtype_palette()$colors[c('C1','C2','C3')];invisible(load_study_palette());mm<-function(x)unit(x,'mm')
txt<-function(t,x,y,sz=6,face='plain',just='left',col='#25282B')grid.text(t,mm(x),mm(y),just=just,gp=gpar(fontfamily='Arial',fontsize=sz,fontface=face,col=col))
rect<-function(x,y,w,h,col)grid.rect(mm(x),mm(y),mm(w),mm(h),gp=gpar(fill=col,col=NA))
line<-function(x,y,col='#CCD2D8',lwd=.5,lty=1)grid.lines(mm(x),mm(y),gp=gpar(col=col,lwd=lwd,lty=lty))
emit<-function(n,draw,w,h){f<-file.path(out,'figures',n);cairo_pdf(paste0(f,'.pdf'),width=w/25.4,height=h/25.4,family='Arial');draw();dev.off();svglite::svglite(paste0(f,'.svg'),width=w/25.4,height=h/25.4);draw();dev.off();ragg::agg_png(paste0(f,'.png'),width=w,height=h,units='mm',res=600);draw();dev.off();ragg::agg_png(file.path(out,'review',paste0(n,'.png')),width=w,height=h,units='mm',res=180);draw();dev.off()}
# Exact frozen GSEA results; no selection or FDR is rerun.
groot<-'8-Visualization/refinements/fig2_selection_20261031_v32';full<-fread('1-Clustering/publication/tables/step1_ranked_gsea_main_only_study_adjusted.csv');sel<-fread(file.path(groot,'tables/B_selection18.csv'));ext<-fread(file.path(groot,'tables/SFigure2B_extended72.csv'));sens<-fread(file.path(groot,'tables/B_sensitivity_membership.csv'))
setnames(sel,'Pathway','pathway');stopifnot(nrow(ext)==72,uniqueN(sel$pathway)==18)
summary<-full[,.(Tested=.N,Up=sum(padj<.05&NES>0),Down=sum(padj<.05&NES<0)),by=.(Collection,Subtype)]
fwrite(summary,file.path(out,'source_data/GSEA_all_tested_summary.csv'));file.copy('1-Clustering/publication/tables/step1_ranked_gsea_main_only_study_adjusted.csv',file.path(out,'source_data/GSEA_all_tested_original.csv'),overwrite=TRUE)
for(n in c('B_selection18.csv','B_selected54.csv','B_selection_decisions.csv','B_all_set_disposition.csv','B_sensitivity_summary.csv','B_sensitivity_membership.csv','SFigure2B_extended72.csv'))file.copy(file.path(groot,'tables',n),file.path(out,'source_data',n),overwrite=TRUE)
collections<-c('GO_BP','Reactome','Hallmark','KEGG');collection_labels<-c('GO BP','Reactome','Hallmark','KEGG MEDICUS')
draw_bg<-function(){grid.newpage();txt('Complete GSEA testing and display selection',2,52,7,'bold');txt('Collection',2,43,6,'bold');txt('Tested / contrast',47,43,6,'bold','center');for(j in 1:3){txt(names(pal)[j],86+(j-1)*28,46,6.5,'bold','center');txt('NES+ / NES-',86+(j-1)*28,42,5.7,just='center')};txt('Main sets',164,43,6,'bold','center')
 for(i in 1:4){y<-35-(i-1)*6;if(i%%2)rect(87,y,172,5.5,'#F3F5F7');cc<-collections[i];z<-summary[Collection==cc];stopifnot(length(unique(z$Tested))==1);txt(collection_labels[i],2,y,6);txt(unique(z$Tested),47,y,6,just='center');for(j in 1:3){q<-z[Subtype==names(pal)[j]];txt(paste(q$Up,q$Down,sep=' / '),86+(j-1)*28,y,6,just='center')};txt(sum(sel$Collection==cc),164,y,6,just='center')}
 txt('Counts: original FDR < 0.05 within subtype and collection; no pooling of testing families',2,10,5.7)
 txt('Six subtype x direction queues; FDR then |NES| ranking; measured-gene Jaccard < 0.50; 3 rounds -> 18 sets',2,5.3,5.6)
 txt('All tested NES, FDR and leading-edge genes accompany the figures; extensions are display sensitivity',2,1.8,5.4,col='#5A656F')}
emit('GSEA_background',draw_bg,174,56)
# Reuse only the frozen extended-heatmap drawer, redirecting all output to this revision.
code<-readLines('8-Visualization/scripts/render_fig2_selection_v32.R');env<-new.env();eval(parse(text=code[c(2:3,6:13,16:35)]),envir=env)
# Function boundaries checked below if historical source changes.
emit('GSEA_extended24',env$make_b(ext,100,114),100,114)
sc<-unique(sens[,.(Threshold,Budget)]);setorder(sc,Threshold,Budget);mainids<-sel[order(SelectionOrder),pathway];freq<-fread(file.path(groot,'tables/B_selection_frequency.csv'));ids<-c(mainids,freq[!Pathway%in%mainids][order(-SelectedInScenarios,Pathway),Pathway]);stopifnot(length(ids)==31,nrow(sc)==9)
label<-env$label
sdraw<-function(){grid.newpage();xs<-seq(101,168,length.out=9);txt('GSEA display sensitivity',2,131,7,'bold');txt('Overlap',75,124,5.6,just='right');txt('Count',75,120,5.6,just='right')
 for(k in 1:9){txt(sprintf('%.2f',sc$Threshold[k]),xs[k],124,5.5,just='center');txt(sc$Budget[k],xs[k],120,5.5,just='center')}
 for(i in seq_along(ids)){y<-114-(i-1)*3.45;txt(label(ids[i]),2,y,5.7,if(ids[i]%in%mainids)'bold'else'plain');for(k in 1:9){yes<-nrow(sens[Pathway==ids[i]&Threshold==sc$Threshold[k]&Budget==sc$Budget[k]])>0;rect(xs[k],y,5.5,2.6,if(yes)'#506D82'else'#EEF1F3')}}
 txt('Bold: main18; filled: selected. Nested display budgets are not independent validation.',2,3,5.5,col='#5A656F')}
emit('GSEA_selection_sensitivity',sdraw,174,135)
bgs<-full[,.(Tested=uniqueN(pathway),SignificantAny=uniqueN(pathway[padj<.05])),by=Collection];bgs[,Main:=vapply(Collection,function(cc)sum(sel$Collection==cc),integer(1))];fwrite(bgs,file.path(out,'source_data/GSEA_background_unique_sets.csv'))
bdraw<-function(){grid.newpage();txt('GSEA selection background',1,109,7,'bold');txt('Gene sets',1,99,5.8,'bold');txt('Tested',32,99,5.6,'bold','center');txt('Any sig.',48,99,5.6,'bold','center');txt('Main',63,99,5.6,'bold','center')
 for(i in 1:4){y<-91-(i-1)*7;q<-bgs[Collection==collections[i]];if(i%%2)rect(35,y,69,6.5,'#F3F5F7');txt(collection_labels[i],1,y,5.8);txt(q$Tested,32,y,5.8,just='center');txt(q$SignificantAny,48,y,5.8,just='center');txt(q$Main,63,y,5.8,just='center')}
 txt('Any sig.: original FDR < 0.05 in >=1 subtype',1,61,5.3);txt('Testing families: subtype x collection',1,56,5.3)
 for(i in 1:6)txt(c('Selection rule','1  Six subtype x NES-direction queues','2  Rank by FDR, then absolute NES','3  Exclude measured-set Jaccard >=0.50','4  Three rounds: 18 main sets','5  Fourth round: 24-set extension')[i],1,47-(i-1)*6,if(i==1)6 else 5.5,if(i==1)'bold'else'plain')
 txt('Complete original NES/FDR and leading edges',1,10,5.2);txt('are provided for all tested sets in Source Data.',1,6,5.2)
}
emit('GSEA_background_compact',bdraw,70,114)
# Reorganize the same99TME effects by scientific question; no cross-method averaging.
em<-fread('8-Visualization/refinements/fig2_layout_audit_20261031_v31/E/tables/E_selection_map.csv');ez<-fread('8-Visualization/refinements/fig2_layout_audit_20261031_v31/E/tables/E_selected_contrasts.csv')
labels<-c('Immune score','Stromal score','Purity proxy','CD8 T cells','NK cells','B cells','Monocytes / lineage','Macrophages','Dendritic cells','Neutrophils','Fibroblasts / CAFs','Endothelial cells')
ys<-c(63,60.1,57.2,50.5,47.6,44.7,41.8,38.9,36,33.1,26.4,23.5);em[,DisplayY:=ys[match(Label,labels)]];ez[,DisplayY:=ys[match(Label,labels)]];stopifnot(!anyNA(ez$DisplayY),nrow(ez)==99)
sh<-c('MCP-counter'=21,EPIC=22,xCell=24,CIBERSORT=23,ESTIMATE=25);off<-c('MCP-counter'=.8,EPIC=.267,xCell=-.267,CIBERSORT=-.8,ESTIMATE=0)
draw_e<-function(){grid.newpage();cx<-c(65,86,107);txt('Biological question',2,73,5.5,'bold');txt('Readout',25,73,5.6,'bold');for(j in 1:3){txt(names(pal)[j],cx[j],73,6.3,'bold','center');rect(cx[j],70,18,1.1,pal[j]);line(rep(cx[j],2),c(21,65),'#CDD3D8',.5,2)}
 blocks<-list(list(l=labels[1:3],title='Tissue\nbackground'),list(l=labels[4:10],title='Immune-cell\nsignals'),list(l=labels[11:12],title='Stromal /\nvascular signals'))
 for(i in 1:3){yy<-ys[match(blocks[[i]]$l,labels)];lo<-min(yy)-1.45;hi<-max(yy)+1.45;if(i!=2)rect(60,(lo+hi)/2,119,hi-lo,'#F3F5F7');txt(blocks[[i]]$title,2,(lo+hi)/2,5.3,'bold')}
 for(i in seq_along(labels))txt(labels[i],25,ys[i],5.3)
 dot<-function(x,y,method,fill,col)grid.points(mm(x),mm(y),pch=sh[method],size=mm(1.05),gp=gpar(col=col,fill=fill,lwd=.6))
 for(i in seq_len(nrow(ez))){r<-ez[i];j<-match(r$Subtype,names(pal));y<-r$DisplayY+off[r$Method];xx<-function(v)cx[j]+4.4*v;stopifnot(r$Lower>=-2,r$Upper<=2,r$Estimable);line(xx(c(r$Lower,r$Upper)),c(y,y),adjustcolor(pal[j],alpha.f=.5),.6);dot(xx(r$Effect),y,r$Method,if(r$FDR<.05)pal[j]else'white',pal[j])}
 for(j in 1:3){line(cx[j]+c(-8.8,8.8),c(20.5,20.5),'grey50');for(v in c(-2,0,2))txt(v,cx[j]+v*4.4,18,5,just='center')};txt('Study-adjusted difference (within-study SD)',86,14.5,5.2,just='center')
 xx<-c(3,30,49,70,98);for(i in 1:5){dot(xx[i],9,names(sh)[i],'#444C53','#444C53');txt(names(sh)[i],xx[i]+2,9,5.2)};dot(3,4,'MCP-counter','#444C53','#444C53');txt('FDR <0.05',5,4,5.2);dot(32,4,'MCP-counter','white','#444C53');txt('FDR >=0.05',34,4,5.2);line(c(65,71),c(4,4),'grey50');txt('95% CI',74,4,5.2)
}
emit('Figure2E_question_order',draw_e,120,76);fwrite(ez,file.path(out,'source_data/TME_all99_question_order.csv'));fwrite(em,file.path(out,'source_data/TME_selection_question_order.csv'))
