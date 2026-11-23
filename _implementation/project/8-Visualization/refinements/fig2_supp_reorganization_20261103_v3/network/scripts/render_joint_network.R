# Display-only union of frozen meta-effects and network-preservation statistics.
suppressPackageStartupMessages({library(data.table);library(grid);library(svglite);library(ragg);library(jsonlite)})
source('8-Visualization/scripts/fig_style.R')
out<-'8-Visualization/refinements/fig2_supp_reorganization_20261103_v3/network'
for(d in c('scripts','figures','source_data','review'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
v1<-'8-Visualization/refinements/fig2_global_revision_20261103_v1/principal_support/tables/'
v2<-'8-Visualization/refinements/fig2_global_revision_20261103_v2/network_support/source_data/'
inputs<-c(effects=paste0(v1,'network_principal_effects.csv'),order=paste0(v1,'module_order.csv'),
  preservation=paste0(v2,'preservation_ordered.csv'),meta=paste0(v2,'meta_step2_wgcna_subtype_random_effects_meta.csv'),
  canonical_preservation=paste0(v2,'preservation_step2_wgcna_module_preservation.csv'),
  contract=paste0(v2,'contract_step2_wgcna_analysis_contract.csv'),gate=paste0(v2,'gate_step2_wgcna_release_gate_status.csv'),
  sensitivity=paste0(v2,'sensitivity_step2_wgcna_reference_outlier_network_sensitivity.csv'),
  subtype_palette=SUBTYPE_PALETTE_CSV,study_palette=STUDY_PALETTE_CSV,
  style='8-Visualization/scripts/fig_style.R',style_config='8-Visualization/scripts/figure_config.json')
hashes<-unname(tools::md5sum(inputs));stopifnot(!anyNA(hashes))
read<-function(k)fread(inputs[[k]])
effects<-read('effects');ord<-read('order');pres<-read('preservation');meta<-read('meta');cp<-read('canonical_preservation');contract<-read('contract');gate<-read('gate')
pal<-load_subtype_palette()$colors;study_pal<-load_study_palette()$colors
setorder(ord,RowOrder);setorder(pres,RowOrder)
stopifnot(nrow(ord)==24,nrow(effects)==72,nrow(pres)==24,sum(ord$MainFigure)==11,
  identical(ord$Module,pres$Module),identical(ord$ModuleSize,pres$ModuleSize),
  gate[Check=='reference_connectivity_outlier_check',Status]=='QUALIFIED',
  contract[Parameter=='subtype_estimand',Value]=='target subtype minus equal-weighted mean of the other two')
joined<-merge(effects,meta,by=c('Module','TargetSubtype'),suffixes=c('','.canonical'))
for(f in c('RandomEffect','FDR','CI_low','CI_high','I2Percent'))stopifnot(max(abs(joined[[f]]-joined[[paste0(f,'.canonical')]]))<1e-12)
joinedp<-merge(pres,cp,by='Module',suffixes=c('','.canonical'))
for(f in c('ZsummaryPreservation','ZdensityPreservation','ZconnectivityPreservation'))stopifnot(max(abs(joinedp[[f]]-joinedp[[paste0(f,'.canonical')]]))<1e-12)
for(k in names(inputs))file.copy(inputs[k],file.path(out,'source_data',paste0(k,'_',basename(inputs[k]))),overwrite=TRUE)
fwrite(data.table(Source=unname(inputs),MD5=hashes),file.path(out,'source_data/input_hashes.csv'))
source_joint<-merge(effects,pres[,.(Module,PreservationEvaluationSize,ZdensityPreservation,ZconnectivityPreservation,MedianRankPreservation)],by='Module',sort=FALSE)
source_joint[,SubtypeOrder:=match(TargetSubtype,c('C1','C2','C3'))];setorder(source_joint,RowOrder,SubtypeOrder)
fwrite(source_joint,file.path(out,'source_data/all72_effects_and_preservation.csv'))
mm<-function(x)unit(x,'mm');texts<-list();effect_marks<-list();pres_marks<-list()
txt<-function(s,x,y,size=6,face='plain',just='left',col='#28333B'){grid.text(s,mm(x),mm(y),just=just,gp=gpar(fontfamily='Arial',fontsize=size,fontface=face,col=col));texts[[length(texts)+1L]]<<-data.table(Text=as.character(s),X=x,Y=y,Size=size)}
line<-function(x,y,col='#CAD2D8',lwd=.4,lty=1)grid.lines(mm(x),mm(y),gp=gpar(col=col,lwd=lwd,lty=lty))
rect<-function(x,y,w,h,fill,border=NA)grid.rect(mm(x),mm(y),mm(w),mm(h),just=c('left','bottom'),gp=gpar(fill=fill,col=border,lwd=.35))
point<-function(x,y,size=.65,pch=16,col='#273D4E')grid.points(mm(x),mm(y),pch=pch,size=mm(size),gp=gpar(col=col,lwd=.6))
ramp<-colorRamp(c(SCALE_EFFECT_LOW,SCALE_EFFECT_MID,SCALE_EFFECT_HIGH));effect_color<-function(v)rgb(ramp((pmax(-2,pmin(2,v))+2)/4),maxColorValue=255)
num<-function(v)sprintf('%.2f',ifelse(abs(v)<.005,0,v))
draw<-function(){grid.newpage();texts<<-list();effect_marks<<-list();pres_marks<<-list()
  txt('Co-expression associations and network preservation',2,116.5,7.5,'bold')
  txt('Meta: 6 eligible studies (n = 484)   |   Preservation: reference n = 220; pooled held-out n = 274',2,112.4,6)
  txt('Module',2,106.8,6.2,'bold');txt('Genes',28.8,106.8,5.8,'bold','right')
  xs<-c(39,53,67);cw<-12.5
  for(j in 1:3){txt(paste0('C',j),xs[j],106.8,6.3,'bold','center');rect(xs[j]-cw/2,103.4,cw,1,pal[paste0('C',j)])}
  txt('Zsummary',104,106.8,6.2,'bold','center');txt('Zdensity',147.2,106.8,6.2,'bold','center');txt('Zconnectivity',166.4,106.8,6.2,'bold','center')
  fx<-function(v)84+39*v/45
  line(fx(c(0,45)),c(104.2,104.2),'#78858F',.4)
  for(v in c(0,20,40)){line(rep(fx(v),2),c(104.2,103.7),'#78858F',.4);txt(v,fx(v),102.1,5.5,just='center')}
  for(v in c(2,10))line(rep(fx(v),2),c(19.2,100.3),if(v==10)'#B5C1CA'else'#D9E0E5',.45,2)
  for(i in 1:24){m<-ord[i];p<-pres[i];y<-98.9-(i-1)*3.4
    if(i>1 && m$GroupOrder!=ord$GroupOrder[i-1])line(c(2,172),rep(y+1.7,2),'#E2E7EA',.35)
    rect(2,y-1.2,1.2,2.4,m$ModuleColor);txt(m$ModuleColor,4,y,6,if(m$MainFigure)'bold'else'plain');txt(m$ModuleSize,28.8,y,5.8,just='right');if(m$MainFigure)point(30.2,y,.6)
    for(j in 1:3){e<-effects[Module==m$Module & TargetSubtype==paste0('C',j)];stopifnot(nrow(e)==1);color<-if(abs(e$RandomEffect)>1.3)'white'else'#24313B'
      rect(xs[j]-cw/2,y-1.55,cw,3.1,effect_color(e$RandomEffect));txt(num(e$RandomEffect),xs[j]-.65,y,5.9,just='center',col=color);if(e$FDR<.05)point(xs[j]+4.6,y,.62,col=color)
      effect_marks[[length(effect_marks)+1L]]<<-data.table(Module=m$Module,TargetSubtype=e$TargetSubtype,RandomEffect=e$RandomEffect,FDR=e$FDR,Significant=e$FDR<.05,X=xs[j],Y=y)
    }
    point(fx(p$ZsummaryPreservation),y,1.1,pch=if(p$ZsummaryPreservation<10)1 else 16)
    txt(num(p$ZsummaryPreservation),132.1,y,5.9,just='right');txt(num(p$ZdensityPreservation),147.2,y,5.9,just='center');txt(num(p$ZconnectivityPreservation),166.4,y,5.9,just='center')
    pres_marks[[length(pres_marks)+1L]]<<-data.table(Module=m$Module,ModuleSize=m$ModuleSize,MainFigure=m$MainFigure,Zsummary=p$ZsummaryPreservation,Zdensity=p$ZdensityPreservation,Zconnectivity=p$ZconnectivityPreservation,EvaluationSize=p$PreservationEvaluationSize,X=fx(p$ZsummaryPreservation),Y=y)
  }
  txt('Effect (SD)',2,14.8,5.7)
  for(k in 1:128)rect(22+(k-1)*24/128,13.8,24/128+.01,1.35,effect_color(-2+4*(k-1)/127))
  for(v in c(-2,0,2))txt(v,22+(v+2)*6,11.6,5.4,just='center')
  point(53,14.6,.62);txt('Tile dot: original FDR < 0.05',55,14.6,5.7)
  point(53,9.8,.6);txt('Gene-count dot / bold name: main Fig. 2C module',55,9.8,5.7)
  txt('Zsummary guides: 2 and 10',116,14.6,5.7)
  txt('200 permutations; all module genes evaluated',116,9.8,5.7)
  txt('Reference QC gate: QUALIFIED; pooled structural support, not independent validation.',2,4.6,5.9)
}
name<-'WGCNA_all24_joint_association_preservation';p<-file.path(out,'figures',name)
svglite(paste0(p,'.svg'),width=174/25.4,height=120/25.4,bg='white');draw();dev.off()
cairo_pdf(paste0(p,'.pdf'),width=174/25.4,height=120/25.4,family='Arial',bg='white');draw();dev.off()
agg_png(paste0(p,'.png'),width=174,height=120,units='mm',res=600,background='white');draw();dev.off()
agg_png(file.path(out,'review',paste0(name,'_preview.png')),width=174,height=120,units='mm',res=300,background='white');draw();dev.off()
fwrite(rbindlist(texts),file.path(out,'review/text_geometry.csv'));fwrite(rbindlist(effect_marks),file.path(out,'source_data/plotted_effects.csv'));fwrite(rbindlist(pres_marks),file.path(out,'source_data/plotted_preservation.csv'))
stopifnot(identical(hashes,unname(tools::md5sum(inputs))))
write_json(list(Status='RENDERED_PENDING_QA',WidthMM=174,HeightMM=120,Modules=24,MainModules=11,MetaEffects=72,PreservationStatistics=72,SourceInputsUnchanged=TRUE,Refitted=FALSE),file.path(out,'review/render_summary.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'review/sessionInfo.txt'))
