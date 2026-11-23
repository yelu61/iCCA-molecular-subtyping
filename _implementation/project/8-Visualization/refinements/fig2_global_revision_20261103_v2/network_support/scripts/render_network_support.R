# Scientific support from frozen results only; no network, effect or meta-analysis fitting.
suppressPackageStartupMessages({library(data.table);library(grid);library(svglite);library(ragg);library(jsonlite)})
source('8-Visualization/scripts/fig_style.R')
out <- '8-Visualization/refinements/fig2_global_revision_20261103_v2/network_support'
for(d in c('figures','source_data','review','scripts'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
base <- '2-TranscriptomicCharacterization/publication/'
inputs <- c(main='8-Visualization/refinements/fig2_selection_20261031_v32/tables/C_all11_modules.csv',
  order='8-Visualization/refinements/sfig2_compact_20261031_v35/tables/module_order.csv',
  study=paste0(base,'tables/step2_wgcna_subtype_study_effects.csv'),meta=paste0(base,'tables/step2_wgcna_subtype_random_effects_meta.csv'),
  preservation=paste0(base,'tables/step2_wgcna_module_preservation.csv'),counts=paste0(base,'tables/step2_wgcna_study_subtype_counts.csv'),
  contract=paste0(base,'tables/step2_wgcna_analysis_contract.csv'),gate=paste0(base,'tables/step2_wgcna_release_gate_status.csv'),
  sensitivity=paste0(base,'tables/step2_wgcna_reference_outlier_network_sensitivity.csv'),
  model=paste0(base,'models/step2_wgcna_network.rds'),master='governance/data/subtype_master_table.csv',
  analysis_script='2-TranscriptomicCharacterization/scripts/run_step2_wgcna_refresh.R',
  subtype_palette=SUBTYPE_PALETTE_CSV,study_palette=STUDY_PALETTE_CSV,
  style='8-Visualization/scripts/fig_style.R',style_config='8-Visualization/scripts/figure_config.json')
hashes<-unname(tools::md5sum(inputs));stopifnot(!anyNA(hashes))
read<-function(k)fread(inputs[[k]])
main<-read('main');setorder(main,DisplayOrder);ord<-read('order');stud<-read('study');meta<-read('meta');pres<-read('preservation');counts<-read('counts');contract<-read('contract');gate<-read('gate');sensitivity<-read('sensitivity')
pal<-load_subtype_palette()$colors;study_pal<-load_study_palette()$colors
model<-readRDS(inputs['model']);master<-read('master')[SubtypeClass=='Main',.(Sample=as.character(Sample),FinalSubtype,Study)]
stopifnot(nrow(main)==11,nrow(ord)==24,nrow(pres)==24,nrow(counts)==7,nrow(master)==497,
  nrow(gate[Blocking==TRUE & !Status%in%c('PASS','QUALIFIED')])==0,
  contract[Parameter=='subtype_estimand',Value]=='target subtype minus equal-weighted mean of the other two',
  all(pres$ModuleSize==pres$PreservationEvaluationSize),length(model$reference_samples)==220,length(model$validation_samples)==274)
for(k in setdiff(names(inputs),c('model','master')))file.copy(inputs[k],file.path(out,'source_data',paste0(k,'_',basename(inputs[k]))),overwrite=TRUE)
fwrite(master,file.path(out,'source_data/main_sample_identity.csv'))
fwrite(data.table(Source=unname(inputs),MD5=hashes),file.path(out,'source_data/input_hashes.csv'))

# Algebraic read-back of frozen scores verifies equal-weighted other-type means.
# This is a source-integrity audit; source EffectSD/SE/CI are used unchanged for plotting.
scoremat<-model$fixed_loading_projected_study_standardized_scores
ids<-match(rownames(scoremat),master$Sample);stopifnot(!anyNA(ids));scoremeta<-master[ids]
assigned<-merge(stud,main[,.(Module,AssignedSubtypeOrPair,DisplayOrder,ModuleSize)],by='Module',sort=FALSE)[TargetSubtype==AssignedSubtypeOrPair]
assigned_meta<-merge(meta,main[,.(Module,AssignedSubtypeOrPair,DisplayOrder,ModuleSize)],by='Module',sort=FALSE)[TargetSubtype==AssignedSubtypeOrPair]
study_order<-names(study_pal)[names(study_pal)%in%counts$Study]
stopifnot(length(study_order)==7)
assigned[,StudyOrder:=match(Study,study_order)];setorder(assigned,DisplayOrder,StudyOrder);setorder(assigned_meta,DisplayOrder)
stopifnot(nrow(assigned)==77,nrow(assigned[EligibleForInference==TRUE])==66,all(assigned[EligibleForInference==FALSE,Study]=='2023Lee'),all(assigned_meta$DirectionEligibleStudies==6))
audit<-rbindlist(lapply(seq_len(nrow(assigned)),function(i){r<-assigned[i];take<-which(scoremeta$Study==r$Study);g<-scoremeta$FinalSubtype[take];v<-scoremat[take,r$Module];means<-tapply(v,g,mean);nn<-table(factor(g,levels=c('C1','C2','C3')));other<-setdiff(c('C1','C2','C3'),r$TargetSubtype)
  eq<-means[r$TargetSubtype]-mean(means[other]);pooled<-mean(v[g==r$TargetSubtype])-mean(v[g!=r$TargetSubtype]);se<-if(r$EligibleForInference)sqrt(var(v[g==r$TargetSubtype])/nn[r$TargetSubtype]+sum(vapply(other,function(x)var(v[g==x])/nn[x],numeric(1)))/4)else NA_real_
  stopifnot(abs(eq-r$EffectSD)<1e-12,as.integer(nn[r$TargetSubtype])==r$N_Target,as.integer(nn[other[1]])==r$N_Other1,as.integer(nn[other[2]])==r$N_Other2)
  if(r$EligibleForInference)stopifnot(abs(se-r$SE)<1e-12,abs(r$CI_low-(r$EffectSD-1.96*r$SE))<1e-12,abs(r$CI_high-(r$EffectSD+1.96*r$SE))<1e-12)
  data.table(Module=r$Module,Study=r$Study,TargetSubtype=r$TargetSubtype,SourceEffect=r$EffectSD,EqualWeightedEffect=unname(eq),SamplePooledOthersEffect=pooled,EqualWeightedDifference=unname(eq-r$EffectSD),PooledOthersDifference=pooled-r$EffectSD,SourceSE=r$SE,ReadbackSE=unname(se),EligibleForInference=r$EligibleForInference)
}))
fwrite(audit,file.path(out,'source_data/contrast_definition_audit.csv'));fwrite(assigned,file.path(out,'source_data/assigned_study_effects.csv'));fwrite(assigned_meta,file.path(out,'source_data/assigned_meta_effects.csv'))
pres<-merge(ord[,.(Module,RowOrder,GroupOrder,MainFigure)],pres,by='Module',sort=FALSE);setorder(pres,RowOrder);fwrite(pres,file.path(out,'source_data/preservation_ordered.csv'))
meta_counts<-counts[Study!='2023Lee',lapply(.SD,sum),.SDcols=c('C1','C2','C3','Total')]
stopifnot(meta_counts$Total==484)

mm<-function(x)unit(x,'mm');texts<-list();marks<-list();rows<-list()
txt<-function(s,x,y,size=5.5,face='plain',just='left',col='#28333B',rot=0){grid.text(s,mm(x),mm(y),just=just,rot=rot,gp=gpar(fontfamily='Arial',fontsize=size,fontface=face,col=col));texts[[length(texts)+1L]]<<-data.table(Text=as.character(s),X=x,Y=y,Size=size)}
line<-function(x,y,col='#BCC5CB',lwd=.4,lty=1)grid.lines(mm(x),mm(y),gp=gpar(col=col,lwd=lwd,lty=lty))
rect<-function(x,y,w,h,fill,border=NA)grid.rect(mm(x),mm(y),mm(w),mm(h),just=c('left','bottom'),gp=gpar(fill=fill,col=border,lwd=.4))
point<-function(x,y,size=1,pch=16,col='#273944',fill=col)grid.points(mm(x),mm(y),pch=pch,size=mm(size),gp=gpar(col=col,fill=fill,lwd=.6))
clear<-function(){texts<<-list();marks<<-list();rows<<-list()}
fmt<-function(v)sprintf('%.2f',ifelse(abs(v)<.005,0,v))
ci_fmt<-function(e,l,u)paste0(fmt(e),' [',fmt(l),', ',fmt(u),']')
short_study<-c('2019Ahn'='Ahn','2022Dong'='Dong*','2023Deng'='Deng','2023Jeon'='Jeon','2023Lee'='Lee','2024Fan'='Fan','TCGA-CHOL'='TCGA')
axis_min<--1.25;axis_max<-3

draw_forest<-function(){grid.newpage();clear()
  txt('Assigned-subtype effects across studies',2,247.1,8,'bold')
  txt('Target − equal-weighted mean of the other two subtypes; frozen effects and 95% CIs',2,242.7,6.2)
  for(i in 1:11){m<-main[i];x0<-2+((i-1)%%2)*87;top<-237-floor((i-1)/2)*38.5
    r<-assigned[Module==m$Module];a<-assigned_meta[Module==m$Module];color<-pal[m$AssignedSubtypeOrPair]
    rect(x0,top-3.3,1,2.1,m$ModuleColor);txt(paste0(m$ModuleColor,' · ',m$AssignedSubtypeOrPair,' · ',m$ModuleSize,' genes'),x0+2,top-2.25,6.35,'bold')
    txt('Study',x0,top-5.8,5.5,'bold');txt('n (1/2/3)',x0+22.1,top-5.8,5.5,'bold','center');txt('Effect (SD)',x0+47.0,top-5.8,5.5,'bold','center');txt('Estimate [95% CI]',x0+83,top-5.8,5.5,'bold','right')
    fx<-function(v)x0+34+(v-axis_min)/(axis_max-axis_min)*27
    ys<-top-9.0-(0:6)*2.6;ym<-top-27.6
    line(rep(fx(0),2),c(ym-1.15,ys[1]+1.1),'#CBD2D8',.45,2)
    for(j in 1:7){s<-r[Study==study_order[j]];n<-counts[Study==s$Study];y<-ys[j]
      txt(short_study[s$Study],x0,y,5.5);txt(paste(n$C1,n$C2,n$C3,sep='/'),x0+22.1,y,5.5,just='center')
      if(s$EligibleForInference){line(fx(c(s$CI_low,s$CI_high)),c(y,y),color,.7);point(fx(s$EffectSD),y,size=1.0,pch=if(s$Study=='2022Dong')15 else 16,col=color);txt(ci_fmt(s$EffectSD,s$CI_low,s$CI_high),x0+83,y,5.5,just='right')}
      else {txt('NE',x0+47,y,5.5,just='center',col='#78838C');txt(paste0(fmt(s$EffectSD),' [NE]'),x0+83,y,5.5,just='right',col='#78838C')}
      marks[[length(marks)+1L]]<<-data.table(Kind='Study',Module=m$Module,Study=s$Study,TargetSubtype=m$AssignedSubtypeOrPair,Effect=s$EffectSD,CI_low=s$CI_low,CI_high=s$CI_high,Eligible=s$EligibleForInference,N_C1=n$C1,N_C2=n$C2,N_C3=n$C3,X=if(s$EligibleForInference)fx(s$EffectSD)else NA_real_,CI_X_low=if(s$EligibleForInference)fx(s$CI_low)else NA_real_,CI_X_high=if(s$EligibleForInference)fx(s$CI_high)else NA_real_,Y=y)
    }
    line(c(x0,x0+83),rep(ym+1.4,2),'#E2E7EB',.35)
    txt('Meta',x0,ym,5.5,'bold');txt(paste(meta_counts$C1,meta_counts$C2,meta_counts$C3,sep='/'),x0+22.1,ym,5.5,'bold','center');line(fx(c(a$CI_low,a$CI_high)),rep(ym,2),color,.9);point(fx(a$RandomEffect),ym,size=1.35,pch=18,col=color);txt(ci_fmt(a$RandomEffect,a$CI_low,a$CI_high),x0+83,ym,5.5,'bold','right')
    marks[[length(marks)+1L]]<<-data.table(Kind='Meta',Module=m$Module,Study='Meta',TargetSubtype=m$AssignedSubtypeOrPair,Effect=a$RandomEffect,CI_low=a$CI_low,CI_high=a$CI_high,Eligible=TRUE,N_C1=meta_counts$C1,N_C2=meta_counts$C2,N_C3=meta_counts$C3,X=fx(a$RandomEffect),CI_X_low=fx(a$CI_low),CI_X_high=fx(a$CI_high),Y=ym)
    ay<-top-30.0;line(fx(c(axis_min,axis_max)),c(ay,ay),'#74818B',.4)
    for(v in -1:3){line(rep(fx(v),2),c(ay,ay-.45),'#74818B',.4);txt(v,fx(v),ay-1.65,5.5,just='center')}
    txt(paste0('I² = ',sprintf('%.1f',a$I2Percent),'%;  meta FDR = ',formatC(a$FDR,format='e',digits=1)),x0,top-35.15,5.5)
    rows[[length(rows)+1L]]<<-data.table(Module=m$Module,AssignedSubtype=m$AssignedSubtypeOrPair,Row=floor((i-1)/2)+1,Column=((i-1)%%2)+1,I2Percent=a$I2Percent,FDR=a$FDR,ModuleSize=m$ModuleSize)
  }
  x<-89;top<-44.5
  txt('Reading the forest',x,top-2.25,6.35,'bold')
  note<-c('n (1/2/3) gives C1 / C2 / C3 sample counts.',
    'Circles: study effects; square: Dong reference.',
    'Diamonds: original six-study REML meta-analysis.',
    'Dong*: 223 scored; 220 used to construct the network.',
    'Lee: C2 n = 1; point estimate retained, inference NE.',
    'Meta n = 484; all 7 studies are displayed.',
    'I² is for the assigned subtype contrast in that block.',
    'Across-study support is not independent validation.')
  for(j in seq_along(note))txt(note[j],x,top-7.0-(j-1)*3.3,5.5)
  txt('Same score units and x-axis across all modules; NE does not denote a zero effect.',2,2.6,5.5)
}

draw_preservation<-function(wide=FALSE){grid.newpage();clear()
  txt('Structural preservation of all 24 modules',2,126.6,7,'bold')
  txt('Reference n = 220; pooled held-out n = 274 across 6 studies',2,122.5,5.8)
  txt('Module',2,117.7,5.6,'bold');txt('Genes',25.8,117.7,5.4,'bold','right')
  fields<-c('ZsummaryPreservation','ZdensityPreservation','ZconnectivityPreservation');labs<-c('Zsummary','Zdensity','Zconnectivity');starts<-if(wide)c(35,81,127)else c(30.7,56.7,82.7);numbers<-if(wide)c(77,123,169)else c(53.7,79.7,105.7);pw<-if(wide)33 else 16.8
  for(k in 1:3){txt(labs[k],starts[k]+pw/2,117.7,5.7,'bold','center');line(c(starts[k],starts[k]+pw),rep(114.5,2),'#7D8991',.4)
    for(v in c(0,20,40)){x<-starts[k]+pw*v/45;line(c(x,x),c(114.5,114.0),'#7D8991',.4);txt(v,x,112.5,5.2,just='center')}
    line(rep(starts[k],2),c(28.6,110.0),'#EDF0F2',.35)
  }
  for(v in c(2,10)){x<-starts[1]+pw*v/45;line(rep(x,2),c(28.6,110.0),if(v==10)'#B5C0C8'else'#D6DDE2',.45,2)}
  for(i in 1:24){r<-pres[i];y<-108.3-(i-1)*3.4
    if(i>1 && r$GroupOrder!=pres$GroupOrder[i-1])line(c(2,if(wide)172 else 108),rep(y+1.7,2),'#E1E6EA',.35)
    rect(2,y-1.2,1.1,2.4,r$ModuleColor);txt(r$ModuleColor,4,y,5.6,if(r$MainFigure)'bold'else'plain');txt(r$ModuleSize,25.8,y,5.3,just='right');if(r$MainFigure)point(27.3,y,.6)
    for(k in 1:3){z<-r[[fields[k]]];x<-starts[k]+pw*z/45;point(x,y,1.05,pch=if(k==1 && z<10)1 else 16,col=if(k==1)'#273D4E'else'#657886');txt(sprintf('%.2f',z),numbers[k],y,5.5,just='right')
      marks[[length(marks)+1L]]<<-data.table(Module=r$Module,Metric=fields[k],Value=z,X=x,Y=y,ModuleSize=r$ModuleSize,EvaluationSize=r$PreservationEvaluationSize,MainFigure=r$MainFigure)
    }
  }
  txt('200 permutations; every module member evaluated.',2,24.9,5.6)
  point(2.7,20.7,.6);txt('Main Fig. 2C module; open Zsummary dots: below 10.',4.2,20.7,5.4)
  txt('Dashed Zsummary guides: 2 and 10; no component thresholds.',2,16.7,5.4)
  txt('Pooled structural support; not independent validation.',2,12.7,5.5)
  txt('Reference outlier sensitivity gate: QUALIFIED (see caption).',2,8.7,5.5)
  txt('Z statistics depend on module size and pooled-study composition.',2,4.7,5.4)
}
export<-function(name,draw,w,h){p<-file.path(out,'figures',name)
  svglite(paste0(p,'.svg'),width=w/25.4,height=h/25.4,bg='white');draw();dev.off()
  cairo_pdf(paste0(p,'.pdf'),width=w/25.4,height=h/25.4,family='Arial',bg='white');draw();dev.off()
  agg_png(paste0(p,'.png'),width=w,height=h,units='mm',res=600,background='white');draw();dev.off()
  agg_png(file.path(out,'review',paste0(name,'_preview.png')),width=w,height=h,units='mm',res=300,background='white');draw();dev.off()
  fwrite(rbindlist(texts),file.path(out,'review',paste0(name,'_texts.csv')));fwrite(rbindlist(marks),file.path(out,'source_data',paste0(name,'_marks.csv')))
  if(length(rows))fwrite(rbindlist(rows),file.path(out,'source_data',paste0(name,'_blocks.csv')))
}
export('WGCNA_assigned_contrast_study_forests',draw_forest,174,250)
export('WGCNA_all24_preservation',draw_preservation,110,130)
export('WGCNA_all24_preservation_wide',function()draw_preservation(TRUE),174,130)
stopifnot(identical(hashes,unname(tools::md5sum(inputs))))
write_json(list(Status='RENDERED_PENDING_QA',Forest=list(WidthMM=174,HeightMM=250,Modules=11,StudyRows=77,EligibleStudyRows=66,NERows=11,MetaRows=11,MetaN=484),Preservation=list(WidthMM=110,HeightMM=130,WideWidthMM=174,WideHeightMM=130,Modules=24,Metrics=3,ReferenceN=220,HeldOutN=274,Permutations=200),Contrast='target minus equal-weighted mean of other two',Refitted=FALSE,InputsUnchanged=TRUE,MaxContrastReadbackDifference=max(abs(audit$EqualWeightedDifference))),file.path(out,'review/render_summary.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'review/sessionInfo.txt'))
