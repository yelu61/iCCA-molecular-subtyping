# Figure6 integration: explicit clinical questions, paired identification and drug hypotheses.
# R-only independent PDF masters; old releases and canonical models are read-only.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(grid);library(digest)})
source('8-Visualization/scripts/fig_style.R')
out <- '8-Visualization/refinements/fig6_integrated_20261101_v20'
old <- '8-Visualization/refinements/fig6_storyline_20261101_v19'
joint <- '5-ClassifierValidation/review/20261101_duct_compact_joint_v1'
for(d in c('figures','supplement','review','tables','qa','provenance'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
inputs<-character(); manifest<-list(); checks<-list(); text_bounds<-list(); active<-''; checking<-FALSE
rd<-function(p){inputs<<-unique(c(inputs,p));cols<-names(fread(p,nrows=0));if('Sample'%in%cols)fread(p,colClasses=c(Sample='character'))else fread(p)}
ro<-function(n)rd(file.path(old,'tables',paste0(n,'.csv')))
check<-function(n,b){stopifnot(all(b));checks[[n]]<<-data.table(Check=n,Pass=TRUE)}
pal<-load_subtype_palette()$colors[c('C1','C2','C3')];invisible(load_study_palette());cls<-names(pal)
ink<-'#263E4B';muted<-'#7C8992'; light<-'#E5EBEF'
base<-theme_classic(base_size=6.4,base_family='Arial')+theme(axis.text=element_text(size=6,color=ink),axis.title=element_text(size=6.4),plot.title=element_text(size=7.2,face='bold'),plot.subtitle=element_text(size=6),plot.caption=element_text(size=5.8,hjust=0),plot.title.position='plot',plot.caption.position='plot',legend.position='bottom',legend.title=element_blank(),legend.text=element_text(size=6),strip.background=element_blank(),strip.text=element_text(size=6.3,face='bold'),panel.spacing=unit(1.2,'mm'),plot.margin=margin(1.3,1.3,1.3,1.3,'mm'))
theme_set(base)
start<-function(w,h){W<<-w;H<<-h;grid.newpage();pushViewport(viewport(xscale=c(0,w),yscale=c(h,0)))}
txt<-function(s,x,y,size=6.3,bold=FALSE,col=ink,just='left',italic=FALSE){
 g<-textGrob(s,unit(x,'native'),unit(y,'native'),just=c(just,'center'),gp=gpar(fontfamily='Arial',fontsize=size,fontface=if(italic)'italic' else if(bold)'bold' else 'plain',col=col,lineheight=1.03))
 if(checking){ww<-convertWidth(grobWidth(g),'mm',TRUE);hh<-convertHeight(grobHeight(g),'mm',TRUE);lo<-x-switch(just,left=0,center=ww/2,right=ww)
  text_bounds[[length(text_bounds)+1L]]<<-data.table(Panel=active,Text=s,XMin=lo,XMax=lo+ww,YMin=y-hh/2,YMax=y+hh/2,FontPt=size)
  if(lo< -.05||lo+ww>W+.05||y-hh/2< -.05||y+hh/2>H+.05)stop('Text outside panel: ',active,' / ',s)
 };grid.draw(g)
}
ln<-function(x,y,col=muted,lwd=.5,lty=1)grid.lines(unit(x,'native'),unit(y,'native'),gp=gpar(col=col,lwd=lwd,lty=lty))
dot<-function(x,y,col=ink,pch=16,size=1.25)grid.points(unit(x,'native'),unit(y,'native'),pch=pch,size=unit(size,'mm'),gp=gpar(col=col,fill=col,lwd=.6))
rect<-function(x,y,w,h,fill,col=NA,lwd=.5)grid.rect(unit(x,'native'),unit(y,'native'),unit(w,'native'),unit(h,'native'),gp=gpar(fill=fill,col=col,lwd=lwd))
axisx<-function(xx,ticks,y,labels=as.character(ticks),size=5.7){ln(xx(range(ticks)),c(y,y));for(i in seq_along(ticks)){ln(rep(xx(ticks[i]),2),c(y,y+.6));txt(labels[i],xx(ticks[i]),y+2,size,just='center')}}
interval<-function(r,xx,y,col=ink,pch=16,hr=FALSE){e<-if(hr)r$HR else r$Estimate;stopifnot(nrow(r)==1,is.finite(e),r$Lower95<=e,r$Upper95>=e);v<-xx(c(r$Lower95,e,r$Upper95));stopifnot(all(v>0&v<W));ln(v[c(1,3)],c(y,y),col,.75);dot(v[2],y,col,pch)}
save_panel<-function(id,p,w,h,d,track='Main',more=NULL){
 active<<-id; pdf<-file.path(out,if(track=='Main')'figures' else 'supplement',paste0(id,'.pdf'));png<-file.path(out,'review',paste0(id,'.png'));tb<-file.path(out,'tables',paste0(id,'.csv'))
 draw<-function(){if(is.function(p))p()else print(p)}
 checking<<-TRUE;cairo_pdf(pdf,width=w/25.4,height=h/25.4,family='Arial');draw();dev.off();checking<<-FALSE
 ragg::agg_png(png,width=w,height=h,units='mm',res=450,background='white');draw();dev.off();fwrite(d,tb)
 if(!is.null(more))for(n in names(more))fwrite(more[[n]],file.path(out,'tables',paste0(id,'_',n,'.csv')))
 manifest[[id]]<<-data.table(Panel=id,Track=track,PDF=pdf,PNG=png,SourceTable=tb,WidthMM=w,HeightMM=h,Rows=nrow(d),Reuse='Rendered from traceable source tables');message(id)
}

# A: prediction first; all original covariate/time-specific associations retained below.
ah<-ro('Figure6A_clinical_information_1');a<-ro('Figure6A_clinical_information_2');b<-ro('Figure6B_measured_duct_increment')
drawA<-function(){start(85,76);on.exit(popViewport());txt('Prognostic information',2,3.2,7.3,TRUE);txt('331 patients',83,3.2,6,just='right')
 txt('24-month OOF prediction',2,8.5,6.5,TRUE);txt('AUC',43,13,6.2,TRUE,just='center');txt('Increment',70,13,6.2,TRUE,just='center')
 ax<-function(v)32+(v-.5)/.35*18;dx<-function(v)60+(v+.02)/.20*15
 ln(rep(ax(.5),2),c(16,32),light,.6,2);ln(rep(dx(0),2),c(16,37.5),light,.6,2)
 for(i in 1:3){y<-c(18.5,24.5,30.5)[i];co<-c(muted,pal['C3'],ink)[i];cmp<-c('M0','MC3','MAll')[i];r<-a[Comparison==cmp]
 txt(c('Clinical','+ C3','+ C1/C2/C3')[i],2,y);interval(r,ax,y,co,c(1,16,15)[i]);txt(sprintf('%.3f',r$Estimate),57,y,6.2,just='right')
 if(i==1)txt('Ref.',83,y,6,just='right')else{r<-a[Comparison==paste0(cmp,'-M0')];interval(r,dx,y,co,c(1,16,15)[i]);txt(sprintf('%+.3f',r$Estimate),83,y,6.2,just='right')}
 };r<-a[Comparison=='MAll-MC3'];txt('Full split vs C3',2,36.5,6);interval(r,dx,36.5,muted,18);txt(sprintf('%+.3f',r$Estimate),83,36.5,6,just='right')
 axisx(ax,c(.5,.8),39,c('0.5','0.8'));axisx(dx,c(0,.15),39,c('0','0.15'))
 ln(c(2,83),c(44,44),light);txt('Adjusted associations',2,48,6.5,TRUE);dot(52,48,ink,16);txt('12 mo',54,48,5.8);dot(68,48,ink,1);txt('24 mo',70,48,5.8)
 hx<-function(v)37+(log(v)-log(.25))/(log(4.5)-log(.25))*34;ln(rep(hx(1),2),c(51,70),light,.6,2)
 keys<-c('Age (per 10 years)','Male vs female','Stage III/IV vs I/II','C2 vs C1','C3 vs C1')
 for(i in 1:5){y<-52+(i-1)*4.3;txt(keys[i],2,y,6);co<-if(i<=3)ink else pal[c('C2','C3')[i-3]]
  for(j in 1:2){r<-ah[Contrast==keys[i]&TimeMonths==c(12,24)[j]];yy<-y+c(-.8,.8)[j];interval(r,hx,yy,co,c(16,1)[j],TRUE)}
 };axisx(hx,c(.25,1,4),71.5,c('0.25','1','4'));txt('HR',83,73.5,5.7,just='right')
}
save_panel('Figure6A_clinical_information',drawA,85,76,a,more=list(associations=ah))
check('A retains all six prediction estimates and ten association estimates',nrow(a)==6&nrow(ah)==10)

# B: beyond measured duct is the primary visual question; absolute performance remains visible.
drawB<-function(){start(85,76);on.exit(popViewport());txt('Beyond measured duct',2,3.2,7.3,TRUE);txt('79 patients',83,3.2,6,just='right')
 dot(35,8.5,pal['C3'],16);txt('Primary',37,8.5,5.8);dot(55,8.5,ink,2);txt('Time-varying',57,8.5,5.8)
 txt('Change in 24-month AUC',2,13.5,6.5,TRUE);dx<-function(v)43+(v+.05)/.35*29;ln(rep(dx(0),2),c(16,36),light,.6,2)
 for(i in 1:3){y<-19+(i-1)*7;txt(c('Add C3 to clinical + duct','C3 vs duct (with clinical)','Add duct to clinical + C3')[i],2,y,6.0,bold=i==1)
  for(j in 1:2){r<-b[Comparison==c('M3-M1','M2-M1','M3-M2')[i]&Family==c('RidgePH','RidgeTimeVarying')[j]];yy<-y+c(-1.1,1.1)[j];interval(r,dx,yy,c(pal['C3'],ink)[j],c(16,2)[j]);txt(sprintf('%+.3f',r$Estimate),83,yy,5.8,just='right')}
 };axisx(dx,c(0,.15,.30),37,c('0','0.15','0.30'));ln(c(2,83),c(42,42),light)
 txt('Absolute 24-month AUC',2,46,6.5,TRUE);ax<-function(v)36+(v-.2)/.6*36;ln(rep(ax(.5),2),c(49,69),light,.6,2)
 for(i in 1:4){y<-51+(i-1)*5.7;txt(c('Clinical','Clinical + duct','Clinical + C3','Clinical + duct + C3')[i],2,y,6.0)
  for(j in 1:2){r<-b[Comparison==paste0('M',i-1)&Family==c('RidgePH','RidgeTimeVarying')[j]];yy<-y+c(-.9,.9)[j];interval(r,ax,yy,c(pal['C3'],ink)[j],c(16,2)[j]);if(j==1)txt(sprintf('%.3f',r$Estimate),83,y,6,just='right')}
 };axisx(ax,c(.2,.5,.8),71.5,c('0.2','0.5','0.8'))
}
save_panel('Figure6B_measured_duct_increment',drawB,85,76,b);check('B retains all fourteen estimates including sensitivity',nrow(b)==14)

# C: only measured pathology in the main comparison. Within-column fractions.
csrc<-ro('Figure6C_duct_proportions_1');cp<-csrc[Source=='Pathology'];cp<-melt(cp,id.vars=c('FinalSubtype','N'),measure.vars=c('Small','Large'),variable.name='Duct',value.name='Count');cp[,Fraction:=Count/N]
drawC<-function(){start(48,58);on.exit(popViewport());txt('Morphology and subtype',2,3.4,7.0,TRUE);txt('Measured pathology · n = 93',2,8,5.9)
 for(j in 1:3){xx<-18+(j-1)*10.5;txt(cls[j],xx,14,6.5,TRUE,pal[j],just='center');txt(paste0('n=',c(20,41,32)[j]),xx,18,5.7,just='center')
  for(i in 1:2){yy<-27+(i-1)*13;r<-cp[FinalSubtype==cls[j]&Duct==c('Small','Large')[i]];hi<-i==1&&j==3
   rect(xx,yy,9.8,11.5,adjustcolor(pal[j],alpha.f=if(hi).23 else .09),if(hi)pal[j] else NA,.8)
   txt(r$Count,xx,yy-1.6,8,TRUE,just='center');txt(sprintf('%.1f%%',100*r$Fraction),xx,yy+2.5,5.8,just='center')
  }
 };txt('Small\nduct',2,27,6.2);txt('Large\nduct',2,40,6.2)
 txt('15 of 32 C3 tumors',2,50.5,6.5,TRUE,pal['C3']);txt('have small-duct morphology',2,55,6.0)
}
save_panel('Figure6C_observed_duct_overlap',drawC,48,58,cp);check('C measured-only 93 patients and C3 small15',sum(cp$Count)==93&cp[FinalSubtype=='C3'&Duct=='Small',Count]==15)

# F: all 4473 patient-drug values and aligned distribution summaries.
z<-ro('Figure6F_drug_heterogeneity');cand<-rd('8-Visualization/refinements/fig6_translational_20261101_v17/tables/nine_candidate_inventory.csv');setorder(cand,CandidateSubtype,Drug)
z[,Drug:=factor(Drug,levels=rev(cand$Drug))];z[,CandidateGroup:=factor(paste(CandidateSubtype,'candidates'),levels=paste(cls,'candidates'))]
short<-setNames(sub('_[0-9]+$','',cand$Drug),cand$Drug);short['Sepantronium bromide_1941']<-'Sepantronium bromide'
summ<-z[,.(N=.N,Median=median(Z),Q1=quantile(Z,.25),Q3=quantile(Z,.75)),by=.(Drug,CandidateGroup,FinalSubtype)]
hm<-ggplot(z,aes(PatientOrder,Drug,fill=pmax(-2,pmin(2,Z))))+geom_tile(width=1,height=1)+facet_grid(CandidateGroup~FinalSubtype,scales='free',space='free',labeller=labeller(FinalSubtype=c(C1='C1 · 112',C2='C2 · 223',C3='C3 · 162')))+scale_y_discrete(labels=short,expand=expansion(add=.5))+scale_x_continuous(expand=expansion(add=.5))+scale_fill_gradient2(low='#397CA6',mid='white',high='#B66D7D',limits=c(-2,2),breaks=c(-2,0,2),labels=c('≤−2','0','≥2'),name='Score (SD)',guide=guide_colorbar(title.position='top',barwidth=unit(16,'mm'),barheight=unit(1.5,'mm')))+labs(x=NULL,y=NULL,title='Individual tumors')+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=5.7),strip.text.y=element_blank(),plot.title=element_text(size=6.4),panel.spacing=unit(.7,'mm'))
q<-ggplot(summ,aes(Median,Drug,color=FinalSubtype))+geom_vline(xintercept=0,linetype=2,color='#C0C8CD',linewidth=.25)+geom_errorbar(aes(xmin=Q1,xmax=Q3),orientation='y',position=position_dodge(width=.65),width=0,linewidth=.45)+geom_point(position=position_dodge(width=.65),size=.9)+facet_grid(CandidateGroup~.,scales='free_y',space='free_y')+scale_y_discrete(expand=expansion(add=.5))+scale_x_continuous(limits=c(-2.3,2.3),breaks=c(-2,0,2))+scale_color_manual(values=pal)+labs(x='Score (SD)',y=NULL,title='Median / IQR')+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),strip.text.y=element_blank(),plot.title=element_text(size=6.4),legend.key.width=unit(2,'mm'),legend.spacing.x=unit(0,'mm'),legend.margin=margin(0,0,0,0),panel.spacing=unit(.7,'mm'))
p<-hm+q+plot_layout(widths=c(3.6,1.2))+plot_annotation(title='Predicted drug patterns',subtitle='Nine candidates · GDSC2-trained oncoPredict',theme=theme(plot.title=element_text(size=7.2,face='bold',family='Arial'),plot.subtitle=element_text(size=6,family='Arial'),plot.margin=margin(1.3,1.3,1.3,1.3,'mm')))
save_panel('Figure6F_drug_heterogeneity',p,104,72,z,more=list(distribution_summary=summ));check('F all497 patients nine drugs retained',nrow(z)==4473&nrow(summ)==27)

# G: grouped compound labels save space while preserving all90 assay cells.
g<-ro('Figure6G_molecular_evidence');rows<-unique(g[,.(Drug,Gene,CandidateSubtype)]);rows[,DrugOrder:=match(Drug,cand$Drug)];setorder(rows,DrugOrder,Gene)
rows[,Row:=seq_len(.N)];rows[,Y:=18+(Row-1)*2.7+.65*(DrugOrder-1)];g<-merge(g,rows[,.(Drug,Gene,Y)],by=c('Drug','Gene'))
drawG<-function(){start(66,72);on.exit(popViewport());txt('Molecular context',1,3.4,7.2,TRUE);txt('Nine candidates · 15 associated genes',1,7.5,5.8)
 txt('RNA',39.2,11.5,6.1,TRUE,just='center');txt('Protein',56.8,11.5,6.1,TRUE,just='center')
 for(j in 1:3){txt(cls[j],34.8+(j-1)*4.4,15,5.8,TRUE,pal[j],just='center');txt(cls[j],52.4+(j-1)*4.4,15,5.8,TRUE,pal[j],just='center')}
 for(dr in cand$Drug){r<-rows[Drug==dr];lab<-sub('_[0-9]+$','',dr);lab<-sub('Sepantronium bromide','Sepantronium',lab,fixed=TRUE);col<-pal[r$CandidateSubtype[1]]
  txt(lab,1,mean(r$Y),5.5,FALSE,col);ln(c(.6,.6),range(r$Y)+c(-1.1,1.1),col,1.2)
 }
 for(i in 1:nrow(rows)){r<-rows[i];txt(r$Gene,25.5,r$Y,5.4,FALSE,ink,just='center')
  for(ass in c('RNA','Protein'))for(j in 1:3){d<-g[Drug==r$Drug&Gene==r$Gene&Subtype==cls[j]&Assay==ass];xx<-if(ass=='RNA')34.8+(j-1)*4.4 else 52.4+(j-1)*4.4
   fill<-c(Higher='#BD7588',Lower='#6596AE','Not significant'='white',Unavailable='#ECEFF1')[d$State];rect(xx,r$Y,4,2.4,fill,'#D6DDE1',.3)
   txt(c(Higher='+',Lower='−','Not significant'='·',Unavailable='×')[d$State],xx,r$Y,5.8,FALSE,if(d$State%in%c('Higher','Lower'))'white' else muted,just='center')
  }
 }
 for(i in 1:4){xx<-c(2,17,31,45)[i];rect(xx,66,2,2,c('#BD7588','#6596AE','white','#ECEFF1')[i],'#D6DDE1',.3);txt(c('Higher','Lower','NS','Missing')[i],xx+1.7,66,5.4)}
 txt('Direction at assay-specific FDR < 0.05',1,70,5.5)
}
save_panel('Figure6G_molecular_evidence',drawG,66,72,g);check('G complete90 assay cells15 genes9drugs',nrow(g)==90&uniqueN(g$Gene)==15&uniqueN(g$Drug)==9)

# H: preserve all patients and both genes; no causal or efficacy inference.
hz<-z[as.character(Drug)=='Ulixertinib_1908'];h<-copy(g[Drug=='Ulixertinib_1908']);h[,Y:=ifelse(Gene=='MAPK1',2,1)+c(C1=.21,C2=0,C3=-.21)[Subtype]];h[,Significant:=FDR<.05];h[,Assay:=factor(Assay,levels=c('RNA','Protein'))]
x<-ggplot(hz,aes(FinalSubtype,Z,color=FinalSubtype))+geom_point(position=position_jitter(width=.16,height=0,seed=20261101),size=.32,alpha=.3)+geom_boxplot(width=.45,outlier.shape=NA,fill=NA,linewidth=.35)+scale_color_manual(values=pal,guide='none')+scale_x_discrete(labels=c(C1='C1 · 112',C2='C2 · 223',C3='C3 · 162'))+labs(title='Ulixertinib prediction',x=NULL,y='Relative score (SD)')+theme(plot.title=element_text(size=6.5))
y<-ggplot(h,aes(Effect,Y,color=Subtype,shape=Significant))+geom_vline(xintercept=0,linetype=2,linewidth=.25,color='#ACB8C0')+geom_point(size=1.4)+facet_wrap(~Assay,nrow=1,scales='free_x',labeller=labeller(Assay=c(RNA='RNA · log2 effect',Protein='Protein · SD effect')))+scale_color_manual(values=pal)+scale_shape_manual(values=c('TRUE'=16,'FALSE'=1),labels=c('FALSE'='FDR ≥ 0.05','TRUE'='FDR < 0.05'))+scale_y_continuous(breaks=c(2,1),labels=c('MAPK1','MAPK3'),limits=c(.5,2.5))+scale_x_continuous(expand=expansion(mult=.15))+labs(title='Both mapped ERK genes',x='Study-adjusted expression difference',y=NULL)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.title=element_text(size=6.5),legend.key.width=unit(2.5,'mm'),legend.key.height=unit(2,'mm'),legend.spacing.x=unit(.5,'mm'),legend.margin=margin(0,0,0,0),legend.box='horizontal')
p<-x+y+plot_layout(widths=c(1,1.8))+plot_annotation(title='C3-focused candidate for experimental follow-up',theme=theme(plot.title=element_text(size=7.2,face='bold',family='Arial'),plot.margin=margin(1.2,1.2,1.2,1.2,'mm')))
save_panel('Figure6H_C3_hypothesis',p,174,44,h,more=list(patient_predictions=hz));check('H all497patients and12expressioneffects',nrow(hz)==497&nrow(h)==12)

# D/E: all four prespecified workflows; paired changes reveal the small joint increment.
perf<-rd(file.path(joint,'tables/performance_ci.csv'));diffs<-rd(file.path(joint,'tables/paired_differences.csv'));trans<-rd(file.path(joint,'tables/C3_binary_transition_counts.csv'))
models<-c('Anchor3','CompactWorkflow','JointWorkflow','FullWorkflow');mlabel<-c('RNA anchors','Compact','Joint','Full workflow')
roles<-setNames(c('#85919A','#426C85','#243E50','#6F8C86'),models)
metrics<-c('BalancedAccuracy','C3Sensitivity','C3Specificity');dp<-perf[Model%chin%models&Metric%chin%metrics]
drawD<-function(){start(68,58);on.exit(popViewport());txt('Molecular identification',2,3.4,7.2,TRUE);txt('497 tumors · OOF · 95% intervals',2,8,5.9)
 centers<-c(29,44,59);for(j in 1:3)txt(c('Balanced\naccuracy','C3\nsensitivity','C3\nspecificity')[j],centers[j],14,5.7,TRUE,just='center')
 for(i in 1:4){y<-22+(i-1)*8;txt(mlabel[i],2,y,6.0,bold=i==3)
  for(j in 1:3){r<-dp[Model==models[i]&Metric==metrics[j]];xx<-function(v)centers[j]-5+(100*v-50)/50*10
   ln(xx(c(.5,1)),rep(y+1.3,2),light,.35);interval(r,xx,y+1.3,roles[i]);txt(sprintf('%.1f',100*r$Estimate),centers[j],y-1.8,6.2,bold=i==3,just='center')}
 }
 for(j in 1:3){txt('50',centers[j]-5,51,5.2,just='center');txt('100%',centers[j]+5,51,5.2,just='center')}
 txt('Joint = RNA anchors + adaptive compact',2,55.5,5.7)
}
save_panel('Figure6D_identification_comparison',drawD,68,58,dp,more=list(paired_differences=diffs));check('D twelve estimates fourmodels',nrow(dp)==12)
types<-c('Both correct','Corrected','New error','Both incorrect');tcols<-c('#EDF1F4','#D5E8EF','#F3DDD0','#DDE2E6')
ep<-trans[Comparison%chin%c('CompactWorkflow-Anchor3','JointWorkflow-CompactWorkflow')]
drawE<-function(){start(50,58);on.exit(popViewport());txt('Paired C3 reclassification',1,3.4,7.0,TRUE)
 for(k in 1:2){basey<-if(k==1)0 else 24;cmp<-c('CompactWorkflow-Anchor3','JointWorkflow-CompactWorkflow')[k]
  txt(c('RNA anchors → compact','Compact → joint')[k],1,8.5+basey,6.2,TRUE)
  for(j in 1:4){xx<-17.5+(j-1)*8.8;txt(c('Correct\nboth','Errors\nfixed','New\nerrors','Wrong\nboth')[j],xx,14+basey,5.3,just='center')
   for(i in 1:2){yy<-21+basey+(i-1)*6.1;r<-ep[Comparison==cmp&ReferenceGroup==c('C3','Non-C3')[i]&Transition==types[j]];stopifnot(nrow(r)==1);rect(xx,yy,8.2,5.4,tcols[j]);txt(r$N,xx,yy,6.8,TRUE,just='center')}
  };txt('C3',1,21+basey,5.8,TRUE,pal['C3']);txt('Non-C3',1,27.1+basey,5.8)
 };txt('Reference groups: 162 C3 / 335 non-C3',1,56,5.4)
}
save_panel('Figure6E_paired_reclassification',drawE,50,58,ep)
check('E each comparison retains497',all(ep[,sum(N),by=Comparison]$V1==497))
check('Joint adds one correctedC3 and no newerrors',ep[Comparison=='JointWorkflow-CompactWorkflow'&Transition=='Corrected',sum(N)]==1&ep[Comparison=='JointWorkflow-CompactWorkflow'&Transition=='New error',sum(N)]==0)

# Preserve every existing supplementary panel and its tables, with explicit provenance.
oldmanifest<-rd(file.path(old,'figure_manifest.csv'))
for(i in which(oldmanifest$Track=='Supplement')){
 r<-oldmanifest[i];id<-r$Panel;pdf<-file.path(out,'supplement',paste0(id,'.pdf'));png<-file.path(out,'review',paste0(id,'.png'))
 stopifnot(file.copy(r$PDF,pdf,overwrite=TRUE),file.copy(r$PNG,png,overwrite=TRUE))
 tt<-file.path(old,'tables',paste0(id,'.csv'));dest<-file.path(out,'tables',basename(tt));stopifnot(file.copy(tt,dest,overwrite=TRUE))
 inputs<-unique(c(inputs,r$PDF,r$PNG,tt));manifest[[id]]<-data.table(Panel=id,Track='Supplement',PDF=pdf,PNG=png,SourceTable=dest,WidthMM=r$WidthMM,HeightMM=r$HeightMM,Rows=nrow(fread(dest)),Reuse='Byte-identical existing supplement',OriginPDF=r$PDF)
}

# M: complete model scope (all subtype recalls and C3 binary metrics).
s<-perf[Model%chin%models&Metric%chin%c('BalancedAccuracy','C1Recall','C2Recall','C3Recall','C3Specificity','C3PPV')]
s[,Model:=factor(Model,levels=rev(models))];s[,Metric:=factor(Metric,levels=c('BalancedAccuracy','C1Recall','C2Recall','C3Recall','C3Specificity','C3PPV'),labels=c('Balanced accuracy','C1 recall','C2 recall','C3 sensitivity','C3 specificity','C3 PPV'))]
p<-ggplot(s,aes(100*Estimate,Model,color=Model))+geom_vline(xintercept=c(0,100),color=light,linewidth=.25)+geom_errorbar(aes(xmin=100*Lower95,xmax=100*Upper95),orientation='y',width=.10,linewidth=.35)+geom_point(size=1.4)+geom_text(aes(label=sprintf('%.1f',100*Estimate)),nudge_y=.24,size=6.2/.pt,show.legend=FALSE)+facet_wrap(~Metric,ncol=3)+scale_color_manual(values=roles,guide='none')+scale_y_discrete(labels=setNames(mlabel,models))+scale_x_continuous(limits=c(-5,108),breaks=c(0,50,100))+labs(title='Complete identification performance',subtitle='Same 497 tumors · identical outer partitions · conditional 95% intervals',x='OOF performance (%)',y=NULL)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),panel.spacing=unit(5,'mm'),strip.text=element_text(size=7),plot.title=element_text(size=8),plot.subtitle=element_text(size=7))
save_panel('SFig6M_joint_identification_performance',p,174,116,s,'Supplement')

# N: all four complete confusion matrices; every count including zero is shown.
cm<-rd(file.path(joint,'tables/confusion_matrices.csv'));cm[,Truth:=factor(Truth,levels=rev(cls))];cm[,Model:=factor(Model,levels=models,labels=mlabel)]
p<-ggplot(cm,aes(Prediction,Truth,fill=RowFraction))+geom_tile(color='white',linewidth=.65)+geom_text(aes(label=sprintf('%d\n%.1f%%',Count,100*RowFraction),color=RowFraction>.65),size=7/.pt,lineheight=1.05)+scale_color_manual(values=c('FALSE'=ink,'TRUE'='white'),guide='none')+scale_fill_gradient(low='#F3F6F7',high='#305D75',limits=c(0,1),guide='none')+facet_wrap(~Model,ncol=2)+coord_fixed()+scale_y_discrete(labels=c(C1='C1 ·112',C2='C2 ·223',C3='C3 ·162'))+labs(title='Complete three-class error structure',subtitle='Reference molecular labels versus out-of-fold predictions',x='Predicted subtype',y='Reference subtype · n')+theme(axis.line=element_blank(),axis.ticks=element_blank(),strip.text=element_text(size=8),plot.title=element_text(size=9),plot.subtitle=element_text(size=7),panel.spacing=unit(5,'mm'))
save_panel('SFig6N_joint_confusion_matrices',p,174,158,cm,'Supplement');check('N all four497 predictions',sum(cm$Count)==4*497&nrow(cm)==36)

# O: every paired C3 transition, with zeros explicit and no pooling of reference groups.
tt<-copy(trans);tt[,ReferenceGroup:=factor(ReferenceGroup,levels=c('Non-C3','C3'))];tt[,Transition:=factor(Transition,levels=types,labels=c('Correct in both','Errors corrected','New errors','Wrong in both'))];tt[,Comparison:=factor(Comparison,levels=c('CompactWorkflow-Anchor3','JointWorkflow-Anchor3','JointWorkflow-CompactWorkflow'),labels=c('RNA anchors → compact','RNA anchors → joint','Compact → joint'))];tt[,Denominator:=ifelse(ReferenceGroup=='C3',162L,335L)];tt[,Pct:=100*N/Denominator]
p<-ggplot(tt,aes(Transition,ReferenceGroup,fill=Transition))+geom_tile(width=.94,height=.85)+geom_text(aes(label=sprintf('%d\n%.1f%%',N,Pct)),size=7/.pt,lineheight=1.05,color=ink)+facet_wrap(~Comparison,ncol=1)+scale_fill_manual(values=setNames(tcols,levels(tt$Transition)),guide='none')+scale_y_discrete(labels=c(C3='C3 ·162','Non-C3'='Non-C3 ·335'))+scale_x_discrete(position='top')+labs(title='Paired patient reclassification for C3 detection',subtitle='Every patient contributes once per comparison; row percentages',x=NULL,y=NULL)+theme(axis.line=element_blank(),axis.ticks=element_blank(),strip.text=element_text(size=8),plot.title=element_text(size=9),plot.subtitle=element_text(size=7),axis.text=element_text(size=7),panel.spacing=unit(6,'mm'))
save_panel('SFig6O_all_C3_transitions',p,174,151,tt,'Supplement');check('O three complete paired comparisons',nrow(tt)==24&sum(tt$N)==3*497)

# P: patient-level correspondence explains why the pathology and OOF marginals differ.
duct<-ro('Figure6C_duct_proportions_2');meas<-duct[Measured==TRUE];stopifnot(nrow(meas)==93,uniqueN(meas$Sample)==93)
ac<-meas[,.(N=.N),by=.(FinalSubtype,DuctObserved,Class)];ac<-merge(CJ(FinalSubtype=cls,DuctObserved=c('small duct type','large duct type'),Class=c('small duct type','large duct type')),ac,all.x=TRUE);ac[is.na(N),N:=0L];ac[,ObservedN:=sum(N),by=.(FinalSubtype,DuctObserved)];ac[,Fraction:=N/ObservedN];ac[,Observed:=factor(DuctObserved,levels=c('large duct type','small duct type'),labels=c('Large','Small'))];ac[,Predicted:=factor(Class,levels=c('small duct type','large duct type'),labels=c('Small','Large'))]
p<-ggplot(ac,aes(Predicted,Observed,fill=Fraction))+geom_tile(color='white',linewidth=.65)+geom_text(aes(label=sprintf('%d / %d\n%.1f%%',N,ObservedN,100*Fraction),color=Fraction>.65),size=7/.pt,lineheight=1.1)+scale_color_manual(values=c('FALSE'=ink,'TRUE'='white'),guide='none')+scale_fill_gradient(low='#F3F6F7',high='#596D85',limits=c(0,1),guide='none')+facet_wrap(~FinalSubtype,nrow=1)+coord_fixed()+labs(title='Observed pathology versus paired RNA duct prediction',subtitle='The same 93 patients · counts and within-observed-category percentages',x='OOF RNA duct prediction',y='Measured duct morphology')+theme(axis.line=element_blank(),axis.ticks=element_blank(),strip.text=element_text(size=8),plot.title=element_text(size=9),plot.subtitle=element_text(size=7),axis.text=element_text(size=7),axis.title=element_text(size=7))
save_panel('SFig6P_observed_OOF_duct_agreement',p,174,82,ac,'Supplement',more=list(patient_pairs=meas));check('P patient-level93 observed no projection',sum(ac$N)==93)

# Q: original observed/OOF/projected marginals remain available as a separate supplement.
id<-'SFig6Q_duct_source_proportions';src<-'Figure6C_duct_proportions';pdf<-file.path(out,'supplement',paste0(id,'.pdf'));png<-file.path(out,'review',paste0(id,'.png'));tb<-file.path(out,'tables',paste0(id,'.csv'))
stopifnot(file.copy(file.path(old,'figures',paste0(src,'.pdf')),pdf,overwrite=TRUE),file.copy(file.path(old,'review',paste0(src,'.png')),png,overwrite=TRUE));fwrite(csrc,tb)
manifest[[id]]<-data.table(Panel=id,Track='Supplement',PDF=pdf,PNG=png,SourceTable=tb,WidthMM=144,HeightMM=95,Rows=nrow(csrc),Reuse='Former main C; byte-identical PDF')
inputs<-unique(c(inputs,file.path(old,'figures',paste0(src,'.pdf')),file.path(old,'review',paste0(src,'.png'))))

# R/S: independently computed exploratory duct-stratified survival module.
clin<-'6-ClinicalSurvival/review/20261101_duct_c3_stratified_v1';clman<-rd(file.path(clin,'figure_manifest.csv'))
for(pair in list(c('SFig6R_duct_C3_four_group_survival','Duct_C3_four_group_KM'),c('SFig6S_small_duct_C3_association','Small_duct_C3_association'))){
 id<-pair[1];src<-pair[2];pdf<-file.path(out,'supplement',paste0(id,'.pdf'));png<-file.path(out,'review',paste0(id,'.png'));spdf<-file.path(clin,'figures',paste0(src,'.pdf'));spng<-file.path(clin,'review',paste0(src,'.png'))
 stopifnot(file.copy(spdf,pdf,overwrite=TRUE),file.copy(spng,png,overwrite=TRUE))
 datafile<-file.path(clin,'tables',if(grepl('four_group',id))'km_curves.csv' else 'small_duct_c3_effects.csv');d<-rd(datafile);tb<-file.path(out,'tables',paste0(id,'.csv'));fwrite(d,tb)
 manifest[[id]]<-data.table(Panel=id,Track='Supplement',PDF=pdf,PNG=png,SourceTable=tb,WidthMM=174,HeightMM=if(grepl('four_group',id))113 else 112,Rows=nrow(d),Reuse='New exploratory analysis; imported checked R export',OriginPDF=spdf)
 inputs<-unique(c(inputs,spdf,spng));
}

# Native-size four-row preview for review only; independent PDFs are the masters.
positions<-data.table(Panel=c('Figure6A_clinical_information','Figure6B_measured_duct_increment','Figure6C_observed_duct_overlap','Figure6D_identification_comparison','Figure6E_paired_reclassification','Figure6F_drug_heterogeneity','Figure6G_molecular_evidence','Figure6H_C3_hypothesis'),Letter=LETTERS[1:8],X=c(4,93,4,56,128,4,112,4),Top=c(4,4,84,84,84,146,146,222),Width=c(85,85,48,68,50,104,66,174),Height=c(76,76,58,58,58,72,72,44))
ragg::agg_png(file.path(out,'review/Figure6_layout_preview.png'),width=182,height=270,units='mm',res=300,background='white');start(182,270)
for(i in 1:nrow(positions)){r<-positions[i];im<-png::readPNG(file.path(out,'review',paste0(r$Panel,'.png')));grid.raster(im,x=unit(r$X,'native'),y=unit(r$Top,'native'),width=unit(r$Width,'mm'),height=unit(r$Height,'mm'),just=c('left','top'));txt(r$Letter,r$X-3,r$Top+2,7.5,TRUE)}
popViewport();dev.off();fwrite(positions,file.path(out,'tables/layout_positions_mm.csv'));check('Layout all panel edges within182x270',all(positions$X+positions$Width<=178)&all(positions$Top+positions$Height<=266))

# Complete artifact index and machine-readable audit.
mf<-rbindlist(manifest,fill=TRUE);setorder(mf,Track,Panel);fwrite(mf,file.path(out,'figure_manifest.csv'))
fwrite(rbindlist(manifest,fill=TRUE),file.path(out,'partial_manifest.csv'))
fwrite(rbindlist(text_bounds),file.path(out,'qa/text_bounds.csv'))
fwrite(rbindlist(checks),file.path(out,'qa/numeric_checks.csv'))
fwrite(data.table(Path=inputs,SHA256=vapply(inputs,function(p)digest(file=p,algo='sha256'),character(1))),file.path(out,'tables/input_hashes.csv'))
saveRDS(list(manifest=manifest,checks=checks,inputs=inputs),file.path(out,'provenance/render_state.rds'))
writeLines(capture.output(sessionInfo()),file.path(out,'qa/sessionInfo.txt'))
