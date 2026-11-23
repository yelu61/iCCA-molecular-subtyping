suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(digest)})
source('8-Visualization/scripts/fig_style.R')
out<-'8-Visualization/refinements/fig6_redesign_20261101_v18';old<-'8-Visualization/refinements/fig6_translational_20261101_v17'
for(d in c('figures','supplement','review','tables','qa'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
inputs<-character();rd<-function(p){inputs<<-unique(c(inputs,p));cols<-names(fread(p,nrows=0));if('Sample'%in%cols)fread(p,colClasses=c(Sample='character'))else fread(p)}
br<-'5-ClassifierValidation/review/20261101_duct_compact_bridge_v1/tables';readbr<-function(n)rd(file.path(br,paste0(n,'.csv')))
pal<-load_subtype_palette()$colors[c('C1','C2','C3')];invisible(load_study_palette());cls<-c('C1','C2','C3')
role<-c(Anchor3='#808A91',CompactWorkflow='#305D75',Anchor4='#AE8A58',FullWorkflow='#527F79')
rolelab<-c(Anchor3='Duct RNA',CompactWorkflow='Compact',Anchor4='Duct + MUC6',FullWorkflow='Full')
base<-theme_classic(base_size=8,base_family='Arial')+theme(axis.text=element_text(size=7.5,color='#27333A'),axis.title=element_text(size=8),plot.title=element_text(size=9,face='bold'),plot.subtitle=element_text(size=7.5),plot.caption=element_text(size=7,hjust=0),plot.title.position='plot',plot.caption.position='plot',legend.position='bottom',legend.title=element_blank(),legend.text=element_text(size=7.5),strip.background=element_blank(),strip.text=element_text(size=8,face='bold'),panel.spacing=grid::unit(2,'mm'),plot.margin=margin(2,2,2,2,'mm'))
theme_set(base);manifest<-list();checks<-list()
save<-function(id,p,w,h,d,track='Main'){
 pdf<-file.path(out,if(track=='Main')'figures' else 'supplement',paste0(id,'.pdf'));png<-file.path(out,'review',paste0(id,'.png'));tb<-file.path(out,'tables',paste0(id,'.csv'))
 cairo_pdf(pdf,width=w/25.4,height=h/25.4,family='Arial');print(p);dev.off();ragg::agg_png(png,width=w,height=h,units='mm',res=450,background='white');print(p);dev.off();fwrite(d,tb)
 manifest[[id]]<<-data.table(Panel=id,Track=track,PDF=pdf,PNG=png,SourceTable=tb,WidthMM=w,HeightMM=h,Reuse='Redesigned',Rows=nrow(d));message(id)
}
check<-function(n,b){stopifnot(all(b));checks[[n]]<<-data.table(Check=n,Pass=TRUE)}

# D: absolute estimates, uncertainty and paired gain without duplicated rows.
perf<-readbr('performance_ci');met<-c('BalancedAccuracy','C1Recall','C2Recall','C3Recall');mlab<-c('Balanced\naccuracy','C1 recall','C2 recall','C3 recall')
d<-perf[Model%chin%c('Anchor3','CompactWorkflow')&Metric%chin%met];d[,Y:=5-match(Metric,met)]
wide<-dcast(d,Metric+Y~Model,value.var='Estimate');gain<-readbr('paired_differences')[Comparison=='CompactWorkflow-Anchor3'&Metric=='BalancedAccuracy']
p<-ggplot(d,aes(Estimate*100,Y,color=Model))+geom_vline(xintercept=c(0,50,100),color='#E6E9EB',linewidth=.3)+
 geom_segment(data=wide,aes(x=100*Anchor3,xend=100*CompactWorkflow,y=Y,yend=Y),inherit.aes=FALSE,color='#CCD2D6',linewidth=.7)+
 geom_errorbar(aes(xmin=100*Lower95,xmax=100*Upper95),orientation='y',width=.10,linewidth=.45)+geom_point(size=1.8)+
 geom_text(aes(y=Y+ifelse(Model=='CompactWorkflow',.23,-.23),label=sprintf('%.1f',100*Estimate)),size=7.5/.pt,show.legend=FALSE)+
 scale_color_manual(values=role,breaks=c('Anchor3','CompactWorkflow'),labels=rolelab[c('Anchor3','CompactWorkflow')])+
 scale_x_continuous(limits=c(-3,106),breaks=c(0,50,100),expand=expansion(mult=0))+scale_y_continuous(breaks=4:1,labels=mlab,limits=c(.5,4.5))+
 labs(title='Identification beyond duct RNA',subtitle=sprintf('BA gain +%.1f pp (95%% CI %.1f–%.1f)',gain$Estimate*100,gain$Lower95*100,gain$Upper95*100),x='OOF performance (%)',y=NULL)+
 theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),legend.margin=margin(0,0,0,0),legend.key.height=grid::unit(3,'mm'))
save('Figure6D_identification_gain',p,90,82,d)
check('D preserves all eight estimates and intervals',nrow(d)==8&all(d$Lower95<=d$Estimate)&all(d$Upper95>=d$Estimate))

# E: explicit correctness transitions; zero cells remain visible.
tr<-readbr('patient_transition_counts');types<-c('Both correct','Corrected classification','Lost correct classification','Both incorrect')
e<-merge(CJ(Truth=cls,Transition=types),tr,by=c('Truth','Transition'),all.x=TRUE);e[is.na(N),N:=0L];e[,RowN:=sum(N),by=Truth];e[,Pct:=100*N/RowN]
e[,Category:=factor(Transition,levels=types,labels=c('Correct\nin both','Errors\ncorrected','New\nerrors','Persistent\nerrors'))];e[,Truth:=factor(Truth,levels=rev(cls))]
epal<-setNames(c('#E9EDF0','#D8E9EF','#F4E3D8','#D1D8DE'),levels(e$Category))
p<-ggplot(e,aes(Category,Truth,fill=Category))+geom_tile(width=.94,height=.88)+
 geom_text(aes(label=N),size=10/.pt,fontface='bold',nudge_y=.10,color='#263E4B')+
 geom_text(aes(label=sprintf('%.1f%%',Pct)),size=7.5/.pt,nudge_y=-.20,color='#52616A')+
 scale_fill_manual(values=epal,guide='none')+scale_x_discrete(position='top',expand=expansion(add=.5))+
 scale_y_discrete(labels=c(C1='C1\nn=112',C2='C2\nn=223',C3='C3\nn=162'),expand=expansion(add=.6))+
 labs(title='Patient-level reclassification',subtitle='143 corrected  |  3 new errors',x=NULL,y=NULL,caption='Duct RNA → compact workflow; row percentages')+
 theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(size=7.5),panel.grid=element_blank())
save('Figure6E_patient_corrections',p,90,82,e)
check('E complete paired transitions',sum(e$N)==497&sum(e[Transition=='Corrected classification',N])==143&sum(e[Transition=='Lost correct classification',N])==3&nrow(e)==12)

# F: every patient plus directly aligned subtype distribution summaries.
z<-rd(file.path(old,'tables/Figure6F_patient_drug_patterns.csv'));cand<-rd(file.path(old,'tables/nine_candidate_inventory.csv'));setorder(cand,CandidateSubtype,Drug)
z[,Drug:=factor(Drug,levels=rev(cand$Drug))];z[,CandidateGroup:=factor(paste(CandidateSubtype,'candidates'),levels=paste(cls,'candidates'))]
short<-setNames(sub('_[0-9]+$','',cand$Drug),cand$Drug);clab<-c(C1='C1 · 112',C2='C2 · 223',C3='C3 · 162')
summ<-z[,.(N=.N,Median=median(Z),Q1=quantile(Z,.25),Q3=quantile(Z,.75)),by=.(Drug,CandidateGroup,FinalSubtype)]
fwrite(summ,file.path(out,'tables/Figure6F_distribution_summary.csv'))
hm<-ggplot(z,aes(PatientOrder,Drug,fill=pmax(-2,pmin(2,Z))))+geom_tile(width=1,height=1)+facet_grid(CandidateGroup~FinalSubtype,scales='free',space='free',labeller=labeller(FinalSubtype=clab))+
 scale_y_discrete(labels=short,expand=expansion(add=.5))+scale_x_continuous(expand=expansion(add=.5))+
 scale_fill_gradient2(low='#397CA6',mid='white',high='#B66D7D',limits=c(-2,2),breaks=c(-2,0,2),labels=c('≤−2','0','≥2'),name='Score (SD)',guide=guide_colorbar(title.position='top',barwidth=grid::unit(20,'mm'),barheight=grid::unit(2,'mm')))+
 labs(x=NULL,y=NULL,title='Individual tumors')+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_blank(),strip.text.y=element_blank(),plot.title=element_text(size=8),panel.spacing=grid::unit(1.3,'mm'))
q<-ggplot(summ,aes(Median,Drug,color=FinalSubtype))+geom_vline(xintercept=0,linetype=2,color='#C0C8CD',linewidth=.3)+
 geom_errorbar(aes(xmin=Q1,xmax=Q3),orientation='y',position=position_dodge(width=.65),width=0,linewidth=.6)+geom_point(position=position_dodge(width=.65),size=1.3)+
 facet_grid(CandidateGroup~.,scales='free_y',space='free_y')+scale_y_discrete(expand=expansion(add=.5))+
 scale_x_continuous(limits=c(-2.3,2.3),breaks=c(-2,0,2))+scale_color_manual(values=pal)+
 labs(x='Score (SD)',y=NULL,title='Median / IQR')+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),strip.text.y=element_text(size=7),plot.title=element_text(size=8),legend.key.width=grid::unit(3,'mm'),legend.spacing.x=grid::unit(0,'mm'),panel.spacing=grid::unit(1.3,'mm'))
p<-hm+q+plot_layout(widths=c(3.2,1))+plot_annotation(title='Candidate drug patterns and patient heterogeneity',subtitle='GDSC2-trained predictions · 497 tumors',theme=theme(plot.title=element_text(size=9,face='bold',family='Arial'),plot.subtitle=element_text(size=7.5,family='Arial'),plot.margin=margin(2,2,2,2,'mm')))
save('Figure6F_drug_heterogeneity',p,174,101,z)
check('F all values and subtype distribution cells retained',nrow(z)==4473&nrow(summ)==27&all(summ$Q1<=summ$Median)&all(summ$Q3>=summ$Median)&sum(summ$N)==4473)

# G: paired effects on nine drug rows; separate heterogeneity and selection stability.
g<-rd(file.path(old,'tables/Figure6G_drug_candidate_uncertainty.csv'))
ret<-rd('5-ClassifierValidation/review/20261101_drug_evidence_refresh_v1/tables/original_candidate_retention.csv')
pos<-data.table(Drug=cand$Drug,Y=c(10,9,8,6.5,5.5,4.5,3,2,1));g<-merge(g,pos,by='Drug');g[,YPoint:=Y+ifelse(Other==sort(Other)[1],.16,-.16),by=Drug]
rank<-merge(g[,.(Y=Y[1],Subtype=Subtype[1],MinI2=min(I2_percent),MaxI2=max(I2_percent)),by=Drug],ret,by=c('Drug','Subtype'))
rank[,I2Label:=sprintf('%.0f–%.0f%%',MinI2,MaxI2)];rank[,RankLabel:=paste0(TopThreeRetained,' / ',Omissions)]
fwrite(rank,file.path(out,'tables/Figure6G_stability_summary.csv'))
forest<-ggplot(g,aes(Effect,YPoint,color=Subtype,shape=Other))+geom_vline(xintercept=0,linetype=2,color='#B4BDC3',linewidth=.4)+
 geom_errorbar(aes(xmin=Lower,xmax=Upper),orientation='y',width=.08,linewidth=.45)+geom_point(size=1.7)+
 scale_color_manual(values=pal,name='Candidate')+scale_shape_manual(values=c(C1=16,C2=17,C3=15),labels=c(C1='vs C1',C2='vs C2',C3='vs C3'),name='Comparator')+
 guides(color=guide_legend(order=1,override.aes=list(shape=16)),shape=guide_legend(order=2,override.aes=list(color='#26333B')))+
 scale_y_continuous(breaks=pos$Y,labels=short[pos$Drug],limits=c(.5,10.6),expand=expansion(mult=0))+
 scale_x_continuous(limits=c(min(g$Lower)-.1,.22),breaks=c(-2,-1,0))+
 labs(title='Pairwise effects',x='Predicted score difference (SD; 95% CI)',y=NULL)+theme(plot.title=element_text(size=8),axis.line.y=element_blank(),axis.ticks.y=element_blank())
stab<-ggplot(rank,aes(y=Y))+geom_text(aes(x=1,label=I2Label),size=7.5/.pt,color='#40515C')+geom_text(aes(x=2,label=RankLabel,color=Subtype),size=8/.pt,fontface='bold')+
 scale_y_continuous(limits=c(.5,10.6),expand=expansion(mult=0))+scale_x_continuous(limits=c(.5,2.5),breaks=c(1,2),labels=c('I² range','Top-3\nretained'),position='top')+scale_color_manual(values=pal,guide='none')+
 labs(title='Across studies',x=NULL,y=NULL)+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_blank(),plot.title=element_text(size=8),axis.text.x=element_text(size=7.5),legend.position='none')
p<-forest+stab+plot_layout(widths=c(3.4,1.25))+plot_annotation(title='Candidate effects, uncertainty and selection stability',subtitle='Nine candidates · 18 contrasts · seven study omissions',theme=theme(plot.title=element_text(size=9,face='bold',family='Arial'),plot.subtitle=element_text(size=7.5,family='Arial'),plot.margin=margin(2,2,2,2,'mm')))
save('Figure6G_candidate_evidence',p,174,98,g)
check('G all intervals and crossing-zero result retained',nrow(g)==18&sum(g$Lower<0&g$Upper>0)==1&nrow(rank)==9&all(rank$TopThreeRetained<=7))

# S: full error destinations remain available with common 0–100% cell shading.
cm<-readbr('confusion_matrices')[Model%chin%c('Anchor3','CompactWorkflow')];cm[,Truth:=factor(Truth,levels=rev(cls))];cm[,Model:=factor(Model,levels=c('Anchor3','CompactWorkflow'),labels=c('Duct RNA','Compact workflow'))]
p<-ggplot(cm,aes(Prediction,Truth,fill=RowFraction))+geom_tile(color='white',linewidth=.8)+geom_text(aes(label=sprintf('%d\n%.1f%%',Count,100*RowFraction),color=RowFraction>.65),size=8/.pt,lineheight=1.1)+
 scale_color_manual(values=c('FALSE'='#293D49','TRUE'='white'),guide='none')+scale_fill_gradient(low='#F3F6F7',high='#305D75',limits=c(0,1),guide='none')+
 facet_wrap(~Model,nrow=1)+coord_fixed()+scale_y_discrete(labels=c(C1='C1 · 112',C2='C2 · 223',C3='C3 · 162'))+
 labs(title='Complete classification error structure',x='Predicted subtype',y='Reference subtype · n')+theme(axis.line=element_blank(),axis.ticks=element_blank())
save('Supplement_confusion_matrices',p,174,77,cm,'Supplement')
check('Confusion matrices retain 994 predictions',sum(cm$Count)==994&nrow(cm)==18)

# S: same row-wise distribution grammar for expression and drug scores.
distplot<-function(a,key,order,title,subtitle,xlab,labels=NULL){
 a[,Item:=factor(get(key),levels=rev(order))];cnt<-a[,.(N=.N),by=.(Item,FinalSubtype)]
 if(is.null(labels))labels<-setNames(order,order)
 maxx<-max(a$Z);minx<-min(a$Z);span<-maxx-minx;cnt[,TextX:=maxx+span*.07]
 variable_n<-any(cnt[,uniqueN(N),by=FinalSubtype]$V1>1)
 facetlabels<-if(variable_n)setNames(cls,cls)else setNames(paste0(cls,' · n=',cnt[match(cls,FinalSubtype),N]),cls)
 p<-ggplot(a,aes(Z,Item,color=FinalSubtype))+geom_vline(xintercept=0,color='#DCE1E5',linewidth=.3)+
 geom_point(position=position_jitter(width=0,height=.17,seed=20261101),size=.32,alpha=.23)+
 geom_boxplot(aes(group=Item),orientation='y',width=.45,outlier.shape=NA,fill=NA,linewidth=.45)+
 facet_wrap(~FinalSubtype,nrow=1,labeller=labeller(FinalSubtype=facetlabels))+scale_color_manual(values=pal,guide='none')+scale_y_discrete(labels=labels)+
 scale_x_continuous(limits=c(minx-span*.03,maxx+span*if(variable_n).15 else .03),expand=expansion(mult=0))+
 labs(title=title,subtitle=subtitle,x=xlab,y=NULL)+theme(panel.spacing=grid::unit(4,'mm'),axis.ticks.y=element_blank(),axis.line.y=element_blank())
 if(variable_n)p<-p+geom_text(data=cnt,aes(x=TextX,label=N),color='#6A7780',hjust=.5,size=6.8/.pt)
 list(p=p,data=a,counts=cnt)
}
anch<-rd(file.path(old,'tables/Supplement_duct_anchor_expression.csv'));ao<-c('S100P','MUC5AC','CRP','CDH2','MUC6','TFF1','NCAM1')
x<-distplot(anch,'Gene',ao,'Duct-related RNA markers','Points: patients · boxes: median / IQR · right labels: n','Within-study centered expression (SD)')
save('Supplement_duct_marker_distributions',x$p,174,108,x$data,'Supplement');fwrite(x$counts,file.path(out,'tables/duct_gene_denominators.csv'));check('Seven genes retain all 3287 observations',nrow(anch)==3287)
ex<-rd(file.path(old,'tables/Supplement_compact17_expression.csv'));weights<-rd('5-ClassifierValidation/runs/bulk-rnaseq__20261031_step5_minpanel_v1/tables/full_candidate_weights.csv');setorder(weights,Rank)
x<-distplot(ex,'Gene',weights$Gene,'Expression of the 17-gene development candidate','Weight-rank order · all patients · descriptive expression, not validation','Within-study centered expression (SD)')
save('Supplement_compact17_distributions',x$p,174,190,x$data,'Supplement');check('17 genes retain 8449 observations',nrow(ex)==8449&setequal(ex$Gene,weights$Gene))
x<-distplot(z,'Drug',cand$Drug,'Patient distributions of all nine candidates','All values, including extremes · no additional significance tests','Relative predicted drug score (SD)',short)
save('Supplement_drug_distributions',x$p,174,126,x$data,'Supplement');check('Drug distributions retain all 4473 observations',nrow(x$data)==4473)

s<-perf[Metric%chin%met];s[,Model:=factor(Model,levels=c('Anchor3','Anchor4','CompactWorkflow','FullWorkflow'))];s[,Metric:=factor(Metric,levels=met,labels=c('Balanced accuracy','C1 recall','C2 recall','C3 recall'))]
p<-ggplot(s,aes(Estimate*100,Model,color=Model))+geom_vline(xintercept=c(0,100),color='#E3E7EA',linewidth=.3)+
 geom_errorbar(aes(xmin=Lower95*100,xmax=Upper95*100),orientation='y',width=.10,linewidth=.4)+geom_point(size=1.5)+geom_text(aes(label=sprintf('%.1f',Estimate*100)),nudge_y=.25,size=7/.pt,show.legend=FALSE)+
 facet_wrap(~Metric,nrow=1)+scale_y_discrete(limits=rev(levels(s$Model)),labels=rolelab)+scale_x_continuous(limits=c(-7,112),breaks=c(0,50,100))+
 scale_color_manual(values=role,guide='none')+labs(title='Sensitivity to the RNA reference and classifier scope',subtitle='Same 497 patients · estimates and conditional 95% intervals',x='OOF performance (%)',y=NULL)+
 theme(strip.text=element_text(size=7.5),axis.line.y=element_blank(),axis.ticks.y=element_blank(),panel.spacing=grid::unit(3,'mm'))
save('Supplement_reference_sensitivity',p,174,78,s,'Supplement');check('Sensitivity preserves 16 estimates',nrow(s)==16)
fwrite(rbindlist(manifest),file.path(out,'new_figure_manifest.csv'));fwrite(rbindlist(checks),file.path(out,'qa/numeric_checks.csv'))
fwrite(data.table(Path=inputs,SHA256=vapply(inputs,function(p)digest(file=p,algo='sha256'),character(1))),file.path(out,'tables/input_hashes.csv'))
