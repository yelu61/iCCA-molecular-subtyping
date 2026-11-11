suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(digest)})
source('8-Visualization/scripts/fig_style.R')
out<-'8-Visualization/refinements/fig6_storyline_20261101_v19';old<-'8-Visualization/refinements/fig6_translational_20261101_v17'
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

# G: complete drug-associated expression evidence, all three molecular subtypes.
v18<-'8-Visualization/refinements/fig6_redesign_20261101_v18'
cand<-rd(file.path(old,'tables/nine_candidate_inventory.csv'));setorder(cand,CandidateSubtype,Drug)
ev<-rd('5-ClassifierValidation/review/20261101_drug_evidence_refresh_v1/tables/candidate_target_multiomic_evidence.csv')
a<-ev[Drug%chin%cand$Drug & !is.na(Gene)&Gene!=''];a<-merge(a,cand[,.(Drug,CandidateSubtype)],by='Drug')
a[,DrugOrder:=match(Drug,cand$Drug)];setorder(a,DrugOrder,Gene,Subtype)
rowmap<-unique(a[,.(Drug,Gene,DrugOrder,CandidateSubtype)]);rowmap[,Row:=.N:1];rowmap[,Label:=paste0(CandidateSubtype,' · ',sub('_[0-9]+$','',Drug),'  |  ',Gene)]
a<-merge(a,rowmap[,.(Drug,Gene,Row,Label)],by=c('Drug','Gene'))
g<-rbind(a[,.(Drug,Gene,Subtype,CandidateSubtype,Row,Label,Assay='RNA',Effect=RNAEffect,FDR=RNAFDR)],a[,.(Drug,Gene,Subtype,CandidateSubtype,Row,Label,Assay='Protein',Effect=ProteinEffect,FDR=ProteinFDR)])
g[,State:=fifelse(!is.finite(Effect)|!is.finite(FDR),'Unavailable',fifelse(FDR<.05&Effect>0,'Higher',fifelse(FDR<.05&Effect<0,'Lower','Not significant')))]
g[,State:=factor(State,levels=c('Higher','Lower','Not significant','Unavailable'))];g[,Assay:=factor(Assay,levels=c('RNA','Protein'))]
g[,Symbol:=c(Higher='+',Lower='−','Not significant'='·',Unavailable='×')[as.character(State)]]
p<-ggplot(g,aes(Subtype,Row,fill=State))+geom_tile(width=.89,height=.84,color='#CDD4D9',linewidth=.2)+geom_text(aes(label=Symbol,color=State),size=8/.pt)+
 facet_grid(.~Assay)+scale_y_continuous(breaks=rowmap$Row,labels=rowmap$Label,expand=expansion(add=.6))+
 scale_x_discrete(position='top',expand=expansion(add=.6))+scale_fill_manual(values=c(Higher='#BD7588',Lower='#6596AE','Not significant'='white',Unavailable='#ECEFF1'),drop=FALSE)+
 scale_color_manual(values=c(Higher='white',Lower='white','Not significant'='#71808A',Unavailable='#71808A'),guide='none')+
 labs(title='Molecular context of all nine drug candidates',subtitle='Mechanism-associated genes · subtype versus other tumors',x=NULL,y=NULL,caption='Direction shown at source FDR < 0.05; expression does not establish dependency.')+
 theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_text(size=7.5),legend.position='bottom',legend.key.size=grid::unit(3,'mm'),legend.spacing.x=grid::unit(1,'mm'))
save('Figure6G_molecular_evidence',p,174,139,g)
check('All nine drugs,15 genes and90 assay cells retained',uniqueN(g$Drug)==9 &uniqueN(g$Gene)==15 & nrow(g)==90)
check('Missing differs from non-significant',any(g$State=='Unavailable')&any(g$State=='Not significant'))

# H: exploratory C3 close; retain both mapped genes and all three subtypes.
z<-rd(file.path(v18,'tables/Figure6F_drug_heterogeneity.csv'))[Drug=='Ulixertinib_1908'];z[,FinalSubtype:=factor(FinalSubtype,levels=cls)]
contr<-rd(file.path(v18,'tables/Figure6G_candidate_evidence.csv'))[Drug=='Ulixertinib_1908'];fwrite(contr,file.path(out,'tables/Figure6H_pairwise_prediction_intervals.csv'))
x<-ggplot(z,aes(FinalSubtype,Z,color=FinalSubtype))+geom_point(position=position_jitter(width=.17,height=0,seed=20261101),size=.45,alpha=.32)+geom_boxplot(width=.52,outlier.shape=NA,fill=NA,linewidth=.45)+
 scale_color_manual(values=pal,guide='none')+scale_x_discrete(labels=c(C1='C1\n112',C2='C2\n223',C3='C3\n162'))+
 labs(title='Ulixertinib prediction',subtitle='All tumors; lower relative score',x='Subtype / n',y='Relative predicted score (SD)')+theme(plot.title=element_text(size=8),plot.subtitle=element_text(size=7))
h<-g[Drug=='Ulixertinib_1908'];h[,Y:=ifelse(Gene=='MAPK1',2,1)+c(C1=.21,C2=0,C3=-.21)[Subtype]];h[,Significant:=FDR<.05]
y<-ggplot(h,aes(Effect,Y,color=Subtype,shape=Significant))+geom_vline(xintercept=0,linetype=2,linewidth=.3,color='#ACB8C0')+geom_point(size=2)+
 facet_wrap(~Assay,nrow=1,scales='free_x',labeller=labeller(Assay=c(RNA='RNA · log2 effect',Protein='Protein · SD effect')))+scale_color_manual(values=pal)+scale_shape_manual(values=c('TRUE'=16,'FALSE'=1),labels=c('FALSE'='FDR ≥ 0.05','TRUE'='FDR < 0.05'))+
 scale_y_continuous(breaks=c(2,1),labels=c('MAPK1','MAPK3'),limits=c(.5,2.5))+scale_x_continuous(expand=expansion(mult=.15))+
 labs(title='Both mapped ERK genes',subtitle='Each subtype versus other tumors',x='Study-adjusted expression difference',y=NULL)+
 theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.title=element_text(size=8),plot.subtitle=element_text(size=7),strip.text=element_text(size=7),legend.text=element_text(size=7),legend.key.width=grid::unit(3,'mm'),legend.box='vertical',legend.margin=margin(0,0,0,0))
p<-x+y+plot_layout(widths=c(1,1.5))+plot_annotation(title='A C3-focused hypothesis for experimental follow-up',subtitle='Ulixertinib / MAPK1–MAPK3 · predictions and abundance, not demonstrated efficacy',theme=theme(plot.title=element_text(size=9,face='bold',family='Arial'),plot.subtitle=element_text(size=7.5,family='Arial'),plot.margin=margin(2,2,2,2,'mm')))
save('Figure6H_C3_hypothesis',p,174,86,h);fwrite(z,file.path(out,'tables/Figure6H_patient_predictions.csv'))
check('H retains497 patients and both genes in all subtypes and assays',nrow(z)==497&nrow(h)==12&uniqueN(z$Sample)==497)
check('C3 MAPK3 corroboration and MAPK1 null RNA preserved',h[Subtype=='C3'&Gene=='MAPK3',all(Effect>0&FDR<.05)]&h[Subtype=='C3'&Gene=='MAPK1'&Assay=='RNA',FDR>=.05])
fwrite(rbindlist(manifest),file.path(out,'new_figure_manifest.csv'));fwrite(rbindlist(checks),file.path(out,'qa/numeric_checks.csv'))
fwrite(data.table(Path=inputs,SHA256=vapply(inputs,function(p)digest(file=p,algo='sha256'),character(1))),file.path(out,'tables/input_hashes.csv'))
