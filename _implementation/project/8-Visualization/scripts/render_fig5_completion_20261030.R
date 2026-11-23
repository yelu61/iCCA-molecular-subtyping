# Figure5 completion: display the fixed primary and exploratory comparisons separately.
source("8-Visualization/scripts/render_fig5_consensus_20261030.R")
suppressPackageStartupMessages({library(ComplexHeatmap);library(circlize);library(e1071)})
grDevices::pdfFonts(Arial=grDevices::pdfFonts("Helvetica")[[1]])
FIG5_COMPLETION<-"8-Visualization/refinements/fig5_complete_20261030_v4"
FIG5_STRATEGY_RUN<-"5-ClassifierValidation/runs/bulk-rnaseq__20261030_step5_feature_strategy_v1"
width_mm <- 174
completion_path<-function(...)file.path(STEP8_ROOT,FIG5_COMPLETION,...)
for(d in c("figures","source_data","review/pdf_rendered","review/qa","before"))dir.create(completion_path(d),recursive=TRUE,showWarnings=FALSE)
COMPLETION_NAMES<-c(A="Figure5A_training_schematic",B="Figure5B_complete114_feature_evidence",C="Figure5C_feature_strategy_comparison",D="Figure5D_oof_performance",E="Figure5E_study_transfer",F="Figure5F_external_expression")
COMPLETION_DIMS<-list(A=c(174,68),B=c(174,120),C=c(174,105),D=c(174,94),E=c(174,108),F=c(174,128))

completion_strategy<-function(){
 v<-fromJSON(file.path(FIG5_STRATEGY_RUN,"verification.json"));stopifnot(isTRUE(v$passed),v$primary_probability_error<1e-10)
 p<-fread(file.path(FIG5_STRATEGY_RUN,"tables/performance_ci.csv"))[Metric=="BalancedAccuracy"]
 n<-fread(file.path(FIG5_STRATEGY_RUN,"tables/feature_count_summary.csv"));p<-merge(p,n,by=c("Spec","Budget","Strategy"))
 order<-c("Limma","RF","LASSO","Consensus2","Core3");labels<-c(Limma="limma",RF="Random forest",LASSO="LASSO",Consensus2="≥2 methods",Core3="All three")
 p[,Strategy:=factor(Strategy,levels=rev(order))];p[,Primary:=ifelse(Strategy=="Consensus2","Prespecified consensus","Other strategies")]
 p[,Budget:=factor(Budget,levels=c("Native","CoreMatched"),labels=c("Native feature sets","Same count within each training set"))]
 p[,NText:=sprintf("%g [%d–%d]",MedianFeatures,MinFeatures,MaxFeatures)]
 fwrite(p,completion_path("source_data/5C_strategy_display.csv"))
 q<-ggplot(p,aes(y=Strategy,x=Estimate,colour=Primary))+geom_errorbar(aes(xmin=Lower95,xmax=Upper95),orientation="y",width=.18,linewidth=.5)+geom_point(size=2.1)+geom_text(aes(x=1.018,label=NText),colour=FIG5_INK,size=2.45,hjust=0)+scale_colour_manual(values=c("Prespecified consensus"=FIG5_INK,"Other strategies"="#8196A5"),guide="none")+scale_y_discrete(labels=labels)+scale_x_continuous(limits=c(.925,1.065),breaks=c(.94,.96,.98,1),labels=scales::label_percent(accuracy=1))+facet_wrap(~Budget,nrow=1)+labs(title="Feature-selection strategies under identical nested CV",subtitle="Linear SVM throughout; feature count = median [range] across 25 outer fits",x="OOF balanced accuracy",y=NULL,caption="Exploratory comparison; 95% conditional patient-bootstrap intervals. Same-count budget follows the training-only three-method intersection size.\nWithin-set ordering uses training-only limma F rank. No strategy or optimal gene count was selected from these results.")+fig5_theme()+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(size=7.2),strip.text=element_text(size=7),plot.caption=element_text(size=6.6))
 q
}

completion_transfer<-function(){
 d<-fig5_read("loso_predictions")[Model==FIG5_PRIMARY]
 counts<-d[,.(Correct=sum(Prediction==Truth),Total=.N),by=.(Study,Truth)];counts[,Recall:=Correct/Total]
 summaries<-fig5_read("loso_per_study_performance")[Model==FIG5_PRIMARY & Metric=="BalancedAccuracy",.(Study,N,LOSO=Estimate)]
 oof<-fig5_read("nested_oof_predictions")[Model==FIG5_PRIMARY]
 inside<-oof[,.(OOF=mean(vapply(c("C1","C2","C3"),function(cl)mean(Prediction[Truth==cl]==cl),numeric(1)))),by=Study]
 summaries<-merge(summaries,inside,by="Study");setorder(summaries,Study);summaries[,Label:=sprintf("%s (n = %d)",Study,N)];summaries[,Label:=factor(Label,levels=rev(Label))]
 counts<-merge(counts,summaries[,.(Study,Label)],by="Study");counts[,Truth:=factor(Truth,levels=c("C1","C2","C3"))];counts[,Text:=sprintf("%d/%d\n%.1f%%",Correct,Total,100*Recall)]
 stopifnot(sum(counts$Total)==497,nrow(counts)==21,all(counts$Total>0))
 fwrite(counts,completion_path("source_data/5E_study_class_recall.csv"));fwrite(summaries,completion_path("source_data/5E_study_balanced_accuracy.csv"))
 p<-ggplot(counts,aes(Truth,Label,fill=Recall))+geom_tile(colour="white",linewidth=.8)+geom_text(aes(label=Text,colour=ifelse(Recall>.8,"white",FIG5_INK)),size=2.55,lineheight=1)+scale_colour_identity()+scale_fill_gradient(low="#ECF2F6",high="#4E8DB1",limits=c(.5,1),oob=scales::squish,guide="none")+labs(title="Class recall in the held-out study",x="Reference subtype",y=NULL)+fig5_theme()+theme(axis.line=element_blank(),axis.ticks=element_blank(),plot.title=element_text(size=8),axis.text.y=element_text(size=7))
 long<-melt(summaries,id.vars=c("Study","N","Label"),measure.vars=c("OOF","LOSO"),variable.name="Evaluation",value.name="BalancedAccuracy")
 q<-ggplot(summaries,aes(y=Label))+geom_vline(xintercept=1,colour="grey85",linewidth=.3)+geom_segment(aes(x=OOF,xend=LOSO,yend=Label),colour="grey70",linewidth=.7)+geom_point(data=long,aes(x=BalancedAccuracy,shape=Evaluation,colour=Evaluation),size=2.3)+scale_shape_manual(values=c(OOF=1,LOSO=16),labels=c(OOF="Internal OOF",LOSO="Held-out study"))+scale_colour_manual(values=c(OOF="#9AA9B4",LOSO=FIG5_INK),labels=c(OOF="Internal OOF",LOSO="Held-out study"))+scale_x_continuous(limits=c(.84,1.025),breaks=c(.85,.9,.95,1),labels=scales::label_percent(accuracy=1))+labs(title="Balanced accuracy",x="Balanced accuracy",y=NULL,colour=NULL,shape=NULL)+fig5_theme()+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),plot.title=element_text(size=8),legend.position="bottom",legend.text=element_text(size=6.6))
 (p+q+plot_layout(widths=c(1.7,1)))+plot_annotation(title="Classifier transfer across seven studies",subtitle="Each study excluded from feature selection, model training and tuning",caption="Left: correct / reference-class total. Right: point estimates; internal OOF and held-out study.\nUnseen-study centering uses its unlabeled cohort mean (transductive); this is not an independent diagnostic validation.",theme=fig5_theme()+theme(plot.caption=element_text(size=6.6)))
}

completion_external_data<-function(){
 old<-readRDS("8-Visualization/refinements/fig5_marker_story_20261030_v2/source_data/display_data.rds")
 genes<-old$model$genes;ref<-data.table(Gene=genes)
 avg<-sapply(c("C1","C2","C3"),function(cl)colMeans(old$z[old$meta$Subtype==cl,genes,drop=FALSE]))
 ref[,ReferenceClass:=factor(c("C1","C2","C3")[max.col(avg,ties.method="first")],levels=c("C1","C2","C3"))]
 ref<-merge(ref,old$evidence[,.(Gene,ConsensusFraction,FRank)],by="Gene");setorder(ref,ReferenceClass,-ConsensusFraction,FRank,Gene)
 predictions<-fig5_read("external_sample_predictions")[Model==FIG5_PRIMARY];predictions[,Sample:=as.character(Sample)]
 cohorts<-sort(unique(predictions$Cohort));stopifnot(length(cohorts)==6,nrow(predictions)==506)
 raw<-list();for(path in c("0-Data/qc/0-Data-rnaseq-microarray1.Rdata","0-Data/qc/0-Data-rnaseq-microarray2.Rdata")){env<-new.env(parent=emptyenv());load(path,env);for(n in intersect(ls(env),cohorts))raw[[n]]<-get(n,env)}
 all<-list();wide<-list();qc<-list();maxerror<-0
 for(cohort in cohorts){x<-as.matrix(raw[[cohort]]);storage.mode(x)<-"double"
  if(anyDuplicated(rownames(x))){nn<-table(rownames(x));x<-rowsum(x,rownames(x),reorder=FALSE);x<-x/as.numeric(nn[rownames(x)])}
  z<-predictions[Cohort==cohort];stopifnot(setequal(colnames(x),z$Sample),!anyDuplicated(z$Sample));z[,Prediction:=factor(Prediction,levels=c("C1","C2","C3"))];setorder(z,Prediction,-MaximumProbability,Sample)
  available<-intersect(genes,rownames(x));v<-x[available,z$Sample,drop=FALSE];v[!is.finite(v)]<-NA_real_;mu<-rowMeans(v,na.rm=TRUE);sdv<-apply(v,1,sd,na.rm=TRUE);mu[!is.finite(mu)]<-0;sdv[!is.finite(sdv)|sdv<1e-8]<-1
  v<-sweep(sweep(v,1,mu,"-"),1,sdv,"/");v[!is.finite(v)]<-0
  model_input<-matrix(0,nrow(z),length(genes),dimnames=list(z$Sample,genes));model_input[,available]<-t(v)
  pr<-attr(predict(old$model$model,model_input,probability=TRUE),"probabilities")[,c("C1","C2","C3"),drop=FALSE];pr<-pr/rowSums(pr)
  error<-max(abs(pr-as.matrix(z[,paste0("Probability",c("C1","C2","C3")),with=FALSE])));maxerror<-max(maxerror,error)
  display<-matrix(NA_real_,114,nrow(z),dimnames=list(ref$Gene,z$Sample));display[available,]<-v
  stopifnot(length(available)==unique(z$MatchedFeatures),error<1e-10)
  all[[cohort]]<-list(mat=display,meta=z)
  long<-as.data.table(as.table(display));setnames(long,c("Gene","Sample","StudyZ"));long[,Cohort:=cohort];wide[[cohort]]<-long
  qc[[cohort]]<-data.table(Cohort=cohort,N=nrow(z),MatchedFeatures=length(available),MissingFeatures=length(setdiff(genes,available)),ProbabilityReproductionError=error)
 }
 fwrite(ref,completion_path("source_data/5F_gene_order.csv"));fwrite(rbindlist(wide),completion_path("source_data/5F_external_expression.csv"));fwrite(predictions,completion_path("source_data/5F_sample_assignments.csv"));fwrite(rbindlist(qc),completion_path("source_data/5F_display_QC.csv"))
 saveRDS(list(ref=ref,cohorts=cohorts,data=all),completion_path("source_data/5F_display_data.rds"))
 write_json(list(passed=TRUE,samples=506,panel_gene_count=length(genes),cohort_objects=6,maximum_prediction_reproduction_error=maxerror,unmeasured_display="NA / grey; not zero"),completion_path("review/external_display_verification.json"),auto_unbox=TRUE,pretty=TRUE)
}

completion_external<-function(){function(){
 d<-readRDS(completion_path("source_data/5F_display_data.rds"));ht<-NULL
 display_names<-c(log2dat_2012Andersen="2012Andersen",log2dat_2013Sia="2013Sia",log2dat_2017Jusakual1="2017Jusakual P1",log2dat_2017Jusakual2="2017Jusakual P2",log2dat_2020Job="2020Job",log2dat_TIGER="TIGER")
 for(cohort in d$cohorts){x<-d$data[[cohort]];first<-is.null(ht)
  anno<-HeatmapAnnotation(Projected=x$meta$Prediction,col=list(Projected=FIG5_SUBTYPE),show_annotation_name=first,simple_anno_size=unit(2.4,"mm"),annotation_name_gp=gpar(fontfamily="Arial",fontsize=6.5),show_legend=first,annotation_legend_param=list(title_gp=gpar(fontfamily="Arial",fontsize=7,fontface="bold"),labels_gp=gpar(fontfamily="Arial",fontsize=7)))
  one<-Heatmap(x$mat,name=if(first)"Within-study z" else cohort,col=circlize::colorRamp2(c(-2.5,0,2.5),STEP8_CONFIG$effect_colors),na_col="#B7B7B7",cluster_rows=FALSE,cluster_columns=FALSE,row_split=d$ref$ReferenceClass,row_title_rot=0,row_title_gp=gpar(fontfamily="Arial",fontsize=7,fontface="bold"),show_row_names=FALSE,show_column_names=FALSE,column_title=sprintf("%s\nn=%d; %d/114",display_names[cohort],nrow(x$meta),unique(x$meta$MatchedFeatures)),column_title_gp=gpar(fontfamily="Arial",fontsize=6.7,fontface="bold"),top_annotation=anno,width=unit(21,"mm"),row_gap=unit(1.3,"mm"),show_heatmap_legend=first,use_raster=FALSE,heatmap_legend_param=list(title_gp=gpar(fontfamily="Arial",fontsize=7,fontface="bold"),labels_gp=gpar(fontfamily="Arial",fontsize=7),at=c(-2.5,0,2.5)))
  ht<-if(first)one else ht+one
 }
 draw(ht,heatmap_legend_side="right",annotation_legend_side="right",merge_legends=TRUE,ht_gap=unit(1.6,"mm"),padding=unit(c(17,4,13,3),"mm"))
 grid.text("External microarray projection: classifier-feature expression",x=unit(3,"mm"),y=unit(1,"npc")-unit(3,"mm"),just=c("left","top"),gp=gpar(fontfamily="Arial",fontsize=9,fontface="bold",col=FIG5_INK))
 grid.text("All 506 samples; 114 gene rows ordered by development reference-class expression",x=unit(3,"mm"),y=unit(1,"npc")-unit(8,"mm"),just=c("left","top"),gp=gpar(fontfamily="Arial",fontsize=7,col=FIG5_INK))
 grid.text("Row groups: highest development mean in C1/C2/C3. Headers: samples; measured panel genes.",x=unit(3,"mm"),y=unit(11,"mm"),just="left",gp=gpar(fontfamily="Arial",fontsize=6.8,col=FIG5_INK))
 grid.text("Grey: unmeasured genes. Within-study z colors saturate at ±2.5. Projection only; no independent subtype truth.",x=unit(3,"mm"),y=unit(6,"mm"),just="left",gp=gpar(fontfamily="Arial",fontsize=6.8,col=FIG5_INK))
}}

render_fig5_completion_panel<-function(id){
 stopifnot(id%in%names(COMPLETION_NAMES));fig5_check()
 stem<-COMPLETION_NAMES[[id]];dest<-completion_path("figures",stem)
 if(id%in%c("A","B","D")){
  src<-switch(id,A="8-Visualization/figures/main/Figure5A_training_schematic",B="8-Visualization/refinements/fig5_marker_logic_20261030_v3/figures/Complete114_feature_evidence",D="8-Visualization/refinements/fig5_marker_logic_20261030_v3/figures/OOF_performance_compact")
  for(ext in c("pdf","svg"))if(file.exists(paste0(src,".",ext)))stopifnot(file.copy(paste0(src,".",ext),paste0(dest,".",ext),overwrite=TRUE))
 }else{
  if(id=="F")completion_external_data()
  p<-switch(id,C=completion_strategy(),E=completion_transfer(),F=completion_external());draw<-function(){if(is.function(p))p() else print(p)};sz<-COMPLETION_DIMS[[id]]
  cairo_pdf(paste0(dest,".pdf"),width=sz[1]/25.4,height=sz[2]/25.4,family="Arial");draw();dev.off()
  svglite(paste0(dest,".svg"),width=sz[1]/25.4,height=sz[2]/25.4);draw();dev.off()
  ragg::agg_png(completion_path("review",paste0(stem,".png")),width=sz[1],height=sz[2],units="mm",res=300,background="white");draw();dev.off()
 }
 stopifnot(system2("pdftoppm",c("-png","-singlefile","-r","150",shQuote(paste0(dest,".pdf")),shQuote(completion_path("review/pdf_rendered",stem))))==0)
 cat("Prepared Figure5",id,"\n");invisible(paste0(dest,".pdf"))
}
render_fig5_completion<-function(){for(id in names(COMPLETION_NAMES))render_fig5_completion_panel(id)}

# Active notebook entry; the accepted A is read-only, and B/D retain exact source bytes.
render_fig5_current_panel<-function(id){
 stopifnot(id%in%names(COMPLETION_NAMES))
 manifest_path<-file.path(STEP8_ROOT,"8-Visualization/tables/figure_manifest.csv")
 m<-fread(manifest_path,colClasses="character",na.strings=NULL);ix<-which(m$Figure=="Figure 5" & m$Panel==id);stopifnot(length(ix)==1)
 target<-file.path(STEP8_ROOT,m$TargetOutput[ix])
 if(id=="A"){stopifnot(file.exists(target));return(invisible(target))}
 rendered<-render_fig5_completion_panel(id);stopifnot(file.copy(rendered,target,overwrite=TRUE),file.copy(sub("[.]pdf$",".svg",rendered),sub("[.]pdf$",".svg",target),overwrite=TRUE))
 sources<-strsplit(m$SourceTable[ix],";",fixed=TRUE)[[1]];stopifnot(all(file.exists(file.path(STEP8_ROOT,sources))))
 fingerprint<-paste(paste0(sources,"=",unname(tools::md5sum(file.path(STEP8_ROOT,sources)))),collapse=";")
 m[ix,`:=`(Status="Rendered",SourceMD5=fingerprint,QAFile="",RenderedAt=format(Sys.time(),tz="UTC",usetz=TRUE))]
 tmp<-tempfile(tmpdir=dirname(manifest_path),fileext=".csv");fwrite(m,tmp,na="");stopifnot(file.rename(tmp,manifest_path))
 invisible(target)
}
