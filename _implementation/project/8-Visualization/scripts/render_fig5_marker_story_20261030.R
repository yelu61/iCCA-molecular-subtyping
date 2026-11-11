# Marker-focused alternatives. Descriptive evidence is kept distinct from CV.
suppressPackageStartupMessages({library(ggplot2);library(data.table);library(patchwork);library(ggrepel);library(ComplexHeatmap);library(circlize);library(uwot);library(grid);library(jsonlite);library(svglite);library(ragg)})
source("8-Visualization/scripts/render_fig5_consensus_20261030.R")
MARKER_REVIEW<-"8-Visualization/refinements/fig5_marker_story_20261030_v2"
MARKER_SIZE_RUN<-"5-ClassifierValidation/runs/bulk-rnaseq__20261030_step5_marker_size_v1"
width_mm <- 174
# Metrics for ComplexHeatmap on its temporary PostScript measurement device.
# Cairo and SVG exports resolve the installed Arial family directly.
grDevices::pdfFonts(Arial=grDevices::pdfFonts("Helvetica")[[1]])
marker_path<-function(...)file.path(STEP8_ROOT,MARKER_REVIEW,...)
for(d in c("figures","source_data","review"))dir.create(marker_path(d),recursive=TRUE,showWarnings=FALSE)
marker_prepare<-function(){
 fig5_check();support<-fig5_read("final_feature_method_support");freq<-fig5_read("feature_selection_stability");model<-readRDS(file.path(FIG5_RUN,"models/candidate_models.rds"))$Consensus_SVM
 evidence<-merge(support,freq,by="Gene",suffixes=c("Full","Folds"));evidence[,EvidenceScore:=(LimmaFraction+RFFraction+LASSOFraction)/3];evidence[,FinalSelected:=Gene%in%model$genes];setorder(evidence,-EvidenceScore,FRank,Gene);evidence[,EvidenceRank:=seq_len(.N)];stopifnot(all(evidence$EvidenceRank>0),sum(evidence$FinalSelected)==114)
 fwrite(evidence,marker_path("source_data/marker_priority.csv"))
 regions<-support[Methods>0,.(N=.N),by=.(Limma,RF,LASSO)];fwrite(regions,marker_path("source_data/marker_overlap_regions.csv"))
 meta<-fig5_read("population");meta[,Sample:=as.character(Sample)];meta[,Subtype:=factor(Subtype,levels=c("C1","C2","C3"))]
 env<-new.env(parent=emptyenv());load("1-Clustering/data/combined_expr_before_correction.rds",env);genes<-model$genes;z<-t(env$combined_expr[genes,meta$Sample,drop=FALSE]);rm(env)
 for(st in unique(meta$Study)){ii<-which(meta$Study==st);z[ii,]<-sweep(z[ii,,drop=FALSE],2,model$study_offsets[[st]][genes],"-")}
 z<-sweep(sweep(z,2,model$center[genes],"-"),2,model$scale[genes],"/");stopifnot(nrow(z)==497,ncol(z)==114,all(is.finite(z)))
 core<-support[Methods==3,.(Gene,FRank)];stopifnot(nrow(core)==20,all(core$Gene%in%colnames(z)))
 means<-sapply(c("C1","C2","C3"),function(cl)colMeans(z[meta$Subtype==cl,core$Gene,drop=FALSE]));core[,HighestMeanSubtype:=c("C1","C2","C3")[max.col(means,ties.method="first")]];core[,HighestMeanSubtype:=factor(HighestMeanSubtype,levels=c("C1","C2","C3"))];setorder(core,HighestMeanSubtype,FRank)
 fwrite(core,marker_path("source_data/core20_gene_order.csv"));fwrite(data.table(Sample=meta$Sample,Study=meta$Study,Subtype=meta$Subtype,z),marker_path("source_data/training_standardized_114_expression.csv"))
 umapfile<-marker_path("source_data/marker114_umap_coordinates.csv")
 if(file.exists(umapfile)){um<-fread(umapfile,colClasses=c(Sample="character"));stopifnot(identical(um$Sample,meta$Sample))}else{set.seed(20261030);xy<-uwot::umap(z,n_neighbors=30,min_dist=.3,n_components=2,metric="euclidean",n_threads=1,n_sgd_threads=1,verbose=FALSE);um<-data.table(Sample=meta$Sample,Study=meta$Study,Subtype=meta$Subtype,UMAP1=xy[,1],UMAP2=xy[,2]);fwrite(um,umapfile)}
 saveRDS(list(evidence=evidence,support=support,regions=regions,model=model,meta=meta,z=z,core=core,umap=um),marker_path("source_data/display_data.rds"))
}
marker_overlap<-function(d){function(){
 grid.newpage();pushViewport(viewport(xscale=c(0,1),yscale=c(0,1)))
 grid.text("Three-method marker overlap",.5,.965,gp=gpar(fontfamily="Arial",fontsize=9,fontface="bold",col=FIG5_INK))
 fill<-c("#A8C5D6","#BDD5C4","#D9C3CF");centres<-list(c(.37,.64),c(.63,.64),c(.5,.43))
 for(i in 1:3)grid.circle(centres[[i]][1],centres[[i]][2],r=.268,gp=gpar(fill=adjustcolor(fill[i],alpha.f=.55),col=FIG5_INK,lwd=.8))
 text<-function(s,x,y,bold=FALSE,size=10)grid.text(s,x,y,gp=gpar(fontfamily="Arial",fontsize=size,col=FIG5_INK,fontface=if(bold)"bold" else "plain"))
 get<-function(l,r,a)d$regions[Limma==l & RF==r & LASSO==a,N]
 text("limma\nn = 129",.18,.885,TRUE,8);text("RF\nn = 200",.82,.885,TRUE,8);text("LASSO  n = 95",.5,.13,TRUE,8)
 text(get(TRUE,FALSE,FALSE),.25,.66);text(get(FALSE,TRUE,FALSE),.75,.66);text(get(FALSE,FALSE,TRUE),.5,.28)
 text(get(TRUE,TRUE,FALSE),.5,.76);text(get(TRUE,FALSE,TRUE),.35,.46);text(get(FALSE,TRUE,TRUE),.65,.46);text(get(TRUE,TRUE,TRUE),.5,.56,TRUE,11)
 text("≥2 methods: 114 genes  |  All three: 20 genes",.5,.075,TRUE,8)
 text("Full-development selection; circle areas are schematic.",.5,.035,FALSE,7);popViewport()
}}
marker_priority<-function(d){
 e<-d$evidence;top<-e[FinalSelected==TRUE][1:12]
 p<-ggplot(e,aes(EvidenceRank,EvidenceScore))+geom_line(colour="grey75",linewidth=.35)+geom_point(data=e[FinalSelected==TRUE],colour=FIG5_INK,size=1)+geom_text_repel(data=top[1:3],aes(label=Gene),size=2.7,seed=20261030,max.overlaps=Inf,box.padding=.25,point.padding=.15,min.segment.length=0,segment.colour="grey60",segment.size=.25,nudge_y=.035)+scale_x_log10(breaks=c(1,10,100,1000,10000),labels=c("1","10","100","1,000","10,000"))+scale_y_continuous(limits=c(0,1.14),breaks=c(0,.25,.5,.75,1))+labs(title="Genome-wide evidence ranking",x="Gene rank (high to low support; log scale)",y="Mean fold-wise selection frequency")+fig5_theme()+theme(plot.title=element_text(size=8))
 tab<-melt(top[,.(Gene,LimmaFraction,RFFraction,LASSOFraction)],id.vars="Gene",variable.name="Method",value.name="Frequency");tab[,Gene:=factor(Gene,levels=rev(top$Gene))];tab[,Method:=factor(Method,levels=c("LimmaFraction","RFFraction","LASSOFraction"),labels=c("limma","RF","LASSO"))]
 q<-ggplot(tab,aes(Method,Gene,fill=Frequency))+geom_tile(colour="white",linewidth=.4)+geom_text(aes(label=sprintf("%.2f",Frequency),colour=ifelse(Frequency>.6,"white",FIG5_INK)),size=2.6)+scale_colour_identity()+scale_fill_gradient(low="white",high=FIG5_INK,limits=c(0,1),guide="none")+labs(title="Top panel genes: fold support",x=NULL,y=NULL)+fig5_theme()+theme(axis.line=element_blank(),axis.ticks=element_blank(),plot.title=element_text(size=8),axis.text.y=element_text(size=7.5))
 fig5_prepare((p+q+plot_layout(widths=c(1.65,1)))+plot_annotation(title="Consensus marker prioritization",caption="Score = mean selection frequency across three methods and 25 outer fits. Dark points: final 114-gene panel. Ranking is descriptive.",theme=fig5_theme()))
}
marker_expression<-function(d){function(){
 genes<-d$core$Gene;ix<-order(d$meta$Subtype,d$meta$Study,d$meta$Sample);mat<-t(d$z[ix,genes,drop=FALSE]);meta<-d$meta[ix];study_pal<-load_study_palette()$colors
 annotation<-HeatmapAnnotation(Subtype=meta$Subtype,Study=meta$Study,col=list(Subtype=FIG5_SUBTYPE,Study=study_pal),simple_anno_size=unit(2.5,"mm"),annotation_name_gp=gpar(fontfamily="Arial",fontsize=7),annotation_legend_param=list(Subtype=list(title_gp=gpar(fontsize=7,fontface="bold"),labels_gp=gpar(fontsize=7)),Study=list(title_gp=gpar(fontsize=7,fontface="bold"),labels_gp=gpar(fontsize=7))))
 colors<-circlize::colorRamp2(c(-2.5,0,2.5),STEP8_CONFIG$effect_colors)
 ht<-Heatmap(mat,name="Training z",col=colors,top_annotation=annotation,cluster_rows=FALSE,cluster_columns=FALSE,row_split=d$core$HighestMeanSubtype,column_split=meta$Subtype,row_title_rot=0,row_title_gp=gpar(fontfamily="Arial",fontsize=8,fontface="bold"),column_title_gp=gpar(fontfamily="Arial",fontsize=8,fontface="bold"),show_column_names=FALSE,row_names_gp=gpar(fontfamily="Arial",fontsize=8),row_gap=unit(1.2,"mm"),column_gap=unit(1.2,"mm"),use_raster=FALSE,heatmap_legend_param=list(title_gp=gpar(fontfamily="Arial",fontsize=7,fontface="bold"),labels_gp=gpar(fontfamily="Arial",fontsize=7),at=c(-2.5,0,2.5)))
 draw(ht,heatmap_legend_side="right",annotation_legend_side="right",merge_legends=TRUE,padding=unit(c(13,2,10,2),"mm"))
 grid.text("Expression of the 20 three-method core genes",x=unit(3,"mm"),y=unit(1,"npc")-unit(3,"mm"),just=c("left","top"),gp=gpar(fontfamily="Arial",fontsize=9,fontface="bold",col=FIG5_INK))
 grid.text("497 development tumors; rows grouped by highest reference-class mean. Full classifier: 114 genes.",x=unit(3,"mm"),y=unit(5,"mm"),just="left",gp=gpar(fontfamily="Arial",fontsize=7,col="grey25"))
}}
marker_umap<-function(d){ggplot(d$umap,aes(UMAP1,UMAP2,colour=Subtype))+geom_point(size=1.15,alpha=.85)+scale_colour_manual(values=FIG5_SUBTYPE,labels=c(C1="C1 (n = 112)",C2="C2 (n = 223)",C3="C3 (n = 162)"))+labs(title="Development-sample view of classifier features",subtitle="UMAP of 114 genes; all 497 tumors",x="UMAP 1",y="UMAP 2",colour=NULL,caption="Descriptive layout after supervised feature selection; not independent validation.")+fig5_theme()+theme(legend.position="top")}
marker_performance<-function(){
 met<-fig5_read("nested_performance_ci")[Model==FIG5_PRIMARY];r<-fig5_read("roc_coordinates")[Model==FIG5_PRIMARY];r[,FPR:=1-Specificity]
 auc<-setNames(vapply(c("C1","C2","C3"),function(cl)sprintf("%s  AUC = %.4f",cl,met[Metric==paste0(cl,"_AUROC"),Estimate]),character(1)),c("C1","C2","C3"))
 p<-ggplot(r,aes(FPR,Sensitivity,colour=Class))+geom_abline(slope=1,intercept=0,linetype=2,linewidth=.3,colour="grey75")+geom_step(linewidth=.65)+scale_colour_manual(values=FIG5_SUBTYPE,labels=auc)+scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1))+scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1))+coord_equal()+labs(title="One-versus-rest ROC",x="False-positive rate",y="Sensitivity",colour=NULL)+fig5_theme()+theme(legend.position="bottom",legend.direction="vertical",plot.title=element_text(size=8),legend.key.height=unit(3,"mm"))
 q<-fig5_confusion()$patches$plots[[1]]+labs(title="Reference versus predicted")+theme(plot.title=element_text(size=8))
 a<-met[Metric=="Accuracy"];caption<-sprintf("Accuracy %.2f%% (95%% interval %.2f–%.2f%%); balanced accuracy %.2f%%; macro-F1 %.3f. Conditional OOF intervals.",100*a$Estimate,100*a$Lower95,100*a$Upper95,100*met[Metric=="BalancedAccuracy",Estimate],met[Metric=="MacroF1",Estimate])
 ((p+q+plot_layout(widths=c(1,1)))+plot_annotation(title="Consensus SVM: overall OOF performance",subtitle="Five-times repeated nested CV; 497 reference-labelled tumors",caption=paste(strwrap(caption,width=94),collapse="\n"),theme=fig5_theme()))
}
marker_size_plot<-function(){
 stopifnot(isTRUE(fromJSON(file.path(MARKER_SIZE_RUN,"verification.json"))$passed));d<-fread(file.path(MARKER_SIZE_RUN,"tables/marker_size_performance.csv"));d[,Metric:=factor(Metric,levels=c("MacroAUROC","BalancedAccuracy","MacroF1"),labels=c("Macro-AUROC","Balanced accuracy","Macro-F1"))]
 colors<-c("Macro-AUROC"="#4C6A7B","Balanced accuracy"="#8C6D31","Macro-F1"="#76627F")
 p<-ggplot(d,aes(K,Estimate,colour=Metric,shape=Metric,linetype=Metric))+geom_ribbon(aes(ymin=Lower95,ymax=Upper95,fill=Metric),colour=NA,alpha=.07,show.legend=FALSE)+geom_line(linewidth=.65)+geom_point(size=2)+scale_colour_manual(values=colors)+scale_fill_manual(values=colors)+scale_x_log10(breaks=sort(unique(d$K)))+scale_y_continuous(limits=c(max(0,floor(min(d$Lower95)*20)/20),1.005),breaks=scales::breaks_pretty(5))+labs(title="Classification performance versus marker count",subtitle="Training-fold selection and tuning repeated at each K",x="Number of training-selected consensus genes (K)",y="OOF performance",colour=NULL,shape=NULL,linetype=NULL,caption="Exploratory size sensitivity; 95% conditional patient-bootstrap intervals. No optimal K or replacement model was selected.")+fig5_theme()+theme(legend.position="top")
 fig5_prepare(p)
}
render_marker_story<-function(){
 marker_prepare();d<-readRDS(marker_path("source_data/display_data.rds"))
 plots<-list(Marker_overlap=marker_overlap(d),Marker_priority=marker_priority(d),Marker_core_expression=marker_expression(d),Marker_UMAP=marker_umap(d),Overall_model_performance=marker_performance())
 dims<-list(Marker_overlap=c(114,100),Marker_priority=c(174,100),Marker_core_expression=c(174,112),Marker_UMAP=c(114,100),Overall_model_performance=c(174,100),Marker_number_performance=c(174,95))
 if(file.exists(file.path(MARKER_SIZE_RUN,"verification.json")))plots$Marker_number_performance<-marker_size_plot()
 for(n in names(plots)){p<-plots[[n]];w<-dims[[n]][1];h<-dims[[n]][2];draw<-function(){if(is.function(p))p() else print(p)};f<-marker_path("figures",n)
  cairo_pdf(paste0(f,".pdf"),width=w/25.4,height=h/25.4,family="Arial");draw();dev.off()
  svglite(paste0(f,".svg"),width=w/25.4,height=h/25.4);draw();dev.off()
  ragg::agg_png(marker_path("review",paste0(n,".png")),width=w,height=h,units="mm",res=300,background="white");draw();dev.off();cat("Rendered",n,"\n")
 }
}
