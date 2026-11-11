# Figure5: display-only rendering of the reviewed consensus classifier run.
# No model fitting, new confidence intervals, hypothesis tests or cutoffs.
suppressPackageStartupMessages({library(ggplot2);library(data.table);library(patchwork);library(grid);library(jsonlite);library(svglite);library(ragg);library(digest)})
if(!exists("STEP8_ROOT"))source("8-Visualization/scripts/fig_style.R")
stopifnot(identical(STEP8_FONT$family,"Arial"))
FIG5_WIDTH_MM<-174
width_mm <- 174
FIG5_RUN<-"5-ClassifierValidation/runs/bulk-rnaseq__20261030_step5_consensus_v1"
FIG5_REVIEW<-"8-Visualization/refinements/fig5_consensus_20261030_v1"
FIG5_PRIMARY<-"Consensus_SVM"
FIG5_INK<-"#293F50";FIG5_MID<-"#758B9B";FIG5_LIGHT<-"#EAF0F4"
FIG5_SUBTYPE<-load_subtype_palette()$colors[c("C1","C2","C3")]
FIG5_STUDY<-load_study_palette()$colors
fig5_read<-function(name)fread(file.path(STEP8_ROOT,FIG5_RUN,"tables",paste0(name,".csv")),check.names=FALSE)
fig5_check<-function(){v<-fromJSON(file.path(STEP8_ROOT,FIG5_RUN,"verification.json"));stopifnot(isTRUE(v$passed));s<-fig5_read("final_model_features");stopifnot(s[Model==FIG5_PRIMARY,.N]==114)}
fig5_theme<-function(){theme_classic(base_size=8,base_family=STEP8_FONT$family)+theme(plot.title=element_text(size=9,face="bold",colour=FIG5_INK),plot.subtitle=element_text(size=7.5,colour="grey30"),plot.title.position="plot",plot.caption=element_text(size=7,hjust=0,colour="grey30",lineheight=1.1),plot.caption.position="plot",axis.title=element_text(size=8),axis.text=element_text(size=7.5,colour="grey20"),axis.line=element_line(linewidth=.3,colour="grey45"),axis.ticks=element_line(linewidth=.25),legend.text=element_text(size=7.5),legend.title=element_text(size=7.5),legend.key.size=unit(3,"mm"),legend.position="top",strip.background=element_rect(fill=FIG5_LIGHT,colour=NA),strip.text=element_text(size=8,face="bold"),panel.spacing=unit(4,"mm"),plot.margin=margin(3,3,3,3,"mm"))}
fig5_empty<-function()theme_void(base_family=STEP8_FONT$family)+theme(plot.margin=margin(3,2,3,2,"mm"),plot.title=element_text(size=8,face="bold",colour=FIG5_INK))
fig5_schematic<-function(){
  # Diagram symbols are illustrative; the counts/support rule are source-backed.
  draw<-function(){
    grid.newpage();pushViewport(viewport(xscale=c(0,174),yscale=c(0,68)))
    txt<-function(s,x,y,size=7,bold=FALSE,just="centre",colour=FIG5_INK){grid.text(s,unit(x,"native"),unit(y,"native"),just=just,gp=gpar(fontfamily=STEP8_FONT$family,fontsize=size,fontface=if(bold)"bold" else "plain",col=colour,lineheight=1.05))}
    rect<-function(x,y,w,h,fill="white",col="grey82",lwd=.6){grid.rect(unit(x,"native"),unit(y,"native"),unit(w,"native"),unit(h,"native"),gp=gpar(fill=fill,col=col,lwd=lwd))}
    line<-function(x,y,arrow=FALSE,col=FIG5_INK,lwd=.7){grid.lines(unit(x,"native"),unit(y,"native"),gp=gpar(col=col,lwd=lwd),arrow=if(arrow)grid::arrow(length=unit(1.1,"mm"),type="closed") else NULL)}
    mat<-function(x,y,w,h,nr=7,nc=9){for(i in seq_len(nr))for(j in seq_len(nc))rect(x+(j-.5)*w/nc,y+(i-.5)*h/nr,w/nc,h/nr,fill=c("#CBD6DF","#B1C1CE","#E7EDF2")[(i+2*j)%%3+1],col="white",lwd=.25)}
    starts<-c(3,33,78.5,108,139.5);widths<-c(27.5,43,27,29,31.5);heads<-c("Data","Three-method selection","Consensus genes","Model training","Evaluation")
    for(i in 1:5){rect(starts[i]+widths[i]/2,39,widths[i],46,fill="white",col="grey90",lwd=.4);rect(starts[i]+widths[i]/2,64.5,widths[i],5,fill="#DBE4EB",col=NA);txt(heads[i],starts[i]+widths[i]/2,64.5,7.8,TRUE)}
    txt("497 principal-state\ntumors",16.75,57,7.4);txt("Expression + labels",16.75,52,7)
    mat(5.5,30,22.5,18);for(j in 1:3){rect(5.5+(j-.5)*7.5,49.4,7.5,1.8,FIG5_SUBTYPE[j],col=NA);txt(names(FIG5_SUBTYPE)[j],5.5+(j-.5)*7.5,47.1,7)};txt("Genes × samples",16.75,25.5,7)
    line(c(28.5,31.8),c(39,39));line(c(31.8,31.8),c(24,52))
    ys<-c(52,38,24);methods<-c("limma","Random forest","LASSO");subs<-c("Study-adjusted\ndifferential expression","Permutation\nimportance","Multinomial\nsparse selection")
    for(i in 1:3){line(c(31.8,34),c(ys[i],ys[i]),TRUE);rect(54.3,ys[i],40,11,fill="#F8FAFB",col="#9AAEBE");txt(methods[i],36,ys[i]+2.5,7.5,TRUE,"left");txt(subs[i],36,ys[i]-1.9,7,FALSE,"left")
      if(i<3)for(k in 1:4){rect(69,ys[i]+3.2-k*1.5,8,.75,"#D3DDE5",col=NA);rect(65+(5-k)*.65,ys[i]+3.2-k*1.5,(5-k)*1.3,.75,FIG5_MID,col=NA)}
      if(i==3){line(c(67,67),c(ys[i]-3.4,ys[i]+3.4),col="grey55");for(k in 1:4){yy<-ys[i]+3.4-k*1.4;line(c(67,67+if(k%%2)4 else 0),c(yy,yy),col=FIG5_MID,lwd=1.8)}}
      line(c(74.5,77,77,79),c(ys[i],ys[i],39,39),TRUE)
    }
    txt("Selected by ≥2 methods",92,57,7.2)
    support<-fig5_read("final_feature_method_support");setorder(support,FRank)
    chosen<-rbind(support[Methods==3][1:2],support[Methods==2][1:3],support[Methods==1][1]);stopifnot(nrow(chosen)==6)
    for(j in 1:3){txt(c("limma","RF","LASSO")[j],82.7+(j-1)*7.6,51,7);for(i in 1:6)rect(82.7+(j-1)*7.6,48-i*3.1,5.4,2.4,if(as.logical(chosen[[c("Limma","RF","LASSO")[j]]][i]))FIG5_MID else "#E2E8ED",col="white",lwd=.4)}
    line(c(80,102),c(29.4,29.4),col="#9AA8B2",lwd=.9);txt("×",103.3,29.4,7,colour="#9AA8B2")
    line(c(92,92),c(27.6,25.5),TRUE);rect(92,21.6,21,7.2,FIG5_LIGHT,col="#9AAEBE");txt("114 genes",92,23,7.8,TRUE);txt("Final panel",92,19.8,7)
    line(c(105.7,108.7),c(39,39),TRUE);txt("Linear SVM",122.5,56,8.4,TRUE);txt("Primary model",122.5,51.8,7)
    rect(122.5,38.3,20,18,fill="#F8FAFB",col=FIG5_INK);line(c(113,121,125),c(30,39,47));line(c(121,132),c(39,35))
    xx<-c(115,117,119,125,128,130,121,124,127);yy<-c(38,41,42,44,42,40,33,34,32)
    grid.points(unit(xx,"native"),unit(yy,"native"),pch=16,size=unit(1.2,"mm"),gp=gpar(col=rep(c("#4D6A80","#A9BDCC","#7E96A7"),each=3)))
    txt("RF / kNN comparators",122.5,23,7)
    line(c(137.2,139,139),c(39,39,54));line(c(139,139),c(39,24));ey<-c(54,39,24)
    for(i in 1:3)line(c(139,141.3),c(ey[i],ey[i]),TRUE)
    for(i in 1:3)for(j in 1:3)rect(142.3+(j-.5)*2,51+(i-.5)*2,2,2,fill="white",col=FIG5_MID,lwd=.45)
    for(i in 1:3)rect(142.3+(i-.5)*2.2,39,1.8,5.6,fill=if(i==3)"#B1C1CE" else FIG5_LIGHT,col=FIG5_MID,lwd=.5)
    for(i in 1:3)mat(141.6+(i-1)*2.6,21,2.1,6,nr=4,nc=2)
    txt("OOF\nevaluation",150.5,54,7.3,FALSE,"left");txt("LOSO transfer\nTransductive",150.5,39,7,FALSE,"left");txt("Study-wise\nmicroarray\nprojection",150.5,24,7,FALSE,"left")
    rect(85,7.2,104,9.5,fill=FIG5_LIGHT,col=NA)
    txt("Preprocessing, selection and tuning within training folds",85,9.3,7)
    txt("Outer: 5-fold × 5 repeats; inner: 4-fold",85,5.3,7)
    popViewport()
  };draw
}
fig5_confusion<-function(){
  d<-fig5_read("confusion_matrix")[Model==FIG5_PRIMARY];d[,Truth:=factor(Truth,levels=c("C3","C2","C1"))];d[,Prediction:=factor(Prediction,levels=c("C1","C2","C3"))];d[,Text:=sprintf("%d\n%.1f%%",Freq,100*RowFraction)];d[,TextColour:=ifelse(RowFraction>.55,"white",FIG5_INK)]
  p<-ggplot(d,aes(Prediction,Truth,fill=RowFraction))+geom_tile(colour="white",linewidth=.8)+geom_text(aes(label=Text,colour=TextColour),size=3.15,lineheight=1.1)+scale_colour_identity()+scale_fill_gradient(low="#F3F6F8",high=FIG5_INK,limits=c(0,1),guide="none")+coord_equal(expand=FALSE)+labs(x="Predicted subtype",y="Reference subtype")+fig5_theme()+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text=element_text(size=8),plot.margin=margin(2,1,2,2,"mm"))
  oof<-fig5_read("nested_oof_predictions")[Model==FIG5_PRIMARY];counts<-oof[,.N,by=Truth];met<-fig5_read("nested_performance_ci")[Model==FIG5_PRIMARY]
  rows<-rbindlist(lapply(c("C1","C2","C3"),function(cl){rec<-met[Metric==paste0(cl,"_Sensitivity"),Estimate];prec<-met[Metric==paste0(cl,"_Precision"),Estimate];data.table(Class=cl,y=match(cl,c("C3","C2","C1")),N=counts[Truth==cl,N],Recall=rec,Precision=prec,F1=2*rec*prec/(rec+prec))}))
  table<-ggplot()+geom_rect(data=rows,aes(xmin=.04,xmax=.11,ymin=y-.35,ymax=y+.35,fill=Class),colour=NA)+scale_fill_manual(values=FIG5_SUBTYPE,guide="none")+geom_text(data=rows,aes(x=.15,y=y,label=Class),hjust=0,size=2.8,fontface="bold")+geom_text(data=rows,aes(x=.43,y=y,label=N),size=2.8)+geom_text(data=rows,aes(x=.79,y=y,label=sprintf("%.1f%%",100*Recall)),size=2.8)+geom_text(data=rows,aes(x=1.22,y=y,label=sprintf("%.1f%%",100*Precision)),size=2.8)+geom_text(data=rows,aes(x=1.62,y=y,label=sprintf("%.3f",F1)),size=2.8)+annotate("text",x=c(.43,.79,1.22,1.62),y=3.65,label=c("n","Recall","Precision","F1"),size=2.8,fontface="bold")+coord_cartesian(xlim=c(0,1.85),ylim=c(.35,3.85),clip="off")+fig5_empty()
  (p+table+plot_layout(widths=c(1.05,1.05)))+plot_annotation(title="Out-of-fold subtype assignments",subtitle=sprintf("497 tumors  |  Accuracy %.2f%%  |  Balanced accuracy %.2f%%",100*met[Metric=="Accuracy",Estimate],100*met[Metric=="BalancedAccuracy",Estimate]),caption="Reference labels originate from the discovery expression resource.",theme=fig5_theme())
}
fig5_comparison<-function(include_baseline=FALSE){
  d<-fig5_read("nested_performance_ci")[Metric%in%c("BalancedAccuracy","MacroF1","MulticlassBrier")]
  order<-c("Consensus_SVM","Consensus_kNN","Consensus_RF");if(include_baseline)order<-c(order,"Limma123_SVM")
  d<-d[Model%in%order];labs<-c(Consensus_SVM="Consensus SVM",Consensus_kNN="Consensus kNN",Consensus_RF="Consensus RF",Limma123_SVM="Limma-123 SVM\nHistorical-size baseline")
  d[,ModelLabel:=factor(Model,levels=rev(order),labels=labs[rev(order)])];d[,MetricLabel:=factor(Metric,levels=c("BalancedAccuracy","MacroF1","MulticlassBrier"),labels=c("Balanced accuracy ↑","Macro-F1 ↑","Multiclass Brier ↓"))]
  cols<-setNames(c(FIG5_INK,"#7D929F","#AAB8C1","#C5CCD1")[seq_along(order)],labs[order])
  ggplot(d,aes(Estimate,ModelLabel,colour=ModelLabel))+geom_segment(aes(x=Lower95,xend=Upper95,yend=ModelLabel),linewidth=.65)+geom_point(size=2.4)+scale_colour_manual(values=cols,guide="none")+facet_wrap(~MetricLabel,scales="free_x",nrow=1)+scale_x_continuous(expand=expansion(mult=c(.15,.12)),n.breaks=4)+labs(title=if(include_baseline)"Classifier comparison including the historical-size baseline" else "Classifiers trained on consensus features",x="OOF estimate and 95% interval",y=NULL,caption="Intervals: 1,000 patient bootstraps conditional on OOF predictions.")+fig5_theme()+theme(panel.grid.major.x=element_line(colour="grey93",linewidth=.25),axis.line.y=element_blank(),axis.ticks.y=element_blank(),strip.text=element_text(size=7.8))
}
fig5_calibration<-function(){
  d<-fig5_read("calibration_bins")[Model==FIG5_PRIMARY];m<-fig5_read("calibration_metrics")[Model==FIG5_PRIMARY]
  lev<-c("C1","C2","C3");lab<-setNames(vapply(lev,function(cl){z<-m[Class==cl];sprintf("%s  |  Brier %.4f; slope %.2f",cl,z$Brier,z$CalibrationSlope)},character(1)),lev)
  d[,Class:=factor(Class,levels=lev)];d[,BinIndex:=seq_len(.N),by=Class];d[,BinMid:=(BinIndex-.5)/10]
  p<-ggplot(d[N>0],aes(Predicted,Observed,colour=Class))+geom_abline(slope=1,intercept=0,colour="grey65",linetype=2,linewidth=.35)+geom_line(linewidth=.6)+geom_point(size=2)+scale_colour_manual(values=FIG5_SUBTYPE,guide="none")+facet_wrap(~Class,nrow=1,labeller=as_labeller(lab))+scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1),expand=expansion(mult=.03))+scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1),expand=expansion(mult=.03))+labs(x=NULL,y="Observed class fraction")+fig5_theme()+theme(plot.margin=margin(2,3,0,3,"mm"),strip.text=element_text(size=7.3))
  support<-ggplot(d,aes(BinMid,N,fill=Class))+geom_col(width=.085)+scale_fill_manual(values=FIG5_SUBTYPE,guide="none")+facet_wrap(~Class,nrow=1)+scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1),expand=expansion(mult=.03))+scale_y_continuous(n.breaks=3,expand=expansion(mult=c(0,.1)))+labs(x="Predicted probability",y="Bin n")+fig5_theme()+theme(strip.text=element_blank(),strip.background=element_blank(),axis.title.y=element_text(size=7),plot.margin=margin(0,3,2,3,"mm"))
  (p/support+plot_layout(heights=c(2.5,1)))+plot_annotation(title="OOF probability calibration and bin support",caption="Fixed 0.1-wide bins; empty bins remain visible in the support axis. Probabilities are not deployment-calibrated.",theme=fig5_theme())
}
fig5_loso<-function(){
  d<-fig5_read("loso_per_study_performance")[Model==FIG5_PRIMARY & Metric=="BalancedAccuracy",.(Study,N,LOSO=Estimate)]
  oof<-fig5_read("nested_oof_predictions")[Model==FIG5_PRIMARY]
  internal<-oof[,.(OOF=mean(vapply(c("C1","C2","C3"),function(cl)mean(Prediction[Truth==cl]==cl),numeric(1)))),by=Study]
  d<-merge(d,internal,by="Study");setorder(d,Study);d[,Label:=sprintf("%s  (n = %d)",Study,N)];d[,Label:=factor(Label,levels=rev(Label))]
  long<-melt(d,id.vars=c("Study","N","Label"),measure.vars=c("OOF","LOSO"),variable.name="Evaluation",value.name="BalancedAccuracy")
  ggplot(d,aes(y=Label))+geom_vline(xintercept=1,colour="grey80",linewidth=.3)+geom_segment(aes(x=OOF,xend=LOSO,yend=Label),colour="grey75",linewidth=.65)+geom_point(data=long,aes(x=BalancedAccuracy,shape=Evaluation,colour=Evaluation),size=2.5)+scale_shape_manual(values=c(OOF=1,LOSO=16),labels=c(OOF="Internal OOF",LOSO="Held-out study"))+scale_colour_manual(values=c(OOF="#8E9FAB",LOSO=FIG5_INK),labels=c(OOF="Internal OOF",LOSO="Held-out study"))+scale_x_continuous(limits=c(.8,1.015),breaks=c(.8,.85,.9,.95,1),labels=scales::label_percent(accuracy=1))+labs(title="Transfer to held-out studies",x="Balanced accuracy",y=NULL,shape=NULL,colour=NULL,caption="Point estimates only. Held-out-study centering uses unlabeled cohort data (transductive).")+fig5_theme()+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),legend.position="top",panel.grid.major.x=element_line(colour="grey94",linewidth=.25))
}
fig5_external<-function(){
  d<-fig5_read("external_cohort_summary")[Model==FIG5_PRIMARY];setorder(d,Cohort)
  display<-sub("log2dat_","",d$Cohort);display<-sub("2017Jusakual1","2017Jusakual / platform 1",display);display<-sub("2017Jusakual2","2017Jusakual / platform 2",display)
  d[,Label:=sprintf("%s  (n = %d)",display,N)];d[,Label:=factor(Label,levels=rev(Label))]
  long<-melt(d,id.vars=c("Label","N"),measure.vars=c("C1","C2","C3"),variable.name="Subtype",value.name="Count");long[,Fraction:=Count/N]
  p<-ggplot(long,aes(Fraction,Label,fill=Subtype))+geom_col(width=.62,position=position_stack(reverse=TRUE))+geom_text(aes(label=Count,group=Subtype,colour=ifelse(Subtype=="C1","white",FIG5_INK)),position=position_stack(vjust=.5,reverse=TRUE),size=2.65)+scale_colour_identity()+scale_fill_manual(values=FIG5_SUBTYPE)+scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1),labels=scales::label_percent(accuracy=1),expand=c(0,0))+labs(title="Projected subtype composition",x="Projected fraction",y=NULL,fill=NULL)+fig5_theme()+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.title=element_text(size=8),legend.position="top")
  cover<-ggplot(d,aes(y=Label))+geom_col(aes(x=TotalFeatures),width=.62,fill="#E2E8ED")+geom_col(aes(x=MatchedFeatures),width=.62,fill=FIG5_INK)+geom_text(aes(x=57,label=paste0(MatchedFeatures," / ",TotalFeatures)),colour="white",size=2.7)+scale_x_continuous(limits=c(0,114),breaks=c(0,57,114),expand=c(0,0))+labs(title="Feature coverage",x="Measured panel genes",y=NULL)+fig5_theme()+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),plot.title=element_text(size=8))
  margin<-ggplot(d,aes(MedianMargin,Label))+geom_segment(aes(x=0,xend=MedianMargin,yend=Label),colour="grey80",linewidth=.4)+geom_point(size=2,colour=FIG5_INK)+geom_text(aes(label=sprintf("%.2f",MedianMargin)),hjust=1.35,size=2.7)+scale_x_continuous(limits=c(0,1.03),breaks=c(0,.5,1),expand=c(0,0))+labs(title="Assignment margin",x="Median probability gap",y=NULL)+fig5_theme()+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),plot.title=element_text(size=8))
  (p+cover+margin+plot_layout(widths=c(2.6,1.05,1.05)))+plot_annotation(title="Study-wise external microarray projection",caption="506 samples; 6 study/platform objects. Projection only; probability gaps do not estimate external accuracy.",theme=fig5_theme())
}
fig5_prepare<-function(p){if(is.function(p))return(p);wrap<-function(s)if(is.character(s))paste(strwrap(s,width=110),collapse="\n") else s;if(inherits(p,"patchwork"))p$patches$annotation$caption<-wrap(p$patches$annotation$caption) else p$labels$caption<-wrap(p$labels$caption);p}
build_fig5_consensus<-function(id){fig5_check();fig5_prepare(switch(id,A=fig5_schematic(),B=fig5_confusion(),C=fig5_comparison(),D=fig5_calibration(),E=fig5_loso(),F=fig5_external(),stop("Unknown Figure5 panel")))}
FIG5_NAMES<-c(A="Figure5A_training_schematic",B="Figure5B_oof_confusion",C="Figure5C_model_comparison",D="Figure5D_calibration",E="Figure5E_loso_forest",F="Figure5F_external_projection")
FIG5_HEIGHT<-c(A=68,B=88,C=76,D=96,E=84,F=96)
fig5_supplement<-function(id){
  fig5_check()
  if(id=="S11A"){
    d<-fig5_read("calibration_metrics")[Model==FIG5_PRIMARY];d<-melt(d,id.vars="Class",measure.vars=c("CalibrationIntercept","CalibrationSlope"),variable.name="Diagnostic",value.name="Value");d[,Diagnostic:=factor(Diagnostic,levels=c("CalibrationIntercept","CalibrationSlope"),labels=c("Calibration intercept","Calibration slope"))]
    refs<-data.table(Diagnostic=factor(c("Calibration intercept","Calibration slope"),levels=levels(d$Diagnostic)),Reference=c(0,1))
    return(ggplot(d,aes(Value,factor(Class,levels=c("C3","C2","C1")),colour=Class))+geom_vline(data=refs,aes(xintercept=Reference),inherit.aes=FALSE,linetype=2,colour="grey60",linewidth=.35)+geom_point(size=2.6)+geom_text(aes(label=sprintf("%.2f",Value)),nudge_y=.18,size=2.7)+scale_colour_manual(values=FIG5_SUBTYPE,guide="none")+facet_wrap(~Diagnostic,scales="free_x")+scale_x_continuous(expand=expansion(mult=c(.2,.2)))+labs(title="Calibration diagnostics for the primary model",x="OOF estimate",y=NULL,caption="Point estimates only; C1/C3 fits reported near-boundary probability warnings. Ideal intercept = 0; ideal slope = 1.")+fig5_theme())
  }
  if(id=="S11B"){
    d<-fig5_read("feature_selection_stability");not_selected<-d[Consensus==0,.N];ever<-d[Consensus>0];dist<-ever[,.N,by=Consensus];setorder(dist,Consensus);sizes<-fig5_read("outer_feature_sizes")[Method=="Consensus"]
    p<-ggplot(dist,aes(Consensus,N))+geom_col(width=.8,fill=FIG5_MID)+scale_x_continuous(limits=c(.5,25.5),breaks=c(1,5,10,15,20,25))+labs(title="Recurrence of consensus genes",x="Outer models selecting a gene (of 25)",y="Genes")+fig5_theme()+theme(plot.title=element_text(size=8))
    sizes[,FoldIndex:=seq_len(.N)];q<-ggplot(sizes,aes(FoldIndex,N))+geom_point(size=1.8,colour=FIG5_INK)+scale_x_continuous(breaks=c(1,5,10,15,20,25))+labs(title="Consensus panel size",x="Outer fit (repeat/fold order)",y="Selected genes")+fig5_theme()+theme(plot.title=element_text(size=8))
    return((p+q)+plot_annotation(title="Fold-wise feature-selection stability",caption=sprintf("Recurrence includes %d genes selected at least once; %d never-selected genes are not shown.",nrow(ever),not_selected),theme=fig5_theme()))
  }
  if(id=="S11C")return(fig5_comparison(TRUE))
  if(id=="S11D"){
    d<-fig5_read("external_sample_predictions")[Model==FIG5_PRIMARY];d[,Cohort:=sub("log2dat_","",Cohort)]
    return(ggplot(d,aes(Margin))+stat_ecdf(geom="step",colour=FIG5_INK,linewidth=.6)+facet_wrap(~Cohort,ncol=3)+scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1))+scale_y_continuous(breaks=c(0,.5,1))+labs(title="External projection probability gaps",x="Largest minus second-largest probability",y="Cumulative sample fraction",caption="All 506 microarray samples; continuous uncalibrated output, without a reliability cutoff or external accuracy interpretation.")+fig5_theme())
  }
  if(id=="S11E"){
    all<-fig5_read("nested_oof_predictions")[Model==FIG5_PRIMARY];rep<-fig5_read("nested_predictions_by_repeat")[Model==FIG5_PRIMARY];out<-list();repout<-list()
    coords<-function(z,cl){a<-z[,.(Score=get(paste0("Probability",cl)),Positive=as.integer(Truth==cl))][,.(Positive=sum(Positive),Negative=sum(1-Positive)),by=Score];setorder(a,-Score);a[,TP:=cumsum(Positive)];a[,FP:=cumsum(Negative)];a[,.(Recall=TP/sum(Positive),FPR=FP/sum(Negative),Precision=TP/(TP+FP))]}
    for(cl in c("C1","C2","C3")){a<-coords(all,cl);out[[cl]]<-rbind(data.table(Class=cl,Type="ROC",X=c(0,a$FPR),Y=c(0,a$Recall)),data.table(Class=cl,Type="Precision–recall",X=a$Recall,Y=a$Precision));repout[[cl]]<-rbindlist(lapply(sort(unique(rep$Repeat)),function(r){a<-coords(rep[Repeat==r],cl);data.table(Class=cl,Repeat=r,X=c(0,a$FPR),Y=c(0,a$Recall),Type="ROC")}))}
    d<-rbindlist(out);b<-rbindlist(repout);d[,Type:=factor(Type,levels=c("ROC","Precision–recall"))];b[,Type:=factor(Type,levels=levels(d$Type))]
    return(ggplot(d,aes(X,Y,colour=Class))+geom_step(data=b,aes(group=Repeat),colour="grey78",linewidth=.3)+geom_step(linewidth=.65)+facet_grid(Type~Class)+scale_colour_manual(values=FIG5_SUBTYPE,guide="none")+scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1))+scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1))+labs(title="OOF discrimination curves",x="False-positive rate (ROC) / recall (PR)",y="Sensitivity (ROC) / precision (PR)",caption="Colored curves use mean per-patient OOF probabilities; gray ROC curves show the five repeats. Coordinates are display aggregations.")+fig5_theme())
  }
  stop("Unknown classifier supplement")
}
render_fig5_consensus<-function(){
  fig5_check();out<-file.path(STEP8_ROOT,FIG5_REVIEW);for(d in c("figures","review","source_data"))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
  for(id in names(FIG5_NAMES)){
    p<-build_fig5_consensus(id);draw<-function(){if(is.function(p))p() else print(p)}
    file<-file.path(out,"figures",FIG5_NAMES[[id]])
    cairo_pdf(paste0(file,".pdf"),width=FIG5_WIDTH_MM/25.4,height=FIG5_HEIGHT[[id]]/25.4,family=STEP8_FONT$family);draw();dev.off()
    svglite(paste0(file,".svg"),width=174/25.4,height=FIG5_HEIGHT[[id]]/25.4);draw();dev.off()
    ragg::agg_png(file.path(out,"review",paste0(FIG5_NAMES[[id]],".png")),width=174,height=FIG5_HEIGHT[[id]],units="mm",res=300,background="white");draw();dev.off()
    cat("Rendered",id,"\n")
  }
  # Keep all four-model estimates inspectable; limma-123 is a historical-size baseline.
  for(id in paste0("S11",LETTERS[1:5])){
    p<-fig5_prepare(fig5_supplement(id));h<-if(id=="S11E")115 else if(id=="S11D")100 else 80;file<-file.path(out,"figures",paste0("Supplement_",id))
    cairo_pdf(paste0(file,".pdf"),width=174/25.4,height=h/25.4,family=STEP8_FONT$family);print(p);dev.off()
    svglite(paste0(file,".svg"),width=174/25.4,height=h/25.4);print(p);dev.off()
    ragg::agg_png(file.path(out,"review",paste0("Supplement_",id,".png")),width=174,height=h,units="mm",res=300,background="white");print(p);dev.off();cat("Rendered",id,"\n")
  }
  for(n in c("confusion_matrix","nested_performance_ci","nested_paired_model_differences","nested_oof_predictions","nested_predictions_by_repeat","calibration_bins","calibration_metrics","loso_per_study_performance","external_cohort_summary","external_coverage_qc","external_sample_predictions","outer_feature_sizes","feature_selection_stability","final_model_features","final_feature_method_support"))file.copy(file.path(STEP8_ROOT,FIG5_RUN,"tables",paste0(n,".csv")),file.path(out,"source_data",paste0(n,".csv")),overwrite=TRUE)
  invisible(out)
}

# Dedicated notebook/root entry: do not route these panels through the legacy
# 123-gene renderer. All other figures keep their existing rendering functions.
render_fig5_panel<-function(id){
  fig5_check();supp<-startsWith(id,"S11")
  name<-if(supp)paste0("Supplement_",id) else FIG5_NAMES[[id]]
  height<-if(supp){if(id=="S11E")115 else if(id=="S11D")100 else 80} else FIG5_HEIGHT[[id]]
  p<-if(supp)fig5_prepare(fig5_supplement(id)) else build_fig5_consensus(id)
  out<-file.path(STEP8_ROOT,FIG5_REVIEW);file<-file.path(out,"figures",name)
  draw<-function(){if(is.function(p))p() else print(p)}
  cairo_pdf(paste0(file,".pdf"),width=174/25.4,height=height/25.4,family=STEP8_FONT$family);draw();dev.off()
  svglite(paste0(file,".svg"),width=174/25.4,height=height/25.4);draw();dev.off()
  ragg::agg_png(file.path(out,"review",paste0(name,".png")),width=174,height=height,units="mm",res=300,background="white");draw();dev.off()
  manifest<-fread(file.path(STEP8_ROOT,"8-Visualization/tables/figure_manifest.csv"),colClasses="character",na.strings=NULL)
  ix<-which(manifest$Figure==(if(supp)"Supplement" else "Figure 5") & manifest$Panel==id)
  if(length(ix)!=1)stop("Missing or duplicate Figure5 contract")
  target<-file.path(STEP8_ROOT,manifest$TargetOutput[ix]);dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
  stopifnot(file.copy(paste0(file,".pdf"),target,overwrite=TRUE));stopifnot(file.copy(paste0(file,".svg"),sub("[.]pdf$",".svg",target),overwrite=TRUE))
  src<-strsplit(manifest$SourceTable[ix],";",fixed=TRUE)[[1]];stopifnot(all(file.exists(file.path(STEP8_ROOT,src))))
  fingerprint<-paste(paste0(src,"=",unname(tools::md5sum(file.path(STEP8_ROOT,src)))),collapse=";")
  manifest[ix,`:=`(Status="Rendered",WidthMM="174",HeightMM=as.character(height),SourceMD5=fingerprint,QAFile="",RenderedAt=format(Sys.time(),tz="UTC",usetz=TRUE))]
  tmp<-tempfile(tmpdir=file.path(STEP8_ROOT,"8-Visualization/tables"),fileext=".csv");fwrite(manifest,tmp,na="");stopifnot(file.rename(tmp,file.path(STEP8_ROOT,"8-Visualization/tables/figure_manifest.csv")))
  message("Rendered selected consensus version: ",name)
  invisible(target)
}
