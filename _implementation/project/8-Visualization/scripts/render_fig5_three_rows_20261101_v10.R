# Three-row layout, same source statistics, final physical widths.
if(!exists('fig5_v9_current_panel',inherits=FALSE)){
 source('8-Visualization/scripts/render_fig5_final_20261031_v9.R');fig5_v9_current_panel<-render_fig5_current_panel
}
V10<-'8-Visualization/refinements/fig5_three_rows_20261101_v10';v10_path<-function(...)file.path(STEP8_ROOT,V10,...)
completion_path<-function(...)v10_path(...)
v10_theme<-function(){v9_theme()+theme(plot.title=element_text(size=7,face='bold'),plot.subtitle=element_text(size=6),axis.title=element_text(size=6),axis.text=element_text(size=5.8),legend.text=element_text(size=5.8),legend.title=element_text(size=6),plot.caption=element_text(size=5.8),plot.margin=margin(1.5,1.5,1.5,1.5,'mm'))}
vt<-function(t,x,y,size=6,bold=FALSE,just='left',col=FIG5_INK,italic=FALSE){grid.text(t,unit(x,'native'),unit(y,'native'),just=just,gp=gpar(fontfamily='Arial',fontsize=size,fontface=if(italic)'italic' else if(bold)'bold' else 'plain',col=col,lineheight=1.05))}
vr<-function(x,y,w,h,col,edge=NA){grid.rect(unit(x,'native'),unit(y,'native'),unit(w,'native'),unit(h,'native'),gp=gpar(fill=col,col=edge,lwd=.35))}
vl<-function(x,y,col=FIG5_INK,lwd=.5,arrow=FALSE){grid.lines(unit(x,'native'),unit(y,'native'),gp=gpar(col=col,lwd=lwd),arrow=if(arrow)grid::arrow(length=unit(1,'mm'),type='closed') else NULL)}
vp<-function(x,y,col=FIG5_INK,size=.8,pch=16){grid.points(unit(x,'native'),unit(y,'native'),pch=pch,size=unit(size,'mm'),gp=gpar(col=col))}
canvas<-function(w,h){grid.newpage();pushViewport(viewport(xscale=c(0,w),yscale=c(0,h)))}
v10_A<-function(){function(){canvas(174,38)
 starts<-c(2,32,77,108,140);ww<-c(26,41,27,28,32);heads<-c('Data','Feature selection','Consensus','Model','Evaluation')
 for(j in 1:5){x<-starts[j];w<-ww[j];vr(x+w/2,20,w,28,'#F6F8FA','#DCE4EA');vr(x+w/2,32,w,5,'#DCE5EC');vt(heads[j],x+w/2,32,7,TRUE,'centre');if(j<5)vl(c(x+w+.3,starts[j+1]-.5),c(20,20),arrow=TRUE)}
 vt('497 principal tumors',15,26,6.3,FALSE,'centre');vt('Expression + labels',15,22,6,FALSE,'centre')
 for(j in 1:3){vr(7+(j-1)*8,16,7.5,3,FIG5_SUBTYPE[j]);vt(names(FIG5_SUBTYPE)[j],7+(j-1)*8,12,6,TRUE,'centre')}
 vt('limma',34,26,6.5,TRUE);vt('Study-adjusted DE',48,26,6)
 vt('RF',34,19,6.5,TRUE);vt('Permutation importance',48,19,6)
 vt('LASSO',34,12,6.5,TRUE);vt('Multinomial selection',48,12,6)
 vt('≥2 methods',90.5,25,6.5,FALSE,'centre');vt('114 genes',90.5,17,8,TRUE,'centre');vt('Full-development panel',90.5,10,5.7,FALSE,'centre')
 vt('Linear SVM',122,24,7.5,TRUE,'centre');vt('Primary classifier',122,18,6,FALSE,'centre');vt('RF / kNN comparators',122,11,5.7,FALSE,'centre')
 vt('Repeated nested CV',156,26,6,FALSE,'centre');vt('Study-held-out transfer',156,19,6,FALSE,'centre');vt('Microarray projection',156,12,6,FALSE,'centre')
 vt('Training-fold preprocessing, selection and tuning  ·  Outer 5-fold × 5 repeats; inner 4-fold',87,2.5,6.3,FALSE,'centre');popViewport()
}}
v10_B<-function(){function(){canvas(72,76);d<-fread(v8_path('source_data/B_subtype_gene_order.csv'));vt('Consensus gene evidence',1,73.5,7.5,TRUE)
 x0<-11;cw<-58/114;d[,X:=x0+(DisplayOrder-.5)*cw+(SubtypeGroup-1)*.8];groups<-d[,.(N=.N,L=min(X)-cw/2,R=max(X)+cw/2),by=SubtypeGroup]
 for(j in 1:3){g<-groups[SubtypeGroup==j];vr((g$L+g$R)/2,66.5,g$R-g$L,5,adjustcolor(FIG5_SUBTYPE[j],alpha.f=.15));vt(sprintf('C%d | %d',j,g$N),(g$L+g$R)/2,66.5,6.4,TRUE,'centre')}
 vt('Method support',x0,60,6.4,TRUE)
 for(j in 1:3){y<-56-(j-1)*3;vt(c('limma','RF','LASSO')[j],x0-1,y,6,FALSE,'right');for(i in 1:114)vr(d$X[i],y,cw-.03,2.4,if(d[[c('LimmaFull','RFFull','LASSOFull')[j]]][i])FIG5_INK else '#E7EDF1')}
 vt('Selection frequency (25 fits)',x0,44.5,6.4,TRUE)
 for(v in c(0,.5,1)){y<-26+14*v;vl(c(x0,71),c(y,y),'#E3E8EC');vt(sprintf('%d%%',100*v),x0-1,y,5.8,FALSE,'right')}
 for(i in 1:114){y<-26+14*d$ConsensusFraction[i];vl(rep(d$X[i],2),c(26,y),'#A7BCCB',.4);vp(d$X[i],y,size=.45)}
 vt('Mean expression',x0,21.5,6.4,TRUE);ramp<-circlize::colorRamp2(c(-1.5,0,1.5),STEP8_CONFIG$effect_colors)
 for(j in 1:3){y<-17.5-(j-1)*3;vt(paste0('C',j),x0-1,y,6,TRUE,'right',FIG5_SUBTYPE[j]);for(i in 1:114)vr(d$X[i],y,cw-.03,2.4,ramp(d[[paste0('C',j,'MeanTrainingZ')]][i]))}
 vr(2,4,1.3,1.3,FIG5_INK);vt('Selected',3.5,4,5.5);vr(16,4,1.3,1.3,'#E7EDF1');vt('Not selected',17.5,4,5.5)
 vt('Mean z',39,4,5.5);vt('−1.5',51,4,5.5,FALSE,'right');for(i in 1:30)vr(52+(i-.5)*.42,4,.42,1.8,ramp(-1.5+(i-.5)*3/30));vt('1.5',65,4,5.5);popViewport()
}}
v10_C<-function(){function(){canvas(98,76);d<-fread(file.path(STEP8_ROOT,BC6,'source_data/C_subtype_ranking_all114.csv'));vt('Subtype-representative genes',1,73.5,7.5,TRUE)
 for(j in 1:3){cl<-names(FIG5_SUBTYPE)[j];z<-d[Class==cl];setorder(z,WithinClassRank);top<-z[Labelled==TRUE];x0<-8+(j-1)*31;x<-function(r)x0+(r-1)/113*26;y<-function(v)30+(v+2.5)/5*30
 vt(cl,x0,68.5,7.3,TRUE,col=FIG5_SUBTYPE[j]);vt(sprintf('%d with Δz > 0',sum(z$HigherThanBoth)),x0,64.8,5.8)
 vl(c(x0,x0+26),rep(y(0),2),'#BAC5CE',.4);vl(c(x0,x0),c(30,60),'#7D8A92',.4);vl(c(x0,x0+26),c(30,30),'#7D8A92',.4)
 vl(x(z$WithinClassRank),y(z$ExpressionMargin),'#D0D8DF',.65);vp(x(z$WithinClassRank),y(z$ExpressionMargin),'#D0D8DF',.55);pos<-z[HigherThanBoth==TRUE];vp(x(pos$WithinClassRank),y(pos$ExpressionMargin),FIG5_SUBTYPE[j],.75)
 for(v in c(-2,0,2)){vl(c(x0-.5,x0),rep(y(v),2),'#7D8A92',.4);if(j==1)vt(v,x0-1,y(v),5.5,FALSE,'right')}
 ly<-seq(59,44,length.out=6);for(i in 1:6){vl(c(x(top$WithinClassRank[i]),x0+12),c(y(top$ExpressionMargin[i]),ly[i]),'#7D91A0',.35);vt(top$Gene[i],x0+12.3,ly[i],5.8,italic=TRUE)}
 for(r in c(1,60,114))vt(r,x(r),27.8,5.5,FALSE,'centre');vt('Rank',x0+13,25,5.8,FALSE,'centre')
 t0<-1+(j-1)*32.3;xs<-t0+c(0,14,20,25,30);heads<-c('Gene','Δz','Sel.%','LRA','RF')
 for(k in 1:5)vt(heads[k],xs[k],21,5.5,TRUE,if(k==1)'left' else 'centre');vl(c(t0,t0+30.7),c(19.5,19.5),'#CCD6DD',.4)
 for(i in 1:6){yy<-17.8-(i-1)*2.9;if(i%%2)vr(t0+15.3,yy,30.7,2.7,'#F1F5F7');vt(top$Gene[i],xs[1],yy,5.5);vt(sprintf('%.2f',top$ExpressionMargin[i]),xs[2],yy,5.5,FALSE,'centre');vt(round(100*top$ConsensusFraction[i]),xs[3],yy,5.5,FALSE,'centre');for(k in 1:3)vr(xs[4]+(k-2)*1.25,yy,.95,.95,if(top[[c('LimmaFull','RFFull','LASSOFull')[k]]][i])FIG5_INK else '#DCE4EA');vt(top$RFRank[i],xs[5],yy,5.5,FALSE,'centre')}
 }
 grid.text('Expression margin (Δz)',unit(2,'native'),unit(45,'native'),rot=90,gp=gpar(fontfamily='Arial',fontsize=5.5,col=FIG5_INK));popViewport()
}}
v10_D<-function(){function(){canvas(60,50);cm<-fig5_read('confusion_matrix')[Model==FIG5_PRIMARY];ramp<-circlize::colorRamp2(c(0,1),c('#F0F5F8','#4E8DB1'))
 vt('Internal OOF (n = 497)',1,47.5,7,TRUE);vt('Accuracy 98.79%  |  BA 98.79%',1,43.5,6)
 for(i in 1:3)for(j in 1:3){z<-cm[Truth==paste0('C',i)&Prediction==paste0('C',j)];x<-18+(j-.5)*10;y<-10+(3-i+.5)*10;vr(x,y,9.7,9.7,ramp(z$RowFraction));vt(sprintf('%d\n%.1f%%',z$Freq,100*z$RowFraction),x,y,6.4,FALSE,'centre',if(z$RowFraction>.5)'white' else FIG5_INK)}
 for(j in 1:3){vt(paste0('C',j),18+(j-.5)*10,7.5,6,FALSE,'centre');vt(paste0('C',j),16,10+(3-j+.5)*10,6,FALSE,'right')};vt('Predicted',33,2.8,6,FALSE,'centre');grid.text('Reference',unit(8,'native'),unit(25,'native'),rot=90,gp=gpar(fontfamily='Arial',fontsize=6,col=FIG5_INK));popViewport()
}}
v10_E<-function(){function(){canvas(60,64);s<-fread(v8_path('source_data/microarray_summary.csv'));d<-fread(v8_path('source_data/microarray_composition.csv'));short<-c('Andersen','Sia','Jusakual P1','Jusakual P2','Job','TIGER')
 vt('Microarray projection (n = 506)',1,61.7,7,TRUE);vt('Study (n) · genes / 114',1,57.4,5.7);vt('Prob. gap',50,57.4,5.7,TRUE,'centre')
 for(i in 1:6){yy<-53.5-(i-1)*8;z<-s[CohortOrder==i];q<-d[CohortOrder==i];setorder(q,Prediction);vt(sprintf('%s (%d) · %d/114',short[i],z$N,z$MatchedFeatures),1,yy,5.7)
  cum<-0;for(j in 1:3){a<-q[Prediction==paste0('C',j)];wid<-34*a$Fraction;vr(1+cum+wid/2,yy-3, wid,2.7,FIG5_SUBTYPE[j]);vt(a$N,1+cum+wid/2,yy-3,5.5,FALSE,'centre',if(j==2)FIG5_INK else 'white');cum<-cum+wid}
  vl(c(42,58),rep(yy-3,2),'#E1E7EC',.4);vl(42+16*c(z$Q1,z$Q3),rep(yy-3,2),FIG5_INK,1.3);vp(42+16*z$MedianGap,yy-3,size=1)
 }
 vt('0',1,7.6,5.5);vt('100%',35,7.6,5.5,FALSE,'right');vt('0',42,7.6,5.5,FALSE,'centre');vt('1',58,7.6,5.5,FALSE,'centre')
 for(j in 1:3){xx<-3+(j-1)*10;vr(xx,3.3,1.8,1.8,FIG5_SUBTYPE[j]);vt(paste0('C',j),xx+1.7,3.3,5.5)};vt('Median [IQR]',50,3.3,5.5,FALSE,'centre');popViewport()
}}
v10_F<-function(){
 curves<-opt_read('external_survival_curves');risk<-opt_read('external_survival_at_risk');tests<-opt_read('external_survival_tests')
 ps<-lapply(c('2017Jusakual','2020Job'),function(st){z<-curves[Study==st];rr<-risk[Study==st];tt<-tests[Study==st];xmax<-max(rr$Time);ticks<-seq(0,xmax,48);rr<-rr[Time%in%ticks]
 band<-rbindlist(lapply(split(z,z$Class),function(v){setorder(v,Time);b<-v[rep(seq_len(nrow(v)),each=2)];b[,Time:=as.vector(rbind(v$Time,c(v$Time[-1],tail(v$Time,1))))];b}))
 p<-ggplot(z,aes(Time,Survival,colour=Class))+geom_ribbon(data=band,aes(ymin=Lower95,ymax=Upper95,fill=Class),alpha=.08,colour=NA,show.legend=FALSE,na.rm=TRUE)+geom_step(linewidth=.4)+geom_point(data=z[Censored>0],shape=3,size=.8,stroke=.25,show.legend=FALSE)+scale_colour_manual(values=FIG5_SUBTYPE)+scale_fill_manual(values=FIG5_SUBTYPE)+scale_x_continuous(limits=c(0,xmax),breaks=ticks,expand=expansion(mult=c(.035,.045)))+scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1),expand=expansion(mult=c(.01,.02)))+labs(title=sprintf('%s (n=%d)',if(st=='2017Jusakual')'Jusakual' else 'Job',tt$N),subtitle=sprintf('Holm P = %.3f',tt$P_Holm),x='OS (months)',y='Survival',colour=NULL)+v10_theme()+theme(legend.position=if(st=='2017Jusakual')'top' else 'none',legend.key.width=unit(2,'mm'),legend.spacing.x=unit(1,'mm'),legend.key.height=unit(1.5,'mm'),plot.margin=margin(1,1,0,1,'mm'))
 rr[,Class:=factor(Class,levels=c('C3','C2','C1'))];q<-ggplot(rr,aes(Time,Class,label=AtRisk,colour=Class))+geom_text(size=1.95)+scale_colour_manual(values=FIG5_SUBTYPE,guide='none')+scale_x_continuous(limits=c(0,xmax),breaks=ticks,expand=expansion(mult=c(.035,.045)))+labs(title='At risk',x=NULL,y=NULL)+v10_theme()+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(margin=margin(r=1.5,unit='mm')),plot.title=element_text(size=5.8),plot.margin=margin(0,1,1,1,'mm'))+coord_cartesian(clip='off')
 p/q+plot_layout(heights=c(3.8,1))
 });wrap_plots(ps,ncol=1)+plot_annotation(title='External survival',theme=v10_theme())
}
v10_G<-function(){
 w<-mn('full_candidate_weights');w[,Gene:=factor(Gene,levels=rev(Gene))]
 p<-ggplot(w,aes(100*RelativeWeight,Gene))+geom_col(fill=FIG5_INK,width=.65)+labs(title='17-gene candidate',subtitle='Full-development weights',x='Relative weight (%)',y=NULL)+v10_theme()+theme(axis.text.y=element_text(face='italic',size=5.8),plot.margin=margin(1,1.5,.5,1.5,'mm'))
 d<-mn('performance_ci')[Metric%in%c('BalancedAccuracy','C1Recall','C2Recall','C3Recall')];d[,Y:=match(Metric,c('C3Recall','C2Recall','C1Recall','BalancedAccuracy'))+ifelse(Model=='Full',.13,-.13)]
 q<-ggplot(d,aes(100*Estimate,Y,shape=Model,colour=Model))+geom_errorbar(aes(xmin=100*Lower95,xmax=100*Upper95),orientation='y',width=.12,linewidth=.3)+geom_point(size=1.3)+scale_y_continuous(breaks=1:4,labels=c('C3','C2','C1','BA'),limits=c(.5,4.5))+scale_x_continuous(breaks=c(94,97,100),limits=c(93,100.4))+scale_colour_manual(values=c(Full='#8197A6',Compact=FIG5_INK),labels=c(Full='Full',Compact='Inner-selected'))+scale_shape_manual(values=c(Full=1,Compact=16),labels=c(Full='Full',Compact='Inner-selected'))+labs(title='Outer OOF',subtitle='Selected K: 23 [10–51]',x='Performance (%)',y=NULL,shape=NULL,colour=NULL)+v10_theme()+theme(legend.position='bottom',legend.key.width=unit(2,'mm'),legend.key.height=unit(1.5,'mm'),legend.spacing.x=unit(1,'mm'),axis.line.y=element_blank(),axis.ticks.y=element_blank(),legend.margin=margin(0,0,0,0),legend.box.spacing=unit(.3,'mm'),plot.margin=margin(.5,1.5,.5,1.5,'mm'))
 (p/q)+plot_layout(heights=c(1.45,1))+plot_annotation(title='Feature reduction',caption='BA loss: 1.16 pp > 1 pp target',theme=v10_theme())
}
v10_SG<-function(){r<-fig5_read('roc_coordinates')[Model==FIG5_PRIMARY];r[,FPR:=1-Specificity];m<-fig5_read('nested_performance_ci')[Model==FIG5_PRIMARY];labels<-setNames(vapply(names(FIG5_SUBTYPE),function(cl)sprintf('%s: %.4f',cl,m[Metric==paste0(cl,'_AUROC'),Estimate]),character(1)),names(FIG5_SUBTYPE));ggplot(r,aes(FPR,Sensitivity,colour=Class))+geom_abline(slope=1,intercept=0,colour='grey80',linetype=2,linewidth=.3)+geom_step(linewidth=.5)+scale_colour_manual(values=FIG5_SUBTYPE,labels=labels)+scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1))+scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1))+coord_equal()+labs(title='Internal OOF ROC (n = 497)',x='False-positive rate',y='True-positive rate',colour='AUROC')+v9_theme()+theme(legend.position='inside',legend.position.inside=c(.7,.23),legend.background=element_rect(fill='white',colour=NA))}
v10_SH<-function(){x<-fread(v8_path('source_data/microarray_samples.csv'));s<-fread(v8_path('source_data/microarray_summary.csv'));setorder(x,CohortOrder,Sample);x[,Y:=7-CohortOrder+((seq_len(.N)-1)%%9-4)*.035,by=Cohort];s[,Y:=7-CohortOrder];ggplot(x,aes(Margin,Y))+geom_point(size=.8,alpha=.5,colour='#9DAFBC')+geom_segment(data=s,aes(x=Q1,xend=Q3,y=Y,yend=Y),inherit.aes=FALSE,linewidth=1,colour=FIG5_INK)+geom_point(data=s,aes(MedianGap,Y),inherit.aes=FALSE,size=1.6,colour=FIG5_INK)+scale_y_continuous(breaks=s$Y,labels=paste0(c('Andersen','Sia','Jusakual P1','Jusakual P2','Job','TIGER'),' (',s$N,')'))+scale_x_continuous(limits=c(0,1.02),breaks=c(0,.5,1))+labs(title='Microarray probability gaps',subtitle='All 506 samples; median and IQR',x='Largest − second-largest probability',y=NULL)+v9_theme()+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())}
V10_NAMES<-c(A='Figure5A_workflow',B='Figure5B_consensus_evidence',C='Figure5C_subtype_genes',D='Figure5D_internal_performance',E='Figure5E_microarray_projection',F='Figure5F_external_survival',G='Figure5G_feature_reduction',V9_NAMES[c('SA','SB','SC','SD','SE','SF')],SG='SupplementaryFigure5G_internal_ROC',SH='SupplementaryFigure5H_probability_gaps')
V10_DIMS<-c(list(A=c(174,38),B=c(72,76),C=c(98,76),D=c(60,50),E=c(60,64),F=c(53,118),G=c(53,106)),V9_DIMS[c('SA','SB','SC','SD','SE','SF')],list(SG=c(90,80),SH=c(90,76)))
render_fig5_v10_panel<-function(id){fig5_check();stem<-v10_path('figures',V10_NAMES[[id]]);sz<-V10_DIMS[[id]]
 if(id%in%c('SA','SB','SC','SD','SE','SF')){for(ext in c('pdf','svg','png'))stopifnot(file.copy(v9_path('figures',paste0(V9_NAMES[[id]],'.',ext)),paste0(stem,'.',ext),overwrite=TRUE))}else{
 p<-switch(id,A=v10_A(),B=v10_B(),C=v10_C(),D=v10_D(),E=v10_E(),F=v10_F(),G=v10_G(),SG=v10_SG(),SH=v10_SH());draw<-function(){if(is.function(p))p() else print(p)}
 cairo_pdf(paste0(stem,'.pdf'),width=sz[1]/25.4,height=sz[2]/25.4,family='Arial');draw();dev.off();svglite::svglite(paste0(stem,'.svg'),width=sz[1]/25.4,height=sz[2]/25.4);draw();dev.off();ragg::agg_png(paste0(stem,'.png'),width=sz[1],height=sz[2],units='mm',res=300,background='white');draw();dev.off()
 }
 stopifnot(system2('pdftoppm',c('-png','-singlefile','-r','150',shQuote(paste0(stem,'.pdf')),shQuote(v10_path('review/pdf_rendered',V10_NAMES[[id]]))))==0);cat('Rendered',id,'\n');invisible(paste0(stem,'.pdf'))
}

V10_IDS<-c(A='A',B='B',C='C',D='D',E='E',F='F',G='G',S11N='SG',S11O='SH')
render_fig5_current_panel<-function(id){
 if(!id%in%names(V10_IDS))return(fig5_v9_current_panel(id))
 mp<-file.path(STEP8_ROOT,'8-Visualization/tables/figure_manifest.csv');m<-fread(mp,colClasses='character',na.strings=NULL)
 ix<-which(m$Figure==if(startsWith(id,'S11'))'Supplement' else 'Figure 5');ix<-ix[m$Panel[ix]==id];stopifnot(length(ix)==1)
 target<-file.path(STEP8_ROOT,m$TargetOutput[ix]);src<-render_fig5_v10_panel(V10_IDS[[id]])
 for(ext in c('pdf','svg','png'))stopifnot(file.copy(sub('[.]pdf$',paste0('.',ext),src),sub('[.]pdf$',paste0('.',ext),target),overwrite=TRUE))
 sources<-strsplit(m$SourceTable[ix],';',fixed=TRUE)[[1]];stopifnot(all(file.exists(file.path(STEP8_ROOT,sources))))
 m[ix,`:=`(Status='Rendered',SourceMD5=paste(paste0(sources,'=',unname(tools::md5sum(file.path(STEP8_ROOT,sources)))),collapse=';'),QAFile='',RenderedAt=format(Sys.time(),tz='UTC',usetz=TRUE))]
 tmp<-tempfile(tmpdir=dirname(mp),fileext='.csv');fwrite(m,tmp,na='');stopifnot(file.rename(tmp,mp));invisible(target)
}
# Raster arrangement is a review preview, not the publication master.
v10_layout_preview<-function(){
 stopifnot(requireNamespace('png',quietly=TRUE));ragg::agg_png(v10_path('review/Figure5_three_row_preview.png'),width=182,height=260,units='mm',res=300,background='white');canvas(182,260)
 place<-data.table(ID=c('A','B','C','D','E','F','G'),X=c(4,4,80,4,4,68,125),Top=c(8,54,54,138,192,138,138))
 for(i in 1:nrow(place)){id<-place$ID[i];sz<-V10_DIMS[[id]];im<-png::readPNG(v10_path('figures',paste0(V10_NAMES[[id]],'.png')));grid.raster(im,x=unit(place$X[i],'native'),y=unit(260-place$Top[i],'native'),width=unit(sz[1],'native'),height=unit(sz[2],'native'),just=c('left','top'),interpolate=TRUE);vt(id,place$X[i]-3,260-place$Top[i]-2,8,TRUE)}
 popViewport();dev.off();fwrite(place,v10_path('source_data/layout_positions_mm.csv'))
}
