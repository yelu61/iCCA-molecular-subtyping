#!/usr/bin/env Rscript
# Layout-only reorganization. Existing eligibility and sensitivity membership are unchanged.
options(stringsAsFactors=FALSE,warn=2)
library(grid)
root<-normalizePath(getwd())
base<-file.path(root,"8-Visualization/refinements/fig2_supp_reorganization_20261103_v3/tf")
v2<-file.path(root,"8-Visualization/refinements/fig2_global_revision_20261103_v2/tf_program_bridge")
v33<-file.path(root,"8-Visualization/refinements/fig2D_selection_20261031_v33")
for(p in list.files(file.path(v2,"source_data"),full.names=TRUE))file.copy(p,file.path(base,"source_data/v2_preserved",basename(p)),overwrite=TRUE)
for(p in list.files(file.path(v33,"tables"),full.names=TRUE))file.copy(p,file.path(base,"source_data/v33_preserved",basename(p)),overwrite=TRUE)
file.copy(file.path(v33,"SELECTION_CONTRACT.json"),file.path(base,"source_data/v33_preserved/SELECTION_CONTRACT.json"),overwrite=TRUE)
eg<-read.csv(file.path(v33,"tables/D_all263_eligibility.csv"),check.names=FALSE)
sel<-read.csv(file.path(v33,"tables/D_selection12.csv"),check.names=FALSE)
mem<-read.csv(file.path(v33,"tables/D_sensitivity_membership.csv"),check.names=FALSE)
ss<-read.csv(file.path(v33,"tables/D_sensitivity_summary.csv"),check.names=FALSE)
oldtext<-read.csv(file.path(v33,"tables/SFigure2D_selection_sensitivity_texts.csv"),check.names=FALSE)
for(f in c("CoveragePass","OriginalSignificant","LOSOStable","Eligible"))eg[[f]]<-tolower(as.character(eg[[f]]))=="true"
stopifnot(nrow(eg)==789L,length(unique(eg$TF))==263L,nrow(sel)==12L,nrow(mem)==120L,nrow(ss)==10L)
pred<-list(rep(TRUE,nrow(eg)),eg$OriginalSignificant,eg$OriginalSignificant&eg$CoveragePass,eg$Eligible)
flow<-data.frame(Stage=c("Tested","Original FDR < 0.05","+ >=10 measured signed targets","+ 7/7 LOSO sign retained"),
 TFs=sapply(pred,function(x)length(unique(eg$TF[x]))),Contrasts=sapply(pred,sum))
stopifnot(identical(as.integer(flow$TFs),c(263L,250L,225L,225L)),identical(as.integer(flow$Contrasts),c(789L,585L,526L,524L)))
queues<-expand.grid(Subtype=c("C1","C2","C3"),Direction=c("Positive","Negative"),stringsAsFactors=FALSE)
queues$Label<-paste0(queues$Subtype,ifelse(queues$Direction=="Positive","+","\u2212"))
queues$Eligible<-mapply(function(st,di)sum(eg$Subtype==st&eg$Direction==di&eg$Eligible),queues$Subtype,queues$Direction)
queues$Selected<-mapply(function(st,di)sum(sel$AnchorSubtype==st&sel$AnchorDirection==di),queues$Subtype,queues$Direction)
stopifnot(all(queues$Selected==2L),sum(is.na(eg$LOSOStudies))==39L,length(unique(eg$TF[is.na(eg$LOSOStudies)]))==13L)
tfrow<-oldtext[oldtext$Text%in%unique(mem$TF),];tfrow<-tfrow[order(-tfrow$Y),]
tfs<-tfrow$Text;stopifnot(length(tfs)==19L,!anyDuplicated(tfs),setequal(tfs,mem$TF))
ss<-ss[order(match(ss$Mode,c("signed","unsigned")),ss$Threshold,ss$Budget),]
ss$SettingID<-paste(ss$Mode,sprintf("%.2f",ss$Threshold),ss$Budget,sep="_");ss$ColumnOrder<-seq_len(nrow(ss))
mem$SettingID<-paste(mem$Mode,sprintf("%.2f",mem$Threshold),mem$Budget,sep="_")
cells<-expand.grid(TF=tfs,SettingID=ss$SettingID,stringsAsFactors=FALSE)
cells$Selected<-paste(cells$TF,cells$SettingID)%in%paste(mem$TF,mem$SettingID)
cells$RowOrder<-match(cells$TF,tfs);cells$ColumnOrder<-match(cells$SettingID,ss$SettingID)
cells$Main12<-cells$TF%in%sel$TF
for(f in c("Mode","Threshold","Budget"))cells[[f]]<-ss[[f]][match(cells$SettingID,ss$SettingID)]
stopifnot(nrow(cells)==190L,sum(cells$Selected)==120L,sum(!cells$Selected)==70L)
for(i in seq_len(nrow(ss)))stopifnot(sum(cells$Selected[cells$SettingID==ss$SettingID[i]])==ss$Selected[i])
write.csv(flow,file.path(base,"source_data/G_sequential_eligibility_flow.csv"),row.names=FALSE)
write.csv(queues,file.path(base,"source_data/G_six_anchor_queues.csv"),row.names=FALSE)
write.csv(cells,file.path(base,"source_data/I_complete_19x10_membership.csv"),row.names=FALSE)
write.csv(ss,file.path(base,"source_data/I_ten_settings_original_summary.csv"),row.names=FALSE)
# H is reused byte for byte; it is not rerendered or rescaled here.
for(ext in c("pdf","svg","png"))file.copy(file.path(v2,"figures",paste0("TF_program_correspondence.",ext)),
 file.path(base,"figures",paste0("H_TF_program_correspondence.",ext)),overwrite=TRUE)
file.copy(file.path(v2,"review/TF_program_correspondence_preview.png"),file.path(base,"review/H_TF_program_correspondence_preview.png"),overwrite=TRUE)
preserved<-c(list.files(file.path(v2,"source_data"),full.names=TRUE),list.files(file.path(v33,"tables"),full.names=TRUE),file.path(v33,"SELECTION_CONTRACT.json"))
write.csv(data.frame(Source=sub(paste0(root,"/"),"",preserved,fixed=TRUE),MD5=unname(tools::md5sum(preserved))),
 file.path(base,"review/preserved_input_hashes.csv"),row.names=FALSE)

family<-"Arial";ink<-"#26313A";soft<-"#626B74";rule<-"#D6DCE2";selected<-"#247AA4"
mm<-function(x)unit(x,"mm")
txt<-function(s,x,y,size=6.8,face="plain",col=ink,just="left")grid.text(s,mm(x),mm(y),just=just,
 gp=gpar(fontfamily=family,fontsize=size,fontface=face,col=col,lineheight=1.04))
seg<-function(x1,y1,x2,y2,col=rule,lwd=.45)grid.segments(mm(x1),mm(y1),mm(x2),mm(y2),gp=gpar(col=col,lwd=lwd))
rect<-function(x,y,w,h,fill,col=NA,lwd=.35)grid.rect(mm(x),mm(y),mm(w),mm(h),gp=gpar(fill=fill,col=col,lwd=lwd))
drawG<-function(){
 grid.newpage();pushViewport(viewport(x=0,y=0,width=mm(174),height=mm(37),just=c("left","bottom")));rect(87,18.5,174,37,"white")
 txt("TF eligibility and display budget",1.5,33.8,7.8,"bold")
 txt("Qualifying TFs / contrasts",1.5,29.5,6.5,col=soft)
 xs<-c(12.5,37,63,89)
 labs<-c("Tested","FDR < 0.05","\u226510 signed targets","7/7 LOSO sign")
 for(i in 1:4){txt(labs[i],xs[i],25.8,6.3,just="center");txt(paste0(flow$TFs[i]," / ",flow$Contrasts[i]),xs[i],21.7,7.3,"bold",just="center")}
 for(x in c(24.7,50,76))txt("\u2192",x,21.7,7.4,col=soft,just="center")
 seg(103,18.0,103,31.2)
 txt("Six anchor queues: selected / eligible",139,29.5,6.5,col=soft,just="center")
 qx<-seq(110.5,169,by=11.7)
 for(i in 1:6){txt(queues$Label[i],qx[i],25.8,6.5,"bold",just="center");txt(paste0(queues$Selected[i]," / ",queues$Eligible[i]),qx[i],21.7,6.8,just="center")}
 txt("The display budget is 12 TFs; this is not an additional significance filter. Main D retains all 36 subtype contrasts.",1.5,16.6,6.5,"bold")
 txt("Missing LOSO: 13 TFs / 39 contrasts, never counted as passing. Two significant contrasts with sufficient coverage fail 7/7 sign retention.",1.5,11.5,6.3,col=soft)
 txt("Original TF FDR: 263 factors separately within each subtype contrast. Coverage counts measured signed target genes, not samples.",1.5,7.2,6.3,col=soft)
 txt("LOSO direction retention applies to the selection anchor; it does not guarantee stability of every displayed counterpart contrast.",1.5,2.9,6.3,col=soft)
 popViewport()
}
drawI<-function(){
 grid.newpage();pushViewport(viewport(x=0,y=0,width=mm(174),height=mm(90),just=c("left","bottom")));rect(87,45,174,90,"white")
 txt("Sensitivity of displayed TF membership",1.5,86.7,7.8,"bold")
 cx<-c(29,42,55,73,86,99,117,130,143,165)
 for(i in 1:3){ii<-((i-1)*3+1):(i*3);x<-mean(cx[ii]);txt("Signed",x,82.0,6.5,"bold",just="center");txt(sprintf("%.2f",ss$Threshold[ii[1]]),x,78.7,6.5,just="center")}
 txt("Unsigned U",cx[10],82.0,6.3,"bold",just="center");txt("0.50",cx[10],78.7,6.5,just="center")
 txt("Mode",1.5,82.0,6.3,col=soft);txt("Jaccard cutoff",1.5,78.7,6.3,col=soft);txt("TF / budget",1.5,75.1,6.5,"bold")
 for(j in 1:10)txt(ss$Budget[j],cx[j],75.1,6.8,"bold",just="center")
 seg(1.5,73.0,172,73.0)
 yy<-71.0-(seq_along(tfs)-1)*2.9
 for(i in seq_along(tfs)){
  tf<-tfs[i];txt(tf,2,yy[i],6.8,if(tf%in%sel$TF)"bold"else"plain")
  for(j in 1:10){d<-cells[cells$RowOrder==i&cells$ColumnOrder==j,];rect(cx[j],yy[i],11.3,2.25,if(d$Selected)selected else "white",if(d$Selected)selected else rule,lwd=.35)}
 }
 seg(1.5,mean(yy[12:13]),172,mean(yy[12:13]),col="#BCC6CE",lwd=.5)
 rect(3.0,13.6,2.6,1.5,selected);txt("Selected",5.3,13.6,6.3)
 rect(28.8,13.6,2.6,1.5,"white",rule);txt("Not selected",31.2,13.6,6.3)
 txt("Bold TFs: main 12",64,13.6,6.3)
 txt("U: unsigned target Jaccard",112,13.6,6.3)
 txt("All 19 TFs selected in any of 10 settings are shown. The 12-TF membership is unchanged across tested cutoffs, including U.",1.5,8.5,6.3,col=soft)
 txt("At budget 18, signed 0.65 selects USF1 instead of NFE2L2. Lists at a fixed cutoff are nested by construction.",1.5,4.2,6.3,col=soft)
 popViewport()
}
export<-function(name,width,height,draw){stem<-file.path(base,"figures",name)
 grDevices::cairo_pdf(paste0(stem,".pdf"),width=width/25.4,height=height/25.4,family=family);draw();dev.off()
 svglite::svglite(paste0(stem,".svg"),width=width/25.4,height=height/25.4,bg="white");draw();dev.off()
 ragg::agg_png(paste0(stem,".png"),width=width,height=height,units="mm",res=600,background="white");draw();dev.off()
 ragg::agg_png(file.path(base,"review",paste0(name,"_preview.png")),width=width,height=height,units="mm",res=300,background="white");draw();dev.off()
}
export("G_TF_selection_background_compact",174,37,drawG)
export("I_TF_selection_sensitivity",174,90,drawI)
cat("G: full eligibility chain and six queues; I: 19 x 10 = 190 cells, 120 selected / 70 unselected; H reused byte-for-byte.\n")
