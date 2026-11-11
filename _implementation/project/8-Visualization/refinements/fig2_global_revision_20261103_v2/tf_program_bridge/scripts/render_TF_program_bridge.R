#!/usr/bin/env Rscript
# Descriptive joins from fixed source tables; no activity/DEG/network refit or new significance test.
options(stringsAsFactors = FALSE, warn = 2)
library(grid)
root <- normalizePath(getwd())
base <- file.path(root,"8-Visualization/refinements/fig2_global_revision_20261103_v2/tf_program_bridge")
inputs <- c(selection="8-Visualization/refinements/fig2D_selection_20261031_v33/tables/D_selection12.csv",
 edges="8-Visualization/refinements/fig2D_selection_20261031_v33/tables/D_display_union_regulon_edges.csv",
 eligibility="8-Visualization/refinements/fig2D_selection_20261031_v33/tables/D_all263_eligibility.csv",
 patterns="8-Visualization/refinements/fig2_assembly_20261031_v34/A/tables/row_order.csv",
 pattern_layout="8-Visualization/refinements/fig2_assembly_20261031_v34/A/tables/pattern_layout.csv",
 modules="2-TranscriptomicCharacterization/publication/tables/step2_wgcna_gene_module_membership.csv",
 main_modules="8-Visualization/refinements/fig2_assembly_20261031_v34/tables/C_optional_heterogeneity_audit.csv",
 deg="1-Clustering/publication/tables/step1_main_only_study_adjusted_deg.csv")
src <- lapply(file.path(root,inputs),read.csv,check.names=FALSE); names(src)<-names(inputs)
sel<-src$selection; edges<-src$edges[src$edges$TF %in% sel$TF,]; pat<-src$patterns; mods<-src$modules
eg<-src$eligibility; deg<-src$deg
for(f in c("CoveragePass","OriginalSignificant","LOSOStable","Eligible")) eg[[f]] <- tolower(as.character(eg[[f]]))=="true"
stopifnot(nrow(sel)==12L,!anyDuplicated(sel$TF),nrow(pat)==842L,!anyDuplicated(pat$Gene),
          length(unique(pat$Pattern))==6L,nrow(mods)==6753L,!anyDuplicated(mods$Gene),
          nrow(src$main_modules)==11L,nrow(eg)==789L,length(unique(eg$TF))==263L,
          nrow(deg)==31593L,!anyDuplicated(paste(deg$Gene,deg$Subtype)),
          !anyDuplicated(paste(edges$TF,edges$Gene)),all(edges$Sign==sign(edges$mor)),all(edges$Sign!=0))
# Order matches selection-anchor groups; anchor identity and membership stay unchanged.
sel<-sel[order(match(sel$AnchorSubtype,c("C1","C2","C3")),match(sel$AnchorDirection,c("Positive","Negative")),sel$Round),]
sel$ColumnOrder<-seq_len(nrow(sel));sel$AnchorSign<-sign(sel$AnchorEffect)
sel$AnchorLabel<-paste0(sel$AnchorSubtype,ifelse(sel$AnchorSign>0,"+","\u2212"))
sel$ActualMeasuredTargets<-as.integer(table(factor(edges$TF,levels=sel$TF)))
stopifnot(all(sel$ActualMeasuredTargets==sel$NonzeroSignedTargets))
ix<-match(edges$TF,sel$TF)
edges$AnchorSubtype<-sel$AnchorSubtype[ix];edges$AnchorDirection<-sel$AnchorDirection[ix]
edges$AnchorEffect<-sel$AnchorEffect[ix];edges$AnchorFDR<-sel$AnchorFDR[ix]
edges$ExpectedTargetSign<-sign(edges$AnchorEffect)*edges$Sign
di<-match(paste(edges$Gene,edges$AnchorSubtype),paste(deg$Gene,deg$Subtype))
stopifnot(!anyNA(di))
for(f in c("Contrast","logFC_subtype_minus_others","SE","PValue","FDR","NSubtype","NOthers","NTotal","StudyAdjusted")) edges[[paste0("Gene_",f)]]<-deg[[f]][di]
edges$ObservedTargetSign<-sign(edges$Gene_logFC_subtype_minus_others)
edges$DirectionClass<-ifelse(!is.finite(edges$Gene_FDR)|!is.finite(edges$Gene_logFC_subtype_minus_others),"Unavailable",
  ifelse(edges$Gene_FDR<.05 & edges$ObservedTargetSign!=0,
         ifelse(edges$ObservedTargetSign==edges$ExpectedTargetSign,"Concordant significant","Opposite significant"),"Not significant"))
stopifnot(!any(edges$DirectionClass=="Unavailable"),all(edges$Gene_NTotal==497L))
edges$Pattern<-pat$Pattern[match(edges$Gene,pat$Gene)]
edges$WGCNAModule<-mods$Module[match(edges$Gene,mods$Gene)]
edges$InPatternUniverse<-!is.na(edges$Pattern);edges$InWGCNAUniverse<-!is.na(edges$WGCNAModule)
edges$InMainWGCNAModule<-edges$WGCNAModule %in% src$main_modules$Module
pt<-src$pattern_layout[order(src$pattern_layout$PatternOrder),]
programs<-rbind(data.frame(Family="Patterns",ProgramID=pt$Pattern,ProgramLabel=pt$Pattern,ProgramGenes=pt$Genes,UniverseGenes=842L),
 data.frame(Family="WGCNA",ProgramID=src$main_modules$Module,ProgramLabel=src$main_modules$Program,
            ProgramGenes=src$main_modules$ModuleSize,UniverseGenes=6753L))
programs$ProgramOrder<-seq_len(nrow(programs))
for(i in seq_len(nrow(programs))) {
 pp<-programs[i,];n<-if(pp$Family=="Patterns")sum(pat$Pattern==pp$ProgramID)else sum(mods$Module==pp$ProgramID)
 stopifnot(n==pp$ProgramGenes)
}
counts<-list();target_program<-list();k<-0L
for(i in seq_len(nrow(programs)))for(j in seq_len(nrow(sel))) {
 p<-programs[i,];s<-sel[j,];e<-edges[edges$TF==s$TF,]
 included<-if(p$Family=="Patterns")!is.na(e$Pattern)&e$Pattern==p$ProgramID else !is.na(e$WGCNAModule)&e$WGCNAModule==p$ProgramID
 e<-e[included,];k<-k+1L
 cc<-as.integer(table(factor(e$DirectionClass,levels=c("Concordant significant","Opposite significant","Not significant","Unavailable"))))
 counts[[k]]<-data.frame(p,TF=s$TF,ColumnOrder=j,AnchorSubtype=s$AnchorSubtype,AnchorDirection=s$AnchorDirection,
   AnchorEffect=s$AnchorEffect,AnchorFDR=s$AnchorFDR,MeasuredSignedTargets=s$ActualMeasuredTargets,
   OverlapTargets=nrow(e),FractionOfMeasuredTargets=nrow(e)/s$ActualMeasuredTargets,FractionOfProgramGenes=nrow(e)/p$ProgramGenes,
   ConcordantSignificant=cc[1],OppositeSignificant=cc[2],NotSignificant=cc[3],Unavailable=cc[4],
   ConcordantFractionOfMeasuredTargets=cc[1]/s$ActualMeasuredTargets,OppositeFractionOfMeasuredTargets=cc[2]/s$ActualMeasuredTargets,
   NSFractionOfMeasuredTargets=cc[3]/s$ActualMeasuredTargets,Genes=paste(e$Gene,collapse=";"))
 if(nrow(e)) {e$ProgramFamily<-p$Family;e$ProgramID<-p$ProgramID;e$ProgramLabel<-p$ProgramLabel;target_program[[length(target_program)+1L]]<-e}
}
counts<-do.call(rbind,counts);target_program<-do.call(rbind,target_program)
stopifnot(nrow(counts)==204L,all(rowSums(counts[c("ConcordantSignificant","OppositeSignificant","NotSignificant","Unavailable")])==counts$OverlapTargets))
total<-do.call(rbind,lapply(seq_len(nrow(sel)),function(j){s<-sel[j,];e<-edges[edges$TF==s$TF,];data.frame(TF=s$TF,
 MeasuredSignedTargets=nrow(e),ConcordantSignificant=sum(e$DirectionClass=="Concordant significant"),OppositeSignificant=sum(e$DirectionClass=="Opposite significant"),
 NotSignificant=sum(e$DirectionClass=="Not significant"),PatternTargets=sum(e$InPatternUniverse),MainModuleTargets=sum(e$InMainWGCNAModule),WGCNAUniverseTargets=sum(e$InWGCNAUniverse))}))
write.csv(sel,file.path(base,"source_data/current12_selection_anchors.csv"),row.names=FALSE)
write.csv(edges,file.path(base,"source_data/current12_all_signed_targets_anchor_DEG.csv"),row.names=FALSE,na="")
write.csv(programs,file.path(base,"source_data/fixed17_program_definitions.csv"),row.names=FALSE)
write.csv(counts,file.path(base,"source_data/current12_by_17programs_204_cells.csv"),row.names=FALSE)
write.csv(target_program,file.path(base,"source_data/all_overlapping_target_program_records.csv"),row.names=FALSE,na="")
write.csv(total,file.path(base,"source_data/current12_target_direction_totals.csv"),row.names=FALSE)
# Full background, including original nonsignificant/low-coverage/missing-LOSO rows.
eg$SelectedMainTF<-eg$TF %in% sel$TF
eg$SelectionAnchor<-paste(eg$TF,eg$Subtype)==paste(sel$TF[match(eg$TF,sel$TF)],sel$AnchorSubtype[match(eg$TF,sel$TF)])
eg$SelectionAnchor[is.na(eg$SelectionAnchor)]<-FALSE
write.csv(eg,file.path(base,"source_data/all263_TF_789_original_contrasts.csv"),row.names=FALSE,na="")
bg<-do.call(rbind,lapply(sort(unique(eg$TF)),function(tf){d<-eg[eg$TF==tf,];stopifnot(length(unique(d$NonzeroSignedTargets))==1L)
 z<-data.frame(TF=tf,TargetsAvailable=d$TargetsAvailable[1],MeasuredSignedTargets=d$NonzeroSignedTargets[1],
 CoveragePass=d$CoveragePass[1],AnyOriginalSignificant=any(d$OriginalSignificant),AnyEligible=any(d$Eligible),SelectedMain=tf%in%sel$TF,
 ExistingLOSOStudies=max(c(0,d$LOSOStudies),na.rm=TRUE),SelectionAnchor=if(tf%in%sel$TF)sel$AnchorLabel[match(tf,sel$TF)]else"")
 for(st in c("C1","C2","C3")){r<-d[d$Subtype==st,];for(f in c("ActivityDifferenceSubtypeMinusOthers","FDR","OriginalSignificant","DirectionRetained","LOSOStable","Eligible","ExclusionReason"))z[[paste(st,f,sep="_")]]<-r[[f]]};z}))
stopifnot(all(bg$MeasuredSignedTargets > 0))
write.csv(bg,file.path(base,"source_data/all263_TF_background.csv"),row.names=FALSE,na="")
flow_pred<-list(rep(TRUE,nrow(eg)),eg$OriginalSignificant,eg$OriginalSignificant&eg$CoveragePass,eg$Eligible,eg$SelectionAnchor)
flow<-data.frame(Stage=c("Tested","Original FDR < 0.05","+ >=10 measured signed targets","+ 7/7 LOSO direction retained","Display-rule anchors"),
 TFs=sapply(flow_pred,function(x)length(unique(eg$TF[x]))),Contrasts=sapply(flow_pred,sum))
stopifnot(identical(as.integer(flow$TFs),c(263L,250L,225L,225L,12L)),identical(as.integer(flow$Contrasts),c(789L,585L,526L,524L,12L)))
write.csv(flow,file.path(base,"source_data/all263_sequential_selection_flow.csv"),row.names=FALSE)
queues<-expand.grid(Subtype=c("C1","C2","C3"),Direction=c("Positive","Negative"),stringsAsFactors=FALSE)
queues$Label<-paste0(queues$Subtype,ifelse(queues$Direction=="Positive","+","\u2212"))
queues$Eligible<-mapply(function(st,dir)sum(eg$Subtype==st&eg$Direction==dir&eg$Eligible),queues$Subtype,queues$Direction)
queues$SelectedAnchors<-mapply(function(st,dir)sum(sel$AnchorSubtype==st&sel$AnchorDirection==dir),queues$Subtype,queues$Direction)
write.csv(queues,file.path(base,"source_data/selection_queue_counts.csv"),row.names=FALSE)
for(nm in names(inputs))file.copy(file.path(root,inputs[[nm]]),file.path(base,"source_data",paste0("input_",nm,".csv")),overwrite=TRUE)
write.csv(data.frame(Source=unname(inputs),MD5=unname(tools::md5sum(file.path(root,inputs)))),file.path(base,"review/input_hashes.csv"),row.names=FALSE)

family<-"Arial";ink<-"#26313A";soft<-"#616A73";rule<-"#D9DEE2"
cols<-c(`Concordant significant`="#247AA4",`Opposite significant`="#D56B43",`Not significant`="#B7BEC4",Unavailable="#6B4C9A")
mm<-function(x)unit(x,"mm")
txt<-function(s,x,y,size=6.8,face="plain",col=ink,just="left")grid.text(s,x=mm(x),y=mm(y),just=just,gp=gpar(fontfamily=family,fontsize=size,fontface=face,col=col,lineheight=1.05))
seg<-function(x1,y1,x2,y2,col=rule,lwd=.45)grid.segments(mm(x1),mm(y1),mm(x2),mm(y2),gp=gpar(col=col,lwd=lwd))
rect<-function(x,y,w,h,fill,col=NA,lwd=.35)grid.rect(mm(x),mm(y),mm(w),mm(h),gp=gpar(fill=fill,col=col,lwd=lwd))
draw_bridge<-function(){
 grid.newpage();pushViewport(viewport(x=0,y=0,width=mm(174),height=mm(114),just=c("left","bottom")));rect(87,57,174,114,"white")
 txt("Measured TF targets in fixed expression programs",1.5,110.5,8,"bold")
 txt("Overlap counts and signed target-expression correspondence at each TF's selection anchor",1.5,106.6,6.8,col=soft)
 cx<-seq(56.1,166.1,by=10)
 txt("Program [gene count]",1.8,95.1,6.8,"bold")
 for(j in 1:12){txt(sel$TF[j],cx[j],101.7,6.9,"bold",just="center");txt(sel$AnchorLabel[j],cx[j],98.5,6.6,just="center");txt(paste0("n=",sel$ActualMeasuredTargets[j]),cx[j],95.2,6.2,col=soft,just="center")}
 seg(1.5,93.1,172,93.1)
 y<-90.9
 for(fam in c("Patterns","WGCNA")){
  rect(86.75,y,170.5,2.8,"#F0F3F5")
  txt(if(fam=="Patterns")"P1\u2013P6: fixed 842-gene pattern membership"else"Main WGCNA programs: full membership from the 6,753-gene network",1.8,y,6.7,"bold")
  y<-y-3.7
  for(id in programs$ProgramID[programs$Family==fam]){
   p<-programs[programs$ProgramID==id,];lab<-paste0(p$ProgramLabel," [",p$ProgramGenes,"]")
   txt(lab,2.2,y,if(fam=="Patterns")6.8 else 6.4)
   for(j in 1:12){d<-counts[counts$ProgramID==id&counts$TF==sel$TF[j],];xx<-cx[j];
    txt(d$OverlapTargets,xx,y+.6,6.6,just="center",col=if(d$OverlapTargets==0)soft else ink)
    left<-xx-4.35;rect(xx,y-1.12,8.7,.85,"white",rule,lwd=.3)
    vv<-unlist(d[c("ConcordantSignificant","OppositeSignificant","NotSignificant","Unavailable")])/d$MeasuredSignedTargets
    for(k in seq_along(vv))if(vv[k]>0){w<-8.7*vv[k];rect(left+w/2,y-1.12,w,.85,unname(cols[k]));left<-left+w}
   }
   y<-y-4.1
  }
  y<-y-.1
 }
 for(k in 1:3){x<-c(1.8,58,117)[k];rect(x+1.1,10.8,2.2,1.2,cols[k]);txt(c("Expected direction, gene FDR < 0.05","Opposite direction, gene FDR < 0.05","Not significant")[k],x+3.1,10.8,6.3)}
 txt("Count above bar; full bar width = all measured signed targets for that TF (n in header). Empty portion = targets outside that program.",1.5,6.5,6.2,col=soft)
 txt("Expected gene direction = anchor activity sign \u00d7 signed mode of regulation. Shared expression evidence; no enrichment or causal test.",1.5,2.3,6.2,col=soft)
 popViewport()
}
draw_background<-function(){
 grid.newpage();pushViewport(viewport(x=0,y=0,width=mm(174),height=mm(74),just=c("left","bottom")));rect(87,37,174,74,"white")
 txt("Full TF background and display selection",1.5,70.4,8,"bold")
 txt("263 tested TFs; 789 original contrasts. Selection anchors are distinct from the 36 displayed counterpart contrasts.",1.5,66.3,6.8,col=soft)
 txt("Sequential criteria",1.5,60.8,7,"bold");txt("TFs",66,60.8,6.8,"bold",just="center");txt("Contrasts",82,60.8,6.8,"bold",just="center")
 seg(1.5,58.4,88,58.4)
 labels<-c("Tested","Original FDR < 0.05","+ \u226510 measured signed targets","+ 7/7 LOSO direction retained","Display-rule anchors")
 yy<-c(54.5,48.5,42.5,36.5,28.5)
 for(i in 1:5){txt(labels[i],1.8,yy[i],6.7,if(i==5)"bold"else"plain");txt(flow$TFs[i],66,yy[i],7,if(i==5)"bold"else"plain",just="center");txt(flow$Contrasts[i],82,yy[i],7,if(i==5)"bold"else"plain",just="center")}
 seg(1.5,32.7,88,32.7)
 txt("Missing LOSO: 13 TFs / 39 contrasts; none passed.",1.5,22.7,6.3,col=soft)
 txt("Two original-significant, coverage-qualified contrasts",1.5,18.8,6.3,col=soft)
 txt("fail 7/7 sign retention; all 225 eligible TFs remain.",1.5,15.2,6.3,col=soft)
 txt("Anchor queues",108,60.8,7,"bold",just="center");txt("Shown / eligible",105,56.8,6.2,col=soft,just="center")
 for(i in 1:6){y<-51.6-(i-1)*5.5;txt(queues$Label[i],92,y,6.8);txt(paste0(queues$SelectedAnchors[i]," / ",queues$Eligible[i]),114,y,6.6,just="center")}
 # Every tested TF is visible as a point; rank is a display coordinate, not an inferential quantity.
 txt("Measured signed targets\n(log scale)",150,60.1,6.6,"bold",just="center")
 b<-bg[order(bg$MeasuredSignedTargets,bg$TF),];b$Rank<-seq_len(nrow(b));x0<-133;x1<-171;y0<-22;y1<-55
 xfun<-function(v)x0+(v-1)/262*(x1-x0);yfun<-function(v)y0+(log10(v)-log10(5))/(log10(314)-log10(5))*(y1-y0)
 for(v in c(5,10,50,314)){seg(x0,yfun(v),x1,yfun(v),if(v==10)"#A7ADB3"else"#EBEEF0");txt(v,x0-1.4,yfun(v),6.0,col=soft,just="right")}
 for(i in which(!b$SelectedMain))grid.circle(mm(xfun(b$Rank[i])),mm(yfun(b$MeasuredSignedTargets[i])),r=mm(.28),gp=gpar(fill="#AEB7BF",col=NA))
 for(i in which(b$SelectedMain))grid.circle(mm(xfun(b$Rank[i])),mm(yfun(b$MeasuredSignedTargets[i])),r=mm(.47),gp=gpar(fill="#247AA4",col="white",lwd=.25))
 txt("1",x0,19.5,6.0,col=soft,just="center");txt("263",x1,19.5,6.0,col=soft,just="center")
 txt("TFs ordered by target count",152,15.3,6.1,col=soft,just="center")
 grid.circle(mm(134),mm(11.8),r=mm(.47),gp=gpar(fill="#247AA4",col=NA));txt("Main 12",136,11.8,6.2)
 txt("Original FDR families: 263 TFs within each subtype contrast. Eligibility is not causal validation; 12 is a display budget.",1.5,7.0,6.3,col=soft)
 txt("Target n counts measured genes, not samples. Coverage threshold is 10; full 263-TF and 789-contrast tables retain exclusions.",1.5,2.6,6.3,col=soft)
 popViewport()
}
export<-function(stem,width,height,draw){
 out<-file.path(base,"figures",stem)
 grDevices::cairo_pdf(paste0(out,".pdf"),width=width/25.4,height=height/25.4,family=family);draw();dev.off()
 svglite::svglite(paste0(out,".svg"),width=width/25.4,height=height/25.4,bg="white");draw();dev.off()
 ragg::agg_png(paste0(out,".png"),width=width,height=height,units="mm",res=600,background="white");draw();dev.off()
 ragg::agg_png(file.path(base,"review",paste0(stem,"_preview.png")),width=width,height=height,units="mm",res=300,background="white");draw();dev.off()
}
export("TF_program_correspondence",174,114,draw_bridge)
export("TF_all263_background",174,74,draw_background)
cat(sprintf("TF bridge: %d measured signed edges; %d complete grid cells, %d explicit zero cells; %d overlapping TF-target-program records.\n",nrow(edges),nrow(counts),sum(counts$OverlapTargets==0L),nrow(target_program)))
print(total);print(flow)
