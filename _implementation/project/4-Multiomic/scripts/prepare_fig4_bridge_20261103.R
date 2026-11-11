suppressPackageStartupMessages({library(data.table);library(limma);library(fgsea)})
out <- '4-Multiomic/refinements/fig4_bridge_20261103_v1'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
inputs <- c('governance/data/subtype_master_table.csv',
 '1-Clustering/data/combined_expr_after_correction.rds',
 '1-Clustering/publication/tables/step1_main_only_study_adjusted_deg.csv',
 '4-Multiomic/publication/tables/step4_matched_rna_protein_subtype_contrasts_long.csv',
 '4-Multiomic/publication/tables/step4_mass_spec_study_adjusted_deg_main_subtypes.csv',
 '8-Visualization/refinements/fig4_concordance_20261030_v6/source_data/matched_sample_roster.csv',
 '2-TranscriptomicCharacterization/runs/fig2_deg_patterns_20261026_v7/reference_gene_sets.csv.gz',
 '4-Multiomic/refinements/fig4_matched_four_gsea_20261030/gsea_by_view.csv',
 'manuscript/review/figure4_slides_20261102/Figure4G_84_cells_verified.csv')
hashes<-data.table(Path=inputs,MD5=unname(tools::md5sum(inputs)))
hp<-file.path(out,'input_hashes.csv')
if(file.exists(hp)&&!identical(fread(hp),hashes))stop('Bridge inputs changed; use a new versioned run directory before recomputation.')
fwrite(hashes,hp)
m <- fread(inputs[1])[SubtypeClass=='Main']; stopifnot(nrow(m)==497,!anyDuplicated(m$Sample))
e<-new.env();load(inputs[2],e);expr<-e$corrected_expr_limma
stopifnot(nrow(expr)==10531,all(m$Sample%in%colnames(expr)),all(is.finite(expr)))
old<-fread(inputs[4]);shared<-sort(unique(old$Gene)); stopifnot(length(shared)==4706,nrow(old)==14118)
roster<-fread(inputs[6]);mm<-m[match(roster$Sample,Sample)]
stopifnot(nrow(mm)==248,!anyNA(mm$Sample),identical(as.integer(table(mm$FinalSubtype)),c(51L,126L,71L)))
fit_equal <- function(meta,genes) {
 md<-as.data.frame(meta);md$FinalSubtype<-factor(md$FinalSubtype,c('C1','C2','C3'));md$Study<-factor(md$Study)
 d<-model.matrix(~0+FinalSubtype+Study,md);colnames(d)<-make.names(colnames(d))
 stopifnot(qr(d)$rank==ncol(d))
 cc<-makeContrasts(C1=FinalSubtypeC1-(FinalSubtypeC2+FinalSubtypeC3)/2,
 C2=FinalSubtypeC2-(FinalSubtypeC1+FinalSubtypeC3)/2,
 C3=FinalSubtypeC3-(FinalSubtypeC1+FinalSubtypeC2)/2,levels=d)
 ft<-eBayes(contrasts.fit(lmFit(expr[genes,md$Sample,drop=FALSE],d),cc),robust=TRUE)
 rbindlist(lapply(c('C1','C2','C3'),function(st)data.table(Gene=genes,Subtype=st,T=ft$t[,st],Effect=ft$coefficients[,st])))
}
r0<-fread(inputs[3])[,.(Gene,Subtype,T=t,Effect=logFC_subtype_minus_others)]
message('Fit equal-weight full RNA')
r1<-fit_equal(m,rownames(expr));message('Fit equal-weight matched RNA')
r2<-fit_equal(mm,rownames(expr));r3check<-fit_equal(mm,shared)
chk<-merge(r3check,old,by=c('Gene','Subtype'));err<-max(abs(chk$T-chk$RNAModeratedT))
stopifnot(err<1e-8)
r3<-old[,.(Gene,Subtype,T=RNAModeratedT,Effect=NA_real_)]
p3<-old[,.(Gene,Subtype,T=ProteinModeratedT,Effect=NA_real_)]
p4<-fread(inputs[5])[,.(Gene=Protein,Subtype,T=ModeratedT,Effect=StandardizedDifference)]
ranks<-rbindlist(list(R0=r0,R1=r1,R2=r2,R3=r3,P3=p3,P4=p4),idcol='Stage')
stopifnot(ranks[,all(is.finite(T))],!anyDuplicated(ranks[,.(Stage,Subtype,Gene)]))
fwrite(ranks,file.path(out,'gene_ranks.csv.gz'))
sets<-unique(fread(inputs[7])[,.(Collection,Pathway,Gene)])
cov<-rbindlist(lapply(unique(ranks$Stage),function(s){
 u<-unique(ranks[Stage==s,Gene]);sets[Gene%in%u,.(MeasuredSize=uniqueN(Gene)),by=.(Collection,Pathway)][,Stage:=s]
}))
elig<-cov[,.(Stages=.N,Eligible=all(MeasuredSize>=15&MeasuredSize<=500)),by=.(Collection,Pathway)][Stages==6&Eligible==TRUE]
fwrite(cov,file.path(out,'pathway_coverage.csv'));fwrite(elig,file.path(out,'common_test_family.csv'))
frozen<-fread(inputs[8]);answer<-list();idx<-0L
for(s in c('R0','R1','R2','R3','P3','P4'))for(st in c('C1','C2','C3')) {
 idx<-idx+1L;path<-file.path(out,paste0('gsea_',s,'_',st,'.csv'))
 if(file.exists(path)){answer[[idx]]<-fread(path);next}
 if(s%in%c('R3','P3')) {
  view<-if(s=='R3')'RNA' else 'Protein';fg<-copy(frozen[View==view&Subtype==st&pathway%in%elig$Pathway])
  fg[,OriginalQ:=FDRAllCollections];fg[,SeedUsed:=Seed]
 } else {
  z<-copy(ranks[Stage==s&Subtype==st]);setorder(z,-T,Gene);stat<-setNames(z$T,z$Gene)
  members<-sets[Pathway%in%elig$Pathway&Gene%in%z$Gene];pw<-split(members$Gene,members$Pathway)
  seed<-2026110300L+idx;set.seed(seed);message(s,' / ',st,': ',length(pw),' pathways')
  fg<-as.data.table(fgseaMultilevel(pw,stat,minSize=15,maxSize=500,eps=0,nPermSimple=10000,BPPARAM=BiocParallel::SerialParam()))
  fg[,leadingEdge:=vapply(leadingEdge,paste,collapse=';',FUN.VALUE=character(1))]
  fg[,`:=`(OriginalQ=NA_real_,SeedUsed=seed)]
 }
 fg[,`:=`(Stage=s,Subtype=st,QCommon=p.adjust(pval,'BH'))]
 stopifnot(nrow(fg)==nrow(elig),!anyNA(fg$pval),all(is.finite(fg$NES)))
 fwrite(fg,path);answer[[idx]]<-fg
}
ans<-rbindlist(answer,fill=TRUE);ans<-merge(ans,elig[,.(pathway=Pathway,CollectionCommon=Collection)],by='pathway')
fwrite(ans,file.path(out,'all_gsea_results.csv.gz'))
g<-unique(fread(inputs[9])[,.(Pathway,Program)])
g<-rbind(g,data.table(Pathway=c('HALLMARK_G2M_CHECKPOINT','GOBP_CYTOPLASMIC_TRANSLATION'),Program=c('G2M checkpoint','Cytoplasmic translation')))
g[,ProgramOrder:=seq_len(.N)];stopifnot(nrow(g)==16,all(g$Pathway%in%elig$Pathway))
display<-merge(ans,g,by.x='pathway',by.y='Pathway');fwrite(display,file.path(out,'displayed_programs.csv'))
fwrite(rbindlist(list(m[,.(Sample,Study,FinalSubtype,Population='Full RNA')],mm[,.(Sample,Study,FinalSubtype,Population='Matched')])),file.path(out,'sample_rosters.csv'))
fwrite(data.table(Stage=c('R0','R1','R2','R3','P3','P4'),N=c(497,497,248,248,248,248),Genes=c(10531,10531,10531,4706,4706,5227),Pathways=nrow(elig)),file.path(out,'stage_manifest.csv'))
fwrite(data.table(Check='Matched shared RNA t reproduction',MaxAbsoluteError=err,Pass=err<1e-8),file.path(out,'model_reproduction.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'session_info.txt'))
message('DONE; shared family = ',nrow(elig))
