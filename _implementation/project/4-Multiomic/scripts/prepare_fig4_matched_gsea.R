# Additive four-collection GSEA using frozen matched RNA/protein contrasts.
suppressPackageStartupMessages({library(data.table);library(fgsea)})
prepare_fig4_matched_gsea<-function(root=getwd()) {
  out<-file.path(root,'4-Multiomic/refinements/fig4_matched_four_gsea_20261030')
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  inputs<-c('4-Multiomic/publication/tables/step4_matched_rna_protein_subtype_contrasts_long.csv',
    '4-Multiomic/publication/tables/step4_matched_rna_protein_pathway_concordance.csv',
    '2-TranscriptomicCharacterization/runs/fig2_deg_patterns_20261026_v7/reference_gene_sets.csv.gz',
    '4-Multiomic/scripts/prepare_fig4_matched_gsea.R')
  hashes<-data.table(Path=inputs,MD5=unname(tools::md5sum(file.path(root,inputs))))
  hp<-file.path(out,'input_hashes.csv')
  if(file.exists(hp)&&file.exists(file.path(out,'paired_gsea.csv'))&&identical(fread(hp),hashes))return(invisible(out))
  d<-fread(file.path(root,inputs[1]));sets<-unique(fread(file.path(root,inputs[3])))
  stopifnot(nrow(d)==14118L,all(table(d$Subtype)==4706L),!anyDuplicated(d[,.(Subtype,Gene)]),
    all(is.finite(d$RNAModeratedT)),all(is.finite(d$ProteinModeratedT)))
  universe<-sort(unique(d$Gene))
  members<-unique(sets[Gene%in%universe,.(Collection,Pathway,Gene)])
  coverage<-merge(sets[,.(ReferenceSize=uniqueN(Gene)),by=.(Collection,Pathway)],
    members[,.(MeasuredSize=uniqueN(Gene)),by=.(Collection,Pathway)],all.x=TRUE)
  coverage[is.na(MeasuredSize),MeasuredSize:=0L];coverage[,Eligible:=MeasuredSize>=15&MeasuredSize<=500]
  eligible<-coverage[Eligible==TRUE];members<-members[Pathway%in%eligible$Pathway]
  pathways<-split(members$Gene,members$Pathway)
  answer<-list();i<-0L
  for(st in c('C1','C2','C3'))for(view in c('RNA','Protein')) {
    i<-i+1L;z<-copy(d[Subtype==st]);column<-paste0(view,'ModeratedT')
    setorderv(z,c(column,'Gene'),c(-1,1));ranks<-setNames(z[[column]],z$Gene)
    set.seed(20261100L+i)
    message('Matched four-collection GSEA ',st,' / ',view,': ',length(pathways),' sets')
    fg<-as.data.table(fgseaMultilevel(pathways,ranks,minSize=15,maxSize=500,eps=0,
      nPermSimple=10000,BPPARAM=BiocParallel::SerialParam()))
    fg[,leadingEdge:=vapply(leadingEdge,paste,collapse=';',FUN.VALUE=character(1))]
    fg[,`:=`(Subtype=st,View=view,Seed=20261100L+i)];answer[[i]]<-fg
  }
  ans<-merge(rbindlist(answer),eligible[,.(pathway=Pathway,Collection)],by='pathway')
  setnames(ans,'padj','FDRAllCollections')
  ans[,FDRWithinCollection:=p.adjust(pval,'BH'),by=.(Subtype,View,Collection)]
  stopifnot(nrow(ans)==nrow(eligible)*6L,!anyNA(ans$pval),all(is.finite(ans$NES)),
    ans[,all(abs(FDRAllCollections-p.adjust(pval,'BH'))<1e-12),by=.(Subtype,View)]$V1)
  rna<-ans[View=='RNA',.(Subtype,pathway,Collection,MeasuredSize=size,NES_RNA=NES,
    padj_RNA=FDRAllCollections,pval_RNA=pval,leadingEdge_RNA=leadingEdge)]
  pro<-ans[View=='Protein',.(Subtype,pathway,NES_Protein=NES,padj_Protein=FDRAllCollections,
    pval_Protein=pval,leadingEdge_Protein=leadingEdge)]
  paired<-merge(rna,pro,by=c('Subtype','pathway'))
  old<-fread(file.path(root,inputs[2]));audit<-merge(old,paired[,.(Subtype,pathway,
    NewNES_RNA=NES_RNA,NewNES_Protein=NES_Protein,NewQ_RNA=padj_RNA,NewQ_Protein=padj_Protein)],
    by=c('Subtype','pathway'))
  audit[,`:=`(RNA_SignMatches=sign(NES_RNA)==sign(NewNES_RNA),
    Protein_SignMatches=sign(NES_Protein)==sign(NewNES_Protein))]
  fwrite(ans,file.path(out,'gsea_by_view.csv'));fwrite(paired,file.path(out,'paired_gsea.csv'))
  fwrite(coverage,file.path(out,'coverage.csv'));fwrite(members,file.path(out,'membership.csv.gz'))
  fwrite(audit,file.path(out,'legacy_comparison.csv'));fwrite(hashes,hp)
  fwrite(paired[,.(Pathways=.N,BothFDR05=sum(padj_RNA<.05&padj_Protein<.05),
    SameSign=sum(sign(NES_RNA)==sign(NES_Protein))),by=.(Subtype,Collection)],file.path(out,'summary.csv'))
  writeLines(capture.output(sessionInfo()),file.path(out,'session_info.txt'))
  writeLines(c('Status: EXECUTED_REVIEW_DRAFT; additive matched four-collection enrichment.',
    'Frozen moderated t ranks: same248 matched Main samples, same4706 genes, same Study-adjusted equal-other-subtype contrasts.',
    'MSigDB2026.1.Hs local snapshot: GO BP, Hallmark, Reactome, KEGG MEDICUS; measured set size15-500.',
    'fgseaMultilevel eps=0 nPermSimple=10000 SerialParam; seeds20261101-20261106, ordered C1 RNA/Protein, C2 RNA/Protein, C3 RNA/Protein.',
    'Display q is BH across all eligible pathways from four collections within each subtype and modality.',
    'No per-selected-set FDR or mixing full-RNA497 / full-protein5227-gene GSEA with this matched comparison.',
    'Locked publication models/GSEA tables unchanged. New stochastic NES/q are versioned here.',
    'Pathway overlap, same-cohort selection and bulk composition preclude independent-validation/causal claims.'),file.path(out,'ANALYSIS_CONTRACT.md'))
  invisible(out)
}
if(sys.nframe()==0)prepare_fig4_matched_gsea()
