# Separate, additive enrichment analysis; locked protein models/results are read only.
suppressPackageStartupMessages({library(data.table);library(fgsea)})
prepare_fig4_protein_gsea <- function(root=getwd()) {
  out<-file.path(root,'4-Multiomic/refinements/fig4_protein_gsea_20261030')
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  files<-c('4-Multiomic/publication/tables/step4_mass_spec_study_adjusted_deg_main_subtypes.csv',
    '2-TranscriptomicCharacterization/runs/fig2_deg_patterns_20261026_v7/reference_gene_sets.csv.gz',
    '4-Multiomic/publication/tables/step4_mass_spec_ranked_gsea_main_subtypes.csv',
    '4-Multiomic/scripts/prepare_fig4_protein_gsea.R')
  hashes<-data.table(Path=files,MD5=unname(tools::md5sum(file.path(root,files))))
  hp<-file.path(out,'input_hashes.csv');result<-file.path(out,'all_collections_gsea.csv')
  if(file.exists(hp)&&file.exists(result)&&identical(fread(hp),hashes))return(invisible(out))
  deg<-fread(file.path(root,files[1]));sets<-unique(fread(file.path(root,files[2])))
  stopifnot(nrow(deg)==15681L,all(is.finite(deg$ModeratedT)),!anyDuplicated(deg[,.(Subtype,Protein)]),
    setequal(sets$Collection,c('GO_BP','Hallmark','Reactome','KEGG')),length(unique(sets$DBVersion))==1)
  universe<-sort(unique(deg$Protein));membership<-unique(sets[Gene%in%universe,.(Collection,Pathway,Gene)])
  coverage<-merge(sets[,.(ReferenceSize=uniqueN(Gene)),by=.(Collection,Pathway)],
    membership[,.(MeasuredSize=uniqueN(Gene)),by=.(Collection,Pathway)],all.x=TRUE)
  coverage[is.na(MeasuredSize),MeasuredSize:=0L]
  coverage[,Eligible:=MeasuredSize>=15 & MeasuredSize<=500]
  eligible<-coverage[Eligible==TRUE];membership<-membership[Pathway%in%eligible$Pathway]
  pathways<-split(membership$Gene,membership$Pathway)
  ans<-rbindlist(lapply(seq_along(c('C1','C2','C3')),function(i) {
    st<-c('C1','C2','C3')[i];z<-deg[Subtype==st];setorder(z,-ModeratedT,Protein)
    ranks<-setNames(z$ModeratedT,z$Protein);set.seed(20261030L+i)
    message('Protein four-collection GSEA ',st,': ',length(pathways),' pathways')
    fg<-as.data.table(fgseaMultilevel(pathways,ranks,minSize=15,maxSize=500,eps=0,
      nPermSimple=10000,BPPARAM=BiocParallel::SerialParam()))
    fg[,leadingEdge:=vapply(leadingEdge,paste,collapse=';',FUN.VALUE=character(1))]
    fg[,Subtype:=st];fg
  }))
  ans<-merge(ans,eligible[,.(pathway=Pathway,Collection)],by='pathway')
  setnames(ans,'padj','FDRAllCollections')
  ans[,FDRWithinCollection:=p.adjust(pval,'BH'),by=.(Subtype,Collection)]
  stopifnot(!anyNA(ans$pval),all(is.finite(ans$NES)),all(table(ans$Subtype)==nrow(eligible)),
    ans[,all(abs(FDRAllCollections-p.adjust(pval,'BH'))<1e-12),by=Subtype]$V1)
  old<-fread(file.path(root,files[3]));audit<-merge(old[,.(Subtype,pathway,OldES=ES,OldSize=size,
    OldNES=NES,OldFDR=padj)],ans[,.(Subtype,pathway,ES,size,NES,FDRAllCollections)],by=c('Subtype','pathway'))
  audit[,`:=`(ESMatches=abs(ES-OldES)<1e-8,SizeMatches=size==OldSize)]
  fwrite(ans,file.path(out,'all_collections_gsea.csv'))
  fwrite(coverage,file.path(out,'all_pathway_coverage.csv'))
  fwrite(membership,file.path(out,'measured_membership.csv.gz'))
  fwrite(audit,file.path(out,'locked_reference_comparison.csv'))
  fwrite(ans[,.(Tested=.N,PositiveFDR05=sum(FDRAllCollections<.05&NES>0),
    NegativeFDR05=sum(FDRAllCollections<.05&NES<0)),by=.(Collection,Subtype)],file.path(out,'summary.csv'))
  fwrite(hashes,hp);writeLines(capture.output(sessionInfo()),file.path(out,'session_info.txt'))
  writeLines(c('Status: EXECUTED_REVIEW_DRAFT; additive four-collection GSEA, no protein-model refit.',
    'Input: all 5227 frozen Study-adjusted moderated t values per principal subtype (n=248).',
    'Reference: local frozen MSigDB 2026.1.Hs GO BP/Hallmark/Reactome/KEGG MEDICUS membership snapshot; no RNA statistics are reused.',
    'Measured set size 15-500; fgseaMultilevel eps=0 nPermSimple=10000 SerialParam; seeds 20261031/32/33.',
    'Primary display FDR: BH over all tested pathways from all four collections within each subtype.',
    'FDRWithinCollection is supplied as a separately named secondary summary, not used for display.',
    'Legacy Hallmark/Reactome NES and FDR remain unchanged in publication tables; new estimates/FDR are separately versioned.',
    'KEGG is CP:KEGG_MEDICUS, not KEGG LEGACY. Pathway overlap is not independent biological evidence.'),
    file.path(out,'ANALYSIS_CONTRACT.md'))
  invisible(out)
}
if(sys.nframe()==0)prepare_fig4_protein_gsea()
