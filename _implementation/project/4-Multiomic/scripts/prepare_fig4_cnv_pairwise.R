# Follow-up contrasts for every cytoband/direction, before omnibus display filtering.
prepare_fig4_cnv_pairwise <- function(root) {
  library(data.table)
  base<-'4-Multiomic/refinements/fig4_cnv_tests_20261030/'
  out<-file.path(root,'4-Multiomic/refinements/fig4_cnv_pairwise_20261030')
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  paths<-c(paste0(base,c('cytoband_sample_fractions','cytoband_tests','input_hashes'),'.csv'),
    '4-Multiomic/scripts/prepare_fig4_cnv_pairwise.R')
  hashes<-data.table(Path=paths,MD5=unname(tools::md5sum(file.path(root,paths))))
  if(file.exists(file.path(out,'input_hashes.csv'))&&identical(fread(file.path(out,'input_hashes.csv')),hashes)&&
    file.exists(file.path(out,'pairwise_display.csv')))return(invisible(out))
  d<-fread(file.path(root,paths[1]));global<-fread(file.path(root,paths[2]))
  pairs<-list(c('C2','C1'),c('C3','C1'),c('C3','C2'))
  tests<-rbindlist(lapply(seq_along(pairs),function(k) {
    pair<-pairs[[k]]
    d[FinalSubtype%in%pair,{
      a<-AlteredGeneFraction[FinalSubtype==pair[1]];b<-AlteredGeneFraction[FinalSubtype==pair[2]]
      invariant<-length(unique(c(a,b)))==1L
      fit<-if(!invariant)wilcox.test(a,b,exact=FALSE,correct=TRUE) else NULL
      .(Target=pair[1],Reference=pair[2],Contrast=paste(pair,collapse=' vs '),ContrastOrder=k,
        TargetN=length(a),ReferenceN=length(b),TargetMean=mean(a),ReferenceMean=mean(b),
        DifferencePP=100*(mean(a)-mean(b)),PValue=if(invariant)1 else fit$p.value,
        RankBiserial=if(invariant)0 else 2*unname(fit$statistic)/(length(a)*length(b))-1,
        TestStatus=if(invariant)'Invariant; assigned P=1' else 'Two-sided Wilcoxon; continuity and tie correction')
    },by=.(Cytoband,Direction)]
  }))
  stopifnot(nrow(tests)==4602L,!anyNA(tests$PValue))
  tests[,FDR:=p.adjust(PValue,'BH')]
  tests<-merge(tests,global[,.(Cytoband,Direction,OmnibusP=PValue,OmnibusFDR=FDR,Chromosome,X)],
    by=c('Cytoband','Direction'))
  tests[,Displayed:=OmnibusFDR<.05]
  tests[,DirectionOrder:=ifelse(Direction=='Gain',1L,2L)]
  setorder(tests,DirectionOrder,X,ContrastOrder)
  # Same global18 entries as S10Q; all54 contrasts remain, including nonsignificance.
  display<-tests[Displayed==TRUE]
  stopifnot(nrow(display)==54L,uniqueN(display[,.(Cytoband,Direction)])==18L)
  patterns<-display[,{
    q<-setNames(FDR,Contrast);diff<-setNames(DifferencePP,Contrast)
    c3both<-all(q[c('C3 vs C1','C3 vs C2')]<.05)&&prod(sign(diff[c('C3 vs C1','C3 vs C2')]))>0
    c2both<-all(q[c('C2 vs C1','C3 vs C2')]<.05)&&prod(sign(c(diff['C2 vs C1'],-diff['C3 vs C2'])))>0
    c1both<-all(q[c('C2 vs C1','C3 vs C1')]<.05)&&prod(sign(-diff[c('C2 vs C1','C3 vs C1')]))>0
    .(SignificantContrasts=sum(FDR<.05),C1DiffersFromBoth=c1both,C2DiffersFromBoth=c2both,C3DiffersFromBoth=c3both)
  },by=.(Cytoband,Direction,OmnibusFDR)]
  summary<-data.table(AllPairwiseTests=nrow(tests),SignificantAllPairwise=sum(tests$FDR<.05),
    DisplayedEntries=18L,DisplayedContrasts=nrow(display),SignificantDisplayedContrasts=sum(display$FDR<.05))
  fwrite(tests,file.path(out,'pairwise_all.csv'));fwrite(display,file.path(out,'pairwise_display.csv'))
  fwrite(patterns,file.path(out,'region_patterns.csv'));fwrite(summary,file.path(out,'summary.csv'))
  fwrite(hashes,file.path(out,'input_hashes.csv'))
  print(summary);print(display[FDR<.05,.(Cytoband,Direction,Contrast,DifferencePP,FDR)]);print(patterns)
  invisible(out)
}
