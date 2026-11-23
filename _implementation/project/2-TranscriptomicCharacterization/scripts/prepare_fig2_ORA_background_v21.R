suppressPackageStartupMessages(library(data.table))
prepare_ORA_background_v21<-function(){
  out<-'2-TranscriptomicCharacterization/runs/fig2_ORA_background_20261030_v21'
  paths<-c('2-TranscriptomicCharacterization/runs/fig2_soft_patterns_20261030_v20/K6/gene_pattern_membership.csv',
    '2-TranscriptomicCharacterization/runs/fig2_six_patterns_20261030_v19/reference_membership_used.csv.gz',
    '1-Clustering/publication/tables/step1_main_only_study_adjusted_deg.csv',
    '2-TranscriptomicCharacterization/runs/fig2_soft_patterns_20261030_v20/actual_pattern_ORA_all_tests.csv.gz',
    '2-TranscriptomicCharacterization/runs/fig2_soft_patterns_20261030_v20/K6/pattern_counts.csv',
    file.path(out,'CONTRACT.md'),'2-TranscriptomicCharacterization/scripts/prepare_fig2_ORA_background_v21.R')
  hashes<-unname(tools::md5sum(paths));stopifnot(!anyNA(hashes))
  mem<-fread(paths[1]);sets<-unique(fread(paths[2]));tested<-unique(fread(paths[3])$Gene);baseline<-fread(paths[4])[K==6]
  stopifnot(length(tested)==10531,nrow(mem)==842,all(mem$Gene %in% tested))
  baseline[,Background:='HVG2000_annotated'];baseline[,K:=NULL]
  coverage<-list();res<-list();i<-1L;j<-1L
  for(mode in c('Full_tested_10531','Full_tested_annotated'))for(col in unique(sets$Collection)){
    ss<-sets[Collection==col];univ<-if(mode=='Full_tested_10531')tested else intersect(tested,unique(ss$Gene))
    bg<-ss[Gene %in% univ];sz<-bg[,.(SetSize=uniqueN(Gene)),by=.(Collection,Pathway)][SetSize>=10 & SetSize<=500]
    bg<-merge(bg,sz,by=c('Collection','Pathway'))
    for(p in paste0('P',1:6)){
      fg<-intersect(mem[Pattern==p,Gene],univ)
      hit<-bg[Gene %in% fg,.(Count=uniqueN(Gene),HitGenes=paste(sort(unique(Gene)),collapse='/')),by=.(Collection,Pathway)]
      dd<-merge(sz,hit,by=c('Collection','Pathway'),all.x=TRUE);dd[is.na(Count),`:=`(Count=0L,HitGenes='')]
      dd[,`:=`(Background=mode,Pattern=p,UniverseGenes=length(univ),SelectedGenes=length(fg))]
      dd[,PValue:=phyper(Count-1,SetSize,UniverseGenes-SetSize,SelectedGenes,lower.tail=FALSE)]
      dd[,FoldEnrichment:=(Count/SelectedGenes)/(SetSize/UniverseGenes)]
      res[[i]]<-dd;i<-i+1L
      coverage[[j]]<-data.table(Background=mode,Collection=col,Pattern=p,PatternGenes=mem[Pattern==p,.N],MappedForeground=length(fg),
        UniverseGenes=length(univ),CollectionAnnotatedTestedGenes=uniqueN(ss[Gene %in% tested,Gene]),TestedTerms=nrow(sz));j<-j+1L
    }
  }
  res<-rbindlist(res);res[,FDR:=p.adjust(PValue,'BH'),by=.(Background,Collection)];res[,FamilySize:=.N,by=.(Background,Collection)]
  baseline[,FoldEnrichment:=(Count/SelectedGenes)/(SetSize/UniverseGenes)]
  all<-rbindlist(list(baseline,res),use.names=TRUE,fill=TRUE);setorder(all,Background,Collection,Pattern,FDR,-Count,Pathway)
  fwrite(all,file.path(out,'background_ORA_all_tests.csv.gz'));fwrite(rbindlist(coverage),file.path(out,'background_coverage.csv'))
  counts<-all[,.(Tests=.N,FDR05=sum(FDR<.05),DisplayCandidates=sum(FDR<.05 & Count>=5)),by=.(Background,Collection,Pattern)]
  fwrite(counts,file.path(out,'significant_counts.csv'))
  comparison<-dcast(counts,Collection+Pattern~Background,value.var='DisplayCandidates',fill=0)
  fwrite(comparison,file.path(out,'display_count_comparison.csv'))
  fwrite(all[Background=='Full_tested_10531' & FDR<.05 & Count>=5],file.path(out,'full10531_significant_ORA.csv'))
  key<-data.table(Pattern=c('P3','P3','P5','P5','P5','P6','P6'),Pathway=c('HALLMARK_OXIDATIVE_PHOSPHORYLATION',
    'REACTOME_RESPIRATORY_ELECTRON_TRANSPORT','HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT','REACTOME_CELL_CYCLE',
    'REACTOME_O_LINKED_GLYCOSYLATION_OF_MUCINS','HALLMARK_ESTROGEN_RESPONSE_LATE'))
  fwrite(merge(key,all,by=c('Pattern','Pathway')),file.path(out,'key_function_comparison.csv'))
  totals<-all[,.(Tests=.N,Significant= sum(FDR<.05 & Count>=5)),by=.(Background,Pattern)]
  fwrite(totals,file.path(out,'pattern_total_comparison.csv'));print(dcast(totals,Pattern~Background,value.var='Significant'))
  fwrite(data.table(Source=paths,MD5=hashes),file.path(out,'input_provenance.csv'));stopifnot(identical(hashes,unname(tools::md5sum(paths))))
  writeLines(capture.output(sessionInfo()),file.path(out,'session_info.txt'))
  jsonlite::write_json(list(Genes=842,Patterns=6,TestedGeneUniverse=10531,NewClustering=FALSE,
    SourceHashesUnchanged=TRUE,TotalTestRows=nrow(all),NewTestRows=nrow(res),BaselineReused=TRUE),file.path(out,'analysis_diagnostics.json'),auto_unbox=TRUE,pretty=TRUE)
}
if(sys.nframe()==0)prepare_ORA_background_v21()
