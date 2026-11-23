suppressPackageStartupMessages(library(data.table))
prepare_soft_patterns_v20<-function(){
  out<-'2-TranscriptomicCharacterization/runs/fig2_soft_patterns_20261030_v20'
  old<-'2-TranscriptomicCharacterization/runs/fig2_six_patterns_20261030_v19'
  paths<-c(file.path(out,'benchmark_model.rds'),file.path(old,'reference_membership_used.csv.gz'),
    '1-Clustering/data/top_variable_genes.csv',file.path(out,'CONTRACT.md'),
    '2-TranscriptomicCharacterization/scripts/prepare_fig2_soft_patterns_v20.R')
  hashes<-unname(tools::md5sum(paths));b<-readRDS(paths[1]);source<-b$source_model
  sets<-fread(paths[2]);hvg<-fread(paths[3])$Gene
  oldmeta<-source$membership;z<-source$z;allmeta<-source$all
  tracks<-list();profiles<-list();orall<-list();coverage<-list();i<-1L;j<-1L
  for(k in c(6L,8L)){
    f<-b$models[[paste0('Fuzzy_',k)]];center<-f$centers
    angle<-(atan2(sqrt(3)/2*(center[,2]-center[,3]),center[,1]-(center[,2]+center[,3])/2)+pi/6)%%(2*pi)
    raworder<-order(angle);map<-setNames(paste0('P',seq_len(k)),as.character(raworder))
    rawcluster<-f$assignment
    mem<-oldmeta[,.(Gene,Center,SD,WGCNAModule,ModuleClass,kME)]
    mem[,`:=`(RawCluster=as.integer(rawcluster[Gene]),Pattern=unname(map[as.character(rawcluster[Gene])]))]
    mem[,Affinity:=apply(f$membership[Gene,,drop=FALSE],1,max)]
    mem[,AffinityGap:=apply(f$membership[Gene,,drop=FALSE],1,function(x){x<-sort(x,decreasing=TRUE);x[1]-x[2]})]
    mem[,LowAffinity:=Affinity<.6]
    pcounts<-mem[,.(Genes=.N,MeanAffinity=mean(Affinity),LowAffinityGenes=sum(LowAffinity)),by=.(Pattern,RawCluster)]
    pcounts[,PatternOrder:=as.integer(sub('P','',Pattern))];setorder(pcounts,PatternOrder)
    for(h in 1:nrow(pcounts)){cl<-pcounts$RawCluster[h];c<-center[cl,];names(c)<-c('C1','C2','C3');oo<-order(c,decreasing=TRUE)
      # Hue is only a descriptive centroid orientation; no equivalent-high test.
      context<-if(c[oo[1]]-c[oo[2]]<.5)paste(sort(names(c)[oo[1:2]]),collapse='_')else names(c)[oo[1]]
      pcounts[h,`:=`(ShapeC1=c[1],ShapeC2=c[2],ShapeC3=c[3],Context=context,
        ShapeOrder=paste(names(c)[oo],collapse=' > '))]
    }
    gsorder<-unlist(lapply(pcounts$Pattern,function(p){gs<-mem[Pattern==p,Gene];gs[hclust(dist(z[gs,source$main$Sample,drop=FALSE]),method='ward.D2')$order]}),use.names=FALSE)
    mem<-mem[match(gsorder,Gene)];mem[,RowOrder:=seq_len(.N)]
    counts<-pcounts[,.(Pattern,Genes,RobustGenes=NA_integer_,PatternOrder,Context,MeanAffinity,LowAffinityGenes,ShapeOrder,ShapeC1,ShapeC2,ShapeC3)]
    # Stable matching is a cluster-level diagnostic; no per-gene evidence gets hidden.
    jj<-fread(file.path(out,'benchmark_LOSO_Jaccard.csv'))[Method=='Fuzzy' & K==k]
    counts[,MeanLOSOJaccard:=vapply(Pattern,function(p){cl<-pcounts[Pattern==p,RawCluster];mean(jj[Cluster==cl,Jaccard])},numeric(1))]
    local<-file.path(out,paste0('K',k));dir.create(local,showWarnings=FALSE)
    fwrite(mem,file.path(local,'gene_pattern_membership.csv'));fwrite(counts,file.path(local,'pattern_counts.csv'))
    fwrite(data.table(Gene=rownames(f$membership),f$membership),file.path(local,'fuzzy_membership.csv'))
    savescores<-rbindlist(lapply(counts$Pattern,function(p){gs<-mem[Pattern==p,Gene];data.table(Pattern=p,Genes=length(gs),
      Sample=allmeta$Sample,Study=allmeta$Study,FinalSubtype=allmeta$FinalSubtype,SubtypeClass=allmeta$SubtypeClass,
      Score=as.numeric(colMeans(z[gs,allmeta$Sample,drop=FALSE])))}))
    fwrite(savescores,file.path(local,'pattern_scores_all_states.csv'))
    summary<-savescores[,.(N=.N,Mean=mean(Score),Median=median(Score),Q25=quantile(Score,.25),Q75=quantile(Score,.75)),by=.(Pattern,Genes,FinalSubtype,SubtypeClass)]
    fwrite(summary,file.path(local,'pattern_state_summary.csv'))
    saveRDS(list(z=z,main=source$main,all=allmeta,membership=mem,counts=counts,pattern_order=counts$Pattern),file.path(local,'frozen_display_model.rds'))
    for(cn in unique(sets$Collection)){
      ss<-sets[Collection==cn];univ<-intersect(hvg,unique(ss$Gene));bg<-ss[Gene %in% univ]
      sizes<-bg[,.(SetSize=uniqueN(Gene)),by=.(Collection,Pathway)][SetSize>=10 & SetSize<=500]
      bg<-merge(bg,sizes,by=c('Collection','Pathway'))
      for(p in counts$Pattern){fg<-intersect(mem[Pattern==p,Gene],univ)
        hit<-bg[Gene %in% fg,.(Count=uniqueN(Gene),HitGenes=paste(sort(unique(Gene)),collapse='/')),by=.(Collection,Pathway)]
        dd<-merge(sizes,hit,by=c('Collection','Pathway'),all.x=TRUE);dd[is.na(Count),`:=`(Count=0L,HitGenes='')]
        dd[,`:=`(K=k,Pattern=p,UniverseGenes=length(univ),SelectedGenes=length(fg))]
        dd[,PValue:=phyper(Count-1,SetSize,UniverseGenes-SetSize,SelectedGenes,lower.tail=FALSE)]
        orall[[i]]<-dd;i<-i+1L
        coverage[[j]]<-data.table(K=k,Pattern=p,Collection=cn,Genes=mem[Pattern==p,.N],MappedGenes=length(fg),UniverseGenes=length(univ),Terms=nrow(sizes));j<-j+1L
      }
    }
    tracks[[as.character(k)]]<-list(counts=counts,model=local)
    print(counts)
  }
  ora<-rbindlist(orall);ora[,FDR:=p.adjust(PValue,'BH'),by=.(K,Collection)];ora[,FamilySize:=.N,by=.(K,Collection)]
  setorder(ora,K,Collection,Pattern,FDR,-Count,Pathway)
  fwrite(ora,file.path(out,'actual_pattern_ORA_all_tests.csv.gz'));fwrite(rbindlist(coverage),file.path(out,'ORA_coverage.csv'))
  for(k in c(6,8))fwrite(ora[K==k & FDR<.05 & Count>=5],file.path(out,paste0('K',k,'/primary_significant_ORA.csv')))
  fwrite(data.table(Source=paths,MD5=hashes),file.path(out,'preparation_input_provenance.csv'));stopifnot(identical(hashes,unname(tools::md5sum(paths))))
  writeLines(capture.output(sessionInfo()),file.path(out,'session_info.txt'))
  cat('Prepared fixed6 and benchmark-selected8 profile partitions and actual-memberORA.\n')
}
if(sys.nframe()==0)prepare_soft_patterns_v20()
