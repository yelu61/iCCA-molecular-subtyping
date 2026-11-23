suppressPackageStartupMessages({library(data.table);library(limma);library(cluster);library(e1071)})
benchmark_pattern_shapes_v20<-function(){
  out<-'2-TranscriptomicCharacterization/runs/fig2_soft_patterns_20261030_v20'
  src<-c('2-TranscriptomicCharacterization/runs/fig2_six_patterns_20261030_v19/frozen_display_model.rds',
    '2-TranscriptomicCharacterization/runs/fig2_six_patterns_20261030_v19/gene_pattern_membership.csv',
    '2-TranscriptomicCharacterization/runs/fig2_six_patterns_20261030_v19/CONTRACT.md',
    file.path(out,'CONTRACT.md'),'2-TranscriptomicCharacterization/scripts/benchmark_fig2_pattern_shapes_v20.R')
  hashes<-unname(tools::md5sum(src));model<-readRDS(src[1]);z<-model$z[model$membership$Gene,model$main$Sample,drop=FALSE]
  meta<-model$main
  features<-function(z,m){m<-copy(m);m[,State:=factor(FinalSubtype,levels=c('C1','C2','C3'))];m[,Study:=factor(Study)]
    x<-model.matrix(~0+State+Study,m);colnames(x)<-make.names(colnames(x));f<-lmFit(z,x)
    av<-f$coefficients[,c('StateC1','StateC2','StateC3'),drop=FALSE];colnames(av)<-c('C1','C2','C3')
    centered<-av-rowMeans(av);norm<-sqrt(rowSums(centered^2));stopifnot(all(norm>0),all(is.finite(av)))
    list(unit=centered/norm,centered=centered,coefficients=av)
  }
  feat<-features(z,meta);x<-feat$unit;ds<-dist(x);dz<-dist(z)
  loo<-lapply(sort(unique(meta$Study)),function(st){mm<-meta[Study!=st];features(z[,mm$Sample,drop=FALSE],mm)$unit});names(loo)<-sort(unique(meta$Study))
  ari<-function(a,b){tt<-table(a,b);ch<-function(q)q*(q-1)/2;aa<-sum(ch(rowSums(tt)));bb<-sum(ch(colSums(tt)));e<-aa*bb/ch(length(a));(sum(ch(tt))-e)/((aa+bb)/2-e)}
  fitcluster<-function(x,k,method,seed){set.seed(seed)
    if(method=='Ward'){hc<-hclust(dist(x),method='ward.D2');a<-cutree(hc,k);centers<-t(vapply(seq_len(k),function(j)colMeans(x[a==j,,drop=FALSE]),numeric(ncol(x))));u<-NULL;objective<-sum(vapply(seq_len(nrow(x)),function(i)sum((x[i,]-centers[a[i],])^2),numeric(1)))}
    else if(method=='Kmeans'){f<-kmeans(x,centers=k,nstart=50,iter.max=300);a<-f$cluster;centers<-f$centers;u<-NULL;objective<-f$tot.withinss}
    else {best<-NULL;bestobj<-Inf
      for(s in 1:10){set.seed(seed+s);ff<-cmeans(x,centers=k,iter.max=300,m=2,control=list(reltol=1e-8));ob<-sum(ff$withinerror)
        if(ob<bestobj){best<-ff;bestobj<-ob}}
      a<-best$cluster;centers<-best$centers;u<-best$membership;objective<-bestobj}
    if(is.null(u)){dd<-sapply(seq_len(k),function(j)rowSums((x-matrix(centers[j,],nrow(x),ncol(x),byrow=TRUE))^2));u<-1/pmax(dd,1e-14);u<-u/rowSums(u)}
    names(a)<-rownames(x);rownames(centers)<-seq_len(k);colnames(u)<-seq_len(k);rownames(u)<-rownames(x)
    list(assignment=a,centers=centers,membership=u,objective=objective)
  }
  models<-list();metrics<-list();stab<-list();matches<-list();i<-1L;j<-1L;l<-1L
  for(method in c('Ward','Kmeans','Fuzzy'))for(k in 3:8){id<-paste(method,k,sep='_');f<-fitcluster(x,k,method,19000+k)
    models[[id]]<-f;sv<-silhouette(f$assignment,ds);rawsv<-silhouette(f$assignment,dz)
    for(st in names(loo)){g<-fitcluster(loo[[st]],k,method,19000+k)
      stab[[j]]<-data.table(Method=method,K=k,Study=st,ARI=ari(f$assignment,g$assignment));j<-j+1L
      for(a in seq_len(k)){ga<-names(f$assignment)[f$assignment==a]
        js<-vapply(seq_len(k),function(b){gb<-names(g$assignment)[g$assignment==b];length(intersect(ga,gb))/length(union(ga,gb))},numeric(1))
        matches[[l]]<-data.table(Method=method,K=k,Cluster=a,Study=st,Jaccard=max(js));l<-l+1L}
    }
    metrics[[i]]<-data.table(Method=method,K=k,MinGenes=min(table(f$assignment)),MaxGenes=max(table(f$assignment)),
      ShapeSilhouette=mean(sv[,'sil_width']),FullSampleSilhouette=mean(rawsv[,'sil_width']),
      MeanMaxAffinity=mean(apply(f$membership,1,max)),FractionAffinityBelow06=mean(apply(f$membership,1,max)<.6))
    i<-i+1L;cat('Compared ',id,'\n',sep='')
  }
  metrics<-rbindlist(metrics);stab<-rbindlist(stab);matches<-rbindlist(matches)
  ss<-stab[,.(MedianLOSOARI=median(ARI),MinLOSOARI=min(ARI)),by=.(Method,K)]
  jj<-matches[,.(MeanJ=mean(Jaccard)),by=.(Method,K,Cluster)][,.(MinPatternMeanJaccard=min(MeanJ)),by=.(Method,K)]
  metrics<-merge(merge(metrics,ss,by=c('Method','K')),jj,by=c('Method','K'))
  metrics[,Eligible:=MinGenes>=20 & MedianLOSOARI>=.65 & MinPatternMeanJaccard>=.60]
  if(any(metrics$Eligible)){best<-max(metrics[Eligible==TRUE,ShapeSilhouette]);cand<-metrics[Eligible==TRUE & ShapeSilhouette>=.8*best]
    setorder(cand,-K,-MinPatternMeanJaccard,-ShapeSilhouette,Method);chosen<-cand[1];status<-'Eligible'}
  else {cand<-metrics[MinGenes>=20];setorder(cand,-ShapeSilhouette,K,Method);chosen<-cand[1];status<-'No stability-eligible candidate; provisional'}
  metrics[,Selected:=Method==chosen$Method & K==chosen$K]
  setorder(metrics,K,Method)
  fwrite(metrics,file.path(out,'partition_benchmark.csv'));fwrite(stab,file.path(out,'benchmark_LOSO_ARI.csv'));fwrite(matches,file.path(out,'benchmark_LOSO_Jaccard.csv'))
  fwrite(data.table(Gene=rownames(x),x),file.path(out,'unit_state_profile_features.csv'))
  fwrite(data.table(Gene=rownames(x),feat$centered),file.path(out,'centered_state_effects.csv'))
  chosenid<-paste(chosen$Method,chosen$K,sep='_');primary<-models[[chosenid]]
  saveRDS(list(x=x,features=feat,z=z,meta=meta,loo=loo,models=models,chosen=chosen,chosenid=chosenid,
    primary=primary,status=status,source_model=model),file.path(out,'benchmark_model.rds'))
  old<-fread(src[2]);aa<-as.integer(factor(old$Pattern[match(rownames(x),old$Gene)]))
  oldmetric<-data.table(Method='v19 statistical categories plus unresolved',K=length(unique(aa)),
    ShapeSilhouette=mean(silhouette(aa,ds)[,'sil_width']),FullSampleSilhouette=mean(silhouette(aa,dz)[,'sil_width']),
    UnresolvedGenes=sum(old$Pattern=='Unresolved'))
  fwrite(oldmetric,file.path(out,'v19_partition_reference.csv'))
  fwrite(data.table(Source=src,MD5=hashes),file.path(out,'benchmark_input_provenance.csv'))
  stopifnot(identical(hashes,unname(tools::md5sum(src))))
  print(metrics);cat('SELECTED ',chosenid,' ',status,'\n',sep='')
}
if(sys.nframe()==0)benchmark_pattern_shapes_v20()
