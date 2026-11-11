# Additive subtype comparisons; source populations and threshold calls are frozen.
prepare_fig4_cnv_statistics <- function(root) {
  library(data.table)
  out<-file.path(root,'4-Multiomic/refinements/fig4_cnv_tests_20261030')
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  base<-'4-Multiomic/refinements/fig4_cnv_20261030/'
  paths<-c(paste0(base,c('burden_samples','gene_annotation_audit','cytoband_layout','cytoband_frequency','input_hashes'),'.csv'),
    '0-Data/raw/2022Dong/WES_CNV_GISTIC2_raw.txt',
    '4-Multiomic/publication/tables/step4_driver_cnv_sample_long_2022dong.csv',
    '4-Multiomic/refinements/fig4_genomics_20261026/gene_order.csv',
    '4-Multiomic/scripts/prepare_fig4_cnv_statistics.R')
  hashes<-data.table(Path=paths,MD5=unname(tools::md5sum(file.path(root,paths))))
  if(file.exists(file.path(out,'input_hashes.csv'))&&identical(fread(file.path(out,'input_hashes.csv')),hashes)&&
    file.exists(file.path(out,'selected_drivers.csv')))return(invisible(out))
  d<-fread(file.path(root,paste0(base,'burden_samples.csv')));d[,Sample:=as.character(Sample)]
  ann<-fread(file.path(root,paste0(base,'gene_annotation_audit.csv')))[Eligible==TRUE]
  layout<-fread(file.path(root,paste0(base,'cytoband_layout.csv')))
  raw<-fread(file.path(root,'0-Data/raw/2022Dong/WES_CNV_GISTIC2_raw.txt'))
  ids<-d$Sample;gene<-trimws(sub('\\|.*','',raw[[1]]))
  valid<-!is.na(gene)&gene!=''&!grepl('^[0-9]{1,2}-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$',gene,ignore.case=TRUE)
  x<-as.matrix(raw[valid,..ids]);cnv<-rowsum(x,gene[valid],reorder=FALSE)/rowsum(1L*!is.na(x),gene[valid],reorder=FALSE)
  indices<-match(ann$Gene,rownames(cnv));stopifnot(!anyNA(indices),all(is.finite(cnv)))
  bands<-ann$Cytoband;ng<-table(bands)
  rates<-list();n<-0L
  for(direction in c('Gain','Loss')) {
    calls<-if(direction=='Gain')cnv[indices,,drop=FALSE]>=.3 else cnv[indices,,drop=FALSE]<= -.3
    m<-rowsum(1L*calls,bands,reorder=FALSE);m<-m/as.integer(ng[rownames(m)])
    colnames(m)<-ids
    z<-as.data.table(as.table(m));setnames(z,c('Cytoband','Sample','AlteredGeneFraction'))
    z[,`:=`(Sample=as.character(Sample),Cytoband=as.character(Cytoband),Direction=direction)]
    z<-merge(z,d[,.(Sample,FinalSubtype)],by='Sample')
    n<-n+1L;rates[[n]]<-z
  }
  rates<-rbindlist(rates)
  region<-rates[,{
    invariant<-length(unique(AlteredGeneFraction))==1L
    fit<-if(!invariant)kruskal.test(AlteredGeneFraction~FinalSubtype) else NULL
    means<-tapply(AlteredGeneFraction,FinalSubtype,mean)
    .(N=.N,C1N=sum(FinalSubtype=='C1'),C2N=sum(FinalSubtype=='C2'),C3N=sum(FinalSubtype=='C3'),
      PValue=if(invariant)1 else fit$p.value,Statistic=if(invariant)0 else unname(fit$statistic),
      TestStatus=if(invariant)'Invariant; assigned P=1' else 'Kruskal-Wallis with tie correction',
      C1Mean=unname(means['C1']),C2Mean=unname(means['C2']),C3Mean=unname(means['C3']),
      MaxMeanDifference=max(means)-min(means))
  },by=.(Cytoband,Direction)]
  region[,FDR:=p.adjust(PValue,'BH')]
  region<-merge(region,layout,by='Cytoband');setorder(region,Direction,X)
  # Check that inference uses the same sample-level quantity as the plotted means.
  old<-fread(file.path(root,paste0(base,'cytoband_frequency.csv')))
  means<-rates[,.(Frequency=mean(AlteredGeneFraction)),by=.(Cytoband,Direction,FinalSubtype)]
  expected<-melt(old,id.vars=c('Cytoband','FinalSubtype'),measure.vars=c('GainFrequency','LossFrequency'),
    variable.name='Direction',value.name='Expected');expected[,Direction:=sub('Frequency','',Direction)]
  check<-merge(means,expected,by=c('Cytoband','Direction','FinalSubtype'))
  stopifnot(nrow(region)==1534L,nrow(check)==4602L,max(abs(check$Frequency-check$Expected))<1e-12)
  driver<-fread(file.path(root,'4-Multiomic/publication/tables/step4_driver_cnv_sample_long_2022dong.csv'))
  driver[,Sample:=as.character(Sample)];driver<-driver[Sample%in%ids]
  stopifnot(nrow(driver)==219L*15L,!anyDuplicated(driver[,.(Sample,Gene)]))
  joined<-d[match(driver$Sample,Sample)];stopifnot(identical(driver$FinalSubtype,joined$FinalSubtype))
  drivers<-rbindlist(lapply(c('Gain','Loss'),function(direction) {
    driver[,{
      called<-Call==direction
      a<-table(factor(FinalSubtype,levels=c('C1','C2','C3')),factor(called,levels=c(TRUE,FALSE)))
      fit<-if(any(called)&&!all(called))fisher.test(a,workspace=2e6) else NULL
      .(PValue=if(is.null(fit))1 else fit$p.value,Altered=sum(called),N=.N,
        C1Frequency=mean(called[FinalSubtype=='C1']),C2Frequency=mean(called[FinalSubtype=='C2']),
        C3Frequency=mean(called[FinalSubtype=='C3']))
    },by=Gene][,Direction:=direction]
  }))
  drivers[,FDR:=p.adjust(PValue,'BH')]
  drivers<-merge(drivers,ann[,.(Gene,Cytoband)],by='Gene',all.x=TRUE)
  drivers<-merge(drivers,region[,.(Cytoband,Direction,RegionPValue=PValue,RegionFDR=FDR)],by=c('Cytoband','Direction'),all.x=TRUE)
  drivers[,SupportedLink:=!is.na(RegionFDR)&RegionFDR<.05&FDR<.05]
  order<-fread(file.path(root,'4-Multiomic/refinements/fig4_genomics_20261026/gene_order.csv'))
  selected<-merge(drivers[,.(AnyRecordedCNV=any(Altered>0),GeneSignificant=any(FDR<.05),
    RegionAndGeneSupported=any(SupportedLink),MinGeneFDR=min(FDR),MinRegionFDR=min(RegionFDR,na.rm=TRUE)),by=Gene],order[,.(Gene,DisplayOrder)],by='Gene')
  selected[,Selected:=RegionAndGeneSupported]
  setorder(selected,DisplayOrder)
  freq<-driver[,.(N=.N,Gain=sum(Call=='Gain'),Loss=sum(Call=='Loss')),by=.(Gene,FinalSubtype)]
  freq[,`:=`(GainFraction=Gain/N,LossFraction=Loss/N)]
  write<-function(x,name)fwrite(x,file.path(out,paste0(name,'.csv')),na='NA')
  write(rates,'cytoband_sample_fractions');write(region,'cytoband_tests');write(drivers,'driver_tests')
  write(selected,'selected_drivers');write(freq,'driver_frequency');write(hashes,'input_hashes')
  summary<-data.table(RegionTests=nrow(region),SignificantRegionGain=sum(region$FDR<.05&region$Direction=='Gain'),
    SignificantRegionLoss=sum(region$FDR<.05&region$Direction=='Loss'),DriverTests=nrow(drivers),
    DriversWithAnyCall=sum(selected$AnyRecordedCNV),SelectedDriverGenes=sum(selected$Selected))
  write(summary,'summary');print(summary);print(selected);print(drivers[FDR<.05])
  invisible(out)
}
