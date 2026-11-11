# Figure 4 genomics extension. Defines functions only; never refreshes locked Step4.
# New inference belongs here, not in the Step8 renderer.
suppressPackageStartupMessages(library(data.table))

fig4_pair_test <- function(a, b, strata, min_events=5L) {
  stopifnot(length(a)==length(b), length(a)==length(strata),
            !anyNA(a), !anyNA(b), !anyNA(strata))
  a <- as.logical(a); b <- as.logical(b)
  out <- list(N=length(a), Gene1Altered=sum(a), Gene2Altered=sum(b),
              Both=sum(a & b), Gene1Only=sum(a & !b), Gene2Only=sum(!a & b),
              Neither=sum(!a & !b), InformativeStrata=0L,
              OddsRatio=NA_real_, Lower95=NA_real_, Upper95=NA_real_,
              PValue=NA_real_, Status='Low_event_support')
  # Marginal support rule is fixed before testing. Zero overlap alone is NOT a gate.
  if (min(sum(a),sum(b),sum(!a),sum(!b)) < min_events) return(out)
  tabs <- lapply(split(seq_along(a),strata), function(i)
    table(factor(a[i],levels=c(TRUE,FALSE)),factor(b[i],levels=c(TRUE,FALSE))))
  tabs <- Filter(function(x) all(rowSums(x)>0) && all(colSums(x)>0), tabs)
  out$InformativeStrata <- length(tabs)
  if (!length(tabs)) {out$Status <- 'No_informative_stratum'; return(out)}
  arr <- array(unlist(tabs),dim=c(2L,2L,length(tabs)))
  test <- if (length(tabs)==1L) fisher.test(tabs[[1]],alternative='two.sided',conf.level=.95)
          else mantelhaen.test(arr,exact=TRUE,alternative='two.sided',conf.level=.95)
  out$OddsRatio <- unname(test$estimate)
  out$Lower95 <- unname(test$conf.int[1]); out$Upper95 <- unname(test$conf.int[2])
  out$PValue <- test$p.value; out$Status <- 'Tested'
  out
}

fig4_pair_family <- function(mat, roster, genes, strata_col, min_events=5L) {
  pairs <- combn(genes,2,simplify=FALSE)
  ans <- rbindlist(lapply(c('C1','C2','C3'),function(st) {
    ids <- which(roster$FinalSubtype==st)
    rbindlist(lapply(pairs,function(g)
      c(list(Subtype=st,Gene1=g[1],Gene2=g[2]),
        fig4_pair_test(mat[ids,g[1]],mat[ids,g[2]],roster[[strata_col]][ids],min_events))))
  }))
  # One prespecified family: all 105 gene pairs x 3 subtypes. Untested slots get p=1.
  ans[, FDR := p.adjust(ifelse(is.na(PValue),1,PValue),method='BH')]
  ans[Status!='Tested',FDR:=NA_real_]
  ans[, FamilySize := length(pairs)*3L]
  ans
}

prepare_fig4_genomics <- function(root, force=FALSE) {
  root <- normalizePath(root,mustWork=TRUE)
  out_rel <- '4-Multiomic/refinements/fig4_genomics_20261026'
  out <- file.path(root,out_rel); dir.create(out,recursive=TRUE,showWarnings=FALSE)
  t4 <- '4-Multiomic/publication/tables/'
  inputs <- c('governance/data/subtype_master_table.csv','0-Data/qc/0-Data-combined_maf_data.csv',
    paste0(t4,c('step4_driver_snv_indel_sample_long.csv','step4_driver_snv_indel_frequency_main_subtypes.csv',
      'step4_driver_snv_indel_study_stratified_tests.csv','step4_nonsynonymous_variant_count_by_sample.csv',
      'step4_driver_cnv_sample_long_2022dong.csv','step4_driver_cnv_frequency_main_subtypes_2022dong.csv')),
    '4-Multiomic/scripts/prepare_fig4_genomics.R')
  hashes <- data.table(Path=inputs,MD5=unname(tools::md5sum(file.path(root,inputs))))
  stopifnot(!anyNA(hashes$MD5))
  required <- c('mutation_display_long.csv','sample_order.csv','gene_order.csv',
    'mutation_frequency.csv','cnv_frequency.csv','pairwise_associations.csv',
    'pairwise_burden_sensitivity.csv','source_coverage.csv','technical_checks.csv')
  cache <- file.path(out,'input_hashes.csv')
  if (!force && file.exists(cache) && all(file.exists(file.path(out,required))) &&
      identical(hashes,fread(cache))) {
    saved <- fread(file.path(out,'output_hashes.csv'))
    stopifnot(identical(saved$MD5,unname(tools::md5sum(file.path(out,saved$Path)))))
    message('Figure4 genomics: verified existing owner tables'); return(invisible(out_rel))
  }
  read <- function(p) {
    d <- fread(file.path(root,p))
    if ('Sample' %in% names(d)) d[,Sample:=as.character(Sample)]
    d
  }
  master <- read(inputs[1]); stopifnot(!anyDuplicated(master$Sample))
  grid <- read(inputs[3]); grid <- grid[SubtypeClass=='Main']
  canon <- master[match(grid$Sample,Sample)]
  stopifnot(!anyNA(canon$Sample),identical(grid$FinalSubtype,canon$FinalSubtype),
            identical(grid$Study,canon$Study),!anyDuplicated(grid[,.(Sample,Gene)]),!anyNA(grid$Mutated))
  roster <- unique(grid[,.(Sample,Study,FinalSubtype,SubtypeClass)])
  stopifnot(nrow(roster)==308L,identical(as.integer(table(factor(roster$FinalSubtype,levels=c('C1','C2','C3')))),c(63L,153L,92L)))
  freq <- read(inputs[4]); tests <- read(inputs[5])
  order <- grid[,.(Altered=sum(Mutated)),by=Gene][order(-Altered,Gene)]
  genes <- order$Gene; stopifnot(length(genes)==15L,nrow(grid)==nrow(roster)*length(genes))
  mat <- matrix(FALSE,nrow=nrow(roster),ncol=length(genes),dimnames=list(roster$Sample,genes))
  mat[cbind(match(grid$Sample,roster$Sample),match(grid$Gene,genes))] <- grid$Mutated
  rebuilt <- grid[,.(N=.N,Altered=sum(Mutated),Fraction=mean(Mutated)),by=.(Gene,FinalSubtype)]
  chk <- merge(rebuilt,freq,by=c('Gene','FinalSubtype'),suffixes=c('.new','.locked'))
  stopifnot(nrow(chk)==45L,all(chk$N.new==chk$N.locked),all(chk$Altered.new==chk$Altered.locked))
  freq <- merge(freq,tests[,.(FinalSubtype=Subtype,Gene,CMHOddsRatio,CMHLower95,CMHUpper95,CMHFDR)],
                by=c('Gene','FinalSubtype'),all.x=TRUE)
  burden <- read(inputs[6]); roster <- merge(roster,burden[,.(Sample,NonsynonymousVariantCount)],by='Sample',sort=FALSE)
  mat <- mat[roster$Sample,,drop=FALSE]

  maf <- fread(file.path(root,inputs[2]),colClasses=c(Tumor_Sample_Barcode='character'))
  maf[,Study:=c(Dong2022='2022Dong',Deng2023='2023Deng',TCGA='TCGA-CHOL')[Dataset]]
  maf[,JoinKey:=ifelse(Study=='TCGA-CHOL',substr(Tumor_Sample_Barcode,1,12),Tumor_Sample_Barcode)]
  keys <- copy(roster);keys[,JoinKey:=ifelse(Study=='TCGA-CHOL',substr(Sample,1,12),Sample)]
  stopifnot(!anyDuplicated(keys[,.(Study,JoinKey)]))
  coverage <- maf[,.(SourceGenes=uniqueN(Hugo_Symbol),SourceBarcodes=uniqueN(Tumor_Sample_Barcode),
                    SourceRows=.N),by=Study]
  maf <- merge(maf,keys[,.(Study,JoinKey,Sample)],by=c('Study','JoinKey'))
  ns <- maf[!is.na(Variant_Classification) & Variant_Classification!='Silent' & !is.na(Hugo_Symbol) & nzchar(Hugo_Symbol)]
  ns[,VariantKey:=paste(Hugo_Symbol,Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2,sep=':')]
  re_count <- ns[,.(Count=uniqueN(VariantKey)),by=Sample]
  counts <- merge(roster,re_count,by='Sample',all.x=TRUE);counts[is.na(Count),Count:=0L]
  stopifnot(all(counts$Count==counts$NonsynonymousVariantCount))
  background <- ns[!Hugo_Symbol %in% genes,.(BackgroundVariantCount=uniqueN(VariantKey)),by=Sample]
  roster <- merge(roster,background,by='Sample',all.x=TRUE,sort=FALSE)
  roster[is.na(BackgroundVariantCount),BackgroundVariantCount:=0L]
  mat <- mat[roster$Sample,,drop=FALSE]
  mapping <- c(Missense_Mutation='Missense',Nonsense_Mutation='Nonsense',Frame_Shift_Del='Frameshift deletion',
    Frame_Shift_Ins='Frameshift insertion',In_Frame_Del='In-frame deletion',In_Frame_Ins='In-frame insertion',
    Splice_Site='Splice',Splice_Region='Splice',Nonstop_Mutation='Other non-silent',
    Translation_Start_Site='Other non-silent',Multi_Hit='Multi-hit')
  ev <- ns[Hugo_Symbol %in% genes]
  stopifnot(all(ev$Variant_Classification %in% names(mapping)))
  classes <- ev[,.(MutationClass=if(uniqueN(VariantKey)>1L || any(Variant_Classification=='Multi_Hit')) 'Multi-hit'
    else unname(mapping[Variant_Classification[1]]),
    SourceClasses=paste(sort(unique(Variant_Classification)),collapse=';'),VariantCount=uniqueN(VariantKey)),
    by=.(Sample,Gene=Hugo_Symbol)]
  display <- merge(grid,classes,by=c('Sample','Gene'),all.x=TRUE,sort=FALSE)
  stopifnot(all(!is.na(display$MutationClass)==display$Mutated))
  display[is.na(MutationClass),`:=`(MutationClass='No recorded variant',SourceClasses='',VariantCount=0L)]
  # Lexicographic waterfall within subtype; all assayed samples, including negatives.
  ord <- do.call(base::order,c(list(match(roster$FinalSubtype,c('C1','C2','C3'))),
    lapply(genes,function(g)-as.integer(mat[,g])),list(roster$Study,roster$Sample)))
  roster <- roster[ord];roster[,DisplayOrder:=seq_len(.N)];mat <- mat[roster$Sample,,drop=FALSE]
  order[,DisplayOrder:=seq_len(.N)]
  coverage <- merge(coverage,roster[,.(MainAssayedN=.N,MainWithRecordedVariant=sum(NonsynonymousVariantCount>0)),by=Study],by='Study')
  coverage[,Interpretation:=ifelse(Study=='2023Deng','Restricted curated gene set; excluded from burden sensitivity',
    'Broad curated gene set; gene-level callable territory unavailable')]
  pairs <- fig4_pair_family(mat,roster,genes,'Study')
  # Burden sensitivity: only broad-record sources, burden excludes all 15 tested drivers.
  sens <- roster[Study!='2023Deng']
  sens[,BurdenGroup:=ifelse(BackgroundVariantCount>median(BackgroundVariantCount),'High','Low'),by=.(Study,FinalSubtype)]
  sens[,BurdenStratum:=paste(Study,BurdenGroup,sep=' / ')]
  sp <- fig4_pair_family(mat[sens$Sample,,drop=FALSE],sens,genes,'BurdenStratum')
  cnv <- read(inputs[8]); cnv_grid <- read(inputs[7]);cnv_grid <- cnv_grid[SubtypeClass=='Main']
  canon <- master[match(cnv_grid$Sample,Sample)]
  stopifnot(!anyDuplicated(cnv_grid[,.(Sample,Gene)]),identical(cnv_grid$FinalSubtype,canon$FinalSubtype),
            all(cnv_grid$Study=='2022Dong'),!anyNA(cnv_grid$ContinuousCNV),
            all(cnv_grid$Call==ifelse(cnv_grid$ContinuousCNV>=.3,'Gain',ifelse(cnv_grid$ContinuousCNV<= -.3,'Loss','Neutral'))))
  cf <- cnv_grid[,.(N=.N,Gain=sum(Call=='Gain'),Loss=sum(Call=='Loss')),by=.(Gene,FinalSubtype)]
  chk <- merge(cf,cnv,by=c('Gene','FinalSubtype'),suffixes=c('.new','.locked'))
  stopifnot(nrow(chk)==45L,all(chk$N.new==chk$N.locked),all(chk$Gain.new==chk$Gain.locked),all(chk$Loss.new==chk$Loss.locked))
  write <- function(d,n) fwrite(d,file.path(out,n),na='NA')
  write(display,'mutation_display_long.csv');write(roster,'sample_order.csv');write(order,'gene_order.csv')
  write(freq,'mutation_frequency.csv');write(cnv,'cnv_frequency.csv');write(pairs,'pairwise_associations.csv')
  write(sp,'pairwise_burden_sensitivity.csv');write(sens,'burden_sensitivity_roster.csv');write(coverage,'source_coverage.csv')
  checks <- data.table(Check=c('Canonical labels and unique sample-gene keys','All 308 DNA-assayed principal samples retained',
    'Mutation classification reproduces locked positives','Variant counts reproduce locked counts',
    'Mutation and CNV frequencies reproduce locked tables','CNV thresholds preserved at +/-0.3','315 prespecified tests per pair family'),
    Status='PASS')
  stopifnot(nrow(pairs)==315L,nrow(sp)==315L)
  write(checks,'technical_checks.csv')
  results <- rbindlist(list(pairs[,.(Analysis='Study-stratified primary',Subtype,N,Status,FDR)],
    sp[,.(Analysis='Study + background-burden sensitivity',Subtype,N,Status,FDR)]))
  summary <- results[,.(N=unique(N),TotalPairs=.N,TestedPairs=sum(Status=='Tested'),
    FDR_lt_0_05=sum(FDR<.05,na.rm=TRUE)),by=.(Analysis,Subtype)]
  write(summary,'pairwise_summary.csv');print(summary)
  write(hashes,'input_hashes.csv')
  files <- c(required,'burden_sensitivity_roster.csv','pairwise_summary.csv')
  write(data.table(Path=files,MD5=unname(tools::md5sum(file.path(out,files)))),'output_hashes.csv')
  invisible(out_rel)
}
