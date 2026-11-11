# Scoped owner analysis for Figure4 frequency association sensitivity.
# Uses locked assay-roster calls; leaves the primary estimates and FDR unchanged.
prepare_fig4_frequency_sensitivity <- function(root) {
  out <- file.path(root,'4-Multiomic/refinements/fig4_frequency_20261030')
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  source_paths <- c('4-Multiomic/publication/tables/step4_driver_snv_indel_sample_long.csv',
    '4-Multiomic/publication/tables/step4_driver_snv_indel_study_stratified_tests.csv',
    'governance/data/subtype_master_table.csv','4-Multiomic/scripts/prepare_fig4_frequency_sensitivity.R')
  d <- data.table::fread(file.path(root,source_paths[1]))[SubtypeClass=='Main']
  m <- data.table::fread(file.path(root,source_paths[3]))
  z <- m[match(d$Sample,Sample)]
  stopifnot(!anyNA(z$Sample),identical(d$Study,z$Study),identical(d$FinalSubtype,z$FinalSubtype),
    !anyDuplicated(d[,.(Sample,Gene)]))
  primary <- data.table::fread(file.path(root,source_paths[2]))
  d <- d[Study %in% c('2022Dong','TCGA-CHOL')]
  rows <- data.table::rbindlist(lapply(c('C1','C2','C3'),function(st) {
    data.table::rbindlist(lapply(unique(d$Gene),function(g) {
      x <- d[Gene==g];target <- x$FinalSubtype==st
      arr <- table(factor(target,levels=c(TRUE,FALSE)),factor(x$Mutated,levels=c(TRUE,FALSE)),factor(x$Study))
      fit <- tryCatch(mantelhaen.test(arr,correct=FALSE),error=function(e)NULL)
      data.table::data.table(Subtype=st,Gene=g,TargetN=sum(target),OtherN=sum(!target),
        TargetAltered=sum(x$Mutated[target]),OtherAltered=sum(x$Mutated[!target]),
        CMHOddsRatio=if(is.null(fit))NA_real_ else unname(fit$estimate),
        CMHLower95=if(is.null(fit)||is.null(fit$conf.int))NA_real_ else fit$conf.int[1],
        CMHUpper95=if(is.null(fit)||is.null(fit$conf.int))NA_real_ else fit$conf.int[2],
        CMHPValue=if(is.null(fit))NA_real_ else fit$p.value)
    }))
  }))
  rows[,CMHFDR:=p.adjust(CMHPValue,'BH',n=15),by=Subtype]
  primary <- primary[,.(Subtype,Gene,PrimaryOR=CMHOddsRatio,PrimaryFDR=CMHFDR)]
  rows <- merge(rows,primary,by=c('Subtype','Gene'))
  # Infinite OR is a positive boundary association, not a direction reversal.
  rows[,DirectionConsistent:=!is.na(CMHOddsRatio)&!is.na(PrimaryOR)&sign(CMHOddsRatio-1)==sign(PrimaryOR-1)]
  data.table::fwrite(rows,file.path(out,'dong_tcga_frequency_associations.csv'),na='NA')
  freq <- d[,.(N=.N,Altered=sum(Mutated),Fraction=mean(Mutated)),by=.(Gene,FinalSubtype)]
  data.table::fwrite(freq,file.path(out,'dong_tcga_frequency.csv'))
  data.table::fwrite(data.table::data.table(Path=source_paths,MD5=unname(tools::md5sum(file.path(root,source_paths)))),file.path(out,'input_hashes.csv'))
  summary <- rows[,.(Comparisons=.N,PrimarySignificant=sum(PrimaryFDR<.05,na.rm=TRUE),
    SensitivitySignificant=sum(CMHFDR<.05,na.rm=TRUE),
    PrimarySignificantSameDirection=sum(PrimaryFDR<.05&DirectionConsistent,na.rm=TRUE)),by=Subtype]
  data.table::fwrite(summary,file.path(out,'summary.csv'));print(summary)
  invisible(out)
}
