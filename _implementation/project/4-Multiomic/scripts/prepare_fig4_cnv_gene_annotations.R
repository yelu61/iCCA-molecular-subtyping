# Gene membership annotations only; inferential region selection and statistics do not change.
prepare_fig4_cnv_gene_annotations <- function(root) {
  library(data.table)
  out<-file.path(root,'4-Multiomic/refinements/fig4_cnv_gene_annotations_20261030')
  ref<-file.path(out,'reference')
  stopifnot(file.exists(file.path(ref,'reference_manifest.csv')))
  display<-fread(file.path(root,'4-Multiomic/refinements/fig4_cnv_pairwise_20261030/pairwise_display.csv'))
  bands<-unique(display[,.(Cytoband,Direction,X)])
  features<-fread(file.path(root,'4-Multiomic/refinements/fig4_cnv_20261030/gene_annotation_audit.csv'))
  features<-features[Eligible==TRUE&Cytoband%in%bands$Cytoband,.(SourceGene=Gene,Cytoband)]
  raw<-fread(file.path(root,'0-Data/raw/2022Dong/WES_CNV_GISTIC2_raw.txt'),select=1:3)
  rawmap<-unique(data.table(SourceGene=trimws(sub('\\|.*','',raw[[1]])),SourceGeneID=as.character(raw[[2]])))
  features<-merge(features,rawmap,by='SourceGene')
  stopifnot(nrow(features)==633L,!anyDuplicated(features$SourceGene))
  refs<-fread(file.path(ref,'reference_manifest.csv'))
  stopifnot(identical(unname(tools::md5sum(file.path(ref,refs$File))),refs$MD5))
  records<-unlist(lapply(refs$File,function(f)jsonlite::fromJSON(file.path(ref,f),simplifyVector=FALSE)$reports),recursive=FALSE)
  text<-function(x)if(is.null(x)||!length(x))NA_character_ else paste(unique(unlist(x)),collapse=';')
  meta<-rbindlist(lapply(records,function(r) {
    g<-r$gene;retired<-identical(r$warning$gene_warning_code,'DISCONTINUED_GENE_ID')
    if(!retired)stopifnot(identical(g$tax_id,'9606'))
    maps<-vapply(g$map_locations,function(x)if(identical(x$map_type,'Cytogenetic'))x$map_value else NA_character_,character(1))
    maps<-maps[!is.na(maps)]
    annotations<-g$annotations
    h<-which(vapply(annotations,function(x)grepl('^GRCh38',x$assembly_name),logical(1)))
    a<-if(length(h))annotations[[h[1]]] else NULL
    loc<-if(!is.null(a)&&length(a$genomic_locations))a$genomic_locations[[1]] else NULL
    ids<-unique(c(text(g$gene_id),as.character(unlist(r$query))))
    data.table(SourceGeneID=ids,NCBIRecordID=text(g$gene_id),NCBIRecordStatus=if(retired)'Discontinued' else 'Active human gene',
      CurrentGeneID=if(retired)NA_character_ else text(g$gene_id),CurrentSymbol=text(g$symbol),
      GeneName=text(g$description),GeneType=text(g$type),CurrentCytoband=text(maps),
      HGNCID=text(g$nomenclature_authority$identifier),Aliases=text(g$synonyms),
      EnsemblGeneIDs=text(g$ensembl_gene_ids),ReferenceAssembly=text(a$assembly_name),
      ReferenceAnnotation=text(a$annotation_name),ReferenceAnnotationDate=text(a$annotation_release_date),
      ReferenceAccession=text(loc$genomic_accession_version),
      ReferenceBegin0=text(loc$genomic_range$begin),ReferenceEnd0=text(loc$genomic_range$end),
      GeneURL=paste0('https://www.ncbi.nlm.nih.gov/gene/',g$gene_id))
  }),fill=TRUE)
  meta<-unique(meta,by='SourceGeneID')
  all<-merge(features,meta,by='SourceGeneID',all.x=TRUE)
  all[,AnnotationStatus:=ifelse(as.numeric(SourceGeneID)<=0,'Source placeholder ID; no NCBI ID query',
    ifelse(!is.na(NCBIRecordStatus)&NCBIRecordStatus=='Discontinued','Discontinued NCBI Gene ID',
      ifelse(is.na(CurrentGeneID),'Not returned by NCBI query','NCBI human gene record')))]
  map_status<-function(source,current) {
    if(is.na(current))return('Unresolved')
    entries<-strsplit(current,';',fixed=TRUE)[[1]]
    if(source%in%entries)return('Exact match')
    parse<-function(x) {
      chr<-sub('[pq].*','',x);arm<-sub('^[0-9]+([pq]).*','\\1',x)
      value<-suppressWarnings(as.numeric(sub('^[0-9]+[pq]','',x)))
      list(chr=chr,pos=if(arm=='p')-value else 100+value)
    }
    sourcepos<-parse(source)
    for(entry in entries[grepl('-',entries,fixed=TRUE)]) {
      x<-strsplit(entry,'-',fixed=TRUE)[[1]]
      if(length(x)!=2)next
      if(grepl('^[pq]',x[2]))x[2]<-paste0(sub('[pq].*','',x[1]),x[2])
      lo<-parse(x[1]);hi<-parse(x[2])
      if(sourcepos$chr==lo$chr&&lo$chr==hi$chr&&all(is.finite(c(sourcepos$pos,lo$pos,hi$pos)))&&
         sourcepos$pos>=min(lo$pos,hi$pos)&&sourcepos$pos<=max(lo$pos,hi$pos))return('Compatible range')
    }
    'Different current annotation'
  }
  all[,MapMatch:=mapply(map_status,Cytoband,CurrentCytoband)]
  all[,SymbolChanged:=!is.na(CurrentSymbol)&CurrentSymbol!=SourceGene]
  examples<-fread(file.path(out,'display_examples.csv'))
  examples<-merge(examples,all,by=c('Cytoband','SourceGene'),all.x=TRUE)
  fwrite(examples,file.path(out,'figure_example_annotation_audit.csv'),na='NA')
  stopifnot(nrow(examples)==53L,!anyNA(examples$CurrentGeneID),all(examples$MapMatch%in%c('Exact match','Compatible range')))
  examples[,ExampleLabel:=CurrentSymbol];setorder(examples,Cytoband,ExampleOrder)
  all[,DisplayExample:=SourceGene%in%examples$SourceGene]
  summary<-all[,.(MeasuredGeneFeatures=.N,CurrentGeneIDs=uniqueN(CurrentGeneID[!is.na(CurrentGeneID)]),
    ProteinCoding=sum(GeneType=='PROTEIN_CODING',na.rm=TRUE),
    SourceSymbols=paste(sort(SourceGene),collapse=';'),
    CurrentSymbols=paste(sort(unique(CurrentSymbol[!is.na(CurrentSymbol)])),collapse=';')),
    by=Cytoband]
  labels<-examples[,.(GeneExamples=paste(ExampleLabel,collapse=', ')),by=Cytoband]
  summary<-merge(summary,labels,by='Cytoband');summary<-merge(summary,bands,by='Cytoband')
  setorder(summary,Direction,X)
  paths<-c('4-Multiomic/scripts/prepare_fig4_cnv_gene_annotations.R',
    '4-Multiomic/scripts/fetch_fig4_cnv_gene_reference.py',
    '4-Multiomic/refinements/fig4_cnv_gene_annotations_20261030/display_examples.csv',
    '4-Multiomic/refinements/fig4_cnv_pairwise_20261030/pairwise_display.csv',
    '4-Multiomic/refinements/fig4_cnv_20261030/gene_annotation_audit.csv',
    '0-Data/raw/2022Dong/WES_CNV_GISTIC2_raw.txt',
    paste0('4-Multiomic/refinements/fig4_cnv_gene_annotations_20261030/reference/',c(refs$File,'reference_manifest.csv','retrieval.json')))
  fwrite(all,file.path(out,'all_measured_region_genes.csv'),na='NA')
  fwrite(examples,file.path(out,'figure_gene_examples.csv'),na='NA')
  fwrite(summary,file.path(out,'region_gene_summary.csv'),na='NA')
  fwrite(data.table(Path=paths,MD5=unname(tools::md5sum(file.path(root,paths)))),file.path(out,'input_hashes.csv'))
  print(all[,.N,by=.(AnnotationStatus,MapMatch)]);print(summary[,.(Cytoband,MeasuredGeneFeatures,GeneExamples)])
  invisible(out)
}
