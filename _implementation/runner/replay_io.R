# Read-through copies make the historical relative paths portable without writing upstream.
origin<-Sys.getenv('ICCA_ORIGINAL_ROOT');work<-Sys.getenv('ICCA_REPLAY_ROOT');bundle<-Sys.getenv('ICCA_BUNDLE_ROOT')
base_exists<-base::file.exists;base_normalize<-base::normalizePath;base_source<-base::source
materialize<-function(p){
 if(!is.character(p)||length(p)!=1||!nzchar(p)||grepl('^(https?|ftp)://',p))return(p)
 cat(Sys.getenv('ICCA_JOB_ID'),p,sep='\t',file=file.path(bundle,'review/production_read_access.tsv'),append=TRUE)
 cat('\n',file=file.path(bundle,'review/production_read_access.tsv'),append=TRUE)
 if(base_exists(p))return(p)
 z<-suppressWarnings(system2('python3',c(shQuote(file.path(bundle,'code/materialize_input.py')),shQuote(origin),shQuote(work),shQuote(bundle),shQuote(p)),stdout=TRUE,stderr=FALSE))
 if(length(z)&&is.null(attr(z,'status')))return(tail(z,1));p
}
file.exists<-function(...){p<-c(...);vapply(p,function(q)base_exists(materialize(q)),logical(1))}
normalizePath<-function(path,winslash='/',mustWork=NA){
 p<-vapply(path,materialize,character(1));base_normalize(p,winslash=winslash,mustWork=mustWork)
}
source<-function(file,...){base_source(materialize(file),...)}
sys.source<-function(file,envir=baseenv(),...){base::sys.source(materialize(file),envir=envir,...)}
read.csv<-function(file,...){utils::read.csv(materialize(file),...)}
read.table<-function(file,...){utils::read.table(materialize(file),...)}
read.delim<-function(file,...){utils::read.delim(materialize(file),...)}
readLines<-function(con=stdin(),...){base::readLines(if(is.character(con))materialize(con) else con,...)}
readRDS<-function(file,...){base::readRDS(materialize(file),...)}
gzfile<-function(description,...){base::gzfile(materialize(description),...)}
load<-function(file,envir=parent.frame(),...){base::load(materialize(file),envir=envir,...)}
fread<-function(input='',...){data.table::fread(if(is.character(input)&&length(input)==1&&!grepl('\n',input,fixed=TRUE))materialize(input)else input,...)}
fwrite<-function(x,file='',...){if(nzchar(file))dir.create(dirname(file),recursive=TRUE,showWarnings=FALSE);data.table::fwrite(x,file,...)}
fromJSON<-function(txt,...){jsonlite::fromJSON(if(is.character(txt)&&length(txt)==1&&!grepl('^[{[]',txt))materialize(txt)else txt,...)}
file.copy<-function(from,to,...){from<-vapply(from,materialize,character(1));dir.create(if(length(from)>1||dir.exists(to))to else dirname(to),recursive=TRUE,showWarnings=FALSE);base::file.copy(from,to,...)}
write.csv<-function(x,file='',...){if(is.character(file)&&nzchar(file))dir.create(dirname(file),recursive=TRUE,showWarnings=FALSE);utils::write.csv(x,file,...)}
cairo_pdf<-function(filename,...){dir.create(dirname(filename),recursive=TRUE,showWarnings=FALSE);grDevices::cairo_pdf(filename,...)}
pdf<-function(file,...){dir.create(dirname(file),recursive=TRUE,showWarnings=FALSE);grDevices::pdf(file,...)}
list.files<-function(path='.',...){p<-materialize(path);base::list.files(p,...)}
if(requireNamespace('readr',quietly=TRUE)){
 trace('read_csv',where=asNamespace('readr'),tracer=quote(file<-materialize(file)),print=FALSE)
 trace('read_tsv',where=asNamespace('readr'),tracer=quote(file<-materialize(file)),print=FALSE)
}
trace('md5sum',where=asNamespace('tools'),tracer=quote(files<-vapply(files,materialize,character(1))),print=FALSE)
setwd(work)
cat('Replay working directory:',getwd(),'\n');print(.libPaths())
