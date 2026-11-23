bundle<-Sys.getenv("ICCA_BUNDLE_ROOT");stopifnot(nzchar(bundle));result<-Sys.getenv("ICCA_FINALIZATION_OUTPUT_ROOT");stopifnot(nzchar(result));dir.create(file.path(result,"figures/panels"),recursive=TRUE,showWarnings=FALSE);dir.create(file.path(result,"review"),recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages({library(data.table);library(grid);library(svglite)})
w<-result;s<-file.path(bundle,'source_data/root_clinical_audit');pal<-c(C1='#0072B2',C2='#E69F00',C3='#CC79A7');cls<-names(pal);ink<-'#263E4B';muted<-'#7C8992';light<-'#E5EBEF';checking<-FALSE;active<-'';text_bounds<-list()
ah<-fread(file.path(s,'Figure6A_clinical_information_1.csv'));a<-fread(file.path(s,'Figure6A_clinical_information_2.csv'));b<-fread(file.path(s,'Figure6B_measured_duct_increment.csv'))
stopifnot(nrow(ah)==10,nrow(a)==6,nrow(b)==14)
names<-c('start','txt','ln','dot','rect','axisx','interval','drawA','drawB')
for(e in parse(file.path(s,'render_fig6_integrated_20261101_v20.R')))if(is.call(e)&&is.symbol(e[[1]])&&as.character(e[[1]])=='<-'&&is.symbol(e[[2]])&&as.character(e[[2]])%in%names)eval(e,.GlobalEnv)
for(id in c('A','B')){f<-file.path(w,'figures/panels',paste0('Figure6',id,'_source_redraw'));for(ext in c('pdf','svg')){if(ext=='pdf')cairo_pdf(paste0(f,'.pdf'),width=85/25.4,height=76/25.4,family='Arial')else svglite(paste0(f,'.svg'),width=85/25.4,height=76/25.4);if(id=='A')drawA()else drawB();dev.off()}}
writeLines(capture.output(sessionInfo()),file.path(w,'review/root_clinical_R_session.txt'))
