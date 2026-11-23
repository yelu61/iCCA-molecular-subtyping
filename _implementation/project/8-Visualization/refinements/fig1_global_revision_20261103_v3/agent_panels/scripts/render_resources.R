#!/usr/bin/env Rscript
# R-only resource table; display only, no statistical fitting or source mutation.
options(stringsAsFactors = FALSE, warn = 2)
library(grid)
root <- normalizePath(getwd())
base <- file.path(root, "8-Visualization/refinements/fig1_global_revision_20261103_v3/agent_panels")
source <- file.path(root, "8-Visualization/refinements/fig1_sfig1_typography_20261030/source_data")
assay <- read.csv(file.path(source, "SFig1A.csv"), check.names = FALSE)
expr <- read.csv(file.path(source, "SFig1B.csv"), check.names = FALSE)
study_colors <- read.csv(file.path(root, "0-Data/qc/study_color_dictionary.csv"))
stopifnot(nrow(assay) == 48L, nrow(expr) == 13L, length(unique(expr$Study)) == 12L,
          length(unique(expr$Dataset)) == 13L, all(expr$Study %in% study_colors$Study),
          all(expr$Genes == expr$FinalProcessedFeatures), sum(expr$Samples) == 1092L,
          sum(expr$Samples[expr$AnalysisRole == "Discovery RNA-seq"]) == 586L,
          sum(expr$Samples[expr$AnalysisRole == "Microarray projection"]) == 506L)
for (st in unique(expr$Study)) {
  a <- assay[assay$Study == st & assay$Assay == "Expression", ]
  stopifnot(nrow(a) == 1L, sum(expr$Samples[expr$Study == st]) == a$Samples,
            a$Samples == a$MatchedCuratedExpressionSamples)
}
display <- expr
assay_names <- c(DNA = "DNA-seq", MS = "Protein (mass spectrometry)", RPPA = "Protein (RPPA)")
for (code in names(assay_names)) {
  a <- assay[assay$Assay == assay_names[[code]], ]
  ix <- match(expr$Study, a$Study)
  stopifnot(!anyNA(ix))
  display[[paste0(code, "_matched")]] <- a$MatchedCuratedExpressionSamples[ix]
  display[[paste0(code, "_recorded")]] <- a$Samples[ix]
  display[[paste0(code, "_denominator_definition")]] <- a$DenominatorDefinition[ix]
  display[[paste0(code, "_display")]] <- a$DisplayLabel[ix]
  is_available <- a$Available[ix]
  stopifnot(all(display[[paste0(code, "_display")]][is_available] ==
    paste0(display[[paste0(code, "_matched")]][is_available], "/", display[[paste0(code, "_recorded")]][is_available])),
    all(display[[paste0(code, "_display")]][!is_available] == "\u2014"))
}
write.csv(display, file.path(base, "source_data/resources_display.csv"), row.names = FALSE, na = "NA")
for (f in c("SFig1A.csv", "SFig1B.csv")) file.copy(file.path(source, f), file.path(base, "source_data", f), overwrite = TRUE)
manifest <- data.frame(Source = file.path(source, c("SFig1A.csv", "SFig1B.csv")),
  MD5 = unname(tools::md5sum(file.path(source, c("SFig1A.csv", "SFig1B.csv")))))
write.csv(manifest, file.path(base, "review/resources_input_hashes.csv"), row.names = FALSE)

family <- "Arial"; ink <- "#30343A"; soft <- "#686F77"; rule <- "#D8DCE0"
mm <- function(x) unit(x, "mm")
text_at <- function(s, x, y, fontsize = 6.8, face = "plain", col = ink, just = "left")
  grid.text(s, x = mm(x), y = mm(y), just = just,
    gp = gpar(fontfamily = family, fontsize = fontsize, fontface = face, col = col, lineheight = 1.06))
line_at <- function(y) grid.segments(mm(1.3),mm(y),mm(172.7),mm(y),gp=gpar(col=rule,lwd=.55))
draw <- function() {
  grid.newpage()
  pushViewport(viewport(x=0,y=0,width=mm(174),height=mm(56),just=c("left","bottom")))
  text_at("Expression resources and matched molecular assays", 1.3, 53, fontsize=7.5, face="bold")
  x <- c(1.3, 28, 74, 99, 120, 143, 166)
  text_at("Study / resource", x[1], 47.4, fontsize=6.7, face="bold")
  text_at("Role", x[2], 47.4, fontsize=6.7, face="bold")
  headers <- c("Processed\nfeatures", "Expression\nsamples")
  for (j in 3:4) text_at(headers[j-2],x[j],47.4,fontsize=6.5,face="bold",just="center")
  text_at("Matched / recorded",143,50.4,fontsize=6.5,face="bold",just="center")
  for (j in 5:7) text_at(c("DNA","Protein-MS","RPPA")[j-4],x[j],46.5,fontsize=6.7,face="bold",just="center")
  line_at(43.4)
  yy <- seq(40.6, 8.2, length.out=13)
  for (i in seq_len(nrow(display))) {
    d <- display[i, ]
    text_at(d$Dataset, x[1], yy[i])
    text_at(d$AnalysisRole, x[2], yy[i])
    text_at(format(d$FinalProcessedFeatures,big.mark=",",scientific=FALSE,trim=TRUE),x[3],yy[i],just="center")
    text_at(d$Samples,x[4],yy[i],just="center")
    text_at(d$DNA_display,x[5],yy[i],just="center")
    text_at(d$MS_display,x[6],yy[i],just="center")
    text_at(d$RPPA_display,x[7],yy[i],just="center")
  }
  line_at(mean(yy[7:8]))
  text_at("Jusakual-1 and -2 are two matrices from one study (34 + 49 = 83 samples).",1.3,5.0,fontsize=6.5,col=soft)
  text_at("Matched/recorded = expression-linked samples / samples in that source assay. A dash denotes unavailable assay data.",1.3,1.9,fontsize=6.5,col=soft)
  popViewport()
}
stem <- file.path(base,"figures/SFig1_resources")
grDevices::cairo_pdf(paste0(stem,".pdf"),width=174/25.4,height=56/25.4,family=family);draw();dev.off()
svglite::svglite(paste0(stem,".svg"),width=174/25.4,height=56/25.4,bg="white");draw();dev.off()
ragg::agg_png(paste0(stem,".png"),width=174,height=56,units="mm",res=600,background="white");draw();dev.off()
ragg::agg_png(file.path(base,"review/SFig1_resources_preview.png"),width=174,height=56,units="mm",res=300,background="white");draw();dev.off()
cat("Resources table: 13 matrices, 12 studies, 1092 expression samples; all source assay cells retained.\n")
