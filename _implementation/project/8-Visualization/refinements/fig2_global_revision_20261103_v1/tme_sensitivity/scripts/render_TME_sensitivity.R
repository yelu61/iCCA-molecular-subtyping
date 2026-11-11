#!/usr/bin/env Rscript
# Display stored contrasts only. No refitting, rescaling, or multiple-testing recalculation.
options(stringsAsFactors = FALSE, warn = 2)
library(grid)
root <- normalizePath(getwd())
base <- file.path(root, "8-Visualization/refinements/fig2_global_revision_20261103_v1/tme_sensitivity")
inputs <- c(selection = "8-Visualization/refinements/sfig2_compact_20261031_v35/tables/TME_primary_LOSO_sensitivity.csv",
            primary = "2-TranscriptomicCharacterization/runs/20261027_tme_display_v3/tables/primary_contrasts.csv",
            sensitivity = "2-TranscriptomicCharacterization/runs/20261027_tme_display_v3/tables/sensitivity_contrasts.csv")
src <- lapply(file.path(root, inputs), read.csv, check.names = FALSE)
names(src) <- names(inputs)
sel <- src$selection
key <- function(d, branch = FALSE) do.call(paste, c(d[c("Method", "Score", "Subtype", if(branch) "Sensitivity")], sep = "||"))
stopifnot(nrow(sel) == 99L, !anyDuplicated(key(sel)), !anyDuplicated(key(src$primary)),
          !anyDuplicated(key(src$sensitivity, TRUE)))
pidx <- match(key(sel), key(src$primary))
sidx <- match(key(sel, TRUE), key(src$sensitivity, TRUE))
stopifnot(!anyNA(pidx), !anyNA(sidx))
primary <- src$primary[pidx, ]; sensitivity <- src$sensitivity[sidx, ]
stored_fields <- c("Effect", "SE", "Lower", "Upper", "P", "FDR", "N", "Studies", "Unit", "NTarget", "NOthers", "Estimable", "Reason", "FamilyTests", "FamilySpecified")
for (f in stored_fields) stopifnot(isTRUE(all.equal(primary[[f]], sel[[f]], check.attributes = FALSE)))
stopifnot(isTRUE(all.equal(sensitivity$Effect, sel$SensitivityEffect, check.attributes = FALSE)),
          identical(sensitivity$N, sel$SensitivityN), identical(sensitivity$Estimable, sel$SensitivityEstimable),
          all(primary$Estimable), sum(sensitivity$Estimable) == 63L,
          all(primary$Unit == "within-study SD"), all(sensitivity$Unit == "within-study SD"))
for (role in c("primary", "sensitivity")) {
  d <- get(role)
  d$Analysis <- if(role == "primary") "Primary" else "Sensitivity"
  d$DisplayRow <- sel$RowID
  d$DisplayLabel <- sel$Label
  d$SelectedSensitivity <- sel$Sensitivity
  d$DisplayColumn <- 2L * (match(d$Subtype, c("C1", "C2", "C3")) - 1L) + if(role == "primary") 1L else 2L
  d$FDRSignificant <- !is.na(d$FDR) & d$FDR < 0.05
  d$DisplayText <- ifelse(d$Estimable & is.finite(d$Effect), sprintf("%+.2f", d$Effect),
                        ifelse(d$Reason == "Unit leverage prevents HC3 covariance", "NE\u00b9", "NE\u00b2"))
  tiny <- d$Estimable & is.finite(d$Effect) & abs(d$Effect) < 0.005
  d$DisplayText[tiny] <- sprintf("%+.3f", d$Effect[tiny])
  # Preserve the sign of small negative stored estimates, including C2 dendritic total.
  d$DisplayClipped <- is.finite(d$Effect) & abs(d$Effect) > 2
  assign(role, d)
}
sensitivity$Sensitivity <- NULL
display <- rbind(primary, sensitivity)
display <- display[order(display$DisplayRow, display$DisplayColumn), ]
stopifnot(nrow(display) == 198L, !any(display$DisplayClipped))
rows <- sel[!duplicated(sel$RowID), c("RowID", "Method", "Score", "Label", "Sensitivity")]
rows <- rows[order(rows$RowID), ]
stopifnot(nrow(rows) == 33L, identical(rows$RowID, seq_len(33L)))
for (i in seq_len(nrow(rows))) {
  d <- display[display$DisplayRow == rows$RowID[i], ]
  p <- d[d$Analysis == "Primary", ]; s <- d[d$Analysis == "Sensitivity", ]
  stopifnot(nrow(p) == 3L, nrow(s) == 3L, length(unique(p$N)) == 1L,
            length(unique(s$N)) == 1L, length(unique(p$Studies)) == 1L, length(unique(s$Studies)) == 1L)
  rows$PrimaryN[i] <- p$N[1]; rows$PrimaryStudies[i] <- p$Studies[1]
  rows$SensitivityN[i] <- s$N[1]; rows$SensitivityStudies[i] <- s$Studies[1]
}
write.csv(display, file.path(base, "source_data/TME_primary_sensitivity_198_cells.csv"), row.names = FALSE, na = "")
write.csv(rows, file.path(base, "source_data/TME_readout_denominators_33_rows.csv"), row.names = FALSE, na = "")
for (nm in names(inputs)) file.copy(file.path(root, inputs[[nm]]), file.path(base, "source_data", paste0("input_", nm, ".csv")), overwrite = TRUE)
write.csv(data.frame(Source = unname(inputs), MD5 = unname(tools::md5sum(file.path(root, inputs)))),
          file.path(base, "review/input_hashes.csv"), row.names = FALSE)

family <- "Arial"; ink <- "#252B32"; soft <- "#59626B"; rule <- "#D5DADF"
mm <- function(x) unit(x, "mm")
txt <- function(s, x, y, size = 6.8, face = "plain", col = ink, just = "left")
  grid.text(s, x = mm(x), y = mm(y), just = just,
            gp = gpar(fontfamily = family, fontsize = size, fontface = face, col = col, lineheight = 1.0))
seg <- function(x1, y1, x2, y2, col = rule, lwd = 0.45)
  grid.segments(mm(x1), mm(y1), mm(x2), mm(y2), gp = gpar(col = col, lwd = lwd))
rect <- function(x, y, w, h, fill, col = NA)
  grid.rect(mm(x), mm(y), mm(w), mm(h), gp = gpar(fill = fill, col = col, lwd = 0.35))
ramp <- grDevices::colorRamp(c("#79AFCC", "#FFFFFF", "#D98A79"), space = "Lab")
effect_color <- function(v) if(!is.finite(v)) "#EDF0F2" else grDevices::rgb(ramp((max(-2, min(2, v)) + 2) / 4), maxColorValue = 255)
method_notes <- c(ESTIMATE = "Sensitivity: seven-study native coverage \u226580%",
                  `MCP-counter` = "Sensitivity: seven-study native coverage \u226580%",
                  EPIC = "Sensitivity: fit Pearson r \u22650.3; HC3 covariance not estimable\u00b9",
                  xCell = "Sensitivity: common RNA gene universe; coverage \u226580%",
                  CIBERSORT = "Sensitivity: deconvolution fit P < 0.05; estimable studies")
draw <- function() {
  grid.newpage()
  pushViewport(viewport(x = 0, y = 0, width = mm(174), height = mm(128), just = c("left", "bottom")))
  rect(87, 64, 174, 128, "white")
  txt("Method-specific TME sensitivity", 1.5, 124.6, size = 8, face = "bold")
  txt("Study-adjusted subtype-versus-others effects; 33 readouts \u00d7 3 contrasts \u00d7 2 analyses", 1.5, 120.7, size = 6.8, col = soft)
  centers <- c(95.5,108.5,124,137,152.5,165.5)
  txt("Readout",1.5,112.8,size=6.8,face="bold")
  txt("n / studies",70,117.0,size=6.8,face="bold",just="center")
  txt("Primary",62,112.8,size=6.6,face="bold",just="center")
  txt("Sensitivity",78,112.8,size=6.6,face="bold",just="center")
  for (k in 1:3) {
    js <- c(2*k-1,2*k)
    txt(paste0("C",k),mean(centers[js]),117,size=7.5,face="bold",just="center")
    for (j in js) txt(c("Primary","Sensitivity")[(j+1)%%2+1],centers[j],112.8,size=6.6,face="bold",just="center")
  }
  seg(1.5,110.4,172,110.4)
  y <- 108.9
  for (method in unique(rows$Method)) {
    rect(86.75,y,170.5,2.8,"#F2F4F6")
    txt(method,1.8,y,size=6.8,face="bold")
    txt(method_notes[[method]],29,y,size=6.3,col=soft)
    y <- y-3.0
    ids <- rows$RowID[rows$Method==method]
    for (id in ids) {
      rr <- rows[rows$RowID==id, ]
      txt(rr$Label,2.5,y,size=6.8)
      txt(paste0(rr$PrimaryN," / ",rr$PrimaryStudies),62,y,size=6.5,just="center")
      txt(paste0(rr$SensitivityN," / ",rr$SensitivityStudies),78,y,size=6.5,just="center")
      d <- display[display$DisplayRow==id, ]
      for (j in 1:6) {
        dd <- d[d$DisplayColumn==j, ]
        rect(centers[j],y,12.65,2.5,effect_color(dd$Effect))
        txt(dd$DisplayText,centers[j]-0.35,y,size=6.6,just="center",col=if(dd$Estimable)ink else soft)
        if (dd$FDRSignificant) grid.circle(mm(centers[j]+4.65),mm(y+0.55),r=mm(.23),gp=gpar(fill=ink,col=NA))
      }
      y <- y-2.6
    }
  }
  # Explicit shared effect scale; no estimates exceed the displayed limits.
  bar_x <- seq(106,138,length.out=128)
  for (i in seq_along(bar_x)) rect(bar_x[i],5.4,32/127+0.02,1.6,effect_color(seq(-2,2,length.out=128)[i]))
  txt("Effect (within-study SD)",101,5.4,size=6.4,just="right")
  txt("\u22122",105,5.4,size=6.1,just="right"); txt("+2",139.2,5.4,size=6.1)
  txt("0",122,2.9,size=6.1,just="center")
  grid.circle(mm(2.0),mm(5.4),r=mm(.23),gp=gpar(fill=ink,col=NA))
  txt("Original BH FDR < 0.05",3.5,5.4,size=6.5)
  txt("NE\u00b9: HC3 covariance unavailable.  NE\u00b2: no eligible observations.  FDR families remain method-specific.",1.5,1.8,size=6.2,col=soft)
  popViewport()
}
stem <- file.path(base,"figures/SFigure2G_TME_primary_sensitivity")
grDevices::cairo_pdf(paste0(stem,".pdf"),width=174/25.4,height=128/25.4,family=family);draw();dev.off()
svglite::svglite(paste0(stem,".svg"),width=174/25.4,height=128/25.4,bg="white");draw();dev.off()
ragg::agg_png(paste0(stem,".png"),width=174,height=128,units="mm",res=600,background="white");draw();dev.off()
ragg::agg_png(file.path(base,"review/SFigure2G_TME_primary_sensitivity_preview.png"),width=174,height=128,units="mm",res=300,background="white");draw();dev.off()
cat(sprintf("Exported 198 cells: primary %d/99 estimable; sensitivity %d/99 estimable. Range %.9f to %.9f; no clipping.\n",sum(primary$Estimable),sum(sensitivity$Estimable),min(display$Effect,na.rm=TRUE),max(display$Effect,na.rm=TRUE)))
