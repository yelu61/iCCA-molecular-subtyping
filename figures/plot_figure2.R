#!/usr/bin/env Rscript
# ==============================================================================
# Script: plot_figure2.R
# Purpose: Generate Figure 2 panels (A-E) and full composite matching publication layout.
# Data Inputs: 8-Visualization/data/figure2/
# Outputs: 8-Visualization/figures/figure2/
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(grid)
  library(ragg)
  library(png)
})

# Pinned font handling for macOS / sandboxed environment
if (!"Arial" %in% names(grDevices::pdfFonts())) grDevices::pdfFonts(Arial = grDevices::pdfFonts("Helvetica")[[1]])

# Paths & CLI Arguments
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default_val) {
  idx <- which(args == flag)
  if (length(idx) > 0 && idx < length(args)) args[idx + 1] else default_val
}

data_cands <- c(get_arg("--data-dir", ""), "figures/data/figure2", "8-Visualization/data/figure2", "data/figure2", "../data/figure2")
data_dir <- data_cands[dir.exists(data_cands) & data_cands != ""][1]
if (is.na(data_dir)) data_dir <- "8-Visualization/data/figure2"
data_dir <- normalizePath(data_dir)

out_cands <- c(get_arg("--out-dir", ""), if (dir.exists("8-Visualization/figures")) "8-Visualization/figures/figure2" else "", "figures/outputs/figure2", "outputs/figure2")
out_dir <- out_cands[out_cands != ""][1]
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Color palettes
dict_cands <- c("governance/data/subtype_label_color_dictionary.csv", "data/subtype_label_color_dictionary.csv", "../data/subtype_label_color_dictionary.csv")
dict_file <- dict_cands[file.exists(dict_cands)][1]
if (!is.na(dict_file)) {
  subtype_pal_df <- read.csv(dict_file, stringsAsFactors = FALSE)
  subtype_colors <- setNames(subtype_pal_df$Color, subtype_pal_df$FinalSubtype)
} else {
  subtype_colors <- c(C1 = "#0072B2", C2 = "#E69F00", C3 = "#CC79A7", C1_C2 = "#738859", C1_C3 = "#6675AC", C2_C3 = "#D98C53")
}
pal <- subtype_colors[c("C1", "C2", "C3")]

study_dict_cands <- c("0-Data/qc/study_color_dictionary.csv", "data/study_color_dictionary.csv", "../data/study_color_dictionary.csv")
study_dict_file <- study_dict_cands[file.exists(study_dict_cands)][1]
if (is.na(study_dict_file)) stop("Missing governed study_color_dictionary.csv; run from the repository root.")
study_pal_df <- read.csv(study_dict_file, stringsAsFactors = FALSE)
scol <- if ("ColorHex" %in% names(study_pal_df)) study_pal_df$ColorHex else study_pal_df$Color
study_colors <- setNames(scol, study_pal_df$Study)

scale_effect_low <- "#2166AC"
scale_effect_mid <- "#F7F7F7"
scale_effect_high <- "#B2182B"

mm <- function(x) unit(x, "mm")

txt <- function(s, x, y, size = 5.8, face = "plain", just = "left", col = "#25282B", lineheight = 1) {
  grid.text(s, mm(x), mm(y), just = just,
            gp = gpar(fontfamily = "Arial", fontsize = size, fontface = face, col = col, lineheight = lineheight))
}

rect_box <- function(x, y, w, h, fill, border = NA, lwd = 0.5, just = c("center", "center")) {
  grid.rect(mm(x), mm(y), mm(w), mm(h), just = just, gp = gpar(fill = fill, col = border, lwd = lwd))
}

line_seg <- function(x, y, col = "#CDD3D8", lwd = 0.5, lty = 1) {
  grid.lines(mm(x), mm(y), gp = gpar(col = col, lwd = lwd, lty = lty))
}

point_dot <- function(x, y, fill, col = fill, size = 1.2, pch = 21, lwd = 0.65) {
  grid.points(mm(x), mm(y), pch = pch, size = mm(size), gp = gpar(fill = fill, col = col, lwd = lwd))
}

base_ramp <- colorRamp(c(scale_effect_low, "#FAFAFA", scale_effect_high))
basecolor <- function(v, lim) {
  val <- pmax(0, pmin(1, (v / lim + 1) / 2))
  rgb(base_ramp(val), maxColorValue = 255)
}

clean_label <- function(s) {
  s <- gsub("^(GOBP|HALLMARK|REACTOME|KEGG_MEDICUS)_", "", s)
  s <- tolower(gsub("_", " ", s))
  s <- paste0(toupper(substr(s, 1, 1)), substring(s, 2))
  for (w in c("wnt", "nmd", "g2m", "rna", "dna", "nfkb", "tnfa")) {
    s <- gsub(paste0("\\b", w, "\\b"), toupper(w), s, ignore.case = TRUE)
  }
  s
}

# ------------------------------------------------------------------------------
# Panel A: Six Patterns Heatmap & ORA (with dedicated legend)
# ------------------------------------------------------------------------------
cat("Rendering Panel A...\n")
prep_panelA <- function() {
  model <- readRDS(file.path(data_dir, "A_pattern_model.rds"))
  cls <- fread(file.path(data_dir, "A_wgcna_classification.csv"))
  aa <- fread(file.path(data_dir, "A_annotation_evidence.csv"))
  
  counts <- model$counts
  rr <- copy(model$membership)
  columns <- copy(model$main)
  z <- model$z
  porder <- model$pattern_order
  
  pcols <- setNames(subtype_colors[counts$Context], counts$Pattern)
  modcols <- c(setNames(cls$ModuleColor, cls$Module), grey = "grey", Not_in_network = "#FFFFFF")
  rr[, ModuleColor := unname(modcols[WGCNAModule])]
  
  aa[, Score := FoldEnrichment]
  font <- "Arial"; termfont <- 5.8; lh <- 1.95; width <- 120
  
  wrap_text <- function(s, maxwidth = 45) {
    words <- strsplit(s, " +")[[1]]
    ans <- character()
    line <- ""
    for (w in words) {
      trial <- if (nzchar(line)) paste(line, w) else w
      # approximate width in mm: 5.8pt is ~2.04mm per em, ~0.55 aspect ratio
      if (nchar(trial) * 1.15 > maxwidth && nzchar(line)) {
        ans <- c(ans, line); line <- w
      } else {
        line <- trial
      }
    }
    paste(c(ans, line), collapse = "\n")
  }
  
  aa[, WrappedText := vapply(DisplayLabel, wrap_text, character(1))]
  aa[, Lines := lengths(strsplit(WrappedText, "\n", fixed = TRUE))]
  aa[, RowHeight := Lines * lh + 0.65]
  
  bottom <- 13; gap <- 0.6; left <- 5; heatwidth <- 61; height <- 123
  top <- height - 12
  plan <- copy(counts)
  plan[, PatternOrder := match(Pattern, porder)]
  setorder(plan, PatternOrder)
  
  usable <- top - bottom - gap * (nrow(plan) - 1)
  pos <- top
  for (i in seq_len(nrow(plan))) {
    h <- usable * plan$Genes[i] / 842
    plan[i, `:=`(HeatTop = pos, HeatBottom = pos - h, HeatCenter = pos - h / 2)]
    pos <- pos - h - gap
  }
  
  heights <- vapply(plan$Pattern, function(p) 2.9 + sum(aa[Pattern == p, RowHeight]), numeric(1))
  offsets <- c(0, cumsum(head(heights, -1) + 0.75))
  desired <- plan$HeatCenter + heights / 2 + offsets
  lower <- bottom + tail(heights, 1) + tail(offsets, 1)
  fitted <- pmax(lower, pmin(top, -isoreg(seq_along(desired), -desired)$yf))
  
  for (i in seq_len(nrow(plan))) {
    p <- plan$Pattern[i]
    ix <- which(aa$Pattern == p)
    ht <- heights[i]
    at <- fitted[i] - offsets[i]
    plan[i, `:=`(HeaderY = at - 0.8, AnnotationTop = at, AnnotationBottom = at - ht)]
    cursor <- at - 2.9
    for (j in ix) {
      aa[j, `:=`(TextY = cursor - Lines * lh / 2 + 0.1, BarY = cursor - Lines * lh - 0.6,
                 YTop = cursor, YBottom = cursor - RowHeight)]
      cursor <- cursor - aa$RowHeight[j]
    }
  }
  
  states <- c("C1", "C2", "C3")
  sw <- (heatwidth - 2 * gap) / nrow(columns)
  xx <- numeric(nrow(columns))
  slices <- list()
  pos <- left
  for (st in states) {
    ix <- which(columns$FinalSubtype == st)
    xx[ix] <- pos + (seq_along(ix) - 1) * sw
    slices[[st]] <- data.table(State = st, Start = pos, End = pos + length(ix) * sw, N = length(ix))
    pos <- pos + length(ix) * sw + gap
  }
  slices <- rbindlist(slices)
  columns[, `:=`(XStart = xx, XEnd = xx + sw)]
  
  for (i in seq_len(nrow(plan))) {
    ix <- which(rr$Pattern == plan$Pattern[i])
    h <- (plan$HeatTop[i] - plan$HeatBottom[i]) / length(ix)
    rr[ix, `:=`(YTop = plan$HeatTop[i] - (seq_along(ix) - 1) * h,
                YBottom = plan$HeatTop[i] - seq_along(ix) * h)]
  }
  
  zcols <- colorRampPalette(c(scale_effect_low, scale_effect_mid, scale_effect_high))(256)
  zcolor <- function(a) zcols[pmax(1L, pmin(256L, 1L + floor((pmax(-2, pmin(2, a)) + 2) / 4 * 255)))]
  maxscore <- 30
  
  draw_A <- function(newpage = TRUE) {
    if (newpage) grid.newpage()
    for (k in seq_len(nrow(slices))) {
      sl <- slices[k]
      rect_box(sl$Start, height - 5.3, sl$End - sl$Start, 1.5, subtype_colors[[sl$State]], just = c("left", "bottom"))
    }
    for (j in seq_len(nrow(columns))) {
      rect_box(xx[j], height - 8.3, sw, 1.5, study_colors[[columns$Study[j]]], just = c("left", "bottom"))
    }
    txt("Subtype", 68.3, height - 4.55, 5.5)
    txt("Study", 68.3, height - 7.55, 5.5)
    
    for (i in seq_len(nrow(plan))) {
      p <- plan$Pattern[i]
      gs <- rr[Pattern == p, Gene]
      yr <- 4L
      for (k in seq_len(nrow(slices))) {
        sl <- slices[k]
        ix <- which(columns$FinalSubtype == sl$State)
        a <- z[gs, columns$Sample[ix], drop = FALSE]
        colors_mat <- matrix(zcolor(a), nrow = nrow(a), ncol = ncol(a))
        expanded <- colors_mat[rep(seq_len(nrow(a)), each = yr), rep(seq_len(ncol(a)), each = 5), drop = FALSE]
        grid.raster(as.raster(expanded), x = mm(sl$Start), y = mm(plan$HeatBottom[i]),
                    width = mm(sl$End - sl$Start), height = mm(plan$HeatTop[i] - plan$HeatBottom[i]),
                    just = c("left", "bottom"), interpolate = FALSE)
      }
      a_rr <- rr[Pattern == p]
      for (j in seq_len(nrow(a_rr))) {
        rect_box(2.3, a_rr$YBottom[j], 1.1, a_rr$YTop[j] - a_rr$YBottom[j],
                 a_rr$ModuleColor[j], border = if (a_rr$WGCNAModule[j] == "Not_in_network") "grey60" else NA,
                 just = c("left", "bottom"))
      }
      line_seg(c(66.7, 68.3, 70.5), c(plan$HeatCenter[i], plan$HeatCenter[i], plan$HeaderY[i] - 0.3))
      txt(p, 72, plan$HeaderY[i], 6.1, face = "bold")
      
      a_aa <- aa[Pattern == p]
      for (j in seq_len(nrow(a_aa))) {
        bw <- 45 * a_aa$Score[j] / maxscore
        rect_box(72, a_aa$BarY[j], bw, 0.45, pcols[[p]], just = c("left", "bottom"))
        txt(a_aa$WrappedText[j], 72, a_aa$TextY[j], termfont, lineheight = lh / (termfont * 25.4 / 72))
      }
    }
    line_seg(c(72, 117), c(9.5, 9.5), col = "grey55")
    for (v in c(0, maxscore / 2, maxscore)) {
      x_v <- 72 + 45 * v / maxscore
      line_seg(c(x_v, x_v), c(9.5, 8.9), col = "grey55")
      txt(as.character(v), x_v, 6.9, 5.1, just = "center")
    }
    txt("Fold enrichment", 94.5, 3.0, 5.8, just = "center")
  }
  
  draw_A_legend <- function(newpage = TRUE) {
    if (newpage) grid.newpage()
    key <- function(label, color, x, y, border = NA) {
      rect_box(x + 0.55, y, 1.1, 1.1, color, border = border, just = c("center", "center"))
      txt(label, x + 1.7, y, 5.0)
    }
    txt("Subtype / pattern", 2, 36, 5.0, face = "bold")
    labs <- c("C1", "C2", "C3", "P2", "P4")
    cols_sp <- c(subtype_colors[c("C1", "C2", "C3")], pcols[c("P2", "P4")])
    for (i in 1:5) key(labs[i], cols_sp[i], 18 + (i - 1) * 8, 36)
    
    txt("Study", 2, 32.5, 5.2, face = "bold")
    stud_names <- names(study_colors)[names(study_colors) %in% columns$Study]
    for (i in seq_along(stud_names)) {
      key(stud_names[i], study_colors[stud_names[i]], 2 + ((i - 1) %% 3) * 20, 29.8 - ((i - 1) %/% 3) * 2.5)
    }
    
    txt("WGCNA", 2, 21.8, 5.2, face = "bold")
    mod_names <- names(modcols)[names(modcols) %in% rr$WGCNAModule]
    for (i in seq_along(mod_names)) {
      m <- mod_names[i]
      key(if (m == "Not_in_network") "Outside network" else sub("^ME", "", m),
          modcols[m], 2 + ((i - 1) %% 4) * 14, 19 - ((i - 1) %/% 4) * 2.5,
          border = if (m == "Not_in_network") "grey60" else NA)
    }
    txt("Z-score", 2, 3.2, 5.0)
    for (i in 1:256) {
      rect_box(17 + (i - 0.5) * 25 / 256, 3.2, 25 / 256 + 0.01, 1.2, zcols[i], just = c("center", "center"))
    }
    for (v in c(-2, 0, 2)) {
      txt(v, 17 + (v + 2) / 4 * 25, 0.9, 5.0, just = "center")
    }
  }
  
  list(draw_A = draw_A, draw_A_legend = draw_A_legend)
}

panelA_objs <- prep_panelA()

out_panelA_png <- file.path(out_dir, "Figure_2_PanelA.png")
out_panelA_pdf <- file.path(out_dir, "Figure_2_PanelA.pdf")
agg_png(out_panelA_png, width = 120, height = 123, units = "mm", res = 600, background = "white")
panelA_objs$draw_A(TRUE)
dev.off()
grDevices::cairo_pdf(out_panelA_pdf, width = 120 / 25.4, height = 123 / 25.4, family = "Arial")
panelA_objs$draw_A(TRUE)
dev.off()

out_panelA_leg_png <- file.path(out_dir, "Figure_2_PanelA_legend.png")
out_panelA_leg_pdf <- file.path(out_dir, "Figure_2_PanelA_legend.pdf")
agg_png(out_panelA_leg_png, width = 60, height = 39, units = "mm", res = 600, background = "white")
panelA_objs$draw_A_legend(TRUE)
dev.off()
grDevices::cairo_pdf(out_panelA_leg_pdf, width = 60 / 25.4, height = 39 / 25.4, family = "Arial")
panelA_objs$draw_A_legend(TRUE)
dev.off()

# ------------------------------------------------------------------------------
# Panel B: GSEA NES Panel (18 Pathways)
# ------------------------------------------------------------------------------
cat("Rendering Panel B...\n")
draw_panelB <- function(newpage = TRUE, w = 60, h = 82) {
  z <- fread(file.path(data_dir, "B_gsea_selected54.csv"))
  color <- function(v) basecolor(v, 4.5)
  plan <- unique(z[, .(pathway, AnchorSubtype, AnchorDirection, SelectionOrder, Collection)])
  plan[, DirOrder := match(AnchorDirection, c("Positive", "Negative"))]
  setorder(plan, AnchorSubtype, DirOrder, SelectionOrder)
  plan[, Label := clean_label(pathway)]
  
  step <- 3.2
  if (newpage) grid.newpage()
  xs <- c(w - 14, w - 8.5, w - 3)
  for (j in 1:3) {
    txt(names(pal)[j], xs[j], h - 4, 6.3, face = "bold", just = "center")
    rect_box(xs[j], h - 7, 5.1, 1.1, pal[j], just = c("center", "center"))
  }
  txt("DB", w - 19, h - 4, 5.0, just = "center", col = "grey45")
  txt("+/-", w - 23, h - 4, 5.0, just = "center", col = "grey45")
  cursor <- h - 12
  
  for (s in names(pal)) {
    txt(paste0("Selected via ", s, " contrast"), 1, cursor, 5.4, face = "bold", col = "grey35")
    cursor <- cursor - 3
    for (id in plan[AnchorSubtype == s, pathway]) {
      r <- plan[pathway == id]
      txt(r$Label, 1, cursor, 5.2)
      txt(if (r$AnchorDirection == "Positive") "+" else "-", w - 23, cursor, 5.2, just = "center", col = pal[s])
      txt(c(GO_BP = "G", Reactome = "R", Hallmark = "H", KEGG = "K")[[r$Collection]],
          w - 19, cursor, 5.1, just = "center", col = "grey45")
      for (j in 1:3) {
        a <- z[pathway == id & Subtype == names(pal)[j]]
        rect_box(xs[j], cursor, 5.1, step - 0.3, color(a$NES), just = c("center", "center"))
        if (a$padj < 0.05) {
          point_dot(xs[j], cursor, fill = if (abs(a$NES) > 3) "white" else "grey20",
                    col = if (abs(a$NES) > 3) "white" else "grey20", size = 0.55, pch = 16)
        }
      }
      cursor <- cursor - step
    }
    cursor <- cursor - 0.8
  }
  
  txt("NES", 1, 8, 5.5)
  for (i in 1:100) {
    rect_box(8 + (i - 0.5) * 0.15, 8, 0.16, 1.4, color(-4.5 + (i - 1) / 99 * 9), just = c("center", "center"))
  }
  for (i in 1:3) {
    txt(c("-4.5", "0", "4.5")[i], c(8, 15.5, 23)[i], 5.2, 5.0, just = "center")
  }
  point_dot(30, 8, fill = "grey20", col = "grey20", size = 0.55, pch = 16)
  txt("FDR < 0.05", 32, 8, 5.2)
  txt("G GO BP   R Reactome   H Hallmark", 1, 1.8, 5.0, col = "grey40")
}

out_panelB_png <- file.path(out_dir, "Figure_2_PanelB.png")
out_panelB_pdf <- file.path(out_dir, "Figure_2_PanelB.pdf")
agg_png(out_panelB_png, width = 60, height = 82, units = "mm", res = 600, background = "white")
draw_panelB(TRUE, 60, 82)
dev.off()
grDevices::cairo_pdf(out_panelB_pdf, width = 60 / 25.4, height = 82 / 25.4, family = "Arial")
draw_panelB(TRUE, 60, 82)
dev.off()

# ------------------------------------------------------------------------------
# Panel C: 11 WGCNA Modules Association and Preservation
# ------------------------------------------------------------------------------
cat("Rendering Panel C...\n")
draw_panelC <- function(newpage = TRUE) {
  cc <- fread(file.path(data_dir, "C_all33_effects.csv"))
  cp <- fread(file.path(data_dir, "C_all11_modules.csv"))
  
  if (newpage) grid.newpage()
  xs <- c(77, 108, 139)
  ys <- 40.5 - (0:10) * 3 - c(0, 0, rep(1.5, 4), rep(3, 5))
  scale <- 28 / 5
  
  txt("Co-expression program", 6, 49, 6.2, face = "bold")
  txt("Module | genes", 43, 49, 5.2, face = "bold")
  for (j in 1:3) {
    txt(names(pal)[j], xs[j], 49, 6.3, face = "bold", just = "center")
    rect_box(xs[j], 46.3, 27, 1.1, pal[j], just = c("center", "center"))
    line_seg(rep(xs[j], 2), c(6, 43), col = "#CCD2D8", lwd = 0.5, lty = 2)
  }
  txt("Preservation", 169, 49, 6.0, face = "bold", just = "center")
  txt("Z summary", 169, 46, 5.2, just = "center", col = "grey40")
  for (v in c(2, 10)) line_seg(rep(158 + v / 25 * 20, 2), c(6, 43), col = "#D4D9DD", lwd = 0.5, lty = 2)
  
  txt("Difference (SD), 95% CI", 6, 44, 5.1, col = "grey40")
  point_dot(86, 43.4, fill = "#455B69", col = "#455B69", size = 0.9)
  txt("FDR < 0.05", 88, 43.4, 5.0)
  point_dot(114, 43.4, fill = "white", col = "#455B69", size = 0.9)
  txt("FDR >= 0.05", 116, 43.4, 5.0)
  
  for (i in seq_len(nrow(cp))) {
    y <- ys[i]
    r <- cp[i]
    rect_box(2, y, 1.1, 2.4, r$ModuleColor, just = c("center", "center"))
    txt(r$Program, 4, y, 5.8, face = "bold")
    txt(paste0(r$ModuleColor, " | ", r$ModuleSize), 43, y, 5.0, col = "grey45")
    for (j in 1:3) {
      z <- cc[Module == r$Module & TargetSubtype == names(pal)[j]]
      line_seg(xs[j] + c(z$CI_low, z$CI_high) * scale, c(y, y), col = pal[j], lwd = 0.75)
      point_dot(xs[j] + z$RandomEffect * scale, y,
                fill = if (z$FDR < 0.05) pal[j] else "white", col = pal[j], size = 1.25)
    }
    line_seg(c(158, 158 + r$ZsummaryPreservation / 25 * 20), c(y, y), col = "#AEBAC2", lwd = 0.65)
    point_dot(158 + r$ZsummaryPreservation / 25 * 20, y, fill = "#455B69", col = "#455B69", size = 1.15)
  }
  for (j in 1:3) {
    line_seg(xs[j] + c(-2.5, 2.5) * scale, c(4.2, 4.2), col = "grey50")
    for (v in c(-2, 0, 2)) txt(v, xs[j] + v * scale, 2.0, 5.0, just = "center")
  }
  for (v in c(0, 10, 20)) txt(v, 158 + v / 25 * 20, 2.0, 5.0, just = "center")
}

out_panelC_png <- file.path(out_dir, "Figure_2_PanelC.png")
out_panelC_pdf <- file.path(out_dir, "Figure_2_PanelC.pdf")
agg_png(out_panelC_png, width = 183, height = 53, units = "mm", res = 600, background = "white")
draw_panelC(TRUE)
dev.off()
grDevices::cairo_pdf(out_panelC_pdf, width = 183 / 25.4, height = 53 / 25.4, family = "Arial")
draw_panelC(TRUE)
dev.off()

# ------------------------------------------------------------------------------
# Panel D: decoupleR Transcription Factor Activity
# ------------------------------------------------------------------------------
cat("Rendering Panel D...\n")
draw_panelD <- function(newpage = TRUE, extended = FALSE) {
  z <- fread(file.path(data_dir, "D_tf_selected36.csv"))
  color <- function(v) basecolor(v, 8)
  plan <- unique(z[, .(TF, AnchorSubtype, AnchorDirection, SelectionOrder, NonzeroSignedTargets)])
  plan[, DirOrder := match(AnchorDirection, c("Positive", "Negative"))]
  setorder(plan, AnchorSubtype, DirOrder, SelectionOrder)
  
  H <- 76; xs <- c(27, 39, 51); cw <- 11; step <- 3.35
  if (newpage) grid.newpage()
  txt("TF", 5, H - 4, 6.3, face = "bold")
  txt("+/-", 16, H - 4, 5.4, just = "center", col = "grey40")
  for (j in 1:3) {
    txt(names(pal)[j], xs[j], H - 4, 6.3, face = "bold", just = "center")
    rect_box(xs[j], H - 6.8, cw, 1.1, pal[j], just = c("center", "center"))
  }
  cursor <- H - 12
  for (s in names(pal)) {
    txt(paste0("Selected via ", s, " contrast"), 2, cursor, 5.5, face = "bold", col = "grey40")
    cursor <- cursor - 3
    for (tf in plan[AnchorSubtype == s, TF]) {
      rr <- plan[TF == tf]
      txt(tf, 2, cursor, 6.1, face = "italic")
      txt(if (rr$AnchorDirection == "Positive") "+" else "-", 16, cursor, 5.8, just = "center", col = pal[s])
      for (j in 1:3) {
        r <- z[TF == tf & Subtype == names(pal)[j]]
        v <- r$ActivityDifferenceSubtypeMinusOthers
        rect_box(xs[j], cursor, cw, step - 0.25, color(v), just = c("center", "center"))
        txt(sprintf("%.1f", v), xs[j] - 0.5, cursor, 5.6, just = "center",
            col = if (abs(v) > 5) "white" else "grey15")
        if (r$FDR < 0.05) {
          point_dot(xs[j] + cw / 2 - 1.3, cursor + 0.8,
                    fill = if (abs(v) > 5) "white" else "grey15",
                    col = if (abs(v) > 5) "white" else "grey15", size = 0.48, pch = 16)
        }
      }
      cursor <- cursor - step
    }
    cursor <- cursor - 0.8
  }
  txt("ULM difference", 2, 8, 5.3)
  for (i in 1:100) {
    rect_box(2 + (i - 0.5) * 0.18, 4.7, 0.19, 1.4, color(-8 + (i - 1) / 99 * 16), just = c("center", "center"))
  }
  for (i in 1:3) {
    txt(c("-8", "0", "8")[i], c(2, 11, 20)[i], 2.0, 5.0, just = "center")
  }
  point_dot(33, 4.7, fill = "grey20", col = "grey20", size = 0.5, pch = 16)
  txt("FDR < 0.05", 35, 4.7, 5.2)
}

out_panelD_png <- file.path(out_dir, "Figure_2_PanelD.png")
out_panelD_pdf <- file.path(out_dir, "Figure_2_PanelD.pdf")
agg_png(out_panelD_png, width = 60, height = 76, units = "mm", res = 600, background = "white")
draw_panelD(TRUE)
dev.off()
grDevices::cairo_pdf(out_panelD_pdf, width = 60 / 25.4, height = 76 / 25.4, family = "Arial")
draw_panelD(TRUE)
dev.off()

# ------------------------------------------------------------------------------
# Panel E: 33 TME Readouts
# ------------------------------------------------------------------------------
cat("Rendering Panel E...\n")
draw_panelE <- function(newpage = TRUE) {
  em <- fread(file.path(data_dir, "E_selection_map.csv"))
  ez <- fread(file.path(data_dir, "E_selected_contrasts.csv"))
  em[, Y := 63 - (Y - 14) / 0.72 * 0.39]
  ez[, Y := 63 - (Y - 14) / 0.72 * 0.39]
  
  eshape <- c("MCP-counter" = 21, "EPIC" = 22, "xCell" = 24, "CIBERSORT" = 23, "ESTIMATE" = 25)
  eoffset <- c("MCP-counter" = 0.8, "EPIC" = 0.2667, "xCell" = -0.2667, "CIBERSORT" = -0.8, "ESTIMATE" = 0)
  
  if (newpage) grid.newpage()
  centers <- c(65, 86, 107)
  dot <- function(x, y, method, fill, col, size = 1.05) {
    grid.points(mm(x), mm(y), pch = eshape[method], size = mm(size), gp = gpar(col = col, fill = fill, lwd = 0.6))
  }
  
  txt("Compartment", 5, 71.5, 5.7, face = "bold")
  txt("Readout", 25, 71.5, 5.7, face = "bold")
  for (j in 1:3) {
    txt(names(pal)[j], centers[j], 71.5, 6.3, face = "bold", just = "center")
    rect_box(centers[j], 68.8, 18, 1.1, pal[j], just = c("center", "center"))
  }
  groups <- unique(em$Group)
  for (i in seq_along(groups)) {
    g <- groups[i]
    yy <- em[Group == g, Y]
    lo <- min(yy) - 1.5; hi <- max(yy) + 1.5
    if (i %% 2 == 1) rect_box(60, (lo + hi) / 2, 118, hi - lo, "#F5F7F9", just = c("center", "center"))
    rect_box(1.8, (lo + hi) / 2, 0.6, hi - lo, "#B8C3CC", just = c("center", "center"))
    txt(if (g == "Overall context") "Overall\ncontext" else g, 3, (lo + hi) / 2, 5.3, face = "bold")
  }
  for (j in 1:3) line_seg(rep(centers[j], 2), c(24, 65), col = "#CCD2D8", lwd = 0.5, lty = 2)
  
  labs <- unique(em[, .(Label, Y)])
  for (i in seq_len(nrow(labs))) txt(labs$Label[i], 25, labs$Y[i], 5.5)
  
  for (i in seq_len(nrow(ez))) {
    r <- ez[i]
    j <- match(r$Subtype, names(pal))
    y <- r$Y + eoffset[r$Method]
    xx <- function(v) centers[j] + 4.4 * v
    line_seg(xx(c(r$Lower, r$Upper)), c(y, y), col = adjustcolor(pal[j], alpha.f = 0.5), lwd = 0.6)
    dot(xx(r$Effect), y, r$Method, fill = if (r$FDR < 0.05) pal[j] else "white", col = pal[j])
  }
  for (j in 1:3) {
    line_seg(centers[j] + c(-8.8, 8.8), c(23.7, 23.7), col = "grey50")
    for (v in c(-2, 0, 2)) txt(v, centers[j] + 4.4 * v, 21.3, 5.0, just = "center")
  }
  txt("Study-adjusted difference (within-study SD)", 86, 17.6, 5.3, just = "center")
  methods <- names(eshape)
  lx <- c(3, 30, 49, 70, 98)
  for (i in seq_along(methods)) {
    dot(lx[i], 11.5, methods[i], fill = "#444C53", col = "#444C53", size = 1.2)
    txt(methods[i], lx[i] + 2.2, 11.5, 5.3)
  }
  dot(3, 6, "MCP-counter", fill = "#444C53", col = "#444C53")
  txt("FDR < 0.05", 5.2, 6, 5.2)
  dot(32, 6, "MCP-counter", fill = "white", col = "#444C53")
  txt("FDR >= 0.05", 34.2, 6, 5.2)
  line_seg(c(65, 71), c(6, 6), col = "grey50")
  txt("95% CI", 74, 6, 5.2)
}

out_panelE_png <- file.path(out_dir, "Figure_2_PanelE.png")
out_panelE_pdf <- file.path(out_dir, "Figure_2_PanelE.pdf")
agg_png(out_panelE_png, width = 120, height = 76, units = "mm", res = 600, background = "white")
draw_panelE(TRUE)
dev.off()
grDevices::cairo_pdf(out_panelE_pdf, width = 120 / 25.4, height = 76 / 25.4, family = "Arial")
draw_panelE(TRUE)
dev.off()

# ------------------------------------------------------------------------------
# Full Composite: Figure 2 (180 mm x 270 mm)
# ------------------------------------------------------------------------------
cat("Assembling Full Figure 2 Composite...\n")
composite_png <- file.path(out_dir, "Figure_2_composite.png")
composite_pdf <- file.path(out_dir, "Figure_2_composite.pdf")

render_composite <- function() {
  grid.newpage()
  # 180 mm x 270 mm canvas
  pushViewport(viewport(x = 0, y = 0, width = mm(180), height = mm(270), just = c("left", "bottom")))
  grid.rect(gp = gpar(fill = "white", col = NA))
  
  draw_raster_panel <- function(png_path, x_mm, y_top_mm, w_mm, h_mm) {
    if (file.exists(png_path)) {
      img <- readPNG(png_path)
      y_bottom_mm <- 270 - y_top_mm - h_mm
      grid.raster(img, x = mm(x_mm + w_mm / 2), y = mm(y_bottom_mm + h_mm / 2),
                  width = mm(w_mm), height = mm(h_mm))
    }
  }
  
  draw_panel_letter <- function(letter, x_mm, y_top_mm) {
    y_from_bottom <- 270 - y_top_mm
    grid.text(letter, x = mm(x_mm), y = mm(y_from_bottom), just = c("left", "top"),
              gp = gpar(fontfamily = "Arial", fontsize = 11, fontface = "bold", col = "#111111"))
  }
  
  # Layout from spec.json:
  # A: [3, 13, 114, 116.85]
  draw_raster_panel(out_panelA_png, 3.0, 13.0, 114.0, 116.85)
  draw_panel_letter("A", 1.0, 11.5)
  
  # B: [120, 13, 57, 77.9]
  draw_raster_panel(out_panelB_png, 120.0, 13.0, 57.0, 77.9)
  draw_panel_letter("B", 117.5, 11.5)
  
  # A_legend: [120, 92.5, 57, 37.05]
  draw_raster_panel(out_panelA_leg_png, 120.0, 92.5, 57.0, 37.05)
  
  # C: [3, 136, 174, 50.4]
  draw_raster_panel(out_panelC_png, 3.0, 136.0, 174.0, 50.4)
  draw_panel_letter("C", 1.0, 134.5)
  
  # D: [3, 192, 57, 72.2]
  draw_raster_panel(out_panelD_png, 3.0, 192.0, 57.0, 72.2)
  draw_panel_letter("D", 1.0, 190.5)
  
  # E: [63, 192, 114, 72.2]
  draw_raster_panel(out_panelE_png, 63.0, 192.0, 114.0, 72.2)
  draw_panel_letter("E", 60.5, 190.5)
  
  popViewport()
}

agg_png(composite_png, width = 180, height = 270, units = "mm", res = 300, background = "white")
render_composite()
dev.off()

grDevices::cairo_pdf(composite_pdf, width = 180 / 25.4, height = 270 / 25.4, family = "Arial")
render_composite()
dev.off()

cat("Figure 2 rendered successfully:\n")
cat("  - Individual panels:", out_dir, "\n")
cat("  - Composite PDF:", composite_pdf, "\n")
cat("  - Composite PNG:", composite_png, "\n")
