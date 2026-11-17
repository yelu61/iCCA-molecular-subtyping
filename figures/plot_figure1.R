#!/usr/bin/env Rscript
# ==============================================================================
# Script: plot_figure1.R
# Purpose: Generate Figure 1 panels (A-G) and full composite matching publication layout.
# Data Inputs: 8-Visualization/data/figure1/
# Outputs: 8-Visualization/figures/figure1/
# ==============================================================================

suppressPackageStartupMessages({
  library(grid)
  library(ragg)
  library(ComplexHeatmap)
  library(circlize)
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

data_cands <- c(get_arg("--data-dir", ""), "figures/data/figure1", "8-Visualization/data/figure1", "data/figure1", "../data/figure1")
data_dir <- data_cands[dir.exists(data_cands) & data_cands != ""][1]
if (is.na(data_dir)) data_dir <- "8-Visualization/data/figure1"
data_dir <- normalizePath(data_dir)

out_cands <- c(get_arg("--out-dir", ""), if (dir.exists("8-Visualization/figures")) "8-Visualization/figures/figure1" else "", "figures/outputs/figure1", "outputs/figure1")
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

study_dict_cands <- c("0-Data/qc/study_color_dictionary.csv", "data/study_color_dictionary.csv", "../data/study_color_dictionary.csv")
study_dict_file <- study_dict_cands[file.exists(study_dict_cands)][1]
if (!is.na(study_dict_file)) {
  study_pal_df <- read.csv(study_dict_file, stringsAsFactors = FALSE)
  scol <- if ("ColorHex" %in% names(study_pal_df)) study_pal_df$ColorHex else study_pal_df$Color
  study_colors <- setNames(scol, study_pal_df$Study)
} else {
  study_colors <- c("2019Ahn" = "#88CCEE", "2022Dong" = "#332288", "2023Deng" = "#44AA99", "2023Jeon" = "#AA4499", "2023Lee" = "#DDCC77", "2024Fan" = "#117733", "TCGA-CHOL" = "#CC6677")
}

states <- c("C1", "C2", "C3", "C1_C2", "C1_C3", "C2_C3")
studies <- study_pal_df$Study[study_pal_df$AnalysisRole == "Discovery RNA-seq"]
short_studies <- setNames(c("Ahn", "Dong", "Deng", "Jeon", "Lee", "Fan", "TCGA"), studies)

method_colors <- list(
  NMF = c("N1" = "#D9D9D9", "N2" = "#969696", "N3" = "#404040"),
  KMeans = c("K1" = "#E0DBEB", "K2" = "#A99AC4", "K3" = "#66517F")
)
avail_colors <- list(
  Genome = c("Yes" = "#466E65", "No" = "#DAEAE6"),
  Protein = c("Yes" = "#856C4B", "No" = "#EAE0D2")
)

mm <- function(x) unit(x, "mm")

txt <- function(label, x, y, size = 6.5, col = "#222222", fontface = "plain", just = "centre", rot = 0, lineheight = 1.06) {
  grid.text(label, x = mm(x), y = mm(y), just = just, rot = rot,
            gp = gpar(fontfamily = "Arial", fontsize = size, col = col, fontface = fontface, lineheight = lineheight))
}

ln <- function(x0, y0, x1, y1, col = "#555555", lwd = 0.5, lty = 1) {
  grid.lines(mm(c(x0, x1)), mm(c(y0, y1)), gp = gpar(col = col, lwd = lwd, lty = lty))
}

rect_box <- function(x, y, w, h, fill, col = NA, lwd = 0.4, just = c("left", "bottom")) {
  grid.rect(mm(x), mm(y), mm(w), mm(h), just = just, gp = gpar(fill = fill, col = col, lwd = lwd))
}

# ------------------------------------------------------------------------------
# Panel A: Study Design Schematic
# ------------------------------------------------------------------------------
cat("Rendering Panel A...\n")
panelA_pdf <- file.path(data_dir, "Figure1A_schematic.pdf")
out_panelA_pdf <- file.path(out_dir, "Figure_1_PanelA.pdf")
out_panelA_png <- file.path(out_dir, "Figure_1_PanelA.png")

if (file.exists(panelA_pdf)) {
  file.copy(panelA_pdf, out_panelA_pdf, overwrite = TRUE)
  # Render PNG preview using ragg and pdftools or sips if available
  tryCatch({
    system2("sips", c("-s", "format", "png", out_panelA_pdf, "--out", out_panelA_png), stdout = FALSE, stderr = FALSE)
  }, error = function(e) NULL)
}

# ------------------------------------------------------------------------------
# Panel B: Expression Heatmap (2,000 variable genes, 586 samples)
# ------------------------------------------------------------------------------
cat("Rendering Panel B...\n")
render_panelB <- function() {
  rows <- read.csv(file.path(data_dir, "Figure1B_rows.csv"), stringsAsFactors = FALSE)
  cols <- read.csv(file.path(data_dir, "Figure1B_columns.csv"), stringsAsFactors = FALSE)
  zmat <- as.matrix(read.csv(file.path(data_dir, "Figure1B_row_zscore.csv.gz"), row.names = 1, check.names = FALSE))
  
  # Saturation limit
  lim <- 2
  zmat[zmat > lim] <- lim
  zmat[zmat < -lim] <- -lim
  
  # Ensure column order matches annotation table
  zmat <- zmat[, cols$Sample, drop = FALSE]
  
  state_keys <- states[states %in% cols$FinalSubtype]
  study_keys <- names(study_colors)[names(study_colors) %in% cols$Study]
  
  ann <- ComplexHeatmap::HeatmapAnnotation(
    State = cols$FinalSubtype,
    NMF = paste0("N", cols$NMF),
    KMeans = paste0("K", cols$KMeans),
    Study = cols$Study,
    `Genome data` = cols$GenomeData,
    `Protein data` = cols$ProteinData,
    col = list(
      State = subtype_colors,
      NMF = method_colors$NMF,
      KMeans = method_colors$KMeans,
      Study = study_colors,
      `Genome data` = avail_colors$Genome,
      `Protein data` = avail_colors$Protein
    ),
    simple_anno_size = unit(1.8, "mm"),
    annotation_name_gp = gpar(fontsize = 6, fontfamily = "Arial"),
    show_legend = FALSE
  )
  
  discrete_legend <- function(title, keys, cols_map, nrow = NULL, ncol = NULL, by_row = TRUE) {
    ComplexHeatmap::Legend(
      title = title, at = keys,
      legend_gp = gpar(fill = unname(cols_map[keys])),
      nrow = nrow, ncol = ncol, by_row = by_row,
      grid_height = unit(2.2, "mm"), grid_width = unit(2.2, "mm"),
      title_gp = gpar(fontsize = 6, fontface = "bold", fontfamily = "Arial"),
      labels_gp = gpar(fontsize = 5.5, fontfamily = "Arial")
    )
  }
  
  effect_col <- circlize::colorRamp2(c(-lim, 0, lim), c("#2166AC", "#F7F7F7", "#B2182B"))
  
  method_column <- ComplexHeatmap::packLegend(
    discrete_legend("NMF", paste0("N", 1:3), method_colors$NMF, nrow = 3),
    discrete_legend("KMeans", paste0("K", 1:3), method_colors$KMeans, nrow = 3),
    direction = "horizontal", gap = unit(2, "mm")
  )
  
  avail_column <- ComplexHeatmap::packLegend(
    discrete_legend("Genome", c("Yes", "No"), avail_colors$Genome, nrow = 2),
    discrete_legend("Protein", c("Yes", "No"), avail_colors$Protein, nrow = 2),
    direction = "horizontal", gap = unit(2, "mm")
  )
  
  legends <- ComplexHeatmap::packLegend(
    discrete_legend("State", state_keys, subtype_colors, nrow = 3, ncol = 2, by_row = FALSE),
    discrete_legend("Study", study_keys, study_colors, nrow = length(study_keys)),
    method_column, avail_column,
    ComplexHeatmap::Legend(
      title = "Row z-score", col_fun = effect_col, at = c(-lim, 0, lim), direction = "horizontal",
      legend_width = unit(16, "mm"), grid_height = unit(2.2, "mm"),
      title_gp = gpar(fontsize = 6, fontface = "bold", fontfamily = "Arial"),
      labels_gp = gpar(fontsize = 5.5, fontfamily = "Arial")
    ),
    direction = "vertical", gap = unit(2, "mm")
  )
  
  ht <- ComplexHeatmap::Heatmap(
    zmat,
    name = "Row z-score",
    col = effect_col,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_split = factor(cols$FinalSubtype, levels = states),
    column_title = NULL,
    column_gap = unit(0.35, "mm"),
    show_column_names = FALSE,
    show_row_names = FALSE,
    top_annotation = ann,
    use_raster = TRUE,
    raster_device = "agg_png",
    raster_quality = 4,
    show_heatmap_legend = FALSE
  )
  
  draw(ht, heatmap_legend_list = legends, heatmap_legend_side = "right", padding = unit(c(1, 1, 1, 1), "mm"))
}

out_panelB_png <- file.path(out_dir, "Figure_1_PanelB.png")
out_panelB_pdf <- file.path(out_dir, "Figure_1_PanelB.pdf")
agg_png(out_panelB_png, width = 114, height = 82, units = "mm", res = 600, background = "white")
render_panelB()
dev.off()

pdf(out_panelB_pdf, width = 114 / 25.4, height = 82 / 25.4, family = "Arial")
render_panelB()
dev.off()

# ------------------------------------------------------------------------------
# Panel C: UMAP (497 principal-state samples)
# ------------------------------------------------------------------------------
cat("Rendering Panel C...\n")
draw_panelC <- function() {
  u <- read.csv(file.path(data_dir, "Figure1C.csv"), stringsAsFactors = FALSE)
  grid.newpage()
  pushViewport(viewport(x = 0, y = 0, width = mm(54), height = mm(45.8), just = c("left", "bottom")))
  
  txt("Principal states (n = 497)", 27, 43.5, 7.2, fontface = "bold")
  lx <- c(13.2, 27.3, 41.4)
  for (i in 1:3) {
    rect_box(lx[i] - 5.5, 39.5, 1.6, 1.6, subtype_colors[states[i]])
    txt(paste0(states[i], " (", c(112, 223, 162)[i], ")"), lx[i] - 3.2, 40.3, 6.2, just = "left")
  }
  
  xr <- c(-4.1, 4.7); yr <- c(-3.2, 3.4)
  scale <- 4.2
  x0 <- 10.5; y0 <- 7.8
  xp <- function(x) x0 + (x - xr[1]) * scale
  yp <- function(y) y0 + (y - yr[1]) * scale
  
  grid.points(mm(xp(u$UMAP1)), mm(yp(u$UMAP2)), pch = 16, size = mm(0.65),
              gp = gpar(col = subtype_colors[u$FinalSubtype], alpha = 0.85))
  
  ln(xp(xr[1]), yp(yr[1]), xp(xr[2]), yp(yr[1]), col = "#555555", lwd = 0.45)
  ln(xp(xr[1]), yp(yr[1]), xp(xr[1]), yp(yr[2]), col = "#555555", lwd = 0.45)
  
  for (x in c(-4, -2, 0, 2, 4)) {
    ln(xp(x), y0, xp(x), y0 - 0.8, col = "#555555", lwd = 0.45)
    txt(x, xp(x), y0 - 2.1, 5.8)
  }
  for (y in c(-2, 0, 2)) {
    ln(x0, yp(y), x0 - 0.8, yp(y), col = "#555555", lwd = 0.45)
    txt(y, x0 - 1.5, yp(y), 5.8, just = "right")
  }
  txt("UMAP 1", (xp(xr[1]) + xp(xr[2])) / 2, 1.8, 6.2)
  txt("UMAP 2", 2.8, (yp(yr[1]) + yp(yr[2])) / 2, 6.2, rot = 90)
  popViewport()
}

out_panelC_png <- file.path(out_dir, "Figure_1_PanelC.png")
out_panelC_pdf <- file.path(out_dir, "Figure_1_PanelC.pdf")
agg_png(out_panelC_png, width = 54, height = 45.8, units = "mm", res = 600, background = "white")
draw_panelC()
dev.off()
pdf(out_panelC_pdf, width = 54 / 25.4, height = 45.8 / 25.4, family = "Arial")
draw_panelC()
dev.off()

# ------------------------------------------------------------------------------
# Panel D: Study Composition Barplot
# ------------------------------------------------------------------------------
cat("Rendering Panel D...\n")
draw_panelD <- function() {
  d <- read.csv(file.path(data_dir, "Figure1D.csv"), stringsAsFactors = FALSE)
  grid.newpage()
  pushViewport(viewport(x = 0, y = 0, width = mm(54), height = mm(40.1), just = c("left", "bottom")))
  
  txt("Principal-state composition", 27, 38.3, 7.2, fontface = "bold")
  positions <- c(17, 29, 41)
  for (i in 1:3) {
    rect_box(positions[i] - 3.8, 33.8, 1.8, 1.8, subtype_colors[states[i]])
    txt(states[i], positions[i] - 1.2, 34.7, 6.2, just = "left")
  }
  
  x0 <- 11.5; width <- 32.5
  for (j in seq_along(studies)) {
    z <- d[d$Study == studies[j], ]
    z <- z[match(states[1:3], z$FinalSubtype), ]
    y <- 30.0 - (j - 1) * 3.4
    txt(short_studies[studies[j]], 0.5, y, 6.2, just = "left")
    pos <- c(0, cumsum(z$Fraction))
    for (i in 1:3) {
      if (z$Fraction[i] > 0) {
        rect_box(x0 + width * pos[i], y - 1.3, width * z$Fraction[i], 2.6,
                 subtype_colors[states[i]], col = "#FFFFFF", lwd = 0.35)
        txt(z$N[i], x0 + width * (pos[i] + pos[i + 1]) / 2, y, 5.8,
            col = if (i == 1) "#FFFFFF" else "#111111")
      }
    }
    txt(paste0("n=", z$Total[1]), 53.5, y, 5.8, just = "right")
  }
  ln(x0, 6.8, x0 + width, 6.8, col = "#999999", lwd = 0.45)
  for (f in c(0, 0.5, 1)) {
    ln(x0 + width * f, 6.8, x0 + width * f, 6.1, col = "#999999", lwd = 0.45)
    txt(paste0(100 * f, "%"), x0 + width * f, 4.6, 5.5)
  }
  popViewport()
}

out_panelD_png <- file.path(out_dir, "Figure_1_PanelD.png")
out_panelD_pdf <- file.path(out_dir, "Figure_1_PanelD.pdf")
agg_png(out_panelD_png, width = 54, height = 40.1, units = "mm", res = 600, background = "white")
draw_panelD()
dev.off()
pdf(out_panelD_pdf, width = 54 / 25.4, height = 40.1 / 25.4, family = "Arial")
draw_panelD()
dev.off()

# ------------------------------------------------------------------------------
# Panel E: Published Classifications Correspondence
# ------------------------------------------------------------------------------
cat("Rendering Panel E...\n")
draw_panelE <- function() {
  e <- read.csv(file.path(data_dir, "Figure1E_all_cells.csv"), stringsAsFactors = FALSE)
  grid.newpage()
  pushViewport(viewport(x = 0, y = 0, width = mm(174), height = mm(44), just = c("left", "bottom")))
  
  ids <- c("Figure1E1", "Figure1E3", "Figure1E2", "Figure1E4")
  heads <- c("Ahn", "Dong: RNA", "Dong: protein", "Deng: protein")
  subs <- c("Original A/B classes", "IG1–IG3", "S1–S4", "S-I–S-III")
  shade <- colorRamp(c("#FAFCFD", "#0072A8"))
  
  for (i in seq_along(ids)) {
    d <- e[e$PanelSource == ids[i], ]
    x <- (i - 1) * 44
    txt(heads[i], x, 42.3, 6.8, fontface = "bold")
    txt(paste0("n = ", unique(d$TotalN)), x + 41, 42.3, 6.2, just = "right")
    txt(subs[i], x, 38.8, 6.4)
    txt("Class", x + 0.5, 34.4, 6.1)
    for (j in 1:3) {
      rect_box(x + 9 + (j - 1) * 8, 32.7, 7.9, 3.4, subtype_colors[paste0("C", j)])
      txt(paste0("C", j), x + 12.95 + (j - 1) * 8, 34.4, 6.2, fontface = "bold", just = "centre")
    }
    txt("Row n", x + 38.1, 34.4, 5.8, just = "centre")
    
    rows <- unique(d[order(d$RowOrder), c("OriginalLabel", "RowOrder")])
    cellh <- 23 / nrow(rows)
    for (k in seq_len(nrow(rows))) {
      r <- d[d$OriginalLabel == rows$OriginalLabel[k], ]
      y <- 32 - cellh * (k - 0.5)
      label <- sub("class: ", "", rows$OriginalLabel[k])
      label <- gsub("Not reported", "Not\nreported", label)
      txt(label, x + 4.3, y, 6.0, just = "centre")
      for (j in 1:3) {
        z <- r[r$FinalSubtype == paste0("C", j), ]
        rgbv <- shade(z$RowFraction) / 255
        fill <- if (grepl("Not", label)) "#F0F1F2" else rgb(rgbv[1], rgbv[2], rgbv[3])
        rect_box(x + 9 + (j - 1) * 8, y - cellh / 2 + 0.08, 7.9, cellh - 0.16, fill)
        txt(sprintf("%d(%d%%)", z$N, round(z$RowFraction * 100)),
            x + 12.95 + (j - 1) * 8, y, 6.2,
            col = if (z$RowFraction > 0.65 && !grepl("Not", label)) "white" else "#222222",
            just = "centre")
      }
      txt(unique(r$RowN), x + 38.1, y, 6.2, just = "centre")
    }
    txt("Total", x + 0.3, 6.5, 6.0)
    for (j in 1:3) {
      txt(unique(d$ColumnN[d$FinalSubtype == paste0("C", j)]),
          x + 12.95 + (j - 1) * 8, 6.5, 6.2, fontface = "bold", just = "centre")
    }
    txt(unique(d$TotalN), x + 38.1, 6.5, 6.2, fontface = "bold", just = "centre")
  }
  txt("Cells: n (% within original class). Dong RNA/protein panels share the same 223 samples.",
      0, 2.5, 6.0, col = "#555555")
  popViewport()
}

out_panelE_png <- file.path(out_dir, "Figure_1_PanelE.png")
out_panelE_pdf <- file.path(out_dir, "Figure_1_PanelE.pdf")
agg_png(out_panelE_png, width = 174, height = 44, units = "mm", res = 600, background = "white")
draw_panelE()
dev.off()
pdf(out_panelE_pdf, width = 174 / 25.4, height = 44 / 25.4, family = "Arial")
draw_panelE()
dev.off()

# ------------------------------------------------------------------------------
# Panel F: Clinicopathological Characteristics
# ------------------------------------------------------------------------------
cat("Rendering Panel F...\n")
draw_panelF <- function() {
  d <- read.csv(file.path(data_dir, "Figure1F_compact.csv"), stringsAsFactors = FALSE)
  grid.newpage()
  pushViewport(viewport(x = 0, y = 0, width = mm(90), height = mm(52.1), just = c("left", "bottom")))
  
  percent_col <- circlize::colorRamp2(c(0, 50, 100), c("#F8FAFC", "#D5E7F0", "#246085"))
  x <- c(36.5, 57.2, 78.0)
  ys <- c(38.3, 28.9, 19.5, 10.1)
  w <- 19.8; h <- 8.3
  
  txt("Clinical characteristics", 2.2, 50.2, 8.5, fontface = "bold", just = "left")
  for (j in 1:3) {
    rect_box(x[j] - w/2, 45.3 - 3.9/2, w, 3.9, subtype_colors[paste0("C", j)])
    txt(paste0("C", j), x[j], 45.3, 7.8, fontface = "bold", col = "#111111")
  }
  
  row_titles <- c("Tumor size", "Stage III/IV", "Vascular invasion", "Large-duct type")
  row_notes <- c("Median [IQR], cm", "Observed proportion", "Observed proportion", "Observed pathology")
  
  for (i in 1:4) {
    ss <- d[d$TraitOrder == i, ]
    yy <- ys[i]
    txt(row_titles[i], 2.2, yy + 1.8, 7.2, fontface = "bold", just = "left")
    txt(row_notes[i], 2.2, yy - 0.3, 5.8, just = "left")
    txt(paste0("Available n = ", sum(ss$AvailableN)), 2.2, yy - 2.3, 5.6, col = "#666666", just = "left")
    
    for (j in 1:3) {
      r <- ss[ss$FinalSubtype == paste0("C", j), ]
      fill <- if (i == 1) "#F5F6F8" else percent_col(r$Value)
      rect_box(x[j] - w/2, yy - h/2, w, h, fill)
      txt(r$Primary, x[j], yy + 1.4, 7.8, fontface = "bold")
      txt(r$Secondary, x[j], yy - 0.8, 6.2)
      if (i == 1) txt(r$Tertiary, x[j], yy - 2.8, 5.8, col = "#666666")
    }
  }
  
  txt("Proportion (%)", 2.2, 2.6, 5.8, just = "left")
  for (k in 1:101) {
    rect_box(27.5 + (k - 0.5) * 31.5 / 101, 2.2, 31.5 / 101, 1.2, percent_col(k - 1))
  }
  txt("0", 27.5, 0.9, 5.5)
  txt("50", 43.2, 0.9, 5.5)
  txt("100", 59.0, 0.9, 5.5)
  txt("n/N: positive / available", 88.0, 2.6, 5.8, just = "right")
  popViewport()
}

out_panelF_png <- file.path(out_dir, "Figure_1_PanelF.png")
out_panelF_pdf <- file.path(out_dir, "Figure_1_PanelF.pdf")
agg_png(out_panelF_png, width = 90, height = 52.1, units = "mm", res = 600, background = "white")
draw_panelF()
dev.off()
pdf(out_panelF_pdf, width = 90 / 25.4, height = 52.1 / 25.4, family = "Arial")
draw_panelF()
dev.off()

# ------------------------------------------------------------------------------
# Panel G: Overall Survival KM Curve
# ------------------------------------------------------------------------------
cat("Rendering Panel G...\n")
draw_panelG <- function() {
  curve <- read.csv(file.path(data_dir, "Figure1G_KM_coordinates.csv"), stringsAsFactors = FALSE)
  risk <- read.csv(file.path(data_dir, "Figure1G_risk_counts.csv"), stringsAsFactors = FALSE)
  
  grid.newpage()
  pushViewport(viewport(x = 0, y = 0, width = mm(76), height = mm(53), just = c("left", "bottom")))
  
  txt("Overall survival", 0, 51.2, 7.5, fontface = "bold", just = "left")
  txt("n = 331; events = 153", 0, 47.5, 6.4, just = "left")
  txt("Study-stratified log-rank P < 0.001", 0, 44.1, 6.4, just = "left")
  
  for (j in 1:3) {
    x <- 42 + (j - 1) * 10
    ln(x, 48.1, x + 3, 48.1, subtype_colors[paste0("C", j)], lwd = 1)
    txt(paste0("C", j), x + 3.7, 48.1, 6.3, just = "left")
  }
  
  xx <- function(v) 10 + v / 96 * 63
  yy <- function(v) 21.5 + v * 19
  
  ln(10, 21.5, 73, 21.5, col = "#555555", lwd = 0.5)
  ln(10, 21.5, 10, 40.5, col = "#555555", lwd = 0.5)
  
  for (v in c(0, 0.25, 0.5, 0.75, 1)) {
    ln(9.2, yy(v), 10, yy(v), col = "#555555", lwd = 0.5)
    txt(sprintf("%.2f", v), 8.5, yy(v), 6.0, just = "right")
  }
  
  times <- seq(0, 96, 12)
  for (t in times) {
    ln(xx(t), 21.5, xx(t), 20.8, col = "#555555", lwd = 0.5)
    txt(t, xx(t), 19.1, 6.1, just = "centre")
  }
  
  txt("Survival probability", 1.8, 30.5, 6.3, rot = 90, just = "centre")
  txt("Overall survival (months)", 41.5, 15.6, 6.3, just = "centre")
  
  for (c in paste0("C", 1:3)) {
    q <- curve[curve$Subtype == c, ]
    prev <- c(1, head(q$Survival, -1))
    
    # Step function segments
    step_x <- xx(c(0, rep(q$Month, each = 2)))
    step_y <- yy(c(1, as.vector(rbind(prev, q$Survival))))
    for (idx in seq(1, length(step_x) - 1)) {
      ln(step_x[idx], step_y[idx], step_x[idx + 1], step_y[idx + 1], subtype_colors[c], lwd = 0.85)
    }
    
    # Censored ticks
    cen <- q[q$Censored > 0, ]
    if (nrow(cen) > 0) {
      for (k in seq_len(nrow(cen))) {
        ln(xx(cen$Month[k]), yy(cen$Survival[k]) - 0.3,
           xx(cen$Month[k]), yy(cen$Survival[k]) + 0.3,
           subtype_colors[c], lwd = 0.5)
      }
    }
  }
  
  txt("Number at risk", 0, 12.5, 6.3, just = "left")
  for (j in 1:3) {
    y <- 9.5 - (j - 1) * 3
    c <- paste0("C", j)
    txt(c, 0, y, 6.3, col = subtype_colors[c], just = "left")
    r <- risk[risk$Subtype == c, ]
    for (k in seq_len(nrow(r))) {
      txt(r$AtRisk[k], xx(r$Month[k]), y, 6.2, just = "centre")
    }
  }
  popViewport()
}

out_panelG_png <- file.path(out_dir, "Figure_1_PanelG.png")
out_panelG_pdf <- file.path(out_dir, "Figure_1_PanelG.pdf")
agg_png(out_panelG_png, width = 76, height = 53, units = "mm", res = 600, background = "white")
draw_panelG()
dev.off()
pdf(out_panelG_pdf, width = 76 / 25.4, height = 53 / 25.4, family = "Arial")
draw_panelG()
dev.off()

# ------------------------------------------------------------------------------
# Full Composite: Figure 1 (180 mm x 270 mm)
# ------------------------------------------------------------------------------
cat("Assembling Full Figure 1 Composite...\n")
composite_png <- file.path(out_dir, "Figure_1_composite.png")
composite_pdf <- file.path(out_dir, "Figure_1_composite.pdf")

render_composite <- function() {
  grid.newpage()
  # Total canvas: 180 mm x 270 mm
  pushViewport(viewport(x = 0, y = 0, width = mm(180), height = mm(270), just = c("left", "bottom")))
  
  # White background
  grid.rect(gp = gpar(fill = "white", col = NA))
  
  # Helper to draw a panel raster from saved PNG
  draw_raster_panel <- function(png_path, x_mm, y_top_mm, w_mm, h_mm) {
    if (file.exists(png_path)) {
      img <- readPNG(png_path)
      y_bottom_mm <- 270 - y_top_mm - h_mm
      grid.raster(img, x = mm(x_mm + w_mm / 2), y = mm(y_bottom_mm + h_mm / 2),
                  width = mm(w_mm), height = mm(h_mm))
    }
  }
  
  # Helper to draw panel letter
  draw_panel_letter <- function(letter, x_mm, y_top_mm) {
    y_from_bottom <- 270 - y_top_mm
    grid.text(letter, x = mm(x_mm), y = mm(y_from_bottom), just = c("left", "top"),
              gp = gpar(fontfamily = "Arial", fontsize = 11, fontface = "bold", col = "#111111"))
  }
  
  # Coordinates from spec.json:
  # Panel A: [4.1, 9, 173.22, 64.29] | Letter at [0.3, 8.2]
  draw_raster_panel(out_panelA_png, 4.1, 9, 173.2, 64.3)
  draw_panel_letter("A", 1.5, 8.2)
  
  # Panel B: [4, 77, 114, 81.86] | Letter at [0.3, 76]
  draw_raster_panel(out_panelB_png, 4.0, 77, 114.0, 81.9)
  draw_panel_letter("B", 1.5, 76.0)
  
  # Panel C: [124, 77, 54, 45.83] | Letter at [120, 76]
  draw_raster_panel(out_panelC_png, 124.0, 77, 54.0, 45.8)
  draw_panel_letter("C", 120.0, 76.0)
  
  # Panel D: [124, 124, 54, 40.12] | Letter at [120, 123]
  draw_raster_panel(out_panelD_png, 124.0, 124, 54.0, 40.1)
  draw_panel_letter("D", 120.0, 123.0)
  
  # Panel E: [4, 165.5, 174, 44] | Letter at [0.3, 164.7]
  draw_raster_panel(out_panelE_png, 4.0, 165.5, 174.0, 44.0)
  draw_panel_letter("E", 1.5, 164.7)
  
  # Panel F: [4, 213, 90, 52.1] | Letter at [0.3, 212]
  draw_raster_panel(out_panelF_png, 4.0, 213, 90.0, 52.1)
  draw_panel_letter("F", 1.5, 212.0)
  
  # Panel G: [102, 213, 76, 53] | Letter at [98, 212]
  draw_raster_panel(out_panelG_png, 102.0, 213, 76.0, 53.0)
  draw_panel_letter("G", 98.0, 212.0)
  
  popViewport()
}

agg_png(composite_png, width = 180, height = 270, units = "mm", res = 300, background = "white")
render_composite()
dev.off()

pdf(composite_pdf, width = 180 / 25.4, height = 270 / 25.4, family = "Arial")
render_composite()
dev.off()

cat("Figure 1 rendered successfully:\n")
cat("  - Individual panels:", out_dir, "\n")
cat("  - Composite PDF:", composite_pdf, "\n")
cat("  - Composite PNG:", composite_png, "\n")
