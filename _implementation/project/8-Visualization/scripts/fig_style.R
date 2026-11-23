# Shared Step8 style: governed palettes, final physical dimensions, PDF only.
suppressPackageStartupMessages(library(ggplot2))
find_project_root <- function(start = getwd()) {
  p <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(p, "governance/data/subtype_master_table.csv")) &&
        file.exists(file.path(p, "8-Visualization/scripts/figure_config.json"))) return(p)
    parent <- dirname(p)
    if (identical(parent, p)) stop("Cannot locate the iCCA project root")
    p <- parent
  }
}
STEP8_ROOT <- find_project_root()
STEP8_CONFIG <- jsonlite::fromJSON(file.path(STEP8_ROOT, "8-Visualization/scripts/figure_config.json"))
STEP8_MM <- STEP8_CONFIG$widths_mm
STEP8_FONT <- list(family=STEP8_CONFIG$font_family, panel=STEP8_CONFIG$font_title,
                  axis=STEP8_CONFIG$font_axis, tick=STEP8_CONFIG$font_tick, note=STEP8_CONFIG$font_note)
SUBTYPE_PALETTE_CSV <- "governance/data/subtype_label_color_dictionary.csv"
STUDY_PALETTE_CSV <- "0-Data/qc/study_color_dictionary.csv"
ASSAY_PALETTE_CSV <- "governance/data/assay_type_color_dictionary.csv"
SUBTYPE_TRUTH_CSV <- "governance/data/subtype_master_table.csv"
MISSING_UNKNOWN_FILL <- STEP8_CONFIG$missing$unknown
MISSING_NOT_ASSAYED <- STEP8_CONFIG$missing$not_assayed
MISSING_NOT_ESTIMABLE <- STEP8_CONFIG$missing$not_estimable
SCALE_EFFECT_LOW <- STEP8_CONFIG$effect_colors[1]
SCALE_EFFECT_MID <- STEP8_CONFIG$effect_colors[2]
SCALE_EFFECT_HIGH <- STEP8_CONFIG$effect_colors[3]
# unknown=grey/?; not_assayed=dark grey/x; not_estimable=white/NE.
# Measured zero remains a measured value, never one of these missing states.
missing_encoding <- function(state) {
  states <- c("unknown", "not_assayed", "not_estimable")
  if (any(!state %in% states)) stop("Unknown missing-data state")
  data.frame(state=state, fill=unlist(STEP8_CONFIG$missing)[state],
             mark=c(unknown="?", not_assayed="x", not_estimable="NE")[state], row.names=NULL)
}
load_palette_table <- function(root, rel, key) {
  df <- read.csv(file.path(root, rel), stringsAsFactors=FALSE, check.names=FALSE)
  color_col <- if ("ColorHex" %in% names(df)) "ColorHex" else "Color"
  required <- c(key, color_col)
  if (!nrow(df) || !all(required %in% names(df))) stop("Invalid palette schema: ", rel)
  if (anyNA(df[[key]]) || any(!nzchar(df[[key]])) || anyDuplicated(df[[key]])) stop("Invalid palette keys: ", rel)
  if (anyNA(df[[color_col]]) || any(!grepl("^#[0-9A-Fa-f]{6}$", df[[color_col]]))) stop("Invalid palette colors: ", rel)
  if ("DisplayOrder" %in% names(df)) df <- df[order(df$DisplayOrder), ]
  df$ColorHex <- df[[color_col]]
  df
}
load_subtype_palette <- function(root=STEP8_ROOT) {
  df <- load_palette_table(root, SUBTYPE_PALETTE_CSV, "FinalSubtype")
  stopifnot(all(c("DisplayLabel", "SubtypeClass") %in% names(df)))
  list(colors=setNames(df$ColorHex, df$FinalSubtype), labels=setNames(df$DisplayLabel, df$FinalSubtype),
       class=setNames(df$SubtypeClass, df$FinalSubtype), df=df)
}
load_study_palette <- function(root=STEP8_ROOT) {
  df <- load_palette_table(root, STUDY_PALETTE_CSV, "Study")
  stopifnot(all(c("DisplayOrder", "AnalysisRole") %in% names(df)))
  list(colors=setNames(df$ColorHex, df$Study), roles=setNames(df$AnalysisRole, df$Study), df=df)
}
load_assay_palette <- function(root=STEP8_ROOT) {
  df <- load_palette_table(root, ASSAY_PALETTE_CSV, "AnalysisRole")
  stopifnot(all(c("DisplayLabel") %in% names(df)))
  list(colors=setNames(df$ColorHex, df$AnalysisRole), labels=setNames(df$DisplayLabel, df$AnalysisRole), df=df)
}
scale_fill_subtype_main <- function(pal=load_subtype_palette()) {
  main <- pal$df[pal$df$SubtypeClass == "Main", ]
  scale_fill_manual(values=setNames(main$ColorHex, main$FinalSubtype), limits=main$FinalSubtype,
                    labels=setNames(main$DisplayLabel, main$FinalSubtype), drop=FALSE)
}
scale_color_subtype_main <- function(pal=load_subtype_palette()) {
  main <- pal$df[pal$df$SubtypeClass == "Main", ]
  scale_color_manual(values=setNames(main$ColorHex, main$FinalSubtype), limits=main$FinalSubtype,
                     labels=setNames(main$DisplayLabel, main$FinalSubtype), drop=FALSE)
}
theme_step8 <- function(base_size=STEP8_FONT$tick, base_family=STEP8_FONT$family) {
  theme_classic(base_size=base_size, base_family=base_family) + theme(
    plot.title=element_text(size=STEP8_FONT$panel, face="bold", hjust=0),
    plot.subtitle=element_text(size=STEP8_FONT$note),
    plot.caption=element_text(size=STEP8_FONT$note, hjust=0),
    plot.caption.position="plot", plot.title.position="plot",
    axis.title=element_text(size=STEP8_FONT$axis), axis.text=element_text(size=STEP8_FONT$tick, color="grey20"),
    legend.title=element_text(size=STEP8_FONT$tick), legend.text=element_text(size=STEP8_FONT$tick),
    strip.text=element_text(size=STEP8_FONT$tick, face="bold"),
    axis.line=element_line(color="grey40", linewidth=.3),
    legend.key.size=grid::unit(3,"mm"), plot.margin=margin(3,3,3,3,"mm"))
}
# draw may be a ggplot/grid object or a function drawing on the open PDF device.
save_panel <- function(plot, target_pdf, width_mm, height_mm, root=STEP8_ROOT) {
  if (!grepl("[.]pdf$", target_pdf, ignore.case=TRUE)) stop("Target must end in .pdf")
  if (any(!is.finite(c(width_mm,height_mm))) || min(width_mm,height_mm) <= 0) stop("Invalid PDF dimensions")
  if (inherits(plot,"ggplot")) {
    for (label in c("title","subtitle","caption")) if (is.character(plot$labels[[label]]))
      plot$labels[[label]] <- paste(strwrap(plot$labels[[label]],width=max(35,floor(width_mm/1.8))),collapse="\n")
  }
  out <- if (grepl("^/", target_pdf)) target_pdf else file.path(root,target_pdf)
  dir.create(dirname(out), showWarnings=FALSE, recursive=TRUE)
  temp <- tempfile(".step8-", tmpdir=dirname(out), fileext=".pdf")
  on.exit(unlink(temp), add=TRUE)
  grDevices::cairo_pdf(temp, width=width_mm/25.4, height=height_mm/25.4, family=STEP8_FONT$family)
  device <- grDevices::dev.cur()
  tryCatch({if (is.function(plot)) plot() else if (inherits(plot,"grob")) grid::grid.draw(plot) else print(plot)},
           finally=grDevices::dev.off(device))
  if (!file.exists(temp) || file.info(temp)$size == 0) stop("Empty PDF export")
  if (!file.rename(temp,out)) stop("Cannot publish PDF: ",out)
  invisible(normalizePath(out))
}
