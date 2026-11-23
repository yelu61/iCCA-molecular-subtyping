# Statistical labels; shared configuration, no panel letters or added tests.
fdr_label <- function(q, style=STEP8_CONFIG$fdr_style) {
  if (!style %in% c("stars","exact")) stop("Unknown FDR label style")
  if (any(!is.na(q) & (!is.finite(q) | q < 0 | q > 1))) stop("FDR must be in [0,1] or NA")
  out <- rep("NA",length(q)); keep <- !is.na(q)
  if (style == "exact") out[keep] <- ifelse(q[keep] < .001,formatC(q[keep],format="e",digits=1),formatC(q[keep],format="f",digits=3))
  else out[keep] <- ifelse(q[keep] >= .05,"ns",ifelse(q[keep] >= .01,"*",ifelse(q[keep] >= .001,"**","***")))
  out
}
ci_label <- function(est,lo,hi,digits=2) {
  if (any(!is.na(lo) & !is.na(hi) & lo > hi)) stop("Reversed CI endpoints")
  fmt <- function(x) ifelse(is.na(x),"NA",formatC(x,format="f",digits=digits))
  paste0(fmt(est)," (",fmt(lo)," to ",fmt(hi),")")
}
# Use a caption instead of covering data points at Inf coordinates.
n_tag <- function(text) ggplot2::labs(caption=text)
auc_label <- function(auc,lo,hi,digits=3) ci_label(auc,lo,hi,digits)
# Compatibility helper: intentionally draws no letter for standalone PDFs.
panel_tag <- function(letter) ggplot2::theme(plot.tag=ggplot2::element_blank())
calibration_caption_stats <- function(brier,intercept,slope,digits=3) {
  paste0("Brier = ",formatC(brier,format="f",digits=digits),"; intercept = ",
         formatC(intercept,format="f",digits=digits),"; slope = ",formatC(slope,format="f",digits=digits))
}
