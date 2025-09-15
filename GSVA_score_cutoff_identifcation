#!/usr/bin/env Rscript
# ============================================================
# Disease-severity cutoff discovery + plots from GSVA scores
# Author: Sumeet R. Deshmukh - 2025-09-15
# ------------------------------------------------------------
# Inputs (CSV):
#   1) metadata: must contain columns Sample.ID, Condition
#   2) scores:   GSVA pathway scores.
#                Supported shapes:
#                 a) wide matrix with a 'pathway' column and Samples as columns
#                 b) long/flat table with columns: Sample.ID, GSVA_Score
#
# What it does:
#   - Merges meta + chosen pathway scores
#   - Learns best cutoff from supervised contrasts (Youden) or unsupervised fallback
#   - Saves ranked contrasts + per-contrast density/ROC plots
#   - Produces an all-groups density with arrowed peaks
#
# Usage (edit CONFIG below or pass via env vars):
#   Rscript scripts/disease_severity_cutoff.R
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(pROC)
  library(effsize)
  library(cowplot)
  library(grid)
  # prefer dplyr verbs over similarly-named S4 methods
  if (requireNamespace("conflicted", quietly = TRUE)) {
    conflicted::conflicts_prefer(
      dplyr::filter, dplyr::select, dplyr::slice, dplyr::lag, dplyr::intersect
    )
  }
})

# -----------------------
# CONFIG (edit here)
# -----------------------
META_PATH      <- Sys.getenv("META_PATH",      unset = "data/Colon_metadata_usti.csv")
SCORES_PATH    <- Sys.getenv("SCORES_PATH",    unset = "data/Colon_usti_scores.csv")
PATHWAY_NAME   <- Sys.getenv("PATHWAY_NAME",   unset = "CDi_disease_severity_genes")
OUTDIR         <- Sys.getenv("OUTDIR",         unset = "outputs")
SEED           <- as.integer(Sys.getenv("SEED", unset = "1"))

set.seed(SEED)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "contrast_plots"), showWarnings = FALSE)

cat("\n[INFO] Reading files:\n",
    "  meta   :", META_PATH, "\n",
    "  scores :", SCORES_PATH, "\n",
    "  pathway:", PATHWAY_NAME, "\n\n")

# -----------------------
# READERS
# -----------------------
read_metadata <- function(path) {
  meta <- readr::read_csv(path, show_col_types = FALSE)
  # normalize common column names
  nm <- names(meta)
  nm <- sub("^sample[._ ]?id$", "Sample.ID", nm, ignore.case = TRUE)
  nm <- sub("^condition$",        "Condition", nm, ignore.case = TRUE)
  names(meta) <- nm
  stopifnot(all(c("Sample.ID","Condition") %in% names(meta)))
  meta$Sample.ID <- as.character(meta$Sample.ID)
  meta$Condition <- as.character(meta$Condition)
  meta
}

# scores reader that accepts:
#  (A) wide table with 'pathway' + sample columns
#  (B) already-flat table with Sample.ID + GSVA_Score
read_scores <- function(path, pathway_name) {
  df <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  nm <- names(df)

  # Flat format?
  if (all(c("Sample.ID","GSVA_Score") %in% nm)) {
    df <- df %>%
      mutate(Sample.ID = as.character(Sample.ID),
             GSVA_Score = suppressWarnings(as.numeric(GSVA_Score))) %>%
      filter(is.finite(GSVA_Score))
    return(df)
  }

  # Wide format with 'pathway' row key?
  has_path_col <- any(tolower(nm) == "pathway")
  if (!has_path_col) {
    stop("Scores file doesn't have 'Sample.ID+GSVA_Score' columns or a 'pathway' column.")
  }

  names(df)[tolower(names(df)) == "pathway"] <- "pathway"
  if (!(pathway_name %in% df$pathway)) {
    stop(sprintf("Pathway '%s' not found in scores.", pathway_name))
  }

  # keep only the requested pathway row and pivot to long
  one <- df %>% filter(.data$pathway == pathway_name)
  one$pathway <- NULL
  # now columns are sample IDs; pivot into two columns
  tibble(
    Sample.ID  = names(one),
    GSVA_Score = as.numeric(one[1, ] %>% unlist(use.names = FALSE))
  ) %>% filter(is.finite(GSVA_Score))
}

# -----------------------
# LOAD & MERGE
# -----------------------
meta <- read_metadata(META_PATH)
sev  <- read_scores(SCORES_PATH, PATHWAY_NAME)

df <- sev %>%
  inner_join(meta, by = "Sample.ID") %>%
  mutate(
    Condition = trimws(as.character(Condition)),
    GSVA_Score = suppressWarnings(as.numeric(GSVA_Score))
  )

readr::write_csv(df, file.path(OUTDIR, "data_input_all.csv"))
cat("[OK] Merged data →", file.path(OUTDIR, "data_input_all.csv"), "\n")

stopifnot(all(c("GSVA_Score","Condition") %in% names(df)))

# -----------------------
# LABEL MASKS (edit synonyms here)
# -----------------------
`%in_any%` <- function(x, patterns) {
  p <- paste0("^(?:", paste(patterns, collapse="|"), ")$")
  grepl(p, x, ignore.case = TRUE)
}

is_t_resp      <- df$Condition %in_any% c("Treated_Responder","Responder","Treated responder")
is_t_nonresp   <- df$Condition %in_any% c("Treated_non_responders","Treated_non_responder","Non-responder","Non_Responder")
is_t_uninfl    <- df$Condition %in_any% c("Treated_CD_Uninflamed","Treated Uninflamed")
is_t_infl      <- df$Condition %in_any% c("Treated_CD_Inflamed","Treated Inflamed")
is_nt_uninfl   <- df$Condition %in_any% c("Not_treated_CD_Uninflamed","Untreated_CD_Uninflamed","Not Treated CD Uninflamed","Untreated Uninflamed")
is_nt_infl     <- df$Condition %in_any% c("Not_treated_CD_Inflamed","Untreated_CD_Inflamed","Not Treated CD Inflamed","Untreated Inflamed")

uninf_like_mask <- is_t_uninfl | is_nt_uninfl
infl_like_mask  <- is_t_infl   | is_nt_infl

have_treated_labels <- sum(is_t_resp,    na.rm=TRUE) >= 2 && sum(is_t_nonresp, na.rm=TRUE) >= 2   # A
have_treated_pair   <- sum(is_t_uninfl,  na.rm=TRUE) >= 2 && sum(is_t_infl,    na.rm=TRUE) >= 2   # B2
have_untreated_pair <- sum(is_nt_uninfl, na.rm=TRUE) >= 2 && sum(is_nt_infl,   na.rm=TRUE) >= 2   # B1
have_combined       <- sum(uninf_like_mask, na.rm=TRUE) >= 2 && sum(infl_like_mask, na.rm=TRUE) >= 2  # B*

# -----------------------
# HELPERS
# -----------------------
unsupervised_cut <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4) {
    return(list(cutoff = median(x, na.rm=TRUE),
                auc = NA_real_,
                orient = "lower_score = positive"))
  }
  d <- density(x, na.rm = TRUE)
  y <- d$y; xg <- d$x
  dy <- diff(y)
  minima_idx <- which(dy[-1] > 0 & dy[-length(dy)] < 0) + 1
  if (length(minima_idx) == 0L) {
    cutoff <- median(x, na.rm=TRUE)
  } else {
    q <- quantile(x, c(0.05, 0.95), na.rm = TRUE)
    in_win <- xg[minima_idx] >= q[1] & xg[minima_idx] <= q[2]
    if (!any(in_win)) in_win <- rep(TRUE, length(minima_idx))
    cand <- minima_idx[in_win]
    cutoff <- xg[cand[which.min(y[cand])]]
  }
  list(cutoff = as.numeric(cutoff), auc = NA_real_, orient = "lower_score = positive")
}

eval_contrast <- function(x_pos, x_neg, name, pos_lab, neg_lab) {
  x_pos <- x_pos[is.finite(x_pos)]
  x_neg <- x_neg[is.finite(x_neg)]
  n_pos <- length(x_pos); n_neg <- length(x_neg)

  if (n_pos < 2 || n_neg < 2) {
    return(tibble(
      contrast = name, pos = pos_lab, neg = neg_lab,
      n_pos = n_pos, n_neg = n_neg,
      auc = NA_real_, cutoff = NA_real_, orient = NA_character_,
      youden = NA_real_, ks_D = NA_real_, cliffs_delta = NA_real_
    ))
  }

  # lower GSVA = healthier/positive → negate so "higher=positive" for ROC
  y <- c(rep(1, n_pos), rep(0, n_neg))
  x <- c(-x_pos, -x_neg)

  roc_obj <- pROC::roc(response = y, predictor = x, quiet = TRUE)
  best    <- pROC::coords(roc_obj, "best", best.method = "youden", transpose = TRUE)
  cutoff  <- -as.numeric(best["threshold"])   # back to original scale
  auc     <- as.numeric(pROC::auc(roc_obj))
  youden  <- as.numeric(best["youden"])
  ks      <- suppressWarnings(stats::ks.test(x_pos, x_neg))
  ks_D    <- unname(ks$statistic)
  cd      <- tryCatch(effsize::cliff.delta(x_pos, x_neg)$estimate, error = function(e) NA_real_)

  tibble(
    contrast = name, pos = pos_lab, neg = neg_lab,
    n_pos = n_pos, n_neg = n_neg,
    auc = auc, cutoff = cutoff,
    orient = "lower_score = positive",
    youden = youden, ks_D = ks_D, cliffs_delta = as.numeric(cd)
  )
}

find_peak <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(list(x = NA_real_, y = NA_real_))
  d <- stats::density(x); i <- which.max(d$y)
  list(x = d$x[i], y = d$y[i])
}

# -----------------------
# SUPERVISED: try A/B2/B1/B*
# -----------------------
rows <- list()
add_row <- function(flag, x_pos, x_neg, name, pos_lab, neg_lab) {
  if (isTRUE(flag)) rows[[length(rows)+1]] <<- eval_contrast(x_pos, x_neg, name, pos_lab, neg_lab)
}

add_row(have_treated_labels, df$GSVA_Score[is_t_resp],  df$GSVA_Score[is_t_nonresp],
        "A",  "Treated_Responder",          "Treated_non_responders")
add_row(have_treated_pair,   df$GSVA_Score[is_t_uninfl],df$GSVA_Score[is_t_infl],
        "B2", "Treated_CD_Uninflamed",      "Treated_CD_Inflamed")
add_row(have_untreated_pair, df$GSVA_Score[is_nt_uninfl], df$GSVA_Score[is_nt_infl],
        "B1", "Not_treated_CD_Uninflamed",  "Not_treated_CD_Inflamed")
add_row(have_combined,       df$GSVA_Score[uninf_like_mask], df$GSVA_Score[infl_like_mask],
        "B*", "Healthy-like (Healthy + Uninflamed)", "Inflamed-like")

results <- dplyr::bind_rows(rows)

if (nrow(results) == 0) {
  message("[WARN] No supervised contrasts available; will fall back to unsupervised.")
} else {
  results <- results %>%
    dplyr::mutate(min_n = pmin(n_pos, n_neg),
                  abs_cliff = abs(cliffs_delta)) %>%
    dplyr::arrange(dplyr::desc(auc), dplyr::desc(abs_cliff),
                   dplyr::desc(ks_D), dplyr::desc(min_n))
  print(results)

  # Pick best
  best      <- dplyr::slice_head(results, n = 1)
  cutoff    <- best$cutoff[1]
  auc_val   <- best$auc[1]
  orient    <- best$orient[1]
  learned_on <- sprintf("Path %s: %s vs %s", best$contrast[1], best$pos[1], best$neg[1])

  message(sprintf("[OK] Best contrast: %s | cutoff=%.3f | AUC=%.3f | orient=%s",
                  learned_on, cutoff, auc_val, orient))
}

# -----------------------
# UNSUPERVISED fallback
# -----------------------
if (!exists("cutoff") || !is.finite(cutoff)) {
  message("[INFO] Falling back to unsupervised valley/median.")
  fit <- unsupervised_cut(df$GSVA_Score)
  cutoff    <- fit$cutoff
  auc_val   <- fit$auc
  orient    <- fit$orient
  learned_on <- "Path C: Unsupervised (valley/median)"
}

# -----------------------
# APPLY cutoff to treated samples (binary label)
# -----------------------
treated_mask <- df$Condition %in% c("Treated_CD_Uninflamed", "Treated_CD_Inflamed",
                                    "Treated_Responder", "Treated_non_responders")

if (identical(orient, "higher_score = positive")) {
  df$Predicted_Response <- ifelse(treated_mask & is.finite(df$GSVA_Score) & df$GSVA_Score >= cutoff, 1,
                           ifelse(treated_mask & is.finite(df$GSVA_Score), 0, NA))
} else {
  df$Predicted_Response <- ifelse(treated_mask & is.finite(df$GSVA_Score) & df$GSVA_Score <= cutoff, 1,
                           ifelse(treated_mask & is.finite(df$GSVA_Score), 0, NA))
}
df$Predicted_Response_Label <- ifelse(is.na(df$Predicted_Response), NA,
                               ifelse(df$Predicted_Response == 1, "Responder", "Non-responder"))

readr::write_csv(results, file.path(OUTDIR, "Contrast_ranking_A_B2_B1_Bstar.csv"))
message("[OK] Wrote: ", file.path(OUTDIR, "Contrast_ranking_A_B2_B1_Bstar.csv"))

# -----------------------
# DENSITY: all groups with cutoff & arrows to each peak
# -----------------------
present_lvls <- unique(stats::na.omit(as.character(df$Condition)))
df_all <- df %>%
  dplyr::filter(Condition %in% present_lvls & is.finite(GSVA_Score)) %>%
  dplyr::mutate(Condition = factor(Condition, levels = present_lvls))

pal_distinct <- c(
  "#E41A1C","#377EB8","#4DAF4A","#984EA3",
  "#FF7F00","#A65628","#F781BF","#999999"
)
pal_all <- setNames(rep(pal_distinct, length.out = length(levels(df_all$Condition))),
                    levels(df_all$Condition))

peaks_all <- df_all %>%
  group_by(Condition) %>%
  summarise(.p = list(find_peak(GSVA_Score)), .groups="drop") %>%
  mutate(peak_x = vapply(.p, `[[`, numeric(1), "x"),
         peak_y = vapply(.p, `[[`, numeric(1), "y")) %>%
  dplyr::select(-.p) %>%
  dplyr::filter(is.finite(peak_x), is.finite(peak_y))

ymax_all <- if (nrow(peaks_all)) max(peaks_all$peak_y) else 1

p_all <- ggplot(df_all, aes(x = GSVA_Score, fill = Condition, color = Condition)) +
  geom_density(alpha = 0.35, linewidth = 0.6, adjust = 1) +
  geom_vline(xintercept = cutoff, linetype = "dashed", linewidth = 0.8) +
  geom_curve(
    data = peaks_all,
    aes(x = ifelse(peak_x < cutoff, peak_x - 0.12, peak_x + 0.12),
        y = pmin(peak_y + 0.15 * ymax_all, ymax_all * 1.05),
        xend = peak_x, yend = peak_y, color = Condition),
    inherit.aes = FALSE, curvature = 0.2, linewidth = 0.6,
    arrow = arrow(length = unit(3, "mm"), type = "closed")
  ) +
  geom_label(
    data = peaks_all,
    aes(x = ifelse(peak_x < cutoff, peak_x - 0.12, peak_x + 0.12),
        y = pmin(peak_y + 0.15 * ymax_all, ymax_all * 1.05),
        label = as.character(Condition), fill = Condition),
    inherit.aes = FALSE, label.size = 0, alpha = 0.9, size = 4
  ) +
  scale_fill_manual(values = pal_all, drop = FALSE) +
  scale_color_manual(values = pal_all, guide = "none", drop = FALSE) +
  theme_minimal(base_size = 16) +
  labs(title = "Disease Severity Densities by Condition (arrows at peaks)",
       subtitle = paste0(learned_on,
                         if (!is.na(auc_val)) sprintf(" | AUC = %.3f", auc_val) else ""),
       x = "Disease Severity (GSVA score)", y = "Density") +
  theme(legend.position = "top", legend.title = element_blank())

ggsave(file.path(OUTDIR, "DiseaseSeverity_density_all_groups.png"), p_all, width = 11, height = 8, dpi = 300)
ggsave(file.path(OUTDIR, "DiseaseSeverity_density_all_groups.pdf"),  p_all, width = 11, height = 8)
message("[OK] Saved all-groups density plots.")

# -----------------------
# PER-CONTRAST density + ROC
# -----------------------
contrasts_avail <- list()
add_contrast <- function(flag, name, pos_mask, neg_mask, pos_lab, neg_lab){
  if (isTRUE(flag)) {
    contrasts_avail[[length(contrasts_avail)+1]] <<- list(
      name = name,
      pos_mask = pos_mask, neg_mask = neg_mask,
      pos_lab  = pos_lab,  neg_lab  = neg_lab
    )
  }
}

add_contrast(have_treated_labels, "A: Treated_Responder vs Treated_non_responders",
             is_t_resp, is_t_nonresp, "Treated_Responder", "Treated_non_responders")
add_contrast(have_treated_pair,   "B2: Treated_CD_Uninflamed vs Treated_CD_Inflamed",
             is_t_uninfl, is_t_infl, "Treated_CD_Uninflamed", "Treated_CD_Inflamed")
add_contrast(have_untreated_pair, "B1: Not_treated_CD_Uninflamed vs Not_treated_CD_Inflamed",
             is_nt_uninfl, is_nt_infl, "Not_treated_CD_Uninflamed", "Not_treated_CD_Inflamed")
add_contrast(have_combined,       "B*: Healthy-like vs Inflamed-like",
             uninf_like_mask, infl_like_mask, "Healthy-like", "Inflamed-like")

make_density_plot <- function(df2, cutoff, title, pos_lab, neg_lab) {
  pal <- c(Pos="#1F77B4", Neg="#FF7F0E")
  ggplot(df2, aes(x = GSVA_Score, fill = grp, color = grp)) +
    geom_density(alpha = 0.35, linewidth = 0.7, adjust = 1) +
    geom_vline(xintercept = cutoff, linetype = "dashed", linewidth = 0.9) +
    annotate("label", x = cutoff, y = Inf, vjust = 1.2,
             label = sprintf("Cutoff = %.3f", cutoff),
             size = 4, label.size = 0, fill = "white") +
    scale_fill_manual(values = pal, labels = c(Pos = pos_lab, Neg = neg_lab)) +
    scale_color_manual(values = pal, guide = "none") +
    labs(title = paste0("Density: ", title),
         x = "Disease Severity (GSVA score)", y = "Density", fill = NULL) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"))
}

make_roc_plot <- function(roc_obj, auc_val, title){
  df_roc <- data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities
  )
  ggplot(df_roc, aes(x = fpr, y = tpr)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
    geom_line(size = 1) +
    coord_equal() +
    labs(title = paste0("ROC: ", title, sprintf(" (AUC = %.3f)", auc_val)),
         x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
}

dir.create(file.path(OUTDIR, "contrast_plots/with_ROC"), recursive = TRUE, showWarnings = FALSE)
all_pages <- list()

compute_cutoff_for <- function(x_pos, x_neg){
  x_pos <- x_pos[is.finite(x_pos)]
  x_neg <- x_neg[is.finite(x_neg)]
  stopifnot(length(x_pos) >= 2, length(x_neg) >= 2)
  y <- c(rep(1, length(x_pos)), rep(0, length(x_neg)))
  x <- c(-x_pos, -x_neg)
  roc_obj <- pROC::roc(y, x, quiet = TRUE)
  best    <- pROC::coords(roc_obj, "best", best.method = "youden", transpose = TRUE)
  list(
    cutoff  = -as.numeric(best["threshold"]),
    auc     = as.numeric(pROC::auc(roc_obj)),
    roc_obj = roc_obj
  )
}

for (cc in contrasts_avail) {
  pos_vals <- df$GSVA_Score[cc$pos_mask]; pos_vals <- pos_vals[is.finite(pos_vals)]
  neg_vals <- df$GSVA_Score[cc$neg_mask]; neg_vals <- neg_vals[is.finite(neg_vals)]
  if (length(pos_vals) < 2 || length(neg_vals) < 2) {
    message("Skipping ", cc$name, " (insufficient samples)."); next
  }

  fit        <- compute_cutoff_for(pos_vals, neg_vals)
  cutoff_use <- fit$cutoff; auc_use <- fit$auc; roc_use <- fit$roc_obj

  df2 <- dplyr::bind_rows(
    tibble::tibble(GSVA_Score = pos_vals, grp = factor("Pos", levels = c("Pos","Neg"))),
    tibble::tibble(GSVA_Score = neg_vals, grp = factor("Neg", levels = c("Pos","Neg")))
  )

  p_den <- make_density_plot(df2, cutoff_use, cc$name, cc$pos_lab, cc$neg_lab)
  p_roc <- make_roc_plot(roc_use, auc_use, cc$name)
  combo <- cowplot::plot_grid(p_den, p_roc, ncol = 2, rel_widths = c(1.2, 1))

  base <- gsub("[^A-Za-z0-9_\\-]+", "_", cc$name)
  ggsave(file.path(OUTDIR, "contrast_plots/with_ROC", paste0(base, "_density.png")),
         p_den, width = 7.8, height = 5.6, dpi = 300)
  ggsave(file.path(OUTDIR, "contrast_plots/with_ROC", paste0(base, "_roc.png")),
         p_roc, width = 7.2, height = 5.6, dpi = 300)
  ggsave(file.path(OUTDIR, "contrast_plots/with_ROC", paste0(base, "_combo.png")),
         combo, width = 13.8, height = 5.8, dpi = 300)

  all_pages[[length(all_pages)+1]] <- combo
}

if (length(all_pages)) {
  pdf(file.path(OUTDIR, "contrast_plots/with_ROC", "ALL_contrasts_combo.pdf"),
      width = 13.8, height = 5.8, onefile = TRUE)
  for (g in all_pages) grid::grid.newpage()
  for (g in all_pages) grid::grid.draw(ggplotGrob(g))
  dev.off()
  message("[OK] Wrote per-contrast PNGs and ALL_contrasts_combo.pdf")
} else {
  message("[INFO] No contrast plots were produced.")
}

cat("\n[DONE] Outputs in:", normalizePath(OUTDIR), "\n")
