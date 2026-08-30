# ══════════════════════════════════════════════════════════════════════════════
# Sweep Threshold Dropoff — Query-Coverage and Gene-Fraction
# ══════════════════════════════════════════════════════════════════════════════
#
# Reads combined_gene_detail.csv / combined_results.csv / combined_length_quantiles.csv
# (from combine_sweep_results.py). Each row carries a sample_id; there is no
# workflow / read_subset grouping and no faceting. Line plots are a single panel
# showing the mean across samples (bold) with each sample as a thin grey line.
#
# Heatmaps are drawn as TWO PANELS side by side: the absolute feature count on
# the left and the percentage of baseline retained on the right. The two need
# independent fill scales, which is why patchwork is used to combine them.
#
# The gene-fraction grid is read from combined_results.csv rather than hardcoded,
# so it always matches whatever was actually swept by coverage_threshold_sweep.py.
#
# Usage:
#   Rscript plot_sweep_dropoff.R
#   Rscript plot_sweep_dropoff.R combined_sweep_results/ figures/
# ══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})
setDTthreads(1)

# patchwork is used to place the count and percentage heatmaps side by side with
# independent fill scales. If it is unavailable the script still runs and writes
# the two panels as separate files instead.
HAVE_PATCHWORK <- requireNamespace("patchwork", quietly = TRUE)
if (HAVE_PATCHWORK) suppressPackageStartupMessages(library(patchwork))

# ── arguments ─────────────────────────────────────────────────────────────────
args    <- commandArgs(trailingOnly = TRUE)
in_dir  <- if (length(args) >= 1) args[1] else "combined_sweep_results/"
out_dir <- if (length(args) >= 2) args[2] else file.path(in_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gene_path <- file.path(in_dir, "combined_gene_detail.csv")
if (!file.exists(gene_path)) stop("combined_gene_detail.csv not found in ", in_dir,
                                  " - run combine_sweep_results.py first")

LEVELS_TO_PLOT <- c("gene_accession", "class", "mechanism", "group")

# Fallback only. This is overwritten below with the actual values found in
# combined_results.csv, so it no longer has to be kept in sync by hand.
GF_SWEEP <- c(0, 0.1, 0.25, 0.5, 0.8)

message("Loading data...")
genes <- fread(gene_path, showProgress = FALSE)
genes <- genes[sample_id != "sample_id"]   # drop any embedded header rows

qcov_values     <- sort(unique(as.double(genes$min_query_coverage)))
identity_values <- if ("min_identity" %in% names(genes))
  sort(unique(as.double(genes$min_identity))) else c(0)
match_qcov_values <- if ("min_match_qcov" %in% names(genes))
  sort(unique(as.double(genes$min_match_qcov))) else c(0)

# Priority for "which axis is primary": match_qcov > identity > qcov.
if (length(match_qcov_values) > 1) {
  primary_values <- match_qcov_values
  primary_col    <- "min_match_qcov"
  primary_label  <- "Min match-qcov"
} else if (length(identity_values) > 1) {
  primary_values <- identity_values
  primary_col    <- "min_identity"
  primary_label  <- "Min identity"
} else {
  primary_values <- qcov_values
  primary_col    <- "min_query_coverage"
  primary_label  <- "Min query coverage"
}
# All four sweep axes are now proportions (0-1), so no axis needs special
# scaling. PRIMARY_IS_IDENTITY is retained as FALSE for backward compatibility
# with any local edits that referenced it.
PRIMARY_IS_IDENTITY <- FALSE

# ── Guard: these figures assume ONE read-level axis was swept ─────────────────
# coverage_threshold_sweep.py evaluates the full cartesian product of
# query-coverage x identity x match-qcov x gene-fraction. Every figure below
# plots against a single "primary" axis. If more than one read-level axis was
# swept, the other axes have to be dealt with somehow; this script holds them at
# their most permissive value (see hold_non_primary below) so each figure means
# "detection vs the primary axis, with the other filters off". That is a
# well-defined statement, but it is NOT the same as exploring the joint grid, so
# say so loudly rather than letting the reader assume otherwise.
n_axes_swept <- sum(c(length(qcov_values) > 1,
                      length(identity_values) > 1,
                      length(match_qcov_values) > 1))
if (n_axes_swept > 1) {
  message("")
  message("  ================================================================")
  message("  NOTE: more than one read-level threshold axis was swept.")
  message("    query-coverage values : ", length(qcov_values))
  message("    identity values       : ", length(identity_values))
  message("    match-qcov values     : ", length(match_qcov_values))
  message("  Figures plot against '", primary_col, "' and hold the other axes")
  message("  at their most permissive value. The CSV outputs contain the full")
  message("  grid and remain the source of truth for joint effects.")
  message("  To plot a different axis, sweep that one alone.")
  message("  ================================================================")
  message("")
}

# Restrict a table to the most permissive value of every read-level threshold
# axis other than the one being plotted. Without this, per-level counts computed
# with uniqueN() would union across the ignored axes, and heatmap tiles would be
# drawn multiple times on top of each other.
# Describe any read-level filter that was pinned at a single NON-ZERO value.
# Such a filter is applied at every point in the grid but never appears on an
# axis, so without this the figures would silently omit it.
fixed_filter_note <- function(dt, keep_col) {
  notes <- character(0)
  for (ax in setdiff(c("min_query_coverage", "min_identity", "min_match_qcov"), keep_col)) {
    if (ax %in% names(dt)) {
      vals <- unique(as.double(dt[[ax]]))
      if (length(vals) == 1 && vals[1] > 0) {
        notes <- c(notes, paste0(ax, " fixed at ", vals[1]))
      }
    }
  }
  if (length(notes) == 0) "" else paste0("Also applied at every point: ",
                                         paste(notes, collapse = "; "))
}

hold_non_primary <- function(dt, keep_col, label) {
  for (ax in setdiff(c("min_query_coverage", "min_identity", "min_match_qcov"), keep_col)) {
    if (ax %in% names(dt)) {
      vals <- as.double(dt[[ax]])
      if (uniqueN(vals) > 1) {
        ax_min <- min(vals, na.rm = TRUE)
        message("  [", label, "] holding ", ax, " at its most permissive value (", ax_min, ")")
        dt <- dt[as.double(get(ax)) == ax_min]
      }
    }
  }
  dt
}

message(sprintf("  %s rows | %d samples",
                formatC(nrow(genes), format = "d", big.mark = ","),
                uniqueN(genes$sample_id)))
message(sprintf("  Primary sweep: %s - values: %s",
                primary_col, paste(primary_values, collapse = ", ")))

# ══════════════════════════════════════════════════════════════════════════════
# Quick summary table (per-sample + averaged across samples)
# ══════════════════════════════════════════════════════════════════════════════
results_path <- file.path(in_dir, "combined_results.csv")
if (file.exists(results_path)) {
  message("\nBuilding quick summary table from combined_results.csv...")
  results <- fread(results_path, showProgress = FALSE)
  results <- results[sample_id != "sample_id"]
  results[, n_fragment_hits_in_passing_genes := as.double(n_fragment_hits_in_passing_genes)]

  # Derive the gene-fraction grid from the data instead of hardcoding it, so the
  # figures always match whatever grid was passed to coverage_threshold_sweep.py.
  if ("min_gene_fraction" %in% names(results)) {
    results[, min_gene_fraction := as.double(min_gene_fraction)]
    GF_SWEEP <- sort(unique(results$min_gene_fraction))
    message("  Gene-fraction grid read from data: ", paste(GF_SWEEP, collapse = ", "))
  } else {
    message("  [NOTE] min_gene_fraction column absent; using default grid: ",
            paste(GF_SWEEP, collapse = ", "))
  }

  THRESHOLD_COLS <- intersect(c("min_query_coverage", "min_identity",
                                "min_match_qcov", "min_gene_fraction"),
                              names(results))
  message("  Threshold columns detected: ", paste(THRESHOLD_COLS, collapse = ", "))

  raw_cols  <- c("n_alignments_in_passing_genes", "n_fragment_hits_in_passing_genes",
                "n_genes_passing_both", "n_classes_passing",
                "n_mechanisms_passing", "n_groups_passing")
  new_names <- c("n_classified_alignments", "n_classified_fragment_hits",
                "n_genes", "n_classes", "n_mechanisms", "n_groups")
  present   <- raw_cols %in% names(results)
  raw_cols  <- raw_cols[present]
  new_names <- new_names[present]

  id_cols <- c("sample_id", THRESHOLD_COLS)
  per_sample_summary <- results[, ..id_cols]
  for (i in seq_along(raw_cols)) per_sample_summary[[new_names[i]]] <- results[[raw_cols[i]]]

  fwrite(per_sample_summary, file.path(out_dir, "summary_table_per_sample.csv"))
  message(sprintf("  saved: summary_table_per_sample.csv (%s rows - one per sample x threshold combination)",
                  formatC(nrow(per_sample_summary), format = "d", big.mark = ",")))

  summary_table_avg <- per_sample_summary[, c(
    lapply(.SD, mean), .(n_samples = .N)
  ), by = THRESHOLD_COLS, .SDcols = new_names]
  setnames(summary_table_avg, new_names, paste0("mean_", new_names))
  setorderv(summary_table_avg, THRESHOLD_COLS)

  fwrite(summary_table_avg, file.path(out_dir, "summary_table_averaged.csv"))
  message(sprintf("  saved: summary_table_averaged.csv (%d rows - one per threshold combination, averaged across samples)",
                  nrow(summary_table_avg)))
  message("\nPreview of summary_table_averaged.csv:")
  print(summary_table_avg)
} else {
  message("\n[NOTE] combined_results.csv not found in ", in_dir,
         " - skipping quick summary table.")
  message("  Using default gene-fraction grid: ", paste(GF_SWEEP, collapse = ", "))
}

show_plot <- function(p, stem, w = 9, h = 7) {
  print(p)
  ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = w, height = h, dpi = 150)
  message("  saved: ", stem, ".png")
}

# ── two-panel heatmap helper ───────────────────────────────────────────────────
# The count panel and the percentage panel need independent fill scales, since a
# raw feature count and a 0-100 percentage share no meaningful range. ggplot2
# facets force a single shared scale, so the panels are built separately and
# combined with patchwork. Falls back to two files when patchwork is absent.
save_two_panel <- function(p_count, p_pct, stem, title = NULL, w = 13, h = 6) {
  if (HAVE_PATCHWORK) {
    combined <- p_count + p_pct + patchwork::plot_layout(ncol = 2)
    if (!is.null(title)) {
      combined <- combined + patchwork::plot_annotation(
        title = title,
        theme = theme(plot.title = element_text(size = 13, face = "bold"))
      )
    }
    show_plot(combined, stem, w = w, h = h)
  } else {
    message("  [NOTE] patchwork not installed; writing count and percent panels separately for ", stem)
    show_plot(p_count, paste0(stem, "_count"), w = w / 2 + 1, h = h)
    show_plot(p_pct,   paste0(stem, "_pct"),   w = w / 2 + 1, h = h)
  }
}

# Shared heatmap skeleton so the two panels are visually consistent.
heatmap_panel <- function(data, x_col, y_col, fill_col, label_expr,
                          fill_name, panel_title, x_lab, y_lab,
                          fill_limits = NULL) {
  ggplot(data, aes(factor(.data[[x_col]]), factor(.data[[y_col]]), fill = .data[[fill_col]])) +
    geom_tile(colour = "white") +
    geom_text(aes(label = label_expr), size = 2.7) +
    scale_fill_distiller(palette = "RdYlBu", direction = 1,
                         name = fill_name, limits = fill_limits) +
    labs(title = panel_title, x = x_lab, y = y_lab) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title  = element_text(size = 10, face = "bold"),
          legend.position = "right")
}

# Combination-effect plot (single panel): x = retention, y = one threshold,
# colour = the other threshold. Optional per-sample grey lines under the mean.
make_combo_plot <- function(data, level, x_col, y_col, colour_col, colour_label,
                            x_lab, y_lab, fname_stem, sample_data = NULL) {
  p <- ggplot(data, aes(.data[[x_col]], .data[[y_col]],
                       colour = factor(.data[[colour_col]])))
  if (!is.null(sample_data)) {
    p <- p + geom_line(
      data = sample_data,
      aes(x = .data[[x_col]], y = .data[[y_col]],
          group = interaction(sample_id, .data[[colour_col]])),
      colour = "grey70", alpha = 0.25, linewidth = 0.25, inherit.aes = FALSE
    )
  }
  p <- p +
    geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
    scale_y_continuous(labels = label_percent()) +
    scale_colour_viridis_d(name = colour_label) +
    labs(title    = paste0(toupper(substr(level,1,1)), substr(level,2,nchar(level)),
                           " - combined effect of both thresholds"),
         subtitle = if (!is.null(sample_data))
           "Bold = mean across samples; grey = each individual sample"
         else
           "Each line holds the other threshold fixed at its labelled value",
         x = x_lab, y = y_lab) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
  p <- p + if (x_col == "mean_pct") {
    scale_x_continuous(labels = label_comma(suffix = "%"), limits = c(0, 100))
  } else {
    scale_x_continuous(labels = label_comma())
  }
  show_plot(p, fname_stem)
}

# ══════════════════════════════════════════════════════════════════════════════
# Reads retained (single panel; mean line + per-sample grey lines)
# ══════════════════════════════════════════════════════════════════════════════
if (exists("results")) {
  message("\nBuilding reads-retained figures from combined_results.csv...")

  # Keep one read-level axis varying; hold the others at their loosest setting so
  # each heatmap tile and each line corresponds to exactly one threshold pair.
  results <- hold_non_primary(results, primary_col, "reads")

  plot_reads_retained <- function(fixed_col, fixed_val, x_col, x_lab, fname_stem) {
    sub <- results[get(fixed_col) == fixed_val]
    value_cols <- c("pct_fragment_hits_in_passing_genes_of_baseline",
                    "n_fragment_hits_in_passing_genes")
    sub_long <- melt(sub, id.vars = c("sample_id", x_col), measure.vars = value_cols)
    sub_long[, metric := factor(variable, levels = value_cols,
                                labels = c("Mean % of baseline retained", "Mean absolute reads retained"))]

    mean_data <- sub_long[, .(mean_val = mean(value)), by = c(x_col, "metric")]

    p <- ggplot(mean_data, aes(.data[[x_col]] * 100, mean_val)) +
      geom_line(data = sub_long, aes(x = .data[[x_col]] * 100, y = value, group = sample_id),
               colour = "grey55", alpha = 0.35, linewidth = 0.3, inherit.aes = FALSE) +
      geom_line(colour = "steelblue4", linewidth = 1) +
      geom_point(colour = "steelblue4", size = 1.6) +
      facet_wrap(~metric, scales = "free_y", ncol = 1) +
      scale_x_continuous(breaks = seq(0, 100, 20), labels = label_percent(scale = 1)) +
      scale_y_continuous(labels = label_comma()) +
      labs(title = paste0("Fragment-hits retained vs ", x_lab, " - every sample (grey) + mean (blue)"),
           x = x_lab, y = NULL) +
      theme_bw(base_size = 10) +
      theme(strip.background = element_rect(fill = "grey92"))
    show_plot(p, fname_stem, h = 9)
  }

  plot_reads_retained("min_gene_fraction", min(GF_SWEEP),
                      primary_col, primary_label, "reads_retained_vs_primary")
  plot_reads_retained(primary_col, min(primary_values),
                      "min_gene_fraction", "Min gene fraction", "reads_retained_vs_gf")

  reads_grid <- results[, .(
    mean_pct = mean(pct_fragment_hits_in_passing_genes_of_baseline),
    mean_n   = mean(n_fragment_hits_in_passing_genes)
  ), by = THRESHOLD_COLS]

  fwrite(reads_grid, file.path(out_dir, "reads_retained_combined_grid.csv"))

  reads_per_sample <- copy(results)
  setnames(reads_per_sample,
          c("pct_fragment_hits_in_passing_genes_of_baseline", "n_fragment_hits_in_passing_genes"),
          c("mean_pct", "mean_n"))

  # ── TWO-PANEL heatmap: absolute count | percent of baseline ────────────────
  p_reads_heat_count <- heatmap_panel(
    reads_grid, primary_col, "min_gene_fraction", "mean_n",
    label_expr  = round(reads_grid$mean_n, 0),
    fill_name   = "Fragment\nhits",
    panel_title = "Absolute fragment-hits retained",
    x_lab = primary_label, y_lab = "Min gene fraction"
  )
  p_reads_heat_pct <- heatmap_panel(
    reads_grid, primary_col, "min_gene_fraction", "mean_pct",
    label_expr  = paste0(round(reads_grid$mean_pct, 0), "%"),
    fill_name   = "% of\nbaseline",
    panel_title = "Percent of baseline retained",
    x_lab = primary_label, y_lab = "Min gene fraction",
    fill_limits = c(0, 100)
  )
  save_two_panel(p_reads_heat_count, p_reads_heat_pct,
                 "reads_retained_combined_heatmap",
                 title = "Fragment-hits retained across the joint threshold grid")

  make_combo_plot(reads_grid, "reads", "mean_pct", primary_col, "min_gene_fraction",
                  "Gene\nfraction", "Mean % of baseline retained", primary_label,
                  "reads_retained_combo_by_primary_pct", sample_data = reads_per_sample)
  make_combo_plot(reads_grid, "reads", "mean_n", primary_col, "min_gene_fraction",
                  "Gene\nfraction", "Mean absolute reads retained", primary_label,
                  "reads_retained_combo_by_primary_count", sample_data = reads_per_sample)
  make_combo_plot(reads_grid, "reads", "mean_pct", "min_gene_fraction", primary_col,
                  gsub("\n", " ", primary_label), "Mean % of baseline retained", "Min gene fraction",
                  "reads_retained_combo_by_gf_pct", sample_data = reads_per_sample)
  make_combo_plot(reads_grid, "reads", "mean_n", "min_gene_fraction", primary_col,
                  gsub("\n", " ", primary_label), "Mean absolute reads retained", "Min gene fraction",
                  "reads_retained_combo_by_gf_count", sample_data = reads_per_sample)

  baseline_per_sample <- results[
    get(primary_col) == min(primary_values) & min_gene_fraction == min(GF_SWEEP),
    .(sample_id, baseline_hits = n_fragment_hits_in_passing_genes)
  ]
  results_lost <- merge(results, baseline_per_sample, by = "sample_id")
  results_lost[, n_hits_lost   := baseline_hits - n_fragment_hits_in_passing_genes]
  results_lost[, pct_hits_lost := as.double(100 - pct_fragment_hits_in_passing_genes_of_baseline)]
  results_lost[, n_hits_lost   := as.double(n_hits_lost)]

  plot_reads_lost <- function(fixed_col, fixed_val, x_col, x_lab, fname_stem) {
    sub <- results_lost[get(fixed_col) == fixed_val]
    value_cols <- c("pct_hits_lost", "n_hits_lost")
    sub_long <- melt(sub, id.vars = c("sample_id", x_col), measure.vars = value_cols)
    sub_long[, metric := factor(variable, levels = value_cols,
                                labels = c("Mean % of baseline LOST", "Mean absolute reads LOST"))]

    mean_data <- sub_long[, .(mean_val = mean(value)), by = c(x_col, "metric")]

    p <- ggplot(mean_data, aes(.data[[x_col]] * 100, mean_val)) +
      geom_line(data = sub_long, aes(x = .data[[x_col]] * 100, y = value, group = sample_id),
               colour = "grey55", alpha = 0.35, linewidth = 0.3, inherit.aes = FALSE) +
      geom_line(colour = "firebrick", linewidth = 1) +
      geom_point(colour = "firebrick", size = 1.6) +
      facet_wrap(~metric, scales = "free_y", ncol = 1) +
      scale_x_continuous(breaks = seq(0, 100, 20), labels = label_percent(scale = 1)) +
      scale_y_continuous(labels = label_comma()) +
      labs(title = paste0("Fragment-hits LOST vs ", x_lab, " - every sample (grey) + mean (red)"),
           subtitle = "Complement of the retained figures; baseline = count at no-filter (qcov=0, gf=0)",
           x = x_lab, y = NULL) +
      theme_bw(base_size = 10) +
      theme(strip.background = element_rect(fill = "grey92"))
    show_plot(p, fname_stem, h = 9)
  }

  plot_reads_lost("min_gene_fraction", min(GF_SWEEP),
                  primary_col, primary_label, "reads_lost_vs_primary")
  plot_reads_lost(primary_col, min(primary_values),
                  "min_gene_fraction", "Min gene fraction", "reads_lost_vs_gf")
}

# ══════════════════════════════════════════════════════════════════════════════
# Percent identity + alignment-length distributions (single panel each)
# ══════════════════════════════════════════════════════════════════════════════
quant_path <- file.path(in_dir, "combined_length_quantiles.csv")
if (file.exists(quant_path)) {
  message("\nBuilding identity and alignment-length figures from combined_length_quantiles.csv...")
  quant <- fread(quant_path, showProgress = FALSE)
  quant <- quant[sample_id != "sample_id"]

  id_cols_present <- grep("^pct_identity_p", names(quant), value = TRUE)
  if (length(id_cols_present) == 0) {
    message("  [NOTE] No pct_identity_p* columns found in combined_length_quantiles.csv.")
  } else {
    quant[, (id_cols_present) := lapply(.SD, as.double), .SDcols = id_cols_present]
    quant[, min_query_coverage := as.double(min_query_coverage)]

    quant_mq_vals  <- if ("min_match_qcov" %in% names(quant)) sort(unique(as.double(quant$min_match_qcov))) else c(0)
    quant_id_vals  <- if ("min_identity"   %in% names(quant)) sort(unique(as.double(quant$min_identity)))   else c(0)
    quant_primary_col <- if (length(quant_mq_vals) > 1) "min_match_qcov" else
                         if (length(quant_id_vals) > 1) "min_identity" else
                         "min_query_coverage"
    quant_primary_label <- switch(quant_primary_col,
      min_match_qcov     = "Min match-qcov",
      min_identity        = "Min identity",
      min_query_coverage  = "Min query coverage")
    quant[, quant_primary := as.double(get(quant_primary_col))]

    # ── Figure 1: identity distribution vs the primary threshold ───────────
    id_summary <- quant[, .(
      p10  = mean(as.double(pct_identity_p10)),
      p25  = mean(as.double(pct_identity_p25)),
      p50  = mean(as.double(pct_identity_p50)),
      p75  = mean(as.double(pct_identity_p75)),
      p90  = mean(as.double(pct_identity_p90)),
      p5   = mean(as.double(pct_identity_p5)),
      p95  = mean(as.double(pct_identity_p95))
    ), by = .(quant_primary)]

    sample_median <- quant[, .(median_identity = mean(as.double(pct_identity_p50))),
                          by = .(sample_id, quant_primary)]

    p_id_dist <- ggplot(id_summary, aes(x = factor(round(quant_primary, 0)))) +
      geom_errorbar(aes(ymin = p5, ymax = p95), width = 0.2, colour = "grey40") +
      geom_crossbar(aes(y = p50, ymin = p25, ymax = p75),
                    fill = "steelblue4", alpha = 0.7, colour = "steelblue4", width = 0.5) +
      geom_line(data = sample_median,
               aes(x = factor(round(quant_primary, 0)),
                   y = median_identity, group = sample_id),
               colour = "grey70", alpha = 0.3, linewidth = 0.25, inherit.aes = FALSE) +
      geom_point(aes(y = p50), colour = "white", size = 2) +
      scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),
                         labels = label_percent(scale = 1)) +
      labs(title = paste0("Percent identity distribution vs ", quant_primary_label),
           subtitle = paste0("Box = IQR (p25-p75), whiskers = p5-p95, white dot = mean median,",
                             " grey lines = per-sample median.\n",
                             "Identity = (aln_length - NM) / aln_length; based on NM edit-distance tag."),
           x = quant_primary_label, y = "Percent identity across aligned region") +
      theme_bw(base_size = 10)
    show_plot(p_id_dist, "identity_distribution_vs_primary", h = 7, w = 10)

    # ── Figure 2: alignment length, WITH and WITHOUT accounting for matches ─
    ml_cols_present <- grep("^matched_length_p", names(quant), value = TRUE)
    if (length(ml_cols_present) == 0) {
      message("  [NOTE] No matched_length_p* columns found - skipping length comparison.")
    } else {
      quant[, (ml_cols_present) := lapply(.SD, as.double), .SDcols = ml_cols_present]
      al_cols_present <- grep("^aln_length_p", names(quant), value = TRUE)
      quant[, (al_cols_present) := lapply(.SD, as.double), .SDcols = al_cols_present]

      len_summary <- quant[, .(
        raw_p25 = mean(aln_length_p25), raw_p50 = mean(aln_length_p50), raw_p75 = mean(aln_length_p75),
        raw_p5  = mean(aln_length_p5),  raw_p95 = mean(aln_length_p95),
        match_p25 = mean(matched_length_p25), match_p50 = mean(matched_length_p50), match_p75 = mean(matched_length_p75),
        match_p5  = mean(matched_length_p5),  match_p95 = mean(matched_length_p95)
      ), by = .(quant_primary)]

      len_long <- rbindlist(list(
        len_summary[, .(quant_primary, length_type = "Raw aligned length (incl. mismatches/indels)",
                        p5 = raw_p5, p25 = raw_p25, p50 = raw_p50, p75 = raw_p75, p95 = raw_p95)],
        len_summary[, .(quant_primary, length_type = "Matched-only length (true matches)",
                        p5 = match_p5, p25 = match_p25, p50 = match_p50, p75 = match_p75, p95 = match_p95)]
      ))

      p_len_dist <- ggplot(len_long, aes(x = factor(round(quant_primary, 0)),
                                         colour = length_type, fill = length_type)) +
        geom_errorbar(aes(ymin = p5, ymax = p95), width = 0.2,
                     position = position_dodge(width = 0.6)) +
        geom_crossbar(aes(y = p50, ymin = p25, ymax = p75), alpha = 0.6, width = 0.5,
                     position = position_dodge(width = 0.6)) +
        scale_colour_manual(values = c("Raw aligned length (incl. mismatches/indels)" = "grey40",
                                       "Matched-only length (true matches)" = "steelblue4")) +
        scale_fill_manual(values = c("Raw aligned length (incl. mismatches/indels)" = "grey70",
                                     "Matched-only length (true matches)" = "steelblue4")) +
        scale_y_continuous(labels = label_comma()) +
        labs(title = paste0("Alignment length vs ", quant_primary_label, " - with and without accounting for matches"),
             subtitle = paste0("Box = IQR (p25-p75), whiskers = p5-p95, across samples.\n",
                               "Matched-only length = aln_length - NM (mismatches/indels subtracted via the NM edit-distance tag)."),
             x = quant_primary_label, y = "Length (bp)", colour = NULL, fill = NULL) +
        theme_bw(base_size = 10) + theme(legend.position = "bottom")
      show_plot(p_len_dist, "alignment_length_with_vs_without_matches", h = 7, w = 10)

      fwrite(len_summary, file.path(out_dir, "alignment_length_quantile_summary.csv"))
      message("  saved: alignment_length_quantile_summary.csv")
    }

    # ── Figure 3: identity quantile heatmap ─────────────────────────────────
    quant_long <- melt(id_summary,
                       id.vars     = c("quant_primary"),
                       measure.vars = c("p5", "p10", "p25", "p50", "p75", "p90", "p95"),
                       variable.name = "quantile", value.name = "mean_identity")

    p_id_heat <- ggplot(quant_long,
                        aes(factor(quantile, levels = c("p5","p10","p25","p50","p75","p90","p95")),
                            factor(round(quant_primary, 0)),
                            fill = mean_identity)) +
      geom_tile(colour = "white") +
      geom_text(aes(label = round(mean_identity, 1)), size = 2.6) +
      scale_fill_distiller(palette = "RdYlBu", direction = 1,
                           limits = c(50, 100), name = "Mean\nidentity (%)") +
      labs(title = paste0("Mean percent identity quantiles by ", quant_primary_label),
           subtitle = "Each cell = mean across samples of that quantile of the identity distribution",
           x = "Identity quantile", y = quant_primary_label) +
      theme_bw(base_size = 9)
    show_plot(p_id_heat, "identity_quantile_heatmap", w = 10)

    fwrite(id_summary, file.path(out_dir, "identity_quantile_summary.csv"))
    message("  saved: identity_quantile_summary.csv")

    # ── Figures 4-6: BASELINE (zero-filter) distribution histograms ────────────
    QLEVELS_HIST  <- c(5, 10, 25, 50, 75, 90, 95)
    INTERVAL_MASS <- c(5, 15, 25, 25, 15,  5)
    TRUE_MASS     <- 0.90
    N_PTS         <- 2000

    reconstruct <- function(qvals, n) {
      qvals <- as.numeric(qvals)
      if (any(is.na(qvals)) || length(qvals) != 7) return(numeric(0))
      pts <- list()
      for (i in 1:6) {
        ni <- round(n * INTERVAL_MASS[i] / 90)
        lo <- qvals[i]; hi <- qvals[i + 1]
        if (hi < lo) hi <- lo
        pts[[i]] <- if (ni > 0) runif(ni, lo, hi) else numeric(0)
      }
      unlist(pts)
    }

    build_pts <- function(data, prefix) {
      qcols <- paste0(prefix, "_p", QLEVELS_HIST)
      if (!all(qcols %in% names(data))) return(NULL)
      rbindlist(lapply(seq_len(nrow(data)), function(i) {
        pts <- reconstruct(as.numeric(data[i, ..qcols]), N_PTS)
        if (length(pts) == 0) return(NULL)
        data.table(value = pts)
      }))
    }

    density_hist_plot <- function(pts_dt, x_lab, title, fname, binwidth, xlim = NULL) {
      if (is.null(pts_dt) || nrow(pts_dt) == 0) {
        message("  [SKIP] No data to plot for ", fname); return(invisible(NULL))
      }
      v <- pts_dt$value
      lo_b <- floor(min(v) / binwidth) * binwidth
      hi_b <- ceiling(max(v) / binwidth) * binwidth + binwidth
      brks <- seq(lo_b, hi_b, by = binwidth)
      h <- hist(v, breaks = brks, plot = FALSE)
      hist_dt <- data.table(mid = h$mids, density = h$counts / length(v) * TRUE_MASS / binwidth)
      d <- density(v); dens_dt <- data.table(x = d$x, y = d$y * TRUE_MASS)

      p <- ggplot() +
        geom_col(data = hist_dt, aes(mid, density), width = binwidth,
                fill = "steelblue4", colour = "white", alpha = 0.8) +
        geom_line(data = dens_dt, aes(x, y), colour = "firebrick", linewidth = 0.7) +
        { if (!is.null(xlim)) coord_cartesian(xlim = xlim) } +
        labs(title = title,
             subtitle = paste0("Baseline (zero-filter) distribution reconstructed from per-sample quantiles ",
                               "(p5-p95; shaded area = ", TRUE_MASS,
                               "; outer 5% tails on each side not shown).\n",
                               uniqueN(pts_dt), " pooled sample-rows."),
             x = x_lab,
             y = paste0("Density (integrates to ~", TRUE_MASS,
                        " over shown range; remaining ", 1 - TRUE_MASS, " in unshown tails)")) +
        theme_bw(base_size = 10)
      show_plot(p, fname, h = 6, w = 9)
    }

    quant_thresh_cols <- intersect(c("min_query_coverage", "min_identity",
                                     "min_match_qcov"), names(quant))
    for (col in quant_thresh_cols) quant[, (col) := as.double(get(col))]
    quant_baseline <- quant[Reduce(`&`, lapply(quant_thresh_cols, function(c) get(c) == 0))]
    message("  Baseline (zero-filter) rows for histograms: ", nrow(quant_baseline))

    if (nrow(quant_baseline) > 0) {
      rl_pts <- build_pts(quant_baseline, "read_length")
      rl_range <- if (!is.null(rl_pts)) diff(range(rl_pts$value)) else 0
      rl_bw <- max(1, round(rl_range / 40))
      density_hist_plot(rl_pts, "Read length (bp)",
                        "Baseline read-length distribution",
                        "baseline_read_length_histogram", binwidth = rl_bw)

      density_hist_plot(build_pts(quant_baseline, "aln_length"),
                        "Aligned length (bp)  [clipping-only, includes mismatches/indels]",
                        "Baseline aligned-length distribution",
                        "baseline_aln_length_histogram", binwidth = rl_bw)

      ml_cols_hist <- paste0("matched_length_p", QLEVELS_HIST)
      if (all(ml_cols_hist %in% names(quant_baseline))) {
        density_hist_plot(build_pts(quant_baseline, "matched_length"),
                          "Matched length (bp)  [aln_length - NM; true matches only]",
                          "Baseline matched-length distribution (matches only)",
                          "baseline_matched_length_histogram", binwidth = rl_bw)
      }

      mq_exact <- paste0("match_qcov_pct_p", QLEVELS_HIST)
      if (all(mq_exact %in% names(quant_baseline))) {
        mq_pts <- build_pts(quant_baseline, "match_qcov_pct")
        mq_title_sfx <- ""
      } else if (all(c(ml_cols_hist, paste0("read_length_p", QLEVELS_HIST)) %in% names(quant_baseline))) {
        message("  [APPROX] match_qcov_pct_p* absent - approximating as matched_length_p{q}/read_length_p{q}")
        approx_mq <- copy(quant_baseline)
        mq_out_cols <- paste0("match_qcov_pct_p", QLEVELS_HIST)
        for (q in QLEVELS_HIST) {
          ml_c <- paste0("matched_length_p", q); rl_c <- paste0("read_length_p", q)
          approx_mq[, paste0("match_qcov_pct_p", q) := get(ml_c) / get(rl_c) * 100]
        }
        for (i in seq_len(nrow(approx_mq))) {
          v <- cummax(as.numeric(approx_mq[i, ..mq_out_cols]))
          approx_mq[i, (mq_out_cols) := as.list(v)]
        }
        mq_pts   <- build_pts(approx_mq, "match_qcov_pct")
        mq_title_sfx <- " [approximate]"
      } else {
        mq_pts <- NULL; mq_title_sfx <- ""
      }
      density_hist_plot(mq_pts,
                        "Match-query percentage (%)  [(matched_length / read_length) x 100]",
                        paste0("Baseline match-query-% distribution", mq_title_sfx),
                        "baseline_match_qcov_pct_histogram", binwidth = 2.5, xlim = c(0, 100))
    }
  }
} else {
  message("\n[NOTE] combined_length_quantiles.csv not found in ", in_dir, " - skipping identity/length figures.")
}

# ══════════════════════════════════════════════════════════════════════════════
# Per-level dropoff: derive counts at each (primary x gene-fraction) from raw
# gene_detail, for each annotation level. Heatmaps are two-panel: count | percent.
# ══════════════════════════════════════════════════════════════════════════════
message("\nBuilding per-level dropoff figures from combined_gene_detail.csv...")

# numeric coercions
num_cols <- intersect(c("min_query_coverage", "min_identity", "min_match_qcov",
                        "gene_fraction", "read_count"), names(genes))
for (col in num_cols) genes[, (col) := as.double(get(col))]

# the primary sweep column present in gene_detail
gd_primary_col <- if ("min_match_qcov" %in% names(genes) && length(match_qcov_values) > 1) "min_match_qcov" else
                  if ("min_identity"   %in% names(genes) && length(identity_values)   > 1) "min_identity"   else
                  "min_query_coverage"

level_available <- LEVELS_TO_PLOT[LEVELS_TO_PLOT %in% names(genes)]
if (length(level_available) == 0) {
  message("  [NOTE] none of ", paste(LEVELS_TO_PLOT, collapse = ", "),
          " are columns in gene_detail - skipping per-level dropoff.")
}

# Same treatment as the reads figures. Critical here: the per-level counts use
# uniqueN() over the filtered rows, so leaving a second read-level axis in place
# would union the detected features across all its values and report the count at
# its most permissive setting under a label implying otherwise.
genes <- hold_non_primary(genes, gd_primary_col, "per-level")

for (level in level_available) {
  message("  level: ", level)
  Level <- paste0(toupper(substr(level, 1, 1)), substr(level, 2, nchar(level)))

  # For each (sample, primary, gene_fraction cutoff), count distinct level entities
  # whose gene_fraction >= cutoff. Built across the data-derived GF_SWEEP grid.
  level_rows <- rbindlist(lapply(GF_SWEEP, function(gf) {
    passed <- genes[gene_fraction >= gf]
    cnt <- passed[, .(n_entities = uniqueN(get(level))),
                  by = c("sample_id", gd_primary_col)]
    cnt[, min_gene_fraction := gf]
    cnt
  }))
  if (nrow(level_rows) == 0) next

  # baseline per sample = count at (min primary, min gf)
  base <- level_rows[get(gd_primary_col) == min(primary_values) &
                       min_gene_fraction == min(GF_SWEEP),
                     .(sample_id, baseline_n = n_entities)]
  level_rows <- merge(level_rows, base, by = "sample_id", all.x = TRUE)
  level_rows[, pct_of_baseline := ifelse(baseline_n > 0, n_entities / baseline_n * 100, NA_real_)]

  # mean across samples
  level_mean <- level_rows[, .(mean_n = mean(n_entities),
                               mean_pct = mean(pct_of_baseline, na.rm = TRUE)),
                           by = c(gd_primary_col, "min_gene_fraction")]

  # ── Dropoff curve: x = primary, y = mean count, colour = gene fraction ──
  p1 <- ggplot(level_mean, aes(.data[[gd_primary_col]] * (if (PRIMARY_IS_IDENTITY) 1 else 100),
                               mean_n, colour = factor(min_gene_fraction))) +
    geom_line(data = level_rows,
             aes(x = .data[[gd_primary_col]] * (if (PRIMARY_IS_IDENTITY) 1 else 100),
                 y = n_entities,
                 group = interaction(sample_id, min_gene_fraction)),
             colour = "grey75", alpha = 0.2, linewidth = 0.25, inherit.aes = FALSE) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.5) +
    scale_colour_viridis_d(name = "Min gene\nfraction") +
    scale_x_continuous(labels = if (PRIMARY_IS_IDENTITY) label_comma() else label_percent(scale = 1)) +
    scale_y_continuous(labels = label_comma()) +
    labs(title = paste0(Level, " count vs ", primary_label),
         subtitle = paste(c("Bold = mean across samples; grey = each sample",
                            fixed_filter_note(genes, gd_primary_col)),
                          collapse = "\n"),
         x = primary_label, y = paste0("Distinct ", level, " count")) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
  show_plot(p1, paste0("dropoff_", level, "_count_vs_primary"))

  # ── % of baseline version ──
  p2 <- ggplot(level_mean, aes(.data[[gd_primary_col]] * (if (PRIMARY_IS_IDENTITY) 1 else 100),
                               mean_pct, colour = factor(min_gene_fraction))) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.5) +
    scale_colour_viridis_d(name = "Min gene\nfraction") +
    scale_x_continuous(labels = if (PRIMARY_IS_IDENTITY) label_comma() else label_percent(scale = 1)) +
    scale_y_continuous(labels = label_percent(scale = 1), limits = c(0, NA)) +
    labs(title = paste0(Level, " retained (% of baseline) vs ", primary_label),
         x = primary_label, y = "% of baseline retained") +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
  show_plot(p2, paste0("dropoff_", level, "_pct_vs_primary"))

  # ── TWO-PANEL heatmap: absolute count | percent of baseline ────────────────
  p3_count <- heatmap_panel(
    level_mean, gd_primary_col, "min_gene_fraction", "mean_n",
    label_expr  = round(level_mean$mean_n, 0),
    fill_name   = paste0("Mean\n", level, "\ncount"),
    panel_title = paste0("Absolute ", level, " count"),
    x_lab = primary_label, y_lab = "Min gene fraction"
  )
  p3_pct <- heatmap_panel(
    level_mean, gd_primary_col, "min_gene_fraction", "mean_pct",
    label_expr  = paste0(round(level_mean$mean_pct, 0), "%"),
    fill_name   = "% of\nbaseline",
    panel_title = paste0("Percent of baseline ", level, " retained"),
    x_lab = primary_label, y_lab = "Min gene fraction",
    fill_limits = c(0, 100)
  )
  save_two_panel(p3_count, p3_pct, paste0("dropoff_", level, "_heatmap"),
                 title = paste0(Level, " detected across the joint threshold grid"))

  fwrite(level_mean, file.path(out_dir, paste0("dropoff_", level, "_summary.csv")))
}

message("\nDone. Figures + summary CSVs written to: ", out_dir)