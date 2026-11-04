## max. code width ============================================================

## Title: Arc Institute in vivo pertubation screens
## Author: DKatzlberger
## Date: 21.01.2026

## Reproducibility ============================================================
seed <- 42
set.seed(seed)

## Settings ===================================================================
if (!interactive()) {
  # Commands
  full_args  <- commandArgs(trailingOnly = FALSE)
  trail_args <- commandArgs(trailingOnly = TRUE)

  # Script
  script_path <- sub("^--file=", "", grep("^--file=", full_args, value = TRUE)[1])
  script_name <- tools::file_path_sans_ext(basename(script_path))
  script_id   <- sub(pattern = "-.*", replacement = "", script_name)

  # Paths
  project_root <- trail_args[1]
  input_path   <- gsub(pattern = "\\\\ ", replacement = " ", trail_args[-1])
  output_path  <- "."

  # Params
  target_alpha <- 0.05

  # Output
  output_name <- paste0("_alpha_", target_alpha)

} else {
  # Script
  script_path <- rstudioapi::getActiveDocumentContext()$path
  script_name <- tools::file_path_sans_ext(basename(script_path))
  script_id   <- sub(pattern = "-.*", replacement = "", script_name)

  # Paths
  project_root <- "/home/scc_e_928652/ArcSimulationSynthetic"
  input_path   <- list.files(path = "results/nextflow/ArcSimSyn_offtarget_mixed", pattern = "_dge\\.rds$", recursive = TRUE, full.names = TRUE)
  output_path  <- "results/manual/ArcSimSyn_offtarget_mixed/15-visualize"

  #stopifnot(dirname(dirname(input_path)) == dirname(output_path))
  if (!dir.exists(output_path)) dir.create(path = output_path, recursive = TRUE)

  # Params
  target_alpha <- 0.05

  # Output
  output_name <- paste0("_alpha_", target_alpha)
}

## r-utils ====================================================================
source(file.path(project_root, "R", "utils.R"))

## Environment ================================================================
scc_info <- scc_status()

script_startup(
  root = project_root,
  renv = TRUE,
  bioc = TRUE,
  packages = c(
    "SummarizedExperiment",
    "SingleCellExperiment",
    "patchwork"
  )
)

## Script info ================================================================
showf(fmt = "\n%s\n", line80(paste0(script_name, ".R")), type = "cat")

## Dataset ====================================================================
readDGE <- function(file_paths) {

  # Validate input files
  file_paths <- as.character(file_paths)
  file_paths <- file_paths[file.exists(file_paths)]
  if (length(file_paths) == 0) showf("No DGE files found in the input.", type = "stop")
  
  # Load files
  showf("Loading %s DGE files.", length(file_paths), type = "message")
  loaded_data <- vector("list", length(file_paths))
  bad_files   <- character()

  for (i in seq_along(file_paths)) {
    current_path <- file_paths[i]

    loaded_data[[i]] <- tryCatch(
      expr = as.data.frame(readRDS(current_path)),
      error = function(e) {
        bad_files <<- c(bad_files, current_path)
        showf(
          "Failed to load %s: %s",
          current_path,
          conditionMessage(e),
          type = "message"
        )
        NULL
      }
    )

    if (i %% 1000 == 0) {
      showf(
        "Loaded %s/%s files.",
        i,
        length(file_paths),
        type = "message"
      )
    }
  }

  if (length(bad_files) > 0) {
    showf(
      "%s DGE files could not be loaded.",
      length(bad_files),
      type = "message"
    )
  }

  # Keep successfully loaded data frames
  valid_data <- Filter(Negate(is.null), loaded_data)

  if (length(valid_data) == 0) {
    stop("None of the DGE files could be loaded successfully.")
  }

  # Standardize columns
  showf("Standardizing columns.", type = "message")
  all_col_lists  <- lapply(valid_data, colnames)
  union_cols     <- unique(unlist(all_col_lists))
  intersect_cols <- Reduce(intersect, all_col_lists)

  non_overlapping_cols <- setdiff(union_cols, intersect_cols)

  if (length(non_overlapping_cols) > 0) {
    showf(
      "Removing non-overlapping columns: %s",
      paste(non_overlapping_cols, collapse = ", "),
      type = "message"
    )
  }

  valid_data <- lapply(
    valid_data,
    function(df) df[, intersect_cols, drop = FALSE]
  )

  # Combine data frames
  dge <- do.call(rbind, valid_data)
  rownames(dge) <- NULL
  dge
}

dge <- readDGE(input_path)

## Levels & order ============================================================
method_levels  <- c("limma", "edger", "deseq2", "dream", "nebula", "wilcox", "sceptre", "scmageck")
formula_levels <- c("~genotype", "~genotype + sample", "~genotype + (1|sample)", "~genotype + sample + (1|guide)", "~genotype + (1|sample) + (1|guide)")

dge$seed    <- factor(dge$seed)
dge$method  <- factor(dge$method, levels = method_levels)
dge$formula <- factor(dge$formula, levels = formula_levels)

## Stats ======================================================================
group_vars <- c("sim_mode", "workflow", "method_type",  "method" "gene_qc", "formula", "time_min", "knockout", "celltype",  "type", "seed")

# Confusion matrix
stats <- aggregate(
  cbind(
    P        = as.integer(!is.na(dge$FDR) & dge$FDR < target_alpha),
    n_tested = as.integer(!is.na(dge$FDR))
  ),
  by = dge[group_vars],
  FUN = sum
)

stats$positive_rate <- stats$P / stats$n_tested


## Panel A ====================================================================
## Panel A00: Annotations =====================================================

# Sim. mode
p_sim_mode <- ggplot2::ggplot(
  data = stats,
  mapping = ggplot2::aes(
    x = method,
    y = "Sim. mode",
    fill = sim_mode
  )
) +
ggplot2::geom_tile(
  color = "black",
  linewidth = 0.1,
  show.legend = TRUE
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::scale_fill_manual(
  name = "Sim. mode", 
  values = c(
    "synthetic"   = RColorBrewer::brewer.pal(8, "Dark2")[3], 
    "permutation" = RColorBrewer::brewer.pal(8, "Dark2")[4]
  )
) +
ggplot2::labs(
  x = NULL,
  y = NULL
) +
theme_DK(
  show_facet_x = TRUE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_borders = FALSE,
  show_grid = FALSE,
  show_axis = FALSE
) +
ggplot2::theme(
  strip.text.x = ggplot2::element_text(
    angle = 90,
    hjust = 0,
    vjust = 0.5
  ),
  axis.title.y = ggplot2::element_text(
    angle = 0,
    vjust = 0.5,
    hjust = 0.5
  ),
  plot.margin = ggplot2::margin(0, 0, 0, 0, "pt")
)

# Workflows
p_workflow <- ggplot2::ggplot(
  data = stats,
  mapping = ggplot2::aes(
    x = method,
    y = "Workflow",
    fill = workflow
  )
) +
ggplot2::geom_tile(
  color = "black",
  linewidth = 0.1
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::scale_fill_brewer(
  name = "Workflow",
  palette = "Set2"
) +
ggplot2::labs(
  x = "Method",
  y = NULL
) +
theme_DK(
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_borders = FALSE,
  show_grid = FALSE,
  show_axis = FALSE
) +
ggplot2::theme(
  axis.title.y = ggplot2::element_text(
    angle = 0,
    vjust = 0.5,
    hjust = 0.5
  ),
  plot.margin = ggplot2::margin(0, 0, 0, 0, "pt")
)

# Method type
p_method_type <- ggplot2::ggplot(
  data = stats,
  mapping = ggplot2::aes(
    x = method,
    y = "Method type",
    fill = method_type
  )
) +
ggplot2::geom_tile(
  color = "black",
  linewidth = 0.1
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::scale_fill_brewer(
  name = "Method type",
  palette = "Set1"
) +
ggplot2::labs(
  x = NULL,
  y = NULL
) +
theme_DK(
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_borders = FALSE,
  show_grid = FALSE,
  show_axis = FALSE
) +
ggplot2::theme(
  axis.title.y = ggplot2::element_text(
    angle = 0,
    vjust = 0.5,
    hjust = 0.5
  ),
  plot.margin = ggplot2::margin(0, 0, 0, 0, "pt")
)

# QC: gene filtering
p_gene_qc <- ggplot2::ggplot(
  data = stats,
  mapping = ggplot2::aes(
    x = method,
    y = "Gene filter",
    fill = gene_qc
  )
) +
ggplot2::geom_tile(
  color = "black",
  linewidth = 0.1
) +
ggplot2::scale_fill_manual(
  name = "Gene filter", 
  values = c(
    "TRUE"  = "#D4A017",
    "FALSE" = "darkred"
  ),
  labels = c(
    "TRUE"  = "method-specific",
    "FALSE" = "intersection"
  )
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = NULL,
  y = NULL
) +
theme_DK(
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_borders = FALSE,
  show_grid = FALSE,
  show_axis = FALSE
) +
ggplot2::theme(
  plot.margin = ggplot2::margin(0, 0, 0, 0)
)

# Patchwork
p_panel_A00 <- wrap_plots(
  p_sim_mode,
  p_workflow,
  p_method_type,
  p_gene_qc,
  ncol = 1
)

## Panel A01: Tested ==========================================================
p_panel_A01 <- ggplot2::ggplot(
  data = stats,
  mapping = ggplot2::aes(
    x = method,
    y = n_tested
  )
) +
ggplot2::stat_summary(
  fun = "mean",
  geom = "col",
  color = "black",
  fill = "grey80", 
  linewidth = 0.1,
  width = 0.8,
  position = ggplot2::position_dodge(
    width = 0.6
  )
) +
ggplot2::stat_summary(
  fun.data = "mean_sdl",
  fun.args = list(mult = 1),
  geom = "errorbar",
  color = "black",
  width = 0.2,
  linewidth = 0.2,
  position = ggplot2::position_dodge(
    width = 0.6
  )
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Method",
  y = "Nr. tests"
) +
theme_DK(
  rotate = 45,
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_grid = FALSE
)

## Panel A02: Positive rate ===================================================
p_panel_A02 <- ggplot2::ggplot(
  mapping = ggplot2::aes(
    x = method,
    y = positive_rate
  )
) +
ggplot2::stat_summary(
  data = stats[stats$type == "real", ],
  mapping = ggplot2::aes(
    fill = "Potential true positives (TP)"
  ),
  fun = "mean",
  geom = "col",
  color = "black",
  linewidth = 0.1,
  width = 0.6,
  position = ggplot2::position_dodge(
    width = 0.6
  )
) +
ggplot2::stat_summary(
  data = stats[stats$type == "perm", ],
  mapping = ggplot2::aes(
    fill = "False positives (FP)"
  ),
  fun = "mean",
  geom = "col",
  color = "black",
  linewidth = 0.1,
  width = 0.6,
  position = ggplot2::position_dodge(
    width = 0.6
  )
) +
ggplot2::stat_summary(
  data = stats[stats$type == "perm", ],
  mapping = ggplot2::aes(
  ),
  fun.data = "mean_sdl",
  fun.args = list(mult = 1),
  geom = "errorbar",
  color = "black",
  width = 0.2,
  linewidth = 0.2,
  position = ggplot2::position_dodge(
    width = 0.6
  )
) +
ggplot2::stat_summary(
  data = stats[stats$type == "real", ],
  mapping = ggplot2::aes(
    shape = knockout,
    color = celltype
  ),
  fun = "mean",
  geom = "point",
  size = 0.25,
  stroke = 0.2,
  position = ggplot2::position_dodge(
    width = 0.6
  )
) +
ggplot2::stat_summary(
  data = stats[stats$type == "real", ],
  mapping = ggplot2::aes(
    fill  = "Total positives (TP + FP)"
  ),
  fun = "mean",
  geom = "point",
  shape = 23,
  size = 1,
  color = "black",
  stroke = 0.1,
  position = ggplot2::position_dodge(
    width = 0.6
  )
) +
ggplot2::scale_fill_manual(
  name = "Discovery status",
  values = c(
    "Total positives (TP + FP)"     = "#008cdd",
    "Potential true positives (TP)" = "#ffd478",
    "False positives (FP)"          = "grey80"
  ),
  breaks = c(
    "Total positives (TP + FP)",
    "Potential true positives (TP)", 
    "False positives (FP)"
  )
) +
ggplot2::scale_color_manual(
  name = "Celltype",
  values = c(
    "Ery"         = "#762A83",
    "Imm. B-cell" = "#9970AB",
    "Gran. P"     = "#D01C8B",
    "GMP"         = "#E76289",
    "Gran."       = "#B2182B",
    "MEP"         = "#8C8D27",
    "HSC"         = "#4D9221",
    "MkP"         = "#1B7837"
  ),
  breaks = c(
    "Ery", "Imm. B-cell", "Gran. P", "GMP",
    "Gran.", "MEP", "HSC", "MkP"
  )
) +
ggplot2::scale_shape_manual(
  values = c(
    "Brd9" = 16,
    "Cebpa" = 17,
    "Chd4" = 15,
    "Kmt2d" = 18,
    "Rcor1" = 3,
    "Klf1" = 4
  )
) +
ggplot2::guides(
  color = ggplot2::guide_legend(
    ncol = 2,
    override.aes = list(
      size = 0.5, 
      stroke = 0.5
    )
  ),
  shape = ggplot2::guide_legend(
    ncol = 2,
    override.aes = list(
      size = 0.5,
      stroke = 0.5
    )
  ),
  fill = ggplot2::guide_legend(
    override.aes = list(
      size = c(1.5, NA, NA)
    )
  )
) +
ggplot2::scale_y_continuous(
  trans = scales::pseudo_log_trans(sigma = 0.01),
  limits = function(x) c(min(0, x), max(1, x)),
  breaks = scales::pretty_breaks(n = 3)
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Method", 
  y = "Positive rate",
  shape = "Knockout"
) +
theme_DK(
  rotate = 45,
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_grid = FALSE
)

## Panel A03: Best method =====================================================
p_panel_A03 <- ggplot2::ggplot(
  data = local({

    rates <- aggregate(
      positive_rate ~ method + formula + type,
      data = stats,
      FUN = mean,
      na.rm = TRUE
    )

    real <- rates[rates$type == "real", c("method", "formula", "positive_rate")]
    names(real)[3] <- "real_positive_rate"

    perm <- rates[rates$type == "perm", c("method", "formula", "positive_rate")]
    names(perm)[3] <- "perm_positive_rate"

    scores <- merge(
      real, perm,
      by = c("method", "formula"),
      all = TRUE
    )

    scores$ratio <- with(
      scores,
      ifelse(
        perm_positive_rate > 0,
        real_positive_rate / perm_positive_rate,
        ifelse(
          perm_positive_rate == 0 & real_positive_rate > 0,
          Inf,
          NA_real_
        )
      )
    )

    # Create one row per method + formula + metric
    metrics <- rbind(
      data.frame(
        scores[c("method", "formula")],
        metric = "Highest TP",
        value = scores$real_positive_rate - scores$perm_positive_rate
      ),
      data.frame(
        scores[c("method", "formula")],
        metric = "Lowest FP",
        value = scores$perm_positive_rate
      )
    )

    metrics$metric <- factor(
      metrics$metric,
      levels = c(
        "Highest TP",
        "Lowest FP"
      )
    )

    # Highlight ALL global winners, including exact ties
    metrics$best <- FALSE

    for (m in levels(metrics$metric)) {

      i <- which(metrics$metric == m & !is.na(metrics$value))

      if (length(i) == 0) next

      best_value <- if (m == "Lowest FP") {
        min(metrics$value[i])
      } else {
        max(metrics$value[i])
      }

      metrics$best[i] <- metrics$value[i] == best_value
    }

    metrics
  }),
  mapping = ggplot2::aes(
    x = method,
    y = metric
  )
) +
ggplot2::geom_tile(
  mapping = ggplot2::aes(
    fill = best
  ),
  color = "black",
  linewidth = 0.1
) +
ggplot2::scale_fill_manual(
  values = c(
    "FALSE" = "white",
    "TRUE"  = "pink"
  ),
  breaks = "TRUE",
  labels = "Best method",
  name = NULL
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Method",
  y = NULL
) +
theme_DK(
  rotate = 45,
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_borders = FALSE,
  show_axis = FALSE,
  show_grid = FALSE
) +
ggplot2::theme(
  plot.margin = ggplot2::margin(0, 0, 0, 0)
)

## Panel A04: Time (min) ======================================================
p_panel_A04 <- ggplot2::ggplot(
  data = aggregate(time_min ~ formula + method, data = stats, FUN = mean),
  mapping = ggplot2::aes(
    x = method,
    y = "Time (min)",
    fill = time_min
  )
) +
ggplot2::geom_tile(
  color = "black",
  linewidth = 0.1
) +
ggplot2::scale_fill_gradientn(
  name = "Time (min)",
  colours = c("#14213D", "#007F86", "#4CC9A6", "#F9C74F"),
  limits = c(0, NA),
  breaks = function(x) c(0, max(x)),
  labels = function(x) round(x, 1)
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Method",
  y = NULL
) +
theme_DK(
  rotate = 45,
  show_facet_x = FALSE,
  show_label_x = TRUE,
  show_ticks_x = TRUE,
  show_borders = FALSE,
  show_axis = FALSE,
  show_grid = FALSE
) +
ggplot2::theme(
  plot.margin = ggplot2::margin(0, 0, 0, 0)
)

## Patchwork & save ===========================================================
p_panel_A <-
  p_panel_A00 + 
  p_panel_A01 +
  p_panel_A02 +
  p_panel_A03 +
  p_panel_A04 +
  patchwork::plot_layout(
    ncol = 1,
    heights = c( 0.05, 0.05, 0.05, 0.05, 0.2, 0.4, 0.1, 0.05),
    guides = "collect"
  ) 

# Save
ggplot2::ggsave(
  p_panel_A,
  file = file.path(
    output_path,
    paste0(
      "panel_A",
      output_name,
      ".svg"
    )
  ),
  width = 10,
  height = 10.5,
  units = "cm"
)

saveRDS(
  p_panel_A,
  file = file.path(
    output_path,
    paste0(
      "panel_A",
      file_name,
      ".rds"
    )
  )
)


## Script info ================================================================
showf("\n◷ %s ====================================================\n", date(), type = "cat")
