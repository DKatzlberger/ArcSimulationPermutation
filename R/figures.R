## max. code width ============================================================

## Title: Arc Institute in vivo pertubation screens
## Author: DKatzlberger
## Date: 21.01.2026

## Reproducibility ============================================================
seed <- 42
set.seed(seed)

## Settings ===================================================================
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  project_root <- args[1]
  input_path   <- args[2]
  output_path  <- "."

} else {
  project_root <- "/home/scc_e_928652/ArcSimulationPermutation"
  input_path   <- "results/nextflow/ArcSimPerm_GSE213511"
  output_path  <- "results/nextflow/ArcSimPerm_GSE213511/figures"

  target_alpha <- 0.05

  file_name <- paste0("_alpha_", target_alpha)
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

## Dataset ====================================================================
readDGE <- function(
  path,
  pattern
){

  # Find DGE statistics as RDS
  file_paths <- list.files(
    path = path, 
    pattern = pattern, 
    full.names = TRUE, 
    recursive = TRUE
  )

  # Check for duplication
  any_duplicated_paths <- any(duplicated(file_paths))
  if (any_duplicated_paths) {
    n_dup_paths <- sum(duplicated(file_paths))
    showf("Note: There are %s duplicated DGE files.", n_dup_paths)
  }

  # Pre-allocate a list to store the data frames
  loaded_data <- vector("list", length(file_paths))
  bad_files <- c()

  # Load the files using a loop to track progress
  for (i in seq_along(file_paths)) {
    current_path <- file_paths[i]
    
    loaded_data[[i]] <- tryCatch({
      as.data.frame(readRDS(current_path))
    },
    error = function(e) {
      bad_files <<- c(bad_files, current_path)
      NULL
    })
    
    # Print progress every 500 files
    if (i %% 1000 == 0) {
      txt <- sprintf("%s/%s", i, length(file_paths))
      showf("Loaded %s", txt, type = "message")
    }
  }

  if (length(bad_files) > 0) {
    showf(
      "NOTE: Not all files could be loaded. %s failed file(s).",
      length(bad_files),
      type = "message"
    )
  }

  # Column standardization
  showf("Standardizing columns.", type = "message")
  
  valid_data <- loaded_data[!vapply(loaded_data, is.null, logical(1))]
  if (length(valid_data) > 0) {

    all_col_lists  <- lapply(valid_data, colnames)
    union_cols     <- unique(unlist(all_col_lists))
    intersect_cols <- Reduce(intersect, all_col_lists)
    
    non_overlapping_cols <- setdiff(union_cols, intersect_cols)
    
    if (length(non_overlapping_cols) > 0) {
      showf(
        "Removing '%s' non-overlapping column(s).", 
        paste(non_overlapping_cols, collapse = ", "), 
        type = "message"
      )
      
      loaded_data <- lapply(loaded_data, function(df) {
        if (is.null(df)) return(NULL)
        return(df[, intersect_cols, drop = FALSE])
      })
    } else {
      showf("All files have matching columns.", type = "message")
    }
  }

  # Combine into a dataframe
  out <- do.call(rbind, loaded_data)
  showf("Finished.", type = "message")
  out
}

dge <- readDGE(path = input_path, pattern = "~genotype.*_dge\\.rds$")

## Levels & order ============================================================
method_levels  <- c("limma", "edger", "deseq2", "dream", "nebula", "wilcox", "sceptre", "scmageck")
formula_levels <- c("~genotype", "~genotype + sample", "~genotype + (1|sample)", "~genotype + sample + (1|guide)", "~genotype + (1|sample) + (1|guide)")


## Panel Z ==================================================================
d_panel_Z <- dplyr::distinct(dge, formula, method, type, seed, knockout, celltype)

d_panel_Z$seed    <- factor(d_panel_Z$seed)
d_panel_Z$method  <- factor(d_panel_Z$method, levels = method_levels)
d_panel_Z$formula <- factor(d_panel_Z$formula, levels = formula_levels)

## Panel Z01: Nr. of seeds ==================================================
p_panel_Z01 <- ggplot2::ggplot(
  data = d_panel_Z,
  mapping = ggplot2::aes(
    x = method,
    group = interaction(celltype, knockout),
  )
) +
ggplot2::geom_bar(
  position = "dodge"
) +
ggplot2::scale_y_continuous(
  breaks = function(x) {c(0, max(x))},
  labels = scales::label_number(accuracy = 1)
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  rows = ggplot2::vars(type),
  scales = "free",
  space = "free"
) +
ggplot2::labs(
  x = "Method",
  y = "Number of seeds"
) +
theme_DK(
  rotate = 45,
  show_grid = FALSE,
  show_axis = FALSE
) +
ggplot2::theme(
  strip.text.x = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5),
  plot.margin = ggplot2::margin(0, 0, 0, 0)
)


## Panel A ====================================================================
group_vars <- c("workflow", "method_type",  "gene_qc", "method", "knockout", "celltype", "formula", "type", "seed", "time_min")

d_panel_A <- aggregate(
  cbind(
    P        = as.integer(!is.na(dge$FDR) & dge$FDR < target_alpha),
    n_tested = as.integer(!is.na(dge$FDR))
  ),
  by = dge[group_vars],
  FUN = sum
)

d_panel_A$positive_rate <- d_panel_A$P / d_panel_A$n_tested

# Factors
d_panel_A$seed    <- factor(d_panel_A$seed)
d_panel_A$method  <- factor(d_panel_A$method, levels = method_levels)
d_panel_A$formula <- factor(d_panel_A$formula, levels = formula_levels)

## Panel A00: Annotations =====================================================

# Workflows
p_workflow <- ggplot2::ggplot(
  data = d_panel_A,
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
  show_facet_x = TRUE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_borders = FALSE,
  show_grid = FALSE,
  show_axis = FALSE
) +
ggplot2::theme(
  strip.text.x = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5),
  axis.title.y = ggplot2::element_text(angle = 0, vjust = 0.5, hjust = 0.5),
  plot.margin = ggplot2::margin(0, 0, 0, 0)
)

# Method type
p_method_type <- ggplot2::ggplot(
  data = d_panel_A,
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
  plot.margin = ggplot2::margin(0, 0, 0, 0)
)

# QC: gene filtering
p_gene_qc <- ggplot2::ggplot(
  data = d_panel_A,
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
  p_workflow,
  p_method_type,
  p_gene_qc,
  ncol = 1
)

## Panel A01: Tested genes ====================================================
p_panel_A01 <- ggplot2::ggplot(
  data = d_panel_A,
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
  width = 0.6,
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
  y = "Tested genes"
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
  data = d_panel_A[d_panel_A$type == "real", ],
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
  data = d_panel_A[d_panel_A$type == "perm", ],
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
  data = d_panel_A[d_panel_A$type == "perm", ],
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
  data = d_panel_A[d_panel_A$type == "real", ],
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
  data = d_panel_A[d_panel_A$type == "real", ],
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
  # labels = scales::percent_format(
  #   accuracy = 1
  # )
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
      data = d_panel_A,
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
  data = aggregate(time_min ~ formula + method, data = d_panel_A, FUN = mean),
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
    heights = c( 0.05, 0.05, 0.05, 0.2, 0.4, 0.1, 0.05),
    guides = "collect"
  ) 

# Save
ggplot2::ggsave(
  p_panel_A,
  file = file.path(
    output_path,
    paste0(
      "panel_A",
      file_name,
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


## Panel B ====================================================================
## Panel B00: Annotations =====================================================

# Workflows
p_workflow <- ggplot2::ggplot(
  data = d_panel_A,
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
  show_facet_x = TRUE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_borders = FALSE,
  show_grid = FALSE,
  show_axis = FALSE
) +
ggplot2::theme(
  strip.text.x = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5),
  axis.title.y = ggplot2::element_text(angle = 0, vjust = 0.5, hjust = 0.5),
  plot.margin = ggplot2::margin(0, 0, 0, 0)
)

# Method type
p_method_type <- ggplot2::ggplot(
  data = d_panel_A,
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
  plot.margin = ggplot2::margin(0, 0, 0, 0)
)

# QC: gene filtering
p_gene_qc <- ggplot2::ggplot(
  data = d_panel_A,
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
p_panel_B00 <- wrap_plots(
  p_workflow,
  p_method_type,
  p_gene_qc,
  ncol = 1
)

## Panel B01: Tested genes ====================================================
p_panel_B01 <- ggplot2::ggplot(
  data = d_panel_A,
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
  width = 0.6,
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
  y = "Tested genes"
) +
theme_DK(
  rotate = 45,
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_grid = FALSE
)

## Panel B02: Nr. positives in perm. data =====================================
p_panel_B02 <- ggplot2::ggplot(
  data = d_panel_A[d_panel_A$type == "perm", ],
  mapping = ggplot2::aes(
    x = method,
    y = P
  )
) +
ggplot2::stat_summary(
  fun = "mean",
  geom = "col",
  color = "black",
  fill = "grey80",
  linewidth = 0.1,
  width = 0.6,
  position = ggplot2::position_dodge(
    width = 0.6
  )
) +
ggplot2::stat_summary(
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
ggplot2::scale_y_continuous(
  limits = function(x) c(min(0, x), max(1, x))
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Method",
  y = "Nr. positives\nin permutations"
) +
theme_DK(
  rotate = 45,
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_grid = FALSE
)

## Panel B03: Nr. perm. with at least one DEG =================================
p_panel_B03 <- ggplot2::ggplot(
  data = aggregate(
    P ~ method + formula + knockout + celltype,
    data = subset(d_panel_A, type == "perm"),
    FUN = function(x) sum(x > 0)
  ),
  mapping = ggplot2::aes(
    x = method,
    y = P
  )
) +
ggplot2::stat_summary(
  fun = "mean",
  geom = "col",
  color = "black",
  fill = "grey80",
  linewidth = 0.1,
  width = 0.6,
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
ggplot2::scale_y_continuous(
  limits = function(x) c(min(0, x), max(1, x))
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Method",
  y = "Nr. perm. with at\nleast one positive"
) +
theme_DK(
  rotate = 45,
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_grid = FALSE
)

## Panel B04: Nr. permutations ================================================
p_panel_B04 <- ggplot2::ggplot(
  data = aggregate(
    seed ~ method + formula,
    data = aggregate(
      seed ~ method + formula + knockout + celltype,
      data = unique(
        subset(d_panel_A, type == "perm")[
          , c("method", "formula", "knockout", "celltype", "seed")
        ]
      ),
      FUN = length
    ),
    FUN = mean
  ),
  mapping = ggplot2::aes(
    x = method,
    y = "Nr. permutations",
    fill = seed
  )
) +
ggplot2::geom_tile(
  color = "black",
  linewidth = 0.1
) +
ggplot2::scale_fill_gradient(
  name = "Nr. permutations",
  low = "white",
  high = "#8C6D31",
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
  show_borders = FALSE,
  show_grid = FALSE,
  show_axis = FALSE
) 

## Patchwork & save ===========================================================
p_panel_B <-
  p_panel_B00 + 
  p_panel_B01 +
  p_panel_B02 +
  p_panel_B03 +
  p_panel_B04 +
  patchwork::plot_layout(
    ncol = 1,
    heights = c( 0.05, 0.05, 0.05, 0.2, 0.4, 0.4, 0.05),
    guides = "collect"
  ) 

# Save
ggplot2::ggsave(
  p_panel_B,
  file = file.path(
    output_path,
    paste0(
      "panel_B",
      file_name,
      ".svg"
    )
  ),
  width = 10,
  height = 10.5,
  units = "cm"
)

saveRDS(
  p_panel_B,
  file = file.path(
    output_path,
    paste0(
      "panel_B",
      file_name,
      ".rds"
    )
  )
)



# d_panel_B <- dge[dge$type == "real" & !is.na(dge$pvalue) & !is.na(dge$logFC), ]
# d_panel_B <- d_panel_B[d_panel_B$celltype == "Ery" & d_panel_B$knockout == "Brd9", ]

# nrow(dge[is.na(dge$seed), ])
# head((d_panel_B[is.na(d_panel_B$method), ]))

# # Correaltion groups
# group_vars <- c("knockout", "celltype", "formula", "method")
# d_panel_B$group_id <- interaction(d_panel_B$knockout, d_panel_B$celltype, d_panel_B$formula, d_panel_B$method, sep = " | ", drop = TRUE, lex.order = TRUE)

# group_meta  <- unique(d_panel_B[c("group_id", group_vars)])
# group_meta  <- group_meta[with(group_meta, order(formula, celltype, knockout, method)), , drop = FALSE]
# group_order <- as.character(group_meta$group_id)

# wide <- reshape(
#   d_panel_B[c("gene", "group_id", "logFC")],
#   idvar = "gene",
#   timevar = "group_id",
#   direction = "wide"
# )

# mat <- as.matrix(wide[, match(group_order, colnames(wide)), drop = FALSE])
# colnames(mat) <- group_order
# rownames(mat) <- wide$gene

# cor_mat <- cor(mat, use = "pairwise.complete.obs", method = "pearson")
# diag(cor_mat) <- NA_real_

# cor_plot_data <- as.data.frame(
#   as.table(cor_mat),
#   stringsAsFactors = FALSE
# )

# names(cor_plot_data) <- c(
#   "group1", "group2", "correlation"
# )

# cor_plot_data$group1 <- as.character(cor_plot_data$group1)
# cor_plot_data$group2 <- as.character(cor_plot_data$group2)

# # Attach metadata without changing row order
# idx1 <- match(cor_plot_data$group1, group_meta$group_id)
# idx2 <- match(cor_plot_data$group2, group_meta$group_id)

# for (v in group_vars) {
#   cor_plot_data[[paste0(v, "1")]] <- group_meta[[v]][idx1]
#   cor_plot_data[[paste0(v, "2")]] <- group_meta[[v]][idx2]
# }

# cor_plot_data$group1 <- factor(
#   cor_plot_data$group1,
#   levels = group_order
# )

# cor_plot_data$group2 <- factor(
#   cor_plot_data$group2,
#   levels = rev(group_order)
# )

# ## Panel B00: Annotations =====================================================

# # Formula
# p_formula <- ggplot2::ggplot(
#   data = cor_plot_data,
#   mapping = ggplot2::aes(
#     x = group1,
#     y = "Formula",
#     fill = formula1
#   )
# ) +
# ggplot2::geom_tile(
#   color = NA,
#   linewidth = 0.1
# ) +
# ggplot2::facet_grid(
#   cols = ggplot2::vars(formula1),
#   scales = "free_x",
#   space = "free"
# ) +
# ggplot2::labs(
#   x = NULL,
#   y = NULL
# ) +
# theme_DK(
#   show_facet_x = TRUE,
#   show_label_x = FALSE,
#   show_ticks_x = FALSE,
#   show_borders = FALSE,
#   show_grid = FALSE,
#   show_axis = FALSE
# ) +
# ggplot2::theme(
#   axis.title.y = ggplot2::element_text(angle = 0, vjust = 0.5, hjust = 0.5),
#   plot.margin = ggplot2::margin(0, 0, 0, 0)
# )

# # Method
# p_method <- ggplot2::ggplot(
#   data = cor_plot_data,
#   mapping = ggplot2::aes(
#     x = group1,
#     y = "Method",
#     fill = method1
#   )
# ) +
# ggplot2::geom_tile(
#   color = NA,
#   linewidth = 0.1
# ) +
# ggplot2::facet_grid(
#   cols = ggplot2::vars(formula1),
#   scales = "free_x",
#   space = "free"
# ) +
# ggplot2::labs(
#   x = NULL,
#   y = NULL
# ) +
# theme_DK(
#   show_facet_x = TRUE,
#   show_label_x = FALSE,
#   show_ticks_x = FALSE,
#   show_borders = FALSE,
#   show_grid = FALSE,
#   show_axis = FALSE
# ) +
# ggplot2::theme(
#   strip.text.x = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5),
#   axis.title.y = ggplot2::element_text(angle = 0, vjust = 0.5, hjust = 0.5),
#   plot.margin = ggplot2::margin(0, 0, 0, 0)
# )

# # Patchwork
# p_fig_1_panel_B00 <- wrap_plots(
#   # p_formula,
#   p_method,
#   ncol = 1
# )

# ## Panel B01: Correlation =====================================================
# p_fig_1_panel_B01 <- ggplot2::ggplot(
#   data = cor_plot_data,
#   mapping = ggplot2::aes(
#     x = group1, 
#     y = group2, 
#     fill = correlation
#   )
# ) +
# ggplot2::geom_tile() +
# ggplot2::scale_fill_gradient2(
#   high = "#B2182B",
#   mid = "white",
#   low = "#2166AC",
#   midpoint = 0,
#   limits = c(-1, 1),
#   name = "Pearson R",
#   na.value = "grey85"
# ) +
# ggplot2::facet_grid(
#   cols = ggplot2::vars(formula1),
#   scales = "free_x",
#   space = "free"
# ) +
# ggplot2::labs(
#   x = NULL,
#   y = NULL
# ) +
# theme_DK(
#   rotate = 90,
#   show_borders = FALSE,
#   show_grid = FALSE,
#   show_axis = FALSE,
#   show_facet_x = FALSE,
#   show_label_y = FALSE,
#   show_ticks_y = FALSE,
#   show_ticks_x = FALSE
# )

# ## Panel B99: Annotations =====================================================
# # Knockout
# p_knockout <- ggplot2::ggplot(
#   data = cor_plot_data,
#   mapping = ggplot2::aes(
#     x = group1,
#     y = "Knockout",
#     fill = knockout1
#   )
# ) +
# ggplot2::geom_tile(
#   color = NA,
#   linewidth = 0.1
# ) +
# ggplot2::facet_grid(
#   cols = ggplot2::vars(formula1),
#   scales = "free_x",
#   space = "free"
# ) +
# ggplot2::labs(
#   x = NULL,
#   y = NULL
# ) +
# theme_DK(
#   show_facet_x = FALSE,
#   show_label_x = FALSE,
#   show_ticks_x = FALSE,
#   show_borders = FALSE,
#   show_grid = FALSE,
#   show_axis = FALSE
# ) +
# ggplot2::theme(
#   axis.title.y = ggplot2::element_text(angle = 0, vjust = 0.5, hjust = 0.5),
#   plot.margin = ggplot2::margin(0, 0, 0, 0)
# )

# # Celltype
# p_celltype <- ggplot2::ggplot(
#   data = cor_plot_data,
#   mapping = ggplot2::aes(
#     x = group1,
#     y = "Celltype",
#     fill = celltype1
#   )
# ) +
# ggplot2::geom_tile(
#   color = NA,
#   linewidth = 0.1
# ) +
# ggplot2::scale_fill_manual(
#   name = "Celltype",
#   values = c(
#     "Ery"         = "#762A83",
#     "Imm. B-cell" = "#9970AB",
#     "Gran. P"     = "#D01C8B",
#     "GMP"         = "#E76289",
#     "Gran."       = "#B2182B",
#     "MEP"         = "#8C8D27",
#     "HSC"         = "#4D9221",
#     "MkP"         = "#1B7837"
#   ),
#   breaks = c(
#     "Ery", "Imm. B-cell", "Gran. P", "GMP",
#     "Gran.", "MEP", "HSC", "MkP"
#   )
# ) +
# ggplot2::facet_grid(
#   cols = ggplot2::vars(formula1),
#   scales = "free_x",
#   space = "free"
# ) +
# ggplot2::labs(
#   x = NULL,
#   y = NULL
# ) +
# theme_DK(
#   rotate = 90,
#   show_facet_x = FALSE,
#   show_label_x = FALSE,
#   show_ticks_x = FALSE,
#   show_borders = FALSE,
#   show_grid = FALSE,
#   show_axis = FALSE
# ) +
# ggplot2::theme(
#   axis.title.y = ggplot2::element_text(angle = 0, vjust = 0.5, hjust = 0.5),
#   plot.margin = ggplot2::margin(0, 0, 0, 0)
# )

# # Patchwork
# p_fig_1_panel_B99 <- wrap_plots(
#   p_knockout,
#   p_celltype,
#   ncol = 1
# )

# ## Patchwork & save ===========================================================
# p_fig_1_panel_B <-
#   p_fig_1_panel_B00 + 
#   p_fig_1_panel_B01 +
#   p_fig_1_panel_B99 + 
#   patchwork::plot_layout(
#     ncol = 1,
#     heights = c(0.05, 1, 0.1),
#     guides = "collect"
#   ) 

# # Save
# ggplot2::ggsave(
#   p_fig_1_panel_B,
#   file = file.path(
#     output_path,
#     paste0(
#       "panel_B",
#       file_name,
#       ".svg"
#     )
#   ),
#   width = 20,
#   height = 11.5,
#   units = "cm"
# )

# saveRDS(
#   p_fig_1_panel_B,
#   file = file.path(
#     output_path,
#     paste0(
#       "panel_B",
#       file_name,
#       ".rds"
#     )
#   )
# )


