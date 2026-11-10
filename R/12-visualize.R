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
  project_root <- "/home/scc_e_928652/ArcSimulationPermutation"
  input_path   <- list.files(path = "results/nextflow/ArcSimPerm_GSE213511", pattern = "_dge\\.rds$", recursive = TRUE, full.names = TRUE)
  output_path  <- "results/manual/ArcSimPerm_GSE213511/12-visualize"

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

.readDGE <- fun_with_ram(readDGE)
dge <- .readDGE(input_path)

## Levels & order =============================================================
method_levels  <- c("limma", "edger", "deseq2", "dream", "nebula", "wilcox", "sceptre", "scmageck")
formula_levels <- c("~genotype", "~genotype + sample", "~genotype + (1|sample)", "~genotype + sample + (1|guide)", "~genotype + (1|sample) + (1|guide)")

dge$seed    <- factor(dge$seed)
dge$method  <- factor(dge$method, levels = method_levels)
dge$formula <- factor(dge$formula, levels = formula_levels)



## Data Panel A & B ===========================================================
pAB_vars_1 <- c("sim_mode", "workflow", "method_type",  "method", "gene_qc", "formula", "time_min")
pAB_vars_2 <- c("celltype", "knockout", "type", "seed")
pAB_group_vars <- c(pAB_vars_1, pAB_vars_2)

# Confusion matrix
pAB_stats <- stats::aggregate(
  cbind(
    P        = as.integer(!is.na(dge$FDR) & dge$FDR < target_alpha),
    n_tested = as.integer(!is.na(dge$FDR))
  ),
  by = dge[pAB_group_vars],
  FUN = sum
)

pAB_stats$positive_rate <- pAB_stats$P / pAB_stats$n_tested

## Data Panel C ===============================================================
pC_dat_perm <- dge[dge$type == "perm", ]
pC_n_perm   <- length(unique(pC_dat_perm$seed))

# Count significant detections per gene
pC_counts <- stats::aggregate(
  FDR ~ method + formula + celltype + knockout + gene,
  data = pC_dat_perm,
  FUN = \(x) sum(x < target_alpha, na.rm = TRUE)
)

names(pC_counts)[6] <- "n_perm"

# Calculate number of identified genes for each permutation count
pC_groups <- split(
  pC_counts,
  interaction(
    pC_counts[c("method", "formula", "celltype", "knockout")],
    drop = TRUE
  )
)

pC_cum_counts <- do.call(rbind, lapply(pC_groups, \(d) {
  data.frame(
    method = d$method[1],
    formula = d$formula[1],
    celltype = d$celltype[1],
    knockout = d$knockout[1],
    n_perm = seq_len(pC_n_perm),
    n_genes = vapply(
      seq_len(pC_n_perm),
      \(x) sum(d$n_perm >= x),
      numeric(1)
    )
  )
}))

rownames(pC_cum_counts) <- NULL

# Calculate mean and SD across knockout-celltype combinations
pC_groups <- split(
  pC_cum_counts,
  interaction(
    pC_cum_counts[c("formula", "method", "n_perm")],
    drop = TRUE
  )
)

pC_dat_summary <- do.call(rbind, lapply(pC_groups, \(d) {
  data.frame(
    formula = d$formula[1],
    method = d$method[1],
    n_perm = d$n_perm[1],
    n_seeds = pC_n_perm,
    pct = 100 * d$n_perm[1] / pC_n_perm,
    mean_n_genes = mean(d$n_genes, na.rm = TRUE),
    sd_n_genes = if (nrow(d) > 1) sd(d$n_genes, na.rm = TRUE) else 0
  )
}))

rownames(pC_dat_summary) <- NULL

# Bin permutation counts and summarise
pC_dat_summary$bin <- cut(
  pC_dat_summary$n_perm,
  breaks = seq(0, pC_n_perm, length.out = 11),
  include.lowest = TRUE,
  labels = FALSE
)

pC_groups <- split(
  pC_dat_summary,
  interaction(
    pC_dat_summary[c("formula", "method", "bin")],
    drop = TRUE
  )
)

pC_dat_binned <- do.call(rbind, lapply(pC_groups, \(d) {
  data.frame(
    formula = d$formula[1],
    method = d$method[1],
    bin = d$bin[1],
    n_perm = mean(d$n_perm),
    n_seeds = d$n_seeds[1],
    pct = mean(d$pct),
    mean_n_genes = mean(d$mean_n_genes, na.rm = TRUE),
    sd_n_genes = mean(d$sd_n_genes, na.rm = TRUE)
  )
}))

rownames(pC_dat_binned) <- NULL

# pC_perm <- dge[dge$type == "perm", ]

# # Count significant detections per gene across permutation datasets
# pC_gene_counts <- stats::aggregate(
#   FDR ~ method + formula + celltype + knockout + gene,
#   data = pC_perm,
#   FUN = \(x) sum(x < target_alpha, na.rm = TRUE)
# )

# names(pC_gene_counts)[6] <- "n_perm"

# # Total number of permutation datasets
# pC_n_seeds <- length(unique(pC_perm$seed))

# # Calculate number of identified genes for each permutation count
# pC_groups <- split(
#   pC_gene_counts,
#   interaction(
#     pC_gene_counts$method,
#     pC_gene_counts$formula,
#     pC_gene_counts$knockout,
#     pC_gene_counts$celltype,
#     drop = TRUE
#   )
# )

# pC_dat_perm <- do.call(rbind, lapply(pC_groups, \(d) {

#   if (max(d$n_perm) == 0) return(NULL)

#   data.frame(
#     method = d$method[1],
#     formula = d$formula[1],
#     knockout = d$knockout[1],
#     celltype = d$celltype[1],
#     n_seeds = pC_n_seeds,
#     n_perm = seq_len(pC_n_seeds),
#     n_genes = vapply(
#       seq_len(pC_n_seeds),
#       \(x) sum(d$n_perm >= x),
#       numeric(1)
#     )
#   )
# }))

# rownames(pC_dat_perm) <- NULL

# # Calculate percentage of permutation datasets
# pC_dat_perm$pct <- 100 * pC_dat_perm$n_perm / pC_dat_perm$n_seeds

# # Calculate mean and SD across knockout-celltype combinations
# pC_groups <- split(
#   pC_dat_perm,
#   interaction(
#     pC_dat_perm$formula,
#     pC_dat_perm$method,
#     pC_dat_perm$n_perm,
#     drop = TRUE
#   )
# )

# pC_dat_summary <- do.call(rbind, lapply(pC_groups, \(d) {
#   data.frame(
#     formula = d$formula[1],
#     method = d$method[1],
#     n_perm = d$n_perm[1],
#     n_seeds = d$n_seeds[1],
#     pct = d$pct[1],
#     mean_n_genes = mean(d$n_genes, na.rm = TRUE),
#     sd_n_genes = if (sum(!is.na(d$n_genes)) > 1) {
#       sd(d$n_genes, na.rm = TRUE)
#     } else {
#       0
#     }
#   )
# }))

# rownames(pC_dat_summary) <- NULL

# # Bin permutation counts (10 bins)
# pC_bin_breaks <- seq(
#   0,
#   max(pC_dat_summary$n_perm),
#   length.out = 11
# )

# pC_dat_summary$bin <- cut(
#   pC_dat_summary$n_perm,
#   breaks = pC_bin_breaks,
#   include.lowest = TRUE,
#   labels = FALSE
# )

# # Summarise each bin
# pC_groups <- split(
#   pC_dat_summary,
#   interaction(
#     pC_dat_summary$formula,
#     pC_dat_summary$method,
#     pC_dat_summary$bin,
#     drop = TRUE
#   )
# )

# pC_dat_binned <- do.call(rbind, lapply(pC_groups, \(d) {
#   data.frame(
#     formula = d$formula[1],
#     method = d$method[1],
#     bin = d$bin[1],
#     n_perm = mean(d$n_perm),
#     n_seeds = d$n_seeds[1],
#     pct = mean(d$pct),
#     mean_n_genes = mean(d$mean_n_genes, na.rm = TRUE),
#     sd_n_genes = mean(d$sd_n_genes, na.rm = TRUE)
#   )
# }))

# rownames(pC_dat_binned) <- NULL

## Data Panel D ===============================================================
pD_vars <- c("formula", "method", "celltype", "knockout")

# Separate original and permutation results
pD_real <- dge[dge$type == "real", ]
pD_perm <- dge[dge$type == "perm", ]

# Identify significant original DEGs and permutation detections
pD_real_sig <- pD_real[!is.na(pD_real$FDR) & pD_real$FDR < target_alpha, ]
pD_perm_sig <- pD_perm[!is.na(pD_perm$FDR) & pD_perm$FDR < target_alpha, ]

# Count permutations
pD_n_seeds <- stats::aggregate(
  seed ~ formula + method + celltype + knockout,
  data = pD_perm,
  FUN = \(x) length(unique(x))
)
names(pD_n_seeds)[5] <- "n_seeds"

# Count detection per gene across permutations
pD_n_detected <- stats::aggregate(
  seed ~ formula + method + celltype + knockout + gene,
  data = unique(pD_perm_sig[c(pD_vars, "gene", "seed")]),
  FUN = length
)
names(pD_n_detected)[6] <- "n_detected"

# Combine original DEGs with permutation detection counts
pD_dat <- merge(
  pD_real_sig,
  pD_n_detected,
  by = c(pD_vars, "gene"),
  all.x = TRUE,
  sort = FALSE
)

pD_dat <- merge(
  pD_dat,
  pD_n_seeds,
  by = pD_vars,
  all.x = TRUE,
  sort = FALSE
)

# Calculate absolute logFC and percentage of permutation detections
pD_dat$n_detected[is.na(pD_dat$n_detected)] <- 0
pD_dat$abs_logFC    <- abs(pD_dat$logFC)
pD_dat$pct_detected <- 100 * pD_dat$n_detected /pD_dat$n_seeds

# Rank genes by absolute logFC within each comparison
pD_groups <- split(pD_dat, interaction(pD_dat[pD_vars], drop = TRUE))

pD_dat <- do.call(rbind, lapply(pD_groups, \(d) {
  d <- d[order(-d$abs_logFC), ]
  d$rank <- seq_len(nrow(d))
  d
}))

rownames(pD_dat) <- NULL

# Correlation across celltypes and knockouts
pD_groups   <- split(seq_len(nrow(pD_dat)), interaction(pD_dat[c("formula", "method")], drop = TRUE))
pD_dat$corr <- NA_real_

for (i in pD_groups) {
  d <- pD_dat[i, ]
  d <- d[complete.cases(d[c("rank", "pct_detected")]), ]

  if (nrow(d) < 2 || sd(d$rank) == 0 ||sd(d$pct_detected) == 0) {
    next
  }

  pD_dat$corr[i] <- cor(
    d$rank,
    d$pct_detected,
    method = "spearman"
  )
}

# Add correlation annotation
pD_dat$anno <- paste0("Spearman = ", round(pD_dat$corr, 2))

## Data Panel E ===============================================================
pE_vars  <- c("formula", "method", "celltype", "knockout")

pE_dat    <- dge[dge$type == "real" & !is.na(dge$logFC), ]
pE_dat$id <- interaction(pE_dat[pE_vars], drop = TRUE, sep = "_")

pE_mat <- stats::xtabs(logFC ~ gene + id, data = pE_dat)
pE_mat <- as.matrix(pE_mat)

pE_cor <- stats::cor(pE_mat, method = "spearman", use = "pairwise.complete.obs")

## Data Panel F ===============================================================
pF_vars <- c("formula", "method", "celltype", "knockout")

pF_dat    <- dge[dge$type == "real" & !is.na(dge$logFC) & is.finite(dge$logFC), ]
pF_groups <- split(pF_dat, interaction(pF_dat[pF_vars], drop = TRUE, sep = "_"))

# Enrichmet DB
pF_hallmark <- msigdbr::msigdbr(db_species = "MM", species = "Mus musculus", collection = "MH")
TERM2GENE   <- unique(data.frame(term = pF_hallmark$gs_name, gene = pF_hallmark$gene_symbol))
pathways    <- split(TERM2GENE$gene, TERM2GENE$term)

pF_fgsea <- do.call(rbind, lapply(pF_groups, \(d) {

  # Require one value per gene and comparison
  if (anyDuplicated(d$gene)) {
    showf("Duplicate genes found in comparison %s", d$id[1], type = "stop")
  }

  # Rank all genes by signed logFC
  pF_stats <- stats::setNames(d$logFC, d$gene)
  pF_stats <- sort(pF_stats, decreasing = TRUE)

  # Run fgsea
  res <- fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = pF_stats,
    minSize = 15,
    maxSize = 500,
    eps = 0
  )

  res <- as.data.frame(res)

  if (nrow(res) == 0) return(NULL)

  # Add comparison metadata
  res$sim_mode    <- d$sim_mode[1]
  res$workflow    <- d$workflow[1]
  res$method_type <- d$method_type[1]
  res$gene_qc     <- d$gene_qc[1]

  res$formula  <- d$formula[1]
  res$method   <- d$method[1]
  res$celltype <- d$celltype[1]
  res$knockout <- d$knockout[1]

  res$database <- "MSigDB Hallmark (mouse)"

  res
}))

rownames(pF_fgsea) <- NULL


## Panel A ====================================================================
## Panel A00: Annotation ------------------------------------------------------

# Sim. mode
p_sim_mode <- ggplot2::ggplot(
  data = pAB_stats,
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
  data = pAB_stats,
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
ggplot2::scale_fill_manual(
  name = "Workflow",
  values = c(
    "singlecell" = "grey10",
    "pseudobulk" = "grey60"
  )
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
  data = pAB_stats,
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
  data = pAB_stats,
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

## Panel A01: Tested genes ----------------------------------------------------
p_panel_A01 <- ggplot2::ggplot(
  data = pAB_stats,
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
  y = "Nr. tests"
) +
theme_DK(
  rotate = 45,
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_grid = FALSE
)

## Panel A02: LogFC -----------------------------------------------------------
p_panel_A02 <- ggplot2::ggplot(
   data = rbind(
    transform(
      dge[dge$type == "real" & !is.na(dge$FDR), ],
      gene_set = factor("Tested genes", levels = c("Tested genes", "Positives"))
    ),
    transform(
      dge[dge$type == "real" & !is.na(dge$FDR) & dge$FDR < target_alpha, ],
      gene_set = factor("Positives", levels = c("Tested genes", "Positives"))
    )
  ),
  mapping = ggplot2::aes(
    x = method,
    y = logFC
  )
) +
ggplot2::geom_violin(
  mapping = ggplot2::aes(
    fill = gene_set
  ),
  color = "black",
  scale = "width",
  linewidth = 0.1,
  show.legend = FALSE
) +
ggplot2::scale_fill_manual(
  name = NULL,
  values = c(
    "Tested genes" = "white",
    "Positives"    = "#008cdd"
  )
) +
ggplot2::scale_y_continuous(
  breaks = scales::pretty_breaks(n = 3)
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free",
  space = "free"
) +
theme_DK(
  rotate = 45,
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_grid = FALSE
)

## Panel A03: Positive rate ---------------------------------------------------
p_panel_A03 <- ggplot2::ggplot(
  mapping = ggplot2::aes(
    x = method,
    y = positive_rate
  )
) +
ggplot2::stat_summary(
  data = pAB_stats[pAB_stats$type == "real", ],
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
  data = pAB_stats[pAB_stats$type == "perm", ],
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
  data = pAB_stats[pAB_stats$type == "perm", ],
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
  data = pAB_stats[pAB_stats$type == "real", ],
  mapping = ggplot2::aes(
    color = celltype,
    shape = knockout
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
  data = pAB_stats[pAB_stats$type == "real", ],
  mapping = ggplot2::aes(
    fill  = "Total positives (TP + FP)"
  ),
  fun = "mean",
  geom = "point",
  shape = 23,
  size = 0.8,
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
  values = colors_celltype
) +
ggplot2::scale_shape_manual(
  name = "Knockout",
  values = shapes_knockout
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
  # trans = scales::pseudo_log_trans(sigma = 0.01),
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
  y = "Positive rate"
) +
theme_DK(
  rotate = 45,
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_grid = FALSE
)

## Panel A04: Best method -----------------------------------------------------
p_panel_A04 <- ggplot2::ggplot(
  data = local({

    rates <- aggregate(
      positive_rate ~ method + formula + type,
      data = pAB_stats,
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

## Panel A05: Time (min) ------------------------------------------------------
p_panel_A05 <- ggplot2::ggplot(
  data = aggregate(time_min ~ formula + method, data = pAB_stats, FUN = mean),
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

## Patchwork & save -----------------------------------------------------------
panel_dir <- file.path(output_path, "panel_A")
if (!dir.exists(panel_dir)) dir.create(panel_dir, recursive = TRUE)

# Patchwork
p_panel_A <-
  p_panel_A00 + 
  p_panel_A01 +
  p_panel_A02 +
  p_panel_A03 +
  p_panel_A04 +
  p_panel_A05 +
  patchwork::plot_layout(
    ncol = 1,
    heights = c( 0.05, 0.05, 0.05, 0.05, 0.2, 0.4, 0.4, 0.1, 0.05),
    guides = "collect"
  ) 

# Save
purrr::walk(
  c("svg", "pdf", "png"),
  \(fmt) ggplot2::ggsave(
    plot = p_panel_A,
    file = file.path(
      panel_dir,
      paste0(
        "panel_A", 
        output_name, 
        ".", fmt
      )
    ),
    width = 10,
    height = 12,
    units = "cm",
    dpi = 300
  )
)

saveRDS(
  object = p_panel_A,
  file = file.path(
    panel_dir,
    paste0(
      "panel_A",
      output_name,
      ".rds"
    )
  )
)


## Panel B ====================================================================
## Panel A00: Annotation ------------------------------------------------------

# Sim. mode
p_sim_mode <- ggplot2::ggplot(
  data = pAB_stats,
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
  data = pAB_stats,
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
ggplot2::scale_fill_manual(
  name = "Workflow",
  values = c(
    "singlecell" = "grey10",
    "pseudobulk" = "grey60"
  )
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
  data = pAB_stats,
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
  data = pAB_stats,
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
  p_sim_mode,
  p_workflow,
  p_method_type,
  p_gene_qc,
  ncol = 1
)

## Panel B01: Tested genes ----------------------------------------------------
p_panel_B01 <- ggplot2::ggplot(
  data = pAB_stats,
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
  y = "Nr. tests"
) +
theme_DK(
  rotate = 45,
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_grid = FALSE
)

## Panel B02: Nr. positives in perm. data -------------------------------------
p_panel_B02 <- ggplot2::ggplot(
  data = pAB_stats[pAB_stats$type == "perm", ],
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

## Panel B03: Nr. perm. with at least one DEG ---------------------------------
p_panel_B03 <- ggplot2::ggplot(
  data = aggregate(
    P ~ method + formula + knockout + celltype,
    data = subset(pAB_stats, type == "perm"),
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
  limits = function(x) c(min(0, x), max(1, x)),
  breaks = scales::breaks_pretty()
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

## Panel B04: Nr. permutations ------------------------------------------------
p_panel_B04 <- ggplot2::ggplot(
  data = aggregate(
    seed ~ method + formula,
    data = aggregate(
      seed ~ method + formula + knockout + celltype,
      data = unique(
        subset(pAB_stats, type == "perm")[
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

## Patchwork & save -----------------------------------------------------------
panel_dir <- file.path(output_path, "panel_B")
if (!dir.exists(panel_dir)) dir.create(panel_dir, recursive = TRUE)

# Patchwork
p_panel_B <-
  p_panel_B00 + 
  p_panel_B01 +
  p_panel_B02 +
  p_panel_B03 +
  p_panel_B04 +
  patchwork::plot_layout(
    ncol = 1,
    heights = c( 0.05, 0.05, 0.05, 0.05, 0.2, 0.4, 0.4, 0.05),
    guides = "collect"
  ) 

# Save
purrr::walk(
  c("svg", "pdf", "png"),
  \(fmt) ggplot2::ggsave(
    plot = p_panel_B,
    file = file.path(
      panel_dir,
      paste0(
        "panel_B", 
        output_name, 
        ".", fmt
      )
    ),
    width = 10,
    height = 10.5,
    units = "cm",
    dpi = 300
  )
)

saveRDS(
  object = p_panel_B,
  file = file.path(
    panel_dir,
    paste0(
      "panel_B",
      output_name,
      ".rds"
    )
  )
)


## Panel C ====================================================================
## Panel C01 ------------------------------------------------------------------
p_list <- purrr::map(unique(pC_dat_binned$formula), \(f) {

  dat <- pC_dat_binned[pC_dat_binned$formula == f, ]

  ggplot2::ggplot(
    data = dat,
    mapping = ggplot2::aes(
      x = n_perm,
      y = mean_n_genes
    )
  ) +
  ggplot2::geom_col(
    fill = "grey80",
    color = "black",
    linewidth = 0.1
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = mean_n_genes - sd_n_genes,
      ymax = mean_n_genes + sd_n_genes
    ),
    width = 1,
    linewidth = 0.2
  ) +
  ggplot2::facet_wrap(
    ~method,
    scales = "free_y"
  ) +
  ggplot2::scale_x_continuous(
  breaks = \(x) {
      c(1, scales::breaks_pretty(n = 3)(x)[
        scales::breaks_pretty(n = 3)(x) > 1
      ])
    },
    labels = \(x) {
      paste0(
        round(x),
        "\n(",
        round(100 * x / dat$n_seeds[1], digits = 0),
        "%)"
      )
    }
  ) +
  ggplot2::scale_y_continuous(
    breaks = scales::breaks_pretty()
  ) +
  ggplot2::labs(
    title = f,
    x = "Nr. permutations where a gene was wrongly identified positive",
    y = "Nr. identified positives"
  ) +
  theme_DK(show_grid = FALSE)
})

## Patchwork & save -----------------------------------------------------------
panel <- "panel_C"

# Directory
panel_dir <- file.path(output_path, panel)
if (!dir.exists(panel_dir)) dir.create(panel_dir, recursive = TRUE)

# Formula and formats
purrr::walk2(
  p_list,
  unique(pC_dat_binned$formula),
  \(p, f) {
    
    # Formula dir
    formula_dir <- file.path(panel_dir, f)
    if (!dir.exists(formula_dir)) dir.create(formula_dir, recursive = TRUE)

    # Save
    purrr::walk(
      c("svg", "pdf", "png"),
      \(fmt) ggplot2::ggsave(
        plot = p,
        filename = file.path(
          formula_dir,
          paste0(panel, "_", f, output_name, ".", fmt)
        ),
        width = 8,
        height = 8,
        units = "cm",
        dpi = 300
      )
    )

    saveRDS(
      object = p,
      file = file.path(
        formula_dir,
        paste0(panel, "_", f, output_name, ".rds")
      )
    )
  }
)


## Panel D ====================================================================
## Panel D01 ------------------------------------------------------------------
p_list <- purrr::map(unique(pD_dat$formula), \(f) {

  dat <- pD_dat[pD_dat$formula == f, ]

  ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = rank,
      y = pct_detected
    )
  ) +
  ggplot2::geom_point(
    mapping = ggplot2::aes(
      color = celltype,
      shape = knockout
    ),
    size = 0.2
  ) +
  # ggplot2::geom_smooth(
  #   method = "loess",
  #   formula = y ~ x,
  #   se = TRUE,
  #   level = 0.95,
  #   color = "steelblue",
  #   fill = "grey80",
  #   linewidth = 0.5
  # ) +
  ggplot2::geom_text(
    data = unique(dat[c("method", "anno")]),
    mapping = ggplot2::aes(
      x = -Inf,
      y = 90,
      label = anno
    ),
    inherit.aes = FALSE,
    size = 1.7,
    hjust = -0.05,
    vjust = -0.5
  ) +
  ggplot2::scale_color_manual(
    name = "Celltype",
    values = colors_celltype
  ) +
  ggplot2::scale_shape_manual(
    name = "Knockout",
    values = shapes_knockout
  ) +
  ggplot2::facet_wrap(
    ~ method,
    scales = "free"
  ) +
  ggplot2::scale_x_continuous(
    breaks = scales::breaks_pretty(n = 3),
    labels = \(x) {
      vapply(x, \(z) {
        if (length(dat$rank) == 0 || all(is.na(dat$rank))) {
          return("")
        }

        idx <- which.min(abs(dat$rank - z))

        paste0(
          dat$rank[idx],
          "\n(",
          round(dat$abs_logFC[idx], 1),
          ")"
        )
      }, character(1))
    }
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 110),
    breaks = seq(0, 100, 25)
  ) +
  ggplot2::labs(
    title = f,
    x = "Gene rank (absolute log2FC) from original data",
    y = "% of permutions wrongly identifying genes as positives"
  ) +
  theme_DK(show_grid = FALSE)
})

## Patchwork & save -----------------------------------------------------------
panel <- "panel_D"

# Directory
panel_dir <- file.path(output_path, panel)
if (!dir.exists(panel_dir)) dir.create(panel_dir, recursive = TRUE)

# Formula and formats
purrr::walk2(
  p_list,
  unique(pD_dat$formula),
  \(p, f) {
    
    # Formula dir
    formula_dir <- file.path(panel_dir, f)
    if (!dir.exists(formula_dir)) dir.create(formula_dir, recursive = TRUE)

    # Save
    purrr::walk(
      c("svg", "pdf", "png"),
      \(fmt) ggplot2::ggsave(
        plot = p,
        filename = file.path(
          formula_dir,
          paste0(panel, "_", f, output_name, ".", fmt)
        ),
        width = 9,
        height = 8,
        units = "cm",
        dpi = 300
      )
    )

    saveRDS(
      object = p,
      file = file.path(
        formula_dir,
        paste0(panel, "_", f, output_name, ".rds")
      )
    )
  }
)


## Panel E ====================================================================
## Panel E01: Correlation heatmap ---------------------------------------------
pE_anno <- unique(pE_dat[c("id", "workflow", pF_vars)])
pE_anno <- pE_anno[match(colnames(pE_cor), pE_anno$id), ]
rownames(pE_anno) <- pE_anno$id

pE_groups <- split(pE_anno, interaction(pE_anno[c("celltype", "knockout")], drop = TRUE, sep = "_"))

p_list <- stats::setNames(
  Map(\(d, group) {

    pE_ids <- d$id
    pE_mat <- pE_cor[pE_ids, pE_ids, drop = FALSE]

    ha <- cmplx_anno(
      Formula = d$formula,
      Method = d$method,
      Workflow = d$workflow,
      col = list(
        Formula = colors_formula, 
        Method = colors_method,
        Workflow = c("singlecell" = "grey10", "pseudobulk" = "grey60")
      ),
      annotation_name_side = "left",
      annotation_legend_param = list(
        Formula = list(ncol = 1),
        Method = list(ncol = 3),
        Workflow = list(ncol = 2)
      ),
      which = "column"
    )

    cmplx_ht(
      pE_mat,
      name = "Spearman",
      col = circlize::colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
      heatmap_legend_param = list(at = c(-1, 0, 1), labels = c("-1", "0", "1")),
      top_annotation = ha,
      column_title = paste("Celltype:", d$celltype[1], "Knockout:", d$knockout[1]),
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_row_names = FALSE,
      show_column_names = FALSE,
      base_size = 5,
      cell_size = 2
    )
  }, pE_groups, names(pE_groups)),
  nm = names(pE_groups)
)

## Patchwork & save -----------------------------------------------------------
panel_dir <- file.path(output_path, "panel_E")
if (!dir.exists(panel_dir)) dir.create(panel_dir, recursive = TRUE)

purrr::iwalk(
  p_list,
  \(p, nm) {

    purrr::walk(
      c("svg", "pdf", "png"),
      \(fmt) cmplx_save(
        plot = p,
        file = file.path(
          panel_dir,
          paste0(
            "panel_E_",
            nm,
            output_name,
            ".",
            fmt
          )
        ),
        width = 4,
        height = 2.5
      )
    )
  }
)


## Panel F ====================================================================
## Panel F00 Annotation -------------------------------------------------------

# Sim. mode
p_sim_mode <- ggplot2::ggplot(
  data = pF_fgsea,
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
  data = pF_fgsea,
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
ggplot2::scale_fill_manual(
  name = "Workflow",
  values = c(
    "singlecell" = "grey10",
    "pseudobulk" = "grey60"
  )
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
  data = pF_fgsea,
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
  data = pF_fgsea,
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

# Enrich database
p_enrich_db <- ggplot2::ggplot(
  data = pF_fgsea,
  mapping = ggplot2::aes(
    x = method,
    y = "Enriched db.",
    fill = database
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
ggplot2::labs(
  x = NULL,
  y = NULL,
  fill = "Enriched database"
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
p_panel_F00 <- wrap_plots(
  p_sim_mode,
  p_workflow,
  p_method_type,
  p_gene_qc,
  p_enrich_db,
  ncol = 1
)

## Panel F01 ------------------------------------------------------------------
p_panel_F01 <- ggplot2::ggplot(
  data = stats::aggregate(
    I(!is.na(padj) & padj < target_alpha) ~ method + formula + celltype + knockout,
    data = pF_fgsea,
    FUN = sum
  ),
  mapping = ggplot2::aes(
    x = method,
    y = `I(!is.na(padj) & padj < target_alpha)`
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
  breaks = scales::pretty_breaks(n = 3)
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Method",
  y = "Nr. sig. pathways"
) +
theme_DK(
  rotate = 45,
  show_facet_x = FALSE,
  show_label_x = FALSE,
  show_ticks_x = FALSE,
  show_grid = FALSE
)

## Panel F02 ------------------------------------------------------------------
p_panel_F02 <- ggplot2::ggplot(
  data = pF_fgsea,
  mapping = ggplot2::aes(
    x = method,
    y = NES
  )
) +
ggplot2::geom_violin(
  color = "black",
  scale = "width",
  linewidth = 0.1
) +
ggplot2::scale_y_continuous(
  breaks = scales::pretty_breaks(n = 3)
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
 x = "Method",
 y = "Norm. enrichment scores"
) +
theme_DK(
  rotate = 45,
  show_facet_x = FALSE,
  show_label_x = TRUE,
  show_ticks_x = TRUE,
  show_grid = FALSE
)

## Patchwork & save -----------------------------------------------------------
panel_dir <- file.path(output_path, "panel_F")
if (!dir.exists(panel_dir)) dir.create(panel_dir, recursive = TRUE)

# Patchwork
p_panel_F <-
  p_panel_F00 +
  p_panel_F01 +
  p_panel_F02 +
  patchwork::plot_layout(
    ncol = 1,
    heights = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.4, 0.4),
    guides = "collect"
  ) 

# Save
purrr::walk(
  c("svg", "pdf", "png"),
  \(fmt) ggplot2::ggsave(
    plot = p_panel_F,
    file = file.path(
      panel_dir,
      paste0(
        "panel_F", 
        output_name, 
        ".", fmt
      )
    ),
    width = 10,
    height = 10.5,
    units = "cm",
    dpi = 300
  )
)

saveRDS(
  object = p_panel_F,
  file = file.path(
    panel_dir,
    paste0(
      "panel_F",
      output_name,
      ".rds"
    )
  )
)



## Script info ================================================================
showf("\n◷ %s ====================================================\n", date(), type = "cat")
