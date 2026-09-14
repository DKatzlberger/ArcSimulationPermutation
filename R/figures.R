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
  project_root <- "/home/scc_e_928652/ArcSimulationSynthetic"
  input_path   <- "results/nextflow/ArcSimPerm_JB22_Bcells"
  output_path  <- "results/nextflow/ArcSimPerm_JB22_Bcells/figures"

  alpha <- 0.05

  file_name <- paste0("_alpha_", alpha)
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

  # Found files
  exp_found <- list.files(
    path = path, 
    pattern = "_dge\\.rds$", 
    full.names = TRUE, 
    recursive = TRUE
  )

  # Check for duplication
  any_duplicated_paths <- any(duplicated(exp_found))
  if (any_duplicated_paths) {
    n_dup_paths <- sum(duplicated(exp_found))
    showf("Note: There are %s duplicated DGE files.", n_dup_paths)
  }

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
      "NOTE: Not all files could be loaded. %s failed files.",
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
  do.call(rbind, loaded_data)
}

dge <- readDGE(
  path = input_path,
  pattern = "~genotype.*_dge\\.rds$"
)

dge$formula  <- factor(
  dge$formula, 
  levels = c(
    "~genotype", 
    "~genotype + orig.ident"
  )
)

## Stats ======================================================================
group_vars <- c("workflow", "method_type", "method", "genotype",  "formula", "type", "seed")

stats <- aggregate(
  cbind(
    P = dge$pvalue < alpha
  ),
  by  = dge[group_vars],
  FUN = function(x) sum(x, na.rm = TRUE)
)


## Panel A ====================================================================
## Panel A0: Annotations ======================================================

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

# Patchwork
p_fig_1_panel_A0 <- wrap_plots(
  p_workflow,
  p_method_type,
  ncol = 1
)

## Panel A1: Positives ========================================================
p_fig_1_panel_A1 <- ggplot2::ggplot(
  mapping = ggplot2::aes(
    x = method,
    y = P
  )
) + 
ggplot2::stat_summary(
  data = stats[stats$type == "real", ],
  mapping = ggplot2::aes(
    group = genotype,
    fill = "Potential true positives (TP)"
  ),
  fun = "mean",
  geom = "col",
  position = ggplot2::position_dodge(
    width = 0.6
  ),
  width = 0.6,
  color = "black",
  linewidth = 0.1
) +
ggplot2::stat_summary(
  data = stats[stats$type == "real", ],
  mapping = ggplot2::aes(
    group = genotype,
    shape = genotype,
    fill  = "Total positives (TP + FP)"
  ),
  fun = "mean",
  geom = "point",
  position = ggplot2::position_dodge(
    width = 0.6
  ),
  size = 1
) +
ggplot2::stat_summary(
  data = stats[stats$type == "perm", ],
  mapping = ggplot2::aes(
    group = genotype,
    fill = "False positives (FP)"
  ),
  fun = "mean",
  geom = "col",
  position = ggplot2::position_dodge(
    width = 0.6
  ),
  width = 0.6,
  color = "black",
  linewidth = 0.1
) +
ggplot2::stat_summary(
  data = stats[stats$type == "perm", ],
  mapping = ggplot2::aes(
    group = genotype
  ),
  fun.data = "mean_sdl",
  fun.args = list(mult = 1),
  geom = "errorbar",
  position = ggplot2::position_dodge(
    width = 0.6
  ),
  color = "black",
  width = 0.2,
  linewidth = 0.2
) +
ggplot2::scale_fill_manual(
  name = NULL,
  values = c(
    "Total positives (TP + FP)"     = "blue",
    "Potential true positives (TP)" = "orange",
    "False positives (FP)"          = "grey80"
  ),
  breaks = c(
    "Total positives (TP + FP)",
    "Potential true positives (TP)", 
    "False positives (FP)"
  )
) +
ggplot2::facet_grid(
  cols = ggplot2::vars(formula),
  scales = "free_x",
  space = "free"
) +
ggplot2::labs(
  x = "Method", 
  y = "Nr. DEGs",
  shape = "Knockout"
) +
theme_DK(
  rotate = 45,
  show_facet_x = FALSE,
  show_grid = FALSE
)

## Patchwork & save ===========================================================
p_fig_1_panel_A <-
  p_fig_1_panel_A0 + 
  p_fig_1_panel_A1 + 
  patchwork::plot_layout(
    ncol = 1,
    heights = c(0.05, 0.05, 0.3),
    guides = "collect"
  ) 

# Save
ggplot2::ggsave(
  p_fig_1_panel_A,
  file = file.path(
    output_path,
    paste0(
      "panel_A",
      file_name,
      ".svg"
    )
  ),
  width = 8,
  height = 9.5,
  units = "cm"
)

saveRDS(
  p_fig_1_panel_A,
  file = file.path(
    output_path,
    paste0(
      "panel_A",
      file_name,
      ".rds"
    )
  ),
)


