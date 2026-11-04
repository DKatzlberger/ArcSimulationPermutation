## max. code width ============================================================

## Title: Arc Institute in vivo pertubation screens
## Author: DKatzlberger
## Date: 21.01.2026

## Shuffling knockout identification (meta data) to destroy any
## true associations between gene expression and the knockout.

## The data is saved in its original form if the seed is '1'. For other seeds
## the meta data columns 'genotype', 'guide' and 'on_guide' are shuffled.


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
  input_path   <- gsub(pattern = "\\\\ ", replacement = " ", trail_args[2])
  output_path  <- "."

  # Params
  seed <- as.numeric(trail_args[3])
  
  # Output
  input_name  <- tools::file_path_sans_ext(basename(input_path))
  output_name <- paste0(script_id, "-", input_name)
  output_name <- paste0(output_name, "_seed_", seed)

} else {
  # Script
  script_path <- rstudioapi::getActiveDocumentContext()$path
  script_name <- tools::file_path_sans_ext(basename(script_path))
  script_id   <- sub(pattern = "-.*", replacement = "", script_name)

  # Paths
  project_root <- "/home/scc_e_928652/ArcSimulationSynthetic"
  input_path   <- "results/manual/ArcSimPerm_GSE213511/01-prepare/01-celltype_Gran. P_knockout_Cebpa.h5ad"
  output_path  <- file.path("results/manual/ArcSimPerm_GSE213511", script_name)

  stopifnot(dirname(dirname(input_path)) == dirname(output_path))
  if (!dir.exists(output_path)) dir.create(path = output_path, recursive = TRUE)

  # Params
  seed <- 42
  
  # Output
  input_name  <- tools::file_path_sans_ext(basename(input_path))
  output_name <- paste0(script_id, "-", input_name)
  output_name <- paste0(output_name, "_seed_", seed)
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
    "SingleCellExperiment"
  )
)

## Script info ================================================================
showf(fmt = "\n%s\n", line80(paste0(script_name, ".R")), type = "cat")

## Dataset ====================================================================
.read_rds  <- fun_with_ram(base::readRDS)
.read_h5ad <- fun_with_ram(anndataR::read_h5ad)
.check_obj <- fun_with_ram(check_obj)

ext <- tools::file_ext(input_path)
txt <- sprintf("%s as 'SingleCellExperiment'", input_path)

sce_obj <- switch(EXPR = ext,
  rds  = .read_rds(input_path, info = txt),
  h5ad = .read_h5ad(input_path, as = "SingleCellExperiment", info = txt),
  showf("Unsupported file extension: %s", ext, type = "stop")
)

# Check input
txt <- "checking input obj."
.check_obj(sce_obj, info = txt)

# Gene truth (pipeline)
sce_rd <- SingleCellExperiment::rowData(sce_obj)
sce_rd <- as.data.frame(sce_rd)

## Shuffle ====================================================================
if (identical(seed, 1)) {

  sce_rd$seed <- as.character(seed)
  sce_rd$type <- "real"
  
} else {
  set.seed(seed)

  idx <- sample(seq_len(ncol(sce_obj)))

  sce_obj$genotype <- sce_obj$genotype[idx]
  sce_obj$guide    <- sce_obj$guide[idx]
  sce_obj$on_guide <- sce_obj$on_guide[idx]

  sce_rd$seed <- as.character(seed)
  sce_rd$type <- "perm"
}

# Add the rowData
SummarizedExperiment::rowData(sce_obj) <- NULL
SummarizedExperiment::rowData(sce_obj) <- sce_rd

## Save =======================================================================
.write_h5ad <- fun_with_ram(anndataR::write_h5ad)
txt         <- sprintf("saved as %s", paste0(output_name, ".h5ad"))

# Integer check
stopifnot(all(SingleCellExperiment::counts(sce_obj) %% 1 == 0))

# Check output
txt <- "checking output obj."
.check_obj(sce_obj, info = txt)

.write_h5ad(
  object = sce_obj, 
  path = file.path(
    output_path,
    paste0(output_name, ".h5ad")
  ),
  mode = "w",
  info = txt
)

## Script info ================================================================
showf("\n◷ %s ====================================================\n", date(), type = "cat")
