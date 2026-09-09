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

  # Meta
  knockout <- args[3]
  seed     <- as.numeric(args[4])
  
  file_name <- paste0("shuffle_knockout_", knockout, "_seed_", seed)

} else {
  project_root <- "/home/scc_e_928652/ArcSimulationPermutation"
  input_path   <- "results/manual/prepare/prepare_knockout_Cysltr2.h5ad"
  output_path  <- "results/manual/shuffle"

  # Meta
  knockout  <- "Cysltr2"
  seed      <- 1
  
  file_name <- paste0("shuffle_knockout_", knockout, "_seed_", seed)
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

## Dataset ====================================================================
.read_h5ad <- fun_with_ram(anndataR::read_h5ad)

txt     <- sprintf("%s as 'SingleCellExperiment'", input_path)
sce_obj <- .read_h5ad(input_path, as = "SingleCellExperiment", info = txt)

# Create rowData
sce_rd <- S4Vectors::DataFrame(gene = rownames(sce_obj))

## Shuffle ====================================================================
if (identical(seed, 1)) {

  sce_rd$seed <- seed
  sce_rd$type <- "real"
  
} else {
  set.seed(seed)
  sce_obj$on_guide <- sample(sce_obj$on_guide)
  sce_obj$genotype <- sample(sce_obj$genotype)

  sce_rd$seed <- seed
  sce_rd$type <- "perm"
}

# Add the rowData
SummarizedExperiment::rowData(sce_obj) <- NULL
SummarizedExperiment::rowData(sce_obj) <- sce_rd

## Save =======================================================================
.write_h5ad <- fun_with_ram(anndataR::write_h5ad)
txt         <- sprintf("saved as %s", paste0(file_name, ".h5ad"))

# Integer check
stopifnot(all(SingleCellExperiment::counts(sce_obj) %% 1 == 0))

.write_h5ad(
  object = sce_obj, 
  path = file.path(
    output_path,
    paste0(file_name, ".h5ad")
  ),
  mode = "w",
  info = txt
)