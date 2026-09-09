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
  seed <- args[3]
  
  file_name <- paste0("shuffled_seed_", seed)

} else {
  project_root <- "/home/scc_e_928652/ArcSimulationPermutation"
  input_path   <- "results/manual/bpcells/bpcells.rds"
  output_path  <- "results/manual/shuffle"

  # Meta
  seed     <- 42
  genotype <- "Cysltr2"
  
  file_name <- paste0("shuffled_seed_", seed)
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
    "Seurat",
    "BPCells"
  )
)

## Dataset ====================================================================
.readRDS <- fun_with_ram(base::readRDS)

txt      <- sprintf("'%s' as 'BPCells Seurat'", input_path)
seu_obj  <- .readRDS(input_path, info = txt)

## Shuffle ====================================================================
set.seed(seed)

seu_obj$on_guide <- sample(seu_obj$on_guide)
seu_obj$genotype <- sample(seu_obj$genotype)

## Save =======================================================================
.saveRDS <- fun_with_ram(base::saveRDS)
txt      <- sprintf("saved as %s", paste0(file_name, ".rds"))

.saveRDS(
  sce_obj, 
  file.path(
    output_path,
    paste0(file_name, ".rds")
  ),
  info = txt
)