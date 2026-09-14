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
  seed     <- args[4]

  file_name <- paste0("pseudobulk_knockout_", knockout, "_seed_", seed)

} else {
  project_root <- "/home/scc_e_928652/ArcSimulationSynthetic"
  input_path   <- "results/nextflow/shuffle/shuffle_knockout_Brd7_seed_34.h5ad"
  output_path  <- "results/manual/pseudobulk"

#   stopifnot(dirname(dirname(input_path)) == dirname(output_path))

  # Meta
  knockout <- "Brd7"
  seed     <- 34

  file_name <- paste0("pseudobulk_knockout_", knockout, "_seed_", seed)
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

## Pseudobulk =================================================================
.aggregate_expression <- fun_with_ram(scuttle::aggregateAcrossCells)

groups  <- c("orig.ident", "guide")
txt     <- sprintf("pseudobulk by %s", paste(groups, collapse = ", "))
sce_obj <- .aggregate_expression(
  sce_obj, 
  ids = colData(sce_obj)[, groups], 
  statistics = "sum",
  use.assay.type = "counts",
  info = txt
)

# Check for no NAs
stopifnot(knockout %in% sce_obj$genotype)
stopifnot("NTC" %in% sce_obj$genotype)

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
