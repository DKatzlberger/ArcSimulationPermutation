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
  knockout  <- args[3]
  
  file_name <- paste0("prepare_knockout_", knockout)

} else {
  project_root <- "/home/scc_e_928652/ArcSimulationPermutation"
  input_path   <- "data/JB22_Spleen/rna_integrated.rds"
  output_path  <- "results/manual/prepare"

  # Meta
  knockout  <- "Cysltr2"

  file_name <- paste0("prepare_knockout_", knockout)
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
    "Seurat"
  )
)

## Dataset ====================================================================
.readRDS <- fun_with_ram(base::readRDS)

txt      <- sprintf("'%s' as 'BPCells Seurat'", input_path)
seu_obj  <- .readRDS(input_path, info = txt)

## Genotype & on_guide ========================================================
.assign_guide_by_max_count <- fun_with_ram(assign_guide_by_max_count)

seu_obj <- .assign_guide_by_max_count(seu_obj)
seu_obj <- subset(seu_obj, subset = genotype %in% c("NTC", knockout))

## Convert to SingleCellExperiment ============================================
matr <- SeuratObject::LayerData(seu_obj, assay = "RNA", layer = "counts")
meta <- seu_obj@meta.data

sce_obj <- SingleCellExperiment::SingleCellExperiment(
  assays  = list(counts = matr),
  colData = meta
)

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
