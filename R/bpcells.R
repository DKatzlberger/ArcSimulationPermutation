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
  experiment <- args[3]
  file_name  <- paste0("bpcells_", experiment)

} else {
  project_root <- "/home/scc_e_928652/ArcSimulationPermutation"
  input_path   <- "data/JB22_Spleen/rna_integrated.rds"
  output_path  <- "results/manual/bpcells"

  # Meta
  experiment <- basename(dirname(input_path))
  file_name  <- paste0("bpcells_", experiment)
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

txt      <- sprintf("'%s' as 'Seurat'", input_path)
seu_obj  <- .readRDS(input_path, info = txt)

rna_counts    <- seu_obj[["RNA"]]$counts
crispr_counts <- seu_obj[["CRISPR"]]$counts

stopifnot(all(rna_counts@x %% 1 == 0), all(crispr_counts@x %% 1 == 0))

## Convert to disk readable ===================================================
.write_bpcells <- fun_with_ram(BPCells::write_matrix_dir)
.convert_type  <- fun_with_ram(BPCells::convert_matrix_type)

txt        <- "converting RNA assay"
rna_counts <- .convert_type(rna_counts, "uint32_t", info = txt)
rna_counts <- .write_bpcells(
  mat = rna_counts,
  dir = file.path(
    output_path,
    paste0(file_name, "_RNA")
  ),
  info = txt
)

txt           <- "converting CRISPR assay"
crispr_counts <- .convert_type(crispr_counts, "uint32_t", info = txt)
crispr_counts <- .write_bpcells(
  mat = crispr_counts,
  dir = file.path(
    output_path,
    paste0(file_name, "_CRISPR")
  ),
  info = txt
)

## BPcells seurat obj. ========================================================
seu_obj <- Seurat::CreateSeuratObject(
  counts = rna_counts,
  assay = "RNA",
  meta.data = seu_obj@meta.data
)

# Add CRISPR assay
seu_obj[["CRISPR"]] <- SeuratObject::CreateAssay5Object(
  counts = crispr_counts
)

## Save =======================================================================
.saveRDS <- fun_with_ram(base::saveRDS)
txt      <- sprintf("saved as %s", paste0(file_name, ".rds"))

.saveRDS(
  seu_obj,
  file.path(
    output_path,
    paste0(file_name, ".rds")
  ),
  info = txt
)
