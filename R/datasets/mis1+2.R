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
  output_path  <- "data_generated/MIS1+2"

  # Count data
  count_path <- "data_raw/MIS1+2/rna_integrated.rds"

  # Meta data
  meta_path <- "data_raw/MIS1+2/cell_metadata.rds"
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
### Cell meta -----------------------------------------------------------------
.read_rds <- fun_with_ram(base::readRDS)

cell_meta <- .read_rds(meta_path)

### Counts --------------------------------------------------------------------
seu_obj <- .read_rds(count_path)

### Make names ----------------------------------------------------------------
colnames(seu_obj@meta.data)[colnames(seu_obj@meta.data) == "Guide"]  <- "guide"
colnames(seu_obj@meta.data)[colnames(seu_obj@meta.data) == "cell_type_immgen_medium"] <- "celltype"
colnames(seu_obj@meta.data)[colnames(seu_obj@meta.data) == "Sample"]  <- "sample"


## Save =======================================================================
.check_obj <- fun_with_ram(check_obj)
.save_rds  <- fun_with_ram(base::saveRDS)

# Check output
txt <- "checking output obj."
.check_obj(seu_obj, info = txt)

# Cell meta
txt <- sprintf("saved as %s", "cell_metadata.rds")
.save_rds(
  object = cell_meta, 
  file = file.path(
    output_path,
    "cell_metadata.rds"
  ),
  info = txt
)

# RNA integrated counts
txt <- sprintf("saved as %s", "rna_integrated.rds")
.save_rds(
  object = seu_obj, 
  file = file.path(
    output_path,
    "rna_integrated.rds"
  ),
  info = txt
)