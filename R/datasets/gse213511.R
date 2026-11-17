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

  data_raw_path <- file.path("data_raw", "GSE213511")

} else {
  project_root <- "/home/scc_e_928652/ArcSimulationPermutation"
  output_path  <- "data_generated/GSE213511"

  data_raw_path <- file.path("data_raw", "GSE213511")

  # Count data
  count_url <- c(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213511/suppl/GSE213511%5FinVivo%5FOP1%5Fckit%5F14d%5F1%2Eh5",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213511/suppl/GSE213511%5FinVivo%5FOP1%5Flin%2D%5F14d%5F1%2Eh5",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213511/suppl/GSE213511%5FinVivo%5FOP1%5Flin%2D%5F14d%5F2%2Eh5",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213511/suppl/GSE213511%5FinVivo%5FOP1%5Flin%2D%5F28d%5F1%2Eh5",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213511/suppl/GSE213511%5FinVivo%5FOP1%5Flin%2D%5F28d%5F2%2Eh5",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213511/suppl/GSE213511%5FinVivo%5FOP2%5Fckit%5F14d%5F1%2Eh5",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213511/suppl/GSE213511%5FinVivo%5FOP2%5Flin%2D%5F14d%5F1%2Eh5",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213511/suppl/GSE213511%5FinVivo%5FOP3%5Fckit%5F14d%5F1%2Eh5",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213511/suppl/GSE213511%5FinVivo%5FOP3%5Flin%2D%5F14d%5F1%2Eh5",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213511/suppl/GSE213511%5FinVivo%5FOP4%5Fckit%5F14d%5F1%2Eh5",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213511/suppl/GSE213511%5FinVivo%5FOP4%5Flin%2D%5F14d%5F1%2Eh5"
  )

  # Meta data
  meta_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213511/suppl/GSE213511%5FCellAnnotation%5Finvivo%2Etsv%2Egz"
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
dest_file <- file.path(data_raw_path, basename(URLdecode(meta_url)))

if (!file.exists(dest_file)) download.file(meta_url, dest_file, mode = "wb", method = "libcurl")

cell_meta <- read.delim(gzfile(dest_file), stringsAsFactors = FALSE, check.names = FALSE)

### Counts --------------------------------------------------------------------
seu_list <- list()
for (url in count_url) {
  file_name   <- basename(URLdecode(url))
  sample_name <- sub("^GSE213511_", "", sub("\\.h5$", "", file_name))
  dest_file   <- file.path(data_raw_path, file_name)

  if (!file.exists(dest_file)) download.file(url, dest_file, mode = "wb", method = "libcurl")
  
  # Seurat object
  f <- Seurat::Read10X_h5(dest_file)

  seu_obj <- Seurat::CreateSeuratObject(counts = f[["Gene Expression"]])
  seu_obj[["CRISPR"]] <- SeuratObject::CreateAssay5Object(counts = f[["CRISPR Guide Capture"]])

  # Add meta
  meta <- cell_meta[cell_meta$Sample == sample_name, , drop = FALSE]
  idx  <- match(colnames(seu_obj), meta$Barcode)

  if (anyDuplicated(meta$Barcode)) showf("Duplicate barcodes in meta for sample: %s", sample_name, type = "stop")
  if (all(is.na(idx))) showf("No matching barcodes for sample: %s", sample_name, type = "stop")

  if (anyNA(idx)) {
    showf("Warning: %s cells in %s have no annotation and will be removed.\n", sum(is.na(idx)), file_name, type = "cat")
    seu_obj <- subset(seu_obj, cells = colnames(seu_obj)[!is.na(idx)])
    idx <- idx[!is.na(idx)]
  }

  # Match meta for remaining cells
  matched_meta <- meta[idx, , drop = FALSE]
  rownames(matched_meta) <- colnames(seu_obj)

  seu_obj <- Seurat::AddMetaData(
    seu_obj,
    metadata = matched_meta[, setdiff(colnames(matched_meta), "Barcode"), drop = FALSE]
  )

  # Append
  seu_list[[sample_name]] <- seu_obj
}

# Merge samples
seu_obj <- merge(
  x = seu_list[[1]],
  y = seu_list[-1],
  merge.data = FALSE
)

seu_obj[["RNA"]]    <- SeuratObject::JoinLayers(seu_obj[["RNA"]])
seu_obj[["CRISPR"]] <- SeuratObject::JoinLayers(seu_obj[["CRISPR"]])

### Make names ----------------------------------------------------------------
colnames(seu_obj@meta.data)[colnames(seu_obj@meta.data) == "Guide"]    <- "guide"
colnames(seu_obj@meta.data)[colnames(seu_obj@meta.data) == "Clusters"] <- "celltype"
colnames(seu_obj@meta.data)[colnames(seu_obj@meta.data) == "Sample"]   <- "sample"



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