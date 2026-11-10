## max. code width ============================================================

## Title: Arc Institute in vivo pertubation screens
## Author: DKatzlberger
## Date: 21.01.2026

## Log/Printing ===============================================================
line80 <- function(x, char = "=") {
  paste(x, strrep(char, times = max(1, 78 - nchar(x))))
}

showf <- function(
  fmt, 
  ..., 
  type = c("message", "print", "cat", "warning", "stop")
) {
  type   <- match.arg(type)
  txt    <- sprintf(fmt, ...)
  caller <- sys.call(-1L)

  switch(
    type,
    message = message(txt),
    print   = print(txt),
    cat     = cat(txt),
    warning = warning(simpleWarning(txt, call = caller)),
    stop    = stop(simpleError(txt, call = caller))
  )

  invisible(txt)
}

fun_with_ram <- function(fun) {

  # Helper: extend '-' to line 80
  .line80 <- function(x) {
    paste(x, strrep("-", max(1, 78 - nchar(x))))
  }

  # Helper: job id  
  .scc_jobid <- function() {
    cg  <- readLines("/proc/self/cgroup", warn = FALSE)
    hit <- grep("job_[0-9]+", cg, value = TRUE)
    if (length(hit) == 0) return(NA_character_)
    sub(".*job_([0-9]+).*", "\\1", hit[1])
  }

  # Helper: Get system cgroup memory
  .get_system_mem <- function() {
    job_id  <- .scc_jobid()
    limit_gb <- NA_real_
    used_gb  <- NA_real_
    
    if (!is.na(job_id) && job_id != "") {
      x <- tryCatch(
        system(sprintf("scontrol show job %s", job_id), intern = TRUE),
        error = function(e) NULL
      )
      if (!is.null(x)) {
        line <- grep("AllocTRES=", x, value = TRUE)
        if (length(line) > 0) {
          limit_str <- sub(".*mem=([^,\\s]+).*", "\\1", line, perl = TRUE)
          
          num_val <- as.numeric(gsub("[^0-9.]", "", limit_str))
          unit    <- toupper(gsub("[0-9.]", "", limit_str))
          
          if (grepl("G", unit)) {
            limit_gb <- num_val          
          } else if (grepl("M", unit)) {
            limit_gb <- num_val / 1024   
          } else if (grepl("T", unit)) {
            limit_gb <- num_val * 1024  
          } else {
            limit_gb <- num_val         
          }
        }
      }
    }
    
    if (file.exists("/proc/self/status")) {
      status <- readLines("/proc/self/status", warn = FALSE)
      vsize  <- grep("^VmRSS:", status, value = TRUE)
      if (length(vsize) > 0) {
        kb <- as.numeric(gsub("[^0-9]", "", vsize))
        used_gb <- kb / (1024^2) # KB -> MB -> GB
      }
    }
    
    list(limit = limit_gb, usage = used_gb)
  }

  # Extract fun. name
  fun_name <- deparse(substitute(fun))
  fun_name <- paste("fun.:", fun_name)

  force(fun)
  force(fun_name)

  function(...) {

    args_list <- list(...)
    
    if (!is.null(names(args_list)) && "info" %in% names(args_list)) {
      info <- args_list[["info"]]
      args_list[["info"]] <- NULL
    } else {
      info <- NULL
    }

    cnd <- rlang::catch_cnd({
      t0 <- Sys.time()
      messages <- character()
      warnings <- character()
      std_logs <- character()  

      mem <- peakRAM::peakRAM({
        std_logs <- capture.output({
          visible_output <- withCallingHandlers(
            withVisible(do.call(fun, args_list, envir = parent.frame())),

            message = function(m) {
              messages <<- c(messages, conditionMessage(m))
              invokeRestart("muffleMessage")
            },

            warning = function(w) {
              warnings <<- c(warnings, conditionMessage(w))
              invokeRestart("muffleWarning")
            }
          )
        }, type = "output")
      })
      t1 <- Sys.time()

      output_obj <- visible_output$value

      # Stats
      time <- as.numeric(difftime(t1, t0, units = "secs"))
      size <- as.numeric(object.size(output_obj)) / 1024^2
      peak <- as.numeric(mem[[grep("^Peak", names(mem), value = TRUE)[1]]][1])

      # Read allocated system memory
      sys_mem <- .get_system_mem()

      # Report
      showf("\n\n", type = "cat")
      showf("%s\n", .line80(fun_name), type = "cat")
      if (!is.null(info)) showf("info: %s\n", info, type = "cat")
      showf("- time: %.2f seconds\n", time, type = "cat")
      showf("- size: %.2f MB\n", size, type = "cat")
      showf("- peak: %.2f MB\n", peak, type = "cat")
      showf("- used: %.2f GB (%.2f GB)\n", sys_mem$usage, sys_mem$limit, type = "cat")
      
      # Std. out
      std_logs <- std_logs[nchar(trimws(std_logs)) > 0]
      if (length(std_logs)) {
        showf("\noutput:\n", type = "cat")
        showf("%s\n", paste(std_logs, collapse = "\n"), type = "cat")
      }

      # Messages
      if (length(messages)) {
        showf("\nmessages:\n", type = "cat")
        showf("%s", paste(messages, collapse = ""), type = "cat")
      }

      # Warnings
      if (length(warnings)) {
        showf("\nwarnings:\n", type = "cat")
        showf("%s\n", paste(warnings, collapse = "\n"), type = "cat")
      }

      showf("%s\n", .line80("✓ executed"), type = "cat")

      # Return
      if (visible_output$visible) {
        return(output_obj)
      } else {
        return(invisible(output_obj))
      }
    }, classes = "error")

    if (!is.null(cnd)) {
      showf("\n\n", type = "cat")
      showf("%s\n", .line80(fun_name), type = "cat")
      if (!is.null(info)) showf("info: %s\n", info, type = "cat")
      
      tb <- rlang::trace_back(cnd)
      showf("\nerror:\n", type = "cat")
      showf("%s\n", conditionMessage(cnd), type = "cat")

      showf("\ntraceback:\n", type = "cat")
      tb     <- rlang::trace_back(cnd)
      std_tb <- capture.output(print(tb))
      showf("%s\n", paste(std_tb, collapse = "\n"), type = "cat")
      
      showf("%s\n", .line80("✗ failed"), type = "cat")

      # Script stop
      old_ops <- options(show.error.messages = FALSE)
      on.exit(options(old_ops))
      stop(structure(list(message = ""), class = c("simpleError", "error", "condition")))
    }
  }
}

## Environment ================================================================
scc_status <- function(verbose = TRUE) {

  .scc_jobid <- function() {
    cg <- readLines("/proc/self/cgroup", warn = FALSE)
    hit <- grep("job_[0-9]+", cg, value = TRUE)
    if (length(hit) == 0) return(NA_character_)
    sub(".*job_([0-9]+).*", "\\1", hit[1])
  }

  .slurm_mem <- function(jobid) {
    if (is.na(jobid)) return(NA_character_)
    x <- system(
      sprintf("scontrol show job %s", jobid),
      intern = TRUE
    )
    line <- grep("AllocTRES=", x, value = TRUE)
    sub(".*mem=([^,]+).*", "\\1", line)
  }

  # User & node
  job_id   <- .scc_jobid()
  username <- Sys.getenv("USER")
  node     <- Sys.info()[["nodename"]]
  n_cores  <- as.integer(system("nproc", intern = TRUE))
  memory   <- .slurm_mem(job_id)

  if (verbose){
    cat(sprintf("%s@%s\n", username, node))
    cat(sprintf("- scc job id: %s\n", job_id))
    cat(sprintf("- cores available: %d\n", n_cores))
    cat(sprintf("- mem available: %s\n", memory))
  }

  invisible(
    list(
      job_id = job_id,
      n_cores = n_cores,
      mem = .slurm_mem(job_id)
    )
  )
}

script_startup <- function(
  root,
  renv = FALSE,
  bioc = FALSE,
  packages = NULL
){

  # renv
  t0 <- Sys.time()
  if (renv) {

    # Environment variables
    csbg_mnt <- "/mnt/ceph/groups_hdd/SCCGroup/computational_systems_biology_group/dkatzlberger"
    Sys.setenv(RENV_PATHS_RENV = ".renv")
    Sys.setenv(RENV_PATHS_ROOT = file.path(csbg_mnt, ".renv"))
    Sys.setenv(RENV_PATHS_CACHE = file.path(csbg_mnt, ".renv/cache"))
    Sys.setenv(RENV_PATHS_SANDBOX = file.path(csbg_mnt, ".renv/sandbox"))

    Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
    Sys.setenv(RENV_CONFIG_CACHE_SYMLINKS = "TRUE")
    Sys.setenv(RENV_CONFIG_INSTALL_TRANSACTIONAL = "FALSE")
    Sys.setenv(RENV_CONFIG_AUTO_SNAPSHOT = "FALSE")
    Sys.setenv(RENV_CONFIG_PKG_TYPES = "source")

    # Activate renv
    renv::load(project = root)

    project <- tryCatch(renv::project(), error = function(e) NULL)
    version <- as.character(packageVersion("renv"))

    # BiocManager 
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        renv::install("BiocManager", prompt = FALSE)
        # BiocManager repositories 
        suppressMessages(options(repos = BiocManager::repositories()))
        renv::install("BiocVersion", prompt = FALSE)
      }
    }
  }
  t1 <- Sys.time()

  # Libraries (only if provided)
  if (!is.null(packages)) {

    if (!is.character(packages)) showf("`packages` must be a character vector\n", type = "stop")

    suppressWarnings(
      suppressPackageStartupMessages({
        for (pkg in unique(packages)) {
          library(pkg, character.only = TRUE, quietly = TRUE)
        }
      })
    )
  }

  t2 <- Sys.time()
  
  # Report
  renv_time <- as.numeric(difftime(t1, t0, units = "secs"))
  lib_time  <- as.numeric(difftime(t2, t1, units = "secs"))
  showf("- renv activation: %.2f seconds\n", renv_time, type = "cat")
  showf("- lib. attachment: %.2f seconds\n", lib_time, type = "cat")

  # Info
  showf("NOTE: R runs through pixi. Make sure all variables are set!\n", type = "cat")
  showf("NOTE: renv sandbox is deactivated!\n", type = "cat")
  
  
  # Return
  invisible(TRUE)
}

## CRISPR =====================================================================
check_obj <- function(obj) {

  class <- paste(class(obj), collapse = ", ")
  dims  <- dim(obj)

  showf(fmt = "Obj.: %s", class, type = "message")
  if (!is.null(dims) && length(dims) == 2) {
    showf(fmt = "Dim.: %s genes x %s cells", dims[1], dims[2], type = "message")
  }
}

assign_guide_by_max_count <- function(
  obj,
  assay = "CRISPR",
  ntc_pattern = "^NTCms|^NonTargeting|^BRDN|^NTC",
  verbose = TRUE
) {

  # Extract crispr counts
  counts <- obj@assays[[assay]]$counts

  # Assign guides
  obj$guide <- apply(counts, 2, function(x) {
    if (max(x) == 0) return(NA)
    rownames(counts)[which.max(x)]
  })

  # Find cells with NTC pattern
  is_ntc <- grepl(ntc_pattern, obj$guide, ignore.case = TRUE)

  # Assign new columns
  obj$genotype <- ifelse(is_ntc, "NTC", sub("-.*$", "", obj$guide))
  obj$on_guide <- ifelse(is_ntc, "NTC", obj$guide)

  # Validate guide assignment 
  valid_cells <- !is.na(obj$guide)
  genotypes_per_guide <- tapply(obj$genotype[valid_cells], obj$guide[valid_cells], function(x) length(unique(x)))
  problematic_guides  <- names(genotypes_per_guide[genotypes_per_guide > 1])

  if (length(problematic_guides) > 0) {
    showf("Identical guides were assigned to different genotypes", type = "stop")
  }
  
  # Verbose
  n_valid_cells <- sum(!is.na(obj$guide))
  n_na_cells    <- sum(is.na(obj$guide))
  n_genotypes   <- length(unique(obj$genotype)) 

  if (verbose) {
    showf("Assigned %s unique genotypes to %s cells.", n_genotypes, n_valid_cells, type = "message")
    if (n_na_cells > 0) showf("NOTE: %s cells couldn't be assigned.", n_na_cells, type = "message")
  }

  # Return
  obj
}

make_guide_target_df <- function(
  data,
  guide_col,
  target_col, 
  ref,
  verbose = TRUE
) {

  data <- as.data.frame(data)
  
  # Unique grna targets
  obj <- unique(data[, c(guide_col, target_col)])
  colnames(obj) <- c("guide_id", "guide_target")
  
  # Targeting vs. non-targeting
  obj$guide_target <- ifelse(
    obj$guide_target == ref, 
    "non-targeting", 
    sprintf("%s%s", target_col, obj$guide_target)
  )

  # Verbose
  n_guides  <- length(unique(obj$guide_id))
  n_targets <- length(unique(obj$guide_target))

  if (verbose) {
    showf("Found %d unique guides mapped to %d unique genotypes.", n_guides, n_targets, type = "message")
    showf("Recoding reference value '%s' to 'non-targeting'.", ref, type = "message")
  }

  # Return
  rownames(obj) <- NULL
  obj
}

make_pseudo_guide_matr <- function(
  data,
  guide_col,
  pseudocount = 1,
  verbose = TRUE
){

  data <- as.data.frame(data)

  cell_ids  <- rownames(data)
  guide_ids <- unique(data[[guide_col]])

  # Create matrix
  obj <- sapply(guide_ids, function(g) {
    as.numeric(data[[guide_col]] == g)
  })

  colnames(obj) <- guide_ids
  rownames(obj) <- cell_ids

  # Add pseudocount
  obj <- obj * pseudocount
  obj <- t(obj)

  # Verbose
  n_cells  <- length(cell_ids)
  n_guides <- length(guide_ids)

  if (verbose) {
    showf("Created pseudocounts for %d guides and %d cells.", n_guides, n_cells, type = "message")
    showf("Added a pseudocount of %d to assigned cells.", pseudocount, type = "message")
  }

  # Return
  obj
}

make_extra_cov_df <- function(
  data,
  target_col
){

  data <- as.data.frame(data)

  # Create extra covariate df
  obj <- data[, target_col, drop = FALSE]
  
  # Return
  obj
}

## DGE analysis ===============================================================
parse_model_formula <- function(
  formula,
  verbose = TRUE
) {

  # From string to formula
  f <- as.formula(formula)

  # Random effects
  bars <- reformulas::findbars(f)
  if (is.null(bars)) {
    rand_vars <- NULL
  } else {
    rand_vars <- unique(unname(unlist(lapply(bars, function(b) all.vars(b[[3]])))))
  }
  
  # Fixed effects
  fixed_vars <- labels(terms(reformulas::nobars(f)))
  if (length(fixed_vars) > 0) {
    target_covariate <- fixed_vars[1]
    extra_covariates <- if (length(fixed_vars) > 1) fixed_vars[-1] else NULL
  } else {
    target_covariate <- NULL
    extra_covariates <- NULL
  }

  # Verbose
  if (verbose) {
    showf("Formula: %s", formula, type = "message")
    showf("Target covariate: %s", target_covariate, type = "message")
    showf("Extra covariates: %s", paste(extra_covariates, collapse = ", "), type = "message")
    showf("Random covariate: %s", paste(rand_vars, collapse = ", "), type = "message")
  }
  
  # Return
  list(
    target_covariate  = target_covariate,
    extra_covariates  = extra_covariates,
    random_covariates = rand_vars
  )
}

validate_dge <- function(
  x
){
  
  # Check
  if (nrow(x) == 0) {
    showf("Statistic has zero rows. Maybe wrong extraction fun.", type = "stop")
  }

  # Universal req. columns
  required_cols <- c("coef", "gene", "logFC", "se", "pvalue")
  missing_cols  <- setdiff(required_cols, colnames(x))

  if (length(missing_cols) > 0) {
    showf(
      paste0("Statistic is missing req. column(s): ", 
      paste(missing_cols, collapse = ", ")), 
      type = "stop"
    )
  }

  # Check for NA columns
  required_cols <- c("pvalue")
  na_cols       <- colSums(is.na(x[required_cols]))

  if (any(na_cols > 0)) {
    showf(
      paste0("NOTE: Statistic has NA in column(s): ", 
      paste(names(na_cols)[na_cols > 0], collapse = ", ")), 
      type = "message"
    )
  }

  # Return
  invisible(TRUE)
}

extract_nebula <- function(
  fit,
  intercept = FALSE,
  adjust_p = TRUE,
  adjust_m = "BH"
){

  stats <- fit$summary

  # Identify coefficients
  coef_names <- unique(sub(
      "^(logFC|se|p)_",
      "",
      grep("^(logFC|se|p)_", colnames(stats), value = TRUE)
  ))

  # Build long-format table
  out <- do.call(
    rbind,
    lapply(coef_names, function(coef){

      data.frame(
        coef = coef,
        gene = stats$gene,
        logFC = stats[[paste0("logFC_", coef)]],
        se = stats[[paste0("se_", coef)]],
        pvalue = stats[[paste0("p_", coef)]],
        stringsAsFactors = FALSE
      )
    })
  )

  # Optionally remove intercept
  if (!intercept) out <- out[out$coef != "(Intercept)", ]

  # Optionally adjust p-values
  if (adjust_p) {
    out$FDR <- ave(
      out$pvalue,
      out$coef,
      FUN = function(z){ p.adjust(z, method = adjust_m)}
    )
  }

  # Validate
  rownames(out) <- NULL
  validate_dge(out)

  # Return
  out
}

extract_sceptre <- function(
  fit,
  intercept = FALSE,
  adjust_p = TRUE,
  adjust_m = "BH"
) {

  # Extract the results table
  stats <- if (inherits(fit, "sceptre_object")) {
    sceptre::get_result(fit, analysis = "run_discovery_analysis")
  } else {
    showf("fit must be a 'sceptre_object'", type = "stop")
  }
  
  out <- data.frame(
    coef = stats$grna_target,
    gene = stats$response_id,
    logFC = stats$log_2_fold_change,
    se = stats$se_fold_change,
    pvalue = stats$p_value,
    stringsAsFactors = FALSE
  )
  
  if (adjust_p) {
    out$FDR <- ave(
      out$pvalue,
      out$coef,
      FUN = function(z) { p.adjust(z, method = adjust_m) }
    )
  }
  
  # Validate
  rownames(out) <- NULL
  validate_dge(out)

  # Return
  out
}

extract_limma <- function(
  fit,
  intercept = FALSE,
  adjust_p = TRUE,
  adjust_m = "BH"
){

  coefs <- colnames(fit$coefficients)

  if (!intercept) {
    coefs <- setdiff(coefs, "(Intercept)")
  }

  out <- do.call(
    rbind,
    lapply(coefs, function(coef){

      idx <- match(coef, colnames(fit$coefficients))

      res <- data.frame(
        coef = coef,
        gene = rownames(fit$coefficients),
        logFC = fit$coefficients[, idx],
        se = fit$stdev.unscaled[, idx] * sqrt(fit$s2.post),
        pvalue = fit$p.value[, idx],
        stringsAsFactors = FALSE
      )

      res
    })
  )

  # Optionally adjust p-values
  if (adjust_p) {
    out$FDR <- ave(
      out$pvalue,
      out$coef,
      FUN = function(z){ p.adjust(z, method = adjust_m)}
    )
  }

  # Validate
  rownames(out) <- NULL
  validate_dge(out)

  # Return
  out
}

extract_edger <- function(
  fit,
  intercept = FALSE,
  adjust_p = TRUE,
  adjust_m = "BH"
){

  coefs <- colnames(fit$design)

  if (!intercept) {
    coefs <- setdiff(coefs, "(Intercept)")
  }

  out <- do.call(
    rbind,
    lapply(coefs, function(coef){

      qlf <- edgeR::glmQLFTest(fit, coef = coef)
      
      res_table <- qlf$table
      coef_idx  <- match(coef, colnames(fit$design))

      res <- data.frame(
        coef = coef,
        gene = rownames(qlf),
        logFC = res_table$logFC,
        se = NA_real_,
        pvalue = res_table$PValue,
        stringsAsFactors = FALSE
      )

      res
    })
  )

  # Optionally adjust p-values
  if (adjust_p) {
    out$FDR <- ave(
      out$pvalue,
      out$coef,
      FUN = function(z){ p.adjust(z, method = adjust_m)}
    )
  }

  # Validate
  rownames(out) <- NULL
  validate_dge(out)

  # Return
  out
}

extract_deseq2 <- function(
  fit,
  gene_qc = FALSE,
  intercept = FALSE,
  adjust_p = TRUE,
  adjust_m = "BH"
){
  coefs <- DESeq2::resultsNames(fit)

  if (!intercept) {
    coefs <- setdiff(coefs, "Intercept")
  }

  out <- do.call(
    rbind,
    lapply(coefs, function(coef){
      
      res_obj <- DESeq2::results(
        fit, 
        name = coef,
        cooksCutoff = gene_qc, 
        independentFiltering = gene_qc
      )

      # Rename to limma format
      clean_coef <- coef
      if (clean_coef == "Intercept") {
        clean_coef <- "(Intercept)"
      } else {
        clean_coef <- gsub("_vs_.*$", "", clean_coef)
        clean_coef <- gsub("_", "", clean_coef)
      }
      
      res <- data.frame(
        coef = clean_coef,
        gene = rownames(res_obj),
        logFC = res_obj$log2FoldChange,
        se = res_obj$lfcSE,
        pvalue = res_obj$pvalue,
        stringsAsFactors = FALSE
      )

      res
    })
  )

  # Optionally adjust p-values
  if (adjust_p) {
    out$FDR <- ave(
      out$pvalue,
      out$coef,
      FUN = function(z){ p.adjust(z, method = adjust_m)}
    )
  }

  # Validate
  rownames(out) <- NULL
  validate_dge(out)

  # Return
  out
}

extract_seurat <- function(
  fit,
  adjust_p = TRUE,
  adjust_m = "BH"
){

  out <- data.frame(
    coef = fit$coef,
    gene = rownames(fit),
    logFC = fit$avg_log2FC,
    se = NA_real_,
    pvalue = fit$p_val,
    stringsAsFactors = FALSE
  )

  # Optionally adjust p-values
  if (adjust_p) {
    out$FDR <- ave(
      out$pvalue,
      out$coef,
      FUN = function(z){ p.adjust(z, method = adjust_m)}
    )
  }

  # Validate
  rownames(out) <- NULL
  validate_dge(out)

  # Return
  out
}

extract_scmageck <- function(
  fit,
  input_ids,
  intercept = FALSE,
  adjust_p = TRUE,
  adjust_m = "BH"
) {

  .match_scmageck_ids <- function(ids, reference) {

    unmatched <- !ids %in% reference
    if (!any(unmatched)) return(ids)

    # Reverse observed make.names() conversions
    candidates <- make.names(reference, unique = FALSE)
    n_matches  <- vapply(ids[unmatched], function(x) sum(candidates == x), integer(1))
    if (any(n_matches > 1L)) showf(fmt = "Ambiguous mapping for %s scMAGeCK IDs.", sum(n_matches > 1L), type = "stop")

    # Match scMAGeCK IDs against converted input IDs
    idx <- match(ids[unmatched], candidates)
    if (anyNA(idx)) showf(fmt = "Unable to map %s scMAGeCK IDs to input IDs.", sum(is.na(idx)), type = "stop")
    
    # Restore original input IDs
    ids[unmatched] <- reference[idx]
    if (anyDuplicated(ids)) showf(fmt = "ID mapping produced duplicated input IDs.", type = "stop")

    # Return
    showf(fmt = "Successfully mapped all IDs.", type = "message")
    ids
  }

  score <- fit$score
  pval  <- fit$pvalue

  coefs <- score$Perturbedgene
  if (!intercept) coefs <- setdiff(coefs, "NegCtrl")
  
  score_genes <- setdiff(colnames(score), "Perturbedgene")
  pval_genes  <- setdiff(colnames(pval), "Perturbedgene")

  # Ensure score and p-value matrices have identical feature names
  if (!identical(score_genes, pval_genes)) showf("IDs differ between score and p-value matrices.", type = "stop")

  # Restore input identifiers only when necessary
  genes <- .match_scmageck_ids(ids = score_genes, reference = input_ids)

  out <- do.call(
    rbind,
    lapply(coefs, function(coef) {

      idx <- match(coef, score$Perturbedgene)

      data.frame(
        coef = coef,
        gene = genes,
        logFC = as.numeric(score[idx, score_genes]),
        se = NA_real_,
        pvalue = as.numeric(pval[idx, pval_genes]),
        stringsAsFactors = FALSE
      )
    })
  )

  if (adjust_p) {
    out$FDR <- ave(
      out$pvalue,
      out$coef,
      FUN = function(z) p.adjust(z, method = adjust_m)
    )
  }

  rownames(out) <- NULL
  validate_dge(out)
  # Return
  out
}

## ggplot2 ====================================================================
theme_DK <- function(
  base_size = 5,
  lgd_key_prop = 1,
  lgd_position = "right",
  make_void = FALSE,
  show_axis = TRUE,
  show_grid = TRUE,
  show_facet_x = TRUE,
  show_facet_y = TRUE,
  show_label_x = TRUE,
  show_label_y = TRUE,
  show_ticks_x = TRUE,
  show_ticks_y = TRUE,
  show_borders = TRUE,
  zero_margins = TRUE,
  show_identity = FALSE,
  xintercept = NULL,
  yintercept = NULL,
  rotate = NULL,
  base_family = "sans",
  ...
) {

  # Unit conversion helpers 
  mm_to_pt <- function(mm) mm / 0.352778
  pt_to_mm <- function(pt) pt * 0.352778

  # Set internal params. for void
  if (make_void) {
    show_axis <- FALSE
    show_grid <- FALSE
    show_label_x <- FALSE
    show_label_y <- FALSE
    show_ticks_x <- FALSE
    show_ticks_y <- FALSE
    show_borders <- FALSE
  }

  # Base theme
  p <- ggplot2::theme(
    plot.background  = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    strip.background = ggplot2::element_rect(fill = "white", color = NA),

    # Facets
    strip.placement = "inside",
    strip.text.x = ggplot2::element_text(size = base_size, family = base_family, margin = ggplot2::margin(b = 1, t = 0, unit = "mm")),
    strip.text.y = ggplot2::element_text(size = base_size, family = base_family, margin = ggplot2::margin(b = 1, t = 0, unit = "mm")),
    panel.spacing.x = grid::unit(1, "mm"),
    panel.spacing.y = grid::unit(1, "mm"),

    axis.line  = ggplot2::element_line(color = "black", linewidth = 0.3),
    axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = ggplot2::unit(1, "mm"),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = if (show_grid) ggplot2::element_line(color = "grey80", linewidth = 0.3) else ggplot2::element_blank(),

    # Text
    axis.text  = ggplot2::element_text(size = base_size, family = base_family),
    axis.title = ggplot2::element_text(size = base_size, family = base_family),
    plot.title = ggplot2::element_text(size = base_size, hjust = 0.5, family = base_family),
    plot.subtitle = ggplot2::element_text(size = base_size, hjust = 0.5, family = base_family),
    plot.caption = ggplot2::element_text(size = base_size, hjust = 0, family = base_family),
    ...
  )

  # Legend sizing
  lgd_key_size  <- base_size * lgd_key_prop

  p <- p + ggplot2::theme(
    legend.title = ggplot2::element_text(size = base_size, family = base_family, margin = ggplot2::margin(b = 1, unit = "mm")),
    legend.text = ggplot2::element_text(size = base_size, family = base_family),

    legend.key.height = ggplot2::unit(pt_to_mm(lgd_key_size), "mm"),
    legend.key.width = ggplot2::unit(pt_to_mm(base_size), "mm"),
    legend.key.spacing = ggplot2::unit(0.5, "mm"),
    legend.box.spacing = ggplot2::unit(1, "mm"),
    legend.background = ggplot2::element_rect(fill = NA, color = NA)
  )

  # Legend position
  theme_lgd_pos    <- lgd_position
  theme_lgd_inside <- NULL
  theme_lgd_just   <- NULL

  if (lgd_position %in% c("inside.top.left", "inside.top.right", "inside.bottom.left", "inside.bottom.right")) {
    theme_lgd_pos <- "inside"
    
    if (lgd_position == "inside.top.left") {
      theme_lgd_inside <- c(0.02, 0.98)
      theme_lgd_just   <- c(0, 1)
    } else if (lgd_position == "inside.top.right") {
      theme_lgd_inside <- c(0.98, 0.98)
      theme_lgd_just   <- c(1, 1)
    } else if (lgd_position == "inside.bottom.left") {
      theme_lgd_inside <- c(0.02, 0.02)
      theme_lgd_just   <- c(0, 0)
    } else if (lgd_position == "inside.bottom.right") {
      theme_lgd_inside <- c(0.98, 0.02)
      theme_lgd_just   <- c(1, 0)
    }
  }

  p <- p + ggplot2::theme(
    legend.position = theme_lgd_pos,
    legend.position.inside = theme_lgd_inside,
    legend.justification = theme_lgd_just
  )

  # Adaptive margins (zero margins)
  if (zero_margins) {
    p <- p + ggplot2::theme(
      legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
      plot.margin = ggplot2::margin(1, 1, 1, 1, unit = "mm")
    )
  }

  # Adaptive facets (x and y facet)
  if (!show_facet_x) p <- p + ggplot2::theme(strip.text.x = ggplot2::element_blank())
  if (!show_facet_y) p <- p + ggplot2::theme(strip.text.y = ggplot2::element_blank())

  # Adaptive boarders
  if (show_borders) { p <- p + ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.line = ggplot2::element_blank()
    )
  }

  # Adaptive axis
  if (!show_axis) {
    p <- p + ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(), 
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 0.5, b = 0, unit = "mm")), # pull x labels closer
      axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 0.5, l = 0, unit = "mm"))  # pull y labels closer
    )
  }

  # Adaptive rotation
  if (!is.null(rotate)) {
    angle <- rotate %% 360
    # Determine justification heuristically
    if (angle == 0) {
      hjust <- 0.5; vjust <- 0.5
    } else if (angle == 90) {
      hjust <- 1; vjust <- 0.5
    } else if (angle == 270) {
      hjust <- 0; vjust <- 0.5
    } else if (angle > 0 && angle < 90) {
      hjust <- 1; vjust <- 1
    } else if (angle > 90 && angle < 180) {
      hjust <- 1; vjust <- 0
    } else if (angle > 180 && angle < 270) {
      hjust <- 0; vjust <- 0
    } else { # e.g., 315°
      hjust <- 0; vjust <- 1
    }

    p <- p + ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = angle,
        hjust = hjust,
        vjust = vjust,
        family = base_family
      )
    )
  }

  if (!show_label_x) {
    p <- p + ggplot2::theme(
      axis.title.x = ggplot2::element_blank()
    )
  }

  if (!show_label_y) {
    p <- p + ggplot2::theme(
      axis.title.y = ggplot2::element_blank()
    )
  }

  if (!show_ticks_x) {
    p <- p + ggplot2::theme(
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank()
    )
  }

  if (!show_ticks_y) {
    p <- p + ggplot2::theme(
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank()
    )
  }

  # Extra layers
  layers <- list()

  # Identity line y = x
  if (show_identity) {
    layers <- append(
      layers,
      list(
        ggplot2::geom_abline(
          slope = 1,
          intercept = 0,
          linetype = "solid",
          linewidth = 0.3,
          color = "black"
        )
      )
    )
  }

   # Vertical intercept(s)
  if (!is.null(xintercept)) {
    layers <- append(
      layers,
      list(
        ggplot2::geom_vline(
          xintercept = xintercept,
          linetype = "solid",
          linewidth = 0.3,
          color = "black"
        )
      )
    )
  }

  # Horizontal intercept(s)
  if (!is.null(yintercept)) {
    layers <- append(
      layers,
      list(
        ggplot2::geom_hline(
          yintercept = yintercept,
          linetype = "solid",
          linewidth = 0.3,
          color = "black"
        )
      )
    )
  }

  # Return theme + layers
  return(c(list(p), layers))
}

add_corner_axis <- function(p, base_size = 5) {

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  xlab <- p$labels$x %||% "x"
  ylab <- p$labels$y %||% "y"

  # Apply coord_fixed
  build <- ggplot2::ggplot_build(p)

  xr <- build$layout$panel_params[[1]]$x.range
  yr <- build$layout$panel_params[[1]]$y.range

  dx <- diff(xr)
  dy <- diff(yr)

  x0 <- xr[1]
  y0 <- yr[1]

  # proportional arrow lengths
  arrow_x <- 0.1 * dx
  arrow_y <- 0.1 * dy

  # proportional label offsets
  offset_x <- 0.03 * dy
  offset_y <- 0.03 * dx

  p <- p + ggplot2::annotate(
    "segment",
    x = x0, 
    xend = x0 + arrow_x,
    y = y0, 
    yend = y0,
    arrow = grid::arrow(
      length = grid::unit(1, "mm"),
      type = "closed"
    ),
    linewidth = 0.3
  ) +
  ggplot2::annotate(
    "segment",
    x = x0, 
    xend = x0,
    y = y0, 
    yend = y0 + arrow_y,
    arrow = grid::arrow(
      length = grid::unit(1, "mm"),
      type = "closed"
    ),
    linewidth = 0.3
  ) +
  ggplot2::annotate(
    "text",
    x = x0 + arrow_x / 2,
    y = y0 - offset_x,
    label = xlab,
    hjust = 0.5,
    vjust = 1,
    size = base_size / 3
  ) +
  ggplot2::annotate(
    "text",
    x = x0 - offset_y,
    y = y0 + arrow_y / 2,
    label = ylab,
    angle = 90,
    hjust = 0.5,
    vjust = 0,
    size = base_size / 3
  )

  return(p)
}

## ComplexHeatmap =============================================================
cmplx_gp <- function(
  base_size = 5, 
  base_family = "sans", 
  ...
) {
  grid::gpar(
    fontsize = base_size,
    fontfamily = base_family,
    fontface = "plain",
    ...
  )
}

cmplx_bar <- function(
  x,
  base_size = 5,
  lgd_key_prop = 1,
  base_family = "sans",
  which = c("row", "column"),
  axis_side = NULL,
  flip = FALSE,
  ...
) {

  which <- match.arg(which)
  
  if (is.null(axis_side)) {
    axis_side <- if (which == "row") "top" else "left"
  } else {
    if (which == "row") {
      axis_side <- match.arg(axis_side, c("top", "bottom"))
    } else {
      axis_side <- match.arg(axis_side, c("left", "right"))
    }
  }

  gp_base <- cmplx_gp(base_size, base_family, lwd = 0.5)
  anno_height <- grid::unit(base_size * 2, "mm")

  max_val <- max(x, na.rm = TRUE)
  tick_positions <- seq(from = 0, to = max_val, length.out = 3) 

  tick_labels <- scales::label_number(
    scale_cut = scales::cut_short_scale(),
    accuracy = 1
  )(tick_positions)

  fun <- ComplexHeatmap::anno_barplot(
    x,
    which = which,
    border = FALSE,
    axis = TRUE,
    gp = grid::gpar(
      fill = "grey40", 
      col = NA
    ),
    axis_param = list(
      gp = gp_base,
      labels_rot = 90,
      side = axis_side,
      at = tick_positions,
      labels = tick_labels,
      direction = if (flip) "reverse" else "normal" 
    ),
    width = if (which == "row") anno_height else NULL,
    height = if (which == "column") anno_height else NULL
  )

  # Return
  fun
}

cmplx_anno <- function(
  ...,
  base_size = 5,
  lgd_key_prop = 1,
  base_family = "sans",
  annotation_legend_param = list(),
  which = c("column", "row")
) {

  which <- match.arg(which)

  dots <- list(...)
  name_side <- dots$annotation_name_side

  if (!is.null(name_side)) {
    if (which == "column" && !name_side %in% c("left", "right")) {
      stop("For column annotations, annotation_name_side must be 'left' or 'right'")
    }
    if (which == "row" && !name_side %in% c("top", "bottom")) {
      stop("For row annotations, annotation_name_side must be 'top' or 'bottom'")
    }
  }

  pt_to_mm <- function(pt) pt * 0.352778

  gp_base <- cmplx_gp(base_size, base_family)
  anno_height <- grid::unit(pt_to_mm(base_size), "mm")
  lgd_unit <- grid::unit(pt_to_mm(base_size * lgd_key_prop), "mm")

  default_legend_param <- list(
    title_position = "topleft",
    title_gp = gp_base,
    labels_gp = gp_base,
    gap = grid::unit(0.5, "mm"),
    grid_width = lgd_unit,
    grid_height = lgd_unit,
    legend_height = lgd_unit
  )

  n_anno <- sum(!names(dots) %in% c(
    "col", "na_col", "gp", "border",
    "annotation_name_side", "annotation_name_offset",
    "annotation_name_rot", "annotation_name_gp",
    "annotation_height", "annotation_width",
    "simple_anno_size", "show_annotation_name",
    "show_legend", "which", "gap",
    "annotation_legend_param"
  ))

  if (length(annotation_legend_param) > 0) {
    annotation_legend_param <- lapply(
      annotation_legend_param,
      \(x) utils::modifyList(default_legend_param, x)
    )
  }

  anno <- ComplexHeatmap::HeatmapAnnotation(
    ...,
    which = which,
    gap = grid::unit(0.7, "mm"),
    simple_anno_size = anno_height,
    annotation_name_gp = gp_base,
    annotation_legend_param = if (length(annotation_legend_param) == 0) {
      rep(list(default_legend_param), n_anno)
    } else {
      annotation_legend_param
    }
  )

  anno
}

cmplx_legend <- function(
  title,
  labels,
  base_size = 5,
  lgd_key_prop = 1,
  base_family = "sans",
  ...
) {

  mm_to_pt <- function(mm) mm / 0.352778
  pt_to_mm <- function(pt) pt * 0.352778

  gp_base <- cmplx_gp(base_size, base_family)

  lgd_key_size <- base_size * lgd_key_prop
  lgd_unit <- grid::unit(pt_to_mm(lgd_key_size), "mm")

  ComplexHeatmap::Legend(
    title = title,
    labels = labels,
    title_gp = gp_base,
    labels_gp = gp_base,
    grid_width = lgd_unit,
    legend_height = lgd_unit,
    ...
  )
}

cmplx_ht <- function(
  mat,
  heatmap_legend_param = list(),
  ...,
  base_size = 5,
  lgd_key_prop = 1,
  base_family = "sans",
  cell_size = NULL
) {

  stopifnot(is.matrix(mat) || is.data.frame(mat))
  mat <- as.matrix(mat)

  # Dimensions
  n_row <- nrow(mat)
  n_col <- ncol(mat)

  # Unit conversion helpers 
  mm_to_pt <- function(mm) mm / 0.352778
  pt_to_mm <- function(pt) pt * 0.352778

  # Annotation scaling
  gp_base <- cmplx_gp(base_size, base_family)
  anno_height <- grid::unit(pt_to_mm(base_size), "mm")

  # Legend scaling
  lgd_key_size <- base_size * lgd_key_prop
  lgd_unit <- grid::unit(pt_to_mm(lgd_key_size), "mm")

  # Conditionally enforce square tiles
  size_args <- list()
  if (!is.null(cell_size)) {
    if (is.numeric(cell_size)) {
      cell_size <- grid::unit(cell_size, "mm")
    }
    size_args$width  <- cell_size * n_col
    size_args$height <- cell_size * n_row
  }

  ht <- do.call(ComplexHeatmap::Heatmap, c(list(
    mat,
    ...,

    # Heatmap text
    row_names_gp = gp_base,
    column_names_gp = gp_base,
    column_title_gp = gp_base,
    row_title_gp = gp_base,

    # Heatmap legend
    heatmap_legend_param = utils::modifyList(
      list(
        title_gp = gp_base,
        labels_gp = gp_base,
        grid_width = lgd_unit,
        legend_height = lgd_unit
      ),
      val = heatmap_legend_param
    ),

    # Dendrogram styling
    row_dend_gp = grid::gpar(lwd = 0.3),
    column_dend_gp = grid::gpar(lwd = 0.3)
  ), size_args))

  # Return
  return(ht)
}

cmplx_save <- function(
  plot,
  file,
  width,
  height,
  units = "in",
  device = NULL,
  res = 300,
  anno_lgd_list = NULL,
  ...
) {

  to_inches <- function(x, units) {
    switch(
      EXPR = units,
      "in" = x,
      "cm" = x / 2.54,
      "mm" = x / 25.4,
      "px" = x / res,
      stop("Unsupported units: ", units)
    )
  }

  width_in  <- to_inches(width, units)
  height_in <- to_inches(height, units)

  if (is.null(device)) {
    ext <- tools::file_ext(file)
    device <- switch(
      EXPR = tolower(ext),
      pdf = grDevices::pdf,
      svg = grDevices::svg,
      png = grDevices::png,
      jpeg = grDevices::jpeg,
      jpg = grDevices::jpeg,
      tiff = grDevices::tiff,
      stop("Unsupported file type: ", ext)
    )
  }

  if (identical(device, grDevices::png) ||
      identical(device, grDevices::jpeg) ||
      identical(device, grDevices::tiff)) {

    device(
      file,
      width = width_in,
      height = height_in,
      units = "in",
      res = res
    )

  } else {

    device(
      file,
      width = width_in,
      height = height_in
    )
  }

  ComplexHeatmap::ht_opt(legend_gap = grid::unit(2.5, "mm"))

  # grid::grid.newpage()
  ComplexHeatmap::draw(
    plot,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    annotation_legend_list = anno_lgd_list,
    legend_grouping = "original",
    align_heatmap_legend = "heatmap_top",
    align_annotation_legend = "heatmap_top",
    merge_legend = FALSE,
    ...
  )

  invisible(grDevices::dev.off())
}


## Aestethics =================================================================
colors_formula <- c(
  "~genotype"                           = "#ECCBAE",
  "~genotype + sample"                  = "#046C9A",
  "~genotype + (1|sample)"              = "#D69C4E",
  "~genotype + sample + (1|guide)"      = "#ABDDDE",
  "~genotype + (1|sample) + (1|guide)"  = "#000000"
)

colors_method <- c(
  "limma"    = "#66C2A5",
  "edger"    = "#FC8D62",
  "deseq2"   = "#8DA0CB",
  "dream"    = "#E78AC3",
  "nebula"   = "#A6D854",
  "wilcox"   = "#FFD92F",
  "sceptre"  = "#E5C494",
  "scmageck" = "#B3B3B3"
)

colors_celltype <- c(
  "Ery"         = "#762A83",
  "Imm. B-cell" = "#9970AB",
  "Gran. P"     = "#D01C8B",
  "GMP"         = "#E76289",
  "Gran."       = "#B2182B",
  "MEP"         = "#8C8D27",
  "HSC"         = "#4D9221",
  "MkP"         = "#1B7837"
)

shapes_knockout <- c(
  "Brd9" = 16,
  "Cebpa" = 17,
  "Chd4" = 15,
  "Kmt2d" = 18,
  "Rcor1" = 3,
  "Klf1" = 4
)