#!/usr/bin/env Rscript

## max. code width ============================================================

## Title: Arc Institute in vivo pertubation screens
## Author: DKatzlberger
## Date: 21.01.2026

## log ========================================================================

showf <- function(fmt, ..., type = c(
  "message", "print", "cat", "warning", "stop"
)) {
  type <- match.arg(type)
  txt <- sprintf(fmt, ...)
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
    cg <- readLines("/proc/self/cgroup", warn = FALSE)
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
        renv::install("BiocVersion", prompt = FALSE)
      }

      # BiocManager repositories 
      suppressMessages(options(repos = BiocManager::repositories()))
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
assign_guide_by_max_count <- function(
  obj,
  assay = "CRISPR",
  ntc_pattern = "NTC|NonTargeting|BRDN",
  verbose = TRUE
) {

  # Extract guide counts
  counts <- obj@assays[[assay]]$counts

  # Assign guide
  obj$guide <- apply(counts, 2, function(x) {
    rownames(counts)[which.max(x)]
  })

  # Assign genotype
  obj$genotype <- ifelse(
    grepl(ntc_pattern, obj$guide, ignore.case = TRUE),
    "NTC",
    sub("-.*$", "", obj$guide)
  )

  # Create on_guide column for NTC comparison
  obj$on_guide <- ifelse(
    grepl(ntc_pattern, obj$guide, ignore.case = TRUE),
    "NTC",
    obj$guide
  )

  if (verbose) {
    showf("Assigned guides to %s cells", ncol(obj), type = "message")
    showf("- %s unique knockouts/genotypes", length(unique(obj$genotype)), type = "message")
  }

  # Return
  obj
}
