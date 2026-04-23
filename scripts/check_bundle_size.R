#!/usr/bin/env Rscript

script_flag <- "--file="
script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub(script_flag, "", script_arg[grep(script_flag, script_arg)][1])
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(project_root)

args <- commandArgs(trailingOnly = TRUE)
warn_only <- "--warn-only" %in% args
limit_arg <- args[startsWith(args, "--limit-gb=")][1]
limit_gb <- if (length(limit_arg) && !is.na(limit_arg)) {
    suppressWarnings(as.numeric(sub("^--limit-gb=", "", limit_arg)))
} else {
    suppressWarnings(as.numeric(Sys.getenv("ATLAS_BUNDLE_SIZE_LIMIT_GB", "5")))
}

if (is.na(limit_gb) || limit_gb <= 0) {
    stop("Bundle size limit must be a positive number of GB.", call. = FALSE)
}

bundle_paths <- c(
    "app.R",
    "R",
    "www",
    "scripts",
    "metadata",
    "annotations",
    "orthogroups",
    "within_species_integrated_datasets",
    "app_ready_integration",
    "DESCRIPTION",
    "LICENSE",
    "renv.lock"
)

path_size <- function(path) {
    if (!file.exists(path)) return(0)
    files <- if (dir.exists(path)) {
        list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    } else {
        path
    }
    files <- files[file.exists(files) & !dir.exists(files)]
    if (!length(files)) return(0)
    sum(file.info(files)$size, na.rm = TRUE)
}

sizes <- data.frame(
    path = bundle_paths,
    bytes = vapply(bundle_paths, path_size, numeric(1)),
    stringsAsFactors = FALSE
)
sizes <- sizes[sizes$bytes > 0, , drop = FALSE]
sizes$gb <- sizes$bytes / 1024^3

total_gb <- sum(sizes$gb)
cat("Bundle size estimate:\n")
for (i in seq_len(nrow(sizes))) {
    cat(sprintf("  %-36s %7.3f GB\n", sizes$path[[i]], sizes$gb[[i]]))
}
cat(sprintf("Total: %.3f GB | Limit: %.3f GB\n", total_gb, limit_gb))

if (total_gb > limit_gb) {
    msg <- sprintf(
        "Estimated ShinyApps.io bundle size %.3f GB exceeds configured limit %.3f GB. Use remote public-read data artifacts or reduce the bundle before deployment.",
        total_gb,
        limit_gb
    )
    if (warn_only) {
        warning(msg, call. = FALSE)
    } else {
        stop(msg, call. = FALSE)
    }
}

cat("Bundle size check passed.\n")
