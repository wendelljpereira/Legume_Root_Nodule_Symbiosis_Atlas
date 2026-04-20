#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Seurat))

script_flag <- "--file="
script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub(script_flag, "", script_arg[grep(script_flag, script_arg)][1])
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(project_root)

source("scripts/atlas_dataset_utils.R")

cross_meta_columns <- c(
    "sample_name",
    "sample",
    "Sample",
    "orig.ident",
    "species",
    "study",
    "condition",
    "Group",
    "time_point",
    "time",
    "cell_class",
    "saturn_label",
    "saturn_ref_label",
    "Rank_1st_label",
    "Rank_1st",
    "Rank_2nd",
    "Rank_3rd",
    "Rank_4th",
    "Rank_5th"
)

cross_specs <- list(
    camex = list(
        label = "Camex",
        path = "app_ready_integration/camex/clustered_dataset.rds",
        out_path = "app_ready_integration/camex/clustered_dataset_app_slim.rds"
    ),
    saturn = list(
        label = "SATURN",
        path = "app_ready_integration/saturn/clustered_dataset.rds",
        out_path = "app_ready_integration/saturn/clustered_dataset_app_slim.rds"
    )
)

build_cross_slim <- function(spec) {
    cat("\n== Building slim cross dataset:", spec$label, "==\n")
    obj <- readRDS(spec$path)
    print_object_audit(paste(spec$label, "before"), spec$path, obj)

    keep_reductions <- slim_reduction_names(obj, prefer_pca_fallback = FALSE)

    slim_obj <- slim_seurat_object(
        obj = obj,
        keep_meta_cols = cross_meta_columns,
        keep_reductions = keep_reductions,
        keep_misc = FALSE
    )

    saveRDS(slim_obj, spec$out_path, compress = "xz")
    print_object_audit(paste(spec$label, "after"), spec$out_path, slim_obj)

    original_size <- file.info(spec$path)$size
    slim_size <- file.info(spec$out_path)$size
    ratio <- slim_size / original_size
    cat(
        sprintf(
            "[%s] compression ratio: %.1f%% (%s MB -> %s MB)\n",
            spec$label,
            ratio * 100,
            round(original_size / 1024^2, 1),
            round(slim_size / 1024^2, 1)
        )
    )

    invisible(gc())
}

invisible(lapply(cross_specs, build_cross_slim))
