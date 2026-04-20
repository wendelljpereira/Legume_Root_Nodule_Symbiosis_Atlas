#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Seurat))

script_flag <- "--file="
script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub(script_flag, "", script_arg[grep(script_flag, script_arg)][1])
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(project_root)

source("scripts/atlas_dataset_utils.R")

within_meta_columns <- c(
    "orig.ident",
    "sample_name",
    "sample",
    "Sample",
    "Group",
    "condition",
    "batch",
    "study",
    "time_point",
    "integrated_dataset",
    "integration_method",
    "Rank_1st",
    "Rank_2nd",
    "Rank_3rd",
    "Rank_4th",
    "Rank_5th"
)

within_specs <- list(
    medicago_combat = list(
        label = "Medicago truncatula (ComBat_BBKNN)",
        path = "within_species_integrated_datasets/ComBat_BBKNN/M_truncatula_clustered_dataset.rds"
    ),
    glycine_combat = list(
        label = "Glycine max (ComBat_BBKNN)",
        path = "within_species_integrated_datasets/ComBat_BBKNN/G_max_clustered_dataset.rds"
    ),
    lotus_combat = list(
        label = "Lotus japonicus (ComBat_BBKNN)",
        path = "within_species_integrated_datasets/ComBat_BBKNN/L_japonicus_clustered_dataset.rds"
    ),
    medicago_seurat = list(
        label = "Medicago truncatula (Seurat)",
        path = "within_species_integrated_datasets/Seurat/M_truncatula_clustered_dataset.rds"
    ),
    glycine_seurat = list(
        label = "Glycine max (Seurat)",
        path = "within_species_integrated_datasets/Seurat/G_max_clustered_dataset.rds"
    ),
    lotus_seurat = list(
        label = "Lotus japonicus (Seurat)",
        path = "within_species_integrated_datasets/Seurat/L_japonicus_clustered_dataset.rds"
    )
)

build_within_slim <- function(spec) {
    out_path <- app_slim_path(spec$path)
    cat("\n== Building slim within-species dataset:", spec$label, "==\n")
    obj <- readRDS(spec$path)
    print_object_audit(paste(spec$label, "before"), spec$path, obj)

    keep_reductions <- slim_reduction_names(obj, prefer_pca_fallback = FALSE)

    slim_obj <- slim_seurat_object(
        obj = obj,
        keep_meta_cols = within_meta_columns,
        keep_reductions = keep_reductions,
        keep_misc = FALSE
    )

    saveRDS(slim_obj, out_path, compress = "xz")
    print_object_audit(paste(spec$label, "after"), out_path, slim_obj)

    original_size <- file.info(spec$path)$size
    slim_size <- file.info(out_path)$size
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

invisible(lapply(within_specs, build_within_slim))
