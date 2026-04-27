#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(Seurat)
})

script_flag <- "--file="
script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub(script_flag, "", script_arg[grep(script_flag, script_arg)][1])
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(project_root)
source("scripts/atlas_dataset_utils.R")
source("R/atlas_annotations.R")

annotations_dir <- atlas_cluster_annotations_dir()
legacy_overrides_dir <- atlas_legacy_celltype_overrides_dir()

within_species_keys <- c("medicago", "glycine", "lotus")
integration_methods <- c("ComBat_BBKNN", "Seurat", "Saturn")

species_registry <- list(
    medicago = list(
        within_paths = c(
            ComBat_BBKNN = "within_species_integrated_datasets/ComBat_BBKNN/M_truncatula_clustered_dataset.rds",
            Seurat = "within_species_integrated_datasets/Seurat/M_truncatula_clustered_dataset.rds",
            Saturn = "within_species_integrated_datasets/Saturn/M_truncatula_clustered_dataset.rds"
        )
    ),
    glycine = list(
        within_paths = c(
            ComBat_BBKNN = "within_species_integrated_datasets/ComBat_BBKNN/G_max_clustered_dataset.rds",
            Seurat = "within_species_integrated_datasets/Seurat/G_max_clustered_dataset.rds",
            Saturn = "within_species_integrated_datasets/Saturn/G_max_clustered_dataset.rds"
        )
    ),
    lotus = list(
        within_paths = c(
            ComBat_BBKNN = "within_species_integrated_datasets/ComBat_BBKNN/L_japonicus_clustered_dataset.rds",
            Seurat = "within_species_integrated_datasets/Seurat/L_japonicus_clustered_dataset.rds",
            Saturn = "within_species_integrated_datasets/Saturn/L_japonicus_clustered_dataset.rds"
        )
    )
)

cross_integration_registry <- list(
    camex = list(
        path = "app_ready_integration/camex/clustered_dataset.rds",
        slim_path = "app_ready_integration/camex/clustered_dataset_app_slim.rds"
    ),
    saturn = list(
        path = "app_ready_integration/saturn/clustered_dataset.rds",
        slim_path = "app_ready_integration/saturn/clustered_dataset_app_slim.rds"
    )
)

dataset_specs <- c(
    lapply(within_species_keys, function(species_key) {
        lapply(integration_methods, function(integration_method) {
            raw_path <- species_registry[[species_key]]$within_paths[[integration_method]]
            list(
                dataset_key = paste(species_key, integration_method, sep = "_"),
                path = pick_first_existing_path(c(app_slim_path(raw_path), raw_path))
            )
        })
    }) |> unlist(recursive = FALSE),
    lapply(names(cross_integration_registry), function(cross_key) {
        cfg <- cross_integration_registry[[cross_key]]
        list(
            dataset_key = cross_key,
            path = pick_first_existing_path(c(cfg$slim_path, cfg$path))
        )
    })
)

for (spec in dataset_specs) {
    dataset_key <- spec$dataset_key
    path <- spec$path
    annotation_path <- cluster_annotation_file_path(dataset_key, annotations_dir = annotations_dir)

    if (!file.exists(path)) {
        stop(sprintf("Dataset '%s' is missing: %s", dataset_key, path), call. = FALSE)
    }

    if (!file.exists(annotation_path)) {
        stop(sprintf("Cluster annotation file is missing for dataset '%s': %s", dataset_key, annotation_path), call. = FALSE)
    }

    obj <- readRDS(path)
    obj <- apply_cluster_annotation_overlay(
        obj,
        dataset_key = dataset_key,
        annotations_dir = annotations_dir,
        legacy_overrides_dir = legacy_overrides_dir
    )

    annotation_tbl <- read_cluster_annotation_file(dataset_key, annotations_dir = annotations_dir, required = TRUE)
    cat(sprintf(
        "OK\t%s\t%s rows\t%s\n",
        dataset_key,
        format(nrow(annotation_tbl), big.mark = ","),
        annotation_path
    ))
}

cat("\nAll cluster annotation files validated successfully.\n")
