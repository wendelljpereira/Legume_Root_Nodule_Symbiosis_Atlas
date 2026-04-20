#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Seurat))

script_flag <- "--file="
script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub(script_flag, "", script_arg[grep(script_flag, script_arg)][1])
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(project_root)

source("scripts/atlas_dataset_utils.R")

app_env <- new.env(parent = globalenv())
sys.source("app.R", envir = app_env)

require_slim_path <- function(path, label) {
    if (!file.exists(path)) {
        stop(sprintf("%s slim file is missing: %s", label, path), call. = FALSE)
    }

    path
}

read_prepared_within <- function(species_key, integration_method) {
    label <- sprintf("%s (%s integration)", app_env$species_registry[[species_key]]$label, integration_method)
    base_path <- app_env$species_registry[[species_key]]$within_paths[[integration_method]]
    slim_path <- app_slim_path(base_path)
    obj <- app_env$safe_read_rds(require_slim_path(slim_path, label), label)
    app_env$prepare_within_object(
        obj,
        override_id = paste(species_key, integration_method, sep = "_")
    )
}

read_prepared_cross <- function(cross_key) {
    cfg <- app_env$cross_integration_registry[[cross_key]]
    label <- sprintf("%s cross-species integration", cfg$label)
    slim_path <- cfg$slim_path %||% app_slim_path(cfg$path)
    obj <- app_env$safe_read_rds(require_slim_path(slim_path, label), label)
    app_env$prepare_cross_object(obj, cross_key)
}

assert_data_layer <- function(obj, label) {
    assay_obj <- obj[[DefaultAssay(obj)]]
    layers <- default_assay_layers(assay_obj)

    if (!("data" %in% layers)) {
        stop(sprintf("%s is missing a data layer after load.", label), call. = FALSE)
    }
}

assert_within_reductions <- function(obj, label) {
    reductions <- Reductions(obj)

    if (!("umap" %in% reductions)) {
        stop(sprintf("%s is missing the 2D UMAP reduction.", label), call. = FALSE)
    }

    has_3d_or_fallback <- any(c(app_env$within_three_d_reduction_name, "umap_3d", "umap3d", "pca") %in% reductions)

    if (!has_3d_or_fallback) {
        stop(sprintf("%s is missing both 3D UMAP and PCA fallback reductions.", label), call. = FALSE)
    }
}

assert_cross_reductions <- function(obj, cross_key) {
    cfg <- app_env$cross_integration_registry[[cross_key]]
    label <- sprintf("%s cross-species integration", cfg$label)
    reductions <- Reductions(obj)

    if (!("umap" %in% reductions)) {
        stop(sprintf("%s is missing the 2D UMAP reduction.", label), call. = FALSE)
    }

    has_3d <- any(c("umap3d", "umap_3d", app_env$within_three_d_reduction_name) %in% reductions) ||
        file.exists(cfg$umap3d_path %||% "")

    if (!has_3d) {
        stop(sprintf("%s has neither an embedded 3D UMAP reduction nor %s.", label, cfg$umap3d_path %||% "<missing path>"), call. = FALSE)
    }
}

report_dataset <- function(label, path, obj) {
    cat(
        sprintf(
            "OK\t%s\tcells=%s\tgenes=%s\tpath=%s\n",
            label,
            format(ncol(obj), big.mark = ","),
            format(nrow(obj), big.mark = ","),
            path
        )
    )
}

within_specs <- expand.grid(
    species_key = app_env$within_species_keys,
    integration_method = unname(app_env$integration_choices),
    stringsAsFactors = FALSE
)

for (i in seq_len(nrow(within_specs))) {
    species_key <- within_specs$species_key[[i]]
    integration_method <- within_specs$integration_method[[i]]
    label <- sprintf("%s (%s integration)", app_env$species_registry[[species_key]]$label, integration_method)
    slim_path <- app_slim_path(app_env$species_registry[[species_key]]$within_paths[[integration_method]])
    obj <- read_prepared_within(species_key, integration_method)

    if (ncol(obj) <= 0 || nrow(obj) <= 0) {
        stop(sprintf("%s loaded but has invalid dimensions.", label), call. = FALSE)
    }

    assert_data_layer(obj, label)
    assert_within_reductions(obj, label)
    report_dataset(label, slim_path, obj)
}

for (cross_key in app_env$cross_integration_keys) {
    cfg <- app_env$cross_integration_registry[[cross_key]]
    label <- sprintf("%s cross-species integration", cfg$label)
    slim_path <- cfg$slim_path %||% app_slim_path(cfg$path)
    obj <- read_prepared_cross(cross_key)

    if (ncol(obj) <= 0 || nrow(obj) <= 0) {
        stop(sprintf("%s loaded but has invalid dimensions.", label), call. = FALSE)
    }

    assert_data_layer(obj, label)
    assert_cross_reductions(obj, cross_key)
    report_dataset(label, slim_path, obj)
}

cat("ALL DATASETS OK\n")
