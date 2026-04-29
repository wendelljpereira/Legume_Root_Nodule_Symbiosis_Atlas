#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(Seurat)
    library(uwot)
})

script_flag <- "--file="
script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub(script_flag, "", script_arg[grep(script_flag, script_arg)][1])
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(project_root)

source("scripts/atlas_dataset_utils.R")

within_three_d_reduction_name <- "umap3d"

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

dataset_specs <- data.frame(
    species = c(
        "medicago", "glycine", "lotus",
        "medicago", "glycine", "lotus",
        "medicago", "glycine", "lotus"
    ),
    integration_method = c(
        "ComBat_BBKNN", "ComBat_BBKNN", "ComBat_BBKNN",
        "Seurat", "Seurat", "Seurat",
        "Saturn", "Saturn", "Saturn"
    ),
    path = c(
        "within_species_integrated_datasets/ComBat_BBKNN/M_truncatula_clustered_dataset.rds",
        "within_species_integrated_datasets/ComBat_BBKNN/G_max_clustered_dataset.rds",
        "within_species_integrated_datasets/ComBat_BBKNN/L_japonicus_clustered_dataset.rds",
        "within_species_integrated_datasets/Seurat/M_truncatula_clustered_dataset.rds",
        "within_species_integrated_datasets/Seurat/G_max_clustered_dataset.rds",
        "within_species_integrated_datasets/Seurat/L_japonicus_clustered_dataset.rds",
        "within_species_integrated_datasets/Saturn/M_truncatula_clustered_dataset.rds",
        "within_species_integrated_datasets/Saturn/G_max_clustered_dataset.rds",
        "within_species_integrated_datasets/Saturn/L_japonicus_clustered_dataset.rds"
    ),
    stringsAsFactors = FALSE
)

compute_umap3d_matrix <- function(obj, dims = 30L, seed = 1234L, reduction = NULL) {
    reduction_candidates <- unique(c(
        reduction,
        "pca",
        "saturn_latent_pca",
        "integrated_pca",
        "harmony",
        "saturn_integration"
    ))
    reduction_candidates <- reduction_candidates[!is.na(reduction_candidates) & nzchar(reduction_candidates)]
    reduction_name <- reduction_candidates[reduction_candidates %in% Reductions(obj)][1]

    if (is.na(reduction_name) || !length(reduction_name)) {
        stop("PCA coordinates are required to compute a 3D UMAP.", call. = FALSE)
    }

    pca_embeddings <- Embeddings(obj, reduction_name)
    dims_use <- seq_len(min(as.integer(dims), ncol(pca_embeddings)))

    if (length(dims_use) < 3L) {
        stop("At least three PCA dimensions are required to compute a 3D UMAP.", call. = FALSE)
    }

    message("  Computing 3D UMAP from reduction '", reduction_name, "' using ", length(dims_use), " dimensions.")
    set.seed(seed)
    embedding <- uwot::umap(
        X = pca_embeddings[, dims_use, drop = FALSE],
        n_components = 3L,
        n_neighbors = 30,
        min_dist = 0.3,
        metric = "cosine",
        init = "spectral",
        ret_model = FALSE,
        verbose = TRUE,
        n_threads = max(1L, min(4L, parallel::detectCores(logical = TRUE) - 1L))
    )

    rownames(embedding) <- rownames(pca_embeddings)
    colnames(embedding) <- paste0("UMAP3D_", seq_len(ncol(embedding)))
    embedding
}

write_updated_object <- function(obj, path) {
    tmp_path <- paste0(path, ".tmp")
    saveRDS(obj, tmp_path, compress = "xz")
    if (!file.rename(tmp_path, path)) {
        unlink(tmp_path)
        stop("Failed to replace ", path, " with the updated object.", call. = FALSE)
    }
}

build_or_read_slim <- function(full_obj, slim_path) {
    if (file.exists(slim_path)) {
        return(readRDS(slim_path))
    }

    slim_seurat_object(
        obj = full_obj,
        keep_meta_cols = within_meta_columns,
        keep_reductions = slim_reduction_names(full_obj, prefer_pca_fallback = FALSE),
        keep_misc = FALSE
    )
}

add_umap3d_to_slim <- function(spec, overwrite = FALSE) {
    slim_path <- app_slim_path(spec$path)
    label <- sprintf("%s %s", spec$species, spec$integration_method)
    message("\n== Processing ", label, " ==")
    message("  Full source: ", spec$path)
    message("  Slim target: ", slim_path)

    full_obj <- readRDS(spec$path)
    slim_obj <- build_or_read_slim(full_obj, slim_path)

    if (!overwrite &&
        within_three_d_reduction_name %in% Reductions(slim_obj) &&
        ncol(Embeddings(slim_obj, within_three_d_reduction_name)) >= 3L) {
        message("  Skipping: ", within_three_d_reduction_name, " already present in the slim file.")
        return(invisible(FALSE))
    }

    start_time <- Sys.time()
    full_embedding <- if (within_three_d_reduction_name %in% Reductions(full_obj) &&
        ncol(Embeddings(full_obj, within_three_d_reduction_name)) >= 3L &&
        !overwrite) {
        message("  Reusing ", within_three_d_reduction_name, " from the full object.")
        Embeddings(full_obj, within_three_d_reduction_name)[, seq_len(3), drop = FALSE]
    } else {
        compute_umap3d_matrix(full_obj)
    }

    missing_cells <- setdiff(colnames(slim_obj), rownames(full_embedding))
    if (length(missing_cells)) {
        stop(
            "The computed 3D UMAP is missing ",
            length(missing_cells),
            " slim-file cells. First missing cell: ",
            missing_cells[[1]],
            call. = FALSE
        )
    }

    embedding <- full_embedding[colnames(slim_obj), seq_len(3), drop = FALSE]
    slim_obj[[within_three_d_reduction_name]] <- CreateDimReducObject(
        embeddings = embedding,
        key = "UMAP3D_",
        assay = DefaultAssay(slim_obj)
    )

    write_updated_object(slim_obj, slim_path)
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    message("  Saved ", within_three_d_reduction_name, " to slim file in ", round(as.numeric(elapsed), 2), " min.")
    invisible(TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
targets <- setdiff(args, "--overwrite")

specs_to_run <- dataset_specs

if (length(targets)) {
    target_values <- tolower(targets)
    specs_to_run <- specs_to_run[
        tolower(specs_to_run$species) %in% target_values |
            tolower(specs_to_run$integration_method) %in% target_values |
            tolower(basename(specs_to_run$path)) %in% target_values,
        ,
        drop = FALSE
    ]
}

if (!nrow(specs_to_run)) {
    stop("No matching datasets found for the supplied arguments.", call. = FALSE)
}

for (idx in seq_len(nrow(specs_to_run))) {
    add_umap3d_to_slim(spec = specs_to_run[idx, , drop = FALSE], overwrite = overwrite)
}
