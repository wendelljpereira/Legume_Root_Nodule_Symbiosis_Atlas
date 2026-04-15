#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(Seurat)
    library(uwot)
})

within_three_d_reduction_name <- "umap3d"

dataset_specs <- data.frame(
    species = c(
        "medicago", "glycine", "lotus",
        "medicago", "glycine", "lotus"
    ),
    integration_method = c(
        "ComBat_BBKNN", "ComBat_BBKNN", "ComBat_BBKNN",
        "Seurat", "Seurat", "Seurat"
    ),
    path = c(
        "within_species_integrated_datasets/ComBat_BBKNN/M_truncatula_clustered_dataset.rds",
        "within_species_integrated_datasets/ComBat_BBKNN/G_max_clustered_dataset.rds",
        "within_species_integrated_datasets/ComBat_BBKNN/L_japonicus_clustered_dataset.rds",
        "within_species_integrated_datasets/Seurat/M_truncatula_clustered_dataset.rds",
        "within_species_integrated_datasets/Seurat/G_max_clustered_dataset.rds",
        "within_species_integrated_datasets/Seurat/L_japonicus_clustered_dataset.rds"
    ),
    stringsAsFactors = FALSE
)

compute_umap3d_matrix <- function(obj, dims = 30L, seed = 1234L) {
    if (!("pca" %in% Reductions(obj))) {
        stop("PCA coordinates are required to compute a 3D UMAP.")
    }

    pca_embeddings <- Embeddings(obj, "pca")
    dims_use <- seq_len(min(as.integer(dims), ncol(pca_embeddings)))

    if (length(dims_use) < 3L) {
        stop("At least three PCA dimensions are required to compute a 3D UMAP.")
    }

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
    saveRDS(obj, tmp_path)
    if (!file.rename(tmp_path, path)) {
        unlink(tmp_path)
        stop("Failed to replace ", path, " with the updated object.")
    }
}

add_umap3d_reduction <- function(path, overwrite = FALSE) {
    message("Processing ", path)
    obj <- readRDS(path)

    if (!overwrite &&
        within_three_d_reduction_name %in% Reductions(obj) &&
        ncol(Embeddings(obj, within_three_d_reduction_name)) >= 3L) {
        message("  Skipping: ", within_three_d_reduction_name, " already present.")
        return(invisible(FALSE))
    }

    start_time <- Sys.time()
    embedding <- compute_umap3d_matrix(obj)
    obj[[within_three_d_reduction_name]] <- CreateDimReducObject(
        embeddings = embedding,
        key = "UMAP3D_",
        assay = DefaultAssay(obj)
    )
    write_updated_object(obj, path)
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    message("  Saved ", within_three_d_reduction_name, " in ", round(as.numeric(elapsed), 2), " min")
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
    stop("No matching datasets found for the supplied arguments.")
}

for (path in specs_to_run$path) {
    add_umap3d_reduction(path = path, overwrite = overwrite)
}
