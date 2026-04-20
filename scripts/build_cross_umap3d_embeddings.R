#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(irlba)
    library(uwot)
})

cross_specs <- data.frame(
    cross_key = c("camex", "saturn"),
    path = c(
        "app_ready_integration/camex/clustered_dataset.rds",
        "app_ready_integration/saturn/clustered_dataset.rds"
    ),
    output_path = c(
        "app_ready_integration/camex/umap3d_embeddings.rds",
        "app_ready_integration/saturn/umap3d_embeddings.rds"
    ),
    stringsAsFactors = FALSE
)

select_variable_features <- function(mat, n_features = 2000L, min_cells = 25L) {
    detected <- Matrix::rowSums(mat > 0) >= as.integer(min_cells)
    mat <- mat[detected, , drop = FALSE]

    if (!nrow(mat)) {
        stop("No features passed the minimum detection threshold.")
    }

    means <- Matrix::rowMeans(mat)
    mat_sq <- mat
    mat_sq@x <- mat_sq@x ^ 2
    variances <- (Matrix::rowMeans(mat_sq) - means ^ 2) * ncol(mat) / max(ncol(mat) - 1, 1)
    feature_order <- order(variances, decreasing = TRUE)
    feature_count <- min(as.integer(n_features), length(feature_order))

    rownames(mat)[feature_order[seq_len(feature_count)]]
}

compute_cross_umap3d <- function(path, n_features = 2000L, pca_dims = 30L, seed = 1234L) {
    obj <- readRDS(path)
    mat <- GetAssayData(obj, layer = "data")
    features <- select_variable_features(mat, n_features = n_features)
    pca_dims_use <- min(as.integer(pca_dims), max(length(features) - 1L, 2L))

    if (pca_dims_use < 3L) {
        stop("At least three informative features are required to compute a 3D UMAP.")
    }

    message("  Selected ", length(features), " variable features")
    pca_fit <- irlba::prcomp_irlba(
        Matrix::t(mat[features, , drop = FALSE]),
        n = pca_dims_use,
        center = TRUE,
        scale. = FALSE
    )

    message("  PCA complete: ", nrow(pca_fit$x), " cells x ", ncol(pca_fit$x), " PCs")
    set.seed(seed)
    embedding <- uwot::umap(
        X = pca_fit$x[, seq_len(pca_dims_use), drop = FALSE],
        n_components = 3L,
        n_neighbors = 30,
        min_dist = 0.3,
        metric = "cosine",
        init = "spectral",
        ret_model = FALSE,
        verbose = TRUE,
        n_threads = max(1L, min(4L, parallel::detectCores(logical = TRUE) - 1L))
    )

    rownames(embedding) <- colnames(mat)
    colnames(embedding) <- paste0("UMAP3D_", seq_len(ncol(embedding)))
    embedding
}

write_embedding <- function(embedding, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp_path <- paste0(path, ".tmp")
    saveRDS(embedding, tmp_path)

    if (!file.rename(tmp_path, path)) {
        unlink(tmp_path)
        stop("Failed to replace ", path, " with the updated embedding.")
    }
}

args <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
targets <- tolower(setdiff(args, "--overwrite"))

specs_to_run <- cross_specs

if (length(targets)) {
    specs_to_run <- cross_specs[
        tolower(cross_specs$cross_key) %in% targets |
            tolower(basename(cross_specs$path)) %in% targets,
        ,
        drop = FALSE
    ]
}

if (!nrow(specs_to_run)) {
    stop("No matching cross-species datasets found for the supplied arguments.")
}

for (idx in seq_len(nrow(specs_to_run))) {
    spec <- specs_to_run[idx, , drop = FALSE]
    output_path <- spec$output_path[[1]]

    if (!overwrite && file.exists(output_path)) {
        message("Skipping ", spec$cross_key[[1]], ": embedding already present at ", output_path)
        next
    }

    message("Processing ", spec$cross_key[[1]], " from ", spec$path[[1]])
    start_time <- Sys.time()
    embedding <- compute_cross_umap3d(spec$path[[1]])
    write_embedding(embedding, output_path)
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    message("  Saved 3D embedding to ", output_path, " in ", round(as.numeric(elapsed), 2), " min")
}
