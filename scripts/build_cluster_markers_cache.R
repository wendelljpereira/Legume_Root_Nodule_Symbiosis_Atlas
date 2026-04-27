#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(tibble)
    library(purrr)
})

script_flag <- "--file="
script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub(script_flag, "", script_arg[grep(script_flag, script_arg)][1])
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(project_root)
source("scripts/atlas_dataset_utils.R")

within_species_keys <- c("medicago", "glycine", "lotus")
integration_methods <- c("ComBat_BBKNN", "Seurat", "Saturn")
cluster_markers_dir <- "metadata/cluster_markers"

species_registry <- list(
    medicago = list(
        label = "Medicago truncatula",
        within_paths = c(
            ComBat_BBKNN = "within_species_integrated_datasets/ComBat_BBKNN/M_truncatula_clustered_dataset.rds",
            Seurat = "within_species_integrated_datasets/Seurat/M_truncatula_clustered_dataset.rds",
            Saturn = "within_species_integrated_datasets/Saturn/M_truncatula_clustered_dataset.rds"
        )
    ),
    glycine = list(
        label = "Glycine max",
        within_paths = c(
            ComBat_BBKNN = "within_species_integrated_datasets/ComBat_BBKNN/G_max_clustered_dataset.rds",
            Seurat = "within_species_integrated_datasets/Seurat/G_max_clustered_dataset.rds",
            Saturn = "within_species_integrated_datasets/Saturn/G_max_clustered_dataset.rds"
        )
    ),
    lotus = list(
        label = "Lotus japonicus",
        within_paths = c(
            ComBat_BBKNN = "within_species_integrated_datasets/ComBat_BBKNN/L_japonicus_clustered_dataset.rds",
            Seurat = "within_species_integrated_datasets/Seurat/L_japonicus_clustered_dataset.rds",
            Saturn = "within_species_integrated_datasets/Saturn/L_japonicus_clustered_dataset.rds"
        )
    )
)

cross_integration_registry <- list(
    camex = list(
        label = "Camex",
        path = "app_ready_integration/camex/clustered_dataset.rds"
    ),
    saturn = list(
        label = "SATURN",
        path = "app_ready_integration/saturn/clustered_dataset.rds"
    )
)

cluster_markers_cache_path <- function(dataset_key, top_n = FALSE) {
    suffix <- if (isTRUE(top_n)) "_top10" else ""
    file.path(cluster_markers_dir, paste0(dataset_key, suffix, ".tsv"))
}

normalize_marker_columns <- function(marker_tbl) {
    if (!nrow(marker_tbl)) {
        return(tibble(
            cluster = character(0),
            gene = character(0),
            avg_log2FC = numeric(0),
            pct.1 = numeric(0),
            pct.2 = numeric(0),
            p_val_adj = numeric(0)
        ))
    }

    gene_values <- marker_tbl$gene %||% rownames(marker_tbl)
    logfc_col <- c("avg_log2FC", "avg_logFC", "avg_log10FC")[c("avg_log2FC", "avg_logFC", "avg_log10FC") %in% colnames(marker_tbl)][1]

    if (is.na(logfc_col) || !length(logfc_col)) {
        stop("FindAllMarkers output is missing an average log-fold-change column.", call. = FALSE)
    }

    marker_tbl %>%
        mutate(
            cluster = as.character(cluster),
            gene = as.character(gene_values),
            avg_log2FC = as.numeric(.data[[logfc_col]]),
            pct.1 = as.numeric(.data[["pct.1"]]),
            pct.2 = as.numeric(.data[["pct.2"]]),
            p_val_adj = as.numeric(.data[["p_val_adj"]])
        ) %>%
        select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
        filter(!is.na(cluster) & nzchar(cluster), !is.na(gene) & nzchar(gene))
}

top_markers_per_cluster <- function(marker_tbl, n = 10L) {
    if (!nrow(marker_tbl)) {
        return(marker_tbl)
    }

    marker_tbl %>%
        arrange(cluster, p_val_adj, desc(avg_log2FC), desc(pct.1), gene) %>%
        group_by(cluster) %>%
        slice_head(n = n) %>%
        ungroup()
}

write_marker_table <- function(marker_tbl, path) {
    write.table(
        marker_tbl,
        file = path,
        sep = "\t",
        row.names = FALSE,
        quote = TRUE,
        na = ""
    )
}

dataset_registry <- bind_rows(
    bind_rows(lapply(within_species_keys, function(species_key) {
        bind_rows(lapply(integration_methods, function(integration_method) {
            raw_path <- species_registry[[species_key]]$within_paths[[integration_method]]
            tibble(
                dataset_key = paste(species_key, integration_method, sep = "_"),
                dataset_label = sprintf("%s (%s)", species_registry[[species_key]]$label, integration_method),
                dataset_path = pick_first_existing_path(c(app_slim_path(raw_path), raw_path))
            )
        }))
    })),
    bind_rows(lapply(names(cross_integration_registry), function(cross_key) {
        raw_path <- cross_integration_registry[[cross_key]]$path
        tibble(
            dataset_key = cross_key,
            dataset_label = cross_integration_registry[[cross_key]]$label,
            dataset_path = pick_first_existing_path(c(app_slim_path(raw_path), raw_path))
        )
    }))
)

dir.create(cluster_markers_dir, recursive = TRUE, showWarnings = FALSE)

walk(seq_len(nrow(dataset_registry)), function(row_idx) {
    dataset_row <- dataset_registry[row_idx, ]
    dataset_key <- dataset_row$dataset_key[[1]]
    dataset_label <- dataset_row$dataset_label[[1]]
    dataset_path <- dataset_row$dataset_path[[1]]

    cat(sprintf("\n[%s] Loading %s\n", dataset_key, dataset_path))
    obj <- readRDS(dataset_path)
    cluster_ids <- as.character(Idents(obj))
    cluster_ids <- cluster_ids[!is.na(cluster_ids) & nzchar(cluster_ids)]
    cluster_n <- length(unique(cluster_ids))

    cat(sprintf(
        "[%s] %s | %s cells | %s genes | %s clusters\n",
        dataset_key,
        dataset_label,
        format(ncol(obj), big.mark = ","),
        format(nrow(obj), big.mark = ","),
        format(cluster_n, big.mark = ",")
    ))

    if (cluster_n < 2L) {
        stop(sprintf("[%s] Need at least two identity classes to compute markers.", dataset_key), call. = FALSE)
    }

    set.seed(123L)
    raw_markers <- suppressMessages(FindAllMarkers(
        object = obj,
        only.pos = TRUE,
        min.pct = 0.25,
        logfc.threshold = 0.25,
        max.cells.per.ident = 2000,
        verbose = FALSE
    ))

    marker_tbl <- normalize_marker_columns(raw_markers)
    top10_tbl <- top_markers_per_cluster(marker_tbl, n = 10L)

    full_path <- cluster_markers_cache_path(dataset_key, top_n = FALSE)
    top10_path <- cluster_markers_cache_path(dataset_key, top_n = TRUE)

    write_marker_table(marker_tbl, full_path)
    write_marker_table(top10_tbl, top10_path)

    cat(sprintf(
        "[%s] Wrote %s (%s rows) and %s (%s rows)\n",
        dataset_key,
        full_path,
        format(nrow(marker_tbl), big.mark = ","),
        top10_path,
        format(nrow(top10_tbl), big.mark = ",")
    ))
})

cat("\nAll cluster marker caches built successfully.\n")
