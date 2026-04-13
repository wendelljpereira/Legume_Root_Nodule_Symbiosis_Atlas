#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Seurat))

script_flag <- "--file="
script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub(script_flag, "", script_arg[grep(script_flag, script_arg)][1])
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(project_root)

within_species_keys <- c("medicago", "glycine", "lotus")

species_registry <- list(
    medicago = list(
        within_paths = c(
            ComBat_BBKNN = "within_species_integrated_datasets/ComBat_BBKNN/M_truncatula_clustered_dataset.rds",
            Seurat = "within_species_integrated_datasets/Seurat/M_truncatula_clustered_dataset.rds"
        )
    ),
    glycine = list(
        within_paths = c(
            ComBat_BBKNN = "within_species_integrated_datasets/ComBat_BBKNN/G_max_clustered_dataset.rds",
            Seurat = "within_species_integrated_datasets/Seurat/G_max_clustered_dataset.rds"
        )
    ),
    lotus = list(
        within_paths = c(
            ComBat_BBKNN = "within_species_integrated_datasets/ComBat_BBKNN/L_japonicus_clustered_dataset.rds",
            Seurat = "within_species_integrated_datasets/Seurat/L_japonicus_clustered_dataset.rds"
        )
    )
)

integration_choices <- c(
    "ComBat/BBKNN" = "ComBat_BBKNN",
    "Seurat" = "Seurat"
)

cross_species_path <- "cross_species_integrated_datasets/Camex/clustered_dataset.rds"
atlas_summary_path <- "metadata/atlas_summary.tsv"

compact_value_list <- function(values, limit = 4, empty_label = "NA") {
    values <- unique(trimws(as.character(values)))
    values <- values[!is.na(values) & nzchar(values)]

    if (!length(values)) {
        return(empty_label)
    }

    if (length(values) <= limit) {
        return(paste(values, collapse = ", "))
    }

    paste0(
        paste(values[seq_len(limit)], collapse = ", "),
        " (+", length(values) - limit, " more)"
    )
}

pick_first_existing_col <- function(df, candidates) {
    match_col <- intersect(candidates, colnames(df))[1]
    if (is.na(match_col)) NA_character_ else match_col
}

pick_sample_column <- function(df) {
    pick_first_existing_col(df, c("sample_name", "sample", "Sample", "orig.ident"))
}

build_within_summary_row <- function(species_key, integration_method, obj) {
    md <- obj@meta.data
    sample_col <- pick_sample_column(md)
    sample_values <- if (!is.na(sample_col)) md[[sample_col]] else character(0)

    data.frame(
        dataset_scope = "within",
        species_key = species_key,
        integration_method = integration_method,
        integration_label = names(integration_choices)[match(integration_method, integration_choices)],
        cells = ncol(obj),
        genes = nrow(obj),
        sample_n = length(unique(sample_values[!is.na(sample_values) & nzchar(sample_values)])),
        group_label = NA_character_,
        group_n = NA_integer_,
        group_preview = NA_character_,
        species_n = NA_integer_,
        species_preview = NA_character_,
        time_n = NA_integer_,
        time_preview = NA_character_,
        stringsAsFactors = FALSE
    )
}

build_cross_summary_row <- function(obj) {
    md <- obj@meta.data
    sample_col <- pick_sample_column(md)
    sample_values <- if (!is.na(sample_col)) md[[sample_col]] else character(0)

    data.frame(
        dataset_scope = "cross",
        species_key = "cross",
        integration_method = NA_character_,
        integration_label = "Camex",
        cells = ncol(obj),
        genes = nrow(obj),
        sample_n = length(unique(sample_values[!is.na(sample_values) & nzchar(sample_values)])),
        group_label = NA_character_,
        group_n = NA_integer_,
        group_preview = NA_character_,
        species_n = NA_integer_,
        species_preview = NA_character_,
        time_n = NA_integer_,
        time_preview = NA_character_,
        stringsAsFactors = FALSE
    )
}

within_rows <- do.call(
    rbind,
    lapply(within_species_keys, function(species_key) {
        do.call(
            rbind,
            lapply(unname(integration_choices), function(integration_method) {
                obj <- readRDS(species_registry[[species_key]]$within_paths[[integration_method]])
                build_within_summary_row(species_key, integration_method, obj)
            })
        )
    })
)

cross_row <- build_cross_summary_row(readRDS(cross_species_path))
summary_df <- rbind(within_rows, cross_row)

dir.create(dirname(atlas_summary_path), recursive = TRUE, showWarnings = FALSE)
write.table(
    summary_df,
    file = atlas_summary_path,
    sep = "\t",
    row.names = FALSE,
    quote = TRUE,
    na = ""
)

cat("Wrote atlas summary cache to", atlas_summary_path, "\n")
