#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Seurat))

script_flag <- "--file="
script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub(script_flag, "", script_arg[grep(script_flag, script_arg)][1])
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(project_root)
source("scripts/atlas_dataset_utils.R")

within_species_keys <- c("medicago", "glycine", "lotus")

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

integration_choices <- c(
    "ComBat/BBKNN" = "ComBat_BBKNN",
    "Seurat" = "Seurat",
    "SATURN" = "Saturn"
)

cross_integration_keys <- c("camex", "saturn")

cross_integration_registry <- list(
    camex = list(
        label = "Camex",
        path = "app_ready_integration/camex/clustered_dataset.rds",
        slim_path = "app_ready_integration/camex/clustered_dataset_app_slim.rds"
    ),
    saturn = list(
        label = "SATURN",
        path = "app_ready_integration/saturn/clustered_dataset.rds",
        slim_path = "app_ready_integration/saturn/clustered_dataset_app_slim.rds"
    )
)

atlas_summary_path <- "metadata/atlas_summary.tsv"

get_cross_dataset_path <- function(cross_key) {
    cfg <- cross_integration_registry[[cross_key]]
    pick_first_existing_path(c(cfg$slim_path, cfg$path))
}

get_within_dataset_path <- function(species_key, integration_method) {
    path <- species_registry[[species_key]]$within_paths[[integration_method]]
    pick_first_existing_path(c(app_slim_path(path), path))
}

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

species_labels <- c(
    medicago = "Medicago truncatula",
    glycine = "Glycine max",
    lotus = "Lotus japonicus"
)

build_metadata_axis_summary <- function(values, label, limit = 4) {
    clean_values <- unique(trimws(as.character(values)))
    clean_values <- clean_values[!is.na(clean_values) & nzchar(clean_values)]

    list(
        label = label,
        n = length(clean_values),
        preview = if (length(clean_values)) compact_value_list(clean_values, limit = limit) else NA_character_
    )
}

pick_metadata_axis_summary <- function(df, candidates, limit = 4) {
    for (candidate in candidates) {
        column_name <- pick_first_existing_col(df, candidate$columns)
        if (!is.na(column_name)) {
            return(build_metadata_axis_summary(df[[column_name]], candidate$label, limit = limit))
        }
    }

    list(label = NA_character_, n = NA_integer_, preview = NA_character_)
}

build_within_summary_row <- function(species_key, integration_method, obj) {
    md <- obj@meta.data
    sample_col <- pick_sample_column(md)
    sample_values <- if (!is.na(sample_col)) as.character(md[[sample_col]]) else character(0)
    group_summary <- pick_metadata_axis_summary(
        md,
        list(
            list(columns = c("Group", "condition"), label = "conditions"),
            list(columns = c("study"), label = "studies")
        )
    )
    time_summary <- pick_metadata_axis_summary(
        md,
        list(
            list(columns = c("time_point", "time", "Time"), label = "time points")
        )
    )
    species_summary <- build_metadata_axis_summary(species_labels[[species_key]], label = "species", limit = 1)

    data.frame(
        dataset_scope = "within",
        species_key = species_key,
        integration_method = integration_method,
        integration_label = names(integration_choices)[match(integration_method, integration_choices)],
        cells = ncol(obj),
        genes = nrow(obj),
        sample_n = length(unique(sample_values[!is.na(sample_values) & nzchar(sample_values)])),
        group_label = group_summary$label,
        group_n = group_summary$n,
        group_preview = group_summary$preview,
        species_n = species_summary$n,
        species_preview = species_summary$preview,
        time_n = time_summary$n,
        time_preview = time_summary$preview,
        stringsAsFactors = FALSE
    )
}

build_cross_summary_row <- function(cross_key, obj) {
    md <- obj@meta.data
    sample_col <- pick_sample_column(md)
    sample_values <- if (!is.na(sample_col)) as.character(md[[sample_col]]) else character(0)
    group_summary <- pick_metadata_axis_summary(
        md,
        list(
            list(columns = c("condition"), label = "conditions"),
            list(columns = c("cell_class"), label = "cell classes"),
            list(columns = c("saturn_ref_label", "saturn_label"), label = "reference labels"),
            list(columns = c("study"), label = "studies")
        )
    )
    species_summary <- pick_metadata_axis_summary(
        md,
        list(
            list(columns = c("species"), label = "species")
        )
    )
    time_summary <- pick_metadata_axis_summary(
        md,
        list(
            list(columns = c("time_point", "time", "Time"), label = "time points")
        )
    )

    data.frame(
        dataset_scope = "cross",
        species_key = cross_key,
        integration_method = cross_key,
        integration_label = cross_integration_registry[[cross_key]]$label,
        cells = ncol(obj),
        genes = nrow(obj),
        sample_n = length(unique(sample_values[!is.na(sample_values) & nzchar(sample_values)])),
        group_label = group_summary$label,
        group_n = group_summary$n,
        group_preview = group_summary$preview,
        species_n = species_summary$n,
        species_preview = species_summary$preview,
        time_n = time_summary$n,
        time_preview = time_summary$preview,
        stringsAsFactors = FALSE
    )
}

within_rows <- do.call(
    rbind,
    lapply(within_species_keys, function(species_key) {
        do.call(
            rbind,
            lapply(unname(integration_choices), function(integration_method) {
                obj <- readRDS(get_within_dataset_path(species_key, integration_method))
                build_within_summary_row(species_key, integration_method, obj)
            })
        )
    })
)

cross_rows <- do.call(
    rbind,
    lapply(cross_integration_keys, function(cross_key) {
        build_cross_summary_row(
            cross_key = cross_key,
            obj = readRDS(get_cross_dataset_path(cross_key))
        )
    })
)

summary_df <- rbind(within_rows, cross_rows)

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
