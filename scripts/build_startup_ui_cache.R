#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
})

script_flag <- "--file="
script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub(script_flag, "", script_arg[grep(script_flag, script_arg)][1])
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(project_root)
source("scripts/atlas_dataset_utils.R")
source("R/atlas_annotations.R")

within_species_keys <- c("medicago", "glycine", "lotus")
integration_methods <- c("ComBat_BBKNN", "Seurat", "Saturn")
distribution_cluster_columns <- c("Rank_1st", "Rank_2nd", "Rank_3rd", "Rank_4th", "Rank_5th", "cluster_label")
cluster_annotations_dir <- atlas_cluster_annotations_dir()
celltype_overrides_dir <- atlas_legacy_celltype_overrides_dir()

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

ui_choice_cache_path <- "metadata/ui_choice_cache.tsv"
ui_cluster_lookup_cache_path <- "metadata/ui_cluster_lookup_cache.tsv"

pick_first_existing_col <- function(df, candidates) {
    match_col <- intersect(candidates, colnames(df))[1]
    if (is.na(match_col)) NA_character_ else match_col
}

cluster_value_levels <- function(values) {
    values_chr <- unique(as.character(values))
    values_chr <- values_chr[!is.na(values_chr) & nzchar(values_chr)]

    if (!length(values_chr)) {
        return(character(0))
    }

    numeric_values <- suppressWarnings(as.numeric(values_chr))

    if (!anyNA(numeric_values)) {
        return(values_chr[order(numeric_values)])
    }

    sort(values_chr)
}

cluster_label_lookup_local <- function(obj, cluster_column = "__idents__") {
    if (!identical(cluster_column, "__idents__") && cluster_column %in% colnames(obj@meta.data)) {
        cluster_ids <- as.character(obj@meta.data[[cluster_column]])
        label_column <- paste0(cluster_column, "_label")
        cluster_labels <- if (label_column %in% colnames(obj@meta.data)) {
            as.character(obj@meta.data[[label_column]])
        } else {
            cluster_ids
        }
    } else {
        cluster_ids <- as.character(Idents(obj))
        cluster_labels <- if ("cluster_label" %in% colnames(obj@meta.data)) {
            as.character(obj$cluster_label)
        } else {
            cluster_ids
        }
    }

    lookup <- tibble(
        cluster = cluster_ids,
        cluster_label = cluster_labels
    ) %>%
        mutate(
            cluster = ifelse(is.na(cluster) | !nzchar(cluster), NA_character_, cluster),
            cluster_label = ifelse(is.na(cluster_label) | !nzchar(cluster_label), cluster, cluster_label)
        ) %>%
        filter(!is.na(cluster)) %>%
        distinct(cluster, .keep_all = TRUE)

    if (!nrow(lookup)) {
        return(tibble(
            cluster = character(0),
            cluster_label = character(0),
            choice_label = character(0)
        ))
    }

    lookup <- lookup %>%
        mutate(cluster = factor(cluster, levels = cluster_value_levels(cluster))) %>%
        arrange(cluster) %>%
        mutate(cluster = as.character(cluster))

    duplicate_labels <- duplicated(lookup$cluster_label) | duplicated(lookup$cluster_label, fromLast = TRUE)

    lookup %>%
        mutate(
            choice_label = ifelse(
                duplicate_labels,
                paste0(cluster_label, " [", cluster, "]"),
                cluster_label
            )
        )
}

within_group_choices_local <- function(obj) {
    available_cols <- colnames(obj@meta.data)
    choices <- character(0)

    maybe_add <- function(label, column) {
        if (column %in% available_cols && !(column %in% unname(choices))) {
            choices <<- c(choices, setNames(column, label))
        }
    }

    maybe_add("Clustering opt 1", "Rank_1st")
    maybe_add("Clustering opt 2", "Rank_2nd")
    maybe_add("Clustering opt 3", "Rank_3rd")
    maybe_add("Clustering opt 4", "Rank_4th")
    maybe_add("Clustering opt 5", "Rank_5th")
    maybe_add("Time point", "time_point")
    maybe_add("Time point", "Time point")
    maybe_add("Samples", "sample_name")
    maybe_add("Samples", "Sample")
    maybe_add("Samples", "sample")
    choices
}

within_distribution_split_choices_local <- function(obj) {
    available_cols <- colnames(obj@meta.data)
    choices <- c("No split" = "none")

    if ("time_point" %in% available_cols) {
        choices <- c(choices, "Time point" = "time_point")
    } else if ("Time point" %in% available_cols) {
        choices <- c(choices, "Time point" = "Time point")
    } else if ("Group" %in% available_cols) {
        choices <- c(choices, "Time point" = "Group")
    } else if ("condition" %in% available_cols) {
        choices <- c(choices, "Time point" = "condition")
    }

    if ("sample_name" %in% available_cols) {
        choices <- c(choices, "Samples" = "sample_name")
    } else if ("Sample" %in% available_cols) {
        choices <- c(choices, "Samples" = "Sample")
    } else if ("sample" %in% available_cols) {
        choices <- c(choices, "Samples" = "sample")
    }

    choices
}

within_feature_split_choices_local <- function(obj) {
    within_distribution_split_choices_local(obj)
}

within_composition_choices_local <- function(obj) {
    available_cols <- colnames(obj@meta.data)
    choices <- character(0)

    time_col <- pick_first_existing_col(obj@meta.data, c("time_point", "Time point"))
    sample_col <- pick_first_existing_col(obj@meta.data, c("sample_name", "Sample", "sample"))

    if (!is.na(time_col) && time_col %in% available_cols) {
        choices <- c(choices, "Time point" = time_col)
    }

    if (!is.na(sample_col) && sample_col %in% available_cols) {
        choices <- c(choices, "Samples" = sample_col)
    }

    choices
}

cross_composition_choices_local <- function(obj) {
    available_cols <- colnames(obj@meta.data)
    choices <- character(0)

    maybe_add <- function(label, column) {
        if (column %in% available_cols && !(column %in% unname(choices))) {
            choices <<- c(choices, setNames(column, label))
        }
    }

    maybe_add("Species", "species")
    maybe_add("Time point", "time_point")
    maybe_add("Sample", "sample_name")
    maybe_add("Sample", "Sample")
    maybe_add("Sample", "sample")
    choices
}

cross_group_choices_local <- function(obj) {
    available_cols <- colnames(obj@meta.data)
    choices <- character(0)

    maybe_add <- function(label, column) {
        if (column %in% available_cols && !(column %in% unname(choices))) {
            choices <<- c(choices, setNames(column, label))
        }
    }

    maybe_add("Species + label", "species_cell_class")
    maybe_add("Cell class", "cell_class")
    maybe_add("SATURN label", "saturn_label")
    maybe_add("SATURN ref label", "saturn_ref_label")
    maybe_add("Clustering opt 1 label", "Rank_1st_label")
    maybe_add("Clustering opt 2 label", "Rank_2nd_label")
    maybe_add("Clustering opt 3 label", "Rank_3rd_label")
    maybe_add("Clustering opt 4 label", "Rank_4th_label")
    maybe_add("Clustering opt 5 label", "Rank_5th_label")
    maybe_add("Species", "species")
    maybe_add("Cluster", "cluster_label")
    maybe_add("Clustering opt 1", "Rank_1st")
    maybe_add("Clustering opt 2", "Rank_2nd")
    maybe_add("Clustering opt 3", "Rank_3rd")
    maybe_add("Clustering opt 4", "Rank_4th")
    maybe_add("Clustering opt 5", "Rank_5th")
    maybe_add("Condition", "condition")
    maybe_add("Study", "study")
    maybe_add("Time point", "time_point")
    choices
}

cross_distribution_group_choices_local <- function(obj) {
    available_cols <- colnames(obj@meta.data)
    choices <- character(0)

    maybe_add <- function(label, column) {
        if (column %in% available_cols && !(column %in% unname(choices))) {
            choices <<- c(choices, setNames(column, label))
        }
    }

    maybe_add("Time point", "time_point")
    maybe_add("Species", "species")
    maybe_add("Clustering opt 1 label", "Rank_1st_label")
    maybe_add("Clustering opt 2 label", "Rank_2nd_label")
    maybe_add("Clustering opt 3 label", "Rank_3rd_label")
    maybe_add("Clustering opt 4 label", "Rank_4th_label")
    maybe_add("Clustering opt 5 label", "Rank_5th_label")
    maybe_add("Clustering opt 1", "Rank_1st")
    maybe_add("Clustering opt 2", "Rank_2nd")
    maybe_add("Clustering opt 3", "Rank_3rd")
    maybe_add("Clustering opt 4", "Rank_4th")
    maybe_add("Clustering opt 5", "Rank_5th")
    choices
}

choice_rows <- function(dataset_key, choice_family, choices) {
    if (!length(choices)) {
        return(tibble(
            dataset_key = character(0),
            choice_family = character(0),
            label = character(0),
            value = character(0),
            sort_order = integer(0)
        ))
    }

    tibble(
        dataset_key = dataset_key,
        choice_family = choice_family,
        label = names(choices),
        value = unname(choices),
        sort_order = seq_along(choices)
    )
}

cluster_lookup_rows <- function(dataset_key, obj, cluster_sources) {
    bind_rows(lapply(cluster_sources, function(cluster_source) {
        lookup <- cluster_label_lookup_local(obj, cluster_source)

        if (!nrow(lookup)) {
            return(NULL)
        }

        lookup %>%
            mutate(
                dataset_key = dataset_key,
                cluster_source = cluster_source,
                sort_order = seq_len(n())
            ) %>%
            select(dataset_key, cluster_source, cluster, cluster_label, choice_label, sort_order)
    }))
}

get_within_dataset_path <- function(species_key, integration_method) {
    path <- species_registry[[species_key]]$within_paths[[integration_method]]
    pick_first_existing_path(c(app_slim_path(path), path))
}

get_cross_dataset_path <- function(cross_key) {
    cfg <- cross_integration_registry[[cross_key]]
    pick_first_existing_path(c(cfg$slim_path, cfg$path))
}

build_within_cache <- function(species_key, integration_method) {
    dataset_key <- paste(species_key, integration_method, sep = "_")
    obj <- readRDS(get_within_dataset_path(species_key, integration_method))
    obj <- apply_cluster_annotation_overlay(
        obj,
        dataset_key = dataset_key,
        annotations_dir = cluster_annotations_dir,
        legacy_overrides_dir = celltype_overrides_dir
    )

    choice_tbl <- bind_rows(
        choice_rows(dataset_key, "within_group", within_group_choices_local(obj)),
        choice_rows(dataset_key, "within_distribution_split", within_distribution_split_choices_local(obj)),
        choice_rows(dataset_key, "within_feature_split", within_feature_split_choices_local(obj)),
        choice_rows(dataset_key, "within_composition", within_composition_choices_local(obj))
    )

    cluster_sources <- unique(c(
        "__idents__",
        distribution_cluster_columns[distribution_cluster_columns %in% colnames(obj@meta.data)]
    ))

    list(
        choice_tbl = choice_tbl,
        cluster_tbl = cluster_lookup_rows(dataset_key, obj, cluster_sources)
    )
}

build_cross_cache <- function(cross_key) {
    obj <- readRDS(get_cross_dataset_path(cross_key))
    obj <- apply_cluster_annotation_overlay(
        obj,
        dataset_key = cross_key,
        annotations_dir = cluster_annotations_dir,
        legacy_overrides_dir = celltype_overrides_dir
    )

    choice_tbl <- bind_rows(
        choice_rows(cross_key, "cross_group", cross_group_choices_local(obj)),
        choice_rows(cross_key, "cross_distribution_group", cross_distribution_group_choices_local(obj)),
        choice_rows(cross_key, "cross_composition", cross_composition_choices_local(obj))
    )

    cluster_sources <- unique(c(
        "__idents__",
        c("species_cell_class", "cell_class", "cluster_label", "Rank_1st_label", "Rank_1st", "Rank_2nd", "Rank_3rd", "Rank_4th", "Rank_5th", "saturn_label", "saturn_ref_label")[c("species_cell_class", "cell_class", "cluster_label", "Rank_1st_label", "Rank_1st", "Rank_2nd", "Rank_3rd", "Rank_4th", "Rank_5th", "saturn_label", "saturn_ref_label") %in% colnames(obj@meta.data)]
    ))

    list(
        choice_tbl = choice_tbl,
        cluster_tbl = cluster_lookup_rows(cross_key, obj, cluster_sources)
    )
}

within_cache <- lapply(within_species_keys, function(species_key) {
    lapply(integration_methods, function(integration_method) {
        build_within_cache(species_key, integration_method)
    })
})

within_choice_tbl <- bind_rows(lapply(within_cache, function(x) bind_rows(lapply(x, `[[`, "choice_tbl"))))
within_cluster_tbl <- bind_rows(lapply(within_cache, function(x) bind_rows(lapply(x, `[[`, "cluster_tbl"))))

cross_cache <- lapply(names(cross_integration_registry), build_cross_cache)
cross_choice_tbl <- bind_rows(lapply(cross_cache, `[[`, "choice_tbl"))
cross_cluster_tbl <- bind_rows(lapply(cross_cache, `[[`, "cluster_tbl"))

choice_tbl <- bind_rows(within_choice_tbl, cross_choice_tbl) %>%
    arrange(dataset_key, choice_family, sort_order)

cluster_tbl <- bind_rows(within_cluster_tbl, cross_cluster_tbl) %>%
    arrange(dataset_key, cluster_source, sort_order)

dir.create(dirname(ui_choice_cache_path), recursive = TRUE, showWarnings = FALSE)

write.table(
    choice_tbl,
    file = ui_choice_cache_path,
    sep = "\t",
    row.names = FALSE,
    quote = TRUE,
    na = ""
)

write.table(
    cluster_tbl,
    file = ui_cluster_lookup_cache_path,
    sep = "\t",
    row.names = FALSE,
    quote = TRUE,
    na = ""
)

cat("Wrote", ui_choice_cache_path, "\n")
cat("Wrote", ui_cluster_lookup_cache_path, "\n")
