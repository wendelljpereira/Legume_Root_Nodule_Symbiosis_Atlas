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

within_species_keys <- c("medicago", "glycine", "lotus")

species_registry <- list(
    medicago = list(
        within_paths = c(
            ComBat_BBKNN = "within_species_integrated_datasets/ComBat_BBKNN/M_truncatula_clustered_dataset.rds",
            Seurat = "within_species_integrated_datasets/Seurat/M_truncatula_clustered_dataset.rds"
        ),
        canonicalize = function(x) trimws(as.character(x))
    ),
    glycine = list(
        within_paths = c(
            ComBat_BBKNN = "within_species_integrated_datasets/ComBat_BBKNN/G_max_clustered_dataset.rds",
            Seurat = "within_species_integrated_datasets/Seurat/G_max_clustered_dataset.rds"
        ),
        canonicalize = function(x) trimws(as.character(x))
    ),
    lotus = list(
        within_paths = c(
            ComBat_BBKNN = "within_species_integrated_datasets/ComBat_BBKNN/L_japonicus_clustered_dataset.rds",
            Seurat = "within_species_integrated_datasets/Seurat/L_japonicus_clustered_dataset.rds"
        ),
        canonicalize = function(x) {
            x <- trimws(as.character(x))
            gsub("_LC$", "", x)
        }
    )
)

annotation_paths <- c(
    medicago = "annotations/medicago_gene_annotations.tsv",
    glycine = "annotations/glycine_gene_annotations.tsv",
    lotus = "annotations/lotus_gene_annotations.tsv"
)

integration_methods <- c("ComBat_BBKNN", "Seurat")
gene_catalog_cache_dir <- "metadata/gene_catalogs"
cross_species_path <- "app_ready_integration/camex/clustered_dataset.rds"
cross_species_slim_path <- "app_ready_integration/camex/clustered_dataset_app_slim.rds"
cross_feature_lookup_path <- "metadata/cross_feature_lookup.tsv"

get_within_dataset_path <- function(species_key, integration_method) {
    path <- species_registry[[species_key]]$within_paths[[integration_method]]
    pick_first_existing_path(c(app_slim_path(path), path))
}

canonicalize_gene_ids <- function(species_key, ids) {
    ids <- as.character(ids)
    out <- species_registry[[species_key]]$canonicalize(ids)
    out[is.na(ids) | !nzchar(trimws(ids))] <- NA_character_
    out
}

first_nonempty <- function(x) {
    x <- trimws(as.character(x))
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x)) x[[1]] else NA_character_
}

pick_first_existing_col <- function(df, candidates) {
    match_col <- intersect(candidates, colnames(df))[1]
    if (is.na(match_col)) NA_character_ else match_col
}

read_gene_annotations <- function(species_key) {
    path <- annotation_paths[[species_key]]

    if (is.na(path) || !file.exists(path)) {
        return(data.frame(
            canonical_gene_id = character(0),
            annotation_gene_id = character(0),
            common_name = character(0),
            synonyms = character(0),
            description = character(0),
            stringsAsFactors = FALSE
        ))
    }

    raw_df <- read.delim(
        path,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    id_col <- pick_first_existing_col(raw_df, c("gene_id", "id", "locus_tag"))
    common_name_col <- pick_first_existing_col(raw_df, c("common_name", "gene_name", "symbol", "acronym", "name"))
    synonyms_col <- pick_first_existing_col(raw_df, c("synonyms", "aliases", "alias"))
    description_col <- pick_first_existing_col(raw_df, c("description", "geneProduct", "product", "genePublication"))

    if (is.na(id_col)) {
        return(data.frame(
            canonical_gene_id = character(0),
            annotation_gene_id = character(0),
            common_name = character(0),
            synonyms = character(0),
            description = character(0),
            stringsAsFactors = FALSE
        ))
    }

    data.frame(
        annotation_gene_id = as.character(raw_df[[id_col]]),
        canonical_gene_id = canonicalize_gene_ids(species_key, raw_df[[id_col]]),
        common_name = if (!is.na(common_name_col)) trimws(as.character(raw_df[[common_name_col]])) else NA_character_,
        synonyms = if (!is.na(synonyms_col)) trimws(as.character(raw_df[[synonyms_col]])) else NA_character_,
        description = if (!is.na(description_col)) trimws(as.character(raw_df[[description_col]])) else NA_character_,
        stringsAsFactors = FALSE
    ) %>%
        filter(!is.na(canonical_gene_id) & nzchar(canonical_gene_id)) %>%
        group_by(canonical_gene_id) %>%
        summarise(
            annotation_gene_id = first_nonempty(annotation_gene_id),
            common_name = first_nonempty(common_name),
            synonyms = paste(unique(synonyms[!is.na(synonyms) & nzchar(synonyms)]), collapse = "; "),
            description = first_nonempty(description),
            .groups = "drop"
        ) %>%
        mutate(
            synonyms = ifelse(synonyms == "", NA_character_, synonyms),
            description = ifelse(description == "", NA_character_, description)
        )
}

build_gene_catalog <- function(species_key, integration_method) {
    obj <- readRDS(get_within_dataset_path(species_key, integration_method))
    feature_ids <- rownames(obj)
    annotations <- read_gene_annotations(species_key)

    data.frame(
        feature_id = feature_ids,
        canonical_gene_id = canonicalize_gene_ids(species_key, feature_ids),
        stringsAsFactors = FALSE
    ) %>%
        left_join(annotations, by = "canonical_gene_id") %>%
        mutate(
            display_label = ifelse(
                !is.na(common_name) & nzchar(common_name) & common_name != feature_id,
                paste0(common_name, " (", feature_id, ")"),
                feature_id
            ),
            search_tokens = vapply(
                seq_len(n()),
                function(i) {
                    values <- c(feature_id[i], common_name[i], synonyms[i])
                    values <- unique(trimws(as.character(values)))
                    values <- values[!is.na(values) & nzchar(values)]
                    paste(values, collapse = " ")
                },
                character(1)
            )
        )
}

dir.create(gene_catalog_cache_dir, recursive = TRUE, showWarnings = FALSE)

legacy_catalogs <- list.files(
    gene_catalog_cache_dir,
    pattern = "_(ComBat_BBKNN|Seurat)\\.tsv$",
    full.names = TRUE
)

if (length(legacy_catalogs)) {
    file.remove(legacy_catalogs)
}

for (species_key in within_species_keys) {
    catalog <- build_gene_catalog(species_key, integration_methods[[1]])
    out_path <- file.path(gene_catalog_cache_dir, paste0(species_key, ".tsv"))
    write.table(
        catalog,
        file = out_path,
        sep = "\t",
        row.names = FALSE,
        quote = TRUE,
        na = ""
    )
    cat("Wrote", out_path, "\n")
}

cross_obj <- readRDS(pick_first_existing_path(c(cross_species_slim_path, cross_species_path)))
cross_lookup <- data.frame(
    feature_id = rownames(cross_obj),
    canonical_gene_id = canonicalize_gene_ids("medicago", rownames(cross_obj)),
    stringsAsFactors = FALSE
) %>%
    filter(!is.na(canonical_gene_id) & nzchar(canonical_gene_id)) %>%
    distinct(canonical_gene_id, feature_id)

write.table(
    cross_lookup,
    file = cross_feature_lookup_path,
    sep = "\t",
    row.names = FALSE,
    quote = TRUE,
    na = ""
)

cat("Wrote", cross_feature_lookup_path, "\n")
