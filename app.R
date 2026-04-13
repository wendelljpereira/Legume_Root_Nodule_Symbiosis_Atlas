library(shiny)
library(Seurat)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(patchwork)
library(svglite)
library(scCustomize)

# ==============================================================================
# 1. GLOBAL SETTINGS, DATA REGISTRY, AND HELPERS
# ==============================================================================

`%||%` <- function(x, y) {
    if (is.null(x) || !length(x) || all(is.na(x))) y else x
}

within_species_keys <- c("medicago", "glycine", "lotus")

species_registry <- list(
    medicago = list(
        label = "Medicago truncatula",
        orthogroup_col = "medicago.fa",
        within_paths = c(
            ComBat_BBKNN = "within_species_integrated_datasets/ComBat_BBKNN/M_truncatula_clustered_dataset.rds",
            Seurat = "within_species_integrated_datasets/Seurat/M_truncatula_clustered_dataset.rds"
        ),
        canonicalize = function(x) trimws(as.character(x))
    ),
    glycine = list(
        label = "Glycine max",
        orthogroup_col = "glycine.fa",
        within_paths = c(
            ComBat_BBKNN = "within_species_integrated_datasets/ComBat_BBKNN/G_max_clustered_dataset.rds",
            Seurat = "within_species_integrated_datasets/Seurat/G_max_clustered_dataset.rds"
        ),
        canonicalize = function(x) trimws(as.character(x))
    ),
    lotus = list(
        label = "Lotus japonicus",
        orthogroup_col = "lotus.fa",
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

cross_species_path <- "cross_species_integrated_datasets/Camex/clustered_dataset.rds"
atlas_summary_path <- "metadata/atlas_summary.tsv"
gene_catalog_cache_dir <- "metadata/gene_catalogs"
cross_feature_lookup_path <- "metadata/cross_feature_lookup.tsv"

species_choices <- c(
    "Medicago truncatula" = "medicago",
    "Glycine max" = "glycine",
    "Lotus japonicus" = "lotus"
)

integration_choices <- c(
    "ComBat/BBKNN" = "ComBat_BBKNN",
    "Seurat" = "Seurat"
)

annotation_paths <- c(
    medicago = "annotations/medicago_gene_annotations.tsv",
    glycine = "annotations/glycine_gene_annotations.tsv",
    lotus = "annotations/lotus_gene_annotations.tsv"
)

app_palette <- c(
    green_dark = "#194430",
    green = "#2f6b4f",
    green_soft = "#dce9df",
    warm = "#b36f22",
    warm_soft = "#fff0dd",
    text = "#193526",
    muted = "#5a7263",
    border = "#c9d6c4",
    slate = "#8aa59a",
    red_soft = "#f9ece8"
)

species_label <- function(species_key) {
    species_registry[[species_key]]$label
}

species_label_tag <- function(species_key) {
    tags$span(class = "species-binomial", species_label(species_key))
}

canonicalize_gene_ids <- function(species_key, ids) {
    ids <- as.character(ids)
    out <- species_registry[[species_key]]$canonicalize(ids)
    out[is.na(ids) | !nzchar(trimws(ids))] <- NA_character_
    out
}

format_stat_value <- function(value) {
    format(as.numeric(value), big.mark = ",", scientific = FALSE, trim = TRUE)
}

first_nonempty <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x <- x[!is.na(x) & nzchar(x)]

    if (length(x)) x[[1]] else NA_character_
}

pick_first_existing_col <- function(df, candidates) {
    match_col <- intersect(candidates, colnames(df))[1]

    if (is.na(match_col)) {
        NA_character_
    } else {
        match_col
    }
}

split_orthogroup_genes <- function(value) {
    if (length(value) == 0 || is.null(value) || all(is.na(value))) {
        return(character(0))
    }

    value <- trimws(as.character(value[[1]]))

    if (!nzchar(value)) {
        return(character(0))
    }

    genes <- trimws(unlist(strsplit(value, ",", fixed = TRUE)))
    genes[nzchar(genes)]
}

compact_gene_list <- function(values, limit = 5, empty_label = "None") {
    values <- unique(as.character(values))
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

display_gene_labels <- function(species_key, gene_ids) {
    gene_ids <- as.character(gene_ids)

    annotations <- read_gene_annotations(species_key)

    tibble(
        gene_id = gene_ids,
        canonical_gene_id = canonicalize_gene_ids(species_key, gene_ids)
    ) %>%
        left_join(annotations, by = "canonical_gene_id") %>%
        mutate(
            display_label = if_else(
                !is.na(common_name) &
                    nzchar(common_name) &
                    common_name != gene_id,
                paste0(common_name, " (", gene_id, ")"),
                gene_id
            )
        ) %>%
        pull(display_label)
}

metric_tile <- function(value, label) {
    div(
        class = "metric-tile",
        div(class = "metric-value", value),
        div(class = "metric-label", label)
    )
}

plot_download_button <- function(output_id) {
    downloadButton(
        outputId = output_id,
        label = tagList(icon("download"), "Download"),
        class = "btn btn-default btn-sm plot-download-btn"
    )
}

spinning_plot_output <- function(output_id, proxy_height = "360px") {
    plotOutput(output_id, height = "auto")
}

notice_card <- function(title, body, tone = c("info", "warning")) {
    tone <- match.arg(tone)

    div(
        class = paste("alert-card", tone),
        div(class = "alert-title", title),
        div(class = "alert-body", body)
    )
}

html_summary_table <- function(df) {
    if (!nrow(df)) {
        return(div(class = "summary-placeholder", "No rows are available for the current selection."))
    }

    display_df <- df %>%
        mutate(across(everything(), ~ ifelse(is.na(.x), "", as.character(.x))))

    header <- tags$tr(lapply(names(display_df), tags$th))
    body_rows <- lapply(seq_len(nrow(display_df)), function(i) {
        tags$tr(
            lapply(seq_len(ncol(display_df)), function(j) {
                tags$td(display_df[[j]][i])
            })
        )
    })

    div(
        class = "table-scroll",
        tags$table(
            class = "summary-table",
            tags$thead(header),
            tags$tbody(body_rows)
        )
    )
}

app_plot_theme <- function(base_size = 13) {
    theme_minimal(base_size = base_size) +
        theme(
            plot.title = element_text(face = "bold", colour = app_palette["text"]),
            plot.subtitle = element_text(colour = app_palette["muted"]),
            axis.title = element_text(face = "bold", colour = app_palette["text"]),
            axis.text = element_text(colour = app_palette["muted"]),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            legend.position = "top",
            legend.title = element_text(face = "bold"),
            strip.text = element_text(face = "bold", colour = app_palette["text"]),
            strip.background = element_rect(fill = app_palette["green_soft"], colour = NA),
            plot.margin = margin(10, 16, 10, 10)
        )
}

plot_title_annotation <- function(title) {
    plot_annotation(
        title = title,
        theme = theme(
            plot.title = element_text(
                face = "bold",
                colour = app_palette["text"],
                size = 16,
                hjust = 0
            )
        )
    )
}

wrap_titled_plot <- function(plot_obj, title) {
    wrap_elements(full = plot_obj + plot_title_annotation(title))
}

resolve_feature_grid_cols <- function(requested_cols, feature_labels, split_by = "none") {
    requested_cols <- max(1L, as.integer(requested_cols %||% 1L))
    feature_labels <- as.character(feature_labels %||% character(0))

    if (!identical(split_by, "none")) {
        return(1L)
    }

    if (length(feature_labels) > 4L || any(nchar(feature_labels) > 34L, na.rm = TRUE)) {
        return(1L)
    }

    requested_cols
}

reorder_within <- function(x, by, within, fun = mean, sep = "___") {
    grouped <- paste(x, within, sep = sep)
    stats <- tapply(by, grouped, fun)
    factor(grouped, levels = names(sort(stats)))
}

scale_y_reordered <- function(sep = "___") {
    scale_y_discrete(labels = function(x) sub(paste0(sep, ".*$"), "", x))
}

resolve_choice <- function(value, choices, default = NULL) {
    choice_values <- unname(choices)
    fallback <- default %||% choice_values[[1]]

    if (length(value) && !is.null(value) && value %in% choice_values) {
        value
    } else {
        fallback
    }
}

gene_catalog_cache_path <- function(species_key, integration_method) {
    file.path(
        gene_catalog_cache_dir,
        paste0(species_key, "_", integration_method, ".tsv")
    )
}

read_tsv_cache <- function(path) {
    if (!file.exists(path)) {
        return(NULL)
    }

    read.delim(
        path,
        sep = "\t",
        stringsAsFactors = FALSE,
        na.strings = c("", "NA"),
        check.names = FALSE
    ) %>%
        as_tibble()
}

split_panel_count <- function(obj, split_by) {
    if (is.null(split_by) || !length(split_by) || identical(split_by, "none")) {
        return(1L)
    }

    values <- obj@meta.data[[split_by]]

    if (is.null(values)) {
        return(1L)
    }

    values <- as.character(values)
    values <- values[!is.na(values) & nzchar(values)]
    max(1L, length(unique(values)))
}

prepare_within_object <- function(obj) {
    obj$cluster_label <- as.character(Idents(obj))
    obj
}

prepare_cross_object <- function(obj) {
    obj$cluster_label <- as.character(Idents(obj))
    obj$species_cell_class <- paste(
        as.character(obj$species %||% "Unknown"),
        as.character(obj$cell_class %||% "Unknown"),
        sep = " | "
    )
    obj
}

dataset_cache <- new.env(parent = emptyenv())

cache_get <- function(key, builder) {
    if (!exists(key, envir = dataset_cache, inherits = FALSE)) {
        assign(key, builder(), envir = dataset_cache)
    }

    get(key, envir = dataset_cache, inherits = FALSE)
}

get_within_object <- function(species_key, integration_method) {
    cache_key <- paste("within", species_key, integration_method, sep = "::")

    cache_get(cache_key, function() {
        readRDS(species_registry[[species_key]]$within_paths[[integration_method]]) %>%
            prepare_within_object()
    })
}

get_cross_object <- function() {
    cache_get("cross::object", function() {
        readRDS(cross_species_path) %>%
            prepare_cross_object()
    })
}

get_within_feature_lookup <- function(species_key, integration_method) {
    cache_key <- paste("within_lookup", species_key, integration_method, sep = "::")

    cache_get(cache_key, function() {
        cache_path <- gene_catalog_cache_path(species_key, integration_method)
        cached_lookup <- read_tsv_cache(cache_path)

        if (!is.null(cached_lookup) &&
            all(c("feature_id", "canonical_gene_id") %in% colnames(cached_lookup))) {
            return(
                cached_lookup %>%
                    select(feature_id, canonical_gene_id) %>%
                    filter(!is.na(feature_id) & nzchar(feature_id)) %>%
                    filter(!is.na(canonical_gene_id) & nzchar(canonical_gene_id)) %>%
                    distinct(canonical_gene_id, feature_id)
            )
        }

        feature_ids <- rownames(get_within_object(species_key, integration_method))

        tibble(
            feature_id = feature_ids,
            canonical_gene_id = canonicalize_gene_ids(species_key, feature_ids)
        ) %>%
            distinct(canonical_gene_id, feature_id)
    })
}

get_cross_feature_lookup <- function() {
    cache_get("cross::lookup", function() {
        cached_lookup <- read_tsv_cache(cross_feature_lookup_path)

        if (!is.null(cached_lookup) &&
            all(c("feature_id", "canonical_gene_id") %in% colnames(cached_lookup))) {
            return(
                cached_lookup %>%
                    select(feature_id, canonical_gene_id) %>%
                    filter(!is.na(feature_id) & nzchar(feature_id)) %>%
                    filter(!is.na(canonical_gene_id) & nzchar(canonical_gene_id)) %>%
                    distinct(canonical_gene_id, feature_id)
            )
        }

        feature_ids <- rownames(get_cross_object())

        tibble(
            feature_id = feature_ids,
            canonical_gene_id = canonicalize_gene_ids("medicago", feature_ids)
        ) %>%
            distinct(canonical_gene_id, feature_id)
    })
}

read_gene_annotations <- function(species_key) {
    cache_key <- paste("annotation", species_key, sep = "::")

    cache_get(cache_key, function() {
        path <- annotation_paths[[species_key]]

        if (is.na(path) || !file.exists(path)) {
            return(tibble(
                canonical_gene_id = character(0),
                annotation_gene_id = character(0),
                common_name = character(0),
                synonyms = character(0),
                description = character(0)
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
            return(tibble(
                canonical_gene_id = character(0),
                annotation_gene_id = character(0),
                common_name = character(0),
                synonyms = character(0),
                description = character(0)
            ))
        }

        tibble(
            annotation_gene_id = as.character(raw_df[[id_col]]),
            canonical_gene_id = canonicalize_gene_ids(species_key, raw_df[[id_col]]),
            common_name = if (!is.na(common_name_col)) trimws(as.character(raw_df[[common_name_col]])) else NA_character_,
            synonyms = if (!is.na(synonyms_col)) trimws(as.character(raw_df[[synonyms_col]])) else NA_character_,
            description = if (!is.na(description_col)) trimws(as.character(raw_df[[description_col]])) else NA_character_
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
                synonyms = if_else(synonyms == "", NA_character_, synonyms),
                description = if_else(description == "", NA_character_, description)
            )
    })
}

build_gene_catalog <- function(species_key, integration_method) {
    cache_key <- paste("gene_catalog", species_key, integration_method, sep = "::")

    cache_get(cache_key, function() {
        cache_path <- gene_catalog_cache_path(species_key, integration_method)
        cached_catalog <- read_tsv_cache(cache_path)

        if (!is.null(cached_catalog) &&
            all(c("feature_id", "display_label", "search_tokens") %in% colnames(cached_catalog))) {
            return(cached_catalog %>% mutate(across(everything(), as.character)))
        }

        annotations <- read_gene_annotations(species_key)
        obj <- get_within_object(species_key, integration_method)
        feature_ids <- rownames(obj)

        tibble(
            feature_id = feature_ids,
            canonical_gene_id = canonicalize_gene_ids(species_key, feature_ids)
        ) %>%
            left_join(annotations, by = "canonical_gene_id") %>%
            mutate(
                display_label = display_gene_labels(species_key, feature_id),
                search_tokens = pmap_chr(
                    list(feature_id, common_name, synonyms),
                    function(feature_id, common_name, synonyms) {
                        values <- c(feature_id, common_name, synonyms)
                        values <- unique(trimws(as.character(values)))
                        values <- values[!is.na(values) & nzchar(values)]
                        paste(values, collapse = " ")
                    }
                )
            )
    })
}

build_gene_choices <- function(species_key, integration_method) {
    cache_key <- paste("gene_choices", species_key, integration_method, sep = "::")

    cache_get(cache_key, function() {
        catalog <- build_gene_catalog(species_key, integration_method) %>%
            arrange(display_label)

        list(
            choices = setNames(catalog$feature_id, catalog$display_label),
            tokens = catalog$search_tokens,
            feature_ids = catalog$feature_id
        )
    })
}

pick_sample_column <- function(df) {
    pick_first_existing_col(df, c("sample_name", "sample", "Sample", "orig.ident"))
}

same_gene_selection <- function(x, y) {
    identical(sort(unique(as.character(x))), sort(unique(as.character(y))))
}

summarise_value_preview <- function(values, limit = 4, empty_label = "NA") {
    values <- unique(as.character(values))
    values <- values[!is.na(values) & nzchar(values)]

    compact_gene_list(values, limit = limit, empty_label = empty_label)
}

build_within_dataset_summary_row <- function(species_key, integration_method, obj) {
    md <- obj@meta.data
    sample_col <- pick_sample_column(md)
    sample_values <- if (!is.na(sample_col)) md[[sample_col]] else character(0)

    tibble(
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
        time_preview = NA_character_
    )
}

build_cross_dataset_summary_row <- function(obj) {
    md <- obj@meta.data
    sample_col <- pick_sample_column(md)
    sample_values <- if (!is.na(sample_col)) md[[sample_col]] else character(0)

    tibble(
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
        time_preview = NA_character_
    )
}

compute_atlas_summary_table <- function() {
    within_rows <- bind_rows(lapply(within_species_keys, function(species_key) {
        bind_rows(lapply(unname(integration_choices), function(integration_method) {
            build_within_dataset_summary_row(
                species_key = species_key,
                integration_method = integration_method,
                obj = get_within_object(species_key, integration_method)
            )
        }))
    }))

    bind_rows(
        within_rows,
        build_cross_dataset_summary_row(get_cross_object())
    )
}

get_atlas_summary_table <- function() {
    cache_get("atlas::summary_table", function() {
        if (file.exists(atlas_summary_path)) {
            summary_df <- read.delim(
                atlas_summary_path,
                sep = "\t",
                stringsAsFactors = FALSE,
                na.strings = c("", "NA"),
                check.names = FALSE
            ) %>%
                as_tibble()

            if (nrow(summary_df) && "sample_n" %in% colnames(summary_df)) {
                return(summary_df)
            }
        }

        compute_atlas_summary_table()
    })
}

get_within_dataset_summary <- function(species_key, integration_method) {
    summary_row <- get_atlas_summary_table() %>%
        filter(
            dataset_scope == "within",
            species_key == !!species_key,
            integration_method == !!integration_method
        ) %>%
        slice_head(n = 1)

    if (!nrow(summary_row)) {
        summary_row <- build_within_dataset_summary_row(
            species_key = species_key,
            integration_method = integration_method,
            obj = get_within_object(species_key, integration_method)
        )
    }

    as.list(summary_row)
}

get_cross_dataset_summary <- function() {
    summary_row <- get_atlas_summary_table() %>%
        filter(dataset_scope == "cross") %>%
        slice_head(n = 1)

    if (!nrow(summary_row)) {
        summary_row <- build_cross_dataset_summary_row(get_cross_object())
    }

    as.list(summary_row)
}

orthogroup_raw <- read.delim(
    "orthogroups/joint_orthogroups.tsv",
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

orthogroup_members <- tibble(
    orthogroup = orthogroup_raw$Orthogroup,
    medicago = lapply(orthogroup_raw[[species_registry$medicago$orthogroup_col]], split_orthogroup_genes),
    glycine = lapply(orthogroup_raw[[species_registry$glycine$orthogroup_col]], split_orthogroup_genes),
    lotus = lapply(orthogroup_raw[[species_registry$lotus$orthogroup_col]], split_orthogroup_genes)
)

orthogroup_long <- bind_rows(lapply(within_species_keys, function(species_key) {
    tibble(
        orthogroup = orthogroup_members$orthogroup,
        gene_id = orthogroup_members[[species_key]]
    ) %>%
        unnest_longer(gene_id, values_to = "gene_id") %>%
        mutate(
            species = species_key,
            canonical_gene_id = canonicalize_gene_ids(species_key, gene_id)
        )
})) %>%
    distinct(species, canonical_gene_id, orthogroup, gene_id)

orthogroup_index <- orthogroup_long %>%
    select(species, canonical_gene_id, orthogroup) %>%
    distinct()

get_orthogroup_members <- function(orthogroup, target_species) {
    if (is.na(orthogroup) || !length(orthogroup)) {
        return(character(0))
    }

    row_index <- match(orthogroup, orthogroup_members$orthogroup)

    if (is.na(row_index)) {
        return(character(0))
    }

    orthogroup_members[[target_species]][[row_index]] %||% character(0)
}

resolve_source_orthogroups <- function(source_species, source_genes) {
    source_genes <- unique(as.character(source_genes))

    if (!length(source_genes)) {
        return(tibble(
            source_gene = character(0),
            source_canonical = character(0),
            orthogroup = character(0)
        ))
    }

    tibble(
        source_gene = source_genes,
        source_canonical = canonicalize_gene_ids(source_species, source_genes)
    ) %>%
        left_join(
            orthogroup_index %>% filter(species == source_species),
            by = c("source_canonical" = "canonical_gene_id")
        )
}

match_target_features <- function(target_species, target_genes, integration_method = NULL, cross_space = FALSE) {
    target_genes <- unique(as.character(target_genes))
    target_genes <- target_genes[!is.na(target_genes) & nzchar(target_genes)]

    if (!length(target_genes)) {
        return(character(0))
    }

    lookup <- if (isTRUE(cross_space)) {
        get_cross_feature_lookup()
    } else {
        get_within_feature_lookup(target_species, integration_method)
    }

    effective_species <- if (isTRUE(cross_space)) "medicago" else target_species

    tibble(
        canonical_gene_id = canonicalize_gene_ids(effective_species, target_genes)
    ) %>%
        inner_join(lookup, by = "canonical_gene_id") %>%
        pull(feature_id) %>%
        unique()
}

resolve_target_mapping <- function(source_species, source_genes, target_species, integration_method, cross_space = FALSE) {
    source_genes <- unique(as.character(source_genes))

    if (!length(source_genes)) {
        return(list(
            mapping = tibble(),
            plot_table = tibble(),
            plot_features = character(0),
            label_map = c(),
            no_orthogroup = character(0),
            no_target_members = character(0),
            missing_features = character(0),
            multiplicity = tibble()
        ))
    }

    effective_target_species <- if (isTRUE(cross_space)) "medicago" else target_species
    target_lookup <- if (isTRUE(cross_space)) {
        get_cross_feature_lookup()
    } else {
        get_within_feature_lookup(target_species, integration_method)
    }

    if (!isTRUE(cross_space) && identical(source_species, target_species)) {
        mapped_tbl <- tibble(
            source_gene = source_genes,
            source_canonical = canonicalize_gene_ids(source_species, source_genes),
            orthogroup = NA_character_,
            target_gene_original = source_genes,
            target_canonical = canonicalize_gene_ids(target_species, source_genes)
        ) %>%
            left_join(target_lookup, by = c("target_canonical" = "canonical_gene_id")) %>%
            rename(target_feature_id = feature_id)

        no_orthogroup <- character(0)
        no_target_members <- character(0)
    } else {
        source_orthogroups <- resolve_source_orthogroups(source_species, source_genes)

        genes_with_orthogroups <- source_orthogroups %>%
            filter(!is.na(orthogroup)) %>%
            distinct(source_gene) %>%
            pull(source_gene)

        no_orthogroup <- setdiff(source_genes, genes_with_orthogroups)

        candidate_tbl <- source_orthogroups %>%
            mutate(
                target_gene_original = lapply(
                    orthogroup,
                    function(og) get_orthogroup_members(og, effective_target_species)
                )
            ) %>%
            unnest_longer(target_gene_original, keep_empty = TRUE) %>%
            mutate(
                target_gene_original = ifelse(
                    is.na(target_gene_original),
                    NA_character_,
                    as.character(target_gene_original)
                ),
                target_canonical = canonicalize_gene_ids(
                    effective_target_species,
                    target_gene_original
                )
            )

        genes_with_target_members <- candidate_tbl %>%
            filter(!is.na(target_gene_original) & nzchar(target_gene_original)) %>%
            distinct(source_gene) %>%
            pull(source_gene)

        no_target_members <- setdiff(genes_with_orthogroups, genes_with_target_members)

        mapped_tbl <- candidate_tbl %>%
            left_join(target_lookup, by = c("target_canonical" = "canonical_gene_id")) %>%
            rename(target_feature_id = feature_id)
    }

    genes_with_candidate_features <- mapped_tbl %>%
        filter(!is.na(target_gene_original) & nzchar(target_gene_original)) %>%
        distinct(source_gene) %>%
        pull(source_gene)

    genes_with_mapped_features <- mapped_tbl %>%
        filter(!is.na(target_feature_id)) %>%
        distinct(source_gene) %>%
        pull(source_gene)

    missing_features <- setdiff(genes_with_candidate_features, genes_with_mapped_features)

    plot_table <- mapped_tbl %>%
        filter(!is.na(target_feature_id)) %>%
        distinct(source_gene, orthogroup, target_gene_original, target_feature_id) %>%
        group_by(target_feature_id) %>%
        summarise(
            source_gene_label = paste(unique(source_gene), collapse = "; "),
            orthogroup_label = paste(unique(na.omit(orthogroup)), collapse = "; "),
            target_gene_label = paste(unique(na.omit(target_gene_original)), collapse = "; "),
            .groups = "drop"
        ) %>%
        mutate(
            source_gene_display = map_chr(
                source_gene_label,
                function(label_text) {
                    genes <- trimws(unlist(strsplit(label_text, ";", fixed = TRUE)))
                    genes <- genes[nzchar(genes)]
                    paste(unique(display_gene_labels(source_species, genes)), collapse = "; ")
                }
            )
        )

    target_display_lookup <- setNames(
        display_gene_labels(effective_target_species, plot_table$target_feature_id),
        plot_table$target_feature_id
    )

    plot_table <- plot_table %>%
        mutate(
            plot_label = if (!isTRUE(cross_space) && identical(source_species, target_species)) {
                unname(target_display_lookup[target_feature_id])
            } else {
                paste0(source_gene_display, " -> ", unname(target_display_lookup[target_feature_id]))
            }
        )

    multiplicity <- mapped_tbl %>%
        filter(!is.na(target_feature_id)) %>%
        distinct(source_gene, target_feature_id) %>%
        count(source_gene, name = "mapped_gene_count") %>%
        filter(mapped_gene_count > 1)

    list(
        mapping = mapped_tbl,
        plot_table = plot_table,
        plot_features = plot_table$target_feature_id,
        label_map = setNames(plot_table$plot_label, plot_table$target_feature_id),
        no_orthogroup = sort(unique(no_orthogroup)),
        no_target_members = sort(unique(no_target_members)),
        missing_features = sort(unique(missing_features)),
        multiplicity = multiplicity
    )
}

within_group_choices <- function(obj) {
    available_cols <- colnames(obj@meta.data)
    choices <- c("Cluster" = "cluster_label")

    maybe_add <- function(label, column) {
        if (column %in% available_cols && !(column %in% unname(choices))) {
            choices <<- c(choices, setNames(column, label))
        }
    }

    maybe_add("Condition / group", "Group")
    maybe_add("Condition", "condition")
    maybe_add("Sample label", "Sample")
    maybe_add("Sample id", "sample")
    maybe_add("Batch", "batch")
    maybe_add("Integrated dataset", "integrated_dataset")
    maybe_add("Study", "study")
    maybe_add("Time point", "time_point")

    choices
}

within_split_choices <- function(obj) {
    available_cols <- colnames(obj@meta.data)
    choices <- c("No split" = "none")

    maybe_add <- function(label, column) {
        if (column %in% available_cols && !(column %in% unname(choices))) {
            choices <<- c(choices, setNames(column, label))
        }
    }

    maybe_add("Condition / group", "Group")
    maybe_add("Condition", "condition")
    maybe_add("Sample label", "Sample")
    maybe_add("Sample id", "sample")
    maybe_add("Batch", "batch")
    maybe_add("Integrated dataset", "integrated_dataset")
    maybe_add("Study", "study")
    maybe_add("Time point", "time_point")

    choices
}

cross_group_choices <- function(obj) {
    available_cols <- colnames(obj@meta.data)
    choices <- c(
        "Species + cell class" = "species_cell_class",
        "Cell class" = "cell_class",
        "Species" = "species",
        "Cluster" = "cluster_label"
    )

    maybe_add <- function(label, column) {
        if (column %in% available_cols && !(column %in% unname(choices))) {
            choices <<- c(choices, setNames(column, label))
        }
    }

    maybe_add("Condition", "condition")
    maybe_add("Study", "study")
    maybe_add("Time point", "time_point")

    choices
}

cross_split_choices <- function(obj) {
    available_cols <- colnames(obj@meta.data)
    choices <- c("No split" = "none")

    maybe_add <- function(label, column) {
        if (column %in% available_cols && !(column %in% unname(choices))) {
            choices <<- c(choices, setNames(column, label))
        }
    }

    maybe_add("Species", "species")
    maybe_add("Cell class", "cell_class")
    maybe_add("Condition", "condition")
    maybe_add("Study", "study")
    maybe_add("Time point", "time_point")

    choices
}

species_tab_ui <- function(species_key) {
    label <- species_label(species_key)

    tabPanel(
        title = species_label_tag(species_key),
        value = species_key,
        div(
            class = "section-card",
            div(
                class = "section-header",
                div(class = "section-eyebrow", species_label_tag(species_key)),
                h2(tagList(species_label_tag(species_key), " within-species explorer")),
                p("Selected source genes appear directly in the source-species tab and map automatically to orthologs when you view the other species.")
            ),
            fluidRow(
                column(
                    width = 4,
                    div(
                        class = "option-group",
                        selectInput(
                            inputId = paste0(species_key, "_integration"),
                            label = "Integration",
                            choices = integration_choices,
                            selected = "ComBat_BBKNN"
                        )
                    )
                ),
                column(
                    width = 4,
                    div(
                        class = "option-group",
                        uiOutput(paste0(species_key, "_distribution_group_by_ui"))
                    )
                ),
                column(
                    width = 4,
                    div(
                        class = "option-group",
                        uiOutput(paste0(species_key, "_distribution_split_by_ui"))
                    )
                )
            ),
            div(
                class = "subsection-header",
                h3("Cell distribution UMAP"),
                p("Inspect clusters, sample mixing, and metadata splits independently of gene expression.")
            ),
            fluidRow(
                column(
                    width = 4,
                    div(
                        class = "option-group",
                        sliderInput(
                            inputId = paste0(species_key, "_distribution_pt_size"),
                            label = "Cell distribution point size",
                            min = 0.1,
                            max = 2.5,
                            value = 0.55,
                            step = 0.05
                        )
                    )
                )
            ),
            fluidRow(
                column(
                    width = 12,
                    div(
                        class = "plot-card",
                        div(
                            class = "plot-card-header",
                            div(class = "plot-card-title", "Cell distribution"),
                            plot_download_button(paste0("dl_", species_key, "_distribution_umap"))
                        ),
                        spinning_plot_output(paste0(species_key, "_distribution_umap_plot"), proxy_height = "520px")
                    )
                )
            ),
            div(
                class = "subsection-header",
                h3("Gene expression"),
                p("Plot the selected source genes directly or through ortholog mappings in this species atlas.")
            ),
            fluidRow(
                column(
                    width = 4,
                    div(
                        class = "option-group",
                        uiOutput(paste0(species_key, "_split_by_ui"))
                    )
                ),
                column(
                    width = 4,
                    div(
                        class = "option-group",
                        uiOutput(paste0(species_key, "_group_by_ui"))
                    )
                ),
                column(
                    width = 4,
                    div(
                        class = "option-group",
                        sliderInput(
                            inputId = paste0(species_key, "_pt_size"),
                            label = "Feature UMAP point size",
                            min = 0.1,
                            max = 2.5,
                            value = 0.55,
                            step = 0.05
                        )
                    )
                )
            ),
            fluidRow(
                column(
                    width = 4,
                    div(
                        class = "option-group",
                        selectInput(
                            inputId = paste0(species_key, "_umap_columns"),
                            label = "Feature UMAP panels per row",
                            choices = c("1" = 1, "2" = 2, "3" = 3, "4" = 4),
                            selected = 1
                        )
                    )
                )
            ),
            uiOutput(paste0(species_key, "_notice_ui")),
            fluidRow(
                column(
                    width = 12,
                    div(
                        class = "plot-card",
                        div(
                            class = "plot-card-header",
                            div(class = "plot-card-title", "Feature UMAPs"),
                            plot_download_button(paste0("dl_", species_key, "_umap"))
                        ),
                        spinning_plot_output(paste0(species_key, "_umap_plot"), proxy_height = "520px")
                    )
                )
            ),
            fluidRow(
                column(
                    width = 7,
                    div(
                        class = "plot-card",
                        div(
                            class = "plot-card-header",
                            div(class = "plot-card-title", "Expression violin plots"),
                            plot_download_button(paste0("dl_", species_key, "_violin"))
                        ),
                        spinning_plot_output(paste0(species_key, "_violin_plot"), proxy_height = "520px")
                    )
                ),
                column(
                    width = 5,
                    div(
                        class = "plot-card",
                        div(class = "plot-card-title", "Top groups across mapped genes"),
                        spinning_plot_output(paste0(species_key, "_rank_bubble_plot"), proxy_height = "360px")
                    )
                )
            ),
            fluidRow(
                column(
                    width = 8,
                    div(
                        class = "plot-card",
                        div(
                            class = "plot-card-header",
                            div(class = "plot-card-title", "Multi-gene dot plot"),
                            plot_download_button(paste0("dl_", species_key, "_dot"))
                        ),
                        spinning_plot_output(paste0(species_key, "_dot_plot"), proxy_height = "420px")
                    )
                ),
                column(
                    width = 4,
                    div(
                        class = "table-card",
                        div(class = "plot-card-title", "Top groups per mapped gene"),
                        uiOutput(paste0(species_key, "_rank_table_ui"))
                    )
                )
            )
        )
    )
}

# ==============================================================================
# 2. USER INTERFACE
# ==============================================================================

ui <- fluidPage(
    tags$head(
        tags$link(rel = "preconnect", href = "https://fonts.googleapis.com"),
        tags$link(rel = "preconnect", href = "https://fonts.gstatic.com", crossorigin = "anonymous"),
        tags$link(
            rel = "stylesheet",
            href = "https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600;9..144,700&family=Montserrat:wght@400;500;600;700;800&display=swap"
        ),
        tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
    ),
        div(
            class = "app-shell",
            div(
                class = "busy-indicator",
                div(class = "busy-indicator-card",
                    div(class = "busy-indicator-spinner"),
                    div(class = "busy-indicator-copy",
                        div(class = "busy-indicator-text", "Updating plots"),
                        div(class = "busy-indicator-subtext", "Please wait while the atlas refreshes.")
                    )
                )
            ),
        div(
            class = "app-brand-shell",
            div(
                class = "app-brand-lockup",
                div(class = "app-brand-glyph", icon("leaf")),
                div(
                    class = "app-brand-text-wrap",
                    div(class = "app-brand-title", "Legume Root Nodule Symbiosis Atlas")
                )
            )
        ),
        div(
            class = "hero-shell",
                div(
                    class = "hero-copy",
                    div(class = "hero-eyebrow", "Single-cell explorer"),
                    h1("An scRNA-seq atlas for determinated and indeterminated nodules, including cross species comparison"),
                    p(
                        class = "hero-text",
                        HTML("Select a source species, search by gene ID or common name, and compare matched expression patterns across <em>Medicago truncatula</em>, <em>Glycine max</em>, <em>Lotus japonicus</em>, and the cross-species integration.")
                    )
                ),
            div(
                class = "summary-panel",
                div(class = "summary-panel-kicker", "Atlas summary"),
                uiOutput("atlas_summary_ui")
            )
        ),
        div(
            class = "section-card",
            div(
                class = "section-header",
                div(class = "section-eyebrow", "1. Configure the explorer"),
                h2("Choose the source gene panel"),
                p("The selected genes always stay synchronized across the species-specific tabs and the cross-species integration tab.")
            ),
            fluidRow(
                column(
                    width = 4,
                    div(
                        class = "option-group source-species-picker",
                        shinyWidgets::pickerInput(
                            inputId = "source_species",
                            label = "Source species",
                            choices = species_choices,
                            selected = "medicago",
                            multiple = FALSE
                        )
                    )
                ),
                column(
                    width = 6,
                    div(
                        class = "option-group",
                        shinyWidgets::pickerInput(
                            inputId = "selected_genes",
                            label = "Search by gene ID or common name",
                            choices = NULL,
                            multiple = TRUE,
                            options = list(
                                "actions-box" = TRUE,
                                "live-search" = TRUE,
                                "live-search-normalize" = TRUE,
                                "max-options" = 50,
                                "selected-text-format" = "count > 3",
                                "virtualScroll" = 30
                            )
                        ),
                        actionButton(
                            inputId = "apply_gene_selection",
                            label = "Apply gene panel",
                            icon = icon("play"),
                            class = "btn btn-default apply-selection-btn"
                        ),
                        div(class = "selection-meta", textOutput("gene_selection_status"))
                    )
                ),
                column(
                    width = 2,
                    div(
                        class = "option-group",
                        selectInput(
                            inputId = "dl_format",
                            label = "Download format",
                            choices = c("SVG", "PNG"),
                            selected = "PNG"
                        )
                    )
                )
            )
        ),
        tabsetPanel(
            id = "main_tabs",
            tabPanel(
                title = "Overview",
                value = "overview",
                div(
                    class = "section-card",
                    div(
                        class = "section-header",
                        div(class = "section-eyebrow", "2. Atlas overview"),
                        h2("Check the selected genes and their ortholog routing"),
                        p("Use this overview to see which selected genes have orthogroups and which targets are missing in each species atlas.")
                    ),
                    uiOutput("overview_alerts_ui"),
                    div(
                        class = "table-card",
                        div(class = "plot-card-title", "Ortholog mapping summary"),
                        uiOutput("overview_mapping_table_ui")
                    )
                )
            ),
            species_tab_ui("medicago"),
            species_tab_ui("glycine"),
            species_tab_ui("lotus"),
            tabPanel(
                title = "Cross-species integration",
                value = "cross_species",
                div(
                    class = "section-card",
                    div(
                        class = "section-header",
                        div(class = "section-eyebrow", "Cross-species integration"),
                        h2("Shared expression space across the three species"),
                        p("This tab uses the Camex integration object. All queries resolve to shared Medicago-space features before plots are generated.")
                    ),
                    fluidRow(
                        column(
                            width = 3,
                            div(
                                class = "option-group",
                                sliderInput(
                                    inputId = "cross_pt_size",
                                    label = "UMAP point size",
                                    min = 0.1,
                                    max = 2.5,
                                    value = 0.45,
                                    step = 0.05
                                )
                            )
                        ),
                        column(
                            width = 3,
                            div(
                                class = "option-group",
                                selectInput(
                                    inputId = "cross_umap_columns",
                                    label = "UMAP panels per row",
                                    choices = c("1" = 1, "2" = 2, "3" = 3, "4" = 4),
                                    selected = 1
                                )
                            )
                        ),
                        column(
                            width = 3,
                            div(
                                class = "option-group",
                                uiOutput("cross_split_by_ui")
                            )
                        ),
                        column(
                            width = 3,
                            div(
                                class = "option-group",
                                uiOutput("cross_group_by_ui")
                            )
                        )
                    ),
                    uiOutput("cross_notice_ui"),
                    fluidRow(
                        column(
                            width = 12,
                            div(
                                class = "plot-card",
                                div(
                                    class = "plot-card-header",
                                    div(class = "plot-card-title", "Integrated feature UMAPs"),
                                    plot_download_button("dl_cross_umap")
                                ),
                                spinning_plot_output("cross_umap_plot", proxy_height = "520px")
                            )
                        )
                    ),
                    fluidRow(
                        column(
                            width = 7,
                            div(
                                class = "plot-card",
                                div(
                                    class = "plot-card-header",
                                    div(class = "plot-card-title", "Cross-species dot plot"),
                                    plot_download_button("dl_cross_dot")
                                ),
                                spinning_plot_output("cross_dot_plot", proxy_height = "420px")
                            )
                        ),
                        column(
                            width = 5,
                            div(
                                class = "table-card",
                                div(class = "plot-card-title", "Cross-space mapped features"),
                                uiOutput("cross_mapping_table_ui")
                            )
                        )
                    )
                )
            )
        ),
        tags$footer(
            "Single-cell atlas explorer for legume root nodule symbiosis.",
            strong(" 2026."),
            class = "app-footer"
        )
    )
)

# ==============================================================================
# 3. SERVER LOGIC
# ==============================================================================

server <- function(input, output, session) {
    selection_notice <- reactiveVal(NULL)
    applied_selected_genes <- reactiveVal(character(0))

    current_species_integration <- function(species_key) {
        input[[paste0(species_key, "_integration")]] %||% "ComBat_BBKNN"
    }

    staged_source_genes <- reactive({
        genes <- input$selected_genes %||% character(0)
        unique(as.character(genes))
    })

    selected_source_genes <- reactive({
        applied_selected_genes()
    })

    source_integration <- reactive({
        current_species_integration(input$source_species %||% "medicago")
    })

    source_gene_catalog <- reactive({
        build_gene_catalog(
            input$source_species %||% "medicago",
            source_integration()
        )
    })

    observeEvent(
        source_gene_catalog(),
        {
            choice_bundle <- build_gene_choices(
                input$source_species %||% "medicago",
                source_integration()
            )

            current_selection <- isolate(staged_source_genes())
            current_applied_selection <- isolate(applied_selected_genes())
            valid_selection <- intersect(current_selection, choice_bundle$feature_ids)
            valid_applied_selection <- intersect(current_applied_selection, choice_bundle$feature_ids)
            dropped_genes <- unique(c(
                setdiff(current_selection, valid_selection),
                setdiff(current_applied_selection, valid_applied_selection)
            ))

            if (length(dropped_genes)) {
                selection_notice(
                    paste0(
                        "Dropped genes not present in the new source atlas: ",
                        compact_gene_list(dropped_genes, limit = 6)
                    )
                )
            } else {
                selection_notice(NULL)
            }

            applied_selected_genes(valid_applied_selection)

            shinyWidgets::updatePickerInput(
                session = session,
                inputId = "selected_genes",
                choices = choice_bundle$choices,
                selected = valid_selection
            )
        },
        ignoreInit = FALSE
    )

    observeEvent(
        input$apply_gene_selection,
        {
            applied_selected_genes(staged_source_genes())
        }
    )

    source_orthogroup_status <- reactive({
        genes <- selected_source_genes()

        if (!length(genes)) {
            return(list(
                with_orthogroup = character(0),
                without_orthogroup = character(0)
            ))
        }

        source_orthogroups <- resolve_source_orthogroups(
            input$source_species %||% "medicago",
            genes
        )

        genes_with_orthogroup <- source_orthogroups %>%
            filter(!is.na(orthogroup)) %>%
            distinct(source_gene) %>%
            pull(source_gene)

        list(
            with_orthogroup = genes_with_orthogroup,
            without_orthogroup = setdiff(genes, genes_with_orthogroup)
        )
    })

    cross_resolution <- reactive({
        resolve_target_mapping(
            source_species = input$source_species %||% "medicago",
            source_genes = selected_source_genes(),
            target_species = "medicago",
            integration_method = current_species_integration("medicago"),
            cross_space = TRUE
        )
    })

    output$gene_selection_status <- renderText({
        staged_genes <- staged_source_genes()
        applied_genes <- selected_source_genes()
        n_genes <- length(applied_genes)
        note <- selection_notice()

        if (!same_gene_selection(staged_genes, applied_genes)) {
            staged_message <- if (length(staged_genes)) {
                paste0(
                    length(staged_genes),
                    " gene(s) staged. Click Apply gene panel to refresh the plots."
                )
            } else {
                "No genes staged. Click Apply gene panel to clear the plots."
            }

            if (!is.null(note)) {
                return(paste(staged_message, note))
            }

            return(staged_message)
        }

        if (n_genes == 0) {
            if (!is.null(note)) {
                return(paste("No genes selected yet.", note))
            }

            return("No genes selected yet.")
        }

        base_message <- paste0(
            n_genes,
            " source gene(s) selected. Orthologs sync automatically across the other tabs."
        )

        if (!is.null(note)) {
            paste(base_message, note)
        } else {
            base_message
        }
    })

    output$atlas_summary_ui <- renderUI({
        summary_cards <- lapply(within_species_keys, function(species_key) {
            summary <- get_within_dataset_summary(
                species_key,
                current_species_integration(species_key)
            )

            div(
                class = "dataset-tile",
                div(class = "dataset-title", species_label_tag(species_key)),
                div(class = "dataset-value", paste(format_stat_value(summary$cells), "cells")),
                div(class = "dataset-note", paste(format_stat_value(summary$genes), "expressed genes")),
                div(class = "dataset-note", paste(format_stat_value(summary$sample_n), "samples"))
            )
        })

        cross_summary <- get_cross_dataset_summary()
        cross_card <- div(
            class = "dataset-tile",
            div(class = "dataset-title", "Cross-species integration"),
            div(class = "dataset-value", paste(format_stat_value(cross_summary$cells), "cells")),
            div(class = "dataset-note", paste(format_stat_value(cross_summary$genes), "expressed genes")),
            div(class = "dataset-note", paste(format_stat_value(cross_summary$sample_n), "samples"))
        )

        div(
            class = "dataset-grid atlas-summary-grid",
            tagList(summary_cards, cross_card)
        )
    })

    overview_mapping_table <- reactive({
        genes <- selected_source_genes()

        if (!length(genes)) {
            return(tibble())
        }

        source_species <- input$source_species %||% "medicago"
        source_orthogroups <- resolve_source_orthogroups(source_species, genes) %>%
            distinct(source_gene, orthogroup, .keep_all = TRUE) %>%
            arrange(match(source_gene, genes))

        source_orthogroups %>%
            mutate(
                `Source species` = species_label(source_species),
                `Source gene` = display_gene_labels(source_species, source_gene),
                Orthogroup = ifelse(is.na(orthogroup), "No orthogroup", orthogroup),
                Medicago = map_chr(orthogroup, ~ compact_gene_list(get_orthogroup_members(.x, "medicago"))),
                `Glycine max` = map_chr(orthogroup, ~ compact_gene_list(get_orthogroup_members(.x, "glycine"))),
                `Lotus japonicus` = map_chr(orthogroup, ~ compact_gene_list(get_orthogroup_members(.x, "lotus")))
            ) %>%
            select(
                `Source species`,
                `Source gene`,
                Orthogroup,
                Medicago,
                `Glycine max`,
                `Lotus japonicus`
            )
    })

    output$overview_mapping_table_ui <- renderUI({
        mapping_tbl <- overview_mapping_table()

        if (!nrow(mapping_tbl)) {
            return(div(
                class = "summary-placeholder",
                "Select one or more source-species genes to inspect the ortholog mapping summary."
            ))
        }

        html_summary_table(mapping_tbl)
    })

    output$overview_alerts_ui <- renderUI({
        genes <- selected_source_genes()

        if (!length(genes)) {
            return(
                div(
                    class = "alert-stack",
                    notice_card(
                        title = "Ready for selection",
                        body = "Choose a source species and add genes from that species to populate the within-species and cross-species panels.",
                        tone = "info"
                    )
                )
            )
        }

        source_species <- input$source_species %||% "medicago"
        orthogroup_status <- source_orthogroup_status()
        other_species <- setdiff(within_species_keys, source_species)

        cards <- list()

        if (length(orthogroup_status$without_orthogroup)) {
            cards <- append(cards, list(
                notice_card(
                    title = "Selected genes without orthogroups",
                    body = compact_gene_list(orthogroup_status$without_orthogroup, limit = 8),
                    tone = "warning"
                )
            ))
        }

        for (target_species in other_species) {
            resolution <- resolve_target_mapping(
                source_species = source_species,
                source_genes = genes,
                target_species = target_species,
                integration_method = current_species_integration(target_species),
                cross_space = FALSE
            )

            if (length(resolution$no_target_members)) {
                cards <- append(cards, list(
                    notice_card(
                        title = paste(species_label(target_species), "orthogroup gaps"),
                        body = paste(
                            "Orthogroups were found, but they do not contain target-species members for:",
                            compact_gene_list(resolution$no_target_members, limit = 8)
                        ),
                        tone = "warning"
                    )
                ))
            }

            if (length(resolution$missing_features)) {
                cards <- append(cards, list(
                    notice_card(
                        title = paste(species_label(target_species), "feature gaps"),
                        body = paste(
                            "Orthologs were found, but no matching features are present in the selected atlas for:",
                            compact_gene_list(resolution$missing_features, limit = 8)
                        ),
                        tone = "warning"
                    )
                ))
            }

        }

        cross_res <- cross_resolution()

        if (length(cross_res$no_target_members)) {
            cards <- append(cards, list(
                notice_card(
                    title = "Cross-species Medicago-space gaps",
                    body = paste(
                        "The orthogroups below do not contain Medicago members for the cross-species space:",
                        compact_gene_list(cross_res$no_target_members, limit = 8)
                    ),
                    tone = "warning"
                )
            ))
        }

        if (!length(cards)) {
            cards <- list(
                notice_card(
                    title = "Mappings look consistent",
                    body = "The current gene panel resolves cleanly across the species-specific atlases and the cross-species Medicago-space feature set.",
                    tone = "info"
                )
            )
        }

        div(class = "alert-stack", tagList(cards))
    })

    get_ext <- function() {
        tolower(input$dl_format %||% "svg")
    }

    save_ggplot <- function(file, plot_obj, width, height) {
        ggsave(
            filename = file,
            plot = plot_obj,
            device = get_ext(),
            width = width,
            height = height,
            dpi = 300,
            limitsize = FALSE
        )
    }

    register_species_tab <- function(species_key) {
        local({
            prefix <- species_key
            tab_label <- species_label(species_key)

            tab_object <- reactive({
                get_within_object(
                    species_key,
                    current_species_integration(species_key)
                )
            })

            tab_resolution <- reactive({
                resolve_target_mapping(
                    source_species = input$source_species %||% "medicago",
                    source_genes = selected_source_genes(),
                    target_species = species_key,
                    integration_method = current_species_integration(species_key),
                    cross_space = FALSE
                )
            })

            tab_group_choices <- reactive(within_group_choices(tab_object()))
            tab_split_choices <- reactive(within_split_choices(tab_object()))

            output[[paste0(prefix, "_distribution_group_by_ui")]] <- renderUI({
                choices <- tab_group_choices()

                selectInput(
                    inputId = paste0(prefix, "_distribution_group_by"),
                    label = "Color cells by",
                    choices = choices,
                    selected = resolve_choice(
                        input[[paste0(prefix, "_distribution_group_by")]],
                        choices,
                        default = "cluster_label"
                    )
                )
            })

            output[[paste0(prefix, "_distribution_split_by_ui")]] <- renderUI({
                choices <- tab_split_choices()

                selectInput(
                    inputId = paste0(prefix, "_distribution_split_by"),
                    label = "Split UMAP by",
                    choices = choices,
                    selected = resolve_choice(
                        input[[paste0(prefix, "_distribution_split_by")]],
                        choices,
                        default = "none"
                    )
                )
            })

            output[[paste0(prefix, "_group_by_ui")]] <- renderUI({
                choices <- tab_group_choices()

                selectInput(
                    inputId = paste0(prefix, "_group_by"),
                    label = "Summarize expression by",
                    choices = choices,
                    selected = resolve_choice(
                        input[[paste0(prefix, "_group_by")]],
                        choices,
                        default = "cluster_label"
                    )
                )
            })

            output[[paste0(prefix, "_split_by_ui")]] <- renderUI({
                choices <- tab_split_choices()

                selectInput(
                    inputId = paste0(prefix, "_split_by"),
                    label = "Split feature UMAP by",
                    choices = choices,
                    selected = resolve_choice(
                        input[[paste0(prefix, "_split_by")]],
                        choices,
                        default = "none"
                    )
                )
            })

            output[[paste0(prefix, "_notice_ui")]] <- renderUI({
                genes <- selected_source_genes()
                resolution <- tab_resolution()
                source_species <- input$source_species %||% "medicago"
                cards <- list()

                if (!length(genes)) {
                    return(
                        div(
                            class = "alert-stack",
                            notice_card(
                                title = paste("No source genes selected for", tab_label, "expression"),
                                body = "Add one or more source-species genes to populate the expression panels in this tab. The cell distribution UMAP above remains available without gene selection.",
                                tone = "info"
                            )
                        )
                    )
                }

                if (identical(source_species, species_key)) {
                    cards <- append(cards, list(
                        notice_card(
                            title = "Source-species view",
                            body = sprintf("Showing the selected source genes directly in the expression panels for %s.", tab_label),
                            tone = "info"
                        )
                    ))
                } else {
                    cards <- append(cards, list(
                        notice_card(
                            title = "Ortholog-mapped view",
                            body = sprintf(
                                "Using orthogroup mappings from %s into the expression panels for %s.",
                                species_label(source_species),
                                tab_label
                            ),
                            tone = "info"
                        )
                    ))
                }

                if (length(resolution$no_orthogroup)) {
                    cards <- append(cards, list(
                        notice_card(
                            title = "Selected genes without orthogroups",
                            body = compact_gene_list(resolution$no_orthogroup, limit = 8),
                            tone = "warning"
                        )
                    ))
                }

                if (length(resolution$no_target_members)) {
                    cards <- append(cards, list(
                        notice_card(
                            title = paste("No", tab_label, "members in the orthogroups"),
                            body = compact_gene_list(resolution$no_target_members, limit = 8),
                            tone = "warning"
                        )
                    ))
                }

                if (length(resolution$missing_features)) {
                    cards <- append(cards, list(
                        notice_card(
                            title = "Mapped orthologs missing from the selected atlas",
                            body = compact_gene_list(resolution$missing_features, limit = 8),
                            tone = "warning"
                        )
                    ))
                }

                if (nrow(resolution$multiplicity)) {
                    multiplicity_text <- resolution$multiplicity %>%
                        mutate(label = paste0(source_gene, " (", mapped_gene_count, " orthologs)")) %>%
                        pull(label)

                    cards <- append(cards, list(
                        notice_card(
                            title = "One-to-many orthologs",
                            body = compact_gene_list(multiplicity_text, limit = 6),
                            tone = "info"
                        )
                    ))
                }

                div(class = "alert-stack", tagList(cards))
            })

            tab_distribution_group_by <- reactive({
                resolve_choice(
                    input[[paste0(prefix, "_distribution_group_by")]],
                    tab_group_choices(),
                    default = "cluster_label"
                )
            })

            tab_distribution_split_by <- reactive({
                resolve_choice(
                    input[[paste0(prefix, "_distribution_split_by")]],
                    tab_split_choices(),
                    default = "none"
                )
            })

            tab_group_by <- reactive({
                resolve_choice(
                    input[[paste0(prefix, "_group_by")]],
                    tab_group_choices(),
                    default = "cluster_label"
                )
            })

            tab_split_by <- reactive({
                resolve_choice(
                    input[[paste0(prefix, "_split_by")]],
                    tab_split_choices(),
                    default = "none"
                )
            })

            distribution_umap_plot_obj <- reactive({
                obj <- tab_object()
                group_by <- tab_distribution_group_by()
                split_by <- tab_distribution_split_by()
                pt_size <- as.numeric(input[[paste0(prefix, "_distribution_pt_size")]] %||% 0.55)
                split_panels <- if (identical(split_by, "none")) 1L else split_panel_count(obj, split_by)
                split_columns <- if (identical(split_by, "none")) {
                    NULL
                } else if (split_panels <= 4L) {
                    2L
                } else {
                    3L
                }

                distribution_plot <- scCustomize::DimPlot_scCustom(
                    seurat_object = obj,
                    group.by = group_by,
                    split.by = if (identical(split_by, "none")) NULL else split_by,
                    pt.size = pt_size,
                    label = identical(group_by, "cluster_label") && identical(split_by, "none"),
                    repel = TRUE,
                    raster = TRUE,
                    num_columns = split_columns
                )

                distribution_plot &
                    labs(title = NULL, color = NULL) &
                    app_plot_theme() &
                    theme(
                        legend.title = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        axis.line = element_blank(),
                        plot.margin = margin(8, 14, 10, 10)
                    )
            })

            rank_data <- reactive({
                resolution <- tab_resolution()
                obj <- tab_object()

                validate(
                    need(
                        length(resolution$plot_features) > 0,
                        paste("No mapped genes are available for", tab_label, "in the selected atlas.")
                    )
                )

                dot_data <- Seurat::DotPlot(
                    object = obj,
                    features = resolution$plot_features,
                    group.by = tab_group_by()
                )$data

                label_map <- resolution$label_map

                as_tibble(dot_data) %>%
                    transmute(
                        gene_id = features.plot,
                        Gene = unname(label_map[features.plot]),
                        Group = id,
                        `Pct. expressing` = round(pct.exp, 1),
                        `Scaled expression` = round(avg.exp.scaled, 2),
                        `Average expression` = round(avg.exp, 3)
                    ) %>%
                    group_by(Gene) %>%
                    slice_max(order_by = `Scaled expression`, n = 5, with_ties = FALSE) %>%
                    ungroup()
            })

            umap_plot_obj <- reactive({
                resolution <- tab_resolution()
                obj <- tab_object()
                split_by <- tab_split_by()
                requested_cols <- as.integer(input[[paste0(prefix, "_umap_columns")]] %||% 1)
                pt_size <- as.numeric(input[[paste0(prefix, "_pt_size")]] %||% 0.55)
                feature_grid_cols <- resolve_feature_grid_cols(
                    requested_cols = requested_cols,
                    feature_labels = unname(resolution$label_map[resolution$plot_features]),
                    split_by = split_by
                )

                validate(
                    need(
                        length(resolution$plot_features) > 0,
                        paste("No mapped genes are available for", tab_label, "in the selected atlas.")
                    )
                )

                plot_list <- lapply(resolution$plot_features, function(feature_id) {
                    feature_plot <- scCustomize::FeaturePlot_scCustom(
                        seurat_object = obj,
                        features = feature_id,
                        split.by = if (identical(split_by, "none")) NULL else split_by,
                        pt.size = pt_size,
                        order = TRUE,
                        raster = TRUE,
                        num_columns = if (identical(split_by, "none")) NULL else requested_cols,
                        combine = TRUE
                    )

                    feature_plot <- feature_plot &
                        labs(title = NULL, color = NULL) &
                        app_plot_theme() &
                        theme(
                            plot.title = element_blank(),
                            legend.title = element_blank(),
                            axis.title = element_blank(),
                            axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            axis.line = element_blank(),
                            plot.margin = margin(8, 14, 10, 10)
                        )

                    wrap_titled_plot(
                        plot_obj = feature_plot,
                        title = unname(resolution$label_map[feature_id])
                    )
                })

                wrap_plots(plotlist = plot_list, ncol = feature_grid_cols)
            })

            violin_plot_obj <- reactive({
                resolution <- tab_resolution()
                obj <- tab_object()

                validate(
                    need(
                        length(resolution$plot_features) > 0,
                        paste("No mapped genes are available for", tab_label, "in the selected atlas.")
                    )
                )

                label_map <- resolution$label_map
                n_groups <- dplyr::n_distinct(obj@meta.data[[tab_group_by()]])

                plot_list <- lapply(resolution$plot_features, function(feature_id) {
                    violin_plot <- scCustomize::VlnPlot_scCustom(
                        seurat_object = obj,
                        features = feature_id,
                        group.by = tab_group_by(),
                        colors_use = rep(unname(app_palette["warm"]), n_groups),
                        pt.size = 0.02,
                        num_columns = 1,
                        raster = TRUE
                    )

                    violin_plot <- violin_plot &
                        labs(title = NULL, x = NULL, y = "Normalized expression") &
                        app_plot_theme() &
                        theme(
                            plot.title = element_blank(),
                            legend.position = "none",
                            axis.text.x = element_text(angle = 35, hjust = 1),
                            panel.grid.major.x = element_blank()
                        )

                    wrap_titled_plot(
                        plot_obj = violin_plot,
                        title = unname(label_map[feature_id])
                    )
                })

                wrap_plots(plotlist = plot_list, ncol = 1)
            })

            dot_plot_obj <- reactive({
                resolution <- tab_resolution()
                obj <- tab_object()

                validate(
                    need(
                        length(resolution$plot_features) > 0,
                        paste("No mapped genes are available for", tab_label, "in the selected atlas.")
                    )
                )

                scCustomize::DotPlot_scCustom(
                    seurat_object = obj,
                    features = resolution$plot_features,
                    group.by = tab_group_by()
                ) +
                    scale_x_discrete(labels = resolution$label_map) +
                    app_plot_theme() +
                    labs(x = NULL, y = NULL, color = "Scaled expression", size = "% expressing") +
                    theme(
                        panel.grid.major.y = element_blank(),
                        axis.text.x = element_text(angle = 35, hjust = 1)
                    )
            })

            output[[paste0(prefix, "_umap_plot")]] <- renderPlot(
                {
                    umap_plot_obj()
                },
                height = function() {
                    feature_n <- tryCatch(length(tab_resolution()$plot_features), error = function(e) 0L)
                    split_by <- tryCatch(tab_split_by(), error = function(e) "none")
                    feature_cols <- tryCatch(
                        resolve_feature_grid_cols(
                            requested_cols = as.integer(input[[paste0(prefix, "_umap_columns")]] %||% 1),
                            feature_labels = unname(tab_resolution()$label_map[tab_resolution()$plot_features]),
                            split_by = split_by
                        ),
                        error = function(e) 1L
                    )

                    if (identical(split_by, "none")) {
                        feature_rows <- ceiling(feature_n / max(1L, feature_cols))
                        max(520, feature_rows * 360)
                    } else {
                        panels_per_gene <- tryCatch(
                            split_panel_count(tab_object(), split_by),
                            error = function(e) 1L
                        )
                        rows_per_gene <- ceiling(panels_per_gene / max(1L, feature_cols))

                        max(520, feature_n * (rows_per_gene * 290 + 90))
                    }
                },
                res = 110
            )

            output[[paste0(prefix, "_distribution_umap_plot")]] <- renderPlot(
                {
                    distribution_umap_plot_obj()
                },
                height = function() {
                    split_by <- tryCatch(tab_distribution_split_by(), error = function(e) "none")

                    if (identical(split_by, "none")) {
                        return(520)
                    }

                    panel_n <- tryCatch(
                        split_panel_count(tab_object(), split_by),
                        error = function(e) 1L
                    )
                    split_cols <- if (panel_n <= 4L) 2L else 3L
                    split_rows <- ceiling(panel_n / split_cols)

                    max(520, split_rows * 320)
                },
                res = 110
            )

            output[[paste0(prefix, "_violin_plot")]] <- renderPlot(
                {
                    violin_plot_obj()
                },
                height = function() {
                    feature_n <- tryCatch(length(tab_resolution()$plot_features), error = function(e) 0L)
                    max(420, 280 * feature_n)
                },
                res = 110
            )

            output[[paste0(prefix, "_rank_bubble_plot")]] <- renderPlot(
                {
                    data <- rank_data()

                    ggplot(
                        data,
                        aes(
                            x = `Scaled expression`,
                            y = reorder_within(Group, `Scaled expression`, Gene),
                            size = `Pct. expressing`,
                            colour = `Scaled expression`
                        )
                    ) +
                        geom_point(alpha = 0.92) +
                        facet_wrap(~Gene, scales = "free_y", ncol = 1) +
                        scale_y_reordered() +
                        scale_size_continuous(range = c(3.5, 12)) +
                        scale_colour_gradient2(
                            low = app_palette["green_soft"],
                            mid = app_palette["warm_soft"],
                            high = app_palette["green_dark"],
                            midpoint = 0
                        ) +
                        labs(
                            x = "Scaled average expression",
                            y = NULL,
                            size = "% expressing",
                            colour = "Scaled expression"
                        ) +
                        app_plot_theme() +
                        theme(panel.grid.major.y = element_blank())
                },
                height = function() {
                    feature_n <- tryCatch(length(tab_resolution()$plot_features), error = function(e) 0L)
                    max(360, 150 * feature_n)
                },
                res = 110
            )

            output[[paste0(prefix, "_rank_table_ui")]] <- renderUI({
                genes <- selected_source_genes()

                if (!length(genes)) {
                    return(div(
                        class = "summary-placeholder",
                        sprintf("Add source-species genes to see the top groups in %s.", tab_label)
                    ))
                }

                html_summary_table(rank_data())
            })

            output[[paste0(prefix, "_dot_plot")]] <- renderPlot(
                {
                    dot_plot_obj()
                },
                height = function() {
                    feature_n <- tryCatch(length(tab_resolution()$plot_features), error = function(e) 0L)
                    max(420, 95 * feature_n + 120)
                },
                res = 110
            )

            output[[paste0("dl_", prefix, "_umap")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_umap.", get_ext())
                },
                content = function(file) {
                    feature_n <- length(tab_resolution()$plot_features)
                    split_by <- tab_split_by()
                    feature_cols <- resolve_feature_grid_cols(
                        requested_cols = as.integer(input[[paste0(prefix, "_umap_columns")]] %||% 1),
                        feature_labels = unname(tab_resolution()$label_map[tab_resolution()$plot_features]),
                        split_by = split_by
                    )

                    plot_height <- if (identical(split_by, "none")) {
                        feature_rows <- ceiling(feature_n / max(1L, feature_cols))
                        max(7, feature_rows * 3.8)
                    } else {
                        panels_per_gene <- split_panel_count(tab_object(), split_by)
                        rows_per_gene <- ceiling(panels_per_gene / max(1L, feature_cols))
                        max(7, feature_n * (rows_per_gene * 2.9 + 1))
                    }

                    save_ggplot(
                        file = file,
                        plot_obj = umap_plot_obj(),
                        width = 15,
                        height = plot_height
                    )
                }
            )

            output[[paste0("dl_", prefix, "_distribution_umap")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_distribution_umap.", get_ext())
                },
                content = function(file) {
                    split_by <- tab_distribution_split_by()

                    plot_height <- if (identical(split_by, "none")) {
                        7
                    } else {
                        panel_n <- split_panel_count(tab_object(), split_by)
                        split_cols <- if (panel_n <= 4L) 2L else 3L
                        split_rows <- ceiling(panel_n / split_cols)
                        max(7, split_rows * 3.4)
                    }

                    save_ggplot(
                        file = file,
                        plot_obj = distribution_umap_plot_obj(),
                        width = 14,
                        height = plot_height
                    )
                }
            )

            output[[paste0("dl_", prefix, "_violin")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_violins.", get_ext())
                },
                content = function(file) {
                    feature_n <- length(tab_resolution()$plot_features)

                    save_ggplot(
                        file = file,
                        plot_obj = violin_plot_obj(),
                        width = 10,
                        height = max(6, feature_n * 3.2)
                    )
                }
            )

            output[[paste0("dl_", prefix, "_dot")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_dotplot.", get_ext())
                },
                content = function(file) {
                    feature_n <- length(tab_resolution()$plot_features)

                    save_ggplot(
                        file = file,
                        plot_obj = dot_plot_obj(),
                        width = 11,
                        height = max(6, feature_n * 0.85 + 2)
                    )
                }
            )
        })
    }

    walk(within_species_keys, register_species_tab)

    cross_object <- reactive({
        get_cross_object()
    })

    output$cross_group_by_ui <- renderUI({
        choices <- cross_group_choices(cross_object())

        selectInput(
            inputId = "cross_group_by",
            label = "Summarize expression by",
            choices = choices,
            selected = resolve_choice(
                input$cross_group_by,
                choices,
                default = "species_cell_class"
            )
        )
    })

    output$cross_split_by_ui <- renderUI({
        choices <- cross_split_choices(cross_object())

        selectInput(
            inputId = "cross_split_by",
            label = "Split UMAP by",
            choices = choices,
            selected = resolve_choice(
                input$cross_split_by,
                choices,
                default = "species"
            )
        )
    })

    cross_group_by <- reactive({
        resolve_choice(
            input$cross_group_by,
            cross_group_choices(cross_object()),
            default = "species_cell_class"
        )
    })

    cross_split_by <- reactive({
        resolve_choice(
            input$cross_split_by,
            cross_split_choices(cross_object()),
            default = "species"
        )
    })

    output$cross_notice_ui <- renderUI({
        genes <- selected_source_genes()
        resolution <- cross_resolution()
        cards <- list(
            notice_card(
                title = "Shared Medicago-space features",
                body = "The Camex object stores a shared feature space represented with Medicago identifiers. Soybean and Lotus selections are therefore projected to Medicago orthologs in this tab.",
                tone = "info"
            )
        )

        if (!length(genes)) {
            cards <- append(cards, list(
                notice_card(
                    title = "No source genes selected",
                    body = "Add one or more source-species genes to populate the cross-species integration plots.",
                    tone = "info"
                )
            ))
        }

        if (length(resolution$no_orthogroup)) {
            cards <- append(cards, list(
                notice_card(
                    title = "Selected genes without orthogroups",
                    body = compact_gene_list(resolution$no_orthogroup, limit = 8),
                    tone = "warning"
                )
            ))
        }

        if (length(resolution$no_target_members)) {
            cards <- append(cards, list(
                notice_card(
                    title = "Orthogroups without Medicago members",
                    body = compact_gene_list(resolution$no_target_members, limit = 8),
                    tone = "warning"
                )
            ))
        }

        if (length(resolution$missing_features)) {
            cards <- append(cards, list(
                notice_card(
                    title = "Mapped orthologs missing from the Camex feature set",
                    body = compact_gene_list(resolution$missing_features, limit = 8),
                    tone = "warning"
                )
            ))
        }

        if (nrow(resolution$multiplicity)) {
            multiplicity_text <- resolution$multiplicity %>%
                mutate(label = paste0(source_gene, " (", mapped_gene_count, " Medicago-space features)")) %>%
                pull(label)

            cards <- append(cards, list(
                notice_card(
                    title = "One-to-many cross-space mappings",
                    body = compact_gene_list(multiplicity_text, limit = 6),
                    tone = "info"
                )
            ))
        }

        div(class = "alert-stack", tagList(cards))
    })

    cross_mapping_table <- reactive({
        resolution <- cross_resolution()

        if (!nrow(resolution$plot_table)) {
            return(tibble())
        }

        resolution$plot_table %>%
            transmute(
                `Source gene(s)` = source_gene_label,
                `Medicago-space feature` = target_feature_id,
                Orthogroup = ifelse(nzchar(orthogroup_label), orthogroup_label, "NA")
            )
    })

    output$cross_mapping_table_ui <- renderUI({
        genes <- selected_source_genes()

        if (!length(genes)) {
            return(div(
                class = "summary-placeholder",
                "Add source-species genes to see which Medicago-space features are used in the Camex tab."
            ))
        }

        mapping_tbl <- cross_mapping_table()

        if (!nrow(mapping_tbl)) {
            return(div(
                class = "summary-placeholder",
                "No selected genes resolve to the shared Medicago-space feature set."
            ))
        }

        html_summary_table(mapping_tbl)
    })

    cross_umap_plot_obj <- reactive({
        resolution <- cross_resolution()
        obj <- cross_object()
        split_by <- cross_split_by()
        requested_cols <- as.integer(input$cross_umap_columns %||% 1)
        pt_size <- as.numeric(input$cross_pt_size %||% 0.45)
        feature_grid_cols <- resolve_feature_grid_cols(
            requested_cols = requested_cols,
            feature_labels = unname(resolution$label_map[resolution$plot_features]),
            split_by = split_by
        )

        validate(
            need(
                length(resolution$plot_features) > 0,
                "No selected genes resolve to Medicago-space features in the Camex integration."
            )
        )

        plot_list <- lapply(resolution$plot_features, function(feature_id) {
            feature_plot <- scCustomize::FeaturePlot_scCustom(
                seurat_object = obj,
                features = feature_id,
                split.by = if (identical(split_by, "none")) NULL else split_by,
                pt.size = pt_size,
                order = TRUE,
                raster = TRUE,
                num_columns = if (identical(split_by, "none")) NULL else requested_cols,
                combine = TRUE
            )

            feature_plot <- feature_plot &
                labs(title = NULL, color = NULL) &
                app_plot_theme() &
                theme(
                    plot.title = element_blank(),
                    legend.title = element_blank(),
                    axis.title = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.line = element_blank(),
                    plot.margin = margin(8, 14, 10, 10)
                )

            wrap_titled_plot(
                plot_obj = feature_plot,
                title = unname(resolution$label_map[feature_id])
            )
        })

        wrap_plots(plotlist = plot_list, ncol = feature_grid_cols)
    })

    cross_dot_plot_obj <- reactive({
        resolution <- cross_resolution()
        obj <- cross_object()

        validate(
            need(
                length(resolution$plot_features) > 0,
                "No selected genes resolve to Medicago-space features in the Camex integration."
            )
        )

        scCustomize::DotPlot_scCustom(
            seurat_object = obj,
            features = resolution$plot_features,
            group.by = cross_group_by()
        ) +
            scale_x_discrete(labels = resolution$label_map) +
            app_plot_theme() +
            labs(x = NULL, y = NULL, color = "Scaled expression", size = "% expressing") +
            theme(
                panel.grid.major.y = element_blank(),
                axis.text.x = element_text(angle = 35, hjust = 1)
            )
    })

    output$cross_umap_plot <- renderPlot(
        {
            cross_umap_plot_obj()
        },
        height = function() {
            feature_n <- tryCatch(length(cross_resolution()$plot_features), error = function(e) 0L)
            split_by <- tryCatch(cross_split_by(), error = function(e) "none")
            feature_cols <- tryCatch(
                resolve_feature_grid_cols(
                    requested_cols = as.integer(input$cross_umap_columns %||% 1),
                    feature_labels = unname(cross_resolution()$label_map[cross_resolution()$plot_features]),
                    split_by = split_by
                ),
                error = function(e) 1L
            )

            if (identical(split_by, "none")) {
                feature_rows <- ceiling(feature_n / max(1L, feature_cols))
                max(520, feature_rows * 360)
            } else {
                panels_per_gene <- tryCatch(
                    split_panel_count(cross_object(), split_by),
                    error = function(e) 1L
                )
                rows_per_gene <- ceiling(panels_per_gene / max(1L, feature_cols))

                max(520, feature_n * (rows_per_gene * 290 + 90))
            }
        },
        res = 110
    )

    output$cross_dot_plot <- renderPlot(
        {
            cross_dot_plot_obj()
        },
        height = function() {
            feature_n <- tryCatch(length(cross_resolution()$plot_features), error = function(e) 0L)
            max(420, 95 * feature_n + 120)
        },
        res = 110
    )

    output$dl_cross_umap <- downloadHandler(
        filename = function() {
            paste0("cross_species_umap.", get_ext())
        },
        content = function(file) {
            feature_n <- length(cross_resolution()$plot_features)
            split_by <- cross_split_by()
            feature_cols <- resolve_feature_grid_cols(
                requested_cols = as.integer(input$cross_umap_columns %||% 1),
                feature_labels = unname(cross_resolution()$label_map[cross_resolution()$plot_features]),
                split_by = split_by
            )

            plot_height <- if (identical(split_by, "none")) {
                feature_rows <- ceiling(feature_n / max(1L, feature_cols))
                max(7, feature_rows * 3.8)
            } else {
                panels_per_gene <- split_panel_count(cross_object(), split_by)
                rows_per_gene <- ceiling(panels_per_gene / max(1L, feature_cols))
                max(7, feature_n * (rows_per_gene * 2.9 + 1))
            }

            save_ggplot(
                file = file,
                plot_obj = cross_umap_plot_obj(),
                width = 15,
                height = plot_height
            )
        }
    )

    output$dl_cross_dot <- downloadHandler(
        filename = function() {
            paste0("cross_species_dotplot.", get_ext())
        },
        content = function(file) {
            feature_n <- length(cross_resolution()$plot_features)

            save_ggplot(
                file = file,
                plot_obj = cross_dot_plot_obj(),
                width = 12,
                height = max(6, feature_n * 0.85 + 2)
            )
        }
    )
}

shinyApp(ui = ui, server = server)
