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

atlas_summary_path <- "metadata/atlas_summary.tsv"
gene_catalog_cache_dir <- "metadata/gene_catalogs"
within_three_d_reduction_name <- "umap3d"
cross_feature_lookup_path <- "metadata/cross_feature_lookup.tsv"

cross_integration_keys <- c("camex", "saturn")

cross_integration_registry <- list(
    camex = list(
        label = "Camex",
        tab_title = "Cross-species integration (Camex)",
        eyebrow = "Cross-species integration",
        section_title = "Shared expression space across the three species",
        description = "This tab uses the Camex integration object. All queries resolve to shared Medicago-space features before plots are generated.",
        path = "app_ready_integration/camex/clustered_dataset.rds",
        slim_path = "app_ready_integration/camex/clustered_dataset_app_slim.rds",
        feature_mode = "medicago_space",
        default_group_by = "species_cell_class",
        default_split_by = "species"
    ),
    saturn = list(
        label = "SATURN",
        tab_title = "Cross-species integration (SATURN)",
        eyebrow = "Cross-species integration",
        section_title = "Shared expression space across the three species",
        description = "This tab uses the SATURN integration object. Selected source genes resolve to ortholog features from all three species in the integrated SATURN feature space.",
        path = "app_ready_integration/saturn/clustered_dataset.rds",
        slim_path = "app_ready_integration/saturn/clustered_dataset_app_slim.rds",
        feature_mode = "species_prefixed",
        default_group_by = "species",
        default_split_by = "species"
    )
)

pick_first_existing_path <- function(paths) {
    existing_path <- paths[file.exists(paths)][1]

    if (is.na(existing_path) || !length(existing_path)) {
        return(paths[[1]])
    }

    existing_path
}

get_cross_dataset_path <- function(cross_key) {
    cfg <- cross_integration_registry[[cross_key]]
    pick_first_existing_path(c(cfg$slim_path, cfg$path))
}

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

cross_integration_label <- function(cross_key) {
    cross_integration_registry[[cross_key]]$label
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

display_gene_labels <- function(species_key, gene_ids, include_gene_id_with_common = TRUE) {
    gene_ids <- as.character(gene_ids)

    annotations <- read_gene_annotations(species_key)

    tibble(
        gene_id = gene_ids,
        canonical_gene_id = canonicalize_gene_ids(species_key, gene_ids)
    ) %>%
        left_join(annotations, by = "canonical_gene_id") %>%
        mutate(
            has_common_name = !is.na(common_name) &
                nzchar(common_name) &
                common_name != gene_id,
            display_label = case_when(
                has_common_name & isTRUE(include_gene_id_with_common) ~ paste0(common_name, " (", gene_id, ")"),
                has_common_name ~ common_name,
                TRUE ~ gene_id
            )
        ) %>%
        pull(display_label)
}

set_heatmap_row_labels <- function(heatmap_obj, labels) {
    labels <- as.character(labels)

    heatmap_obj@row_names_param$labels <- labels

    row_name_anno <- heatmap_obj@row_names_param$anno
    if (inherits(row_name_anno, "AnnotationFunction") &&
        exists("value", envir = row_name_anno@var_env, inherits = FALSE)) {
        assign("value", labels, envir = row_name_anno@var_env)
    }

    heatmap_obj
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
        label = "Download",
        class = "btn btn-default btn-sm plot-download-btn"
    )
}

spinning_plot_output <- function(output_id, proxy_height = "360px", shell_class = NULL) {
    plot_tag <- plotOutput(output_id, height = "auto")

    if (is.null(shell_class) || !nzchar(shell_class)) {
        plot_tag
    } else {
        div(class = shell_class, plot_tag)
    }
}

spinning_plotly_output <- function(output_id, proxy_height = "360px", shell_class = NULL) {
    plot_tag <- plotly::plotlyOutput(output_id, height = proxy_height)

    if (is.null(shell_class) || !nzchar(shell_class)) {
        plot_tag
    } else {
        div(class = shell_class, plot_tag)
    }
}

umap_plot_shell_class <- function(split_by) {
    if (is.null(split_by) || !length(split_by) || identical(split_by, "none")) {
        "umap-plot-shell is-narrow"
    } else {
        "umap-plot-shell"
    }
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

format_within_feature_panel_title <- function(title, source_species, target_species) {
    title <- as.character(title %||% "")

    if (
        target_species %in% c("glycine", "lotus") &&
        !identical(source_species, target_species) &&
        grepl(" -> ", title, fixed = TRUE)
    ) {
        return(sub(" -> ", "\n", title, fixed = TRUE))
    }

    title
}

compact_feature_legend_guides <- function() {
    guides(
        colour = guide_colourbar(
            title = NULL,
            barwidth = grid::unit(1.45, "cm"),
            barheight = grid::unit(0.18, "cm"),
            ticks = FALSE
        ),
        fill = guide_colourbar(
            title = NULL,
            barwidth = grid::unit(1.45, "cm"),
            barheight = grid::unit(0.18, "cm"),
            ticks = FALSE
        )
    )
}

compact_feature_legend_theme <- function() {
    theme(
        legend.title = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.text = element_text(size = 8),
        legend.key.height = grid::unit(0.2, "cm"),
        legend.key.width = grid::unit(0.9, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.spacing.x = grid::unit(0.08, "cm"),
        legend.spacing.y = grid::unit(0.02, "cm")
    )
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

within_feature_grid_cols <- function(feature_n, split_by = "none") {
    feature_n <- max(1L, as.integer(feature_n %||% 1L))

    if (!identical(split_by, "none")) {
        return(1L)
    }

    if (feature_n <= 2L) {
        return(feature_n)
    }

    if (feature_n <= 4L) {
        return(2L)
    }

    3L
}

feature_umap_height_px <- function(feature_n, feature_cols, split_by = "none", panels_per_gene = 1L) {
    feature_n <- max(0L, as.integer(feature_n %||% 0L))
    feature_cols <- max(1L, as.integer(feature_cols %||% 1L))

    if (identical(split_by, "none")) {
        feature_rows <- ceiling(feature_n / feature_cols)
        return(max(560L, as.integer(ceiling(feature_rows * 400L))))
    }

    rows_per_gene <- ceiling(max(1L, as.integer(panels_per_gene %||% 1L)) / feature_cols)
    max(700L, as.integer(ceiling(feature_n * (rows_per_gene * 280 + 120))))
}

feature_umap_height_inches <- function(feature_n, feature_cols, split_by = "none", panels_per_gene = 1L) {
    feature_umap_height_px(
        feature_n = feature_n,
        feature_cols = feature_cols,
        split_by = split_by,
        panels_per_gene = panels_per_gene
    ) / 95
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

condition_level_order <- c(
    "Roots",
    "roots",
    "0.5h",
    "6h",
    "24h",
    "48h",
    "96h",
    "5dpi",
    "10dpi",
    "12dpi",
    "14dpi",
    "15dpi",
    "14d",
    "21pdi",
    "28dpi"
)

condition_root_palette <- c(
    "Roots" = "#1B4332",
    "roots" = "#1B4332"
)

condition_progression_palette <- c(
    "#C84C09",
    "#2C7FB8",
    "#8E44AD",
    "#16A085",
    "#D81B60",
    "#B8860B",
    "#457B9D",
    "#C06C84",
    "#E9C46A",
    "#8D99AE",
    "#4D9078",
    "#E07A5F",
    "#B7A58F"
)

distribution_sample_palette <- c(
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
    "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
    "#8C564B", "#E377C2", "#17BECF", "#BCBD22", "#7F7F7F"
)

metadata_level_order <- function(values, column_name) {
    values_chr <- as.character(values)
    present_levels <- unique(values_chr[!is.na(values_chr) & nzchar(values_chr)])

    if (!length(present_levels)) {
        return(character(0))
    }

    if (column_name %in% c("Group", "condition")) {
        return(c(
            condition_level_order[condition_level_order %in% present_levels],
            sort(setdiff(present_levels, condition_level_order))
        ))
    }

    if (column_name %in% c("Sample", "sample")) {
        return(sort(present_levels))
    }

    if (column_name %in% c("cluster_label", "Rank_1st", "Rank_2nd", "Rank_3rd", "Rank_4th", "Rank_5th")) {
        return(cluster_value_levels(present_levels))
    }

    if (is.factor(values)) {
        return(levels(values)[levels(values) %in% present_levels])
    }

    sort(present_levels)
}

order_metadata_values <- function(values, column_name) {
    values_chr <- as.character(values)
    ordered_levels <- metadata_level_order(values, column_name)

    factor(values_chr, levels = ordered_levels, ordered = TRUE)
}

condition_color_map <- function(values, column_name) {
    ordered_levels <- metadata_level_order(values, column_name)

    if (!length(ordered_levels)) {
        return(setNames(character(0), character(0)))
    }

    palette_values <- setNames(rep(NA_character_, length(ordered_levels)), ordered_levels)
    root_levels <- ordered_levels[ordered_levels %in% names(condition_root_palette)]

    if (length(root_levels)) {
        palette_values[root_levels] <- unname(condition_root_palette[root_levels])
    }

    non_root_levels <- ordered_levels[is.na(palette_values)]

    if (length(non_root_levels)) {
        progression_values <- condition_progression_palette[seq_len(min(length(non_root_levels), length(condition_progression_palette)))]

        if (length(non_root_levels) > length(condition_progression_palette)) {
            progression_values <- c(
                progression_values,
                grDevices::hcl.colors(
                    length(non_root_levels) - length(condition_progression_palette),
                    palette = "Dynamic"
                )
            )
        }

        palette_values[non_root_levels] <- progression_values
    }

    palette_values[ordered_levels]
}

metadata_colors_use <- function(values, column_name) {
    ordered_levels <- metadata_level_order(values, column_name)

    if (!length(ordered_levels)) {
        return(NULL)
    }

    if (column_name %in% c("Group", "condition")) {
        return(condition_color_map(values, column_name))
    }

    if (column_name %in% c("Sample", "sample")) {
        palette_values <- distribution_sample_palette[seq_len(min(length(ordered_levels), length(distribution_sample_palette)))]

        if (length(ordered_levels) > length(distribution_sample_palette)) {
            palette_values <- c(
                palette_values,
                grDevices::hcl.colors(
                    length(ordered_levels) - length(distribution_sample_palette),
                    palette = "Dynamic"
                )
            )
        }

        names(palette_values) <- ordered_levels
        return(palette_values)
    }

    NULL
}

composition_colors_use <- function(values, column_name) {
    ordered_values <- order_metadata_values(values, column_name)
    ordered_levels <- levels(ordered_values) %||% unique(as.character(ordered_values))
    ordered_levels <- ordered_levels[!is.na(ordered_levels) & nzchar(ordered_levels)]

    if (!length(ordered_levels)) {
        return(NULL)
    }

    if (column_name %in% c("Group", "condition")) {
        return(condition_color_map(values, column_name))
    }

    if (column_name %in% c("Sample", "sample")) {
        return(metadata_colors_use(values, column_name))
    }

    metadata_colors_use(values, column_name)
}

apply_metadata_display_order <- function(obj, columns) {
    valid_columns <- unique(columns[columns %in% colnames(obj@meta.data)])

    if (!length(valid_columns)) {
        return(obj)
    }

    ordered_obj <- obj

    for (column_name in valid_columns) {
        ordered_obj@meta.data[[column_name]] <- order_metadata_values(
            ordered_obj@meta.data[[column_name]],
            column_name
        )
    }

    ordered_obj
}

distribution_colors_use <- function(obj, group_by) {
    if (!(group_by %in% colnames(obj@meta.data))) {
        return(NULL)
    }

    unname(distribution_color_map(obj@meta.data[[group_by]], group_by))
}

distribution_color_map <- function(values, column_name) {
    ordered_levels <- metadata_level_order(values, column_name)

    if (!length(ordered_levels)) {
        return(setNames(character(0), character(0)))
    }

    palette_values <- metadata_colors_use(values, column_name)

    if (is.null(palette_values)) {
        palette_values <- scCustomize::scCustomize_Palette(
            num_groups = length(ordered_levels),
            color_seed = 123
        )
        names(palette_values) <- ordered_levels
    }

    palette_values[ordered_levels]
}

stratified_point_sample <- function(df, group_col, max_points = 30000L, seed = 123L) {
    if (!nrow(df) || nrow(df) <= max_points || !(group_col %in% colnames(df))) {
        return(df)
    }

    group_values <- as.character(df[[group_col]])
    group_values[is.na(group_values) | !nzchar(group_values)] <- "__missing__"
    split_indices <- split(seq_len(nrow(df)), group_values, drop = TRUE)
    group_sizes <- lengths(split_indices)

    if (!length(group_sizes)) {
        return(df)
    }

    raw_quota <- group_sizes * max_points / sum(group_sizes)
    quota <- pmax(1L, floor(raw_quota))
    quota <- pmin(quota, group_sizes)
    remainder <- max_points - sum(quota)

    if (remainder > 0L) {
        spare_capacity <- group_sizes - quota
        order_idx <- order(raw_quota - quota, decreasing = TRUE)

        for (idx in order_idx) {
            if (remainder <= 0L) {
                break
            }
            if (spare_capacity[[idx]] <= 0L) {
                next
            }
            quota[[idx]] <- quota[[idx]] + 1L
            spare_capacity[[idx]] <- spare_capacity[[idx]] - 1L
            remainder <- remainder - 1L
        }
    }

    set.seed(seed)
    sampled_indices <- purrr::map2(
        split_indices,
        quota,
        function(idx, n_keep) {
            if (length(idx) <= n_keep) {
                idx
            } else {
                sample(idx, size = n_keep)
            }
        }
    ) %>%
        unlist(use.names = FALSE)

    sampled_indices <- sort(sampled_indices)
    df[sampled_indices, , drop = FALSE]
}

metadata_column_label <- function(column_name) {
    dplyr::case_when(
        identical(column_name, "Group") ~ "Condition",
        identical(column_name, "condition") ~ "Condition",
        identical(column_name, "Sample") ~ "Sample",
        identical(column_name, "sample") ~ "Sample",
        identical(column_name, "batch") ~ "Batch",
        identical(column_name, "time_point") ~ "Time point",
        identical(column_name, "cluster_label") ~ "Cluster",
        identical(column_name, "saturn_label") ~ "SATURN label",
        identical(column_name, "saturn_ref_label") ~ "SATURN ref label",
        identical(column_name, "Rank_1st_label") ~ "Clustering opt 1 label",
        identical(column_name, "Rank_1st") ~ "Clustering opt 1",
        identical(column_name, "Rank_2nd") ~ "Clustering opt 2",
        identical(column_name, "Rank_3rd") ~ "Clustering opt 3",
        identical(column_name, "Rank_4th") ~ "Clustering opt 4",
        identical(column_name, "Rank_5th") ~ "Clustering opt 5",
        TRUE ~ column_name
    )
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

prepare_within_object <- function(obj) {
    obj$cluster_label <- as.character(Idents(obj))
    obj
}

prepare_cross_object <- function(obj, cross_key) {
    obj$cluster_label <- as.character(Idents(obj))

    label_column <- dplyr::case_when(
        "cell_class" %in% colnames(obj@meta.data) ~ "cell_class",
        "Rank_1st_label" %in% colnames(obj@meta.data) ~ "Rank_1st_label",
        "saturn_label" %in% colnames(obj@meta.data) ~ "saturn_label",
        TRUE ~ "cluster_label"
    )

    obj$species_cell_class <- paste(
        as.character(obj$species %||% "Unknown"),
        as.character(obj@meta.data[[label_column]] %||% "Unknown"),
        sep = " | "
    )
    obj$cross_integration <- cross_key
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
        verbose = FALSE,
        n_threads = max(1L, min(4L, parallel::detectCores(logical = TRUE) - 1L))
    )

    rownames(embedding) <- rownames(pca_embeddings)
    colnames(embedding) <- paste0("UMAP3D_", seq_len(ncol(embedding)))
    embedding
}

get_within_umap3d <- function(species_key, integration_method) {
    cache_key <- paste("within_umap3d", species_key, integration_method, sep = "::")

    cache_get(cache_key, function() {
        obj <- get_within_object(species_key, integration_method)
        available_reductions <- Reductions(obj)
        reduction_name <- c(within_three_d_reduction_name, "umap_3d", "umap3d")[c(within_three_d_reduction_name, "umap_3d", "umap3d") %in% available_reductions][1]

        if (!is.na(reduction_name) && length(reduction_name) && ncol(Embeddings(obj, reduction_name)) >= 3L) {
            return(Embeddings(obj, reduction_name)[, seq_len(3), drop = FALSE])
        }

        compute_umap3d_matrix(obj)
    })
}

get_cross_object <- function(cross_key) {
    cache_key <- paste("cross::object", cross_key, sep = "::")

    cache_get(cache_key, function() {
        readRDS(get_cross_dataset_path(cross_key)) %>%
            prepare_cross_object(cross_key)
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

get_cross_feature_lookup <- function(cross_key) {
    cache_key <- paste("cross::lookup", cross_key, sep = "::")

    cache_get(cache_key, function() {
        feature_mode <- cross_integration_registry[[cross_key]]$feature_mode

        if (identical(feature_mode, "medicago_space")) {
            feature_ids <- rownames(get_cross_object(cross_key))
            object_lookup <- tibble(
                feature_id = feature_ids,
                canonical_gene_id = canonicalize_gene_ids("medicago", feature_ids),
                feature_species = "medicago"
            ) %>%
                filter(!is.na(feature_id) & nzchar(feature_id)) %>%
                filter(!is.na(canonical_gene_id) & nzchar(canonical_gene_id)) %>%
                distinct(feature_species, canonical_gene_id, feature_id)

            cached_lookup <- read_tsv_cache(cross_feature_lookup_path)

            if (!is.null(cached_lookup) &&
                all(c("feature_id", "canonical_gene_id") %in% colnames(cached_lookup))) {
                cached_lookup <- cached_lookup %>%
                    select(feature_id, canonical_gene_id) %>%
                    mutate(feature_species = "medicago") %>%
                    filter(!is.na(feature_id) & nzchar(feature_id)) %>%
                    filter(!is.na(canonical_gene_id) & nzchar(canonical_gene_id)) %>%
                    distinct(feature_species, canonical_gene_id, feature_id)

                if (nrow(cached_lookup) == nrow(object_lookup) &&
                    setequal(cached_lookup$feature_id, object_lookup$feature_id)) {
                    return(cached_lookup)
                }
            }

            return(object_lookup)
        }

        feature_ids <- rownames(get_cross_object(cross_key))

        tibble(feature_id = feature_ids) %>%
            mutate(
                feature_species = sub("^([^:]+)::.*$", "\\1", feature_id),
                feature_gene_id = sub("^[^:]+::", "", feature_id)
            ) %>%
            filter(!is.na(feature_species) & !is.na(feature_gene_id)) %>%
            filter(feature_species %in% within_species_keys) %>%
            mutate(
                canonical_gene_id = purrr::map2_chr(
                    feature_species,
                    feature_gene_id,
                    function(feature_species, feature_gene_id) {
                        canonicalize_gene_ids(feature_species, feature_gene_id)
                    }
                )
            ) %>%
            filter(!is.na(canonical_gene_id) & nzchar(canonical_gene_id)) %>%
            distinct(feature_species, canonical_gene_id, feature_id, feature_gene_id)
    })
}

resolve_cross_integration_mapping <- function(source_species, source_genes, cross_key) {
    feature_mode <- cross_integration_registry[[cross_key]]$feature_mode

    if (identical(feature_mode, "medicago_space")) {
        return(resolve_target_mapping(
            source_species = source_species,
            source_genes = source_genes,
            target_species = "medicago",
            integration_method = NULL,
            cross_space = TRUE
        ))
    }

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

    source_orthogroups <- resolve_source_orthogroups(source_species, source_genes)
    genes_with_orthogroups <- source_orthogroups %>%
        filter(!is.na(orthogroup)) %>%
        distinct(source_gene) %>%
        pull(source_gene)
    no_orthogroup <- setdiff(source_genes, genes_with_orthogroups)

    candidate_tbl <- source_orthogroups %>%
        filter(!is.na(orthogroup)) %>%
        mutate(
            target_members = purrr::map(
                orthogroup,
                function(og) {
                    bind_rows(lapply(within_species_keys, function(target_species) {
                        tibble(
                            target_species = target_species,
                            target_gene_original = get_orthogroup_members(og, target_species)
                        )
                    }))
                }
            )
        ) %>%
        tidyr::unnest(target_members) %>%
        mutate(
            target_gene_original = as.character(target_gene_original),
            target_canonical = purrr::map2_chr(
                target_species,
                target_gene_original,
                function(target_species, target_gene_original) {
                    canonicalize_gene_ids(target_species, target_gene_original)
                }
            )
        )

    genes_with_target_members <- candidate_tbl %>%
        filter(!is.na(target_gene_original) & nzchar(target_gene_original)) %>%
        distinct(source_gene) %>%
        pull(source_gene)
    no_target_members <- setdiff(genes_with_orthogroups, genes_with_target_members)

    lookup <- get_cross_feature_lookup(cross_key)

    mapped_tbl <- candidate_tbl %>%
        left_join(
            lookup,
            by = c("target_species" = "feature_species", "target_canonical" = "canonical_gene_id")
        )

    genes_with_candidate_features <- mapped_tbl %>%
        filter(!is.na(target_gene_original) & nzchar(target_gene_original)) %>%
        distinct(source_gene) %>%
        pull(source_gene)
    genes_with_mapped_features <- mapped_tbl %>%
        filter(!is.na(feature_id)) %>%
        distinct(source_gene) %>%
        pull(source_gene)
    missing_features <- setdiff(genes_with_candidate_features, genes_with_mapped_features)

    plot_table <- mapped_tbl %>%
        filter(!is.na(feature_id)) %>%
        distinct(source_gene, orthogroup, target_species, target_gene_original, feature_id) %>%
        group_by(feature_id, target_species) %>%
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
                    paste(
                        unique(display_gene_labels(
                            source_species,
                            genes,
                            include_gene_id_with_common = FALSE
                        )),
                        collapse = "; "
                    )
                }
            ),
            target_species_label = vapply(target_species, species_label, character(1)),
            target_display = purrr::map2_chr(
                target_species,
                target_gene_label,
                function(target_species, target_gene_label) {
                    target_gene <- trimws(unlist(strsplit(target_gene_label, ";", fixed = TRUE)))[1]
                    display_gene_labels(
                        target_species,
                        target_gene,
                        include_gene_id_with_common = FALSE
                    )[[1]]
                }
            ),
            plot_label = paste0(source_gene_display, " -> ", target_display, " [", target_species_label, "]")
        )

    multiplicity <- mapped_tbl %>%
        filter(!is.na(feature_id)) %>%
        distinct(source_gene, feature_id) %>%
        count(source_gene, name = "mapped_gene_count") %>%
        filter(mapped_gene_count > 1)

    list(
        mapping = mapped_tbl,
        plot_table = plot_table,
        plot_features = plot_table$feature_id,
        label_map = setNames(plot_table$plot_label, plot_table$feature_id),
        no_orthogroup = sort(unique(no_orthogroup)),
        no_target_members = sort(unique(no_target_members)),
        missing_features = sort(unique(missing_features)),
        multiplicity = multiplicity
    )
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

build_cross_dataset_summary_row <- function(cross_key, obj) {
    md <- obj@meta.data
    sample_col <- pick_sample_column(md)
    sample_values <- if (!is.na(sample_col)) md[[sample_col]] else character(0)

    tibble(
        dataset_scope = "cross",
        species_key = cross_key,
        integration_method = cross_key,
        integration_label = cross_integration_label(cross_key),
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

    cross_rows <- bind_rows(lapply(cross_integration_keys, function(cross_key) {
        build_cross_dataset_summary_row(
            cross_key = cross_key,
            obj = get_cross_object(cross_key)
        )
    }))

    bind_rows(
        within_rows,
        cross_rows
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

            required_cross_keys <- cross_integration_keys
            cached_cross_keys <- summary_df %>%
                filter(dataset_scope == "cross") %>%
                pull(integration_method) %>%
                unique()

            if (nrow(summary_df) && "sample_n" %in% colnames(summary_df) &&
                all(required_cross_keys %in% cached_cross_keys)) {
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

get_cross_dataset_summary <- function(cross_key) {
    summary_row <- get_atlas_summary_table() %>%
        filter(
            dataset_scope == "cross",
            integration_method == !!cross_key
        ) %>%
        slice_head(n = 1)

    if (!nrow(summary_row)) {
        summary_row <- build_cross_dataset_summary_row(
            cross_key = cross_key,
            obj = get_cross_object(cross_key)
        )
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
        get_cross_feature_lookup("camex")
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
        get_cross_feature_lookup("camex")
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
                    paste(
                        unique(display_gene_labels(
                            source_species,
                            genes,
                            include_gene_id_with_common = FALSE
                        )),
                        collapse = "; "
                    )
                }
            )
        )

    target_display_lookup <- setNames(
        display_gene_labels(
            effective_target_species,
            plot_table$target_feature_id,
            include_gene_id_with_common = FALSE
        ),
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
    maybe_add("Condition", "Group")
    maybe_add("Sample", "Sample")

    choices
}

within_distribution_split_choices <- function(obj) {
    available_cols <- colnames(obj@meta.data)
    choices <- c("No split" = "none")

    if ("Group" %in% available_cols) {
        choices <- c(choices, "Condition" = "Group")
    } else if ("condition" %in% available_cols) {
        choices <- c(choices, "Condition" = "condition")
    }

    if ("Sample" %in% available_cols) {
        choices <- c(choices, "Sample" = "Sample")
    } else if ("sample" %in% available_cols) {
        choices <- c(choices, "Sample" = "sample")
    }

    choices
}

within_feature_split_choices <- function(obj) {
    available_cols <- colnames(obj@meta.data)
    choices <- c("No split" = "none")

    if ("Group" %in% available_cols) {
        choices <- c(choices, "Condition" = "Group")
    } else if ("condition" %in% available_cols) {
        choices <- c(choices, "Condition" = "condition")
    }

    if ("Sample" %in% available_cols) {
        choices <- c(choices, "Sample" = "Sample")
    } else if ("sample" %in% available_cols) {
        choices <- c(choices, "Sample" = "sample")
    }

    choices
}

within_composition_choices <- function(obj) {
    available_cols <- colnames(obj@meta.data)
    choices <- character(0)

    condition_col <- pick_first_existing_col(obj@meta.data, c("Group", "condition"))
    sample_col <- pick_first_existing_col(obj@meta.data, c("Sample", "sample"))

    if (!is.na(condition_col) && condition_col %in% available_cols) {
        choices <- c(choices, "Condition" = condition_col)
    }

    if (!is.na(sample_col) && sample_col %in% available_cols) {
        choices <- c(choices, "Sample" = sample_col)
    }

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
    maybe_add("SATURN label", "saturn_label")
    maybe_add("Clustering opt 1 label", "Rank_1st_label")
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
                            value = 1.1,
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
                        uiOutput(paste0(species_key, "_distribution_umap_plot_ui"))
                    )
                )
            ),
            div(
                class = "subsection-header",
                h3("Cluster composition"),
                p("See what fraction of cells in each cluster comes from each condition or sample. Cluster definitions follow the active clustering choice above when available.")
            ),
            fluidRow(
                column(
                    width = 4,
                    div(
                        class = "option-group",
                        uiOutput(paste0(species_key, "_composition_by_ui"))
                    )
                )
            ),
            fluidRow(
                column(
                    width = 12,
                    div(
                        class = "plot-card",
                        div(class = "plot-card-title", "Cells per cluster"),
                        uiOutput(paste0(species_key, "_composition_plot_ui"))
                    )
                )
            ),
            div(
                class = "subsection-header",
                h3("Gene expression"),
                p("Plot the selected source genes directly or through ortholog mappings in this species atlas.")
            ),
            uiOutput(paste0(species_key, "_notice_ui")),
            fluidRow(
                column(
                    width = 12,
                    div(
                        class = "option-group",
                        uiOutput(paste0(species_key, "_split_by_ui"))
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
                            div(class = "plot-card-title", "Feature UMAPs"),
                            plot_download_button(paste0("dl_", species_key, "_umap"))
                        ),
                        uiOutput(paste0(species_key, "_umap_plot_ui"))
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
                            div(class = "plot-card-title", "Expression violin plots"),
                            plot_download_button(paste0("dl_", species_key, "_violin"))
                        ),
                        spinning_plot_output(paste0(species_key, "_violin_plot"), proxy_height = "520px")
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
                            div(class = "plot-card-title", "Multi-gene dot plot"),
                            plot_download_button(paste0("dl_", species_key, "_dot"))
                        ),
                        spinning_plot_output(paste0(species_key, "_dot_plot"), proxy_height = "420px")
                    )
                )
            )
        )
    )
}

cross_tab_ui <- function(cross_key) {
    integration_cfg <- cross_integration_registry[[cross_key]]
    prefix <- paste0("cross_", cross_key)

    tabPanel(
        title = integration_cfg$tab_title,
        value = prefix,
        div(
            class = "section-card",
            div(
                class = "section-header",
                div(class = "section-eyebrow", integration_cfg$eyebrow),
                h2(integration_cfg$section_title),
                p(integration_cfg$description)
            ),
            fluidRow(
                column(
                    width = 3,
                    div(
                        class = "option-group",
                        sliderInput(
                            inputId = paste0(prefix, "_pt_size"),
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
                            inputId = paste0(prefix, "_umap_columns"),
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
                        uiOutput(paste0(prefix, "_split_by_ui"))
                    )
                ),
                column(
                    width = 3,
                    div(
                        class = "option-group",
                        uiOutput(paste0(prefix, "_group_by_ui"))
                    )
                )
            ),
            uiOutput(paste0(prefix, "_notice_ui")),
            fluidRow(
                column(
                    width = 12,
                    div(
                        class = "plot-card",
                        div(
                            class = "plot-card-header",
                            div(class = "plot-card-title", "Integrated feature UMAPs"),
                            plot_download_button(paste0("dl_", prefix, "_umap"))
                        ),
                        uiOutput(paste0(prefix, "_umap_plot_ui"))
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
                            plot_download_button(paste0("dl_", prefix, "_dot"))
                        ),
                        spinning_plot_output(paste0(prefix, "_dot_plot"), proxy_height = "420px")
                    )
                ),
                column(
                    width = 5,
                    div(
                        class = "table-card",
                        div(class = "plot-card-title", paste(cross_integration_label(cross_key), "mapped features")),
                        uiOutput(paste0(prefix, "_mapping_table_ui"))
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
        do.call(
            tabsetPanel,
            c(
                list(
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
                    species_tab_ui("lotus")
                ),
                lapply(cross_integration_keys, cross_tab_ui)
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

    overview_cross_resolution <- reactive({
        resolve_cross_integration_mapping(
            source_species = input$source_species %||% "medicago",
            source_genes = selected_source_genes(),
            cross_key = "camex"
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
        within_cards <- lapply(within_species_keys, function(species_key) {
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

        cross_cards <- lapply(cross_integration_keys, function(cross_key) {
            summary <- get_cross_dataset_summary(cross_key)

            div(
                class = "dataset-tile",
                div(class = "dataset-title", paste("Cross-species", cross_integration_label(cross_key))),
                div(class = "dataset-value", paste(format_stat_value(summary$cells), "cells")),
                div(class = "dataset-note", paste(format_stat_value(summary$genes), "expressed genes")),
                div(class = "dataset-note", paste(format_stat_value(summary$sample_n), "samples"))
            )
        })

        div(
            class = "dataset-grid atlas-summary-grid",
            tagList(within_cards, cross_cards)
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

        cross_res <- overview_cross_resolution()

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
            clustering_columns <- c("Rank_1st", "Rank_2nd", "Rank_3rd", "Rank_4th", "Rank_5th", "cluster_label")

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

            expression_prompt <- sprintf(
                "Select a gene first to populate the expression panels for %s.",
                tab_label
            )

            expression_resolution <- reactive({
                validate(
                    need(
                        length(selected_source_genes()) > 0,
                        expression_prompt
                    )
                )

                tab_resolution()
            })

            tab_group_choices <- reactive(within_group_choices(tab_object()))
            tab_distribution_split_choices <- reactive(within_distribution_split_choices(tab_object()))
            tab_feature_split_choices <- reactive(within_feature_split_choices(tab_object()))
            tab_composition_choices <- reactive(within_composition_choices(tab_object()))

            output[[paste0(prefix, "_distribution_group_by_ui")]] <- renderUI({
                choices <- tab_group_choices()

                selectInput(
                    inputId = paste0(prefix, "_distribution_group_by"),
                    label = "Color cells by",
                    choices = choices,
                    selected = resolve_choice(
                        input[[paste0(prefix, "_distribution_group_by")]],
                        choices,
                        default = "Rank_1st"
                    )
                )
            })

            output[[paste0(prefix, "_distribution_split_by_ui")]] <- renderUI({
                choices <- tab_distribution_split_choices()

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

            output[[paste0(prefix, "_composition_by_ui")]] <- renderUI({
                choices <- tab_composition_choices()

                if (!length(choices)) {
                    return(NULL)
                }

                selectInput(
                    inputId = paste0(prefix, "_composition_by"),
                    label = "Show percentages by",
                    choices = choices,
                    selected = resolve_choice(
                        input[[paste0(prefix, "_composition_by")]],
                        choices,
                        default = unname(choices[[1]])
                    )
                )
            })

            output[[paste0(prefix, "_split_by_ui")]] <- renderUI({
                choices <- tab_feature_split_choices()

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

                resolution <- tab_resolution()

                if (!identical(source_species, species_key)) {
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

                if (!length(cards) &&
                    !length(resolution$no_orthogroup) &&
                    !length(resolution$no_target_members) &&
                    !length(resolution$missing_features) &&
                    !nrow(resolution$multiplicity)) {
                    return(NULL)
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
                    default = "Rank_1st"
                )
            })

            tab_distribution_split_by <- reactive({
                resolve_choice(
                    input[[paste0(prefix, "_distribution_split_by")]],
                    tab_distribution_split_choices(),
                    default = "none"
                )
            })

            tab_composition_by <- reactive({
                choices <- tab_composition_choices()
                default_choice <- if (length(choices)) unname(choices[[1]]) else NULL

                resolve_choice(
                    input[[paste0(prefix, "_composition_by")]],
                    choices,
                    default = default_choice
                )
            })

            tab_composition_cluster_by <- reactive({
                available_clusters <- clustering_columns[clustering_columns %in% colnames(tab_object()@meta.data)]
                preferred_cluster <- tab_distribution_group_by()

                if (length(preferred_cluster) && preferred_cluster %in% available_clusters) {
                    return(preferred_cluster)
                }

                if (length(available_clusters)) {
                    return(available_clusters[[1]])
                }

                NA_character_
            })

            tab_group_by <- reactive({
                resolve_choice(
                    tab_distribution_group_by(),
                    tab_group_choices(),
                    default = "Rank_1st"
                )
            })

            tab_split_by <- reactive({
                resolve_choice(
                    input[[paste0(prefix, "_split_by")]],
                    tab_feature_split_choices(),
                    default = "none"
                )
            })

            output[[paste0(prefix, "_distribution_umap_plot_ui")]] <- renderUI({
                split_by <- tab_distribution_split_by()

                if (identical(split_by, "none")) {
                    return(
                        div(
                            class = "distribution-split-layout",
                            div(
                                class = "distribution-split-panel",
                                div(class = "distribution-view-title", "2D UMAP"),
                                spinning_plot_output(
                                    paste0(prefix, "_distribution_umap_plot"),
                                    proxy_height = "540px",
                                    shell_class = "umap-plot-shell"
                                )
                            ),
                            div(
                                class = "distribution-split-panel",
                                div(class = "distribution-view-title", "3D UMAP"),
                                spinning_plotly_output(
                                    paste0(prefix, "_distribution_umap3d_plot"),
                                    proxy_height = "540px",
                                    shell_class = "plotly-plot-shell"
                                )
                            )
                        )
                    )
                }

                spinning_plot_output(
                    paste0(prefix, "_distribution_umap_plot"),
                    proxy_height = "520px",
                    shell_class = umap_plot_shell_class(split_by)
                )
            })

            output[[paste0(prefix, "_composition_plot_ui")]] <- renderUI({
                spinning_plot_output(
                    paste0(prefix, "_composition_plot"),
                    proxy_height = "480px"
                )
            })

            output[[paste0(prefix, "_umap_plot_ui")]] <- renderUI({
                spinning_plot_output(
                    paste0(prefix, "_umap_plot"),
                    proxy_height = "520px",
                    shell_class = "umap-plot-shell"
                )
            })

            distribution_umap_plot_obj <- reactive({
                obj <- tab_object()
                group_by <- tab_distribution_group_by()
                split_by <- tab_distribution_split_by()
                pt_size <- as.numeric(input[[paste0(prefix, "_distribution_pt_size")]] %||% 1.1)
                split_enabled <- !identical(split_by, "none")
                plot_obj <- apply_metadata_display_order(obj, c(group_by, split_by))
                color_map <- distribution_color_map(plot_obj@meta.data[[group_by]], group_by)
                colors_use <- unname(color_map)
                show_cluster_labels <- group_by %in% clustering_columns && identical(split_by, "none")
                split_panels <- if (split_enabled) split_panel_count(plot_obj, split_by) else 1L
                split_columns <- if (identical(split_by, "none")) {
                    NULL
                } else {
                    min(4L, split_panels)
                }

                distribution_plot <- scCustomize::DimPlot_scCustom(
                    seurat_object = plot_obj,
                    colors_use = colors_use,
                    group.by = group_by,
                    split.by = if (split_enabled) split_by else NULL,
                    pt.size = pt_size,
                    label = show_cluster_labels,
                    repel = TRUE,
                    raster = TRUE,
                    num_columns = split_columns
                )

                distribution_plot <- distribution_plot &
                    app_plot_theme() &
                    theme(
                        legend.title = element_blank(),
                        legend.position = if (split_enabled || group_by %in% clustering_columns) "none" else "top",
                        panel.grid = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        axis.line = element_blank(),
                        plot.margin = margin(8, 14, 10, 10)
                    )

                if (split_enabled) {
                    distribution_plot &
                        labs(color = NULL) &
                        theme(
                            plot.title = element_text(
                                face = "bold",
                                colour = app_palette["text"],
                                size = 13,
                                hjust = 0.5
                            )
                        )
                } else {
                    distribution_plot &
                        labs(title = NULL, color = NULL)
                }
            })

            distribution_umap3d_plot_data <- reactive({
                validate(
                    need(
                        identical(tab_distribution_split_by(), "none"),
                        "3D UMAP is available only when the distribution view is not split."
                    )
                )

                group_by <- tab_distribution_group_by()
                obj <- apply_metadata_display_order(tab_object(), group_by)
                embedding <- get_within_umap3d(
                    species_key = species_key,
                    integration_method = current_species_integration(species_key)
                )
                cell_ids <- intersect(rownames(obj@meta.data), rownames(embedding))

                validate(
                    need(length(cell_ids) > 0, "No 3D UMAP coordinates are available for this atlas.")
                )

                distribution_df <- tibble(
                    cell_id = cell_ids,
                    group_value = as.character(obj@meta.data[cell_ids, group_by, drop = TRUE]),
                    umap3d_1 = embedding[cell_ids, 1],
                    umap3d_2 = embedding[cell_ids, 2],
                    umap3d_3 = embedding[cell_ids, 3]
                ) %>%
                    filter(
                        !is.na(group_value) & nzchar(group_value),
                        !is.na(umap3d_1),
                        !is.na(umap3d_2),
                        !is.na(umap3d_3)
                    )

                validate(
                    need(nrow(distribution_df) > 0, "No cells are available for the 3D UMAP view.")
                )

                distribution_df <- distribution_df %>%
                    mutate(group_value = order_metadata_values(group_value, group_by))

                distribution_df <- stratified_point_sample(
                    distribution_df,
                    group_col = "group_value",
                    max_points = 30000L
                )

                full_color_map <- distribution_color_map(obj@meta.data[[group_by]], group_by)
                present_levels <- names(full_color_map)[names(full_color_map) %in% as.character(distribution_df$group_value)]

                list(
                    data = distribution_df,
                    color_map = full_color_map[present_levels],
                    group_by = group_by
                )
            })

            composition_plot_obj <- reactive({
                obj <- tab_object()
                composition_by <- tab_composition_by()
                cluster_by <- tab_composition_cluster_by()
                plot_obj <- apply_metadata_display_order(obj, composition_by)

                validate(
                    need(!is.null(composition_by) && nzchar(composition_by), "No composition metadata are available for this atlas."),
                    need(!is.na(cluster_by) && nzchar(cluster_by), "No clustering metadata are available for this atlas.")
                )

                md <- plot_obj@meta.data
                composition_df <- tibble(
                    cluster = as.character(md[[cluster_by]]),
                    composition = as.character(md[[composition_by]])
                ) %>%
                    filter(
                        !is.na(cluster) & nzchar(cluster),
                        !is.na(composition) & nzchar(composition)
                    ) %>%
                    count(cluster, composition, name = "cell_count")

                validate(
                    need(nrow(composition_df) > 0, "No cluster composition data are available for the selected metadata.")
                )

                composition_df <- composition_df %>%
                    mutate(
                        cluster = factor(cluster, levels = cluster_value_levels(cluster)),
                        composition = order_metadata_values(composition, composition_by)
                    )

                fill_values <- composition_colors_use(composition_df$composition, composition_by)
                legend_rows <- max(1L, min(3L, ceiling(dplyr::n_distinct(composition_df$composition) / 8L)))

                ggplot(
                    composition_df,
                    aes(x = cluster, y = cell_count, fill = composition)
                ) +
                    geom_col(
                        position = "fill",
                        width = 0.82,
                        colour = "#FFFFFF",
                        linewidth = 0.18
                    ) +
                    scale_y_continuous(
                        labels = scales::label_percent(accuracy = 1),
                        expand = expansion(mult = c(0, 0.02))
                    ) +
                    {if (!is.null(fill_values)) scale_fill_manual(values = fill_values) else scale_fill_discrete()} +
                    guides(fill = guide_legend(nrow = legend_rows, byrow = TRUE)) +
                    labs(
                        x = metadata_column_label(cluster_by),
                        y = "Percent of cells",
                        fill = NULL
                    ) +
                    app_plot_theme() +
                    theme(
                        panel.grid.major.x = element_blank(),
                        axis.text.x = element_text(size = 11),
                        legend.position = "top",
                        legend.text = element_text(size = 9),
                        plot.margin = margin(8, 12, 8, 10)
                    )
            })

            rank_data <- reactive({
                resolution <- expression_resolution()
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
	                resolution <- expression_resolution()
	                obj <- tab_object()
	                source_species <- input$source_species %||% "medicago"
	                reference_group_by <- tab_composition_cluster_by()
	                split_by <- tab_split_by()
	                pt_size <- max(2.0, as.numeric(input[[paste0(prefix, "_distribution_pt_size")]] %||% 1.1) * 1.8)
	                split_enabled <- !identical(split_by, "none")
	                feature_panel_n <- length(resolution$plot_features) + 1L
	                plot_obj <- apply_metadata_display_order(obj, c(split_by, reference_group_by))
	                split_panels <- if (split_enabled) split_panel_count(plot_obj, split_by) else 1L
	                split_columns <- if (split_enabled) min(4L, split_panels) else NULL
	                feature_grid_cols <- within_feature_grid_cols(
	                    feature_n = feature_panel_n,
	                    split_by = split_by
	                )

	                validate(
	                    need(
	                        length(resolution$plot_features) > 0,
	                        paste("No mapped genes are available for", tab_label, "in the selected atlas.")
	                    ),
	                    need(
	                        !is.na(reference_group_by) && nzchar(reference_group_by),
	                        paste("No clustering metadata are available for", tab_label, ".")
	                    )
	                )

	                reference_color_map <- distribution_color_map(
	                    plot_obj@meta.data[[reference_group_by]],
	                    reference_group_by
	                )

	                reference_plot <- scCustomize::DimPlot_scCustom(
	                    seurat_object = plot_obj,
	                    colors_use = unname(reference_color_map),
	                    group.by = reference_group_by,
	                    split.by = if (split_enabled) split_by else NULL,
	                    pt.size = pt_size,
	                    label = identical(split_by, "none"),
	                    repel = TRUE,
	                    raster = TRUE,
	                    num_columns = split_columns
	                )

	                reference_plot <- reference_plot &
	                    labs(color = NULL) &
	                    app_plot_theme() &
	                    theme(
	                        legend.title = element_blank(),
	                        legend.position = "none",
	                        panel.grid = element_blank(),
	                        axis.title = element_blank(),
	                        axis.text = element_blank(),
	                        axis.ticks = element_blank(),
	                        axis.line = element_blank(),
	                        plot.margin = margin(4, 8, 6, 6)
	                    )

	                reference_plot <- if (split_enabled) {
	                    reference_plot &
	                        theme(
	                            plot.title = element_text(
	                                face = "bold",
	                                colour = app_palette["text"],
	                                size = 13,
	                                hjust = 0.5
	                            )
	                        )
	                } else {
	                    reference_plot &
	                        labs(title = NULL) &
	                        theme(plot.title = element_blank())
	                }

	                plot_list <- c(list(
	                    wrap_titled_plot(
	                        plot_obj = reference_plot,
	                        title = metadata_column_label(reference_group_by)
	                    )
	                ), lapply(resolution$plot_features, function(feature_id) {
	                    feature_plot <- scCustomize::FeaturePlot_scCustom(
	                        seurat_object = plot_obj,
	                        features = feature_id,
                        split.by = if (split_enabled) split_by else NULL,
                        pt.size = pt_size,
                        order = TRUE,
                        raster = TRUE,
                        num_columns = split_columns,
                        label = FALSE,
                        combine = TRUE
                    )

                    feature_plot <- feature_plot &
                        labs(color = NULL) &
                        app_plot_theme() &
                        compact_feature_legend_guides() &
                        compact_feature_legend_theme() &
                        scale_x_continuous(expand = expansion(mult = 0.01)) &
                        scale_y_continuous(expand = expansion(mult = 0.01)) &
                        theme(
                            panel.grid = element_blank(),
                            axis.title = element_blank(),
                            axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            axis.line = element_blank(),
                            plot.margin = margin(4, 8, 6, 6)
                        )

                    feature_plot <- if (split_enabled) {
                        feature_plot &
                            theme(
                                plot.title = element_text(
                                    face = "bold",
                                    colour = app_palette["text"],
                                    size = 13,
                                    hjust = 0.5
                                )
                            )
                    } else {
                        feature_plot &
                            labs(title = NULL) &
                            theme(plot.title = element_blank())
                    }

	                    wrap_titled_plot(
	                        plot_obj = feature_plot,
	                        title = format_within_feature_panel_title(
	                            title = unname(resolution$label_map[feature_id]),
	                            source_species = source_species,
	                            target_species = species_key
	                        )
	                    )
	                }))

	                wrap_plots(plotlist = plot_list, ncol = feature_grid_cols)
	            })

            violin_plot_obj <- reactive({
                resolution <- expression_resolution()
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
                        pt.size = 0.12,
                        num_columns = 1,
                        raster = FALSE
                    )

                    if (length(violin_plot$layers) >= 2) {
                        violin_plot$layers[[2]]$aes_params$colour <- unname(app_palette["green_dark"])
                        violin_plot$layers[[2]]$aes_params$alpha <- 0.42
                        violin_plot$layers[[2]]$aes_params$size <- 0.34
                    }

                    violin_plot <- violin_plot &
                        labs(
                            title = unname(label_map[feature_id]),
                            x = NULL,
                            y = "Normalized expression"
                        ) &
                        app_plot_theme() &
                        theme(
                            plot.title = element_text(
                                face = "bold",
                                colour = app_palette["text"],
                                size = 16,
                                hjust = 0
                            ),
                            legend.position = "none",
                            axis.text.x = element_text(angle = 35, hjust = 1),
                            panel.grid.major.x = element_blank(),
                            plot.margin = margin(8, 14, 12, 10)
                        )

                    violin_plot
                })

                if (length(plot_list) == 1) {
                    plot_list[[1]]
                } else {
                    wrap_plots(plotlist = plot_list, ncol = 1)
                }
            })

	            dot_plot_obj <- reactive({
	                resolution <- expression_resolution()
	                obj <- tab_object()

                validate(
                    need(
                        length(resolution$plot_features) > 0,
                        paste("No mapped genes are available for", tab_label, "in the selected atlas.")
                    )
                )

	                label_map <- unname(resolution$label_map[resolution$plot_features])
	                cluster_feature_enabled <- length(resolution$plot_features) > 1

	                dot_plot_precheck <- tryCatch(
	                    scCustomize::DotPlot_scCustom(
	                        seurat_object = obj,
	                        features = resolution$plot_features,
	                        group.by = tab_group_by()
	                    )$data,
	                    error = function(e) NULL
	                )

	                if (!is.null(dot_plot_precheck) && "avg.exp.scaled" %in% colnames(dot_plot_precheck)) {
	                    has_nonfinite_feature <- any(
	                        !is.finite(dot_plot_precheck$avg.exp.scaled) &
	                            !is.na(dot_plot_precheck$features.plot)
	                    )

	                    if (isTRUE(has_nonfinite_feature)) {
	                        cluster_feature_enabled <- FALSE
	                    }
	                }

	                build_clustered_dot_plot <- function(cluster_feature) {
	                    scCustomize::Clustered_DotPlot(
	                        seurat_object = obj,
	                        features = resolution$plot_features,
	                        group.by = tab_group_by(),
	                        cluster_feature = cluster_feature,
	                        cluster_ident = TRUE,
	                        flip = FALSE,
	                        row_names_side = "left",
	                        column_names_side = "bottom",
	                        show_ident_colors = FALSE,
	                        show_ident_legend = FALSE,
	                        grid_color = NA,
	                        row_label_size = 11,
	                        column_label_size = 11,
	                        legend_position = "right",
	                        legend_label_size = 10,
	                        legend_title_size = 10,
	                        x_lab_rotate = TRUE,
	                        raster = FALSE,
	                        plot_km_elbow = FALSE
	                    )
	                }

	                plot_obj <- tryCatch(
	                    build_clustered_dot_plot(cluster_feature_enabled),
	                    error = function(e) {
	                        if (isTRUE(cluster_feature_enabled)) {
	                            build_clustered_dot_plot(FALSE)
	                        } else {
	                            stop(e)
	                        }
	                    }
	                )

	                clustered_dot_radius_mm <- 3.2
	                legend_dot_radius_mm <- 2
	                heatmap_obj <- plot_obj@ht_list[[1]]
	                percent_scale_max <- NA_real_

	                if (!is.null(heatmap_obj@matrix_param$cell_fun)) {
	                    cell_fun_env <- environment(heatmap_obj@matrix_param$cell_fun)
	                    percent_mat <- get("percent_mat", envir = cell_fun_env)
	                    percent_scale_max <- max(percent_mat, na.rm = TRUE)
	                    heatmap_obj@matrix_param$cell_fun <- eval(
	                        bquote(
	                            function(j, i, x, y, w, h, fill) {
	                                grid::grid.rect(
	                                    x = x,
	                                    y = y,
	                                    width = w,
	                                    height = h,
	                                    gp = grid::gpar(col = grid_color, fill = NA)
	                                )
	                                grid::grid.circle(
	                                    x = x,
	                                    y = y,
	                                    r = sqrt(percent_mat[i, j] / .(max(percent_scale_max, 1e-08))) * grid::unit(.(clustered_dot_radius_mm), "mm"),
	                                    gp = grid::gpar(fill = col_fun(exp_mat[i, j]), col = NA)
	                                )
	                            }
	                        ),
	                        envir = cell_fun_env
	                    )
	                }

	                if (!is.null(heatmap_obj@matrix_param$layer_fun)) {
	                    layer_fun_env <- environment(heatmap_obj@matrix_param$layer_fun)
	                    if (!is.finite(percent_scale_max) || percent_scale_max <= 0) {
	                        percent_scale_max <- max(get("percent_mat", envir = layer_fun_env), na.rm = TRUE)
	                    }
	                    heatmap_obj@matrix_param$layer_fun <- eval(
	                        bquote(
	                            function(j, i, x, y, w, h, fill) {
	                                grid::grid.rect(
	                                    x = x,
	                                    y = y,
	                                    width = w,
	                                    height = h,
	                                    gp = grid::gpar(col = grid_color, fill = NA)
	                                )
	                                grid::grid.circle(
	                                    x = x,
	                                    y = y,
	                                    r = sqrt(ComplexHeatmap::pindex(percent_mat, i, j) / .(max(percent_scale_max, 1e-08))) * grid::unit(.(clustered_dot_radius_mm), "mm"),
	                                    gp = grid::gpar(
	                                        fill = col_fun(ComplexHeatmap::pindex(exp_mat, i, j)),
	                                        col = NA
	                                    )
	                                )
	                            }
	                        ),
	                        envir = layer_fun_env
	                    )
	                }

	                if (!is.finite(percent_scale_max) || percent_scale_max <= 0) {
	                    percent_scale_max <- 1
	                }

	                legend_breaks <- pretty(c(0, percent_scale_max), n = 4)
	                legend_breaks <- legend_breaks[legend_breaks > 0 & legend_breaks < percent_scale_max]
	                legend_breaks <- unique(c(legend_breaks, percent_scale_max))
	                legend_breaks <- legend_breaks[legend_breaks > 0]

	                if (!length(legend_breaks)) {
	                    legend_breaks <- percent_scale_max
	                }

	                legend_digits <- if (percent_scale_max < 5) {
	                    2L
	                } else if (percent_scale_max < 20) {
	                    1L
	                } else {
	                    0L
	                }

	                legend_labels <- sub(
	                    "\\.?0+$",
	                    "",
	                    formatC(legend_breaks, format = "f", digits = legend_digits)
	                )

	                size_legend <- ComplexHeatmap::Legend(
	                    title = "Percent Expressing",
	                    at = legend_labels,
	                    graphics = lapply(legend_breaks, function(value) {
	                        force(value)
	                        function(x, y, w, h) {
	                            grid::grid.circle(
	                                x = x,
	                                y = y,
	                                r = sqrt(value / percent_scale_max) * grid::unit(legend_dot_radius_mm, "mm"),
	                                gp = grid::gpar(fill = "black", col = NA)
	                            )
	                        }
	                    }),
	                    labels_gp = grid::gpar(fontsize = 10),
	                    title_gp = grid::gpar(fontsize = 10, fontface = "bold")
	                )

	                row_feature_ids <- rownames(heatmap_obj@matrix)
	                if (is.null(row_feature_ids) || !length(row_feature_ids)) {
	                    row_feature_ids <- resolution$plot_features
	                }

	                dot_plot_row_labels <- unname(resolution$label_map[row_feature_ids])
	                missing_label_idx <- which(is.na(dot_plot_row_labels) | !nzchar(dot_plot_row_labels))
	                if (length(missing_label_idx)) {
	                    dot_plot_row_labels[missing_label_idx] <- row_feature_ids[missing_label_idx]
	                }

	                plot_obj@ht_list[[1]] <- set_heatmap_row_labels(
	                    heatmap_obj,
	                    dot_plot_row_labels
	                )

	                if (length(plot_obj@heatmap_legend_param$list)) {
	                    plot_obj@heatmap_legend_param$list[[1]] <- size_legend
	                } else {
	                    plot_obj@heatmap_legend_param$list <- list(size_legend)
	                }

	                plot_obj
	            })

            draw_clustered_dot_plot <- function(plot_obj) {
                ComplexHeatmap::draw(
                    plot_obj,
                    heatmap_legend_side = "right",
                    annotation_legend_side = "right",
                    merge_legends = TRUE
                )
            }

            clustered_dot_plot_height <- function() {
                feature_n <- tryCatch(length(expression_resolution()$plot_features), error = function(e) 0L)
                group_n <- tryCatch(dplyr::n_distinct(tab_object()@meta.data[[tab_group_by()]]), error = function(e) 0L)

                max(
                    460,
                    170 + feature_n * 42 + group_n * 11
                )
            }

	            output[[paste0(prefix, "_umap_plot")]] <- renderPlot(
	                {
	                    umap_plot_obj()
	                },
	                height = function() {
	                    feature_n <- tryCatch(length(expression_resolution()$plot_features) + 1L, error = function(e) 0L)
	                    split_by <- tryCatch(tab_split_by(), error = function(e) "none")
	                    feature_cols <- if (identical(split_by, "none")) {
	                        within_feature_grid_cols(feature_n = feature_n, split_by = split_by)
                    } else {
                        4L
                    }

                    if (identical(split_by, "none")) {
                        feature_umap_height_px(
                            feature_n = feature_n,
                            feature_cols = feature_cols,
                            split_by = split_by
                        )
                    } else {
                        panels_per_gene <- tryCatch(
                            split_panel_count(tab_object(), split_by),
                            error = function(e) 1L
                        )
                        feature_umap_height_px(
                            feature_n = feature_n,
                            feature_cols = feature_cols,
                            split_by = split_by,
                            panels_per_gene = panels_per_gene
                        )
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
                        return(540)
                    }

                    panel_n <- tryCatch(
                        split_panel_count(tab_object(), split_by),
                        error = function(e) 1L
                    )
                    split_cols <- min(4L, max(1L, panel_n))
                    split_rows <- ceiling(panel_n / split_cols)

                    max(520, split_rows * 320)
                },
                res = 110
            )

            output[[paste0(prefix, "_distribution_umap3d_plot")]] <- plotly::renderPlotly({
                plot_data <- distribution_umap3d_plot_data()
                df <- plot_data$data
                color_map <- plot_data$color_map
                group_by <- plot_data$group_by
                marker_size <- max(1.8, as.numeric(input[[paste0(prefix, "_distribution_pt_size")]] %||% 1.1) * 2.2)
                group_levels <- names(color_map)

                p <- plotly::plot_ly(
                    type = "scatter3d",
                    mode = "markers"
                )

                for (group_level in group_levels) {
                    group_df <- df %>%
                        filter(as.character(group_value) == !!group_level)

                    if (!nrow(group_df)) {
                        next
                    }

                    p <- p %>%
                        plotly::add_trace(
                            data = group_df,
                            x = ~umap3d_1,
                            y = ~umap3d_2,
                            z = ~umap3d_3,
                            name = group_level,
                            marker = list(
                                size = marker_size,
                                color = unname(color_map[group_level]),
                                opacity = 0.8
                            ),
                            text = ~paste0(
                                metadata_column_label(group_by), ": ", group_value,
                                "<br>Cell: ", cell_id
                            ),
                            hovertemplate = "%{text}<extra></extra>"
                        )
                }

                p %>%
                    plotly::layout(
                        margin = list(l = 0, r = 0, b = 0, t = 10),
                        legend = list(
                            orientation = "h",
                            x = 0,
                            y = 1.08,
                            font = list(
                                size = 14,
                                color = unname(app_palette["text"])
                            ),
                            itemsizing = "constant",
                            itemwidth = 40,
                            bgcolor = "rgba(255,255,255,0.82)",
                            bordercolor = "rgba(201,214,196,0.92)",
                            borderwidth = 1
                        ),
                        scene = list(
                            aspectmode = "data",
                            xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                            yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                            zaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                            camera = list(
                                eye = list(x = 1.22, y = 1.06, z = 0.75)
                            )
                        )
                    ) %>%
                    plotly::config(displaylogo = FALSE)
            })

            output[[paste0(prefix, "_composition_plot")]] <- renderPlot(
                {
                    composition_plot_obj()
                },
                height = function() {
                    composition_by <- tryCatch(tab_composition_by(), error = function(e) NULL)
                    obj <- tryCatch(tab_object(), error = function(e) NULL)

                    if (is.null(obj) || is.null(composition_by) || !nzchar(composition_by)) {
                        return(460)
                    }

                    level_count <- tryCatch({
                        values <- obj@meta.data[[composition_by]]
                        values <- as.character(values)
                        length(unique(values[!is.na(values) & nzchar(values)]))
                    }, error = function(e) 1L)

                    legend_rows <- max(1L, min(3L, ceiling(level_count / 8L)))
                    max(460, 380 + legend_rows * 44)
                },
                res = 110
            )

            output[[paste0(prefix, "_violin_plot")]] <- renderPlot(
                {
                    violin_plot_obj()
                },
                height = function() {
                    feature_n <- tryCatch(length(expression_resolution()$plot_features), error = function(e) 0L)
                    max(320, 250 * feature_n)
                },
                res = 110
            )

            output[[paste0(prefix, "_dot_plot")]] <- renderPlot(
                {
                    draw_clustered_dot_plot(dot_plot_obj())
                },
                height = function() {
                    clustered_dot_plot_height()
                },
                res = 110
            )

	            output[[paste0("dl_", prefix, "_umap")]] <- downloadHandler(
	                filename = function() {
	                    paste0(prefix, "_umap.", get_ext())
	                },
	                content = function(file) {
	                    resolution <- expression_resolution()
	                    feature_n <- length(resolution$plot_features) + 1L
	                    split_by <- tab_split_by()
	                    feature_cols <- if (identical(split_by, "none")) {
	                        within_feature_grid_cols(feature_n = feature_n, split_by = split_by)
                    } else {
                        4L
                    }

                    plot_height <- feature_umap_height_inches(
                        feature_n = feature_n,
                        feature_cols = feature_cols,
                        split_by = split_by,
                        panels_per_gene = if (identical(split_by, "none")) {
                            1L
                        } else {
                            split_panel_count(tab_object(), split_by)
                        }
                    )

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
                        split_cols <- min(4L, max(1L, panel_n))
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
                    feature_n <- length(expression_resolution()$plot_features)

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
                    plot_obj <- dot_plot_obj()
                    ext <- get_ext()
                    plot_height <- clustered_dot_plot_height() / 95

                    if (identical(ext, "png")) {
                        png(filename = file, width = 1900, height = max(1200, clustered_dot_plot_height() * 2), res = 220)
                    } else if (identical(ext, "svg")) {
                        svglite::svglite(file = file, width = 14, height = plot_height)
                    } else {
                        pdf(file = file, width = 14, height = plot_height)
                    }

                    draw_clustered_dot_plot(plot_obj)
                    grDevices::dev.off()
                }
            )
        })
    }

    walk(within_species_keys, register_species_tab)

    register_cross_tab <- function(cross_key) {
        local({
            integration_cfg <- cross_integration_registry[[cross_key]]
            prefix <- paste0("cross_", cross_key)

            cross_object <- reactive({
                get_cross_object(cross_key)
            })

            cross_resolution <- reactive({
                resolve_cross_integration_mapping(
                    source_species = input$source_species %||% "medicago",
                    source_genes = selected_source_genes(),
                    cross_key = cross_key
                )
            })

            output[[paste0(prefix, "_group_by_ui")]] <- renderUI({
                choices <- cross_group_choices(cross_object())

                selectInput(
                    inputId = paste0(prefix, "_group_by"),
                    label = "Summarize expression by",
                    choices = choices,
                    selected = resolve_choice(
                        input[[paste0(prefix, "_group_by")]],
                        choices,
                        default = integration_cfg$default_group_by
                    )
                )
            })

            output[[paste0(prefix, "_split_by_ui")]] <- renderUI({
                choices <- cross_split_choices(cross_object())

                selectInput(
                    inputId = paste0(prefix, "_split_by"),
                    label = "Split UMAP by",
                    choices = choices,
                    selected = resolve_choice(
                        input[[paste0(prefix, "_split_by")]],
                        choices,
                        default = integration_cfg$default_split_by
                    )
                )
            })

            cross_group_by <- reactive({
                resolve_choice(
                    input[[paste0(prefix, "_group_by")]],
                    cross_group_choices(cross_object()),
                    default = integration_cfg$default_group_by
                )
            })

            cross_split_by <- reactive({
                resolve_choice(
                    input[[paste0(prefix, "_split_by")]],
                    cross_split_choices(cross_object()),
                    default = integration_cfg$default_split_by
                )
            })

            output[[paste0(prefix, "_umap_plot_ui")]] <- renderUI({
                spinning_plot_output(
                    paste0(prefix, "_umap_plot"),
                    proxy_height = "520px",
                    shell_class = umap_plot_shell_class(cross_split_by())
                )
            })

            output[[paste0(prefix, "_notice_ui")]] <- renderUI({
                genes <- selected_source_genes()
                resolution <- cross_resolution()
                feature_mode <- integration_cfg$feature_mode

                intro_card <- if (identical(feature_mode, "medicago_space")) {
                    notice_card(
                        title = "Shared Medicago-space features",
                        body = paste(
                            cross_integration_label(cross_key),
                            "stores a shared feature space represented with Medicago identifiers. Soybean and Lotus selections are projected to Medicago orthologs in this tab."
                        ),
                        tone = "info"
                    )
                } else {
                    notice_card(
                        title = paste(cross_integration_label(cross_key), "ortholog feature space"),
                        body = paste(
                            cross_integration_label(cross_key),
                            "stores species-prefixed features. Source genes are resolved to ortholog features from Medicago, Glycine, and Lotus before the integrated plots are drawn."
                        ),
                        tone = "info"
                    )
                }

                cards <- list(intro_card)

                if (!length(genes)) {
                    cards <- append(cards, list(
                        notice_card(
                            title = "No source genes selected",
                            body = paste(
                                "Add one or more source-species genes to populate the",
                                cross_integration_label(cross_key),
                                "plots."
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
                            title = if (identical(feature_mode, "medicago_space")) {
                                "Orthogroups without Medicago members"
                            } else {
                                "Orthogroups without mapped members"
                            },
                            body = if (identical(feature_mode, "medicago_space")) {
                                compact_gene_list(resolution$no_target_members, limit = 8)
                            } else {
                                paste(
                                    "No ortholog members from Medicago, Glycine, or Lotus were found in the mapped orthogroups for:",
                                    compact_gene_list(resolution$no_target_members, limit = 8)
                                )
                            },
                            tone = "warning"
                        )
                    ))
                }

                if (length(resolution$missing_features)) {
                    cards <- append(cards, list(
                        notice_card(
                            title = paste("Mapped orthologs missing from the", cross_integration_label(cross_key), "feature set"),
                            body = compact_gene_list(resolution$missing_features, limit = 8),
                            tone = "warning"
                        )
                    ))
                }

                if (nrow(resolution$multiplicity)) {
                    multiplicity_text <- resolution$multiplicity %>%
                        mutate(label = paste0(source_gene, " (", mapped_gene_count, " mapped features)")) %>%
                        pull(label)

                    cards <- append(cards, list(
                        notice_card(
                            title = paste("One-to-many mappings in", cross_integration_label(cross_key)),
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

                if (identical(integration_cfg$feature_mode, "medicago_space")) {
                    resolution$plot_table %>%
                        transmute(
                            `Source gene(s)` = source_gene_display,
                            `Medicago-space feature` = display_gene_labels(
                                "medicago",
                                target_feature_id,
                                include_gene_id_with_common = FALSE
                            ),
                            Orthogroup = ifelse(nzchar(orthogroup_label), orthogroup_label, "NA")
                        )
                } else {
                    resolution$plot_table %>%
                        transmute(
                            `Source gene(s)` = source_gene_display,
                            Species = target_species_label,
                            `SATURN feature` = target_display,
                            Orthogroup = ifelse(nzchar(orthogroup_label), orthogroup_label, "NA")
                        )
                }
            })

            output[[paste0(prefix, "_mapping_table_ui")]] <- renderUI({
                genes <- selected_source_genes()

                if (!length(genes)) {
                    return(div(
                        class = "summary-placeholder",
                        paste(
                            "Add source-species genes to see which features are used in the",
                            cross_integration_label(cross_key),
                            "tab."
                        )
                    ))
                }

                mapping_tbl <- cross_mapping_table()

                if (!nrow(mapping_tbl)) {
                    return(div(
                        class = "summary-placeholder",
                        paste(
                            "No selected genes resolve to the",
                            cross_integration_label(cross_key),
                            "feature set."
                        )
                    ))
                }

                html_summary_table(mapping_tbl)
            })

            cross_umap_plot_obj <- reactive({
                resolution <- cross_resolution()
                obj <- cross_object()
                split_by <- cross_split_by()
                requested_cols <- as.integer(input[[paste0(prefix, "_umap_columns")]] %||% 1)
                pt_size <- max(0.65, as.numeric(input[[paste0(prefix, "_pt_size")]] %||% 0.45))
                feature_grid_cols <- resolve_feature_grid_cols(
                    requested_cols = requested_cols,
                    feature_labels = unname(resolution$label_map[resolution$plot_features]),
                    split_by = split_by
                )

                validate(
                    need(
                        length(resolution$plot_features) > 0,
                        paste(
                            "No selected genes resolve to features in the",
                            cross_integration_label(cross_key),
                            "integration."
                        )
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
                        compact_feature_legend_guides() &
                        compact_feature_legend_theme() &
                        scale_x_continuous(expand = expansion(mult = 0.01)) &
                        scale_y_continuous(expand = expansion(mult = 0.01)) &
                        theme(
                            plot.title = element_blank(),
                            panel.grid = element_blank(),
                            axis.title = element_blank(),
                            axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            axis.line = element_blank(),
                            plot.margin = margin(4, 8, 6, 6)
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
                        paste(
                            "No selected genes resolve to features in the",
                            cross_integration_label(cross_key),
                            "integration."
                        )
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

            output[[paste0(prefix, "_umap_plot")]] <- renderPlot(
                {
                    cross_umap_plot_obj()
                },
                height = function() {
                    feature_n <- tryCatch(length(cross_resolution()$plot_features), error = function(e) 0L)
                    split_by <- tryCatch(cross_split_by(), error = function(e) "none")
                    feature_cols <- tryCatch(
                        resolve_feature_grid_cols(
                            requested_cols = as.integer(input[[paste0(prefix, "_umap_columns")]] %||% 1),
                            feature_labels = unname(cross_resolution()$label_map[cross_resolution()$plot_features]),
                            split_by = split_by
                        ),
                        error = function(e) 1L
                    )

                    if (identical(split_by, "none")) {
                        feature_umap_height_px(
                            feature_n = feature_n,
                            feature_cols = feature_cols,
                            split_by = split_by
                        )
                    } else {
                        panels_per_gene <- tryCatch(
                            split_panel_count(cross_object(), split_by),
                            error = function(e) 1L
                        )
                        feature_umap_height_px(
                            feature_n = feature_n,
                            feature_cols = feature_cols,
                            split_by = split_by,
                            panels_per_gene = panels_per_gene
                        )
                    }
                },
                res = 110
            )

            output[[paste0(prefix, "_dot_plot")]] <- renderPlot(
                {
                    cross_dot_plot_obj()
                },
                height = function() {
                    feature_n <- tryCatch(length(cross_resolution()$plot_features), error = function(e) 0L)
                    max(420, 95 * feature_n + 120)
                },
                res = 110
            )

            output[[paste0("dl_", prefix, "_umap")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_umap.", get_ext())
                },
                content = function(file) {
                    feature_n <- length(cross_resolution()$plot_features)
                    split_by <- cross_split_by()
                    feature_cols <- resolve_feature_grid_cols(
                        requested_cols = as.integer(input[[paste0(prefix, "_umap_columns")]] %||% 1),
                        feature_labels = unname(cross_resolution()$label_map[cross_resolution()$plot_features]),
                        split_by = split_by
                    )

                    plot_height <- feature_umap_height_inches(
                        feature_n = feature_n,
                        feature_cols = feature_cols,
                        split_by = split_by,
                        panels_per_gene = if (identical(split_by, "none")) {
                            1L
                        } else {
                            split_panel_count(cross_object(), split_by)
                        }
                    )

                    save_ggplot(
                        file = file,
                        plot_obj = cross_umap_plot_obj(),
                        width = 15,
                        height = plot_height
                    )
                }
            )

            output[[paste0("dl_", prefix, "_dot")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_dotplot.", get_ext())
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
        })
    }

    walk(cross_integration_keys, register_cross_tab)
}

shinyApp(ui = ui, server = server)
