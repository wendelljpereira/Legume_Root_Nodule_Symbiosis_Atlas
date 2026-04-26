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

pick_first_existing_path <- function(paths) {
    existing_path <- paths[file.exists(paths)][1]

    if (is.na(existing_path) || !length(existing_path)) {
        return(paths[[1]])
    }

    existing_path
}

app_slim_path <- function(path) {
    sub("\\.rds$", "_app_slim.rds", path, ignore.case = TRUE)
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
cluster_markers_cache_dir <- "metadata/cluster_markers"
ui_choice_cache_path <- "metadata/ui_choice_cache.tsv"
ui_cluster_lookup_cache_path <- "metadata/ui_cluster_lookup_cache.tsv"
within_three_d_reduction_name <- "umap3d"
cross_feature_lookup_path <- "metadata/cross_feature_lookup.tsv"

# User-overridable paths. See README for Docker usage; set these env vars to
# point at a locally-maintained orthogroups table or celltype override folder
# without editing the source.
resolve_env_path <- function(env_var, default) {
    value <- Sys.getenv(env_var, unset = "")
    if (!nzchar(value)) default else value
}

source(file.path("R", "atlas_annotations.R"), local = TRUE)

orthogroups_path <- resolve_env_path("ATLAS_ORTHOLOGS_PATH", "orthogroups/joint_orthogroups.tsv")
cluster_annotations_dir <- resolve_env_path("ATLAS_CLUSTER_ANNOTATIONS_DIR", "annotations/cluster_annotations")
celltype_overrides_dir <- resolve_env_path("ATLAS_CELLTYPE_OVERRIDES_DIR", "celltype_overrides")

atlas_version <- resolve_env_path("ATLAS_VERSION", "1.0")
atlas_last_updated <- resolve_env_path("ATLAS_LAST_UPDATED", format(Sys.Date(), "%Y-%m-%d"))
atlas_citation_text <- resolve_env_path(
    "ATLAS_CITATION",
    "Pereira W. et al. A cross-species single-cell atlas of legume root nodule symbiosis. Legume Root Nodule Symbiosis Atlas, version 1.0."
)
atlas_cache_max_entries <- max(
    8L,
    suppressWarnings(as.integer(Sys.getenv("ATLAS_CACHE_MAX_ENTRIES", "64")))
)
atlas_expression_cache_max_entries <- max(
    8L,
    suppressWarnings(as.integer(Sys.getenv("ATLAS_EXPRESSION_CACHE_MAX_ENTRIES", "96")))
)
atlas_plot_feature_limit <- max(
    1L,
    suppressWarnings(as.integer(Sys.getenv("ATLAS_PLOT_FEATURE_LIMIT", "80")))
)

# Optional access gate. Set ATLAS_ACCESS_PASSWORD on the server
# (ShinyApps.io > Advanced > Environment variables) when a deployment
# should require a password before the app renders.
atlas_access_password <- Sys.getenv("ATLAS_ACCESS_PASSWORD", unset = "")
atlas_access_required <- nzchar(atlas_access_password)

cross_integration_keys <- c("camex", "saturn")
visible_cross_integration_keys <- setdiff(cross_integration_keys, "camex")

cross_integration_registry <- list(
    camex = list(
        label = "Camex",
        tab_title = "Cross-species integration (Camex)",
        eyebrow = "Cross-species integration",
        section_title = "Shared expression space across the three species",
        description = "This tab uses the Camex integration object. All queries resolve to shared Medicago-space features before plots are generated.",
        path = "app_ready_integration/camex/clustered_dataset.rds",
        slim_path = "app_ready_integration/camex/clustered_dataset_app_slim.rds",
        umap3d_path = "app_ready_integration/camex/umap3d_embeddings.rds",
        feature_mode = "medicago_space",
        default_group_by = "species",
        default_split_by = "species"
    ),
    saturn = list(
        label = "SATURN",
        tab_title = "Cross-species integration (SATURN)",
        eyebrow = "Cross-species integration",
        section_title = "Shared expression space across the three species",
        description = "Selected source genes resolve to ortholog features from all three species in the integrated SATURN feature space.",
        path = "app_ready_integration/saturn/clustered_dataset.rds",
        slim_path = "app_ready_integration/saturn/clustered_dataset_app_slim.rds",
        umap3d_path = "app_ready_integration/saturn/umap3d_embeddings.rds",
        feature_mode = "species_prefixed",
        default_group_by = "species",
        default_split_by = "species"
    )
)

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
    accent_dark = "#3f5343",
    accent = "#5f7a64",
    accent_soft = "#eef1ec",
    accent_subtle = "#e8ece7",
    text = "#1a1f26",
    muted = "#556170",
    subtle = "#7b8591",
    border = "#e4e7ec",
    border_strong = "#cfd4dc",
    warn = "#a16207",
    warn_soft = "#fefce8",
    # legacy aliases kept for compatibility with existing references
    green_dark = "#3f5343",
    green = "#5f7a64",
    green_soft = "#eef1ec",
    warm = "#a16207",
    warm_soft = "#fefce8",
    slate = "#7b8591",
    red_soft = "#fef2f2"
)

species_palette <- c(
    "Glycine max" = "#4E79A7",
    "Lotus japonicus" = "#59A14F",
    "Medicago truncatula" = "#E15759"
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

format_gene_ids_for_display <- function(species_key, ids) {
    ids <- as.character(ids)
    display_ids <- ids

    if (identical(species_key, "glycine")) {
        display_ids <- sub("\\.Wm82\\.a[0-9]+\\.v1$", "", display_ids, perl = TRUE)
    }

    display_ids[is.na(ids)] <- NA_character_
    display_ids
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

compact_display_gene_list <- function(species_key, values, limit = 5, empty_label = "None", include_gene_id_with_common = TRUE) {
    display_values <- display_gene_labels(
        species_key,
        values,
        include_gene_id_with_common = include_gene_id_with_common
    )
    compact_gene_list(display_values, limit = limit, empty_label = empty_label)
}

gene_search_tokens <- function(species_key, feature_id, display_label, common_name = NULL, synonyms = NULL, existing_tokens = NULL) {
    values <- c(
        feature_id,
        format_gene_ids_for_display(species_key, feature_id),
        display_label,
        common_name,
        synonyms,
        existing_tokens
    )
    values <- unique(trimws(as.character(values)))
    values <- values[!is.na(values) & nzchar(values)]
    paste(values, collapse = " ")
}

display_gene_labels <- function(species_key, gene_ids, include_gene_id_with_common = TRUE) {
    gene_ids <- as.character(gene_ids)

    annotations <- read_gene_annotations(species_key)

    tibble(
        gene_id = gene_ids,
        display_gene_id = format_gene_ids_for_display(species_key, gene_ids),
        canonical_gene_id = canonicalize_gene_ids(species_key, gene_ids)
    ) %>%
        left_join(annotations, by = "canonical_gene_id") %>%
        mutate(
            has_common_name = !is.na(common_name) &
                nzchar(common_name) &
                common_name != display_gene_id,
            display_label = case_when(
                has_common_name & isTRUE(include_gene_id_with_common) ~ paste0(common_name, " (", display_gene_id, ")"),
                has_common_name ~ common_name,
                TRUE ~ display_gene_id
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

figure_preset_config <- function(preset = c("exploratory", "presentation", "publication")) {
    preset <- match.arg(preset)

    switch(
        preset,
        exploratory = list(
            base_size = 13,
            title_size = 13.5,
            axis_title_size = 12.5,
            axis_text_size = 11,
            legend_title_size = 11.5,
            legend_text_size = 10.5,
            strip_text_size = 11.5,
            legend_position = "top",
            width_scale = 1,
            height_scale = 1
        ),
        presentation = list(
            base_size = 16,
            title_size = 17,
            axis_title_size = 15,
            axis_text_size = 13,
            legend_title_size = 13.5,
            legend_text_size = 12.5,
            strip_text_size = 13.5,
            legend_position = "bottom",
            width_scale = 1.08,
            height_scale = 1.08
        ),
        publication = list(
            base_size = 12,
            title_size = 12.5,
            axis_title_size = 11.5,
            axis_text_size = 10,
            legend_title_size = 10.5,
            legend_text_size = 9.5,
            strip_text_size = 10.5,
            legend_position = "top",
            width_scale = 0.96,
            height_scale = 0.96
        )
    )
}

app_plot_theme <- function(base_size = NULL, preset = c("exploratory", "presentation", "publication")) {
    cfg <- figure_preset_config(preset)

    if (is.null(base_size) || !is.finite(base_size)) {
        base_size <- cfg$base_size
    }

    theme_minimal(base_size = base_size) +
        theme(
            plot.title = element_text(face = "bold", colour = app_palette["text"], size = cfg$title_size),
            plot.subtitle = element_text(colour = app_palette["muted"], size = cfg$axis_text_size),
            axis.title = element_text(face = "bold", colour = app_palette["text"], size = cfg$axis_title_size),
            axis.text = element_text(colour = app_palette["muted"], size = cfg$axis_text_size),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            legend.position = cfg$legend_position,
            legend.title = element_text(face = "bold", size = cfg$legend_title_size),
            legend.text = element_text(size = cfg$legend_text_size),
            strip.text = element_text(face = "bold", colour = app_palette["text"], size = cfg$strip_text_size),
            strip.background = element_rect(fill = app_palette["green_soft"], colour = NA),
            plot.margin = margin(10, 16, 10, 10)
        )
}

plot_title_annotation <- function(title, size = 16, hjust = 0) {
    plot_annotation(
        title = title,
        theme = theme(
            plot.title = element_text(
                face = "bold",
                colour = app_palette["text"],
                size = size,
                hjust = hjust
            )
        )
    )
}

wrap_titled_plot <- function(plot_obj, title, size = 16, hjust = 0) {
    wrap_elements(full = plot_obj + plot_title_annotation(title, size = size, hjust = hjust))
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

sc_customize_feature_palette <- function(colorblind_safe = FALSE) {
    if (isTRUE(colorblind_safe)) {
        return(viridisLite::cividis(100, end = 0.95))
    }

    palette_ref <- formals(scCustomize::FeaturePlot_scCustom)$colors_use
    if (is.symbol(palette_ref)) {
        return(get(as.character(palette_ref), envir = asNamespace("scCustomize")))
    }

    as.character(palette_ref)
}

sc_customize_feature_na_color <- function() {
    as.character(formals(scCustomize::FeaturePlot_scCustom)$na_color %||% "lightgray")
}

feature_umap_point_size <- function(pt_size) {
    pt_size <- suppressWarnings(as.numeric(pt_size)[[1]])
    if (!is.finite(pt_size)) {
        pt_size <- 1
    }

    # Match the feature UMAP point size to the cell-distribution UMAP control so
    # expression panels retain the same visual density as the reference map.
    max(0.015, pt_size)
}

# Build a feature plot where expressing cells are drawn on top while
# keeping the same point size as the background cells. This preserves
# the species-restricted SATURN layout while matching scCustomize's
# default feature color scale.
emphasized_feature_plot <- function(
    obj,
    feature_id,
    reduction = "umap",
    split_by = NULL,
    expression_values = NULL,
    cell_ids = NULL,
    fixed_max = NULL,
    pt_size = 0.45,
    colors_use = NULL,
    zero_color = NULL,
    expressing_size_boost = 1.8,
    background_alpha = 0.22,
    colorblind_safe = FALSE
) {
    available_reductions <- Reductions(obj)
    if (!length(available_reductions)) {
        stop("Seurat object has no reductions to plot.")
    }
    if (!(reduction %in% available_reductions)) {
        reduction <- available_reductions[[1]]
    }

    embeddings <- Embeddings(obj, reduction)[, 1:2, drop = FALSE]
    coord_df <- tibble::as_tibble(embeddings, rownames = "cell_id")
    colnames(coord_df)[2:3] <- c(".dim1", ".dim2")

    if (!is.null(expression_values)) {
        if (is.null(names(expression_values))) {
            expr_values <- as.numeric(expression_values)
            if (length(expr_values) != nrow(coord_df)) {
                expr_values <- rep(NA_real_, nrow(coord_df))
            }
        } else {
            expr_values <- as.numeric(expression_values[coord_df$cell_id])
        }
    } else {
        expr <- tryCatch(
            FetchData(obj, vars = feature_id),
            error = function(e) NULL
        )

        if (is.null(expr) || !nrow(expr)) {
            expr_values <- rep(NA_real_, nrow(coord_df))
        } else {
            expr_values <- as.numeric(expr[[1]])[match(coord_df$cell_id, rownames(expr))]
        }
    }

    coord_df$.expression <- expr_values
    coord_df$.expression[is.na(coord_df$.expression)] <- 0

    if (!is.null(cell_ids)) {
        cell_ids <- unique(as.character(cell_ids))
        coord_df <- coord_df %>%
            filter(.data$cell_id %in% cell_ids)
    }

    use_split <- !is.null(split_by) && !identical(split_by, "none") && nzchar(split_by)
    if (use_split) {
        if (split_by %in% colnames(obj@meta.data)) {
            meta_vec <- obj[[split_by, drop = TRUE]]
            coord_df$.split <- meta_vec[match(coord_df$cell_id, colnames(obj))]
        } else {
            use_split <- FALSE
        }
    }

    max_expr <- fixed_max %||% suppressWarnings(max(coord_df$.expression, na.rm = TRUE))
    if (!is.finite(max_expr) || max_expr <= 0) max_expr <- 1
    colors_use <- colors_use %||% sc_customize_feature_palette(colorblind_safe)
    zero_color <- zero_color %||% sc_customize_feature_na_color()

    coord_df <- coord_df %>%
        dplyr::mutate(
            .is_expressing = .data$.expression > 0,
            .alpha = ifelse(.data$.is_expressing, 0.95, background_alpha)
        ) %>%
        dplyr::arrange(.data$.expression)

    p <- ggplot(coord_df, aes(x = .data$.dim1, y = .data$.dim2)) +
        geom_point(
            data = dplyr::filter(coord_df, !.data$.is_expressing),
            aes(alpha = .data$.alpha),
            color = zero_color,
            size = pt_size,
            stroke = 0,
            shape = 16
        ) +
        geom_point(
            data = dplyr::filter(coord_df, .data$.is_expressing),
            aes(color = .data$.expression, alpha = .data$.alpha),
            size = pt_size,
            stroke = 0,
            shape = 16
        ) +
        scale_alpha_identity() +
        coord_fixed()

    p <- p + scale_color_gradientn(
        colors = colors_use,
        limits = c(1e-09, max_expr),
        na.value = zero_color
    )

    if (use_split) {
        p <- p + facet_wrap(~ .split)
    }

    p
}

empty_umap_message_plot <- function(message) {
    message <- paste(strwrap(as.character(message), width = 68), collapse = "\n")

    ggplot() +
        annotate(
            "text",
            x = 0.5,
            y = 0.5,
            label = message,
            colour = app_palette["muted"],
            size = 4.2
        ) +
        coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
        theme_void() +
        theme(
            plot.background = element_rect(fill = "white", colour = NA),
            plot.margin = margin(10, 10, 10, 10)
        )
}

aggregate_feature_expression <- function(obj, feature_ids) {
    get_cached_feature_expression(obj, feature_ids)
}

cross_species_cell_ids <- function(obj, species_key) {
    if (!("species" %in% colnames(obj@meta.data))) {
        return(character(0))
    }

    species_values <- as.character(obj@meta.data$species)
    rownames(obj@meta.data)[species_values == species_label(species_key)]
}

panel_feature_entries <- function(feature_ids, feature_labels = NULL) {
    feature_ids <- as.character(feature_ids)
    feature_labels <- as.character(feature_labels %||% feature_ids)

    tibble(
        feature_id = feature_ids,
        feature_label = feature_labels
    ) %>%
        filter(!is.na(feature_id) & nzchar(feature_id)) %>%
        mutate(
            feature_label = ifelse(is.na(feature_label) | !nzchar(feature_label), feature_id, feature_label)
        ) %>%
        distinct(feature_id, .keep_all = TRUE)
}

panel_feature_count <- function(entries) {
    if (is.null(entries) || !nrow(entries)) 0L else nrow(entries)
}

panel_feature_entry_at <- function(entries, idx, repeat_singleton = TRUE) {
    if (is.null(entries) || !nrow(entries)) {
        return(panel_feature_entries(character(0)))
    }

    if (nrow(entries) == 1L && isTRUE(repeat_singleton)) {
        return(entries[1, , drop = FALSE])
    }

    if (idx > nrow(entries)) {
        return(entries[0, , drop = FALSE])
    }

    entries[idx, , drop = FALSE]
}

cross_comparison_panel_title <- function(source_gene_display, slot_entries) {
    if (is.null(slot_entries) || !length(slot_entries)) {
        return(source_gene_display)
    }

    slot_labels <- unique(unlist(lapply(slot_entries, function(entries) {
        if (is.null(entries) || !nrow(entries)) {
            return(character(0))
        }

        entries$feature_label[[1]]
    })))
    slot_labels <- slot_labels[!is.na(slot_labels) & nzchar(slot_labels)]

    if (length(slot_labels) == 1L) {
        return(paste0(source_gene_display, " -> ", slot_labels[[1]]))
    }

    slot_species_labels <- unique(unlist(lapply(names(slot_entries), function(species_key) {
        entries <- slot_entries[[species_key]]

        if (is.null(entries) || !nrow(entries)) {
            return(character(0))
        }

        feature_label <- entries$feature_label[[1]]
        if (is.na(feature_label) || !nzchar(feature_label)) {
            return(character(0))
        }

        paste0(species_label(species_key), ": ", feature_label)
    })))
    slot_species_labels <- slot_species_labels[!is.na(slot_species_labels) & nzchar(slot_species_labels)]

    if (!length(slot_species_labels)) {
        return(source_gene_display)
    }

    paste0(source_gene_display, " -> ", paste(slot_species_labels, collapse = " | "))
}

cross_comparison_panel_specs <- function(resolution, cross_key) {
    plot_tbl <- resolution$plot_table

    if (is.null(plot_tbl) || !nrow(plot_tbl)) {
        return(tibble())
    }

    feature_mode <- cross_integration_registry[[cross_key]]$feature_mode

    if (identical(feature_mode, "medicago_space")) {
        grouped_tbl <- plot_tbl %>%
            mutate(
                panel_order = dplyr::row_number(),
                orthogroup_label = ifelse(is.na(orthogroup_label), "", orthogroup_label),
                target_feature_id = as.character(target_feature_id)
            ) %>%
            group_by(source_gene_display, orthogroup_label) %>%
            summarise(
                panel_order = min(panel_order),
                medicago_panels = list(panel_feature_entries(
                    feature_ids = unique(target_feature_id),
                    feature_labels = display_gene_labels("medicago", unique(target_feature_id), include_gene_id_with_common = FALSE)
                )),
                glycine_panels = list(panel_feature_entries(
                    feature_ids = unique(target_feature_id),
                    feature_labels = display_gene_labels("medicago", unique(target_feature_id), include_gene_id_with_common = FALSE)
                )),
                lotus_panels = list(panel_feature_entries(
                    feature_ids = unique(target_feature_id),
                    feature_labels = display_gene_labels("medicago", unique(target_feature_id), include_gene_id_with_common = FALSE)
                )),
                .groups = "drop"
            ) %>%
            arrange(panel_order)
    } else {
        grouped_tbl <- plot_tbl %>%
            mutate(
                panel_order = dplyr::row_number(),
                orthogroup_label = ifelse(is.na(orthogroup_label), "", orthogroup_label),
                feature_id = as.character(feature_id)
            ) %>%
            group_by(source_gene_display, orthogroup_label) %>%
            summarise(
                panel_order = min(panel_order),
                medicago_panels = list(panel_feature_entries(
                    feature_ids = feature_id[target_species == "medicago"],
                    feature_labels = target_display[target_species == "medicago"]
                )),
                glycine_panels = list(panel_feature_entries(
                    feature_ids = feature_id[target_species == "glycine"],
                    feature_labels = target_display[target_species == "glycine"]
                )),
                lotus_panels = list(panel_feature_entries(
                    feature_ids = feature_id[target_species == "lotus"],
                    feature_labels = target_display[target_species == "lotus"]
                )),
                .groups = "drop"
            ) %>%
            arrange(panel_order)
    }

    grouped_tbl <- bind_rows(lapply(seq_len(nrow(grouped_tbl)), function(row_idx) {
        grouped_row <- grouped_tbl[row_idx, , drop = FALSE]
        slot_count <- max(
            panel_feature_count(grouped_row$medicago_panels[[1]]),
            panel_feature_count(grouped_row$glycine_panels[[1]]),
            panel_feature_count(grouped_row$lotus_panels[[1]])
        )

        if (!is.finite(slot_count) || slot_count < 1L) {
            slot_count <- 1L
        }

        bind_rows(lapply(seq_len(slot_count), function(slot_idx) {
            slot_entries <- list(
                medicago = panel_feature_entry_at(grouped_row$medicago_panels[[1]], slot_idx, repeat_singleton = TRUE),
                glycine = panel_feature_entry_at(grouped_row$glycine_panels[[1]], slot_idx, repeat_singleton = TRUE),
                lotus = panel_feature_entry_at(grouped_row$lotus_panels[[1]], slot_idx, repeat_singleton = TRUE)
            )

            tibble(
                panel_order = grouped_row$panel_order[[1]],
                orthogroup_label = grouped_row$orthogroup_label[[1]],
                slot_index = slot_idx,
                source_gene_display = grouped_row$source_gene_display[[1]],
                medicago_panels = list(slot_entries$medicago),
                glycine_panels = list(slot_entries$glycine),
                lotus_panels = list(slot_entries$lotus),
                title = cross_comparison_panel_title(
                    grouped_row$source_gene_display[[1]],
                    slot_entries
                )
            )
        }))
    })) %>%
        arrange(panel_order, slot_index)

    titles <- grouped_tbl$title
    duplicated_titles <- duplicated(titles) | duplicated(titles, fromLast = TRUE)

    if (any(duplicated_titles)) {
        titles[duplicated_titles & nzchar(grouped_tbl$orthogroup_label)] <- paste0(
            titles[duplicated_titles & nzchar(grouped_tbl$orthogroup_label)],
            " (",
            grouped_tbl$orthogroup_label[duplicated_titles & nzchar(grouped_tbl$orthogroup_label)],
            ")"
        )
        titles <- make.unique(titles, sep = " • ")
    }

grouped_tbl %>%
        mutate(title = titles)
}

cross_render_load_summary <- function(resolution, source_species, cross_key) {
    panel_specs <- cross_comparison_panel_specs(resolution, cross_key)
    multiplicity <- resolution$multiplicity %>%
        arrange(desc(mapped_gene_count), source_gene)

    expansion_labels <- if (nrow(multiplicity)) {
        paste0(
            display_gene_labels(
                source_species,
                multiplicity$source_gene,
                include_gene_id_with_common = FALSE
            ),
            " (",
            multiplicity$mapped_gene_count,
            " mapped ortholog plots)"
        )
    } else {
        character(0)
    }

    list(
        panel_count = nrow(panel_specs),
        mapped_feature_count = nrow(resolution$plot_table),
        expanded_source_count = nrow(multiplicity),
        expansion_labels = expansion_labels
    )
}

enforce_plot_feature_limit <- function(resolution, limit = atlas_plot_feature_limit) {
    limit <- suppressWarnings(as.integer(limit))
    if (!is.finite(limit) || limit <= 0L) {
        limit <- atlas_plot_feature_limit
    }

    all_features <- unique(as.character(resolution$plot_features %||% character(0)))
    all_features <- all_features[!is.na(all_features) & nzchar(all_features)]
    total_features <- length(all_features)

    resolution$feature_limit <- limit
    resolution$feature_count_before_limit <- total_features
    resolution$feature_limit_exceeded <- total_features > limit
    resolution$features_omitted <- character(0)

    if (!isTRUE(resolution$feature_limit_exceeded)) {
        return(resolution)
    }

    kept_features <- all_features[seq_len(limit)]
    omitted_features <- setdiff(all_features, kept_features)
    resolution$features_omitted <- omitted_features
    resolution$plot_features <- as.character(resolution$plot_features)
    resolution$plot_features <- resolution$plot_features[resolution$plot_features %in% kept_features]

    if (!is.null(resolution$plot_table) && nrow(resolution$plot_table)) {
        feature_col <- intersect(c("feature_id", "target_feature_id"), colnames(resolution$plot_table))[1]
        if (!is.na(feature_col)) {
            resolution$plot_table <- resolution$plot_table %>%
                filter(.data[[feature_col]] %in% kept_features)
        }
    }

    if (length(resolution$label_map)) {
        resolution$label_map <- resolution$label_map[names(resolution$label_map) %in% kept_features]
    }

    resolution
}

cross_comparison_height_units <- function(panel_specs) {
    if (is.null(panel_specs) || !nrow(panel_specs)) {
        return(integer(0))
    }

    pmax(
        1L,
        vapply(panel_specs$medicago_panels, panel_feature_count, integer(1)),
        vapply(panel_specs$glycine_panels, panel_feature_count, integer(1)),
        vapply(panel_specs$lotus_panels, panel_feature_count, integer(1))
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

gene_catalog_cache_path <- function(species_key, integration_method = NULL) {
    file.path(
        gene_catalog_cache_dir,
        paste0(species_key, ".tsv")
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

read_delimited_cache <- function(path) {
    if (!file.exists(path)) {
        return(NULL)
    }

    if (grepl("\\.csv$", path, ignore.case = TRUE)) {
        return(
            read.csv(
                path,
                stringsAsFactors = FALSE,
                na.strings = c("", "NA"),
                check.names = FALSE
            ) %>%
                as_tibble()
        )
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
    "Roots" = "#1B4332"
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

label_like_metadata_columns <- c(
    "cell_class",
    "cluster_label",
    "species_cell_class",
    "saturn_label",
    "saturn_ref_label",
    "Rank_1st_label",
    "Rank_2nd_label",
    "Rank_3rd_label",
    "Rank_4th_label",
    "Rank_5th_label"
)

distribution_cluster_columns <- c("Rank_1st", "Rank_2nd", "Rank_3rd", "Rank_4th", "Rank_5th", "cluster_label")

is_cluster_distribution_group <- function(column_name) {
    length(column_name) == 1L &&
        !is.na(column_name) &&
        column_name %in% c(distribution_cluster_columns, paste0(distribution_cluster_columns, "_label"))
}

stable_label_hash <- function(value) {
    chars <- utf8ToInt(enc2utf8(as.character(value)))

    if (!length(chars)) {
        return(0)
    }

    sum(chars * (seq_along(chars) + 17))
}

stable_label_palette <- function(values) {
    ordered_levels <- unique(as.character(values))
    ordered_levels <- ordered_levels[!is.na(ordered_levels) & nzchar(ordered_levels)]

    if (!length(ordered_levels)) {
        return(setNames(character(0), character(0)))
    }

    hashes <- vapply(ordered_levels, stable_label_hash, numeric(1))
    hues <- (hashes * 137.508) %% 360
    chroma_idx <- (hashes %% 4) + 1L
    luminance_idx <- ((hashes %/% 4) %% 3) + 1L
    chroma_values <- c(58, 66, 74, 62)[chroma_idx]
    luminance_values <- c(50, 58, 64)[luminance_idx]
    palette_values <- grDevices::hcl(
        h = hues,
        c = chroma_values,
        l = luminance_values,
        fixup = TRUE
    )

    setNames(palette_values, ordered_levels)
}

metadata_level_order <- function(values, column_name) {
    values_chr <- as.character(values)
    if (column_name %in% c("Group", "condition", "time_point")) {
        values_chr[values_chr == "roots"] <- "Roots"
    }
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

    if (identical(column_name, "species")) {
        return(c(
            vapply(within_species_keys, species_label, character(1)),
            sort(setdiff(present_levels, vapply(within_species_keys, species_label, character(1))))
        ))
    }

    if (column_name %in% c(
        "cluster_label",
        "Rank_1st", "Rank_2nd", "Rank_3rd", "Rank_4th", "Rank_5th",
        "Rank_1st_label", "Rank_2nd_label", "Rank_3rd_label", "Rank_4th_label", "Rank_5th_label"
    )) {
        return(cluster_value_levels(present_levels))
    }

    if (is.factor(values)) {
        return(levels(values)[levels(values) %in% present_levels])
    }

    sort(present_levels)
}

order_metadata_values <- function(values, column_name) {
    values_chr <- as.character(values)
    if (column_name %in% c("Group", "condition", "time_point")) {
        values_chr[values_chr == "roots"] <- "Roots"
    }
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

    if (column_name %in% c("species", "species_label")) {
        species_values <- species_palette[ordered_levels]
        if (all(!is.na(species_values))) {
            return(species_values)
        }

        if (any(!is.na(species_values))) {
            fallback_values <- stable_label_palette(ordered_levels[is.na(species_values)])
            species_values[is.na(species_values)] <- fallback_values[names(species_values)[is.na(species_values)]]
            return(species_values)
        }
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

    if (column_name %in% label_like_metadata_columns) {
        return(stable_label_palette(ordered_levels))
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
        palette_values <- metadata_colors_use(values, column_name)
        if (!is.null(palette_values)) return(palette_values)
    }

    palette_values <- scCustomize::scCustomize_Palette(
        num_groups = length(ordered_levels),
        color_seed = 123
    )
    names(palette_values) <- ordered_levels
    palette_values
}

apply_metadata_display_order <- function(obj, columns) {
    valid_columns <- unique(columns[columns %in% colnames(obj@meta.data)])

    if (!length(valid_columns)) {
        return(obj)
    }

    ordered_obj <- obj

    for (column_name in valid_columns) {
        current_values <- ordered_obj@meta.data[[column_name]]
        ordered_values <- order_metadata_values(current_values, column_name)

        if (identical(current_values, ordered_values)) {
            next
        }

        ordered_obj@meta.data[[column_name]] <- ordered_values
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

group_palette_for_toggle <- function(values, column_name, colorblind_safe = FALSE) {
    ordered_levels <- metadata_level_order(values, column_name)

    if (!length(ordered_levels)) {
        return(setNames(character(0), character(0)))
    }

    if (isTRUE(colorblind_safe)) {
        palette_values <- viridisLite::viridis(length(ordered_levels), end = 0.95)
        names(palette_values) <- ordered_levels
        return(palette_values)
    }

    distribution_color_map(values, column_name)
}

expression_heatmap_z_palette <- function(colorblind_safe = FALSE) {
    if (isTRUE(colorblind_safe)) {
        return(c(low = "#315a8c", mid = "#fff8df", high = "#b45f1f"))
    }

    c(low = "#315f72", mid = "#fff7d6", high = "#9b3f52")
}

row_scale_average_expression_matrix <- function(avg_expr_mat) {
    avg_expr_mat <- as.matrix(avg_expr_mat)
    scaled_mat <- t(apply(avg_expr_mat, 1, function(values) {
        missing_values <- is.na(values)

        if (all(missing_values)) {
            return(rep(NA_real_, length(values)))
        }

        gene_mean <- mean(values, na.rm = TRUE)
        gene_sd <- stats::sd(values, na.rm = TRUE)

        if (!is.finite(gene_sd) || gene_sd == 0) {
            scaled_values <- rep(0, length(values))
            scaled_values[missing_values] <- NA_real_
            return(scaled_values)
        }

        scaled_values <- (values - gene_mean) / gene_sd
        scaled_values[missing_values] <- NA_real_
        scaled_values
    }))

    rownames(scaled_mat) <- rownames(avg_expr_mat)
    colnames(scaled_mat) <- colnames(avg_expr_mat)
    scaled_mat
}

sample_cell_ids_by_group <- function(obj, group_by, max_cells = 2000L, seed = 123L) {
    dataset_key <- object_dataset_key(obj)
    cache_key <- if (!is.na(dataset_key) && nzchar(dataset_key)) {
        paste("sample_cells", group_by, max_cells, seed, sep = "::")
    } else {
        NA_character_
    }

    if (!is.na(cache_key) && nzchar(cache_key)) {
        entry <- get_dataset_cache_entry(dataset_key, function() annotate_dataset_key(obj, dataset_key))

        if (expression_cache_exists(entry, cache_key)) {
            return(expression_cache_get(entry, cache_key))
        }
    }

    cell_ids <- colnames(obj)

    if (!length(cell_ids)) {
        return(character(0))
    }

    if (!(group_by %in% colnames(obj@meta.data))) {
        return(cell_ids[seq_len(min(length(cell_ids), max_cells))])
    }

    sampling_df <- tibble(
        cell_id = cell_ids,
        group_value = as.character(obj@meta.data[[group_by]])
    ) %>%
        filter(!is.na(group_value) & nzchar(group_value))

    if (!nrow(sampling_df)) {
        return(cell_ids[seq_len(min(length(cell_ids), max_cells))])
    }

    sampled_df <- stratified_point_sample(
        sampling_df,
        group_col = "group_value",
        max_points = max_cells,
        seed = seed
    ) %>%
        mutate(group_value = order_metadata_values(group_value, group_by)) %>%
        arrange(group_value, cell_id)

    sampled_ids <- sampled_df$cell_id

    if (!is.na(cache_key) && nzchar(cache_key)) {
        entry <- get_dataset_cache_entry(dataset_key, function() annotate_dataset_key(obj, dataset_key))
        entry <- expression_cache_set(entry, cache_key, sampled_ids)
        set_dataset_cache_entry(dataset_key, entry)
    }

    sampled_ids
}

get_cached_sampled_object <- function(obj, group_by, max_cells = 2000L, seed = 123L) {
    sampled_cell_ids <- sample_cell_ids_by_group(
        obj,
        group_by = group_by,
        max_cells = max_cells,
        seed = seed
    )

    if (!length(sampled_cell_ids)) {
        return(list(
            object = obj[, sampled_cell_ids],
            cell_ids = sampled_cell_ids
        ))
    }

    dataset_key <- object_dataset_key(obj)

    if (is.na(dataset_key) || !nzchar(dataset_key)) {
        return(list(
            object = obj[, sampled_cell_ids],
            cell_ids = sampled_cell_ids
        ))
    }

    cache_key <- paste("sampled_object", group_by, max_cells, seed, sep = "::")
    entry <- get_dataset_cache_entry(dataset_key, function() annotate_dataset_key(obj, dataset_key))

    if (!expression_cache_exists(entry, cache_key)) {
        entry <- expression_cache_set(
            entry,
            cache_key,
            list(
                object = obj[, sampled_cell_ids],
                cell_ids = sampled_cell_ids
            )
        )
        set_dataset_cache_entry(dataset_key, entry)
    }

    expression_cache_get(entry, cache_key)
}

expression_feature_labels <- function(feature_ids, label_map = NULL) {
    feature_ids <- as.character(feature_ids %||% character(0))

    if (!length(feature_ids)) {
        return(character(0))
    }

    labels <- as.character(label_map[feature_ids] %||% feature_ids)
    labels[is.na(labels) | !nzchar(labels)] <- feature_ids[is.na(labels) | !nzchar(labels)]
    stats::setNames(labels, feature_ids)
}

build_expression_heatmap_plot <- function(
    obj,
    feature_ids,
    label_map,
    group_by,
    colorblind_safe = FALSE
) {
    if (!(group_by %in% colnames(obj@meta.data))) {
        stop("The selected grouping column is unavailable for the current heatmap selection.")
    }

    data_matrix <- get_cached_data_matrix(obj)
    feature_ids <- intersect(unique(as.character(feature_ids)), rownames(data_matrix))

    if (!length(feature_ids)) {
        stop("No mapped genes are available for the current heatmap selection.")
    }

    avg_expr_mat <- get_cached_average_expression_matrix(obj, feature_ids, group_by)
    scaled_expr_mat <- row_scale_average_expression_matrix(avg_expr_mat)
    group_levels <- colnames(scaled_expr_mat)

    if (!length(group_levels)) {
        stop("No groups are available for the current heatmap selection.")
    }

    feature_labels <- expression_feature_labels(feature_ids, label_map)
    feature_labels <- stats::setNames(
        make.unique(unname(feature_labels[feature_ids]), sep = " | "),
        feature_ids
    )

    heatmap_df <- as_tibble(scaled_expr_mat, rownames = "feature_id") %>%
        tidyr::pivot_longer(
            cols = -feature_id,
            names_to = "group_value",
            values_to = "scaled_expression"
        ) %>%
        mutate(
            group_value = factor(group_value, levels = group_levels, ordered = TRUE),
            feature_label = factor(
                feature_labels[feature_id],
                levels = rev(unname(feature_labels[feature_ids])),
                ordered = TRUE
            )
        )

    z_limit <- suppressWarnings(max(abs(heatmap_df$scaled_expression), na.rm = TRUE))
    if (!is.finite(z_limit) || z_limit <= 0) {
        z_limit <- 1
    }

    heatmap_palette <- expression_heatmap_z_palette(colorblind_safe)

    ggplot(heatmap_df, aes(x = group_value, y = feature_label, fill = scaled_expression)) +
        geom_tile(color = scales::alpha("white", 0.55), linewidth = 0.35) +
        scale_fill_gradient2(
            low = unname(heatmap_palette[["low"]]),
            mid = unname(heatmap_palette[["mid"]]),
            high = unname(heatmap_palette[["high"]]),
            midpoint = 0,
            limits = c(-z_limit, z_limit),
            oob = scales::squish,
            name = "Row-scaled average\nexpression (z-score)",
            na.value = "#f3f4ef"
        ) +
        labs(
            x = metadata_column_label(group_by),
            y = NULL
        ) +
        app_plot_theme() +
        theme(
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
            axis.ticks = element_blank(),
            legend.title = element_text(face = "bold"),
            plot.margin = margin(8, 14, 12, 10)
        )
}

add_species_caption <- function(plot_obj, species_key, sublabel = NULL) {
    caption_lines <- c(species_label(species_key))

    if (!is.null(sublabel) && !is.na(sublabel) && nzchar(sublabel)) {
        caption_lines <- c(caption_lines, as.character(sublabel))
    }

    plot_obj +
        labs(caption = paste(caption_lines, collapse = "\n")) +
        theme(
            plot.caption = element_text(
                colour = app_palette["muted"],
                size = 11,
                hjust = 0.5,
                margin = margin(t = 6)
            ),
            plot.caption.position = "plot"
        )
}

ortholog_trace_height_px <- function(gene_n) {
    max(
        620L,
        as.integer(220L + max(1L, gene_n) * 420L)
    )
}

within_integration_cell_join_key <- function(obj) {
    cell_ids <- colnames(obj)
    barcode_ids <- sub("^.*_", "", cell_ids)
    barcode_ids <- sub("(-[0-9]+)-[0-9]+$", "\\1", barcode_ids)
    sample_col <- pick_first_existing_col(obj@meta.data, c("Sample", "sample", "sample_name"))

    if (is.na(sample_col)) {
        return(stats::setNames(barcode_ids, cell_ids))
    }

    sample_values <- as.character(obj@meta.data[[sample_col]][match(cell_ids, rownames(obj@meta.data))])
    sample_values[is.na(sample_values) | !nzchar(sample_values)] <- "unknown"
    stats::setNames(paste(sample_values, barcode_ids, sep = "::"), cell_ids)
}

stratified_point_sample <- function(df, group_col, max_points = 30000L, seed = 123L) {
    if (!nrow(df) || nrow(df) <= max_points || !(group_col %in% colnames(df))) {
        return(df)
    }

    group_values <- as.character(df[[group_col]])
    group_values[is.na(group_values) | !nzchar(group_values)] <- "__missing__"
    split_indices <- split(seq_len(nrow(df)), group_values, drop = TRUE)
    group_sizes <- as.numeric(lengths(split_indices))

    if (!length(group_sizes)) {
        return(df)
    }

    raw_quota <- group_sizes * as.numeric(max_points) / sum(group_sizes)
    quota <- pmax(1, floor(raw_quota))
    quota <- pmin(quota, group_sizes)
    remainder <- as.integer(max_points - sum(quota))

    if (remainder > 0L) {
        spare_capacity <- group_sizes - quota
        order_idx <- order(raw_quota - quota, decreasing = TRUE, na.last = TRUE)

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

    quota <- as.integer(quota)

    set.seed(seed)
    sampled_indices <- purrr::map2(
        split_indices,
        quota,
        function(idx, n_keep) {
            n_keep <- as.integer(n_keep)
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
        identical(column_name, "Rank_2nd_label") ~ "Clustering opt 2 label",
        identical(column_name, "Rank_3rd_label") ~ "Clustering opt 3 label",
        identical(column_name, "Rank_4th_label") ~ "Clustering opt 4 label",
        identical(column_name, "Rank_5th_label") ~ "Clustering opt 5 label",
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

read_celltype_override <- function(override_id) {
    read_legacy_celltype_override(override_id, overrides_dir = celltype_overrides_dir)
}

apply_celltype_override <- function(obj, override_id) {
    apply_cluster_annotation_overlay(
        obj,
        dataset_key = override_id,
        annotations_dir = cluster_annotations_dir,
        legacy_overrides_dir = celltype_overrides_dir
    )
}

prepare_within_object <- function(obj, override_id = NULL) {
    obj <- apply_celltype_override(obj, override_id)
    obj <- apply_metadata_display_order(obj, colnames(obj@meta.data))
    obj
}

prepare_cross_object <- function(obj, cross_key) {
    obj <- apply_celltype_override(obj, cross_key)

    label_column <- dplyr::case_when(
        "celltype" %in% colnames(obj@meta.data) ~ "celltype",
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
    obj <- apply_metadata_display_order(obj, colnames(obj@meta.data))
    obj
}

expected_data_files <- function() {
    within_entries <- unlist(lapply(within_species_keys, function(species_key) {
        lapply(names(species_registry[[species_key]]$within_paths), function(integration_method) {
            list(
                path = pick_first_existing_path(c(
                    app_slim_path(species_registry[[species_key]]$within_paths[[integration_method]]),
                    species_registry[[species_key]]$within_paths[[integration_method]]
                )),
                label = sprintf(
                    "%s (%s integration)",
                    species_registry[[species_key]]$label,
                    integration_method
                )
            )
        })
    }), recursive = FALSE)

    cross_entries <- lapply(cross_integration_keys, function(cross_key) {
        list(
            path = pick_first_existing_path(c(
                cross_integration_registry[[cross_key]]$slim_path,
                cross_integration_registry[[cross_key]]$path
            )),
            label = sprintf("%s cross-species integration", cross_integration_registry[[cross_key]]$label)
        )
    })

    metadata_entries <- list(
        list(path = atlas_summary_path, label = "Atlas summary table"),
        list(path = cross_feature_lookup_path, label = "Cross-feature lookup")
    )

    gene_catalog_entries <- lapply(within_species_keys, function(species_key) {
        list(
            path = gene_catalog_cache_path(species_key),
            label = sprintf(
                "%s gene catalog (run `Rscript scripts/build_gene_catalog_cache.R` if missing)",
                species_registry[[species_key]]$label
            )
        )
    })

    c(within_entries, cross_entries, metadata_entries, gene_catalog_entries)
}

check_missing_data_files <- function() {
    entries <- expected_data_files()
    Filter(function(entry) !file.exists(entry$path), entries)
}

dataset_cache <- new.env(parent = emptyenv())
dataset_cache_order <- character(0)

cache_touch_key <- function(key) {
    key <- as.character(key)
    dataset_cache_order <<- c(setdiff(dataset_cache_order, key), key)
    invisible(key)
}

cache_prune_if_needed <- function() {
    max_entries <- atlas_cache_max_entries

    if (!is.finite(max_entries) || max_entries <= 0L) {
        return(invisible(FALSE))
    }

    current_keys <- ls(envir = dataset_cache, all.names = TRUE)
    dataset_cache_order <<- dataset_cache_order[dataset_cache_order %in% current_keys]

    overflow_n <- length(dataset_cache_order) - max_entries
    if (overflow_n <= 0L) {
        return(invisible(FALSE))
    }

    remove_keys <- dataset_cache_order[seq_len(overflow_n)]
    rm(list = remove_keys, envir = dataset_cache)
    dataset_cache_order <<- setdiff(dataset_cache_order, remove_keys)
    invisible(TRUE)
}

cache_assign <- function(key, value) {
    assign(key, value, envir = dataset_cache)
    cache_touch_key(key)
    cache_prune_if_needed()
    invisible(value)
}

cache_get <- function(key, builder) {
    if (!exists(key, envir = dataset_cache, inherits = FALSE)) {
        cache_assign(key, builder())
    } else {
        cache_touch_key(key)
    }

    get(key, envir = dataset_cache, inherits = FALSE)
}

new_dataset_cache_entry <- function(object) {
    list(
        object = object,
        data_matrix = NULL,
        expression_cache = new.env(parent = emptyenv()),
        expression_cache_order = character(0)
    )
}

ensure_dataset_cache_entry <- function(entry) {
    if (is.null(entry$expression_cache) || !is.environment(entry$expression_cache)) {
        entry$expression_cache <- new.env(parent = emptyenv())
    }

    if (is.null(entry$expression_cache_order)) {
        entry$expression_cache_order <- character(0)
    }

    entry
}

get_dataset_cache_entry <- function(key, builder) {
    if (!exists(key, envir = dataset_cache, inherits = FALSE)) {
        cache_assign(key, new_dataset_cache_entry(builder()))
    } else {
        cache_touch_key(key)
    }

    entry <- get(key, envir = dataset_cache, inherits = FALSE)
    entry <- ensure_dataset_cache_entry(entry)

    if (!identical(entry, get(key, envir = dataset_cache, inherits = FALSE))) {
        cache_assign(key, entry)
    }

    entry
}

set_dataset_cache_entry <- function(key, entry) {
    cache_assign(key, ensure_dataset_cache_entry(entry))
    invisible(entry)
}

expression_cache_touch_key <- function(entry, cache_key) {
    entry$expression_cache_order <- c(
        setdiff(entry$expression_cache_order %||% character(0), cache_key),
        cache_key
    )
    entry
}

expression_cache_prune_if_needed <- function(entry) {
    max_entries <- atlas_expression_cache_max_entries

    if (!is.finite(max_entries) || max_entries <= 0L) {
        return(entry)
    }

    current_keys <- ls(envir = entry$expression_cache, all.names = TRUE)
    entry$expression_cache_order <- (entry$expression_cache_order %||% character(0))
    entry$expression_cache_order <- entry$expression_cache_order[entry$expression_cache_order %in% current_keys]

    overflow_n <- length(entry$expression_cache_order) - max_entries
    if (overflow_n <= 0L) {
        return(entry)
    }

    remove_keys <- entry$expression_cache_order[seq_len(overflow_n)]
    rm(list = remove_keys, envir = entry$expression_cache)
    entry$expression_cache_order <- setdiff(entry$expression_cache_order, remove_keys)
    entry
}

expression_cache_exists <- function(entry, cache_key) {
    exists(cache_key, envir = entry$expression_cache, inherits = FALSE)
}

expression_cache_get <- function(entry, cache_key) {
    entry <- expression_cache_touch_key(entry, cache_key)
    get(cache_key, envir = entry$expression_cache, inherits = FALSE)
}

expression_cache_set <- function(entry, cache_key, value) {
    assign(cache_key, value, envir = entry$expression_cache)
    entry <- expression_cache_touch_key(entry, cache_key)
    expression_cache_prune_if_needed(entry)
}

annotate_dataset_key <- function(obj, dataset_key) {
    obj@misc$atlas_dataset_key <- dataset_key
    obj
}

object_dataset_key <- function(obj) {
    key <- obj@misc$atlas_dataset_key %||% NA_character_

    if (!length(key) || is.na(key) || !nzchar(key)) {
        return(NA_character_)
    }

    as.character(key[[1]])
}

hash_cache_values <- function(values, empty_key = "__all__") {
    values <- unique(as.character(values))
    values <- values[!is.na(values) & nzchar(values)]

    if (!length(values)) {
        return(empty_key)
    }

    as.character(rlang::hash(sort(values)))
}

compute_feature_expression_values <- function(data_matrix, feature_ids, cell_ids = NULL) {
    feature_ids <- unique(as.character(feature_ids))
    feature_ids <- feature_ids[!is.na(feature_ids) & nzchar(feature_ids)]
    feature_ids <- intersect(feature_ids, rownames(data_matrix))

    if (!length(feature_ids)) {
        return(setNames(numeric(0), character(0)))
    }

    cell_ids <- if (is.null(cell_ids)) {
        colnames(data_matrix)
    } else {
        intersect(as.character(cell_ids), colnames(data_matrix))
    }

    if (!length(cell_ids)) {
        return(setNames(numeric(0), character(0)))
    }

    expr_mat <- data_matrix[feature_ids, cell_ids, drop = FALSE]
    expr_values <- if (length(feature_ids) == 1L) {
        as.numeric(expr_mat[1, ])
    } else {
        as.numeric(Matrix::colMeans(expr_mat))
    }

    stats::setNames(expr_values, cell_ids)
}

extract_assay_data_matrix <- function(obj) {
    default_assay <- DefaultAssay(obj)
    assay <- obj[[default_assay]]

    if (inherits(assay, "Assay5") && "layers" %in% slotNames(assay) && "data" %in% names(assay@layers)) {
        data_matrix <- assay@layers[["data"]]
        methods::slot(data_matrix, "Dimnames") <- list(rownames(obj), colnames(obj))
        return(data_matrix)
    }

    GetAssayData(obj, layer = "data")
}

get_cached_data_matrix <- function(obj) {
    dataset_key <- object_dataset_key(obj)

    if (is.na(dataset_key) || !nzchar(dataset_key)) {
        return(extract_assay_data_matrix(obj))
    }

    entry <- get_dataset_cache_entry(dataset_key, function() annotate_dataset_key(obj, dataset_key))

    if (is.null(entry$data_matrix)) {
        entry$data_matrix <- extract_assay_data_matrix(entry$object)
        set_dataset_cache_entry(dataset_key, entry)
    }

    entry$data_matrix
}

get_cached_feature_expression <- function(obj, feature_ids, cell_ids = NULL) {
    dataset_key <- object_dataset_key(obj)
    data_matrix <- get_cached_data_matrix(obj)
    feature_ids <- unique(as.character(feature_ids))
    feature_ids <- feature_ids[!is.na(feature_ids) & nzchar(feature_ids)]
    feature_ids <- intersect(feature_ids, rownames(data_matrix))

    if (!length(feature_ids)) {
        return(setNames(numeric(0), character(0)))
    }

    if (is.null(cell_ids)) {
        cell_ids <- colnames(data_matrix)
        cell_hash <- "__all__"
    } else {
        cell_ids <- intersect(as.character(cell_ids), colnames(data_matrix))
        cell_hash <- hash_cache_values(cell_ids, empty_key = "__none__")
    }

    if (!length(cell_ids)) {
        return(setNames(numeric(0), character(0)))
    }

    if (is.na(dataset_key) || !nzchar(dataset_key)) {
        return(compute_feature_expression_values(data_matrix, feature_ids, cell_ids))
    }

    entry <- get_dataset_cache_entry(dataset_key, function() annotate_dataset_key(obj, dataset_key))
    cache_key <- paste(hash_cache_values(feature_ids), cell_hash, sep = "::")

    if (!expression_cache_exists(entry, cache_key)) {
        entry <- expression_cache_set(
            entry,
            cache_key,
            compute_feature_expression_values(data_matrix, feature_ids, cell_ids)
        )
        set_dataset_cache_entry(dataset_key, entry)
    }

    expression_cache_get(entry, cache_key)
}

compute_average_expression_matrix <- function(data_matrix, metadata, feature_ids, group_by) {
    feature_ids <- intersect(unique(as.character(feature_ids)), rownames(data_matrix))

    if (!length(feature_ids)) {
        return(matrix(numeric(0), nrow = 0L, ncol = 0L))
    }

    if (!(group_by %in% colnames(metadata))) {
        stop("The selected grouping column is unavailable for this expression summary.")
    }

    group_values <- order_metadata_values(metadata[[group_by]], group_by)
    valid_idx <- !is.na(group_values) & nzchar(as.character(group_values))
    valid_cells <- colnames(data_matrix)[valid_idx]

    if (!length(valid_cells)) {
        return(matrix(numeric(0), nrow = length(feature_ids), ncol = 0L))
    }

    group_factor <- droplevels(group_values[valid_idx])
    group_levels <- levels(group_factor)

    avg_expr_mat <- vapply(group_levels, function(group_value) {
        group_cell_ids <- valid_cells[group_factor == group_value]

        if (!length(group_cell_ids)) {
            return(rep(NA_real_, length(feature_ids)))
        }

        group_expr_mat <- data_matrix[feature_ids, group_cell_ids, drop = FALSE]

        if (length(group_cell_ids) == 1L) {
            return(as.numeric(group_expr_mat[, 1, drop = TRUE]))
        }

        as.numeric(Matrix::rowMeans(group_expr_mat))
    }, numeric(length(feature_ids)))

    if (is.null(dim(avg_expr_mat))) {
        avg_expr_mat <- matrix(avg_expr_mat, nrow = length(feature_ids), ncol = length(group_levels))
    }

    rownames(avg_expr_mat) <- feature_ids
    colnames(avg_expr_mat) <- group_levels
    avg_expr_mat
}

get_cached_average_expression_matrix <- function(obj, feature_ids, group_by) {
    dataset_key <- object_dataset_key(obj)
    data_matrix <- get_cached_data_matrix(obj)
    feature_ids <- intersect(unique(as.character(feature_ids)), rownames(data_matrix))

    if (!length(feature_ids)) {
        return(matrix(numeric(0), nrow = 0L, ncol = 0L))
    }

    if (is.na(dataset_key) || !nzchar(dataset_key)) {
        return(compute_average_expression_matrix(data_matrix, obj@meta.data, feature_ids, group_by))
    }

    entry <- get_dataset_cache_entry(dataset_key, function() annotate_dataset_key(obj, dataset_key))
    cache_key <- paste("avg_expr", group_by, hash_cache_values(feature_ids), sep = "::")

    if (!expression_cache_exists(entry, cache_key)) {
        entry <- expression_cache_set(
            entry,
            cache_key,
            compute_average_expression_matrix(data_matrix, obj@meta.data, feature_ids, group_by)
        )
        set_dataset_cache_entry(dataset_key, entry)
    }

    expression_cache_get(entry, cache_key)
}

safe_read_rds <- function(path, label) {
    if (!length(path) || is.na(path) || !nzchar(path)) {
        stop(sprintf("No path configured for %s.", label), call. = FALSE)
    }

    if (!file.exists(path)) {
        stop(
            sprintf(
                "%s is not available: the expected file '%s' was not found. If you are running a local copy, confirm the data directory is mounted correctly.",
                label, path
            ),
            call. = FALSE
        )
    }

    tryCatch(
        readRDS(path),
        error = function(e) {
            stop(
                sprintf(
                    "Failed to load %s from '%s'. The file may be corrupt or was saved with an incompatible R/Seurat version. Original error: %s",
                    label, path, conditionMessage(e)
                ),
                call. = FALSE
            )
        }
    )
}

get_within_object <- function(species_key, integration_method) {
    cache_key <- paste("within", species_key, integration_method, sep = "::")

    get_dataset_cache_entry(cache_key, function() {
        label <- sprintf(
            "%s (%s integration)",
            species_registry[[species_key]]$label,
            integration_method
        )
        safe_read_rds(
            pick_first_existing_path(c(
                app_slim_path(species_registry[[species_key]]$within_paths[[integration_method]]),
                species_registry[[species_key]]$within_paths[[integration_method]]
            )),
            label
        ) %>%
            prepare_within_object(override_id = paste(species_key, integration_method, sep = "_")) %>%
            annotate_dataset_key(cache_key)
    })$object
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

read_cached_umap3d <- function(path, label) {
    if (!length(path) || is.na(path) || !nzchar(path) || !file.exists(path)) {
        return(NULL)
    }

    embedding <- tryCatch(
        safe_read_rds(path, label),
        error = function(e) NULL
    )

    if (is.null(embedding)) {
        return(NULL)
    }

    embedding <- as.matrix(embedding)

    if (ncol(embedding) < 3L || is.null(rownames(embedding))) {
        return(NULL)
    }

    embedding[, seq_len(3), drop = FALSE]
}

get_cross_umap3d <- function(cross_key) {
    cache_key <- paste("cross_umap3d", cross_key, sep = "::")

    cache_get(cache_key, function() {
        obj <- get_cross_object(cross_key)
        available_reductions <- Reductions(obj)
        reduction_name <- c("umap3d", "umap_3d", within_three_d_reduction_name)[
            c("umap3d", "umap_3d", within_three_d_reduction_name) %in% available_reductions
        ][1]

        if (!is.na(reduction_name) && length(reduction_name) && ncol(Embeddings(obj, reduction_name)) >= 3L) {
            return(Embeddings(obj, reduction_name)[, seq_len(3), drop = FALSE])
        }

        cached_embedding <- read_cached_umap3d(
            cross_integration_registry[[cross_key]]$umap3d_path %||% "",
            sprintf("%s 3D UMAP embedding", cross_integration_registry[[cross_key]]$label)
        )

        if (!is.null(cached_embedding)) {
            return(cached_embedding)
        }

        tryCatch(compute_umap3d_matrix(obj), error = function(e) NULL)
    })
}

get_cross_object <- function(cross_key) {
    cache_key <- paste("cross::object", cross_key, sep = "::")

    get_dataset_cache_entry(cache_key, function() {
        label <- sprintf(
            "%s cross-species integration",
            cross_integration_registry[[cross_key]]$label
        )
        safe_read_rds(
            pick_first_existing_path(c(
                cross_integration_registry[[cross_key]]$slim_path,
                cross_integration_registry[[cross_key]]$path
            )),
            label
        ) %>%
            prepare_cross_object(cross_key) %>%
            annotate_dataset_key(cache_key)
    })$object
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
            cached_lookup <- read_tsv_cache(cross_feature_lookup_path)

            if (!is.null(cached_lookup) &&
                all(c("feature_id", "canonical_gene_id") %in% colnames(cached_lookup))) {
                return(
                    cached_lookup %>%
                        select(feature_id, canonical_gene_id) %>%
                        mutate(feature_species = "medicago") %>%
                        filter(!is.na(feature_id) & nzchar(feature_id)) %>%
                        filter(!is.na(canonical_gene_id) & nzchar(canonical_gene_id)) %>%
                        distinct(feature_species, canonical_gene_id, feature_id)
                )
            }

            feature_ids <- rownames(get_cross_object(cross_key))

            return(
                tibble(
                    feature_id = feature_ids,
                    canonical_gene_id = canonicalize_gene_ids("medicago", feature_ids),
                    feature_species = "medicago"
                ) %>%
                    distinct(feature_species, canonical_gene_id, feature_id)
            )
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
    source_candidates <- source_orthogroups %>%
        filter(!is.na(orthogroup))

    if (!nrow(source_candidates)) {
        return(list(
            mapping = tibble(),
            plot_table = tibble(),
            plot_features = character(0),
            label_map = c(),
            no_orthogroup = sort(unique(no_orthogroup)),
            no_target_members = character(0),
            missing_features = character(0),
            multiplicity = tibble()
        ))
    }

    candidate_tbl <- source_candidates %>%
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
            cached_catalog <- cached_catalog %>%
                mutate(across(everything(), as.character))

            if (identical(species_key, "glycine")) {
                cached_catalog <- cached_catalog %>%
                    mutate(
                        display_label = display_gene_labels(species_key, feature_id),
                        search_tokens = purrr::pmap_chr(
                            list(feature_id, display_label, search_tokens),
                            function(feature_id, display_label, search_tokens) {
                                gene_search_tokens(
                                    species_key = species_key,
                                    feature_id = feature_id,
                                    display_label = display_label,
                                    existing_tokens = search_tokens
                                )
                            }
                        )
                    )
            }

            return(cached_catalog)
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
                    list(feature_id, display_label, common_name, synonyms),
                    function(feature_id, display_label, common_name, synonyms) {
                        gene_search_tokens(
                            species_key = species_key,
                            feature_id = feature_id,
                            display_label = display_label,
                            common_name = common_name,
                            synonyms = synonyms
                        )
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
            choices = catalog %>%
                transmute(
                    value = feature_id,
                    label = display_label,
                    tokens = search_tokens
                ),
            tokens = catalog$search_tokens,
            feature_ids = catalog$feature_id
        )
    })
}

match_gene_choices <- function(choice_bundle, query) {
    query <- trimws(as.character(query %||% ""))
    if (!nzchar(query)) {
        return(character(0))
    }

    choices <- as_tibble(choice_bundle$choices %||% tibble())
    if (!nrow(choices)) {
        return(character(0))
    }

    query_terms <- unlist(strsplit(query, "\\s+", perl = TRUE), use.names = FALSE)
    query_terms <- trimws(query_terms)
    query_terms <- query_terms[nzchar(query_terms)]

    if (!length(query_terms)) {
        return(character(0))
    }

    haystack <- paste(
        choices$label %||% "",
        choices$value %||% "",
        choices$tokens %||% "",
        sep = " "
    )

    haystack_lower <- tolower(haystack)

    matched <- Reduce(
        `|`,
        lapply(query_terms, function(term) {
            grepl(tolower(term), haystack_lower, fixed = TRUE)
        })
    )

    choices$value[matched]
}

cluster_markers_cache_path <- function(dataset_key, top_n = FALSE) {
    suffix <- if (isTRUE(top_n)) "_top10" else ""
    file.path(cluster_markers_cache_dir, paste0(dataset_key, suffix, ".tsv"))
}

cluster_markers_cache_candidates <- function(dataset_key, top_n = FALSE) {
    suffix <- if (isTRUE(top_n)) "_top10" else ""
    base_path <- file.path(cluster_markers_cache_dir, paste0(dataset_key, suffix))

    c(
        paste0(base_path, ".csv"),
        paste0(base_path, ".tsv")
    )
}

cluster_markers_split_cache_paths <- function(dataset_key, top_n = FALSE) {
    if (!dir.exists(cluster_markers_cache_dir)) {
        return(character(0))
    }

    paths <- list.files(
        cluster_markers_cache_dir,
        full.names = TRUE,
        ignore.case = TRUE
    )

    if (!length(paths)) {
        return(character(0))
    }

    expected_prefix <- paste0(dataset_key, "_markers_")
    paths <- paths[startsWith(basename(paths), expected_prefix)]
    paths <- paths[grepl("\\.(csv|tsv)$", basename(paths), ignore.case = TRUE)]

    if (isTRUE(top_n)) {
        paths <- paths[grepl("top10", basename(paths), ignore.case = TRUE)]
    } else {
        paths <- paths[!grepl("top10", basename(paths), ignore.case = TRUE)]
    }

    sort(paths)
}

marker_cluster_source_from_path <- function(path, dataset_key) {
    file_name <- basename(path)

    prefix <- paste0(dataset_key, "_markers_")
    source_part <- if (startsWith(file_name, prefix)) {
        substr(file_name, nchar(prefix) + 1L, nchar(file_name))
    } else {
        file_name
    }

    source_part <- sub("\\.(csv|tsv)$", "", source_part, ignore.case = TRUE)
    source_part <- sub("_top10$", "", source_part, ignore.case = TRUE)
    source_part <- sub("^_+", "", source_part)

    if (!nzchar(source_part)) {
        return("__idents__")
    }

    normalize_marker_cluster_source(source_part)
}

normalize_marker_cluster_source <- function(values) {
    values_chr <- trimws(as.character(values))
    values_key <- tolower(values_chr)

    dplyr::case_when(
        values_key %in% c("", "na") ~ "__idents__",
        values_key %in% c("__idents__", "saved identities", "saved identity", "idents", "identities", "identity") ~ "__idents__",
        values_key %in% c("rank_1st", "rank 1st", "clustering opt 1", "clustering option 1") ~ "Rank_1st",
        values_key %in% c("rank_1st_label", "rank 1st label", "clustering opt 1 label", "clustering option 1 label") ~ "Rank_1st_label",
        values_key %in% c("rank_2nd", "rank 2nd", "clustering opt 2", "clustering option 2") ~ "Rank_2nd",
        values_key %in% c("rank_3rd", "rank 3rd", "clustering opt 3", "clustering option 3") ~ "Rank_3rd",
        values_key %in% c("rank_4th", "rank 4th", "clustering opt 4", "clustering option 4") ~ "Rank_4th",
        values_key %in% c("rank_5th", "rank 5th", "clustering opt 5", "clustering option 5") ~ "Rank_5th",
        values_key %in% c("cluster", "cluster label") ~ "cluster_label",
        values_key %in% c("species + label", "species_cell_class") ~ "species_cell_class",
        TRUE ~ values_chr
    )
}

normalize_cluster_markers_cache <- function(marker_tbl, source_hint = NULL) {
    if (is.null(marker_tbl) || !nrow(marker_tbl)) {
        return(NULL)
    }

    pick_marker_col <- function(candidates) {
        matches <- candidates[candidates %in% colnames(marker_tbl)]
        if (!length(matches)) {
            return(NA_character_)
        }
        matches[[1]]
    }

    cluster_col <- pick_marker_col(c("cluster", "Cluster", "cluster_id", "ident", "identity"))
    gene_col <- pick_marker_col(c("gene", "Gene ID", "feature", "feature_id"))
    cluster_source_col <- pick_marker_col(c(
        "cluster_source",
        "clustering_level",
        "clustering_column",
        "cluster_level",
        "cluster_scheme",
        "group_by",
        "cluster_by"
    ))
    logfc_col <- pick_marker_col(c("avg_log2FC", "avg_logFC", "avg_log10FC", "Average log2(FoldChange)"))
    pct1_col <- pick_marker_col(c("pct.1", "pct1", "pct_1", "% of cells expressing gene in cluster"))
    pct2_col <- pick_marker_col(c("pct.2", "pct2", "pct_2", "% of cells expressing gene in remaining clusters"))
    pval_adj_col <- pick_marker_col(c("p_val_adj", "padj", "p_adj", "Adjusted p-value"))

    required_cols <- c(cluster_col, gene_col, logfc_col, pct1_col, pct2_col, pval_adj_col)

    if (any(is.na(required_cols))) {
        return(NULL)
    }

    cluster_source_values <- if (!is.na(cluster_source_col)) {
        as.character(marker_tbl[[cluster_source_col]])
    } else if (!is.null(source_hint) && length(source_hint) && !is.na(source_hint) && nzchar(source_hint)) {
        rep(source_hint, nrow(marker_tbl))
    } else {
        rep("__idents__", nrow(marker_tbl))
    }

    marker_tbl %>%
        transmute(
            cluster_source = normalize_marker_cluster_source(cluster_source_values),
            cluster = trimws(as.character(.data[[cluster_col]])),
            gene = trimws(as.character(.data[[gene_col]])),
            avg_log2FC = as.numeric(.data[[logfc_col]]),
            pct.1 = as.numeric(.data[[pct1_col]]),
            pct.2 = as.numeric(.data[[pct2_col]]),
            p_val_adj = as.numeric(.data[[pval_adj_col]])
        ) %>%
        mutate(
            cluster_source = ifelse(is.na(cluster_source) | !nzchar(cluster_source), "__idents__", cluster_source)
        ) %>%
        filter(
            !is.na(cluster_source) & nzchar(cluster_source),
            !is.na(cluster) & nzchar(cluster),
            !is.na(gene) & nzchar(gene)
        )
}

marker_cluster_source_label <- function(cluster_source) {
    if (is.null(cluster_source) || !length(cluster_source) || is.na(cluster_source) || !nzchar(cluster_source)) {
        return("marker clusters")
    }

    if (identical(cluster_source, "__idents__")) {
        return("saved identities")
    }

    metadata_column_label(cluster_source)
}

select_marker_cluster_source <- function(available_sources, preferred_source = NULL, default_source = "Rank_1st") {
    available_sources <- unique(as.character(available_sources))
    available_sources <- available_sources[!is.na(available_sources) & nzchar(available_sources)]

    if (!length(available_sources)) {
        return(NA_character_)
    }

    if (!is.null(preferred_source) && length(preferred_source) && !is.na(preferred_source) && nzchar(preferred_source)) {
        if (preferred_source %in% available_sources) {
            return(preferred_source)
        }

        preferred_label_source <- paste0(preferred_source, "_label")
        if (preferred_label_source %in% available_sources) {
            return(preferred_label_source)
        }
    }

    if (!is.null(default_source) && length(default_source) && !is.na(default_source) && nzchar(default_source)) {
        if (default_source %in% available_sources) {
            return(default_source)
        }

        default_label_source <- paste0(default_source, "_label")
        if (default_label_source %in% available_sources) {
            return(default_label_source)
        }
    }

    if ("__idents__" %in% available_sources) {
        return("__idents__")
    }

    available_sources[[1]]
}

marker_cluster_source_matches <- function(cluster_source, preferred_source) {
    if (is.null(cluster_source) || !length(cluster_source) || is.na(cluster_source) || !nzchar(cluster_source)) {
        return(FALSE)
    }

    if (is.null(preferred_source) || !length(preferred_source) || is.na(preferred_source) || !nzchar(preferred_source)) {
        return(FALSE)
    }

    identical(cluster_source, preferred_source) ||
        identical(cluster_source, paste0(preferred_source, "_label")) ||
        identical(preferred_source, paste0(cluster_source, "_label"))
}

read_cluster_markers_cache <- function(dataset_key, top_n = FALSE) {
    cache_key <- paste("cluster_markers", dataset_key, if (isTRUE(top_n)) "top10" else "full", sep = "::")

    cache_get(cache_key, function() {
        direct_candidates <- cluster_markers_cache_candidates(dataset_key, top_n = top_n)
        direct_path <- direct_candidates[file.exists(direct_candidates)][1]

        if (length(direct_path) && !is.na(direct_path) && nzchar(direct_path)) {
            cached_markers <- read_delimited_cache(direct_path)
            return(normalize_cluster_markers_cache(cached_markers))
        }

        split_paths <- cluster_markers_split_cache_paths(dataset_key, top_n = top_n)

        if (!length(split_paths)) {
            return(NULL)
        }

        marker_tables <- lapply(split_paths, function(path) {
            normalize_cluster_markers_cache(
                read_delimited_cache(path),
                source_hint = marker_cluster_source_from_path(path, dataset_key)
            )
        })

        marker_tables <- Filter(function(tbl) !is.null(tbl) && nrow(tbl), marker_tables)

        if (!length(marker_tables)) {
            return(NULL)
        }

        bind_rows(marker_tables) %>%
            distinct(cluster_source, cluster, gene, .keep_all = TRUE)
    })
}

cluster_label_lookup <- function(obj, cluster_column = "__idents__") {
    if (!is.null(cluster_column) &&
        length(cluster_column) &&
        !identical(cluster_column, "__idents__") &&
        cluster_column %in% colnames(obj@meta.data)) {
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

read_ui_choice_cache <- function() {
    cache_tbl <- read_tsv_cache(ui_choice_cache_path)

    if (is.null(cache_tbl) ||
        !all(c("dataset_key", "choice_family", "label", "value", "sort_order") %in% colnames(cache_tbl))) {
        return(NULL)
    }

    cache_tbl %>%
        transmute(
            dataset_key = as.character(dataset_key),
            choice_family = as.character(choice_family),
            label = as.character(label),
            value = as.character(value),
            sort_order = suppressWarnings(as.integer(sort_order))
        ) %>%
        filter(
            !is.na(dataset_key) & nzchar(dataset_key),
            !is.na(choice_family) & nzchar(choice_family),
            !is.na(label) & nzchar(label),
            !is.na(value) & nzchar(value)
        ) %>%
        mutate(sort_order = ifelse(is.na(sort_order), row_number(), sort_order))
}

get_ui_choice_cache <- function() {
    cache_get("ui::choice_cache", function() {
        read_ui_choice_cache()
    })
}

get_cached_ui_choices <- function(dataset_key, choice_family) {
    cache_tbl <- get_ui_choice_cache()

    if (is.null(cache_tbl) || !nrow(cache_tbl)) {
        return(NULL)
    }

    rows <- cache_tbl %>%
        filter(
            .data$dataset_key == .env$dataset_key,
            .data$choice_family == .env$choice_family
        ) %>%
        arrange(.data$sort_order, .data$label)

    if (!nrow(rows)) {
        return(NULL)
    }

    setNames(rows$value, rows$label)
}

read_ui_cluster_lookup_cache <- function() {
    cache_tbl <- read_tsv_cache(ui_cluster_lookup_cache_path)

    if (is.null(cache_tbl) ||
        !all(c("dataset_key", "cluster_source", "cluster", "cluster_label", "choice_label", "sort_order") %in% colnames(cache_tbl))) {
        return(NULL)
    }

    cache_tbl %>%
        transmute(
            dataset_key = as.character(dataset_key),
            cluster_source = normalize_marker_cluster_source(cluster_source),
            cluster = as.character(cluster),
            cluster_label = as.character(cluster_label),
            choice_label = as.character(choice_label),
            sort_order = suppressWarnings(as.integer(sort_order))
        ) %>%
        filter(
            !is.na(dataset_key) & nzchar(dataset_key),
            !is.na(cluster_source) & nzchar(cluster_source),
            !is.na(cluster) & nzchar(cluster)
        ) %>%
        mutate(
            cluster_label = ifelse(is.na(cluster_label) | !nzchar(cluster_label), cluster, cluster_label),
            choice_label = ifelse(is.na(choice_label) | !nzchar(choice_label), cluster_label, choice_label),
            sort_order = ifelse(is.na(sort_order), row_number(), sort_order)
        )
}

get_ui_cluster_lookup_cache <- function() {
    cache_get("ui::cluster_lookup_cache", function() {
        read_ui_cluster_lookup_cache()
    })
}

get_cached_cluster_lookup <- function(dataset_key, cluster_source) {
    cache_tbl <- get_ui_cluster_lookup_cache()
    cluster_source <- normalize_marker_cluster_source(cluster_source)

    if (is.null(cache_tbl) || !nrow(cache_tbl) ||
        is.null(cluster_source) || !length(cluster_source) || is.na(cluster_source) || !nzchar(cluster_source)) {
        return(NULL)
    }

    rows <- cache_tbl %>%
        filter(
            .data$dataset_key == .env$dataset_key,
            .data$cluster_source == .env$cluster_source
        ) %>%
        arrange(.data$sort_order, .data$cluster)

    if (!nrow(rows)) {
        return(NULL)
    }

    rows %>%
        transmute(
            cluster = .data$cluster,
            cluster_label = .data$cluster_label,
            choice_label = .data$choice_label
        )
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

build_metadata_axis_summary <- function(values, label, limit = 4) {
    clean_values <- unique(trimws(as.character(values)))
    clean_values <- clean_values[!is.na(clean_values) & nzchar(clean_values)]

    list(
        label = label,
        n = length(clean_values),
        preview = if (length(clean_values)) {
            summarise_value_preview(clean_values, limit = limit)
        } else {
            NA_character_
        }
    )
}

pick_metadata_axis_summary <- function(md, candidates, limit = 4) {
    for (candidate in candidates) {
        column_name <- pick_first_existing_col(md, candidate$columns)
        if (!is.na(column_name)) {
            return(build_metadata_axis_summary(
                values = md[[column_name]],
                label = candidate$label,
                limit = limit
            ))
        }
    }

    list(label = NA_character_, n = NA_integer_, preview = NA_character_)
}

summary_value <- function(summary, field) {
    value <- summary[[field]]
    if (!length(value)) return(NA)
    value[[1]]
}

summary_strip_tile <- function(label_tag, summary) {
    div(
        class = "summary-strip-tile",
        div(class = "summary-strip-label", label_tag),
        div(class = "summary-strip-value", format_stat_value(summary_value(summary, "cells"))),
        div(
            class = "summary-strip-note",
            sprintf(
                "%s genes \u00b7 %s samples",
                format_stat_value(summary_value(summary, "genes")),
                format_stat_value(summary_value(summary, "sample_n"))
            )
        )
    )
}

atlas_summary_cache_is_current <- function(summary_df) {
    required_cols <- c(
        "dataset_scope",
        "species_key",
        "integration_method",
        "cells",
        "genes",
        "sample_n"
    )

    nrow(summary_df) > 0 &&
        all(required_cols %in% colnames(summary_df))
}

build_within_dataset_summary_row <- function(species_key, integration_method, obj) {
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
    species_summary <- build_metadata_axis_summary(species_label(species_key), label = "species", limit = 1)

    tibble(
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
        time_preview = time_summary$preview
    )
}

build_cross_dataset_summary_row <- function(cross_key, obj) {
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

    tibble(
        dataset_scope = "cross",
        species_key = cross_key,
        integration_method = cross_key,
        integration_label = cross_integration_label(cross_key),
        cells = ncol(obj),
        genes = nrow(obj),
        sample_n = length(unique(sample_values[!is.na(sample_values) & nzchar(sample_values)])),
        group_label = group_summary$label,
        group_n = group_summary$n,
        group_preview = group_summary$preview,
        species_n = species_summary$n,
        species_preview = species_summary$preview,
        time_n = time_summary$n,
        time_preview = time_summary$preview
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
                all(required_cross_keys %in% cached_cross_keys) &&
                atlas_summary_cache_is_current(summary_df)) {
                return(summary_df)
            }
        }

        compute_atlas_summary_table()
    })
}

get_within_dataset_summary <- function(species_key, integration_method) {
    summary_row <- get_atlas_summary_table() %>%
        filter(
            .data$dataset_scope == "within",
            .data$species_key == .env$species_key,
            .data$integration_method == .env$integration_method
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
            .data$dataset_scope == "cross",
            .data$integration_method == .env$cross_key
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

if (!file.exists(orthogroups_path)) {
    stop(sprintf(
        "Orthogroups table not found at '%s'. Set ATLAS_ORTHOLOGS_PATH to a TSV with columns 'Orthogroup', '%s', '%s', '%s' to use a custom table.",
        orthogroups_path,
        species_registry$medicago$orthogroup_col,
        species_registry$glycine$orthogroup_col,
        species_registry$lotus$orthogroup_col
    ), call. = FALSE)
}

orthogroup_sep <- if (grepl("\\.csv$", orthogroups_path, ignore.case = TRUE)) "," else "\t"

orthogroup_raw <- read.delim(
    orthogroups_path,
    sep = orthogroup_sep,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

required_ortho_cols <- c(
    "Orthogroup",
    species_registry$medicago$orthogroup_col,
    species_registry$glycine$orthogroup_col,
    species_registry$lotus$orthogroup_col
)
missing_ortho_cols <- setdiff(required_ortho_cols, colnames(orthogroup_raw))
if (length(missing_ortho_cols)) {
    stop(sprintf(
        "Orthogroups table at '%s' is missing required columns: %s",
        orthogroups_path,
        paste(missing_ortho_cols, collapse = ", ")
    ), call. = FALSE)
}

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

map_target_genes_to_source_genes <- function(target_species, target_genes, source_species) {
    target_genes <- unique(as.character(target_genes))
    target_genes <- target_genes[!is.na(target_genes) & nzchar(target_genes)]

    if (!length(target_genes) || !(target_species %in% within_species_keys) || !(source_species %in% within_species_keys)) {
        return(character(0))
    }

    if (identical(target_species, source_species)) {
        return(target_genes)
    }

    tibble(
        target_gene = target_genes,
        target_canonical = canonicalize_gene_ids(target_species, target_genes)
    ) %>%
        filter(!is.na(target_canonical) & nzchar(target_canonical)) %>%
        left_join(
            orthogroup_index %>% filter(species == target_species),
            by = c("target_canonical" = "canonical_gene_id")
        ) %>%
        mutate(
            source_members = purrr::map(orthogroup, function(og) {
                if (is.na(og) || !nzchar(og)) {
                    return(character(0))
                }

                get_orthogroup_members(og, source_species)
            }),
            source_gene = purrr::map_chr(source_members, first_nonempty)
        ) %>%
        pull(source_gene) %>%
        unique() %>%
        { .[!is.na(.) & nzchar(.)] }
}

parse_cross_feature_ids <- function(feature_ids) {
    feature_ids <- as.character(feature_ids)

    infer_feature_species <- function(feature_id) {
        if (grepl("^Mtrun", feature_id)) {
            return("medicago")
        }

        if (grepl("^Glyma[.]", feature_id)) {
            return("glycine")
        }

        if (grepl("^LotjaGi", feature_id)) {
            return("lotus")
        }

        for (species_key in within_species_keys) {
            canonical_id <- canonicalize_gene_ids(species_key, feature_id)

            if (!is.na(canonical_id) &&
                orthogroup_long %>%
                    filter(species == species_key, canonical_gene_id == !!canonical_id) %>%
                    nrow() > 0) {
                return(species_key)
            }
        }

        NA_character_
    }

    tibble(feature_id = feature_ids) %>%
        mutate(
            feature_species = ifelse(
                grepl("::", feature_id, fixed = TRUE),
                sub("^([^:]+)::.*$", "\\1", feature_id),
                vapply(feature_id, infer_feature_species, character(1))
            ),
            feature_gene_id = ifelse(
                grepl("::", feature_id, fixed = TRUE),
                sub("^[^:]+::", "", feature_id),
                feature_id
            )
        )
}

display_cross_feature_labels <- function(cross_key, feature_ids) {
    feature_ids <- as.character(feature_ids)

    if (!length(feature_ids)) {
        return(character(0))
    }

    parsed_tbl <- parse_cross_feature_ids(feature_ids)

    purrr::pmap_chr(
        list(parsed_tbl$feature_species, parsed_tbl$feature_gene_id, parsed_tbl$feature_id),
        function(feature_species, feature_gene_id, feature_id) {
            if (!(feature_species %in% within_species_keys) || !nzchar(feature_gene_id)) {
                return(feature_id)
            }

            feature_label <- display_gene_labels(feature_species, feature_gene_id)[[1]]

            if (identical(cross_integration_registry[[cross_key]]$feature_mode, "medicago_space") &&
                identical(feature_species, "medicago")) {
                return(feature_label)
            }

            paste0(feature_label, " [", species_label(feature_species), "]")
        }
    )
}

map_cross_marker_features_to_source_genes <- function(cross_key, feature_ids, source_species) {
    feature_ids <- unique(as.character(feature_ids))
    feature_ids <- feature_ids[!is.na(feature_ids) & nzchar(feature_ids)]

    if (!length(feature_ids)) {
        return(character(0))
    }

    parsed_tbl <- parse_cross_feature_ids(feature_ids) %>%
        filter(feature_species %in% within_species_keys, !is.na(feature_gene_id) & nzchar(feature_gene_id))

    if (!nrow(parsed_tbl)) {
        return(character(0))
    }

    unlist(purrr::pmap(
        list(parsed_tbl$feature_species, parsed_tbl$feature_gene_id),
        function(feature_species, feature_gene_id) {
            map_target_genes_to_source_genes(feature_species, feature_gene_id, source_species)
        }
    ), use.names = FALSE) %>%
        unique() %>%
        { .[!is.na(.) & nzchar(.)] }
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
    maybe_add("Cluster", "cluster_label")
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

cross_composition_choices <- function(obj) {
    available_cols <- colnames(obj@meta.data)
    choices <- character(0)
    maybe_add <- function(label, column) {
        if (column %in% available_cols && !(column %in% unname(choices))) {
            choices <<- c(choices, setNames(column, label))
        }
    }
    maybe_add("Species", "species")
    maybe_add("Time point", "time_point")
    maybe_add("Sample", "Sample")
    maybe_add("Sample", "sample")
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

cross_distribution_group_choices <- function(obj, cross_key) {
    available_cols <- colnames(obj@meta.data)
    choices <- character(0)

    maybe_add <- function(label, column) {
        if (column %in% available_cols && !(column %in% unname(choices))) {
            choices <<- c(choices, setNames(column, label))
        }
    }

    maybe_add("Species", "species")
    maybe_add("Clustering opt 1", "Rank_1st")
    maybe_add("Clustering opt 2", "Rank_2nd")
    maybe_add("Clustering opt 3", "Rank_3rd")
    maybe_add("Clustering opt 4", "Rank_4th")
    maybe_add("Clustering opt 5", "Rank_5th")

    choices
}

cross_distribution_split_choices <- function(obj) {
    available_cols <- colnames(obj@meta.data)
    choices <- c("No split" = "none")

    if ("species" %in% available_cols) {
        choices <- c(choices, "Species" = "species")
    }

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
    maybe_add("Clustering opt 2 label", "Rank_2nd_label")
    maybe_add("Clustering opt 3 label", "Rank_3rd_label")
    maybe_add("Clustering opt 4 label", "Rank_4th_label")
    maybe_add("Clustering opt 5 label", "Rank_5th_label")
    maybe_add("Condition", "condition")
    maybe_add("Study", "study")
    maybe_add("Time point", "time_point")

    choices
}

gene_panel_selectize_options <- function(enable_bulk = FALSE) {
    options <- list(
        plugins = list("remove_button"),
        valueField = "value",
        labelField = "label",
        searchField = c("label", "value", "tokens"),
        searchConjunction = "or",
        maxOptions = 50,
        closeAfterSelect = FALSE,
        hideSelected = TRUE,
        placeholder = "Type a gene ID, common name, or synonym"
    )

    if (isTRUE(enable_bulk)) {
        options$onType <- I("function(query) { window.atlasScheduleBulkGeneOption(this, query); }")
        options$onDropdownOpen <- I("function() { window.atlasScheduleBulkGeneOption(this, this.lastValue || this.$control_input.val()); }")
        options$onItemAdd <- I("function(value) { window.atlasHandleBulkGeneOption(this, value); }")
    }

    options
}
