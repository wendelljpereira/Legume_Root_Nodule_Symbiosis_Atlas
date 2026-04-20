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
within_three_d_reduction_name <- "umap3d"
cross_feature_lookup_path <- "metadata/cross_feature_lookup.tsv"

# User-overridable paths. See README for Docker usage; set these env vars to
# point at a locally-maintained orthogroups table or celltype override folder
# without editing the source.
resolve_env_path <- function(env_var, default) {
    value <- Sys.getenv(env_var, unset = "")
    if (!nzchar(value)) default else value
}

orthogroups_path <- resolve_env_path("ATLAS_ORTHOLOGS_PATH", "orthogroups/joint_orthogroups.tsv")
celltype_overrides_dir <- resolve_env_path("ATLAS_CELLTYPE_OVERRIDES_DIR", "celltype_overrides")

atlas_version <- resolve_env_path("ATLAS_VERSION", "1.0")
atlas_last_updated <- resolve_env_path("ATLAS_LAST_UPDATED", format(Sys.Date(), "%Y-%m-%d"))
atlas_citation_text <- resolve_env_path(
    "ATLAS_CITATION",
    "Pereira W. et al. A cross-species single-cell atlas of legume root nodule symbiosis. (in preparation, 2026). Please cite prior to publication as: 'Legume Root Nodule Symbiosis Atlas, pre-publication release.'"
)

# Pre-publication access gate. Set ATLAS_ACCESS_PASSWORD on the server
# (ShinyApps.io > Advanced > Environment variables) to require a password
# before the app renders. Leave unset during local Docker use.
atlas_access_password <- Sys.getenv("ATLAS_ACCESS_PASSWORD", unset = "")
atlas_access_required <- nzchar(atlas_access_password)

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
        description = "This tab uses the SATURN integration object. Selected source genes resolve to ortholog features from all three species in the integrated SATURN feature space.",
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

expression_export_row_limit <- 50000L

build_expression_export_feature_map <- function(resolution, source_species, source_genes) {
    source_genes <- unique(as.character(source_genes))

    if (!length(source_genes)) {
        return(tibble(
            source_gene = character(0),
            gene_label = character(0),
            column_name = character(0),
            feature_ids = list()
        ))
    }

    base_tbl <- tibble(
        source_gene = source_genes,
        gene_label = display_gene_labels(
            source_species,
            source_genes,
            include_gene_id_with_common = FALSE
        )
    ) %>%
        mutate(
            gene_label = ifelse(is.na(gene_label) | !nzchar(gene_label), source_gene, gene_label)
        )

    mapping_tbl <- resolution$mapping %||% tibble()
    feature_col <- if (nrow(mapping_tbl)) {
        pick_first_existing_col(mapping_tbl, c("target_feature_id", "feature_id"))
    } else {
        NA_character_
    }

    if (is.na(feature_col) || !nrow(mapping_tbl)) {
        base_tbl$feature_ids <- replicate(nrow(base_tbl), character(0), simplify = FALSE)
    } else {
        mapped_tbl <- mapping_tbl %>%
            mutate(feature_id = as.character(.data[[feature_col]])) %>%
            filter(!is.na(feature_id) & nzchar(feature_id)) %>%
            group_by(source_gene) %>%
            summarise(feature_ids = list(unique(feature_id)), .groups = "drop")

        base_tbl <- left_join(base_tbl, mapped_tbl, by = "source_gene")

        if (!("feature_ids" %in% colnames(base_tbl))) {
            base_tbl$feature_ids <- vector("list", nrow(base_tbl))
        }

        for (idx in seq_len(nrow(base_tbl))) {
            if (is.null(base_tbl$feature_ids[[idx]])) {
                base_tbl$feature_ids[[idx]] <- character(0)
            }
        }
    }

    base_tbl %>%
        mutate(column_name = make.unique(gene_label, sep = " | "))
}

collect_expression_export_values <- function(obj, feature_map) {
    if (is.null(feature_map) || !nrow(feature_map)) {
        return(tibble())
    }

    all_features <- unique(unlist(feature_map$feature_ids))
    data_matrix <- get_cached_data_matrix(obj)
    all_features <- intersect(as.character(all_features), rownames(data_matrix))
    expr_data <- if (length(all_features)) {
        data_matrix[all_features, , drop = FALSE]
    } else {
        NULL
    }

    expr_cols <- lapply(seq_len(nrow(feature_map)), function(idx) {
        feature_ids <- intersect(feature_map$feature_ids[[idx]], rownames(data_matrix))

        if (!length(feature_ids)) {
            return(rep(NA_real_, ncol(obj)))
        }

        if (length(feature_ids) == 1L) {
            return(as.numeric(expr_data[feature_ids[[1]], ]))
        }

        as.numeric(Matrix::colMeans(expr_data[feature_ids, , drop = FALSE]))
    })

    tibble::as_tibble(stats::setNames(expr_cols, feature_map$column_name))
}

expression_export_metadata <- function(obj, default_species = NA_character_) {
    md <- obj@meta.data
    cluster_ids <- as.character(Idents(obj))
    cluster_lookup <- cluster_label_lookup(obj)
    cluster_labels <- cluster_lookup$cluster_label[match(cluster_ids, cluster_lookup$cluster)]
    cluster_labels[is.na(cluster_labels) | !nzchar(cluster_labels)] <- cluster_ids[is.na(cluster_labels) | !nzchar(cluster_labels)]

    meta_tbl <- tibble(
        cell_id = rownames(md),
        cluster_id = cluster_ids,
        cluster_label = cluster_labels
    )

    species_col <- pick_first_existing_col(md, c("species"))
    sample_col <- pick_first_existing_col(md, c("Sample", "sample_name", "sample", "orig.ident"))
    condition_col <- pick_first_existing_col(md, c("Group", "group", "condition", "treatment"))
    time_col <- pick_first_existing_col(md, c("time", "time_point"))
    study_col <- pick_first_existing_col(md, c("study"))

    if (!is.na(species_col)) {
        meta_tbl$species <- as.character(md[[species_col]])
    } else if (!is.na(default_species) && nzchar(default_species)) {
        meta_tbl$species <- default_species
    }

    if (!is.na(sample_col)) {
        meta_tbl$sample <- as.character(md[[sample_col]])
    }

    if (!is.na(condition_col)) {
        meta_tbl$condition <- as.character(md[[condition_col]])
    }

    if (!is.na(time_col)) {
        meta_tbl$time <- as.character(md[[time_col]])
    }

    if (!is.na(study_col)) {
        meta_tbl$study <- as.character(md[[study_col]])
    }

    meta_tbl
}

summarise_expression_export_values <- function(values) {
    values <- unique(as.character(values))
    values <- values[!is.na(values) & nzchar(values)]

    if (!length(values)) {
        return(NA_character_)
    }

    paste(sort(values), collapse = "; ")
}

resolve_expression_export_mode <- function(requested_mode, row_count, row_limit = expression_export_row_limit) {
    requested_mode <- requested_mode %||% "per_cell"

    if (identical(requested_mode, "per_cell") && isTRUE(row_count > row_limit)) {
        return("per_cluster")
    }

    requested_mode
}

build_expression_export_table <- function(obj, feature_map, mode = "per_cell", default_species = NA_character_) {
    meta_tbl <- expression_export_metadata(obj, default_species = default_species)
    expr_tbl <- collect_expression_export_values(obj, feature_map)
    cell_tbl <- bind_cols(meta_tbl, expr_tbl)

    if (!identical(mode, "per_cluster")) {
        return(cell_tbl)
    }

    if (!nrow(cell_tbl)) {
        return(cell_tbl)
    }

    cluster_levels <- cluster_value_levels(cell_tbl$cluster_id)
    meta_cols <- intersect(c("species", "sample", "condition", "time", "study"), colnames(cell_tbl))
    expr_cols <- feature_map$column_name
    cluster_indices <- split(
        seq_len(nrow(cell_tbl)),
        factor(cell_tbl$cluster_id, levels = cluster_levels)
    )
    cluster_indices <- cluster_indices[lengths(cluster_indices) > 0]

    cluster_meta_tbl <- cell_tbl %>%
        group_by(cluster_id, cluster_label) %>%
        summarise(
            across(all_of(meta_cols), summarise_expression_export_values),
            .groups = "drop"
        ) %>%
        mutate(cluster_id = factor(cluster_id, levels = cluster_levels)) %>%
        arrange(cluster_id) %>%
        mutate(cluster_id = as.character(cluster_id))

    cluster_expr_tbl <- tibble(cluster_id = names(cluster_indices))

    for (col_name in expr_cols) {
        values <- cell_tbl[[col_name]]
        cluster_expr_tbl[[col_name]] <- vapply(cluster_indices, function(idx) {
            cluster_values <- values[idx]

            if (!length(cluster_values) || all(is.na(cluster_values))) {
                return(NA_real_)
            }

            mean(cluster_values, na.rm = TRUE)
        }, numeric(1))
    }

    cluster_expr_tbl <- cluster_expr_tbl %>%
        mutate(cluster_id = factor(cluster_id, levels = cluster_levels)) %>%
        arrange(cluster_id) %>%
        mutate(cluster_id = as.character(cluster_id))

    cluster_meta_tbl %>%
        left_join(cluster_expr_tbl, by = "cluster_id") %>%
        select(any_of(c("cluster_id", "cluster_label", meta_cols, expr_cols)))
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

label_like_metadata_columns <- c(
    "cell_class",
    "cluster_label",
    "species_cell_class",
    "saturn_label",
    "saturn_ref_label",
    "Rank_1st_label"
)

distribution_cluster_columns <- c("Rank_1st", "Rank_2nd", "Rank_3rd", "Rank_4th", "Rank_5th", "cluster_label")

is_cluster_distribution_group <- function(column_name) {
    length(column_name) == 1L && !is.na(column_name) && column_name %in% distribution_cluster_columns
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

expression_heatmap_palette <- function(colorblind_safe = FALSE) {
    if (isTRUE(colorblind_safe)) {
        return(viridisLite::cividis(100, end = 0.95))
    }

    c(
        "#f6f8f4",
        "#d8e2d6",
        unname(app_palette["accent"]),
        unname(app_palette["accent_dark"])
    )
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

        if (exists(cache_key, envir = entry$expression_cache, inherits = FALSE)) {
            return(get(cache_key, envir = entry$expression_cache, inherits = FALSE))
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
        assign(cache_key, sampled_ids, envir = entry$expression_cache)
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

    if (!exists(cache_key, envir = entry$expression_cache, inherits = FALSE)) {
        assign(
            cache_key,
            list(
                object = obj[, sampled_cell_ids],
                cell_ids = sampled_cell_ids
            ),
            envir = entry$expression_cache
        )
        set_dataset_cache_entry(dataset_key, entry)
    }

    get(cache_key, envir = entry$expression_cache, inherits = FALSE)
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
    colorblind_safe = FALSE,
    max_cells = 2000L
) {
    sampled_selection <- get_cached_sampled_object(
        obj,
        group_by = group_by,
        max_cells = max_cells,
        seed = 123L
    )
    sampled_obj <- sampled_selection$object
    sampled_cell_ids <- sampled_selection$cell_ids

    if (!length(sampled_cell_ids)) {
        stop("No cells are available for the current heatmap selection.")
    }

    feature_ids <- intersect(unique(as.character(feature_ids)), rownames(sampled_obj))

    if (!length(feature_ids)) {
        stop("No mapped genes are available for the current heatmap selection.")
    }

    feature_labels <- expression_feature_labels(feature_ids, label_map)
    group_values <- sampled_obj@meta.data[[group_by]]
    group_colors <- group_palette_for_toggle(
        values = group_values,
        column_name = group_by,
        colorblind_safe = colorblind_safe
    )

    plot_obj <- Seurat::DoHeatmap(
        object = sampled_obj,
        features = feature_ids,
        cells = sampled_cell_ids,
        group.by = group_by,
        group.bar = TRUE,
        group.colors = unname(group_colors),
        disp.min = 0,
        disp.max = NULL,
        slot = "data",
        label = TRUE,
        size = 4.1,
        angle = 0,
        raster = TRUE,
        draw.lines = FALSE
    )

    plot_obj$scales$scales <- Filter(
        function(scale_obj) !("fill" %in% scale_obj$aesthetics),
        plot_obj$scales$scales
    )

    plot_obj +
        scale_fill_gradientn(
            colors = expression_heatmap_palette(colorblind_safe),
            name = "Expression"
        ) +
        scale_y_discrete(labels = feature_labels) +
        labs(
            x = metadata_column_label(group_by),
            y = NULL
        ) +
        app_plot_theme() +
        theme(
            panel.grid = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = margin(8, 14, 12, 10)
        )
}

build_expression_ridge_plot <- function(
    obj,
    feature_ids,
    label_map,
    group_by,
    colorblind_safe = FALSE,
    max_cells = 3000L
) {
    sampled_selection <- get_cached_sampled_object(
        obj,
        group_by = group_by,
        max_cells = max_cells,
        seed = 456L
    )
    sampled_obj <- sampled_selection$object
    sampled_cell_ids <- sampled_selection$cell_ids

    if (!length(sampled_cell_ids)) {
        stop("No cells are available for the current ridge plot selection.")
    }

    feature_ids <- intersect(unique(as.character(feature_ids)), rownames(sampled_obj))

    if (!length(feature_ids)) {
        stop("No mapped genes are available for the current ridge plot selection.")
    }

    feature_labels <- expression_feature_labels(feature_ids, label_map)
    group_values <- sampled_obj@meta.data[[group_by]]
    group_colors <- group_palette_for_toggle(
        values = group_values,
        column_name = group_by,
        colorblind_safe = colorblind_safe
    )

    ridge_plots <- Seurat::RidgePlot(
        object = sampled_obj,
        features = feature_ids,
        group.by = group_by,
        cols = unname(group_colors),
        same.y.lims = TRUE,
        combine = FALSE,
        fill.by = "ident",
        layer = "data"
    )

    ridge_plots <- purrr::imap(ridge_plots, function(plot_obj, idx) {
        feature_id <- feature_ids[[idx]]

        plot_obj +
            labs(
                title = unname(feature_labels[[feature_id]]),
                x = "Normalized expression",
                y = metadata_column_label(group_by)
            ) +
            app_plot_theme() +
            theme(
                plot.title = element_text(
                    face = "bold",
                    colour = app_palette["text"],
                    size = 16,
                    hjust = 0
                ),
                legend.position = "none",
                panel.grid.major.y = element_blank(),
                plot.margin = margin(8, 14, 12, 10)
            )
    })

    if (length(ridge_plots) == 1L) {
        ridge_plots[[1]]
    } else {
        wrap_plots(plotlist = ridge_plots, ncol = 1)
    }
}

expression_ridge_height_px <- function(feature_n, group_n) {
    max(
        360L,
        as.integer(170L + feature_n * 210L + group_n * 10L)
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
        420L,
        as.integer(160L + max(1L, gene_n) * 270L)
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
    if (!nzchar(celltype_overrides_dir) || !dir.exists(celltype_overrides_dir)) {
        return(NULL)
    }

    candidate_paths <- file.path(
        celltype_overrides_dir,
        c(paste0(override_id, ".csv"), paste0(override_id, ".tsv"))
    )
    path <- candidate_paths[file.exists(candidate_paths)][1]
    if (is.na(path) || !length(path)) {
        return(NULL)
    }

    sep <- if (grepl("\\.tsv$", path, ignore.case = TRUE)) "\t" else ","
    df <- tryCatch(
        utils::read.delim(path, sep = sep, stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e) {
            warning(sprintf("Failed to read celltype override '%s': %s", path, conditionMessage(e)), call. = FALSE)
            NULL
        }
    )
    if (is.null(df) || !nrow(df)) return(NULL)

    cols <- tolower(colnames(df))
    id_col <- which(cols %in% c("cluster_id", "cluster", "id"))[1]
    label_col <- which(cols %in% c("label", "celltype", "cell_type", "cluster_label"))[1]
    if (is.na(id_col) || is.na(label_col)) {
        warning(sprintf(
            "Celltype override '%s' must have columns 'cluster_id' and 'label' (received: %s). Skipping.",
            path, paste(colnames(df), collapse = ", ")
        ), call. = FALSE)
        return(NULL)
    }

    data.frame(
        cluster_id = as.character(df[[id_col]]),
        label = as.character(df[[label_col]]),
        stringsAsFactors = FALSE
    )
}

apply_celltype_override <- function(obj, override_id) {
    mapping <- read_celltype_override(override_id)
    obj$cluster_label <- as.character(Idents(obj))
    if (is.null(mapping) || !nrow(mapping)) {
        return(obj)
    }

    lookup <- setNames(mapping$label, mapping$cluster_id)
    current <- obj$cluster_label
    replaced <- lookup[current]
    obj$cluster_label <- ifelse(is.na(replaced), current, unname(replaced))
    obj$celltype <- obj$cluster_label
    obj
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

cache_get <- function(key, builder) {
    if (!exists(key, envir = dataset_cache, inherits = FALSE)) {
        assign(key, builder(), envir = dataset_cache)
    }

    get(key, envir = dataset_cache, inherits = FALSE)
}

get_dataset_cache_entry <- function(key, builder) {
    if (!exists(key, envir = dataset_cache, inherits = FALSE)) {
        assign(
            key,
            list(
                object = builder(),
                data_matrix = NULL,
                expression_cache = new.env(parent = emptyenv())
            ),
            envir = dataset_cache
        )
    }

    entry <- get(key, envir = dataset_cache, inherits = FALSE)

    if (is.null(entry$expression_cache) || !is.environment(entry$expression_cache)) {
        entry$expression_cache <- new.env(parent = emptyenv())
        assign(key, entry, envir = dataset_cache)
    }

    entry
}

set_dataset_cache_entry <- function(key, entry) {
    assign(key, entry, envir = dataset_cache)
    invisible(entry)
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

    if (!exists(cache_key, envir = entry$expression_cache, inherits = FALSE)) {
        assign(
            cache_key,
            compute_feature_expression_values(data_matrix, feature_ids, cell_ids),
            envir = entry$expression_cache
        )
        set_dataset_cache_entry(dataset_key, entry)
    }

    get(cache_key, envir = entry$expression_cache, inherits = FALSE)
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

normalize_cluster_markers_cache <- function(marker_tbl) {
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

    cluster_col <- pick_marker_col(c("cluster", "cluster_id", "ident", "identity"))
    gene_col <- pick_marker_col(c("gene", "feature", "feature_id"))
    cluster_source_col <- pick_marker_col(c(
        "cluster_source",
        "clustering_level",
        "clustering_column",
        "cluster_level",
        "cluster_scheme",
        "group_by",
        "cluster_by"
    ))
    logfc_col <- pick_marker_col(c("avg_log2FC", "avg_logFC", "avg_log10FC"))
    pct1_col <- pick_marker_col(c("pct.1", "pct1", "pct_1"))
    pct2_col <- pick_marker_col(c("pct.2", "pct2", "pct_2"))
    pval_adj_col <- pick_marker_col(c("p_val_adj", "padj", "p_adj"))

    required_cols <- c(cluster_col, gene_col, logfc_col, pct1_col, pct2_col, pval_adj_col)

    if (any(is.na(required_cols))) {
        return(NULL)
    }

    cluster_source_values <- if (!is.na(cluster_source_col)) {
        as.character(marker_tbl[[cluster_source_col]])
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
        cache_path <- pick_first_existing_path(cluster_markers_cache_candidates(dataset_key, top_n = top_n))

        if (!length(cache_path) || is.na(cache_path) || !nzchar(cache_path)) {
            return(NULL)
        }

        cached_markers <- read_delimited_cache(cache_path)
        normalize_cluster_markers_cache(cached_markers)
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
        "sample_n",
        "group_label", "group_n", "group_preview",
        "species_n", "species_preview",
        "time_n", "time_preview"
    )

    if (!nrow(summary_df) || !all(required_cols %in% colnames(summary_df))) {
        return(FALSE)
    }

    valid_preview <- function(n_col, preview_col) {
        counts <- suppressWarnings(as.numeric(summary_df[[n_col]]))
        previews <- trimws(as.character(summary_df[[preview_col]] %||% ""))
        length(counts) == nrow(summary_df) &&
            length(previews) == nrow(summary_df) &&
            all((!is.na(counts) & counts > 0) | !nzchar(previews)) &&
            any(!is.na(counts) & counts > 0 & nzchar(previews))
    }

    valid_preview("group_n", "group_preview") &&
        valid_preview("species_n", "species_preview") &&
        valid_preview("time_n", "time_preview")
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

cross_composition_choices <- function(obj) {
    available_cols <- colnames(obj@meta.data)
    choices <- character(0)
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

cross_distribution_group_choices <- function(obj, cross_key) {
    available_cols <- colnames(obj@meta.data)
    choices <- character(0)

    maybe_add <- function(label, column) {
        if (column %in% available_cols && !(column %in% unname(choices))) {
            choices <<- c(choices, setNames(column, label))
        }
    }

    maybe_add("Time point", "time_point")
    maybe_add("Species", "species")
    maybe_add("Clustering opt 1", "Rank_1st")
    maybe_add("Clustering opt 2", "Rank_2nd")
    maybe_add("Clustering opt 3", "Rank_3rd")
    maybe_add("Clustering opt 4", "Rank_4th")
    maybe_add("Clustering opt 5", "Rank_5th")

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
                h2(tagList(species_label_tag(species_key), " explorer")),
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
                class = "subsection-header permalink-panel",
                `data-permalink-panel` = paste0(species_key, "_distribution"),
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
                class = "subsection-header permalink-panel",
                `data-permalink-panel` = paste0(species_key, "_composition"),
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
                class = "subsection-header permalink-panel",
                `data-permalink-panel` = paste0(species_key, "_markers"),
                h3("Cluster markers"),
                p("Review precomputed positive markers for the active clustering solution when marker tables are available for it, add the top hits to the shared gene panel, or download the current cluster as CSV.")
            ),
            fluidRow(
                column(
                    width = 4,
                    div(
                        class = "option-group",
                        uiOutput(paste0(species_key, "_marker_cluster_ui"))
                    )
                ),
                column(
                    width = 2,
                    div(
                        class = "option-group",
                        numericInput(
                            inputId = paste0(species_key, "_marker_top_n"),
                            label = "Top markers to add",
                            value = 10,
                            min = 1,
                            max = 25,
                            step = 1
                        )
                    )
                ),
                column(
                    width = 6,
                    div(
                        class = "option-group marker-action-group",
                        div(
                            class = "marker-action-row",
                            actionButton(
                                inputId = paste0(species_key, "_add_markers"),
                                label = "Add top N to gene panel",
                                icon = icon("plus"),
                                class = "btn btn-default btn-sm plot-download-btn"
                            ),
                            downloadButton(
                                outputId = paste0("dl_", species_key, "_markers"),
                                label = "Download markers for this cluster (CSV)",
                                class = "btn btn-default btn-sm plot-download-btn"
                            )
                        ),
                        uiOutput(paste0(species_key, "_markers_status_ui"))
                    )
                )
            ),
            fluidRow(
                column(
                    width = 12,
                    div(
                        class = "table-card permalink-panel",
                        `data-permalink-panel` = paste0(species_key, "_markers"),
                        div(class = "plot-card-title", "Positive markers for the selected cluster"),
                        DT::DTOutput(paste0(species_key, "_markers_table"))
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
                    width = 8,
                    div(
                        class = "option-group matrix-export-group",
                        radioButtons(
                            inputId = paste0(species_key, "_matrix_mode"),
                            label = "Expression matrix export",
                            choices = c(
                                "Per cell (long)" = "per_cell",
                                "Per cluster mean (wide)" = "per_cluster"
                            ),
                            selected = "per_cell",
                            inline = TRUE
                        ),
                        uiOutput(paste0(species_key, "_matrix_hint_ui"))
                    )
                ),
                column(
                    width = 4,
                    div(
                        class = "option-group matrix-export-download-group",
                        downloadButton(
                            outputId = paste0("dl_", species_key, "_expression_matrix"),
                            label = "Download expression matrix (CSV)",
                            class = "btn btn-default btn-sm plot-download-btn"
                        )
                    )
                )
            ),
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
                        class = "plot-card permalink-panel",
                        `data-permalink-panel` = paste0(species_key, "_umap"),
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
                        class = "plot-card permalink-panel",
                        `data-permalink-panel` = paste0(species_key, "_violin"),
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
                        class = "plot-card permalink-panel",
                        `data-permalink-panel` = paste0(species_key, "_heatmap"),
                        div(
                            class = "plot-card-header",
                            div(class = "plot-card-title", "Expression heatmap"),
                            plot_download_button(paste0("dl_", species_key, "_heatmap"))
                        ),
                        spinning_plot_output(paste0(species_key, "_heatmap_plot"), proxy_height = "520px")
                    )
                )
            ),
            fluidRow(
                column(
                    width = 12,
                    div(
                        class = "plot-card permalink-panel",
                        `data-permalink-panel` = paste0(species_key, "_ridge"),
                        div(
                            class = "plot-card-header",
                            div(class = "plot-card-title", "Expression ridge plots"),
                            plot_download_button(paste0("dl_", species_key, "_ridge"))
                        ),
                        spinning_plot_output(paste0(species_key, "_ridge_plot"), proxy_height = "520px")
                    )
                )
            ),
            fluidRow(
                column(
                    width = 12,
                    div(
                        class = "plot-card permalink-panel",
                        `data-permalink-panel` = paste0(species_key, "_dot"),
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
            # ── Cell distribution UMAP ──────────────────────────────────────
            div(
                class = "subsection-header permalink-panel",
                `data-permalink-panel` = paste0(prefix, "_distribution"),
                h3("Cell distribution UMAP"),
                p("Inspect cluster structure and species mixing in the integrated embedding.")
            ),
            fluidRow(
                column(
                    width = 4,
                    div(
                        class = "option-group",
                        uiOutput(paste0(prefix, "_dist_group_by_ui"))
                    )
                ),
                column(
                    width = 4,
                    div(
                        class = "option-group",
                        sliderInput(
                            inputId = paste0(prefix, "_dist_pt_size"),
                            label = "Distribution point size",
                            min = 0.1, max = 2.5, value = 0.75, step = 0.05
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
                            plot_download_button(paste0("dl_", prefix, "_dist_umap"))
                        ),
                        uiOutput(paste0(prefix, "_dist_umap_plot_ui"))
                    )
                )
            ),
            # ── Cluster composition ─────────────────────────────────────────
            div(
                class = "subsection-header permalink-panel",
                `data-permalink-panel` = paste0(prefix, "_composition"),
                h3("Cluster composition"),
                p("See what fraction of cells in each cluster comes from each species or condition.")
            ),
            fluidRow(
                column(
                    width = 4,
                    div(
                        class = "option-group",
                        uiOutput(paste0(prefix, "_composition_by_ui"))
                    )
                )
            ),
            fluidRow(
                column(
                    width = 12,
                    div(
                        class = "plot-card",
                        div(class = "plot-card-title", "Cells per cluster"),
                        uiOutput(paste0(prefix, "_composition_plot_ui"))
                    )
                )
            ),
            div(
                class = "subsection-header permalink-panel",
                `data-permalink-panel` = paste0(prefix, "_markers"),
                h3("Cluster markers"),
                p("Review precomputed positive markers for the active clustering solution when marker tables are available for it, add the top hits to the shared source-gene panel, or download the current cluster as CSV.")
            ),
            fluidRow(
                column(
                    width = 4,
                    div(
                        class = "option-group",
                        uiOutput(paste0(prefix, "_marker_cluster_ui"))
                    )
                ),
                column(
                    width = 2,
                    div(
                        class = "option-group",
                        numericInput(
                            inputId = paste0(prefix, "_marker_top_n"),
                            label = "Top markers to add",
                            value = 10,
                            min = 1,
                            max = 25,
                            step = 1
                        )
                    )
                ),
                column(
                    width = 6,
                    div(
                        class = "option-group marker-action-group",
                        div(
                            class = "marker-action-row",
                            actionButton(
                                inputId = paste0(prefix, "_add_markers"),
                                label = "Add top N to gene panel",
                                icon = icon("plus"),
                                class = "btn btn-default btn-sm plot-download-btn"
                            ),
                            downloadButton(
                                outputId = paste0("dl_", prefix, "_markers"),
                                label = "Download markers for this cluster (CSV)",
                                class = "btn btn-default btn-sm plot-download-btn"
                            )
                        ),
                        uiOutput(paste0(prefix, "_markers_status_ui"))
                    )
                )
            ),
            fluidRow(
                column(
                    width = 12,
                    div(
                        class = "table-card permalink-panel",
                        `data-permalink-panel` = paste0(prefix, "_markers"),
                        div(class = "plot-card-title", "Positive markers for the selected cluster"),
                        DT::DTOutput(paste0(prefix, "_markers_table"))
                    )
                )
            ),
            # ── Gene expression ─────────────────────────────────────────────
            div(
                class = "subsection-header",
                h3("Gene expression"),
                p("Compare each selected gene across the shared embedding with a species overview plus aligned Medicago, Glycine, and Lotus expression panels.")
            ),
            fluidRow(
                column(
                    width = 4,
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
                    width = 4,
                    div(
                        class = "option-group",
                        selectInput(
                            inputId = paste0(prefix, "_umap_columns"),
                            label = "Comparison panels per row",
                            choices = c("1" = 1, "2" = 2),
                            selected = 1
                        )
                    )
                ),
                column(
                    width = 4,
                    div(
                        class = "option-group",
                        uiOutput(paste0(prefix, "_group_by_ui"))
                    )
                )
            ),
            uiOutput(paste0(prefix, "_notice_ui")),
            uiOutput(paste0(prefix, "_ortholog_trace_notice_ui")),
            fluidRow(
                column(
                    width = 8,
                    div(
                        class = "option-group matrix-export-group",
                        radioButtons(
                            inputId = paste0(prefix, "_matrix_mode"),
                            label = "Expression matrix export",
                            choices = c(
                                "Per cell (long)" = "per_cell",
                                "Per cluster mean (wide)" = "per_cluster"
                            ),
                            selected = "per_cell",
                            inline = TRUE
                        ),
                        uiOutput(paste0(prefix, "_matrix_hint_ui"))
                    )
                ),
                column(
                    width = 4,
                    div(
                        class = "option-group matrix-export-download-group",
                        downloadButton(
                            outputId = paste0("dl_", prefix, "_expression_matrix"),
                            label = "Download expression matrix (CSV)",
                            class = "btn btn-default btn-sm plot-download-btn"
                        )
                    )
                )
            ),
            fluidRow(
                column(
                    width = 12,
                    div(
                        class = "plot-card permalink-panel",
                        `data-permalink-panel` = paste0(prefix, "_ortholog_trace"),
                        div(
                            class = "plot-card-header",
                            div(class = "plot-card-title", "Ortholog trace"),
                            plot_download_button(paste0("dl_", prefix, "_ortholog_trace"))
                        ),
                        spinning_plot_output(paste0(prefix, "_ortholog_trace_plot"), proxy_height = "760px")
                    )
                )
            ),
            fluidRow(
                column(
                    width = 12,
                    div(
                        class = "plot-card permalink-panel",
                        `data-permalink-panel` = paste0(prefix, "_umap"),
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
                    width = 12,
                    div(
                        class = "plot-card permalink-panel",
                        `data-permalink-panel` = paste0(prefix, "_heatmap"),
                        div(
                            class = "plot-card-header",
                            div(class = "plot-card-title", "Expression heatmap"),
                            plot_download_button(paste0("dl_", prefix, "_heatmap"))
                        ),
                        spinning_plot_output(paste0(prefix, "_heatmap_plot"), proxy_height = "520px")
                    )
                )
            ),
            fluidRow(
                column(
                    width = 12,
                    div(
                        class = "plot-card permalink-panel",
                        `data-permalink-panel` = paste0(prefix, "_ridge"),
                        div(
                            class = "plot-card-header",
                            div(class = "plot-card-title", "Expression ridge plots"),
                            plot_download_button(paste0("dl_", prefix, "_ridge"))
                        ),
                        spinning_plot_output(paste0(prefix, "_ridge_plot"), proxy_height = "520px")
                    )
                )
            ),
            fluidRow(
                column(
                    width = 7,
                    div(
                        class = "plot-card permalink-panel",
                        `data-permalink-panel` = paste0(prefix, "_dot"),
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
            href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=IBM+Plex+Sans:wght@500;600;700&family=IBM+Plex+Mono:wght@400;500&display=swap"
        ),
        tags$link(rel = "icon", type = "image/svg+xml", href = "favicon.svg"),
        tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
        tags$script(HTML(
            "(function() {
                function atlasVisiblePermalinkPanel() {
                    var anchors = Array.prototype.slice.call(document.querySelectorAll('.permalink-panel'));
                    anchors = anchors.filter(function(el) {
                        return el && el.offsetParent !== null;
                    });
                    if (!anchors.length) return null;

                    var targetOffset = 132;
                    var best = null;
                    var bestScore = Infinity;

                    anchors.forEach(function(el) {
                        var rect = el.getBoundingClientRect();
                        var visible = rect.bottom > targetOffset && rect.top < window.innerHeight * 0.8;
                        var score = Math.abs(rect.top - targetOffset);
                        if (!visible) score += 10000;
                        if (score < bestScore) {
                            bestScore = score;
                            best = el;
                        }
                    });

                    return best ? best.getAttribute('data-permalink-panel') : null;
                }

                function atlasFindPermalinkPanel(panelId) {
                    var anchors = Array.prototype.slice.call(document.querySelectorAll('.permalink-panel'));
                    for (var i = 0; i < anchors.length; i += 1) {
                        if (anchors[i].getAttribute('data-permalink-panel') === panelId) {
                            return anchors[i];
                        }
                    }
                    return null;
                }

                function atlasShowToast(message, tone) {
                    var toast = document.getElementById('atlas-toast');
                    if (!toast) {
                        toast = document.createElement('div');
                        toast.id = 'atlas-toast';
                        toast.className = 'atlas-toast';
                        document.body.appendChild(toast);
                    }

                    toast.textContent = message;
                    toast.className = 'atlas-toast is-visible ' + (tone || 'success');

                    if (toast._hideTimer) {
                        window.clearTimeout(toast._hideTimer);
                    }

                    toast._hideTimer = window.setTimeout(function() {
                        toast.className = 'atlas-toast';
                    }, 2200);
                }

                function atlasAppBusy() {
                    return !!document.querySelector('.shiny-busy, .recalculating');
                }

                window.atlasBulkGeneOptionValue = function(query) {
                    return '__bulk__::' + encodeURIComponent(query || '');
                };

                window.atlasUpdateBulkGeneOption = function(selectize, query) {
                    if (!selectize) return;

                    Object.keys(selectize.options || {}).forEach(function(key) {
                        if (key.indexOf('__bulk__::') === 0) {
                            selectize.removeOption(key, true);
                        }
                    });

                    query = (query || '').trim();
                    if (!query.length) {
                        selectize.refreshOptions(false);
                        return;
                    }

                    var optionValue = window.atlasBulkGeneOptionValue(query);
                    selectize.addOption({
                        value: optionValue,
                        label: 'Select all matches for \\\"' + query + '\\\"',
                        tokens: query,
                        bulk_select: true
                    });
                    selectize.refreshOptions(false);

                    window.requestAnimationFrame(function() {
                        var optionNode = selectize.getOption(optionValue);
                        if (optionNode && optionNode.length) {
                            optionNode.addClass('selectize-bulk-option');
                            var parent = optionNode.parent();
                            if (parent && parent.length) {
                                parent.prepend(optionNode);
                            }
                        }
                    });
                };

                window.atlasScheduleBulkGeneOption = function(selectize, query) {
                    window.atlasUpdateBulkGeneOption(selectize, query);
                    [90, 220].forEach(function(delayMs) {
                        window.setTimeout(function() {
                            window.atlasUpdateBulkGeneOption(selectize, query);
                        }, delayMs);
                    });
                };

                window.atlasHandleBulkGeneOption = function(selectize, value) {
                    if (!value || value.indexOf('__bulk__::') !== 0) {
                        return false;
                    }

                    var query = '';
                    try {
                        query = decodeURIComponent(value.replace('__bulk__::', ''));
                    } catch (e) {
                        query = value.replace('__bulk__::', '');
                    }

                    window.setTimeout(function() {
                        selectize.removeItem(value, true);
                        selectize.setTextboxValue('');
                        selectize.close();
                        window.atlasUpdateBulkGeneOption(selectize, '');
                    }, 0);

                    if (window.Shiny && Shiny.setInputValue) {
                        Shiny.setInputValue('selected_genes_bulk_query', {
                            query: query,
                            nonce: Date.now()
                        }, { priority: 'event' });
                    }

                    return true;
                };

                function atlasSetButtonBusy(buttonId, busy, busyLabel) {
                    var button = document.getElementById(buttonId);
                    if (!button) return;

                    if (!button.dataset.atlasOriginalHtml) {
                        button.dataset.atlasOriginalHtml = button.innerHTML;
                    }

                    if (busy) {
                        button.disabled = true;
                        button.classList.add('is-working');
                        button.setAttribute('aria-busy', 'true');
                        button.innerHTML = '<span class=\"atlas-inline-spinner\" aria-hidden=\"true\"></span><span>' + (busyLabel || 'Working...') + '</span>';
                        return;
                    }

                    button.disabled = false;
                    button.classList.remove('is-working');
                    button.removeAttribute('aria-busy');

                    if (button.dataset.atlasOriginalHtml) {
                        button.innerHTML = button.dataset.atlasOriginalHtml;
                    }
                }

                function atlasReportVisiblePermalinkPanel(panelId) {
                    if (!(window.Shiny && Shiny.setInputValue)) {
                        return;
                    }

                    if (!panelId && !window.__atlasRestoringPanelId && atlasAppBusy()) {
                        return;
                    }

                    Shiny.setInputValue('visible_permalink_panel', {
                        panel: panelId || window.__atlasRestoringPanelId || atlasVisiblePermalinkPanel() || '',
                        nonce: Date.now()
                    }, { priority: 'event' });
                }

                function atlasPinboardState() {
                    if (!window.__atlasPinboardState) {
                        window.__atlasPinboardState = {
                            items: [],
                            maxItems: 4,
                            nextId: 1
                        };
                    }
                    return window.__atlasPinboardState;
                }

                function atlasEnsurePinboardUi() {
                    var root = document.getElementById('atlas-pinboard');
                    if (!root) return null;

                    return {
                        root: root,
                        empty: root.querySelector('.atlas-pinboard-empty'),
                        grid: root.querySelector('.atlas-pinboard-grid'),
                        clear: root.querySelector('.atlas-pinboard-clear')
                    };
                }

                function atlasActiveTabLabel() {
                    var active = document.querySelector('.nav-tabs li.active a');
                    return active ? (active.textContent || '').trim() : 'Current tab';
                }

                function atlasSerializeSvg(svgEl) {
                    var svgMarkup = new XMLSerializer().serializeToString(svgEl);
                    return 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svgMarkup);
                }

                function atlasCapturePlotCard(card) {
                    var titleNode = card.querySelector('.plot-card-title');
                    var title = titleNode ? (titleNode.textContent || '').trim() : 'Pinned plot';
                    var context = atlasActiveTabLabel();

                    return new Promise(function(resolve, reject) {
                        var plotlyEl = card.querySelector('.js-plotly-plot');
                        if (plotlyEl && window.Plotly && window.Plotly.toImage) {
                            window.Plotly.toImage(plotlyEl, {
                                format: 'png',
                                width: Math.max(plotlyEl.clientWidth || 0, 880),
                                height: Math.max(plotlyEl.clientHeight || 0, 520)
                            }).then(function(dataUrl) {
                                resolve({ kind: 'img', src: dataUrl, title: title, context: context });
                            }).catch(function() {
                                reject(new Error('plotly-export-failed'));
                            });
                            return;
                        }

                        var imgEl = card.querySelector('img');
                        if (imgEl && (imgEl.currentSrc || imgEl.src)) {
                            resolve({ kind: 'img', src: imgEl.currentSrc || imgEl.src, title: title, context: context });
                            return;
                        }

                        var canvasEl = card.querySelector('canvas');
                        if (canvasEl && canvasEl.toDataURL) {
                            resolve({ kind: 'img', src: canvasEl.toDataURL('image/png'), title: title, context: context });
                            return;
                        }

                        var svgEl = card.querySelector('svg');
                        if (svgEl) {
                            resolve({ kind: 'img', src: atlasSerializeSvg(svgEl), title: title, context: context });
                            return;
                        }

                        reject(new Error('plot-not-ready'));
                    });
                }

                function atlasRenderPinboard() {
                    var ui = atlasEnsurePinboardUi();
                    if (!ui) return;

                    var state = atlasPinboardState();
                    ui.grid.innerHTML = '';

                    if (!state.items.length) {
                        ui.root.classList.remove('has-items');
                        return;
                    }

                    ui.root.classList.add('has-items');

                    state.items.forEach(function(item) {
                        var tile = document.createElement('div');
                        tile.className = 'atlas-pinboard-tile';

                        var head = document.createElement('div');
                        head.className = 'atlas-pinboard-tile-head';

                        var textWrap = document.createElement('div');
                        textWrap.className = 'atlas-pinboard-tile-meta';

                        var title = document.createElement('div');
                        title.className = 'atlas-pinboard-title';
                        title.textContent = item.title;

                        var context = document.createElement('div');
                        context.className = 'atlas-pinboard-context';
                        context.textContent = item.context;

                        var removeBtn = document.createElement('button');
                        removeBtn.type = 'button';
                        removeBtn.className = 'atlas-pinboard-remove';
                        removeBtn.setAttribute('data-pinboard-remove', item.id);
                        removeBtn.setAttribute('aria-label', 'Remove pinned plot');
                        removeBtn.textContent = 'Remove';

                        textWrap.appendChild(title);
                        textWrap.appendChild(context);
                        head.appendChild(textWrap);
                        head.appendChild(removeBtn);

                        var preview = document.createElement('div');
                        preview.className = 'atlas-pinboard-preview';

                        if (item.kind === 'img') {
                            var img = document.createElement('img');
                            img.src = item.src;
                            img.alt = item.title;
                            preview.appendChild(img);
                        }

                        tile.appendChild(head);
                        tile.appendChild(preview);
                        ui.grid.appendChild(tile);
                    });
                }

                function atlasPinPlotCard(card) {
                    if (!card) return;

                    atlasCapturePlotCard(card).then(function(snapshot) {
                        var state = atlasPinboardState();
                        snapshot.id = 'pin-' + state.nextId;
                        state.nextId += 1;

                        if (state.items.length >= state.maxItems) {
                            state.items.shift();
                            atlasShowToast('Scratchpad full. Replaced the oldest plot.', 'success');
                        }

                        state.items.push(snapshot);
                        atlasRenderPinboard();
                        atlasShowToast('Pinned plot', 'success');
                    }).catch(function() {
                        atlasShowToast('Wait for the plot to finish rendering before pinning it.', 'error');
                    });
                }

                function atlasDecoratePlotCards() {
                    var cards = Array.prototype.slice.call(document.querySelectorAll('.plot-card'));

                    cards.forEach(function(card) {
                        if (card.__atlasPinReady) return;

                        var title = card.querySelector(':scope > .plot-card-header .plot-card-title, :scope > .plot-card-title');
                        if (!title) return;

                        var header = card.querySelector(':scope > .plot-card-header');
                        if (!header) {
                            header = document.createElement('div');
                            header.className = 'plot-card-header';
                            title.parentNode.insertBefore(header, title);
                            header.appendChild(title);
                        }

                        var pinBtn = document.createElement('button');
                        pinBtn.type = 'button';
                        pinBtn.className = 'plot-card-pin-btn plot-download-btn';
                        pinBtn.textContent = 'Pin plot';
                        pinBtn.setAttribute('aria-label', 'Pin ' + ((title.textContent || '').trim() || 'plot') + ' to scratchpad');

                        var downloadBtn = header.querySelector('.plot-download-btn, .shiny-download-link');
                        if (downloadBtn) {
                            header.insertBefore(pinBtn, downloadBtn);
                        } else {
                            header.appendChild(pinBtn);
                        }

                        card.__atlasPinReady = true;
                    });
                }

                function atlasTourState() {
                    if (!window.__atlasTourState) {
                        window.__atlasTourState = {
                            active: false,
                            stepIndex: 0,
                            markSeen: true
                        };
                    }
                    return window.__atlasTourState;
                }

                function atlasTourSteps() {
                    return [
                        {
                            title: 'Pick a source species and enter genes',
                            body: 'Start here. Choose the species that supplies your gene IDs, then stage one or more genes in the search box.',
                            target: function() {
                                return document.querySelector('.source-species-picker') || document.getElementById('selected_genes');
                            }
                        },
                        {
                            title: 'Apply the panel',
                            body: 'Click Apply gene panel when the staged list looks right. That panel becomes the shared context for every atlas tab.',
                            target: function() {
                                return document.getElementById('apply_gene_selection');
                            }
                        },
                        {
                            title: 'Browse a species tab',
                            body: 'The Medicago, Glycine, and Lotus tabs show within-species structure, markers, and expression views for the current panel.',
                            target: function() {
                                var tabs = Array.prototype.slice.call(document.querySelectorAll('.nav-tabs li a'));
                                return tabs.find(function(node) {
                                    var text = (node.textContent || '').trim();
                                    return text === 'Medicago truncatula' || text === 'Glycine max' || text === 'Lotus japonicus';
                                }) || document.querySelector('.nav-tabs');
                            }
                        },
                        {
                            title: 'Compare across species',
                            body: 'Use the Camex and SATURN tabs to inspect ortholog mappings and shared cross-species embeddings for the same panel.',
                            target: function() {
                                var tabs = Array.prototype.slice.call(document.querySelectorAll('.nav-tabs li a'));
                                return tabs.find(function(node) {
                                    return /Cross-species/.test(node.textContent || '');
                                }) || document.querySelector('.nav-tabs');
                            }
                        },
                        {
                            title: 'Download plots and matrices',
                            body: 'Use the plot download buttons and matrix export controls to take the current visualizations and gene-by-cell summaries with you.',
                            target: function() {
                                return document.querySelector('.composite-dl-btn') || document.querySelector('.matrix-export-download-group');
                            }
                        }
                    ];
                }

                function atlasClearTourTarget() {
                    var current = document.querySelector('.atlas-tour-target');
                    if (current) current.classList.remove('atlas-tour-target');
                }

                function atlasEnsureTourUi() {
                    var root = document.getElementById('atlas-tour-root');
                    if (root) return root;

                    root = document.createElement('div');
                    root.id = 'atlas-tour-root';
                    root.innerHTML = [
                        '<div class=\"atlas-tour-backdrop\"></div>',
                        '<div class=\"atlas-tour-card\" role=\"dialog\" aria-modal=\"true\" aria-live=\"polite\">',
                        '  <div class=\"atlas-tour-step\"></div>',
                        '  <h3 class=\"atlas-tour-title\"></h3>',
                        '  <p class=\"atlas-tour-body\"></p>',
                        '  <div class=\"atlas-tour-actions\">',
                        '    <button type=\"button\" class=\"atlas-tour-btn atlas-tour-btn-secondary\" data-tour-action=\"back\">Back</button>',
                        '    <button type=\"button\" class=\"atlas-tour-btn atlas-tour-btn-secondary\" data-tour-action=\"skip\">Skip</button>',
                        '    <button type=\"button\" class=\"atlas-tour-btn atlas-tour-btn-primary\" data-tour-action=\"next\">Next</button>',
                        '  </div>',
                        '</div>'
                    ].join('');
                    document.body.appendChild(root);

                    root.addEventListener('click', function(event) {
                        var action = event.target && event.target.getAttribute('data-tour-action');
                        if (!action) return;

                        var state = atlasTourState();
                        var lastStep = atlasTourSteps().length - 1;

                        if (action === 'back') {
                            state.stepIndex = Math.max(0, state.stepIndex - 1);
                            atlasRenderTour();
                            return;
                        }

                        if (action === 'skip') {
                            atlasStopTour(false);
                            return;
                        }

                        if (state.stepIndex >= lastStep) {
                            atlasStopTour(true);
                            return;
                        }

                        state.stepIndex += 1;
                        atlasRenderTour();
                    });

                    return root;
                }

                function atlasRenderTour() {
                    var state = atlasTourState();
                    var steps = atlasTourSteps();
                    var step = steps[state.stepIndex];
                    var root = atlasEnsureTourUi();
                    var card = root.querySelector('.atlas-tour-card');
                    var stepEl = root.querySelector('.atlas-tour-step');
                    var titleEl = root.querySelector('.atlas-tour-title');
                    var bodyEl = root.querySelector('.atlas-tour-body');
                    var backBtn = root.querySelector('[data-tour-action=\"back\"]');
                    var nextBtn = root.querySelector('[data-tour-action=\"next\"]');
                    var target = step.target ? step.target() : null;

                    atlasClearTourTarget();

                    if (target) {
                        target.classList.add('atlas-tour-target');
                        target.scrollIntoView({ behavior: 'smooth', block: 'center' });
                    }

                    stepEl.textContent = 'Step ' + (state.stepIndex + 1) + ' of ' + steps.length;
                    titleEl.textContent = step.title;
                    bodyEl.textContent = step.body;
                    backBtn.disabled = state.stepIndex === 0;
                    nextBtn.textContent = state.stepIndex === steps.length - 1 ? 'Done' : 'Next';

                    root.className = 'atlas-tour-root is-active';

                    var defaultTop = Math.max(88, Math.round(window.innerHeight * 0.2));
                    var defaultLeft = Math.max(20, window.innerWidth - 380);

                    if (!target) {
                        card.style.top = defaultTop + 'px';
                        card.style.left = defaultLeft + 'px';
                        return;
                    }

                    var rect = target.getBoundingClientRect();
                    var cardWidth = 340;
                    var top = Math.max(88, Math.min(window.innerHeight - 220, rect.bottom + 18));
                    var left = rect.left;

                    if (left + cardWidth > window.innerWidth - 20) {
                        left = window.innerWidth - cardWidth - 20;
                    }
                    left = Math.max(20, left);

                    if (top > window.innerHeight - 220) {
                        top = Math.max(88, rect.top - 190);
                    }

                    card.style.top = Math.round(top) + 'px';
                    card.style.left = Math.round(left) + 'px';
                }

                function atlasStopTour(markSeen) {
                    var state = atlasTourState();
                    var root = document.getElementById('atlas-tour-root');

                    atlasClearTourTarget();

                    state.active = false;

                    if (root) {
                        root.className = 'atlas-tour-root';
                    }

                    if (markSeen && state.markSeen) {
                        try {
                            window.localStorage.setItem('atlas_tour_seen_v1', '1');
                        } catch (e) { /* localStorage blocked — ignore */ }
                    }
                }

                function atlasStartTour(markSeen) {
                    var overlay = document.getElementById('app-gate-overlay');
                    if (overlay && overlay.classList.contains('is-active')) {
                        window.__atlasPendingTourStart = !!markSeen;
                        return;
                    }

                    var state = atlasTourState();
                    state.active = true;
                    state.stepIndex = 0;
                    state.markSeen = markSeen !== false;

                    if (state.markSeen) {
                        try {
                            window.localStorage.setItem('atlas_tour_seen_v1', '1');
                        } catch (e) { /* localStorage blocked — ignore */ }
                    }

                    atlasRenderTour();
                }

                function initAtlasClient() {
                    if (!(window.Shiny && Shiny.setInputValue && Shiny.addCustomMessageHandler)) {
                        window.setTimeout(initAtlasClient, 250);
                        return;
                    }

                    if (window.__atlasClientInitialized) {
                        return;
                    }
                    window.__atlasClientInitialized = true;
                    window.__atlasPermalinkReportTimer = null;
                    window.__atlasRestoringPanelId = null;
                    window.__atlasBusyButtons = {};

                    atlasRenderPinboard();
                    atlasDecoratePlotCards();
                    [1500, 6000, 15000, 30000].forEach(function(delayMs) {
                        window.setTimeout(function() {
                            atlasDecoratePlotCards();
                            atlasReportVisiblePermalinkPanel();
                        }, delayMs);
                    });

                    try {
                        var seen = window.localStorage.getItem('atlas_tour_seen_v1');
                        if (!seen) {
                            window.setTimeout(function() {
                                atlasStartTour(true);
                            }, 900);
                        }
                    } catch (e) { /* localStorage blocked — skip tour */ }

                    var queryPayload = {};
                    var searchParams = new URLSearchParams(window.location.search);
                    searchParams.forEach(function(value, key) {
                        queryPayload[key] = value;
                    });
                    [1200, 5000, 15000, 30000].forEach(function(delayMs) {
                        window.setTimeout(function() {
                            var payload = Object.assign({}, queryPayload, { _nonce: Date.now() + delayMs });
                            Shiny.setInputValue('url_query_payload', payload, { priority: 'event' });
                            atlasReportVisiblePermalinkPanel();
                        }, delayMs);
                    });

                    var schedulePermalinkReport = function(delayMs) {
                        if (window.__atlasPermalinkReportTimer) {
                            window.clearTimeout(window.__atlasPermalinkReportTimer);
                        }
                        window.__atlasPermalinkReportTimer = window.setTimeout(function() {
                            atlasReportVisiblePermalinkPanel();
                        }, delayMs || 180);
                    };

                    var copyLink = document.getElementById('copy_permalink');
                    if (copyLink && !copyLink.__atlasBound) {
                        copyLink.__atlasBound = true;
                        copyLink.addEventListener('click', function(event) {
                            event.preventDefault();

                            var panelId = atlasVisiblePermalinkPanel();
                            var url = new URL(window.location.href);

                            if (panelId) {
                                url.searchParams.set('panel', panelId);
                            } else {
                                url.searchParams.delete('panel');
                            }

                            atlasReportVisiblePermalinkPanel(panelId);

                            navigator.clipboard.writeText(url.toString()).then(function() {
                                atlasShowToast('Copied permalink', 'success');
                            }).catch(function() {
                                atlasShowToast('Could not copy permalink', 'error');
                            });
                        });
                    }

                    var helpLink = document.getElementById('open_help_tour');
                    if (helpLink && !helpLink.__atlasTourBound) {
                        helpLink.__atlasTourBound = true;
                        helpLink.addEventListener('click', function(event) {
                            event.preventDefault();
                            atlasStartTour(false);
                        });
                    }

                    document.addEventListener('click', function(event) {
                        var tabLink = event.target && event.target.closest('.nav-tabs a');
                        if (tabLink) {
                            window.setTimeout(function() {
                                atlasReportVisiblePermalinkPanel();
                            }, 260);
                        }

                        var pinBtn = event.target && event.target.closest('.plot-card-pin-btn');
                        if (pinBtn) {
                            event.preventDefault();
                            atlasPinPlotCard(pinBtn.closest('.plot-card'));
                            return;
                        }

                        var removeBtn = event.target && event.target.closest('[data-pinboard-remove]');
                        if (removeBtn) {
                            event.preventDefault();
                            var state = atlasPinboardState();
                            var targetId = removeBtn.getAttribute('data-pinboard-remove');
                            state.items = state.items.filter(function(item) {
                                return item.id !== targetId;
                            });
                            atlasRenderPinboard();
                            return;
                        }

                        var clearBtn = event.target && event.target.closest('.atlas-pinboard-clear');
                        if (clearBtn) {
                            event.preventDefault();
                            atlasPinboardState().items = [];
                            atlasRenderPinboard();
                        }
                    });

                    window.addEventListener('scroll', function() {
                        schedulePermalinkReport(180);
                    }, { passive: true });
                    window.addEventListener('resize', function() {
                        schedulePermalinkReport(180);
                    });

                    document.addEventListener('shiny:idle', function() {
                        var busyButtons = window.__atlasBusyButtons || {};

                        Object.keys(busyButtons).forEach(function(buttonId) {
                            atlasSetButtonBusy(buttonId, false);
                        });

                        window.__atlasBusyButtons = {};
                        Shiny.setInputValue('atlas_client_idle', Date.now(), { priority: 'event' });
                    });

                    Shiny.addCustomMessageHandler('atlas_gate_unlock', function(msg) {
                        var overlay = document.getElementById('app-gate-overlay');
                        if (overlay) overlay.classList.remove('is-active');
                        if (window.__atlasPendingTourStart) {
                            window.__atlasPendingTourStart = false;
                            window.setTimeout(function() {
                                atlasStartTour(true);
                            }, 700);
                        }
                    });
                    Shiny.addCustomMessageHandler('atlas_restore_panel', function(msg) {
                        var panelId = msg && msg.panel ? msg.panel : null;
                        if (!panelId) return;

                        window.__atlasRestoringPanelId = panelId;
                        var attempts = 0;
                        var panelIsVisible = function(target) {
                            if (!target) return false;
                            var rect = target.getBoundingClientRect();
                            return rect.bottom > 80 && rect.top < window.innerHeight * 0.9;
                        };
                        var tryScroll = function() {
                            var target = atlasFindPermalinkPanel(panelId);
                            if (target && target.offsetParent !== null && panelIsVisible(target)) {
                                window.__atlasRestoringPanelId = null;
                                atlasReportVisiblePermalinkPanel(panelId);
                                return;
                            }

                            if (target && target.offsetParent !== null) {
                                target.scrollIntoView({ behavior: 'smooth', block: attempts < 2 ? 'start' : 'center' });
                                window.setTimeout(function() {
                                    var refreshedTarget = atlasFindPermalinkPanel(panelId);
                                    if (refreshedTarget && refreshedTarget.offsetParent !== null && panelIsVisible(refreshedTarget)) {
                                        window.__atlasRestoringPanelId = null;
                                        atlasReportVisiblePermalinkPanel(panelId);
                                        return;
                                    }
                                    if (attempts >= 40) {
                                        window.__atlasRestoringPanelId = null;
                                        atlasReportVisiblePermalinkPanel(panelId);
                                        return;
                                    }
                                    attempts += 1;
                                    window.setTimeout(tryScroll, 400);
                                }, 260);
                                return;
                            }

                            if (attempts >= 40) {
                                window.__atlasRestoringPanelId = null;
                                atlasReportVisiblePermalinkPanel(panelId);
                                return;
                            }
                            attempts += 1;
                            window.setTimeout(tryScroll, 400);
                        };

                        window.setTimeout(function() {
                            var initialTarget = atlasFindPermalinkPanel(panelId);
                            if (initialTarget && initialTarget.offsetParent !== null && panelIsVisible(initialTarget)) {
                                window.__atlasRestoringPanelId = null;
                                atlasReportVisiblePermalinkPanel(panelId);
                            } else {
                                tryScroll();
                            }
                        }, 250);
                    });
                    Shiny.addCustomMessageHandler('atlas_button_busy', function(msg) {
                        if (!msg || !msg.button_id) return;

                        var isBusy = !(msg.busy === false);
                        atlasSetButtonBusy(msg.button_id, isBusy, msg.label);

                        if (isBusy) {
                            window.__atlasBusyButtons[msg.button_id] = true;
                        } else if (window.__atlasBusyButtons) {
                            delete window.__atlasBusyButtons[msg.button_id];
                        }
                    });
                }

                if (document.readyState === 'loading') {
                    document.addEventListener('DOMContentLoaded', initAtlasClient);
                } else {
                    initAtlasClient();
                }
            })();"
        ))
    ),
    uiOutput("app_gate_ui"),
    tags$div(id = "app-gate-overlay", class = if (atlas_access_required) "app-gate-overlay is-active" else "app-gate-overlay"),
        div(
            class = "app-shell",
            div(
                class = "busy-indicator",
                role = "status",
                `aria-live` = "polite",
                `aria-atomic` = "true",
                div(class = "busy-indicator-card",
                    div(class = "busy-indicator-spinner", `aria-hidden` = "true"),
                    div(class = "busy-indicator-copy",
                        div(class = "busy-indicator-text", "Updating...")
                    )
                )
            ),
        uiOutput("startup_data_check"),
        tags$header(
            class = "app-brand-shell",
            div(
                class = "app-brand-lockup",
                div(class = "app-brand-glyph", `aria-hidden` = "true", icon("leaf")),
                div(
                    class = "app-brand-text-wrap",
                    div(class = "app-brand-title", "Legume Root Nodule Symbiosis Atlas")
                )
            )
        ),
        div(
            class = "hero-copy hero-copy-wide",
            div(class = "hero-eyebrow", "Single-cell explorer"),
            h1("An scRNA-seq atlas for determinate and indeterminate nodules, with cross-species comparison"),
            p(
                class = "hero-text",
                HTML("Search by gene ID or common name to compare expression across <em>Medicago truncatula</em>, <em>Glycine max</em>, <em>Lotus japonicus</em>, and the cross-species integration.")
            )
        ),
        uiOutput("atlas_summary_ui"),
        tags$main(
            id = "main-content",
            tabindex = "-1",
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
                        selectizeInput(
                            inputId = "selected_genes",
                            label = "Search by gene ID or common name",
                            choices = NULL,
                            multiple = TRUE,
                            options = list(
                                plugins = list("remove_button"),
                                valueField = "value",
                                labelField = "label",
                                searchField = c("label", "value", "tokens"),
                                searchConjunction = "or",
                                maxOptions = 50,
                                closeAfterSelect = FALSE,
                                hideSelected = TRUE,
                                placeholder = "Type a gene ID, common name, or synonym",
                                onType = I("function(query) { window.atlasScheduleBulkGeneOption(this, query); }"),
                                onDropdownOpen = I("function() { window.atlasScheduleBulkGeneOption(this, this.lastValue || this.$control_input.val()); }"),
                                onItemAdd = I("function(value) { window.atlasHandleBulkGeneOption(this, value); }")
                            )
                        ),
                        div(
                            class = "gene-action-row",
                            actionButton(
                                inputId = "apply_gene_selection",
                                label = "Apply gene panel",
                                icon = icon("play"),
                                class = "btn btn-default apply-selection-btn"
                            ),
                            actionButton(
                                inputId = "open_gene_import",
                                label = "Import list...",
                                icon = icon("file-import"),
                                class = "btn btn-default gene-import-btn"
                            )
                        ),
                        div(
                            class = "selection-meta",
                            `aria-live` = "polite",
                            `aria-atomic` = "true",
                            textOutput("gene_selection_status")
                        )
                    )
                ),
                column(
                    width = 2,
                    div(
                        class = "option-group",
                        selectInput(
                            inputId = "dl_format",
                            label = "Download format",
                            choices = c("SVG", "PNG", "PDF"),
                            selected = "PDF"
                        ),
                        selectInput(
                            inputId = "figure_preset",
                            label = "Figure style",
                            choices = c(
                                "Exploratory" = "exploratory",
                                "Presentation" = "presentation",
                                "Publication" = "publication"
                            ),
                            selected = "publication"
                        ),
                        downloadButton(
                            outputId = "dl_composite",
                            label = "Composite PDF",
                            icon = icon("file-pdf"),
                            class = "btn btn-default btn-sm composite-dl-btn"
                        ),
                        tags$div(
                            class = "option-toggle",
                            checkboxInput(
                                inputId = "colorblind_safe",
                                label = "Colorblind-safe expression colors",
                                value = FALSE
                            )
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
                                div(
                                    class = "plot-card-header",
                                    div(class = "plot-card-title", "Ortholog mapping summary"),
                                    downloadButton("dl_overview_mapping_table", "Download", class = "btn btn-default btn-sm plot-download-btn")
                                ),
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
        )
        ),
        div(
            id = "atlas-pinboard",
            class = "section-card atlas-pinboard",
            div(
                class = "section-header",
                div(class = "section-eyebrow", "Scratchpad"),
                h2("Pinned plots"),
                p("Pin up to four rendered plots from any tab to compare them side-by-side without losing the current view.")
            ),
            div(class = "atlas-pinboard-empty", "No plots pinned yet. Use “Pin plot” in any plot card to add the current view here."),
            div(class = "atlas-pinboard-actions",
                tags$button(type = "button", class = "plot-download-btn atlas-pinboard-clear", "Clear scratchpad")
            ),
            div(class = "atlas-pinboard-grid")
        ),
        tags$footer(
            class = "app-footer",
            div(
                class = "app-footer-row",
                div(
                    class = "app-footer-meta",
                    span(class = "app-footer-version", sprintf("Version %s", atlas_version)),
                    span(class = "app-footer-sep", aria_hidden = "true"),
                    span(class = "app-footer-updated", sprintf("Data updated %s", atlas_last_updated))
                ),
                div(
                    class = "app-footer-actions",
                    actionLink(inputId = "open_citation", label = "Cite this atlas", icon = icon("quote-right")),
                    actionLink(inputId = "copy_permalink", label = "Copy permalink", icon = icon("link")),
                    actionLink(inputId = "open_help_tour", label = "Help & walkthrough", icon = icon("circle-question"))
                )
            ),
            div(
                class = "app-footer-copy",
                "Single-cell atlas explorer for legume root nodule symbiosis.",
                strong(" 2026.")
            )
        )
    )
)

# ==============================================================================
# 3. SERVER LOGIC
# ==============================================================================

server <- function(input, output, session) {
    selection_notice <- reactiveVal(NULL)
    applied_selected_genes <- reactiveVal(character(0))
    gene_panel_busy_message <- reactiveVal(NULL)
    marker_job_messages <- reactiveVal(list())

    app_unlocked <- reactiveVal(!atlas_access_required)
    pending_url_genes <- reactiveVal(NULL)
    pending_url_cluster <- reactiveVal(NULL)
    pending_url_panel <- reactiveVal(NULL)
    restoring_permalink_panel <- reactiveVal(NULL)
    permalink_panel_state <- reactiveVal(NULL)

    output$app_gate_ui <- renderUI({
        if (isTRUE(app_unlocked())) return(NULL)
        div(
            class = "app-gate",
            role = "dialog",
            `aria-modal` = "true",
            `aria-labelledby` = "app-gate-title",
            div(
                class = "app-gate-card",
                div(class = "app-gate-icon", `aria-hidden` = "true", icon("leaf")),
                h2(id = "app-gate-title", class = "app-gate-title", "Pre-publication access"),
                p(
                    class = "app-gate-copy",
                    "This atlas is currently restricted while the accompanying manuscript is under peer review. Enter the access password shared with reviewers to continue."
                ),
                div(
                    class = "app-gate-form",
                    passwordInput("app_gate_password", label = "Password", placeholder = ""),
                    actionButton("app_gate_submit", "Unlock", class = "btn-primary"),
                    uiOutput("app_gate_error")
                )
            )
        )
    })

    observeEvent(input$app_gate_submit, {
        pw <- input$app_gate_password %||% ""
        if (identical(pw, atlas_access_password)) {
            app_unlocked(TRUE)
            session$sendCustomMessage("atlas_gate_unlock", TRUE)
        } else {
            output$app_gate_error <- renderUI({
                div(class = "app-gate-error", "Incorrect password.")
            })
        }
    })

    output$startup_data_check <- renderUI({
        missing <- check_missing_data_files()
        if (!length(missing)) {
            return(NULL)
        }

        bullets <- tags$ul(
            class = "startup-missing-list",
            lapply(missing, function(entry) {
                tags$li(tags$strong(entry$label), tags$code(entry$path))
            })
        )

        div(
            class = "alert-card warning startup-data-banner",
            div(class = "alert-title", "Atlas data is incomplete"),
            div(
                class = "alert-body",
                tagList(
                    tags$p("The following files were not found. Affected tabs will show an error until the data is in place:"),
                    bullets,
                    tags$p(
                        "If you are running a local copy, check that the data directory is mounted and that the ATLAS_DATA_DIR environment variable (if used) points to the right location."
                    )
                )
            )
        )
    })

    current_species_integration <- function(species_key) {
        input[[paste0(species_key, "_integration")]] %||% "ComBat_BBKNN"
    }

    set_marker_job_message <- function(prefix, message = NULL) {
        current <- marker_job_messages()

        if (is.null(message) || !nzchar(message)) {
            current[[prefix]] <- NULL
        } else {
            current[[prefix]] <- message
        }

        marker_job_messages(current)
    }

    update_selected_genes_input <- function(choice_bundle, selected = character(0)) {
        updateSelectizeInput(
            session = session,
            inputId = "selected_genes",
            choices = choice_bundle$choices,
            selected = selected,
            server = TRUE
        )
    }

    staged_source_genes <- reactive({
        genes <- input$selected_genes %||% character(0)
        genes <- unique(as.character(genes))
        genes[!grepl("^__bulk__::", genes)]
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
                        compact_display_gene_list(input$source_species %||% "medicago", dropped_genes, limit = 6)
                    )
                )
                showNotification(
                    ui = tagList(
                        tags$strong("Some genes were removed from the panel."),
                        tags$br(),
                        sprintf(
                            "%d gene(s) are not present in the new source atlas: %s",
                            length(dropped_genes),
                            compact_display_gene_list(input$source_species %||% "medicago", dropped_genes, limit = 6)
                        )
                    ),
                    type = "warning",
                    duration = 10
                )
            } else {
                selection_notice(NULL)
            }

            applied_selected_genes(valid_applied_selection)

            update_selected_genes_input(
                choice_bundle = choice_bundle,
                selected = valid_selection
            )
        },
        ignoreInit = FALSE
    )

    observeEvent(
        input$apply_gene_selection,
        {
            staged_genes <- staged_source_genes()
            gene_panel_busy_message(
                if (length(staged_genes)) {
                    sprintf(
                        "Refreshing plots for %d staged gene(s) across the atlas...",
                        length(staged_genes)
                    )
                } else {
                    "Clearing the current gene panel across the atlas..."
                }
            )
            session$sendCustomMessage("atlas_button_busy", list(
                button_id = "apply_gene_selection",
                label = if (length(staged_genes)) "Applying..." else "Clearing..."
            ))
            applied_selected_genes(staged_genes)
            session$onFlushed(function() {
                session$sendCustomMessage("atlas_button_busy", list(
                    button_id = "apply_gene_selection",
                    busy = FALSE
                ))
            }, once = TRUE)
        }
    )

    parse_pasted_gene_list <- function(raw) {
        if (is.null(raw) || !nzchar(raw)) return(character(0))
        tokens <- unlist(strsplit(raw, "[\\s,;\\t]+", perl = TRUE), use.names = FALSE)
        tokens <- trimws(tokens)
        tokens <- tokens[nzchar(tokens)]
        unique(tokens)
    }

    read_gene_csv <- function(path) {
        if (is.null(path) || !nzchar(path) || !file.exists(path)) return(character(0))
        sep <- if (grepl("\\.tsv$", path, ignore.case = TRUE)) "\t" else ","
        df <- tryCatch(
            utils::read.delim(path, sep = sep, stringsAsFactors = FALSE, check.names = FALSE, header = TRUE),
            error = function(e) NULL
        )
        if (is.null(df) || !ncol(df)) return(character(0))
        first_col <- df[[1]]
        first_col <- trimws(as.character(first_col))
        first_col <- first_col[nzchar(first_col)]
        unique(first_col)
    }

    observeEvent(input$open_citation, {
        showModal(modalDialog(
            title = "Cite this atlas",
            size = "m",
            easyClose = TRUE,
            tags$p("If you use data or visualizations from this explorer in your own work, please cite:"),
            tags$blockquote(class = "citation-block", atlas_citation_text),
            tags$p(
                class = "citation-meta",
                sprintf("App version %s — data snapshot %s.", atlas_version, atlas_last_updated)
            ),
            tags$p(
                class = "citation-meta",
                HTML("<em>A peer-reviewed publication describing this atlas is in preparation. This citation will be updated once the paper is published.</em>")
            ),
            footer = modalButton("Close")
        ))
    })

    url_restoring <- reactiveVal(TRUE)
    url_query_applied <- reactiveVal(FALSE)

    observeEvent(input$visible_permalink_panel, {
        payload <- input$visible_permalink_panel %||% list()
        panel_id <- payload[["panel"]] %||% ""
        panel_id <- trimws(as.character(panel_id[[1]] %||% panel_id))
        restoring_panel <- restoring_permalink_panel()

        if (!nzchar(panel_id)) {
            if (!is.null(restoring_panel) && nzchar(restoring_panel)) {
                return()
            }
            permalink_panel_state(NULL)
            return()
        }

        if (!is.null(restoring_panel) && nzchar(restoring_panel)) {
            if (!identical(panel_id, restoring_panel)) {
                return()
            }
            restoring_permalink_panel(NULL)
        }

        permalink_panel_state(panel_id)
    }, ignoreInit = TRUE)

    current_permalink_cluster <- reactive({
        current_tab <- input$main_tabs %||% "overview"

        if (current_tab %in% within_species_keys) {
            return(input[[paste0(current_tab, "_marker_cluster")]] %||% "")
        }

        if (current_tab %in% paste0("cross_", cross_integration_keys)) {
            return(input[[paste0(current_tab, "_marker_cluster")]] %||% "")
        }

        ""
    })

    build_query_string <- function(panel = permalink_panel_state()) {
        applied <- applied_selected_genes()
        species <- input$source_species %||% "medicago"
        tab <- input$main_tabs %||% "overview"
        cb <- if (isTRUE(input$colorblind_safe)) "1" else "0"

        parts <- c(
            sprintf("species=%s", utils::URLencode(species, reserved = TRUE)),
            sprintf("tab=%s", utils::URLencode(tab, reserved = TRUE)),
            sprintf("cb=%s", cb)
        )

        for (species_key in within_species_keys) {
            integration_method <- current_species_integration(species_key)

            if (!is.null(integration_method) && nzchar(integration_method)) {
                parts <- c(
                    parts,
                    sprintf(
                        "%s_integration=%s",
                        species_key,
                        utils::URLencode(integration_method, reserved = TRUE)
                    )
                )
            }
        }

        cluster_id <- current_permalink_cluster()
        if (!is.null(cluster_id) && nzchar(cluster_id)) {
            parts <- c(parts, sprintf("cluster=%s", utils::URLencode(cluster_id, reserved = TRUE)))
        }

        if (!is.null(panel) && nzchar(panel)) {
            parts <- c(parts, sprintf("panel=%s", utils::URLencode(panel, reserved = TRUE)))
        }

        if (length(applied)) {
            parts <- c(parts, sprintf("genes=%s", utils::URLencode(paste(applied, collapse = ","), reserved = TRUE)))
        }

        paste0("?", paste(parts, collapse = "&"))
    }

    observeEvent(input$url_query_payload, {
        if (isTRUE(url_query_applied())) {
            return()
        }

        query <- input$url_query_payload %||% list()

        if (!length(query)) {
            url_query_applied(TRUE)
            url_restoring(FALSE)
            return()
        }

        species <- query[["species"]]
        if (!is.null(species) && nzchar(species) && (species %in% names(species_choices) || species %in% unname(species_choices))) {
            shinyWidgets::updatePickerInput(session, "source_species", selected = species)
        }

        for (species_key in within_species_keys) {
            integration_param <- query[[paste0(species_key, "_integration")]]
            if (!is.null(integration_param) && nzchar(integration_param) && integration_param %in% unname(integration_choices)) {
                updateSelectInput(
                    session,
                    inputId = paste0(species_key, "_integration"),
                    selected = integration_param
                )
            }
        }

        tab <- query[["tab"]]
        if (!is.null(tab) && nzchar(tab)) {
            updateTabsetPanel(session, "main_tabs", selected = tab)
        }

        genes_param <- query[["genes"]]
        if (!is.null(genes_param) && nzchar(genes_param)) {
            pending_url_genes(strsplit(genes_param, ",", fixed = TRUE)[[1]])
        }

        cb <- query[["cb"]]
        if (!is.null(cb) && nzchar(cb)) {
            updateCheckboxInput(session, "colorblind_safe", value = identical(cb, "1"))
        }

        cluster_param <- query[["cluster"]]
        if (!is.null(cluster_param) && nzchar(cluster_param)) {
            pending_url_cluster(cluster_param)
        }

        panel_param <- query[["panel"]]
        if (!is.null(panel_param) && nzchar(panel_param)) {
            pending_url_panel(panel_param)
            restoring_permalink_panel(panel_param)
            permalink_panel_state(panel_param)
        }

        url_query_applied(TRUE)
        url_restoring(FALSE)
    }, ignoreInit = FALSE)

    observe({
        source_gene_catalog()
        pending <- pending_url_genes()
        if (is.null(pending) || !length(pending)) return()

        bundle <- build_gene_choices(
            input$source_species %||% "medicago",
            source_integration()
        )
        valid_ids <- bundle$feature_ids %||% character(0)
        pending_valid <- intersect(pending, valid_ids)
        if (length(pending_valid)) {
            update_selected_genes_input(
                choice_bundle = bundle,
                selected = pending_valid
            )
            applied_selected_genes(pending_valid)
        }
        pending_url_genes(NULL)
    })

    observeEvent(input$selected_genes_bulk_query, {
        payload <- input$selected_genes_bulk_query %||% list()
        query <- trimws(as.character(payload$query %||% ""))

        if (!nzchar(query)) {
            return()
        }

        source_species <- input$source_species %||% "medicago"
        bundle <- build_gene_choices(source_species, source_integration())
        matched_ids <- match_gene_choices(bundle, query)
        current_selection <- isolate(staged_source_genes())
        new_matches <- setdiff(matched_ids, current_selection)
        updated_selection <- unique(c(current_selection, matched_ids))

        update_selected_genes_input(
            choice_bundle = bundle,
            selected = updated_selection
        )

        if (!length(matched_ids)) {
            showNotification(
                sprintf("No genes in %s match \"%s\".", species_label(source_species), query),
                type = "warning",
                duration = 6
            )
            return()
        }

        if (!length(new_matches)) {
            showNotification(
                sprintf(
                    "All %d gene(s) matching \"%s\" are already staged.",
                    length(matched_ids),
                    query
                ),
                type = "message",
                duration = 6
            )
            return()
        }

        showNotification(
            sprintf(
                "Added %d gene(s) matching \"%s\" to the staged panel.",
                length(new_matches),
                query
            ),
            type = "message",
            duration = 6
        )
    }, ignoreInit = TRUE)

    observeEvent(input$atlas_client_idle, {
        gene_panel_busy_message(NULL)
        marker_job_messages(list())
    }, ignoreInit = TRUE)

    observe({
        if (isTRUE(url_restoring())) return()
        updateQueryString(build_query_string(), mode = "replace")
    })

    observe({
        panel_id <- pending_url_panel()

        if (isTRUE(url_restoring()) || is.null(panel_id) || !nzchar(panel_id) || !isTRUE(app_unlocked())) {
            return()
        }

        session$sendCustomMessage("atlas_restore_panel", list(panel = panel_id))
        pending_url_panel(NULL)
    })

    observeEvent(input$open_gene_import, {
        showModal(modalDialog(
            title = "Import gene list",
            size = "m",
            easyClose = TRUE,
            tags$p(
                class = "modal-hint",
                "Paste or upload a list of gene IDs (one per line, or separated by commas / tabs / semicolons). Genes are matched against the current source species; unrecognized genes are dropped."
            ),
            tags$label(`for` = "gene_import_paste", "Paste gene IDs"),
            tags$textarea(
                id = "gene_import_paste",
                class = "form-control gene-import-textarea",
                rows = 8,
                placeholder = "e.g.\nMedtr7g029290\nMedtr4g104370\n..."
            ),
            tags$div(class = "modal-divider", "or"),
            fileInput(
                inputId = "gene_import_file",
                label = "Upload CSV / TSV (first column is the gene ID)",
                accept = c(".csv", ".tsv", ".txt"),
                buttonLabel = "Choose file",
                placeholder = "No file selected"
            ),
            checkboxInput(
                inputId = "gene_import_replace",
                label = "Replace current selection instead of appending",
                value = FALSE
            ),
            footer = tagList(
                modalButton("Cancel"),
                actionButton("confirm_gene_import", "Import", class = "btn-primary")
            )
        ))
    })

    observeEvent(input$confirm_gene_import, {
        pasted <- parse_pasted_gene_list(input$gene_import_paste)
        file_df <- input$gene_import_file
        from_file <- if (!is.null(file_df) && nrow(file_df)) read_gene_csv(file_df$datapath[1]) else character(0)
        imported <- unique(c(pasted, from_file))

        if (!length(imported)) {
            showNotification("No gene IDs detected in the pasted text or uploaded file.", type = "warning", duration = 6)
            return()
        }

        bundle <- build_gene_choices(
            input$source_species %||% "medicago",
            source_integration()
        )
        valid_ids <- bundle$feature_ids

        if (is.null(valid_ids)) valid_ids <- character(0)

        matches <- imported[imported %in% valid_ids]

        lower_map <- setNames(valid_ids, tolower(valid_ids))
        case_insensitive_hits <- lower_map[tolower(imported[!(imported %in% valid_ids)])]
        case_insensitive_hits <- unname(case_insensitive_hits[!is.na(case_insensitive_hits)])
        matches <- unique(c(matches, case_insensitive_hits))

        unmatched <- setdiff(imported, matches)
        unmatched <- unmatched[!(tolower(unmatched) %in% tolower(matches))]

        current <- isolate(staged_source_genes())
        new_selection <- if (isTRUE(input$gene_import_replace)) matches else unique(c(current, matches))

        update_selected_genes_input(
            choice_bundle = bundle,
            selected = new_selection
        )

        removeModal()

        summary_msg <- sprintf(
            "Imported %d gene(s)%s. %d not recognized in the current source catalog%s.",
            length(matches),
            if (isTRUE(input$gene_import_replace)) " (replaced current selection)" else " (appended)",
            length(unmatched),
            if (length(unmatched)) paste0(": ", compact_display_gene_list(input$source_species %||% "medicago", unmatched, limit = 6)) else ""
        )
        showNotification(summary_msg, type = if (length(unmatched)) "warning" else "message", duration = 8)
    })

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

    overview_cross_resolutions <- reactive({
        source_species <- input$source_species %||% "medicago"
        genes <- selected_source_genes()
        setNames(
            lapply(cross_integration_keys, function(cross_key) {
                resolve_cross_integration_mapping(
                    source_species = source_species,
                    source_genes = genes,
                    cross_key = cross_key
                )
            }),
            cross_integration_keys
        )
    }) %>% bindCache(
        input$source_species %||% "medicago",
        selected_source_genes(),
        cache = "app"
    )

    output$gene_selection_status <- renderText({
        staged_genes <- staged_source_genes()
        applied_genes <- selected_source_genes()
        n_genes <- length(applied_genes)
        note <- selection_notice()
        busy_message <- gene_panel_busy_message()

        if (!is.null(busy_message) && nzchar(busy_message)) {
            return(busy_message)
        }

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

            summary_strip_tile(species_label_tag(species_key), summary)
        })

        cross_summary <- get_cross_dataset_summary(cross_integration_keys[[1]])

        cross_card <- summary_strip_tile(
            tags$span(style = "text-transform: none; letter-spacing: normal;", "Cross-species"),
            cross_summary
        )

        div(
            class = "atlas-summary-strip",
            tagList(within_cards, cross_card)
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
                Medicago = map_chr(orthogroup, ~ compact_display_gene_list("medicago", get_orthogroup_members(.x, "medicago"))),
                `Glycine max` = map_chr(orthogroup, ~ compact_display_gene_list("glycine", get_orthogroup_members(.x, "glycine"))),
                `Lotus japonicus` = map_chr(orthogroup, ~ compact_display_gene_list("lotus", get_orthogroup_members(.x, "lotus")))
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

    output$dl_overview_mapping_table <- downloadHandler(
        filename = function() paste0("ortholog_mapping_", Sys.Date(), ".csv"),
        content = function(file) {
            export_tbl <- add_export_provenance_columns(
                overview_mapping_table(),
                tab_label = "Overview",
                integration_label = "Cross-atlas overview",
                extra = list(atlas_export_type = "ortholog_mapping_table")
            )

            write.csv(export_tbl, file, row.names = FALSE)
        }
    )

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

        cards <- list(
            notice_card(
                title = "How to read ortholog mappings",
                body = "Orthogroup membership is not proof of one-to-one equivalence. One-to-many mappings can reflect paralogs, and missing members or features can come from orthology-table or atlas-coverage limits rather than true biological absence.",
                tone = "info"
            )
        )

        if (length(orthogroup_status$without_orthogroup)) {
            cards <- append(cards, list(
                notice_card(
                    title = "Selected genes without orthogroups",
                    body = paste(
                        compact_display_gene_list(source_species, orthogroup_status$without_orthogroup, limit = 8),
                        "These genes are absent from the current orthogroup table, so cross-species comparisons cannot be interpreted for them here."
                    ),
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
                            compact_display_gene_list(source_species, resolution$no_target_members, limit = 8)
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
                            compact_display_gene_list(source_species, resolution$missing_features, limit = 8)
                        ),
                        tone = "warning"
                    )
                ))
            }

        }

        cross_resolutions <- overview_cross_resolutions()

        for (cross_key in cross_integration_keys) {
            cross_res <- cross_resolutions[[cross_key]]
            cross_label <- cross_integration_label(cross_key)

            if (length(cross_res$no_target_members)) {
                cards <- append(cards, list(
                    notice_card(
                        title = paste(cross_label, "orthogroup gaps"),
                        body = paste(
                            "The orthogroups below lack usable target members in the",
                            cross_label,
                            "space for:",
                            compact_display_gene_list(source_species, cross_res$no_target_members, limit = 8)
                        ),
                        tone = "warning"
                    )
                ))
            }

            if (length(cross_res$missing_features)) {
                cards <- append(cards, list(
                    notice_card(
                        title = paste(cross_label, "feature gaps"),
                        body = paste(
                            "Orthologs were resolved but no matching features exist in the",
                            cross_label,
                            "object for:",
                            compact_display_gene_list(source_species, cross_res$missing_features, limit = 8)
                        ),
                        tone = "warning"
                    )
                ))
            }
        }

        if (!length(cards)) {
            cards <- list(
                notice_card(
                    title = "Mappings look consistent",
                    body = "The current gene panel resolves cleanly across the species-specific atlases and both cross-species integrations.",
                    tone = "info"
                )
            )
        }

        div(class = "alert-stack", tagList(cards))
    })

    get_ext <- function() {
        tolower(input$dl_format %||% "svg")
    }

    figure_preset <- reactive({
        input$figure_preset %||% "publication"
    })

    observeEvent(figure_preset(), {
        recommended_format <- switch(
            figure_preset(),
            exploratory = "PNG",
            presentation = "SVG",
            publication = "PDF",
            "PDF"
        )

        updateSelectInput(session, "dl_format", selected = recommended_format)
    }, ignoreInit = TRUE)

    selected_gene_export_summary <- function(limit = 4L) {
        source_species <- input$source_species %||% "medicago"
        gene_labels <- display_gene_labels(
            source_species,
            selected_source_genes(),
            include_gene_id_with_common = FALSE
        )

        if (!length(gene_labels)) {
            return("none")
        }

        compact_gene_list(gene_labels, limit = limit)
    }

    export_provenance_fields <- function(tab_label, integration_label = NULL, extra = list()) {
        base_fields <- list(
            atlas_version = atlas_version,
            atlas_exported_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
            atlas_source_species = species_label(input$source_species %||% "medicago"),
            atlas_tab = tab_label %||% (input$main_tabs %||% "overview"),
            atlas_integration = integration_label %||% "",
            atlas_selected_genes = selected_gene_export_summary()
        )

        c(base_fields, extra)
    }

    add_export_provenance_columns <- function(tbl, tab_label, integration_label = NULL, extra = list()) {
        tbl <- as_tibble(tbl)
        provenance_tbl <- as_tibble(export_provenance_fields(
            tab_label = tab_label,
            integration_label = integration_label,
            extra = extra
        ))

        if (!nrow(tbl)) {
            return(bind_cols(provenance_tbl[rep(1, 0), , drop = FALSE], tbl))
        }

        bind_cols(provenance_tbl[rep(1, nrow(tbl)), , drop = FALSE], tbl)
    }

    export_plot_caption <- function(tab_label, integration_label = NULL, extra = list()) {
        fields <- export_provenance_fields(tab_label, integration_label, extra = extra)
        extra_fields <- extra

        caption_parts <- c(
            paste0("Atlas ", fields$atlas_version),
            paste0("Exported ", fields$atlas_exported_at),
            paste0("Tab: ", fields$atlas_tab),
            paste0("Source species: ", fields$atlas_source_species),
            if (nzchar(fields$atlas_integration)) paste0("Integration: ", fields$atlas_integration) else NULL,
            paste0("Genes: ", fields$atlas_selected_genes),
            if (length(extra_fields)) {
                paste(
                    vapply(names(extra_fields), function(name) {
                        paste0(gsub("^atlas_", "", name), ": ", extra_fields[[name]])
                    }, character(1)),
                    collapse = " | "
                )
            } else {
                NULL
            }
        )

        paste(caption_parts[nzchar(caption_parts)], collapse = " | ")
    }

    add_plot_export_provenance <- function(plot_obj, tab_label, integration_label = NULL, extra = list()) {
        caption_text <- export_plot_caption(
            tab_label = tab_label,
            integration_label = integration_label,
            extra = extra
        )

        if (!nzchar(caption_text)) {
            return(plot_obj)
        }

        tryCatch(
            plot_obj + labs(caption = caption_text) + theme(
                plot.caption = element_text(
                    size = 8.5,
                    colour = app_palette["muted"],
                    hjust = 0
                )
            ),
            error = function(e) plot_obj
        )
    }

    apply_export_figure_preset <- function(plot_obj) {
        preset <- figure_preset()
        preset_theme <- app_plot_theme(preset = preset)

        tryCatch(
            plot_obj & preset_theme,
            error = function(e) plot_obj + preset_theme
        )
    }

    save_ggplot <- function(file, plot_obj, width, height, tab_label = NULL, integration_label = NULL, extra = list()) {
        preset_cfg <- figure_preset_config(figure_preset())
        plot_obj <- apply_export_figure_preset(plot_obj)
        plot_obj <- add_plot_export_provenance(
            plot_obj,
            tab_label = tab_label,
            integration_label = integration_label,
            extra = extra
        )

        ggsave(
            filename = file,
            plot = plot_obj,
            device = get_ext(),
            width = width * preset_cfg$width_scale,
            height = height * preset_cfg$height_scale,
            dpi = 300,
            limitsize = FALSE
        )
    }

    output$dl_composite <- downloadHandler(
        filename = function() {
            sprintf("nodule_atlas_composite_%s.pdf", format(Sys.time(), "%Y%m%d_%H%M"))
        },
        content = function(file) {
            source_species <- input$source_species %||% "medicago"
            source_genes <- selected_source_genes()

            if (!length(source_genes)) {
                pdf(file = file, width = 8.5, height = 4)
                plot.new()
                title(main = "No genes selected — apply a gene panel before exporting.")
                dev.off()
                return()
            }

            cb <- isTRUE(input$colorblind_safe)

            panels_per_species <- lapply(within_species_keys, function(sp) {
                integration_method <- current_species_integration(sp)
                resolution <- tryCatch(
                    resolve_target_mapping(
                        source_species = source_species,
                        source_genes = source_genes,
                        target_species = sp,
                        integration_method = integration_method,
                        cross_space = FALSE
                    ),
                    error = function(e) NULL
                )
                if (is.null(resolution) || !length(resolution$plot_features)) return(NULL)

                obj <- tryCatch(get_within_object(sp, integration_method), error = function(e) NULL)
                if (is.null(obj)) return(NULL)

                plots <- lapply(resolution$plot_features, function(feature_id) {
                    p <- tryCatch(
                        emphasized_feature_plot(
                            obj = obj,
                            feature_id = feature_id,
                            colorblind_safe = cb
                        ) +
                            app_plot_theme() +
                            labs(
                                title = sprintf(
                                    "%s — %s",
                                    species_registry[[sp]]$label,
                                    unname(resolution$label_map[feature_id]) %||% feature_id
                                ),
                                color = NULL
                            ) +
                            theme(
                                axis.title = element_blank(),
                                axis.text = element_blank(),
                                axis.ticks = element_blank(),
                                plot.title = element_text(face = "bold", size = 11)
                            ),
                        error = function(e) NULL
                    )
                    p
                })
                plots <- Filter(Negate(is.null), plots)
                if (!length(plots)) return(NULL)
                plots
            })

            all_panels <- unlist(panels_per_species, recursive = FALSE)
            all_panels <- Filter(Negate(is.null), all_panels)

            if (!length(all_panels)) {
                pdf(file = file, width = 8.5, height = 4)
                plot.new()
                title(main = "No panels could be generated for the current selection.")
                dev.off()
                return()
            }

            ncol <- min(length(source_genes), 3L)
            if (ncol < 1L) ncol <- 1L
            composite <- patchwork::wrap_plots(all_panels, ncol = ncol)
            composite <- add_plot_export_provenance(
                composite,
                tab_label = "Composite atlas export",
                integration_label = "Cross-atlas composite",
                extra = list(atlas_export_type = "composite_pdf")
            )
            composite <- apply_export_figure_preset(composite)

            page_width <- max(7, ncol * 3.4)
            nrow <- ceiling(length(all_panels) / ncol)
            page_height <- max(5, nrow * 3.0)
            preset_cfg <- figure_preset_config(figure_preset())

            ggsave(
                filename = file,
                plot = composite,
                device = "pdf",
                width = page_width * preset_cfg$width_scale,
                height = page_height * preset_cfg$height_scale,
                limitsize = FALSE
            )
        }
    )

    register_species_tab <- function(species_key) {
        local({
            prefix <- species_key
            tab_label <- species_label(species_key)
            clustering_columns <- distribution_cluster_columns

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
            }) %>% bindCache(
                input$source_species %||% "medicago",
                selected_source_genes(),
                species_key,
                current_species_integration(species_key),
                cache = "app"
            )

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
                                "Using orthogroup mappings from %s into the expression panels for %s. When several target genes fall in the same orthogroup, they are shown separately and should be treated as paralog-aware candidates rather than direct one-to-one equivalents.",
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
                            body = paste(
                                compact_display_gene_list(source_species, resolution$no_orthogroup, limit = 8),
                                "These genes have no orthogroup in the current table, so the absence of mapped panels is an orthology-coverage limit, not a statement about expression."
                            ),
                            tone = "warning"
                        )
                    ))
                }

                if (length(resolution$no_target_members)) {
                    cards <- append(cards, list(
                        notice_card(
                            title = paste("No", tab_label, "members in the orthogroups"),
                            body = paste(
                                compact_display_gene_list(source_species, resolution$no_target_members, limit = 8),
                                "Their orthogroups do not contain genes from this target species, so no cross-species inference can be made for those cases in this atlas."
                            ),
                            tone = "warning"
                        )
                    ))
                }

                if (length(resolution$missing_features)) {
                    cards <- append(cards, list(
                        notice_card(
                            title = "Mapped orthologs missing from the selected atlas",
                            body = paste(
                                compact_display_gene_list(source_species, resolution$missing_features, limit = 8),
                                "Orthologs were found, but the selected atlas lacks matching features for them. Treat this as annotation or feature-space coverage loss rather than proof of no expression."
                            ),
                            tone = "warning"
                        )
                    ))
                }

                if (nrow(resolution$multiplicity)) {
                    multiplicity_text <- resolution$multiplicity %>%
                        mutate(
                            label = paste0(
                                display_gene_labels(input$source_species %||% "medicago", source_gene, include_gene_id_with_common = FALSE),
                                " (", mapped_gene_count, " orthologs)"
                            )
                        ) %>%
                        pull(label)

                    cards <- append(cards, list(
                        notice_card(
                            title = "One-to-many orthologs",
                            body = paste(
                                compact_gene_list(multiplicity_text, limit = 6),
                                "Multiple target genes belong to the same orthogroup. Compare them separately and avoid interpreting any single one as the direct conserved counterpart without external evidence."
                            ),
                            tone = "warning"
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

            tab_expression_export_feature_map <- reactive({
                build_expression_export_feature_map(
                    resolution = tab_resolution(),
                    source_species = input$source_species %||% "medicago",
                    source_genes = selected_source_genes()
                )
            })

            tab_expression_export_requested_mode <- reactive({
                input[[paste0(prefix, "_matrix_mode")]] %||% "per_cell"
            })

            tab_expression_export_effective_mode <- reactive({
                resolve_expression_export_mode(
                    requested_mode = tab_expression_export_requested_mode(),
                    row_count = ncol(tab_object())
                )
            })

            output[[paste0(prefix, "_matrix_hint_ui")]] <- renderUI({
                genes <- selected_source_genes()

                if (!length(genes)) {
                    return(
                        tags$p(
                            class = "marker-status-hint",
                            sprintf("Select one or more source genes to export an expression matrix for %s.", tab_label)
                        )
                    )
                }

                feature_map <- tab_expression_export_feature_map()
                missing_gene_labels <- feature_map$gene_label[lengths(feature_map$feature_ids) == 0]
                hint_lines <- list()

                if (identical(tab_expression_export_requested_mode(), "per_cell") &&
                    identical(tab_expression_export_effective_mode(), "per_cluster")) {
                    hint_lines <- append(hint_lines, list(
                        tags$p(
                            class = "marker-status-hint",
                            sprintf(
                                "Per-cell export is capped at %s rows. %s currently has %s cells, so this download will switch to per-cluster means.",
                                format(expression_export_row_limit, big.mark = ","),
                                tab_label,
                                format(ncol(tab_object()), big.mark = ",")
                            )
                        )
                    ))
                }

                if (length(missing_gene_labels)) {
                    hint_lines <- append(hint_lines, list(
                        tags$p(
                            class = "marker-status-hint",
                            sprintf(
                                "Genes without mapped features in this tab export as blank columns: %s.",
                                compact_gene_list(missing_gene_labels, limit = 4)
                            )
                        )
                    ))
                }

                if (!length(hint_lines)) {
                    return(NULL)
                }

                div(class = "matrix-export-hints", tagList(hint_lines))
            })

            tab_marker_dataset_key <- reactive({
                paste(species_key, current_species_integration(species_key), sep = "_")
            })

            tab_markers_full <- reactive({
                read_cluster_markers_cache(tab_marker_dataset_key(), top_n = FALSE)
            })

            tab_marker_source <- reactive({
                marker_tbl <- tab_markers_full()

                if (is.null(marker_tbl) || !nrow(marker_tbl)) {
                    return(NA_character_)
                }

                available_sources <- unique(as.character(marker_tbl$cluster_source))
                preferred_source <- tab_distribution_group_by()
                if (!(length(preferred_source) && preferred_source %in% distribution_cluster_columns)) {
                    preferred_source <- NULL
                }

                select_marker_cluster_source(
                    available_sources = available_sources,
                    preferred_source = preferred_source,
                    default_source = "Rank_1st"
                )
            })

            tab_marker_lookup <- reactive({
                cluster_label_lookup(tab_object(), tab_marker_source())
            })

            tab_marker_choices <- reactive({
                marker_tbl <- tab_markers_full()

                if (is.null(marker_tbl) || !nrow(marker_tbl)) {
                    return(character(0))
                }

                marker_source <- tab_marker_source()

                marker_tbl <- marker_tbl %>%
                    filter(cluster_source == !!marker_source)

                available_clusters <- unique(as.character(marker_tbl$cluster))
                lookup_tbl <- tab_marker_lookup() %>%
                    filter(cluster %in% available_clusters)

                if (!nrow(lookup_tbl)) {
                    return(setNames(cluster_value_levels(available_clusters), cluster_value_levels(available_clusters)))
                }

                setNames(lookup_tbl$cluster, lookup_tbl$choice_label)
            })

            output[[paste0(prefix, "_marker_cluster_ui")]] <- renderUI({
                choices <- tab_marker_choices()

                if (!length(choices)) {
                    return(
                        tags$p(
                            class = "marker-status-hint",
                            tagList(
                                "Marker cache missing for this dataset. Run ",
                                tags$code("Rscript scripts/build_cluster_markers_cache.R"),
                                "."
                            )
                        )
                    )
                }

                selectInput(
                    inputId = paste0(prefix, "_marker_cluster"),
                    label = "Cluster",
                    choices = choices,
                    selected = resolve_choice(
                        input[[paste0(prefix, "_marker_cluster")]],
                        choices,
                        default = unname(choices[[1]])
                    )
                )
            })

            observeEvent(tab_marker_choices(), {
                pending_cluster <- pending_url_cluster()

                if (is.null(pending_cluster) || !nzchar(pending_cluster) || !identical(input$main_tabs %||% "overview", species_key)) {
                    return()
                }

                choices <- tab_marker_choices()

                if (!length(choices)) {
                    return()
                }

                selected_cluster <- resolve_choice(
                    pending_cluster,
                    choices,
                    default = unname(choices[[1]])
                )

                updateSelectInput(
                    session,
                    inputId = paste0(prefix, "_marker_cluster"),
                    selected = selected_cluster
                )
                pending_url_cluster(NULL)
            }, ignoreInit = FALSE)

            tab_marker_cluster <- reactive({
                choices <- tab_marker_choices()
                default_choice <- if (length(choices)) unname(choices[[1]]) else NULL

                resolve_choice(
                    input[[paste0(prefix, "_marker_cluster")]],
                    choices,
                    default = default_choice
                )
            })

            tab_marker_top_n <- reactive({
                marker_n <- suppressWarnings(as.integer(input[[paste0(prefix, "_marker_top_n")]] %||% 10L))

                if (is.na(marker_n) || marker_n < 1L) {
                    marker_n <- 10L
                }

                min(marker_n, 25L)
            })

            tab_marker_table_raw <- reactive({
                marker_tbl <- tab_markers_full()

                validate(
                    need(
                        !is.null(marker_tbl),
                        tagList(
                            "Cluster markers are not cached for this dataset. Run ",
                            tags$code("Rscript scripts/build_cluster_markers_cache.R"),
                            "."
                        )
                    ),
                    need(nrow(marker_tbl) > 0, "No cluster markers are available for this dataset.")
                )

                cluster_id <- tab_marker_cluster()
                marker_source <- tab_marker_source()
                filtered_tbl <- marker_tbl %>%
                    filter(cluster_source == !!marker_source) %>%
                    filter(cluster == !!cluster_id) %>%
                    arrange(p_val_adj, desc(avg_log2FC), desc(pct.1), gene)

                validate(need(nrow(filtered_tbl) > 0, "No marker rows are available for the selected cluster."))

                filtered_tbl
            })

            output[[paste0(prefix, "_markers_status_ui")]] <- renderUI({
                marker_tbl <- tab_markers_full()
                busy_message <- marker_job_messages()[[prefix]]

                if (is.null(marker_tbl) || !nrow(marker_tbl)) {
                    return(NULL)
                }

                if (!is.null(busy_message) && nzchar(busy_message)) {
                    return(tags$p(class = "marker-status-hint is-working", busy_message))
                }

                source_species <- input$source_species %||% "medicago"
                cluster_choices <- tab_marker_choices()
                marker_source <- tab_marker_source()
                active_group_by <- tab_distribution_group_by()
                current_cluster_label <- names(cluster_choices)[match(tab_marker_cluster(), unname(cluster_choices))]
                current_cluster_label <- current_cluster_label[!is.na(current_cluster_label) & nzchar(current_cluster_label)][1] %||% tab_marker_cluster()
                helper_text <- if (identical(source_species, species_key)) {
                    "Top markers will be added directly into the current source-gene panel."
                } else {
                    sprintf(
                        "Top markers will be mapped back into the %s source-gene panel through orthogroups before they are added.",
                        species_label(source_species)
                    )
                }
                source_sentence <- if (marker_cluster_source_matches(marker_source, active_group_by)) {
                    paste0("using ", marker_cluster_source_label(marker_source), ".")
                } else {
                    paste0(
                        "using ",
                        marker_cluster_source_label(marker_source),
                        " because ",
                        metadata_column_label(active_group_by),
                        " is not a cached clustering solution."
                    )
                }

                tags$p(
                    class = "marker-status-hint",
                    paste(
                        "Showing markers for",
                        current_cluster_label,
                        source_sentence,
                        helper_text
                    )
                )
            })

            output[[paste0(prefix, "_markers_table")]] <- DT::renderDT({
                marker_tbl <- tab_marker_table_raw()

                display_tbl <- marker_tbl %>%
                    transmute(
                        gene = display_gene_labels(species_key, gene),
                        avg_log2FC = avg_log2FC,
                        pct.1 = pct.1,
                        pct.2 = pct.2,
                        p_val_adj = p_val_adj
                    )

                DT::datatable(
                    display_tbl,
                    rownames = FALSE,
                    escape = TRUE,
                    selection = "none",
                    class = "stripe hover order-column compact",
                    options = list(
                        dom = "tip",
                        pageLength = 10,
                        autoWidth = TRUE,
                        order = list(list(1, "desc"))
                    )
                ) %>%
                    DT::formatRound(columns = c("avg_log2FC", "pct.1", "pct.2"), digits = 3) %>%
                    DT::formatSignif(columns = "p_val_adj", digits = 3)
            })

            observeEvent(input[[paste0(prefix, "_add_markers")]], {
                marker_tbl <- tab_marker_table_raw() %>%
                    slice_head(n = tab_marker_top_n())

                source_species <- input$source_species %||% "medicago"
                set_marker_job_message(
                    prefix,
                    sprintf(
                        "Adding the top %d markers from %s into the shared %s source-gene panel...",
                        nrow(marker_tbl),
                        tab_marker_cluster(),
                        species_label(source_species)
                    )
                )
                session$sendCustomMessage("atlas_button_busy", list(
                    button_id = paste0(prefix, "_add_markers"),
                    label = "Adding..."
                ))
                candidate_genes <- if (identical(species_key, source_species)) {
                    marker_tbl$gene
                } else {
                    map_target_genes_to_source_genes(species_key, marker_tbl$gene, source_species)
                }

                valid_source_ids <- build_gene_choices(
                    source_species,
                    source_integration()
                )$feature_ids %||% character(0)

                genes_to_add <- unique(candidate_genes[candidate_genes %in% valid_source_ids])

                if (!length(genes_to_add)) {
                    showNotification(
                        sprintf(
                            "No top markers from this cluster could be added to the %s source-gene panel.",
                            species_label(source_species)
                        ),
                        type = "warning",
                        duration = 8
                    )
                    return()
                }

                updated_selection <- unique(c(staged_source_genes(), genes_to_add))

                update_selected_genes_input(
                    choice_bundle = build_gene_choices(
                        source_species,
                        source_integration()
                    ),
                    selected = updated_selection
                )
                applied_selected_genes(updated_selection)

                skipped_n <- nrow(marker_tbl) - length(genes_to_add)
                showNotification(
                    paste0(
                        "Added ",
                        length(genes_to_add),
                        " marker gene(s) to the shared panel",
                        if (skipped_n > 0) paste0("; ", skipped_n, " could not be mapped into the current source species.") else "."
                    ),
                    type = "message",
                    duration = 6
                )
            }, ignoreInit = TRUE)

            output[[paste0("dl_", prefix, "_markers")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_cluster_", tab_marker_cluster(), "_markers.csv")
                },
                content = function(file) {
                    marker_tbl <- tab_marker_table_raw()
                    lookup_tbl <- tab_marker_lookup()
                    cluster_labels <- lookup_tbl$cluster_label[match(marker_tbl$cluster, lookup_tbl$cluster)]

                    export_tbl <- marker_tbl %>%
                        mutate(
                            cluster_label = ifelse(is.na(cluster_labels) | !nzchar(cluster_labels), cluster, cluster_labels),
                            gene_label = display_gene_labels(species_key, gene)
                        ) %>%
                        select(cluster, cluster_label, gene, gene_label, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
                        add_export_provenance_columns(
                            tab_label = tab_label,
                            integration_label = current_species_integration(species_key),
                            extra = list(atlas_export_type = "cluster_markers")
                        )

                    write.csv(export_tbl, file = file, row.names = FALSE, na = "")
                }
            )

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
            }) %>% bindCache(
                species_key,
                current_species_integration(species_key),
                tab_distribution_group_by(),
                tab_distribution_split_by(),
                as.numeric(input[[paste0(prefix, "_distribution_pt_size")]] %||% 1.1),
                cache = "app"
            )

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
                    {if (!is.null(fill_values)) scale_fill_manual(values = fill_values)} +
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
            }) %>% bindCache(
                species_key,
                current_species_integration(species_key),
                tab_composition_by(),
                tab_composition_cluster_by(),
                cache = "app"
            )

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
	                    feature_plot_args <- list(
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
                    if (isTRUE(input$colorblind_safe)) {
                        feature_plot_args$colors_use <- viridisLite::cividis(100, end = 0.95)
                    }
                    feature_plot <- do.call(scCustomize::FeaturePlot_scCustom, feature_plot_args)

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
	            }) %>% bindCache(
                species_key,
                current_species_integration(species_key),
                input$source_species %||% "medicago",
                selected_source_genes(),
                tab_split_by(),
                tab_composition_cluster_by(),
                as.numeric(input[[paste0(prefix, "_distribution_pt_size")]] %||% 1.1),
                isTRUE(input$colorblind_safe),
                cache = "app"
            )

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
            }) %>% bindCache(
                species_key,
                current_species_integration(species_key),
                input$source_species %||% "medicago",
                selected_source_genes(),
                tab_group_by(),
                cache = "app"
            )

            heatmap_plot_obj <- reactive({
                resolution <- expression_resolution()
                obj <- tab_object()

                validate(
                    need(
                        length(resolution$plot_features) > 0,
                        paste("No mapped genes are available for", tab_label, "in the selected atlas.")
                    )
                )

                build_expression_heatmap_plot(
                    obj = obj,
                    feature_ids = resolution$plot_features,
                    label_map = resolution$label_map,
                    group_by = tab_group_by(),
                    colorblind_safe = isTRUE(input$colorblind_safe)
                )
            }) %>% bindCache(
                species_key,
                current_species_integration(species_key),
                input$source_species %||% "medicago",
                selected_source_genes(),
                tab_group_by(),
                isTRUE(input$colorblind_safe),
                cache = "app"
            )

            ridge_plot_obj <- reactive({
                resolution <- expression_resolution()
                obj <- tab_object()

                validate(
                    need(
                        length(resolution$plot_features) > 0,
                        paste("No mapped genes are available for", tab_label, "in the selected atlas.")
                    )
                )

                build_expression_ridge_plot(
                    obj = obj,
                    feature_ids = resolution$plot_features,
                    label_map = resolution$label_map,
                    group_by = tab_group_by(),
                    colorblind_safe = isTRUE(input$colorblind_safe)
                )
            }) %>% bindCache(
                species_key,
                current_species_integration(species_key),
                input$source_species %||% "medicago",
                selected_source_genes(),
                tab_group_by(),
                isTRUE(input$colorblind_safe),
                cache = "app"
            )

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
	            }) %>% bindCache(
                species_key,
                current_species_integration(species_key),
                input$source_species %||% "medicago",
                selected_source_genes(),
                tab_group_by(),
                cache = "app"
            )

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

            output[[paste0(prefix, "_heatmap_plot")]] <- renderPlot(
                {
                    heatmap_plot_obj()
                },
                height = function() {
                    feature_n <- tryCatch(length(expression_resolution()$plot_features), error = function(e) 0L)
                    max(420, 110 + feature_n * 34)
                },
                res = 110
            )

            output[[paste0(prefix, "_ridge_plot")]] <- renderPlot(
                {
                    ridge_plot_obj()
                },
                height = function() {
                    feature_n <- tryCatch(length(expression_resolution()$plot_features), error = function(e) 0L)
                    group_n <- tryCatch(dplyr::n_distinct(tab_object()@meta.data[[tab_group_by()]]), error = function(e) 0L)
                    expression_ridge_height_px(feature_n = feature_n, group_n = group_n)
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

            output[[paste0("dl_", prefix, "_expression_matrix")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_expression_matrix_", tab_expression_export_effective_mode(), ".csv")
                },
                content = function(file) {
                    if (!length(selected_source_genes())) {
                        stop("Select one or more genes before downloading an expression matrix.", call. = FALSE)
                    }

                    export_tbl <- build_expression_export_table(
                        obj = tab_object(),
                        feature_map = tab_expression_export_feature_map(),
                        mode = tab_expression_export_effective_mode(),
                        default_species = species_label(species_key)
                    )

                    export_tbl <- add_export_provenance_columns(
                        export_tbl,
                        tab_label = tab_label,
                        integration_label = current_species_integration(species_key),
                        extra = list(
                            atlas_export_type = "expression_matrix",
                            atlas_matrix_mode = tab_expression_export_effective_mode()
                        )
                    )

                    utils::write.csv(export_tbl, file, row.names = FALSE, na = "")
                }
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
                        height = plot_height,
                        tab_label = tab_label,
                        integration_label = current_species_integration(species_key),
                        extra = list(atlas_export_type = "umap")
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
                        height = plot_height,
                        tab_label = tab_label,
                        integration_label = current_species_integration(species_key),
                        extra = list(atlas_export_type = "distribution_umap")
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
                        height = max(6, feature_n * 3.2),
                        tab_label = tab_label,
                        integration_label = current_species_integration(species_key),
                        extra = list(atlas_export_type = "violin")
                    )
                }
            )

            output[[paste0("dl_", prefix, "_heatmap")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_heatmap.", get_ext())
                },
                content = function(file) {
                    feature_n <- length(expression_resolution()$plot_features)

                    save_ggplot(
                        file = file,
                        plot_obj = heatmap_plot_obj(),
                        width = 13,
                        height = max(5.5, 1.8 + feature_n * 0.55),
                        tab_label = tab_label,
                        integration_label = current_species_integration(species_key),
                        extra = list(atlas_export_type = "heatmap")
                    )
                }
            )

            output[[paste0("dl_", prefix, "_ridge")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_ridgeplot.", get_ext())
                },
                content = function(file) {
                    feature_n <- length(expression_resolution()$plot_features)
                    group_n <- dplyr::n_distinct(tab_object()@meta.data[[tab_group_by()]])

                    save_ggplot(
                        file = file,
                        plot_obj = ridge_plot_obj(),
                        width = 10,
                        height = max(6, expression_ridge_height_px(feature_n, group_n) / 95),
                        tab_label = tab_label,
                        integration_label = current_species_integration(species_key),
                        extra = list(atlas_export_type = "ridgeplot")
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
                    preset_cfg <- figure_preset_config(figure_preset())
                    export_caption <- export_plot_caption(
                        tab_label = tab_label,
                        integration_label = current_species_integration(species_key),
                        extra = list(atlas_export_type = "dotplot")
                    )

                    if (identical(ext, "png")) {
                        png(
                            filename = file,
                            width = 1900 * preset_cfg$width_scale,
                            height = max(1200, clustered_dot_plot_height() * 2) * preset_cfg$height_scale,
                            res = 220
                        )
                    } else if (identical(ext, "svg")) {
                        svglite::svglite(
                            file = file,
                            width = 14 * preset_cfg$width_scale,
                            height = plot_height * preset_cfg$height_scale
                        )
                    } else {
                        pdf(
                            file = file,
                            width = 14 * preset_cfg$width_scale,
                            height = plot_height * preset_cfg$height_scale
                        )
                    }

                    draw_clustered_dot_plot(plot_obj)
                    grid::grid.text(
                        export_caption,
                        x = grid::unit(0.01, "npc"),
                        y = grid::unit(0.01, "npc"),
                        just = c("left", "bottom"),
                        gp = grid::gpar(fontsize = 8.5, col = unname(app_palette["muted"]))
                    )
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
            }) %>% bindCache(
                input$source_species %||% "medicago",
                selected_source_genes(),
                cross_key,
                cache = "app"
            )

            output[[paste0(prefix, "_group_by_ui")]] <- renderUI({
                choices <- cross_group_choices(cross_object())

                selectInput(
                    inputId = paste0(prefix, "_group_by"),
                    label = "Summarize dot plot by",
                    choices = choices,
                    selected = resolve_choice(
                        input[[paste0(prefix, "_group_by")]],
                        choices,
                        default = integration_cfg$default_group_by
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

            cross_dist_group_by <- reactive({
                resolve_choice(
                    input[[paste0(prefix, "_dist_group_by")]],
                    cross_distribution_group_choices(cross_object(), cross_key),
                    default = integration_cfg$default_group_by
                )
            })

            output[[paste0(prefix, "_dist_group_by_ui")]] <- renderUI({
                choices <- cross_distribution_group_choices(cross_object(), cross_key)
                selectInput(
                    inputId = paste0(prefix, "_dist_group_by"),
                    label = "Color cells by",
                    choices = choices,
                    selected = resolve_choice(
                        input[[paste0(prefix, "_dist_group_by")]],
                        choices,
                        default = integration_cfg$default_group_by
                    )
                )
            })

            output[[paste0(prefix, "_composition_by_ui")]] <- renderUI({
                choices <- cross_composition_choices(cross_object())
                if (!length(choices)) return(NULL)
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

            output[[paste0(prefix, "_dist_umap_plot_ui")]] <- renderUI({
                has_3d <- !is.null(get_cross_umap3d(cross_key))

                if (has_3d) {
                    div(
                        class = "distribution-split-layout",
                        div(
                            class = "distribution-split-panel",
                            div(class = "distribution-view-title", "2D UMAP"),
                            spinning_plot_output(
                                paste0(prefix, "_dist_umap_plot"),
                                proxy_height = "540px",
                                shell_class = "umap-plot-shell"
                            )
                        ),
                        div(
                            class = "distribution-split-panel",
                            div(class = "distribution-view-title", "3D UMAP"),
                            spinning_plotly_output(
                                paste0(prefix, "_dist_umap3d_plot"),
                                proxy_height = "540px",
                                shell_class = "plotly-plot-shell"
                            )
                        )
                    )
                } else {
                    spinning_plot_output(
                        paste0(prefix, "_dist_umap_plot"),
                        proxy_height = "540px",
                        shell_class = "umap-plot-shell"
                    )
                }
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

            output[[paste0(prefix, "_notice_ui")]] <- renderUI({
                genes <- selected_source_genes()
                resolution <- cross_resolution()
                feature_mode <- integration_cfg$feature_mode
                source_species <- input$source_species %||% "medicago"
                render_summary <- cross_render_load_summary(resolution, source_species, cross_key)

                intro_card <- if (identical(feature_mode, "medicago_space")) {
                    notice_card(
                        title = "Shared Medicago-space features",
                        body = paste(
                            cross_integration_label(cross_key),
                            "stores a shared feature space represented with Medicago identifiers. Soybean and Lotus selections are projected to Medicago orthologs in this tab, and neighborhood overlap in the integrated embedding should be treated as comparative structure rather than proof of conserved cell states."
                        ),
                        tone = "info"
                    )
                } else {
                    notice_card(
                        title = paste(cross_integration_label(cross_key), "ortholog feature space"),
                        body = paste(
                            cross_integration_label(cross_key),
                            "stores species-prefixed features. Source genes are resolved to ortholog features from Medicago, Glycine, and Lotus before the integrated plots are drawn, and the shared embedding should be interpreted as a comparative view rather than a one-to-one cell-state map."
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
                            body = paste(
                                compact_display_gene_list(input$source_species %||% "medicago", resolution$no_orthogroup, limit = 8),
                                "These genes are not represented in the current orthogroup table, so the cross-species panels cannot speak to them."
                            ),
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
                                compact_display_gene_list(input$source_species %||% "medicago", resolution$no_target_members, limit = 8)
                            } else {
                                paste(
                                    "No ortholog members from Medicago, Glycine, or Lotus were found in the mapped orthogroups for:",
                                    compact_display_gene_list(input$source_species %||% "medicago", resolution$no_target_members, limit = 8),
                                    "Treat this as an orthogroup-content limit, not as evidence that the biology is absent in the other species."
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
                            body = paste(
                                compact_display_gene_list(input$source_species %||% "medicago", resolution$missing_features, limit = 8),
                                "Orthologs were resolved but the integration feature set does not contain them. Treat this as feature-space coverage loss rather than proof of no expression."
                            ),
                            tone = "warning"
                        )
                    ))
                }

                if (nrow(resolution$multiplicity)) {
                    multiplicity_text <- resolution$multiplicity %>%
                        mutate(
                            label = paste0(
                                display_gene_labels(source_species, source_gene, include_gene_id_with_common = FALSE),
                                " (", mapped_gene_count, " mapped features)"
                            )
                        ) %>%
                        pull(label)

                    cards <- append(cards, list(
                        notice_card(
                            title = paste("One-to-many mappings in", cross_integration_label(cross_key)),
                            body = paste(
                                compact_gene_list(multiplicity_text, limit = 6),
                                "These are ambiguous orthogroup mappings. Each mapped gene is plotted separately, and none should be treated as the definitive conserved counterpart without outside support."
                            ),
                            tone = "warning"
                        )
                    ))
                }

                if (length(genes) &&
                    (render_summary$panel_count > length(unique(genes)) ||
                        render_summary$mapped_feature_count > length(unique(genes)))) {
                    cards <- append(cards, list(
                        notice_card(
                            title = paste("Large orthogroup expansion in", cross_integration_label(cross_key)),
                            body = tagList(
                                tags$p(
                                    sprintf(
                                        "This view will render %s comparison blocks from %s selected source genes and %s mapped ortholog features. Every ortholog member is kept, so larger orthogroups will take longer to load.",
                                        format(render_summary$panel_count, big.mark = ","),
                                        format(length(unique(genes)), big.mark = ","),
                                        format(render_summary$mapped_feature_count, big.mark = ",")
                                    )
                                ),
                                if (length(render_summary$expansion_labels)) {
                                    tags$p(
                                        paste(
                                            "Largest expansions:",
                                            compact_gene_list(render_summary$expansion_labels, limit = 4),
                                            "."
                                        )
                                    )
                                }
                            ),
                            tone = "warning"
                        )
                    ))
                }

                div(class = "alert-stack", tagList(cards))
            })

            cross_expression_export_feature_map <- reactive({
                build_expression_export_feature_map(
                    resolution = cross_resolution(),
                    source_species = input$source_species %||% "medicago",
                    source_genes = selected_source_genes()
                )
            })

            cross_expression_export_requested_mode <- reactive({
                input[[paste0(prefix, "_matrix_mode")]] %||% "per_cell"
            })

            cross_expression_export_effective_mode <- reactive({
                resolve_expression_export_mode(
                    requested_mode = cross_expression_export_requested_mode(),
                    row_count = ncol(cross_object())
                )
            })

            output[[paste0(prefix, "_matrix_hint_ui")]] <- renderUI({
                genes <- selected_source_genes()

                if (!length(genes)) {
                    return(
                        tags$p(
                            class = "marker-status-hint",
                            sprintf("Select one or more source genes to export an expression matrix for %s.", cross_integration_label(cross_key))
                        )
                    )
                }

                feature_map <- cross_expression_export_feature_map()
                missing_gene_labels <- feature_map$gene_label[lengths(feature_map$feature_ids) == 0]
                hint_lines <- list()

                if (identical(cross_expression_export_requested_mode(), "per_cell") &&
                    identical(cross_expression_export_effective_mode(), "per_cluster")) {
                    hint_lines <- append(hint_lines, list(
                        tags$p(
                            class = "marker-status-hint",
                            sprintf(
                                "Per-cell export is capped at %s rows. %s currently has %s cells, so this download will switch to per-cluster means.",
                                format(expression_export_row_limit, big.mark = ","),
                                cross_integration_label(cross_key),
                                format(ncol(cross_object()), big.mark = ",")
                            )
                        )
                    ))
                }

                if (length(missing_gene_labels)) {
                    hint_lines <- append(hint_lines, list(
                        tags$p(
                            class = "marker-status-hint",
                            sprintf(
                                "Genes without mapped features in this tab export as blank columns: %s.",
                                compact_gene_list(missing_gene_labels, limit = 4)
                            )
                        )
                    ))
                }

                if (!length(hint_lines)) {
                    return(NULL)
                }

                div(class = "matrix-export-hints", tagList(hint_lines))
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

            ortholog_trace_source_genes <- reactive({
                genes <- selected_source_genes()

                if (!length(genes)) {
                    return(character(0))
                }

                genes
            })

            output[[paste0(prefix, "_ortholog_trace_notice_ui")]] <- renderUI({
                selected_n <- length(selected_source_genes())

                if (selected_n <= 10L) {
                    return(NULL)
                }

                div(
                    class = "alert-stack",
                    notice_card(
                        title = "Rendering the full ortholog trace",
                        body = sprintf(
                            "Showing all %d selected genes in the ortholog trace. Larger source panels and orthogroups now stay expanded here, so this section can take longer to render.",
                            selected_n
                        ),
                        tone = "warning"
                    )
                )
            })

            ortholog_trace_specs <- reactive({
                source_species <- input$source_species %||% "medicago"
                source_genes <- ortholog_trace_source_genes()

                validate(
                    need(
                        length(source_genes) > 0,
                        "Add source-species genes to generate the ortholog trace panel."
                    )
                )

                lapply(source_genes, function(source_gene) {
                    source_label <- display_gene_labels(source_species, source_gene)[[1]]
                    orthogroup_id <- resolve_source_orthogroups(source_species, source_gene) %>%
                        filter(!is.na(orthogroup) & nzchar(orthogroup)) %>%
                        distinct(orthogroup) %>%
                        slice_head(n = 1) %>%
                        pull(orthogroup) %>%
                        first_nonempty()

                    species_specs <- lapply(within_species_keys, function(species_key) {
                        integration_method <- current_species_integration(species_key)

                        target_gene <- if (identical(species_key, source_species)) {
                            source_gene
                        } else if (is.na(orthogroup_id) || !nzchar(orthogroup_id)) {
                            NA_character_
                        } else {
                            first_nonempty(get_orthogroup_members(orthogroup_id, species_key))
                        }

                        target_label <- if (!is.na(target_gene) && nzchar(target_gene)) {
                            display_gene_labels(species_key, target_gene)[[1]]
                        } else {
                            NA_character_
                        }

                        if (is.na(target_gene) || !nzchar(target_gene)) {
                            return(list(
                                species_key = species_key,
                                integration_method = integration_method,
                                status = "missing_ortholog",
                                note = "No ortholog in this species",
                                target_gene = NA_character_,
                                target_label = NA_character_,
                                feature_id = NA_character_
                            ))
                        }

                        feature_id <- match_target_features(
                            target_species = species_key,
                            target_genes = target_gene,
                            integration_method = integration_method,
                            cross_space = FALSE
                        )[1] %||% NA_character_

                        if (is.na(feature_id) || !nzchar(feature_id)) {
                            return(list(
                                species_key = species_key,
                                integration_method = integration_method,
                                status = "missing_feature",
                                note = if (identical(species_key, source_species)) {
                                    "Gene not present in this atlas"
                                } else {
                                    "Ortholog not present in this atlas"
                                },
                                target_gene = target_gene,
                                target_label = target_label,
                                feature_id = NA_character_
                            ))
                        }

                        list(
                            species_key = species_key,
                            integration_method = integration_method,
                            status = "ok",
                            note = NA_character_,
                            target_gene = target_gene,
                            target_label = target_label,
                            feature_id = feature_id
                        )
                    })

                    list(
                        source_gene = source_gene,
                        source_label = source_label,
                        orthogroup = orthogroup_id,
                        title = if (!is.na(orthogroup_id) && nzchar(orthogroup_id)) {
                            paste0(source_label, " via orthogroup ", orthogroup_id)
                        } else {
                            paste0(source_label, " via orthogroup unavailable")
                        },
                        species_specs = species_specs
                    )
                })
            }) %>% bindCache(
                cross_key,
                input$source_species %||% "medicago",
                ortholog_trace_source_genes(),
                current_species_integration("medicago"),
                current_species_integration("glycine"),
                current_species_integration("lotus"),
                cache = "app"
            )

            ortholog_trace_plot_obj <- reactive({
                row_specs <- ortholog_trace_specs()
                pt_size <- max(0.65, as.numeric(input[[paste0(prefix, "_pt_size")]] %||% 0.45))

                row_plots <- lapply(row_specs, function(row_spec) {
                    expression_max <- max(unlist(lapply(row_spec$species_specs, function(spec) {
                        if (!identical(spec$status, "ok")) {
                            return(NA_real_)
                        }

                        obj <- get_within_object(spec$species_key, spec$integration_method)
                        expr_values <- aggregate_feature_expression(obj, spec$feature_id)

                        if (!length(expr_values)) {
                            return(NA_real_)
                        }

                        max(expr_values, na.rm = TRUE)
                    })), na.rm = TRUE)

                    if (!is.finite(expression_max) || expression_max <= 0) {
                        expression_max <- 1
                    }

                    species_panels <- lapply(row_spec$species_specs, function(spec) {
                        if (!identical(spec$status, "ok")) {
                            empty_panel <- empty_umap_message_plot(spec$note) +
                                theme(plot.margin = margin(4, 8, 6, 6))

                            return(add_species_caption(
                                empty_panel,
                                spec$species_key,
                                sublabel = spec$target_label
                            ))
                        }

                        obj <- get_within_object(spec$species_key, spec$integration_method)
                        feature_plot <- emphasized_feature_plot(
                            obj = obj,
                            feature_id = spec$feature_id,
                            fixed_max = expression_max,
                            pt_size = pt_size,
                            colorblind_safe = isTRUE(input$colorblind_safe),
                            expressing_size_boost = 0
                        ) +
                            labs(title = NULL, color = NULL) +
                            app_plot_theme() +
                            compact_feature_legend_guides() +
                            compact_feature_legend_theme() +
                            scale_x_continuous(expand = expansion(mult = 0.01)) +
                            scale_y_continuous(expand = expansion(mult = 0.01)) +
                            theme(
                                plot.title = element_blank(),
                                panel.grid = element_blank(),
                                axis.title = element_blank(),
                                axis.text = element_blank(),
                                axis.ticks = element_blank(),
                                axis.line = element_blank(),
                                plot.margin = margin(4, 8, 6, 6)
                            )

                        add_species_caption(
                            feature_plot,
                            spec$species_key,
                            sublabel = spec$target_label
                        )
                    })

                    wrap_titled_plot(
                        plot_obj = wrap_plots(
                            plotlist = species_panels,
                            ncol = 3,
                            guides = "collect"
                        ) & theme(legend.position = "top"),
                        title = row_spec$title
                    )
                })

                wrap_plots(plotlist = row_plots, ncol = 1)
            })

            cross_umap_plot_obj <- reactive({
                resolution <- cross_resolution()
                obj <- cross_object()
                panel_specs <- cross_comparison_panel_specs(resolution, cross_key)
                block_cols <- min(2L, max(1L, as.integer(input[[paste0(prefix, "_umap_columns")]] %||% 1L)))
                pt_size <- max(0.65, as.numeric(input[[paste0(prefix, "_pt_size")]] %||% 0.45))
                species_cells <- lapply(within_species_keys, function(species_key) {
                    cross_species_cell_ids(obj, species_key)
                })
                names(species_cells) <- within_species_keys

                validate(
                    need(
                        nrow(panel_specs) > 0,
                        paste(
                            "No selected genes resolve to features in the",
                            cross_integration_label(cross_key),
                            "integration."
                        )
                    )
                )

                reference_plot <- scCustomize::DimPlot_scCustom(
                    seurat_object = apply_metadata_display_order(obj, "species"),
                    colors_use = unname(distribution_color_map(obj@meta.data$species, "species")),
                    group.by = "species",
                    pt.size = max(0.45, pt_size * 0.85),
                    label = FALSE,
                    raster = TRUE
                ) &
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

                reference_panel <- wrap_titled_plot(
                    plot_obj = reference_plot + labs(title = NULL, color = NULL),
                    title = "Species overview"
                )

                comparison_blocks <- lapply(seq_len(nrow(panel_specs)), function(panel_idx) {
                    spec <- panel_specs[panel_idx, , drop = FALSE]
                    species_entries <- setNames(lapply(within_species_keys, function(species_key) {
                        spec[[paste0(species_key, "_panels")]][[1]]
                    }), within_species_keys)
                    species_expr <- setNames(lapply(within_species_keys, function(species_key) {
                        entries <- species_entries[[species_key]]

                        if (!panel_feature_count(entries)) {
                            return(list())
                        }

                        lapply(entries$feature_id, function(feature_id) {
                            aggregate_feature_expression(obj, feature_id)
                        })
                    }), within_species_keys)

                    expression_max <- max(unlist(lapply(within_species_keys, function(species_key) {
                        entries <- species_entries[[species_key]]
                        expr_list <- species_expr[[species_key]]
                        cell_ids <- species_cells[[species_key]]

                        if (!panel_feature_count(entries) || !length(cell_ids)) {
                            return(NA_real_)
                        }

                        vapply(seq_len(nrow(entries)), function(idx) {
                            expr_values <- expr_list[[idx]]

                            if (!length(expr_values)) {
                                return(NA_real_)
                            }

                            max(expr_values[cell_ids], na.rm = TRUE)
                        }, numeric(1))
                    })), na.rm = TRUE)

                    if (!is.finite(expression_max) || expression_max <= 0) {
                        expression_max <- 1
                    }

                    species_columns <- lapply(within_species_keys, function(species_key) {
                        entries <- species_entries[[species_key]]
                        expr_list <- species_expr[[species_key]]

                        if (!panel_feature_count(entries)) {
                            return(wrap_titled_plot(
                                plot_obj = empty_umap_message_plot("No mapped feature"),
                                title = species_label(species_key)
                            ))
                        }

                        feature_panels <- lapply(seq_len(nrow(entries)), function(idx) {
                            entry <- entries[idx, , drop = FALSE]

                            feature_plot <- emphasized_feature_plot(
                                obj = obj,
                                feature_id = entry$feature_id[[1]],
                                expression_values = expr_list[[idx]],
                                cell_ids = species_cells[[species_key]],
                                fixed_max = expression_max,
                                pt_size = pt_size,
                                colorblind_safe = isTRUE(input$colorblind_safe),
                                expressing_size_boost = 0
                            )

                            feature_plot <- feature_plot +
                                labs(title = NULL, color = NULL) +
                                app_plot_theme() +
                                compact_feature_legend_guides() +
                                compact_feature_legend_theme() +
                                scale_x_continuous(expand = expansion(mult = 0.01)) +
                                scale_y_continuous(expand = expansion(mult = 0.01)) +
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
                                title = entry$feature_label[[1]],
                                size = 12,
                                hjust = 0.5
                            )
                        })

                        wrap_plots(plotlist = feature_panels, ncol = 1, guides = "collect") +
                            plot_title_annotation(species_label(species_key), size = 15, hjust = 0.5)
                    })

                    comparison_block <- wrap_plots(
                        plotlist = c(list(reference_panel), species_columns),
                        ncol = 4,
                        widths = c(0.72, 1, 1, 1),
                        guides = "collect"
                    ) &
                        theme(legend.position = "top")

                    comparison_block + plot_title_annotation(as.character(spec$title[[1]]))
                })

                wrap_plots(plotlist = comparison_blocks, ncol = block_cols)
            }) %>% bindCache(
                cross_key,
                input$source_species %||% "medicago",
                selected_source_genes(),
                as.integer(input[[paste0(prefix, "_umap_columns")]] %||% 1L),
                as.numeric(input[[paste0(prefix, "_pt_size")]] %||% 0.45),
                isTRUE(input$colorblind_safe),
                cache = "app"
            )

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
            }) %>% bindCache(
                cross_key,
                input$source_species %||% "medicago",
                selected_source_genes(),
                cross_group_by(),
                cache = "app"
            )

            cross_heatmap_plot_obj <- reactive({
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

                build_expression_heatmap_plot(
                    obj = obj,
                    feature_ids = resolution$plot_features,
                    label_map = resolution$label_map,
                    group_by = cross_group_by(),
                    colorblind_safe = isTRUE(input$colorblind_safe)
                )
            }) %>% bindCache(
                cross_key,
                input$source_species %||% "medicago",
                selected_source_genes(),
                cross_group_by(),
                isTRUE(input$colorblind_safe),
                cache = "app"
            )

            cross_ridge_plot_obj <- reactive({
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

                build_expression_ridge_plot(
                    obj = obj,
                    feature_ids = resolution$plot_features,
                    label_map = resolution$label_map,
                    group_by = cross_group_by(),
                    colorblind_safe = isTRUE(input$colorblind_safe)
                )
            }) %>% bindCache(
                cross_key,
                input$source_species %||% "medicago",
                selected_source_genes(),
                cross_group_by(),
                isTRUE(input$colorblind_safe),
                cache = "app"
            )

            output[[paste0(prefix, "_umap_plot")]] <- renderPlot(
                {
                    cross_umap_plot_obj()
                },
                height = function() {
                    panel_specs <- tryCatch(cross_comparison_panel_specs(cross_resolution(), cross_key), error = function(e) tibble())
                    panel_n <- nrow(panel_specs)
                    block_cols <- min(2L, max(1L, as.integer(input[[paste0(prefix, "_umap_columns")]] %||% 1L)))
                    row_units <- tryCatch({
                        block_units <- cross_comparison_height_units(panel_specs)

                        if (!length(block_units)) {
                            integer(0)
                        } else {
                            vapply(split(block_units, ceiling(seq_along(block_units) / block_cols)), max, integer(1))
                        }
                    }, error = function(e) integer(0))

                    if (!length(row_units)) {
                        return(560L)
                    }

                    max(560L, as.integer(180L + sum(row_units) * 310L))
                },
                res = 110
            )

            output[[paste0(prefix, "_ortholog_trace_plot")]] <- renderPlot(
                {
                    ortholog_trace_plot_obj()
                },
                height = function() {
                    ortholog_trace_height_px(length(ortholog_trace_source_genes()))
                },
                res = 110
            )

            output[[paste0(prefix, "_heatmap_plot")]] <- renderPlot(
                {
                    cross_heatmap_plot_obj()
                },
                height = function() {
                    feature_n <- tryCatch(length(cross_resolution()$plot_features), error = function(e) 0L)
                    max(420, 110 + feature_n * 34)
                },
                res = 110
            )

            output[[paste0(prefix, "_ridge_plot")]] <- renderPlot(
                {
                    cross_ridge_plot_obj()
                },
                height = function() {
                    feature_n <- tryCatch(length(cross_resolution()$plot_features), error = function(e) 0L)
                    group_n <- tryCatch(dplyr::n_distinct(cross_object()@meta.data[[cross_group_by()]]), error = function(e) 0L)
                    expression_ridge_height_px(feature_n = feature_n, group_n = group_n)
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

            output[[paste0("dl_", prefix, "_expression_matrix")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_expression_matrix_", cross_expression_export_effective_mode(), ".csv")
                },
                content = function(file) {
                    if (!length(selected_source_genes())) {
                        stop("Select one or more genes before downloading an expression matrix.", call. = FALSE)
                    }

                    export_tbl <- build_expression_export_table(
                        obj = cross_object(),
                        feature_map = cross_expression_export_feature_map(),
                        mode = cross_expression_export_effective_mode()
                    )

                    export_tbl <- add_export_provenance_columns(
                        export_tbl,
                        tab_label = cross_integration_label(cross_key),
                        integration_label = cross_integration_label(cross_key),
                        extra = list(
                            atlas_export_type = "expression_matrix",
                            atlas_matrix_mode = cross_expression_export_effective_mode()
                        )
                    )

                    utils::write.csv(export_tbl, file, row.names = FALSE, na = "")
                }
            )

            output[[paste0("dl_", prefix, "_umap")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_umap.", get_ext())
                },
                content = function(file) {
                    panel_specs <- cross_comparison_panel_specs(cross_resolution(), cross_key)
                    panel_n <- nrow(panel_specs)
                    block_cols <- min(2L, max(1L, as.integer(input[[paste0(prefix, "_umap_columns")]] %||% 1L)))
                    row_units <- cross_comparison_height_units(panel_specs)

                    if (length(row_units)) {
                        row_units <- vapply(split(row_units, ceiling(seq_along(row_units) / block_cols)), max, integer(1))
                    }

                    plot_height <- if (length(row_units)) {
                        max(6, 1.8 + sum(row_units) * 3.25)
                    } else {
                        6
                    }
                    plot_width <- if (block_cols == 1L) 17 else 30

                    save_ggplot(
                        file = file,
                        plot_obj = cross_umap_plot_obj(),
                        width = plot_width,
                        height = plot_height,
                        tab_label = cross_integration_label(cross_key),
                        integration_label = cross_integration_label(cross_key),
                        extra = list(atlas_export_type = "umap")
                    )
                }
            )

            output[[paste0("dl_", prefix, "_ortholog_trace")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_ortholog_trace.", get_ext())
                },
                content = function(file) {
                    gene_n <- length(ortholog_trace_source_genes())

                    save_ggplot(
                        file = file,
                        plot_obj = ortholog_trace_plot_obj(),
                        width = 14,
                        height = max(7, ortholog_trace_height_px(gene_n) / 95),
                        tab_label = cross_integration_label(cross_key),
                        integration_label = cross_integration_label(cross_key),
                        extra = list(atlas_export_type = "ortholog_trace")
                    )
                }
            )

            output[[paste0("dl_", prefix, "_heatmap")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_heatmap.", get_ext())
                },
                content = function(file) {
                    feature_n <- length(cross_resolution()$plot_features)

                    save_ggplot(
                        file = file,
                        plot_obj = cross_heatmap_plot_obj(),
                        width = 13,
                        height = max(5.5, 1.8 + feature_n * 0.55),
                        tab_label = cross_integration_label(cross_key),
                        integration_label = cross_integration_label(cross_key),
                        extra = list(atlas_export_type = "heatmap")
                    )
                }
            )

            output[[paste0("dl_", prefix, "_ridge")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_ridgeplot.", get_ext())
                },
                content = function(file) {
                    feature_n <- length(cross_resolution()$plot_features)
                    group_n <- dplyr::n_distinct(cross_object()@meta.data[[cross_group_by()]])

                    save_ggplot(
                        file = file,
                        plot_obj = cross_ridge_plot_obj(),
                        width = 10,
                        height = max(6, expression_ridge_height_px(feature_n, group_n) / 95),
                        tab_label = cross_integration_label(cross_key),
                        integration_label = cross_integration_label(cross_key),
                        extra = list(atlas_export_type = "ridgeplot")
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
                        height = max(6, feature_n * 0.85 + 2),
                        tab_label = cross_integration_label(cross_key),
                        integration_label = cross_integration_label(cross_key),
                        extra = list(atlas_export_type = "dotplot")
                    )
                }
            )

            # ── Distribution UMAP ───────────────────────────────────────────
            cross_dist_umap_plot_obj <- reactive({
                obj <- cross_object()
                group_by <- cross_dist_group_by()
                pt_size <- as.numeric(input[[paste0(prefix, "_dist_pt_size")]] %||% 0.75)
                plot_obj <- apply_metadata_display_order(obj, group_by)
                colors_use <- distribution_colors_use(plot_obj, group_by)
                show_cluster_labels <- is_cluster_distribution_group(group_by)

                p <- scCustomize::DimPlot_scCustom(
                    seurat_object = plot_obj,
                    colors_use = colors_use,
                    group.by = group_by,
                    pt.size = pt_size,
                    label = show_cluster_labels,
                    repel = show_cluster_labels,
                    raster = TRUE
                ) &
                    app_plot_theme() &
                    theme(
                        legend.title = element_blank(),
                        legend.position = if (show_cluster_labels) "none" else "top",
                        panel.grid = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        axis.line = element_blank(),
                        plot.margin = margin(8, 14, 10, 10)
                    )

                p & labs(title = NULL, color = NULL)
            }) %>% bindCache(
                cross_key,
                cross_dist_group_by(),
                as.numeric(input[[paste0(prefix, "_dist_pt_size")]] %||% 0.75),
                cache = "app"
            )

            output[[paste0(prefix, "_dist_umap_plot")]] <- renderPlot(
                { cross_dist_umap_plot_obj() },
                height = 540,
                res = 110
            )

            output[[paste0("dl_", prefix, "_dist_umap")]] <- downloadHandler(
                filename = function() paste0(prefix, "_dist_umap.", get_ext()),
                content = function(file) {
                    save_ggplot(
                        file = file,
                        plot_obj = cross_dist_umap_plot_obj(),
                        width = 10,
                        height = 8,
                        tab_label = cross_integration_label(cross_key),
                        integration_label = cross_integration_label(cross_key),
                        extra = list(atlas_export_type = "distribution_umap")
                    )
                }
            )

            cross_dist_umap3d_plot_data <- reactive({
                group_by <- cross_dist_group_by()
                obj <- apply_metadata_display_order(cross_object(), group_by)
                embedding <- get_cross_umap3d(cross_key)

                validate(need(!is.null(embedding), "3D UMAP is not available for this integration (no PCA or pre-computed 3D reduction found)."))

                cell_ids <- intersect(rownames(obj@meta.data), rownames(embedding))

                validate(need(length(cell_ids) > 0, "No 3D UMAP coordinates are available for this integration."))

                distribution_df <- tibble(
                    cell_id = cell_ids,
                    group_value = as.character(obj@meta.data[cell_ids, group_by, drop = TRUE]),
                    umap3d_1 = embedding[cell_ids, 1],
                    umap3d_2 = embedding[cell_ids, 2],
                    umap3d_3 = embedding[cell_ids, 3]
                ) %>%
                    filter(
                        !is.na(group_value) & nzchar(group_value),
                        !is.na(umap3d_1), !is.na(umap3d_2), !is.na(umap3d_3)
                    )

                validate(need(nrow(distribution_df) > 0, "No cells are available for the 3D UMAP view."))

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

            output[[paste0(prefix, "_dist_umap3d_plot")]] <- plotly::renderPlotly({
                plot_data <- cross_dist_umap3d_plot_data()
                df <- plot_data$data
                color_map <- plot_data$color_map
                group_by <- plot_data$group_by
                marker_size <- max(1.8, as.numeric(input[[paste0(prefix, "_dist_pt_size")]] %||% 0.75) * 2.2)
                group_levels <- names(color_map)

                p <- plotly::plot_ly(type = "scatter3d", mode = "markers")

                for (group_level in group_levels) {
                    group_df <- df %>% filter(as.character(group_value) == !!group_level)
                    if (!nrow(group_df)) next
                    p <- p %>%
                        plotly::add_trace(
                            data = group_df,
                            x = ~umap3d_1, y = ~umap3d_2, z = ~umap3d_3,
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
                        margin = list(l = 0, r = 180, b = 0, t = 10),
                        legend = list(
                            x = 1.02, xanchor = "left",
                            y = 1, yanchor = "top",
                            font = list(size = 14, color = unname(app_palette["text"])),
                            itemsizing = "constant",
                            bgcolor = "rgba(255,255,255,0.82)",
                            bordercolor = "rgba(201,214,196,0.92)",
                            borderwidth = 1
                        ),
                        scene = list(
                            aspectmode = "data",
                            xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                            yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                            zaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                            camera = list(eye = list(x = 1.22, y = 1.06, z = 0.75))
                        )
                    ) %>%
                    plotly::config(displaylogo = FALSE)
            })

            # ── Cluster composition ─────────────────────────────────────────
            cross_composition_by <- reactive({
                choices <- cross_composition_choices(cross_object())
                resolve_choice(
                    input[[paste0(prefix, "_composition_by")]],
                    choices,
                    default = if (length(choices)) unname(choices[[1]]) else NA_character_
                )
            })

            cross_composition_cluster_by <- reactive({
                clustering_cols <- c("Rank_1st", "Rank_2nd", "Rank_3rd", "cluster_label", "species_cell_class")
                available <- clustering_cols[clustering_cols %in% colnames(cross_object()@meta.data)]
                preferred <- cross_dist_group_by()
                if (length(preferred) && preferred %in% available) return(preferred)
                if (length(available)) return(available[[1]])
                NA_character_
            })

            cross_marker_dataset_key <- reactive(cross_key)

            cross_markers_full <- reactive({
                read_cluster_markers_cache(cross_marker_dataset_key(), top_n = FALSE)
            })

            cross_marker_source <- reactive({
                marker_tbl <- cross_markers_full()

                if (is.null(marker_tbl) || !nrow(marker_tbl)) {
                    return(NA_character_)
                }

                available_sources <- unique(as.character(marker_tbl$cluster_source))
                preferred_source <- cross_dist_group_by()
                clustering_sources <- c("Rank_1st", "Rank_2nd", "Rank_3rd", "Rank_4th", "Rank_5th", "cluster_label", "species_cell_class")
                if (!(length(preferred_source) && preferred_source %in% clustering_sources)) {
                    preferred_source <- NULL
                }

                select_marker_cluster_source(
                    available_sources = available_sources,
                    preferred_source = preferred_source,
                    default_source = "Rank_1st"
                )
            })

            cross_marker_lookup <- reactive({
                cluster_label_lookup(cross_object(), cross_marker_source())
            })

            cross_marker_choices <- reactive({
                marker_tbl <- cross_markers_full()

                if (is.null(marker_tbl) || !nrow(marker_tbl)) {
                    return(character(0))
                }

                marker_source <- cross_marker_source()

                marker_tbl <- marker_tbl %>%
                    filter(cluster_source == !!marker_source)

                available_clusters <- unique(as.character(marker_tbl$cluster))
                lookup_tbl <- cross_marker_lookup() %>%
                    filter(cluster %in% available_clusters)

                if (!nrow(lookup_tbl)) {
                    return(setNames(cluster_value_levels(available_clusters), cluster_value_levels(available_clusters)))
                }

                setNames(lookup_tbl$cluster, lookup_tbl$choice_label)
            })

            output[[paste0(prefix, "_marker_cluster_ui")]] <- renderUI({
                choices <- cross_marker_choices()

                if (!length(choices)) {
                    return(
                        tags$p(
                            class = "marker-status-hint",
                            tagList(
                                "Marker cache missing for this integration. Run ",
                                tags$code("Rscript scripts/build_cluster_markers_cache.R"),
                                "."
                            )
                        )
                    )
                }

                selectInput(
                    inputId = paste0(prefix, "_marker_cluster"),
                    label = "Cluster",
                    choices = choices,
                    selected = resolve_choice(
                        input[[paste0(prefix, "_marker_cluster")]],
                        choices,
                        default = unname(choices[[1]])
                    )
                )
            })

            observeEvent(cross_marker_choices(), {
                pending_cluster <- pending_url_cluster()

                if (is.null(pending_cluster) || !nzchar(pending_cluster) || !identical(input$main_tabs %||% "overview", prefix)) {
                    return()
                }

                choices <- cross_marker_choices()

                if (!length(choices)) {
                    return()
                }

                selected_cluster <- resolve_choice(
                    pending_cluster,
                    choices,
                    default = unname(choices[[1]])
                )

                updateSelectInput(
                    session,
                    inputId = paste0(prefix, "_marker_cluster"),
                    selected = selected_cluster
                )
                pending_url_cluster(NULL)
            }, ignoreInit = FALSE)

            cross_marker_cluster <- reactive({
                choices <- cross_marker_choices()
                default_choice <- if (length(choices)) unname(choices[[1]]) else NULL

                resolve_choice(
                    input[[paste0(prefix, "_marker_cluster")]],
                    choices,
                    default = default_choice
                )
            })

            cross_marker_top_n <- reactive({
                marker_n <- suppressWarnings(as.integer(input[[paste0(prefix, "_marker_top_n")]] %||% 10L))

                if (is.na(marker_n) || marker_n < 1L) {
                    marker_n <- 10L
                }

                min(marker_n, 25L)
            })

            cross_marker_table_raw <- reactive({
                marker_tbl <- cross_markers_full()

                validate(
                    need(
                        !is.null(marker_tbl),
                        tagList(
                            "Cluster markers are not cached for this integration. Run ",
                            tags$code("Rscript scripts/build_cluster_markers_cache.R"),
                            "."
                        )
                    ),
                    need(nrow(marker_tbl) > 0, "No cluster markers are available for this integration.")
                )

                cluster_id <- cross_marker_cluster()
                marker_source <- cross_marker_source()
                filtered_tbl <- marker_tbl %>%
                    filter(cluster_source == !!marker_source) %>%
                    filter(cluster == !!cluster_id) %>%
                    arrange(p_val_adj, desc(avg_log2FC), desc(pct.1), gene)

                validate(need(nrow(filtered_tbl) > 0, "No marker rows are available for the selected cluster."))

                filtered_tbl
            })

            output[[paste0(prefix, "_markers_status_ui")]] <- renderUI({
                marker_tbl <- cross_markers_full()
                busy_message <- marker_job_messages()[[prefix]]

                if (is.null(marker_tbl) || !nrow(marker_tbl)) {
                    return(NULL)
                }

                if (!is.null(busy_message) && nzchar(busy_message)) {
                    return(tags$p(class = "marker-status-hint is-working", busy_message))
                }

                source_species <- input$source_species %||% "medicago"
                cluster_choices <- cross_marker_choices()
                marker_source <- cross_marker_source()
                active_group_by <- cross_dist_group_by()
                current_cluster_label <- names(cluster_choices)[match(cross_marker_cluster(), unname(cluster_choices))]
                current_cluster_label <- current_cluster_label[!is.na(current_cluster_label) & nzchar(current_cluster_label)][1] %||% cross_marker_cluster()
                source_sentence <- if (marker_cluster_source_matches(marker_source, active_group_by)) {
                    paste0("using ", marker_cluster_source_label(marker_source), ".")
                } else {
                    paste0(
                        "using ",
                        marker_cluster_source_label(marker_source),
                        " because ",
                        metadata_column_label(active_group_by),
                        " is not a cached clustering solution."
                    )
                }
                tags$p(
                    class = "marker-status-hint",
                    sprintf(
                        "Showing markers for %s %s Top markers will be mapped back into the %s source-gene panel before they are added.",
                        current_cluster_label,
                        source_sentence,
                        species_label(source_species)
                    )
                )
            })

            output[[paste0(prefix, "_markers_table")]] <- DT::renderDT({
                marker_tbl <- cross_marker_table_raw()

                display_tbl <- marker_tbl %>%
                    transmute(
                        gene = display_cross_feature_labels(cross_key, gene),
                        avg_log2FC = avg_log2FC,
                        pct.1 = pct.1,
                        pct.2 = pct.2,
                        p_val_adj = p_val_adj
                    )

                DT::datatable(
                    display_tbl,
                    rownames = FALSE,
                    escape = TRUE,
                    selection = "none",
                    class = "stripe hover order-column compact",
                    options = list(
                        dom = "tip",
                        pageLength = 10,
                        autoWidth = TRUE,
                        order = list(list(1, "desc"))
                    )
                ) %>%
                    DT::formatRound(columns = c("avg_log2FC", "pct.1", "pct.2"), digits = 3) %>%
                    DT::formatSignif(columns = "p_val_adj", digits = 3)
            })

            observeEvent(input[[paste0(prefix, "_add_markers")]], {
                marker_tbl <- cross_marker_table_raw() %>%
                    slice_head(n = cross_marker_top_n())

                source_species <- input$source_species %||% "medicago"
                set_marker_job_message(
                    prefix,
                    sprintf(
                        "Adding the top %d markers from the %s cluster in %s into the shared %s source-gene panel...",
                        nrow(marker_tbl),
                        cross_marker_cluster(),
                        cross_integration_label(cross_key),
                        species_label(source_species)
                    )
                )
                session$sendCustomMessage("atlas_button_busy", list(
                    button_id = paste0(prefix, "_add_markers"),
                    label = "Adding..."
                ))
                candidate_genes <- map_cross_marker_features_to_source_genes(
                    cross_key = cross_key,
                    feature_ids = marker_tbl$gene,
                    source_species = source_species
                )

                valid_source_ids <- build_gene_choices(
                    source_species,
                    source_integration()
                )$feature_ids %||% character(0)

                genes_to_add <- unique(candidate_genes[candidate_genes %in% valid_source_ids])

                if (!length(genes_to_add)) {
                    showNotification(
                        sprintf(
                            "No top markers from this %s cluster could be mapped into the %s source-gene panel.",
                            cross_integration_label(cross_key),
                            species_label(source_species)
                        ),
                        type = "warning",
                        duration = 8
                    )
                    return()
                }

                updated_selection <- unique(c(staged_source_genes(), genes_to_add))

                update_selected_genes_input(
                    choice_bundle = build_gene_choices(
                        source_species,
                        source_integration()
                    ),
                    selected = updated_selection
                )
                applied_selected_genes(updated_selection)

                skipped_n <- nrow(marker_tbl) - length(genes_to_add)
                showNotification(
                    paste0(
                        "Added ",
                        length(genes_to_add),
                        " mapped marker gene(s) to the shared panel",
                        if (skipped_n > 0) paste0("; ", skipped_n, " could not be mapped into the current source species.") else "."
                    ),
                    type = "message",
                    duration = 6
                )
            }, ignoreInit = TRUE)

            output[[paste0("dl_", prefix, "_markers")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_cluster_", cross_marker_cluster(), "_markers.csv")
                },
                content = function(file) {
                    marker_tbl <- cross_marker_table_raw()
                    lookup_tbl <- cross_marker_lookup()
                    cluster_labels <- lookup_tbl$cluster_label[match(marker_tbl$cluster, lookup_tbl$cluster)]

                    export_tbl <- marker_tbl %>%
                        mutate(
                            cluster_label = ifelse(is.na(cluster_labels) | !nzchar(cluster_labels), cluster, cluster_labels),
                            gene_label = display_cross_feature_labels(cross_key, gene)
                        ) %>%
                        select(cluster, cluster_label, gene, gene_label, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
                        add_export_provenance_columns(
                            tab_label = cross_integration_label(cross_key),
                            integration_label = cross_integration_label(cross_key),
                            extra = list(atlas_export_type = "cluster_markers")
                        )

                    write.csv(export_tbl, file = file, row.names = FALSE, na = "")
                }
            )

            cross_composition_plot_obj <- reactive({
                obj <- cross_object()
                composition_by <- cross_composition_by()
                cluster_by <- cross_composition_cluster_by()
                plot_obj <- apply_metadata_display_order(obj, composition_by)

                validate(
                    need(!is.null(composition_by) && nzchar(composition_by %||% ""), "No composition metadata available for this integration."),
                    need(!is.na(cluster_by) && nzchar(cluster_by), "No clustering metadata available for this integration.")
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

                validate(need(nrow(composition_df) > 0, "No cluster composition data available."))

                composition_df <- composition_df %>%
                    mutate(
                        cluster = factor(cluster, levels = cluster_value_levels(cluster)),
                        composition = order_metadata_values(composition, composition_by)
                    )

                fill_values <- composition_colors_use(composition_df$composition, composition_by)
                legend_rows <- max(1L, min(3L, ceiling(dplyr::n_distinct(composition_df$composition) / 8L)))

                ggplot(composition_df, aes(x = cluster, y = cell_count, fill = composition)) +
                    geom_col(position = "fill", width = 0.82, colour = "#FFFFFF", linewidth = 0.18) +
                    scale_y_continuous(
                        labels = scales::label_percent(accuracy = 1),
                        expand = expansion(mult = c(0, 0.02))
                    ) +
                    {if (!is.null(fill_values)) scale_fill_manual(values = fill_values)} +
                    guides(fill = guide_legend(nrow = legend_rows, byrow = TRUE)) +
                    labs(x = metadata_column_label(cluster_by), y = "Percent of cells", fill = NULL) +
                    app_plot_theme() +
                    theme(
                        panel.grid.major.x = element_blank(),
                        axis.text.x = element_text(size = 11, angle = 35, hjust = 1),
                        legend.position = "top",
                        legend.text = element_text(size = 9),
                        plot.margin = margin(8, 12, 8, 10)
                    )
            }) %>% bindCache(
                cross_key,
                cross_composition_by(),
                cross_composition_cluster_by(),
                cache = "app"
            )

            output[[paste0(prefix, "_composition_plot")]] <- renderPlot(
                { cross_composition_plot_obj() },
                height = function() {
                    composition_by <- tryCatch(cross_composition_by(), error = function(e) NULL)
                    obj <- tryCatch(cross_object(), error = function(e) NULL)
                    if (is.null(obj) || is.null(composition_by) || !nzchar(composition_by %||% "")) return(460)
                    level_count <- tryCatch({
                        vals <- as.character(obj@meta.data[[composition_by]])
                        length(unique(vals[!is.na(vals) & nzchar(vals)]))
                    }, error = function(e) 1L)
                    legend_rows <- max(1L, min(3L, ceiling(level_count / 8L)))
                    max(460, 380 + legend_rows * 44)
                },
                res = 110
            )
        })
    }

    walk(cross_integration_keys, register_cross_tab)
}

shinyApp(ui = ui, server = server)
