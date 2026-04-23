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
                p("Explore this atlas with native genes and local cluster tools, then promote a local panel into the cross-species comparison workflow when you want to compare orthologs.")
            ),
            div(
                class = "atlas-control-row",
                fluidRow(
                    column(
                        width = 3,
                        div(
                            class = "option-group atlas-control-card",
                            selectInput(
                                inputId = paste0(species_key, "_integration"),
                                label = "Integration",
                                choices = integration_choices,
                                selected = "ComBat_BBKNN"
                            )
                        )
                    ),
                    column(
                        width = 3,
                        div(
                            class = "option-group atlas-control-card",
                            uiOutput(paste0(species_key, "_distribution_group_by_ui"))
                        )
                    ),
                    column(
                        width = 3,
                        div(
                            class = "option-group atlas-control-card",
                            uiOutput(paste0(species_key, "_distribution_split_by_ui"))
                        )
                    ),
                    column(
                        width = 3,
                        div(
                            class = "option-group atlas-control-card atlas-control-card-slider",
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
                p("Review precomputed positive markers for the active clustering solution when marker tables are available for it, add the top hits to this species panel, or download the current cluster as CSV.")
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
                                label = "Add top N genes to the 'Gene expression' panel",
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
                p("Add genes with the picker, or move ahead by using the selected cluster markers to populate this Gene expression panel and generate the plots.")
            ),
            fluidRow(
                column(
                    width = 12,
                    div(
                        class = "option-group gene-panel-card",
                        div(
                            class = "gene-panel-picker",
                            selectizeInput(
                                inputId = paste0(species_key, "_local_selected_genes"),
                                label = paste("Search", label, "genes"),
                                choices = NULL,
                                multiple = TRUE,
                                options = gene_panel_selectize_options(enable_bulk = FALSE)
                            )
                        ),
                        div(
                            class = "gene-action-row",
                            actionButton(
                                inputId = paste0(species_key, "_apply_local_genes"),
                                label = "Generate the expression plots",
                                icon = icon("play"),
                                class = "btn btn-default apply-selection-btn"
                            ),
                            actionButton(
                                inputId = paste0(species_key, "_open_gene_import"),
                                label = "Import list...",
                                icon = icon("file-import"),
                                class = "btn btn-default gene-import-btn"
                            )
                        ),
                        div(
                            class = "selection-meta",
                            `aria-live` = "polite",
                            `aria-atomic` = "true",
                            textOutput(paste0(species_key, "_local_gene_selection_status"))
                        )
                    )
                )
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
                            div(class = "plot-card-title", "Average expression by cluster"),
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
                        `data-permalink-panel` = paste0(species_key, "_dot"),
                        div(
                            class = "plot-card-header",
                            div(class = "plot-card-title", "Multi-gene dot plot"),
                            plot_download_button(paste0("dl_", species_key, "_dot"))
                        ),
                        spinning_plot_output(paste0(species_key, "_dot_plot"), proxy_height = "420px")
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
                        spinning_plot_output(paste0(species_key, "_ridge_plot"), proxy_height = "680px")
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
            uiOutput(paste0(prefix, "_comparison_panel_ui")),
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
                p("Review precomputed positive markers for the active clustering solution when marker tables are available for it, add the top hits to the shared comparison panel, or download the current cluster as CSV.")
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
                                label = "Add top N to comparison panel",
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
                            div(class = "plot-card-title", "Average expression by cluster"),
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
                        div(class = "plot-card-title", "Ortholog mapping summary"),
                        uiOutput(paste0(prefix, "_mapping_table_ui"))
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
                        spinning_plot_output(paste0(prefix, "_ridge_plot"), proxy_height = "680px")
                    )
                )
            )
        )
    )
}
