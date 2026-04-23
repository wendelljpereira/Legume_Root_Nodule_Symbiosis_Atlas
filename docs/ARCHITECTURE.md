# Architecture Notes

The Shiny entrypoint is `app.R`. Reusable foundations are split into `R/`:

- `R/atlas_core.R`: package imports, registries, data loading, gene mapping, plotting helpers, cache utilities, and metadata helpers.
- `R/atlas_annotations.R`: CSV-backed cluster annotation overlays used by the app, cache builders, template generator, and release validator.
- `R/atlas_ui_tabs.R`: within-species and cross-species tab builders.
- `app.R`: app-level JavaScript, global UI shell, server wiring, and final `shinyApp()` call.

Cluster annotations are intentionally kept outside the Seurat RDS files. The
runtime loads `annotations/cluster_annotations/<dataset_key>.csv`, validates
that each `cluster_source`/`cluster_id` exists in the loaded object, and adds
display-only label columns such as `cluster_label` and `Rank_1st_label`.
Raw cluster IDs remain in Seurat identities and `Rank_*` metadata columns.

The current refactor is intentionally incremental. The next safe extraction target is server modules for within-species tabs and cross-species tabs, using Shiny modules with namespaced input/output IDs.
