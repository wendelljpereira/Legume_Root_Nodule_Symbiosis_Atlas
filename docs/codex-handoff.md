# Codex Handoff

Last updated: 2026-04-28

## Project

This repository contains the Shiny app for the Legume Root Nodule Symbiosis Atlas: a single-cell RNA-seq atlas for root nodule symbiosis across `Medicago truncatula`, `Glycine max`, and `Lotus japonicus`.

The app is meant for scientific users who want to:

- Inspect within-species expression and cluster structure.
- Seed gene-expression panels from gene search, pasted/imported lists, or cluster markers.
- Compare ortholog-aware expression across species using CAMEX and SATURN integrated embeddings.
- Export plots with provenance for manuscripts and follow-up analyses.

## Current UX Model

- The Overview tab is a light orientation page.
- Each within-species tab is independent. Local gene panels, marker additions, and expression plots should not leak into other species tabs.
- Cross-species tabs are separate from the within-species tabs. They use a shared cross-species Gene expression panel that can stage genes from all three species.
- Orthology-aware comparisons belong in CAMEX/SATURN, not in the within-species tabs.
- Plot rows in within-species tabs appear only after the user clicks `Generate the expression plots`.

## Recent Changes

Most recent commits to inspect:

- `aeaee9d Refine cross-species marker workflow`
- `cb92cd4 Refine atlas expression workflow`
- `7daf832 Prepare atlas app for release documentation`

Key changes in the current work session:

- Reworked within-species tabs so gene panels are local and independent.
- Restored the gene picker "select all matches" bulk option.
- Moved `Split feature UMAP by` into the Gene expression panel controls.
- Removed empty placeholder plot rows until users generate expression plots.
- Changed heatmaps to row-scaled average expression z-scores across clusters.
- Removed all expression ridge plots.
- Refreshed the CSS palette with restrained botanical contrast.
- Simplified summary tiles and removed low-value informational boxes.
- Restricted CAMEX/SATURN Cell distribution UMAP `Color cells by` choices to `Species` and `Clustering opt 1` through `Clustering opt 5`.
- Added CAMEX/SATURN Cell distribution UMAP `Split UMAP by` with default `No split` and optional `Species`.
- Suppressed the integrated 3D UMAP when the distribution UMAP is split by species.
- Moved the cross-species gene selection panel into the `Gene expression` section after `Cluster markers`.
- Reframed the cross-species selector as the shared `Gene expression panel`.
- Added a `Marker feature species` filter for integrated marker tables.
- Added marker feature species to marker CSV exports when it can be inferred.
- Removed `Condition` from within-species distribution/composition selectors; use explicit `Time point` and `Samples` choices instead.
- Added a standalone precompute script for within-species 3D UMAP reductions so slim app files carry static `umap3d` coordinates and the app never computes them during startup or rendering.
- Harmonized Lotus SATURN root/mock metadata into the `Roots` time point when that information is present in the dataset metadata.
- Expanded the SATURN shared Gene expression panel to full width and increased the cross-species heatmap render height for small gene panels.

## Scientific Notes

Integrated per-cluster markers are useful, but they should be interpreted carefully.

- They are appropriate for annotating integrated clusters and detecting shared or enriched programs.
- They can also reflect species imbalance, integration artifacts, or feature-space differences.
- CAMEX stores a shared Medicago-space feature representation, so marker feature species filtering may often be less informative there.
- SATURN uses species-prefixed integrated features, so the species filter should be more informative there.
- A marker feature species filter is not the same as filtering cells by species. If users ask for species-specific markers within integrated clusters, that is a different analysis and should probably use species-stratified differential expression.

## Important Files

- `app.R`: main server logic and cross-tab wiring.
- `R/atlas_core.R`: registries, plotting helpers, cache readers, orthology mapping helpers, marker helpers.
- `R/atlas_ui_tabs.R`: tab-level UI builders for overview, within-species tabs, and cross-species tabs.
- `www/styles.css`: app visual design and layout styles.
- `metadata/cluster_markers/`: expected cache location for precomputed marker tables.
- `scripts/refresh_app_metadata.R`: refresh static metadata/cache tables after dataset changes.
- `scripts/build_gene_catalog_cache.R`: rebuild gene picker catalogs.
- `scripts/build_atlas_summary_cache.R`: rebuild summary cards.
- `scripts/build_cluster_markers_cache.R`: generate marker caches if needed, though the user may provide marker CSVs instead.
- `scripts/build_startup_ui_cache.R`: rebuild cached UI choices/summaries.
- `scripts/add_umap3d_to_within_species_slim.R`: precompute 3D UMAP coordinates from PCA-like reductions in full within-species objects and write them to slim app files.

## How To Launch And Validate

Source-check the app:

```bash
Rscript -e 'source("app.R", local = new.env())'
```

After changing full within-species objects, rebuild static app data before launch:

```bash
Rscript scripts/build_within_app_slim.R
Rscript scripts/add_umap3d_to_within_species_slim.R --overwrite
Rscript scripts/build_gene_catalog_cache.R
Rscript scripts/build_atlas_summary_cache.R
Rscript scripts/build_startup_ui_cache.R
```

Hard-restart local Shiny on the standard port:

```bash
lsof -ti tcp:3803 | xargs -r kill -9
Rscript -e "shiny::runApp('.', host='127.0.0.1', port=3803, launch.browser=FALSE)"
```

Open:

```text
http://127.0.0.1:3803/
```

The user prefers that after each UI fix the app is hard-restarted and opened in Firefox.

## Current Open Questions

- Continue reviewing CAMEX/SATURN cross-species flow and visual hierarchy.
- Decide whether species-specific marker filtering should remain a simple feature-species filter or become a true species-stratified marker analysis.
- Consider making the marker-to-gene-panel workflow clearer, especially when markers cannot be mapped back to the selected source species.
- Inspect cross-species expression panels for legibility when orthogroups expand into many mapped features.
- Re-check whether any URL/permalink state should be updated only on explicit actions, not continuously.

## Local Workspace Notes

At the time this handoff was written, these local items were intentionally not committed:

- `scripts/fix_medicago_feature_ids.R`
If the next session needs them, inspect locally before deciding whether they should be committed. They will not be available on another machine unless they are added and pushed separately.

## Recommended Prompt For Next Codex Session

Start with:

```text
Read docs/codex-handoff.md, inspect the last few commits, launch the app on 127.0.0.1:3803, and help me continue reviewing the cross-species CAMEX and SATURN tabs.
```
