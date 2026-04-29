# Legume Root Nodule Symbiosis Atlas

Interactive Shiny atlas for exploring single-cell root nodule symbiosis datasets across *Medicago truncatula*, *Glycine max*, and *Lotus japonicus*.

The app lets researchers search genes of interest, inspect expression across published single-cell atlases, compare ortholog-aware expression patterns across species, review cluster markers, and export publication-ready figures with atlas provenance.


## What The App Provides

- Within-species expression exploration for *Medicago truncatula*, *Glycine max*, and *Lotus japonicus*.
- Cross-species comparison through CAMEx and SATURN integrations.
- Gene search using gene IDs and available annotations/synonyms.
- Expression UMAPs, violin plots, averaged expression by cluster, and dot plots.
- Cluster marker tables and marker-to-gene-panel workflows.
- Ortholog tracing for one-to-one, one-to-many, and missing mappings.
- CSV-backed cluster annotations that can be edited without modifying Seurat RDS files.
- Docker/local execution for reproducible use and deployment testing.

## Documentation

User-facing documentation is in `docs/` and is configured for ReadTheDocs.

- Start here: `docs/index.md`
- Local build requirements: `docs/requirements.txt`
- ReadTheDocs config: `.readthedocs.yaml`
- Release checklist: `docs/RELEASE_CHECKLIST.md`
- Architecture notes: `docs/ARCHITECTURE.md`

The documentation includes screenshots from the app, quick-start instructions, usage walkthroughs, strengths and limitations, deployment notes, and customization instructions for orthology and cluster annotation tables.

## Repository Layout

```text
app.R                                  Shiny entrypoint and server wiring
R/                                     reusable app helpers, plotting, annotations, UI builders
www/                                   app CSS and browser-side assets
annotations/                          gene annotations and editable cluster annotation CSVs
metadata/                             generated app metadata caches
orthogroups/                          ortholog/orthogroup table
within_species_integrated_datasets/   within-species Seurat objects and slim app objects
app_ready_integration/                 cross-species integration objects and assets
scripts/                              cache builders, validators, smoke tests, release checks
tests/                                testthat and browser-test scaffolding
docs/                                 ReadTheDocs/Sphinx documentation
```

Large atlas objects are intentionally excluded from normal GitHub tracking unless you use a large-file strategy such as Git LFS or external release artifacts.

## Running Locally

### Docker, Recommended

```bash
docker compose up --build
```

Open <http://localhost:3838>.

The included `docker-compose.yml` mounts atlas data, metadata, orthology, and annotation folders into the container as read-only volumes. You can override settings in a local `.env` file or shell environment.

Manual equivalent:

```bash
docker build -t legume-atlas .

docker run --rm -p 3838:3838 \
    -v "$(pwd)/within_species_integrated_datasets:/srv/atlas/within_species_integrated_datasets:ro" \
    -v "$(pwd)/app_ready_integration:/srv/atlas/app_ready_integration:ro" \
    -v "$(pwd)/orthogroups:/srv/atlas/orthogroups:ro" \
    -v "$(pwd)/annotations:/srv/atlas/annotations:ro" \
    -v "$(pwd)/metadata:/srv/atlas/metadata:ro" \
    -v "$(pwd)/celltype_overrides:/srv/atlas/celltype_overrides:ro" \
    legume-atlas
```

### Directly In R

```r
install.packages(c(
    "shiny", "shinyWidgets", "DT", "dplyr", "tidyr", "purrr", "tibble",
    "ggplot2", "patchwork", "svglite", "plotly", "uwot", "viridisLite",
    "Matrix", "rlang", "scales", "Seurat", "scCustomize",
    "renv", "testthat", "shinytest2"
))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("ComplexHeatmap", update = FALSE, ask = FALSE)
shiny::runApp()
```

## Editable Scientific Inputs

The app reads expression data from Seurat RDS files, but the main scientific tables that users may want to replace are plain files.

| Setting | Default | Purpose |
|---|---|---|
| `ATLAS_ORTHOLOGS_PATH` | `orthogroups/joint_orthogroups.tsv` | Ortholog/orthogroup table. TSV or CSV is auto-detected. |
| `ATLAS_CLUSTER_ANNOTATIONS_DIR` | `annotations/cluster_annotations/` | Preferred folder for shipped or user-edited cluster annotation CSVs. |
| `ATLAS_CELLTYPE_OVERRIDES_DIR` | `celltype_overrides/` | Legacy simple overrides with `cluster_id,label`; still supported. |
| `ATLAS_SAMPLE_NAME_OVERRIDES_PATH` | `metadata/sample_name_overrides.tsv` | Optional table for replacing source sample IDs with atlas display names. |
| `ATLAS_ACCESS_PASSWORD` | unset | Enables an optional access gate for private deployments. |
| `ATLAS_VERSION` | `0.1.0-preview` | Displayed in the app footer and exports. |
| `ATLAS_LAST_UPDATED` | current date | Displayed in the app footer. |
| `ATLAS_CITATION` | atlas citation text | Displayed in the citation modal and export provenance. |

Cluster annotation CSVs use this schema:

```csv
dataset_key,cluster_source,cluster_id,cell_type_label,annotation_status,confidence_score,marker_evidence,reference_source,notes
medicago_Seurat,Rank_1st,0,Cortex,curated,0.95,"MtN21; MtN24","Original paper + marker review",""
```

Generate editable templates from the current datasets:

```bash
Rscript scripts/write_cluster_annotation_templates.R
```

After editing labels, validate and refresh app metadata:

```bash
Rscript scripts/validate_cluster_annotations.R
Rscript scripts/refresh_app_metadata.R
```

## Updating Data And Caches

When source Seurat objects change, rebuild slim app objects and generated caches before release.

```bash
Rscript scripts/build_within_app_slim.R
Rscript scripts/add_umap3d_to_within_species_slim.R --overwrite
Rscript scripts/build_cross_app_slim.R
Rscript scripts/build_gene_catalog_cache.R
Rscript scripts/build_cluster_markers_cache.R
Rscript scripts/write_cluster_annotation_templates.R
Rscript scripts/validate_cluster_annotations.R
Rscript scripts/build_atlas_summary_cache.R
Rscript scripts/build_startup_ui_cache.R
Rscript scripts/smoke_test.R
```

The `add_umap3d_to_within_species_slim.R` step precomputes 3D UMAP coordinates from the PCA-like reductions already stored in the full within-species objects and writes them into the slim app files. Keep this out of app startup so the browser only reads static coordinates.

For the strict pre-release gate:

```bash
Rscript scripts/release_check.R --skip-browser
```

Then launch the app locally and run the Firefox browser smoke test, if Node/npm dependencies are available:

```bash
npm install
ATLAS_APP_URL=http://127.0.0.1:3838 npm run test:e2e:firefox
```

## Deployment Notes

The first deployment target is ShinyApps.io. The full local bundle is currently close to the practical 5 GB threshold, so release checks include a size estimate:

```bash
Rscript scripts/check_bundle_size.R
```

If future releases exceed the bundle limit, use versioned public-read data artifacts or another external storage layer for large RDS/cache files while keeping the Shiny interface unchanged.

Docker or a VM/Posit Connect deployment remains the safest option for full-data public hosting if memory or bundle-size limits become restrictive.

## Strengths And Limitations

Strengths:

- Integrates multiple published legume root nodule symbiosis single-cell datasets.
- Supports both within-species and cross-species gene expression exploration.
- Makes ortholog mapping visible rather than hiding one-to-many or missing mappings.
- Keeps user-editable annotations and orthology inputs outside compiled app code.

Limitations:

- Cross-species plots depend on orthogroup membership and feature availability; they are comparative views, not proof of one-to-one conserved cell states.
- Cluster annotations are curated labels over clustering solutions, not immutable biological truth.
- Very large gene panels and large orthogroup expansions are capped to keep plots responsive.
- Public hosting depends on ShinyApps.io bundle and memory constraints unless external data hosting is configured.

## Citation

Cite the atlas version used for your analysis and the original datasets and methods that support the specific results you interpret.

> Pereira W. et al. A cross-species single-cell atlas of legume root nodule symbiosis. Legume Root Nodule Symbiosis Atlas, version 1.0.

## License

See `LICENSE` for the current terms.
