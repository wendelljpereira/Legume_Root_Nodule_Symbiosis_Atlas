# Legume Root Nodule Symbiosis Atlas

Interactive single-cell atlas of root nodule symbiosis across
*Medicago truncatula*, *Glycine max*, and *Lotus japonicus*, with
two within-species integrations (ComBat/BBKNN, Seurat) and two cross-species
integrations (Camex, SATURN).

The app accompanies an in-preparation publication. Until the paper is
published, the hosted deployment on [ShinyApps.io](https://www.shinyapps.io/)
is password-gated.

## Running locally

### Option 1: Docker (recommended)

```bash
docker compose up --build
```

Open http://localhost:3838.

The included `docker-compose.yml` mounts the local atlas data folders and exposes the same container defaults as the `docker run` example below. Override any of `ATLAS_ACCESS_PASSWORD`, `ATLAS_VERSION`, `ATLAS_LAST_UPDATED`, `ATLAS_ORTHOLOGS_PATH`, or `ATLAS_CELLTYPE_OVERRIDES_DIR` in your shell or a local `.env` file before starting Compose.

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

### Option 2: Run directly in R

```r
install.packages(c(
    "shiny", "shinyWidgets", "dplyr", "tidyr", "purrr", "tibble",
    "ggplot2", "patchwork", "svglite", "plotly", "uwot", "viridisLite",
    "Seurat", "scCustomize"
))
shiny::runApp()
```

If the CAMEx or SATURN `.rds` files are refreshed, regenerate their cached 3D
embeddings with:

```bash
Rscript scripts/build_cross_app_slim.R
Rscript scripts/build_cross_umap3d_embeddings.R
```

If any within-species `.rds` files are refreshed, regenerate their slim deploy
artifacts with:

```bash
Rscript scripts/build_within_app_slim.R
```

## User-supplied data overrides

The app reads its heavy data (Seurat `.rds` files) from fixed folders, but
the two scientifically-editable inputs can be swapped via environment
variables without modifying `app.R`:

| Env var | Default | Purpose |
|---|---|---|
| `ATLAS_ORTHOLOGS_PATH` | `orthogroups/joint_orthogroups.tsv` | Ortholog/orthogroup table. Must have columns `Orthogroup`, `medicago.fa`, `glycine.fa`, `lotus.fa` (same format as OrthoFinder output). TSV or CSV is auto-detected by extension. |
| `ATLAS_CELLTYPE_OVERRIDES_DIR` | `celltype_overrides/` | Folder of per-dataset CSVs. Any file named `<species>_<integration>.csv` or `<cross_key>.csv` (e.g. `medicago_Seurat.csv`, `camex.csv`) is read at startup. The CSV must have columns `cluster_id` and `label`; matching cluster IDs are relabeled. |

Example `celltype_overrides/medicago_Seurat.csv`:

```csv
cluster_id,label
0,Cortex
1,Infected cells
2,Nodule primordium
```

## Other environment variables

| Env var | Default | Purpose |
|---|---|---|
| `ATLAS_ACCESS_PASSWORD` | *(unset)* | If set, app renders a password gate at startup. Leave unset during local Docker use. |
| `ATLAS_VERSION` | `0.1.0-preview` | Shown in footer + citation modal. |
| `ATLAS_LAST_UPDATED` | today's date | Shown as "Data updated" in footer. |
| `ATLAS_CITATION` | preview citation | Text shown inside the **Cite this atlas** modal. |

## ShinyApps.io deployment

1. `ATLAS_ACCESS_PASSWORD` — set under **Advanced > Environment variables** in
   the ShinyApps.io dashboard.
2. `ATLAS_VERSION`, `ATLAS_LAST_UPDATED` — set whenever you push a new data
   snapshot.
3. Large Seurat `.rds` files live in the bundle — make sure the bundle stays
   under ShinyApps.io size limits.

## Before deploying

Run the slim-cache builders and smoke test from the main project directory:

```bash
Rscript scripts/build_within_app_slim.R
Rscript scripts/build_cross_app_slim.R
Rscript scripts/build_gene_catalog_cache.R
Rscript scripts/build_cluster_markers_cache.R
Rscript scripts/build_atlas_summary_cache.R
Rscript scripts/smoke_test.R
```

`Rscript scripts/smoke_test.R` exits non-zero if any required slim RDS file is
missing/corrupt or if a dataset lacks the `data` layer / expected reductions.

## Browser console note

The app should load without 404s or JavaScript errors. In Chromium-based
browsers you may still see a `Canvas2D ... willReadFrequently` warning while
Plotly/canvas export paths render; that warning comes from the browser canvas
stack and is harmless.
