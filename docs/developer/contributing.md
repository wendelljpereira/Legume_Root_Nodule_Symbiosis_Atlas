# Development Notes

## Local Checks

Run parser checks and tests:

```bash
Rscript -e 'testthat::test_dir("tests/testthat", reporter = "summary")'
```

Run the main release gate without browser tests:

```bash
Rscript scripts/release_check.R --skip-browser
```

Run the Firefox E2E test when Node/npm dependencies are available:

```bash
npm install
ATLAS_APP_URL=http://127.0.0.1:3838 npm run test:e2e:firefox
```

## Cache Rebuilds

When raw or slim Seurat objects change, rebuild metadata caches before release:

```bash
Rscript scripts/refresh_app_metadata.R
Rscript scripts/refresh_app_metadata.R --check-only --fail-on-stale
```

## Code Organization

- `app.R` contains the Shiny entrypoint and server wiring.
- `R/atlas_core.R` contains reusable data, plotting, mapping, and cache helpers.
- `R/atlas_annotations.R` contains CSV-backed cluster annotation overlays.
- `R/atlas_ui_tabs.R` contains tab-builder UI helpers.

The next major refactor target is extracting within-species and cross-species server logic into Shiny modules.
