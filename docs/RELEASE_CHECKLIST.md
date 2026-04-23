# Release Checklist

Use this checklist before publishing a new atlas build.

1. Confirm the git worktree is healthy with `git status --short`.
2. Rebuild slim datasets when raw Seurat objects change:
   `Rscript scripts/build_within_app_slim.R` and `Rscript scripts/build_cross_app_slim.R`.
3. Generate or refresh cluster annotation templates when cluster definitions change with `Rscript scripts/write_cluster_annotation_templates.R`; edit `cell_type_label`, status, evidence, and notes locally.
4. Validate shipped annotation CSVs with `Rscript scripts/validate_cluster_annotations.R`.
5. Rebuild derived caches with `Rscript scripts/refresh_app_metadata.R`.
6. Require cache freshness with `Rscript scripts/refresh_app_metadata.R --check-only --fail-on-stale`.
7. Run `Rscript scripts/smoke_test.R` to verify every slim object, assay data layer, and UMAP/3D reduction.
8. Run unit tests with `Rscript -e 'testthat::test_dir("tests/testthat")'`.
9. Confirm deployment target can handle the current data footprint with `Rscript scripts/check_bundle_size.R`.
10. Launch the app and run `ATLAS_APP_URL=http://127.0.0.1:3838 npm run test:e2e:firefox`.
11. Test manually in Firefox: no-gene warning, gene picker, marker staging, Generate buttons, annotation labels, downloads, permalink, and large ortholog fan-out warning.
12. Update `ATLAS_VERSION`, `ATLAS_LAST_UPDATED`, and citation text before deployment.
