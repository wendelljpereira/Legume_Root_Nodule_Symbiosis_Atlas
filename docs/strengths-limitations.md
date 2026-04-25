# Strengths And Limitations

## Strengths

- Integrates multiple published legume root nodule symbiosis single-cell datasets.
- Provides within-species and cross-species expression exploration in one interface.
- Makes ortholog mapping transparent instead of hiding missing or one-to-many mappings.
- Separates editable scientific tables from serialized Seurat objects.
- Supports Docker/local execution for reproducibility and customization.

## Limitations

- Cross-species mappings depend on the orthogroup table and integration feature coverage.
- One-to-many orthogroups can inflate plotted feature counts and complicate interpretation.
- Cluster labels are curated annotations of clustering outputs, not definitive cell-type assignments.
- Sparse single-cell expression may be difficult to interpret in violin plots.
- Public ShinyApps.io hosting may be constrained by bundle size and memory limits.

## Best-Practice Interpretation

Use the atlas as a hypothesis-generation and exploration tool. For publication claims, confirm results against source studies, marker evidence, orthology trace tables, and appropriate downstream analyses.
