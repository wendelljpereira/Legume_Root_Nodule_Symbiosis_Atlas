# Legume Root Nodule Symbiosis Atlas Documentation

```{raw} html
<div class="hero-card">
<strong>Explore root nodule symbiosis across species.</strong><br>
This Shiny atlas helps researchers inspect gene expression across published single-cell legume nodule datasets, compare ortholog-aware expression patterns, and export reproducible figures with atlas provenance.
</div>
```

![Atlas overview page](_static/screenshots/overview.png){.doc-screenshot}

## What You Can Do

- Search genes from *Medicago truncatula*, *Glycine max*, or *Lotus japonicus*.
- View expression in within-species single-cell atlases.
- Compare mapped ortholog features in CAMEx and SATURN cross-species integrations.
- Use cluster markers to seed a gene panel.
- Inspect average expression by cluster, dot plots, UMAP expression, violin plots, and ridge plots.
- Download figures and mapping tables for records, teaching, or manuscript preparation.

## Who This Is For

This documentation is written for plant biologists, single-cell users, and computational researchers who want to explore root nodule symbiosis expression patterns without rebuilding the atlas from raw data.

The developer and release sections are for people deploying the app, refreshing data, or replacing project-level annotation/orthology tables.

```{toctree}
:maxdepth: 2
:caption: Using The App

quickstart
user-guide/overview
user-guide/gene-expression
user-guide/cross-species
user-guide/plots
user-guide/markers-annotations
user-guide/exports-citation
strengths-limitations
```

```{toctree}
:maxdepth: 2
:caption: Administration

admin/customization
admin/deployment
RELEASE_CHECKLIST
```

```{toctree}
:maxdepth: 2
:caption: Development

developer/contributing
ARCHITECTURE
```

## Citation

The atlas accompanies an in-preparation manuscript. Until publication, cite the app as:

> Pereira W. et al. A cross-species single-cell atlas of legume root nodule symbiosis. In preparation. Legume Root Nodule Symbiosis Atlas, pre-publication release.

Please also cite the original datasets and computational methods used in the atlas when interpreting or reusing results.
