# Markers And Cluster Annotations

## Cluster Markers

Marker tables are generated from the current atlas objects and cached for the app. They can be used to inspect cluster-enriched genes or to seed the gene-expression panel.

Marker buttons stage genes only. You still need to click **Generate the expression plots** to render plots.

## Cluster Annotations

Cluster annotations are stored outside the Seurat RDS files in editable CSV tables under `annotations/cluster_annotations/`.

This design keeps the shipped labels transparent and makes it possible for Docker/local users to replace annotations without editing app code or serialized objects.

## Annotation Columns

Required columns are:

```text
dataset_key, cluster_source, cluster_id, cell_type_label,
annotation_status, confidence_score, marker_evidence,
reference_source, notes
```

The app preserves raw cluster IDs and adds display-only label columns such as `cluster_label` or `Rank_1st_label` at runtime.

## Interpreting Annotations

Cell-type labels are curated interpretations of clustering solutions. They should be treated as scientific annotations, not as permanent truth. Marker evidence, source references, and notes should be reviewed when using labels in manuscripts.
