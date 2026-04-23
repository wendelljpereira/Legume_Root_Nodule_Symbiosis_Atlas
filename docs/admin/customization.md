# Customizing Scientific Inputs

The app is designed so users can customize scientific tables without editing Seurat RDS files or app source code.

## Orthology Table

Set `ATLAS_ORTHOLOGS_PATH` to point to a custom OrthoFinder-style TSV or CSV.

Required columns:

```text
Orthogroup, medicago.fa, glycine.fa, lotus.fa
```

The app auto-detects TSV versus CSV by file extension.

## Cluster Annotation Tables

Set `ATLAS_CLUSTER_ANNOTATIONS_DIR` to a folder containing one CSV per dataset.

Expected filenames include:

```text
medicago_ComBat_BBKNN.csv
medicago_Seurat.csv
glycine_ComBat_BBKNN.csv
glycine_Seurat.csv
lotus_ComBat_BBKNN.csv
lotus_Seurat.csv
camex.csv
saturn.csv
```

Required columns:

```text
dataset_key,cluster_source,cluster_id,cell_type_label,annotation_status,confidence_score,marker_evidence,reference_source,notes
```

Generate templates:

```bash
Rscript scripts/write_cluster_annotation_templates.R
```

Validate edited annotations:

```bash
Rscript scripts/validate_cluster_annotations.R
```

Refresh app metadata after annotation changes:

```bash
Rscript scripts/refresh_app_metadata.R
```

## Legacy Cell-Type Overrides

`ATLAS_CELLTYPE_OVERRIDES_DIR` remains supported for simple two-column overrides:

```text
cluster_id,label
0,Cortex
1,Infected cells
```

For release work, prefer `ATLAS_CLUSTER_ANNOTATIONS_DIR` because it preserves evidence, review status, and cluster-source context.
