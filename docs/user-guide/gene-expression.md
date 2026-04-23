# Gene Expression Workflow

## Search For Genes

Use the gene picker in a within-species tab to search by gene ID and available annotation text. Selected genes appear as removable tokens.

If you have a list of genes, use the import/paste workflow when available. The app drops genes that are not present in the active dataset and shows a warning.

## Generate Plots

After selecting at least one gene, click **Generate the expression plots**.

If no gene is selected, the app displays a warning and does not trigger plot rendering.

## Use Marker Genes As A Starting Point

Cluster marker sections can add the top marker genes from a selected cluster to the gene panel. These buttons only stage genes; they do not bypass the Generate button.

This is useful for asking questions such as:

- Which clusters express known nodulation markers?
- Are marker genes spatially or transcriptionally concentrated?
- Do candidate genes follow expected infected-cell, cortex, meristem, or primordium patterns?

## Practical Tips

- Start with 1-5 genes when exploring a new question.
- Use the averaged-expression heatmap and dot plot to compare many genes more compactly.
- Split large gene lists into batches, especially in cross-species tabs where one source gene may map to multiple orthologous features.
- Use downloads for records, but keep atlas version/date in your notes.
