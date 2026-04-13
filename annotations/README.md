Place optional per-species gene annotation tables in this folder so the app can
search genes by common name as well as by locus ID.

Supported filenames:

- `medicago_gene_annotations.tsv`
- `glycine_gene_annotations.tsv`
- `lotus_gene_annotations.tsv`

Accepted columns:

- `gene_id`, `id`, or `locus_tag`
- `common_name`, `gene_name`, `symbol`, `acronym`, or `name`
- `synonyms`, `aliases`, or `alias` (optional)
- `description`, `geneProduct`, `product`, or `genePublication` (optional)

Recommended simple format:

```tsv
gene_id	common_name	synonyms	description
MtrunA17Chr5g0437741	NIN	MtNIN; NODULE INCEPTION	Transcription factor
Glyma.14G001500.Wm82.a6.v1	NIN	GmNIN	Soybean ortholog of NIN
LotjaGi5g1v0348600	NIN	LjNIN	Lotus ortholog of NIN
```

The app matches IDs after species-specific normalization, so Lotus IDs can be
provided either with or without the `_LC` suffix.
