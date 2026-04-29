test_that("species-specific gene canonicalization is stable", {
    expect_equal(atlas_env$canonicalize_gene_ids("lotus", "Lj0g00001_LC"), "Lj0g00001")
    expect_equal(atlas_env$canonicalize_gene_ids("medicago", " MtrunA17Chr1g000001 "), "MtrunA17Chr1g000001")
    expect_true(is.na(atlas_env$canonicalize_gene_ids("glycine", "")))
})

test_that("orthogroup gene lists parse compact OrthoFinder cells", {
    expect_equal(
        atlas_env$split_orthogroup_genes("geneA, geneB,, geneC"),
        c("geneA", "geneB", "geneC")
    )
    expect_equal(atlas_env$split_orthogroup_genes(NA_character_), character(0))
})

test_that("plot feature guardrail trims large mapped panels", {
    resolution <- list(
        plot_features = paste0("gene", 1:5),
        plot_table = tibble::tibble(feature_id = paste0("gene", 1:5), label = letters[1:5]),
        label_map = stats::setNames(LETTERS[1:5], paste0("gene", 1:5))
    )

    trimmed <- atlas_env$enforce_plot_feature_limit(resolution, limit = 3)

    expect_true(trimmed$feature_limit_exceeded)
    expect_equal(trimmed$feature_count_before_limit, 5)
    expect_equal(trimmed$plot_features, paste0("gene", 1:3))
    expect_equal(trimmed$features_omitted, paste0("gene", 4:5))
    expect_equal(nrow(trimmed$plot_table), 3)
    expect_equal(names(trimmed$label_map), paste0("gene", 1:3))
})

test_that("average expression summaries are computed by metadata group", {
    mat <- Matrix::Matrix(
        c(
            1, 3, 5, 7,
            2, 4, 6, 8
        ),
        nrow = 2,
        byrow = TRUE,
        sparse = TRUE,
        dimnames = list(c("geneA", "geneB"), paste0("cell", 1:4))
    )
    metadata <- tibble::tibble(cluster = c("0", "0", "1", "1"))

    avg <- atlas_env$compute_average_expression_matrix(mat, metadata, c("geneA", "geneB"), "cluster")

    expect_equal(rownames(avg), c("geneA", "geneB"))
    expect_equal(colnames(avg), c("0", "1"))
    expect_equal(unname(avg[, "0"]), c(2, 3))
    expect_equal(unname(avg[, "1"]), c(6, 7))
})

test_that("sample name overrides apply to display metadata", {
    old_override_path <- atlas_env$sample_name_overrides_path
    atlas_env$sample_name_overrides_path <- file.path(project_root, "metadata", "sample_name_overrides.tsv")
    on.exit({
        atlas_env$sample_name_overrides_path <- old_override_path
    }, add = TRUE)

    lookup <- atlas_env$sample_name_override_lookup("medicago")

    expect_equal(unname(lookup["SRR20997723_extended"]), "Cervantes_48hpi_rep1")
    expect_equal(unname(lookup["SRR20997723"]), "Cervantes_48hpi_rep1")
    expect_equal(
        unname(atlas_env$sample_name_override_lookup("lotus")["Lotus_WT_Mloti_10dpi_rep1"]),
        "Frank_nodules_10d_rep1"
    )
    expect_equal(
        unname(atlas_env$sample_name_override_lookup("glycine")["GSE226149_GSM7065808"]),
        "Cervantes_roots_rep1"
    )

    counts <- Matrix::sparseMatrix(
        i = c(1, 2, 2),
        j = c(1, 1, 2),
        x = c(1, 2, 1),
        dims = c(2, 2),
        dimnames = list(c("gene1", "gene2"), c("cell1", "cell2"))
    )
    obj <- Seurat::CreateSeuratObject(
        counts = counts,
        meta.data = data.frame(
            sample_name = c("SRR20997723_extended", "CRA007122_CRR513370"),
            sample = c("SRR20997723_extended", "CRA007122_CRR513370"),
            species = c("Medicago truncatula", "Glycine max"),
            row.names = c("cell1", "cell2"),
            check.names = FALSE
        )
    )

    renamed <- atlas_env$apply_sample_name_overrides(obj)

    expect_equal(renamed$sample_name[[1]], "Cervantes_48hpi_rep1")
    expect_equal(renamed$sample[[1]], "Cervantes_48hpi_rep1")
    expect_equal(renamed$sample_name[[2]], "Liu_12dpi")
})

test_that("sample display names are synchronized and ordered for plotting", {
    old_override_path <- atlas_env$sample_name_overrides_path
    atlas_env$sample_name_overrides_path <- file.path(project_root, "metadata", "sample_name_overrides.tsv")
    on.exit({
        atlas_env$sample_name_overrides_path <- old_override_path
    }, add = TRUE)

    counts <- Matrix::sparseMatrix(
        i = c(1, 2, 2),
        j = c(1, 1, 2),
        x = c(1, 2, 1),
        dims = c(2, 2),
        dimnames = list(c("gene1", "gene2"), c("cell1", "cell2"))
    )
    obj <- Seurat::CreateSeuratObject(
        counts = counts,
        meta.data = data.frame(
            sample = c("Pereira_A17_roots", "SRR20997723_extended"),
            Sample = c("Time_0h_A17", "SRR20997723_extended"),
            species = rep("Medicago truncatula", 2),
            row.names = c("cell1", "cell2"),
            check.names = FALSE
        )
    )

    renamed <- atlas_env$apply_sample_name_overrides(obj)

    expect_equal(as.character(renamed$Sample), c("Pereira_A17_roots", "Cervantes_48hpi_rep1"))
    expect_equal(as.character(renamed$sample), c("Pereira_A17_roots", "Cervantes_48hpi_rep1"))
    expect_equal(
        atlas_env$sample_display_timepoint(c(
            "Frank_roots_5d_rep1",
            "Frank_roots_10d_rep1",
            "Pereira_A17_roots",
            "Ye_14dpi_rep1",
            "Liu_0.5hpi_rep1",
            "Pereira_A17_96hpi"
        )),
        c("5dpi", "10dpi", "Roots", "14dpi", "0.5h", "96h")
    )
    expect_equal(
        atlas_env$metadata_level_order(
            c("Cervantes_48hpi_rep1", "Pereira_A17_roots", "Liu_roots_rep1"),
            "Sample"
        ),
        c("Pereira_A17_roots", "Liu_roots_rep1", "Cervantes_48hpi_rep1")
    )
    expect_equal(
        atlas_env$metadata_level_order(
            c("Frank_nodules_10d_rep2", "Frank_roots_10d_rep1", "Frank_roots_5d_rep2"),
            "Sample"
        ),
        c("Frank_roots_5d_rep2", "Frank_roots_10d_rep1", "Frank_nodules_10d_rep2")
    )
    expect_equal(
        atlas_env$metadata_level_order(
            c(
                "Zhang_15dpi_rep2",
                "Liu_roots",
                "Cervantes_roots_rep2",
                "Sun_roots",
                "Liu_12dpi",
                "Cervantes_28dpi_rep1"
            ),
            "Sample"
        ),
        c(
            "Liu_roots",
            "Sun_roots",
            "Cervantes_roots_rep2",
            "Liu_12dpi",
            "Zhang_15dpi_rep2",
            "Cervantes_28dpi_rep1"
        )
    )
})
