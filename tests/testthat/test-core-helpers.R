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
