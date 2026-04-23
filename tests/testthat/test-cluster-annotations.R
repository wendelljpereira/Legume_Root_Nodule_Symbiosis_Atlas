make_annotation_test_object <- function() {
    counts <- Matrix::Matrix(
        c(
            1, 0, 3, 0,
            0, 2, 0, 4
        ),
        nrow = 2,
        sparse = TRUE,
        dimnames = list(c("geneA", "geneB"), paste0("cell", 1:4))
    )
    obj <- Seurat::CreateSeuratObject(counts = counts)
    SeuratObject::Idents(obj) <- factor(c("0", "0", "1", "1"), levels = c("0", "1"))
    obj$Rank_1st <- c("0", "0", "1", "1")
    obj
}

annotation_fixture <- function(...) {
    defaults <- data.frame(
        dataset_key = "toy",
        cluster_source = c("__idents__", "__idents__", "Rank_1st", "Rank_1st"),
        cluster_id = c("0", "1", "0", "1"),
        cell_type_label = c("Infected", "Cortex", "Infected", "Cortex"),
        annotation_status = "curated",
        confidence_score = c(0.9, 0.8, 0.9, 0.8),
        marker_evidence = "fixture markers",
        reference_source = "fixture reference",
        notes = "",
        stringsAsFactors = FALSE
    )
    values <- list(...)
    for (nm in names(values)) defaults[[nm]] <- values[[nm]]
    defaults
}

test_that("cluster annotations overlay labels without changing raw cluster IDs", {
    obj <- make_annotation_test_object()
    annotated <- atlas_env$apply_cluster_annotation_table(
        obj,
        annotation_fixture(),
        dataset_key = "toy",
        path = "fixture.csv"
    )

    expect_equal(as.character(SeuratObject::Idents(annotated)), c("0", "0", "1", "1"))
    expect_equal(as.character(annotated$Rank_1st), c("0", "0", "1", "1"))
    expect_equal(as.character(annotated$cluster_label), c("Infected", "Infected", "Cortex", "Cortex"))
    expect_equal(as.character(annotated$Rank_1st_label), c("Infected", "Infected", "Cortex", "Cortex"))
    expect_equal(as.character(annotated$celltype), as.character(annotated$cluster_label))
})

test_that("cluster annotation validation rejects malformed tables", {
    obj <- make_annotation_test_object()

    expect_error(
        atlas_env$normalize_cluster_annotation_table(annotation_fixture()[, -1], dataset_key = "toy", path = "missing.csv"),
        "missing required column"
    )

    duplicate_tbl <- rbind(annotation_fixture(), annotation_fixture()[1, ])
    expect_error(
        atlas_env$validate_cluster_annotation_table(obj, duplicate_tbl, dataset_key = "toy", path = "duplicate.csv"),
        "duplicate"
    )

    unknown_tbl <- annotation_fixture(cluster_id = c("0", "1", "9", "1"))
    expect_error(
        atlas_env$validate_cluster_annotation_table(obj, unknown_tbl, dataset_key = "toy", path = "unknown.csv"),
        "not present"
    )

    empty_label_tbl <- annotation_fixture(cell_type_label = c("Infected", "", "Infected", "Cortex"))
    expect_error(
        atlas_env$validate_cluster_annotation_table(obj, empty_label_tbl, dataset_key = "toy", path = "empty.csv"),
        "empty cell_type_label"
    )
})

test_that("legacy celltype overrides remain supported", {
    obj <- make_annotation_test_object()
    override_dir <- tempfile("legacy-overrides")
    dir.create(override_dir)
    utils::write.csv(
        data.frame(cluster_id = c("0", "1"), label = c("Legacy A", "Legacy B")),
        file.path(override_dir, "toy.csv"),
        row.names = FALSE
    )

    annotated <- atlas_env$apply_cluster_annotation_overlay(
        obj,
        dataset_key = "toy",
        annotations_dir = tempfile("missing-annotations"),
        legacy_overrides_dir = override_dir
    )

    expect_equal(as.character(annotated$cluster_label), c("Legacy A", "Legacy A", "Legacy B", "Legacy B"))
})
