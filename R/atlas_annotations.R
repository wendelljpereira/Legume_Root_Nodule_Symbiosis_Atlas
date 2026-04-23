# Runtime cluster annotation overlays.
# These helpers keep curated labels in editable CSV files instead of writing them
# permanently into the Seurat RDS objects.

cluster_annotation_required_columns <- c(
    "dataset_key",
    "cluster_source",
    "cluster_id",
    "cell_type_label",
    "annotation_status",
    "confidence_score",
    "marker_evidence",
    "reference_source",
    "notes"
)

atlas_cluster_annotations_dir <- function(default = "annotations/cluster_annotations") {
    value <- Sys.getenv("ATLAS_CLUSTER_ANNOTATIONS_DIR", unset = "")
    if (!nzchar(value)) default else value
}

atlas_legacy_celltype_overrides_dir <- function(default = "celltype_overrides") {
    value <- Sys.getenv("ATLAS_CELLTYPE_OVERRIDES_DIR", unset = "")
    if (!nzchar(value)) default else value
}

normalize_cluster_annotation_source <- function(values) {
    values_chr <- trimws(as.character(values))
    values_lower <- tolower(values_chr)

    values_chr[is.na(values_chr) | !nzchar(values_chr)] <- "__idents__"
    values_chr[values_lower %in% c("idents", "identity", "ident", "saved identities", "__idents__")] <- "__idents__"
    values_chr
}

cluster_annotation_label_column <- function(cluster_source) {
    cluster_source <- normalize_cluster_annotation_source(cluster_source)[1]

    if (is.na(cluster_source) || !nzchar(cluster_source) || identical(cluster_source, "__idents__")) {
        return("cluster_label")
    }

    if (grepl("_label$", cluster_source)) {
        return(cluster_source)
    }

    paste0(cluster_source, "_label")
}

cluster_annotation_file_candidates <- function(dataset_key, annotations_dir = atlas_cluster_annotations_dir()) {
    file.path(
        annotations_dir,
        c(paste0(dataset_key, ".csv"), paste0(dataset_key, ".tsv"))
    )
}

cluster_annotation_file_path <- function(dataset_key, annotations_dir = atlas_cluster_annotations_dir()) {
    candidates <- cluster_annotation_file_candidates(dataset_key, annotations_dir)
    existing <- candidates[file.exists(candidates)][1]

    if (is.na(existing) || !length(existing)) {
        return(candidates[[1]])
    }

    existing
}

read_cluster_annotation_file <- function(dataset_key, annotations_dir = atlas_cluster_annotations_dir(), required = FALSE) {
    if (!nzchar(annotations_dir) || !dir.exists(annotations_dir)) {
        if (isTRUE(required)) {
            stop(sprintf("Cluster annotation directory not found: %s", annotations_dir), call. = FALSE)
        }
        return(NULL)
    }

    path <- cluster_annotation_file_path(dataset_key, annotations_dir)
    if (!file.exists(path)) {
        if (isTRUE(required)) {
            stop(sprintf("Cluster annotation file is missing for dataset '%s': %s", dataset_key, path), call. = FALSE)
        }
        return(NULL)
    }

    sep <- if (grepl("\\.tsv$", path, ignore.case = TRUE)) "\t" else ","
    annotation_tbl <- tryCatch(
        utils::read.delim(path, sep = sep, stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e) {
            stop(sprintf("Failed to read cluster annotation file '%s': %s", path, conditionMessage(e)), call. = FALSE)
        }
    )

    normalize_cluster_annotation_table(annotation_tbl, dataset_key = dataset_key, path = path)
}

normalize_cluster_annotation_table <- function(annotation_tbl, dataset_key, path = "<annotation table>") {
    missing_cols <- setdiff(cluster_annotation_required_columns, colnames(annotation_tbl))
    if (length(missing_cols)) {
        stop(sprintf(
            "Cluster annotation file '%s' is missing required column(s): %s",
            path,
            paste(missing_cols, collapse = ", ")
        ), call. = FALSE)
    }

    annotation_tbl <- annotation_tbl[, cluster_annotation_required_columns, drop = FALSE]
    annotation_tbl$dataset_key <- trimws(as.character(annotation_tbl$dataset_key))
    annotation_tbl$cluster_source <- normalize_cluster_annotation_source(annotation_tbl$cluster_source)
    annotation_tbl$cluster_id <- trimws(as.character(annotation_tbl$cluster_id))
    annotation_tbl$cell_type_label <- trimws(as.character(annotation_tbl$cell_type_label))
    annotation_tbl$annotation_status <- trimws(as.character(annotation_tbl$annotation_status))
    annotation_tbl$confidence_score <- suppressWarnings(as.numeric(annotation_tbl$confidence_score))
    annotation_tbl$marker_evidence <- trimws(as.character(annotation_tbl$marker_evidence))
    annotation_tbl$reference_source <- trimws(as.character(annotation_tbl$reference_source))
    annotation_tbl$notes <- trimws(as.character(annotation_tbl$notes))

    if (any(is.na(annotation_tbl$dataset_key) | !nzchar(annotation_tbl$dataset_key))) {
        stop(sprintf("Cluster annotation file '%s' contains empty dataset_key values.", path), call. = FALSE)
    }

    unexpected_dataset_keys <- setdiff(unique(annotation_tbl$dataset_key), dataset_key)
    if (length(unexpected_dataset_keys)) {
        stop(sprintf(
            "Cluster annotation file '%s' is for dataset_key '%s' but contains: %s",
            path,
            dataset_key,
            paste(unexpected_dataset_keys, collapse = ", ")
        ), call. = FALSE)
    }

    if (any(is.na(annotation_tbl$cluster_source) | !nzchar(annotation_tbl$cluster_source))) {
        stop(sprintf("Cluster annotation file '%s' contains empty cluster_source values.", path), call. = FALSE)
    }

    if (any(is.na(annotation_tbl$cluster_id) | !nzchar(annotation_tbl$cluster_id))) {
        stop(sprintf("Cluster annotation file '%s' contains empty cluster_id values.", path), call. = FALSE)
    }

    if (any(is.na(annotation_tbl$cell_type_label) | !nzchar(annotation_tbl$cell_type_label))) {
        stop(sprintf("Cluster annotation file '%s' contains empty cell_type_label values.", path), call. = FALSE)
    }

    duplicated_rows <- duplicated(annotation_tbl[, c("dataset_key", "cluster_source", "cluster_id")])
    if (any(duplicated_rows)) {
        duplicates <- unique(paste(
            annotation_tbl$cluster_source[duplicated_rows],
            annotation_tbl$cluster_id[duplicated_rows],
            sep = ":"
        ))
        stop(sprintf(
            "Cluster annotation file '%s' contains duplicate dataset/source/cluster rows: %s",
            path,
            paste(utils::head(duplicates, 10), collapse = ", ")
        ), call. = FALSE)
    }

    annotation_tbl
}

cluster_source_values <- function(obj, cluster_source) {
    cluster_source <- normalize_cluster_annotation_source(cluster_source)[1]

    if (identical(cluster_source, "__idents__")) {
        return(as.character(SeuratObject::Idents(obj)))
    }

    if (!(cluster_source %in% colnames(obj@meta.data))) {
        stop(sprintf(
            "Cluster annotation source '%s' is not present in this dataset metadata.",
            cluster_source
        ), call. = FALSE)
    }

    as.character(obj@meta.data[[cluster_source]])
}

validate_cluster_annotation_table <- function(obj, annotation_tbl, dataset_key, path = "<annotation table>") {
    annotation_tbl <- normalize_cluster_annotation_table(annotation_tbl, dataset_key = dataset_key, path = path)

    for (cluster_source in unique(annotation_tbl$cluster_source)) {
        available_values <- unique(cluster_source_values(obj, cluster_source))
        available_values <- available_values[!is.na(available_values) & nzchar(available_values)]
        expected_values <- unique(annotation_tbl$cluster_id[annotation_tbl$cluster_source == cluster_source])
        missing_values <- setdiff(expected_values, available_values)

        if (length(missing_values)) {
            stop(sprintf(
                "Cluster annotation file '%s' contains cluster_id value(s) not present in dataset '%s' for source '%s': %s",
                path,
                dataset_key,
                cluster_source,
                paste(utils::head(missing_values, 12), collapse = ", ")
            ), call. = FALSE)
        }
    }

    annotation_tbl
}

initialize_cluster_annotation_labels <- function(obj) {
    obj$cluster_label <- as.character(SeuratObject::Idents(obj))
    obj$celltype <- obj$cluster_label
    obj
}

apply_cluster_annotation_table <- function(obj, annotation_tbl, dataset_key, path = "<annotation table>") {
    annotation_tbl <- validate_cluster_annotation_table(obj, annotation_tbl, dataset_key = dataset_key, path = path)

    for (cluster_source in unique(annotation_tbl$cluster_source)) {
        source_values <- cluster_source_values(obj, cluster_source)
        label_col <- cluster_annotation_label_column(cluster_source)
        rows <- annotation_tbl[annotation_tbl$cluster_source == cluster_source, , drop = FALSE]
        label_lookup <- stats::setNames(rows$cell_type_label, rows$cluster_id)
        labels <- unname(label_lookup[source_values])
        labels[is.na(labels) | !nzchar(labels)] <- source_values[is.na(labels) | !nzchar(labels)]
        obj@meta.data[[label_col]] <- labels
    }

    if ("cluster_label" %in% colnames(obj@meta.data)) {
        obj$celltype <- obj$cluster_label
    }

    obj
}

read_legacy_celltype_override <- function(override_id, overrides_dir = atlas_legacy_celltype_overrides_dir()) {
    if (!nzchar(overrides_dir) || !dir.exists(overrides_dir)) {
        return(NULL)
    }

    candidate_paths <- file.path(
        overrides_dir,
        c(paste0(override_id, ".csv"), paste0(override_id, ".tsv"))
    )
    path <- candidate_paths[file.exists(candidate_paths)][1]
    if (is.na(path) || !length(path)) {
        return(NULL)
    }

    sep <- if (grepl("\\.tsv$", path, ignore.case = TRUE)) "\t" else ","
    df <- tryCatch(
        utils::read.delim(path, sep = sep, stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e) {
            warning(sprintf("Failed to read celltype override '%s': %s", path, conditionMessage(e)), call. = FALSE)
            NULL
        }
    )
    if (is.null(df) || !nrow(df)) return(NULL)

    cols <- tolower(colnames(df))
    id_col <- which(cols %in% c("cluster_id", "cluster", "id"))[1]
    label_col <- which(cols %in% c("label", "celltype", "cell_type", "cluster_label"))[1]
    if (is.na(id_col) || is.na(label_col)) {
        warning(sprintf(
            "Celltype override '%s' must have columns 'cluster_id' and 'label' (received: %s). Skipping.",
            path, paste(colnames(df), collapse = ", ")
        ), call. = FALSE)
        return(NULL)
    }

    data.frame(
        cluster_id = as.character(df[[id_col]]),
        label = as.character(df[[label_col]]),
        stringsAsFactors = FALSE
    )
}

apply_legacy_celltype_override <- function(obj, override_id, overrides_dir = atlas_legacy_celltype_overrides_dir()) {
    mapping <- read_legacy_celltype_override(override_id, overrides_dir = overrides_dir)
    if (is.null(mapping) || !nrow(mapping)) {
        return(obj)
    }

    current_ids <- as.character(SeuratObject::Idents(obj))
    lookup <- stats::setNames(mapping$label, mapping$cluster_id)
    replaced <- unname(lookup[current_ids])
    obj$cluster_label <- ifelse(is.na(replaced) | !nzchar(replaced), current_ids, replaced)
    obj$celltype <- obj$cluster_label
    obj
}

apply_cluster_annotation_overlay <- function(
    obj,
    dataset_key,
    annotations_dir = atlas_cluster_annotations_dir(),
    legacy_overrides_dir = atlas_legacy_celltype_overrides_dir()
) {
    obj <- initialize_cluster_annotation_labels(obj)

    annotation_path <- cluster_annotation_file_path(dataset_key, annotations_dir = annotations_dir)
    annotation_tbl <- read_cluster_annotation_file(dataset_key, annotations_dir = annotations_dir, required = FALSE)
    if (!is.null(annotation_tbl)) {
        obj <- apply_cluster_annotation_table(obj, annotation_tbl, dataset_key = dataset_key, path = annotation_path)
    }

    apply_legacy_celltype_override(obj, dataset_key, overrides_dir = legacy_overrides_dir)
}

cluster_annotation_sources_for_object <- function(obj) {
    metadata_cols <- colnames(obj@meta.data)
    unique(c(
        "__idents__",
        c("Rank_1st", "Rank_2nd", "Rank_3rd", "Rank_4th", "Rank_5th")[c("Rank_1st", "Rank_2nd", "Rank_3rd", "Rank_4th", "Rank_5th") %in% metadata_cols]
    ))
}

cluster_annotation_template_rows <- function(obj, dataset_key, cluster_sources = cluster_annotation_sources_for_object(obj)) {
    rows <- lapply(cluster_sources, function(cluster_source) {
        values <- unique(cluster_source_values(obj, cluster_source))
        values <- values[!is.na(values) & nzchar(values)]
        values <- sort(values)

        data.frame(
            dataset_key = dataset_key,
            cluster_source = cluster_source,
            cluster_id = values,
            cell_type_label = values,
            annotation_status = "unreviewed",
            confidence_score = NA_real_,
            marker_evidence = "",
            reference_source = "",
            notes = "",
            stringsAsFactors = FALSE
        )
    })

    do.call(rbind, rows)
}
