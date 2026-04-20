`%||%` <- function(x, y) {
    if (is.null(x) || !length(x) || all(is.na(x))) y else x
}

pick_first_existing_path <- function(paths) {
    existing_path <- paths[file.exists(paths)][1]

    if (is.na(existing_path) || !length(existing_path)) {
        return(paths[[1]])
    }

    existing_path
}

app_slim_path <- function(path) {
    sub("\\.rds$", "_app_slim.rds", path, ignore.case = TRUE)
}

default_assay_layers <- function(assay_obj) {
    if (inherits(assay_obj, "Assay5")) {
        return(Layers(assay_obj))
    }

    slots <- c("counts", "data", "scale.data")
    slots[vapply(slots, function(slot_name) length(methods::slot(assay_obj, slot_name)) > 0, logical(1))]
}

ensure_assay_data_dgc <- function(assay_obj) {
    layers <- default_assay_layers(assay_obj)

    if (!("data" %in% layers)) {
        stop("The default assay must contain a 'data' layer.", call. = FALSE)
    }

    if (inherits(assay_obj, "Assay5")) {
        data_mat <- LayerData(assay_obj, layer = "data")
        if (!inherits(data_mat, "dgCMatrix")) {
            LayerData(assay_obj, layer = "data") <- methods::as(data_mat, "dgCMatrix")
        }
        return(assay_obj)
    }

    if (!inherits(assay_obj@data, "dgCMatrix")) {
        assay_obj@data <- methods::as(assay_obj@data, "dgCMatrix")
    }

    assay_obj
}

clear_assay_feature_metadata <- function(assay_obj) {
    if ("meta.data" %in% slotNames(assay_obj) && ncol(assay_obj@meta.data)) {
        assay_obj@meta.data <- assay_obj@meta.data[, 0, drop = FALSE]
    }

    if ("meta.features" %in% slotNames(assay_obj) && ncol(assay_obj@meta.features)) {
        assay_obj@meta.features <- assay_obj@meta.features[, 0, drop = FALSE]
    }

    if ("var.features" %in% slotNames(assay_obj) && length(assay_obj@var.features)) {
        assay_obj@var.features <- character(0)
    }

    assay_obj
}

clear_seurat_runtime_baggage <- function(obj, keep_misc = FALSE) {
    if ("commands" %in% slotNames(obj)) {
        obj@commands <- list()
    }
    if ("graphs" %in% slotNames(obj)) {
        obj@graphs <- list()
    }
    if ("neighbors" %in% slotNames(obj)) {
        obj@neighbors <- list()
    }
    if ("images" %in% slotNames(obj)) {
        obj@images <- list()
    }
    if ("tools" %in% slotNames(obj)) {
        obj@tools <- list()
    }
    if (!isTRUE(keep_misc) && "misc" %in% slotNames(obj)) {
        obj@misc <- list()
    }

    obj
}

slim_reduction_names <- function(obj, prefer_pca_fallback = FALSE) {
    available <- Reductions(obj)
    keep <- intersect(c("umap", "umap3d", "umap_3d"), available)

    if (!length(intersect(c("umap3d", "umap_3d"), keep)) &&
        isTRUE(prefer_pca_fallback) &&
        "pca" %in% available) {
        keep <- unique(c(keep, "pca"))
    }

    keep
}

slim_seurat_object <- function(obj, keep_meta_cols = NULL, keep_reductions = NULL, keep_misc = FALSE) {
    assay_name <- DefaultAssay(obj)
    assay_obj <- obj[[assay_name]]
    layers <- default_assay_layers(assay_obj)

    if (!("data" %in% layers)) {
        stop(sprintf("Default assay '%s' does not contain a data layer.", assay_name), call. = FALSE)
    }

    keep_reductions <- keep_reductions %||% slim_reduction_names(obj, prefer_pca_fallback = TRUE)

    if (inherits(assay_obj, "Assay5")) {
        slim_obj <- DietSeurat(
            object = obj,
            assays = assay_name,
            layers = "data",
            dimreducs = intersect(keep_reductions, Reductions(obj)),
            graphs = NULL,
            misc = keep_misc
        )
    } else {
        slim_obj <- suppressWarnings(DietSeurat(
            object = obj,
            assays = assay_name,
            dimreducs = intersect(keep_reductions, Reductions(obj)),
            graphs = NULL,
            misc = keep_misc,
            counts = FALSE,
            data = TRUE,
            scale.data = FALSE
        ))
    }

    slim_obj <- clear_seurat_runtime_baggage(slim_obj, keep_misc = keep_misc)
    slim_obj[[assay_name]] <- ensure_assay_data_dgc(slim_obj[[assay_name]])
    slim_obj[[assay_name]] <- clear_assay_feature_metadata(slim_obj[[assay_name]])

    if (length(VariableFeatures(slim_obj))) {
        VariableFeatures(slim_obj) <- character(0)
    }

    if (!is.null(keep_meta_cols)) {
        keep_meta_cols <- intersect(unique(keep_meta_cols), colnames(slim_obj@meta.data))
        slim_obj@meta.data <- slim_obj@meta.data[, keep_meta_cols, drop = FALSE]
    }

    slim_obj
}

object_audit_row <- function(label, path, obj) {
    assay_name <- DefaultAssay(obj)
    assay_obj <- obj[[assay_name]]
    layers <- default_assay_layers(assay_obj)

    data.frame(
        label = label,
        path = path,
        size_mb = round(file.info(path)$size / 1024^2, 1),
        cells = ncol(obj),
        genes = nrow(obj),
        meta_cols = ncol(obj@meta.data),
        reductions = paste(Reductions(obj), collapse = ","),
        layers = paste(layers, collapse = ","),
        stringsAsFactors = FALSE
    )
}

print_object_audit <- function(label, path, obj) {
    audit <- object_audit_row(label, path, obj)
    cat(
        sprintf(
            "[%s] %s | %.1f MB | %s cells | %s genes | %s meta cols | reductions: %s | layers: %s\n",
            audit$label,
            audit$path,
            audit$size_mb,
            format(audit$cells, big.mark = ","),
            format(audit$genes, big.mark = ","),
            audit$meta_cols,
            audit$reductions %||% "<none>",
            audit$layers %||% "<none>"
        )
    )
}
