#!/usr/bin/env Rscript

script_flag <- "--file="
script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub(script_flag, "", script_arg[grep(script_flag, script_arg)][1])
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(project_root)
source("scripts/atlas_dataset_utils.R")
source("R/atlas_annotations.R")

args <- commandArgs(trailingOnly = TRUE)
force_refresh <- "--force" %in% args
check_only <- "--check-only" %in% args
fail_on_stale <- "--fail-on-stale" %in% args

within_species_keys <- c("medicago", "glycine", "lotus")
integration_methods <- c("ComBat_BBKNN", "Seurat")

species_registry <- list(
    medicago = list(
        within_paths = c(
            ComBat_BBKNN = "within_species_integrated_datasets/ComBat_BBKNN/M_truncatula_clustered_dataset.rds",
            Seurat = "within_species_integrated_datasets/Seurat/M_truncatula_clustered_dataset.rds"
        )
    ),
    glycine = list(
        within_paths = c(
            ComBat_BBKNN = "within_species_integrated_datasets/ComBat_BBKNN/G_max_clustered_dataset.rds",
            Seurat = "within_species_integrated_datasets/Seurat/G_max_clustered_dataset.rds"
        )
    ),
    lotus = list(
        within_paths = c(
            ComBat_BBKNN = "within_species_integrated_datasets/ComBat_BBKNN/L_japonicus_clustered_dataset.rds",
            Seurat = "within_species_integrated_datasets/Seurat/L_japonicus_clustered_dataset.rds"
        )
    )
)

cross_integration_registry <- list(
    camex = list(
        path = "app_ready_integration/camex/clustered_dataset.rds",
        slim_path = "app_ready_integration/camex/clustered_dataset_app_slim.rds"
    ),
    saturn = list(
        path = "app_ready_integration/saturn/clustered_dataset.rds",
        slim_path = "app_ready_integration/saturn/clustered_dataset_app_slim.rds"
    )
)

annotation_paths <- c(
    medicago = "annotations/medicago_gene_annotations.tsv",
    glycine = "annotations/glycine_gene_annotations.tsv",
    lotus = "annotations/lotus_gene_annotations.tsv"
)

cluster_annotation_paths <- if (dir.exists(atlas_cluster_annotations_dir())) {
    list.files(
        atlas_cluster_annotations_dir(),
        pattern = "\\.(csv|tsv)$",
        full.names = TRUE
    )
} else {
    character(0)
}

resolved_within_paths <- unlist(lapply(within_species_keys, function(species_key) {
    vapply(integration_methods, function(integration_method) {
        path <- species_registry[[species_key]]$within_paths[[integration_method]]
        pick_first_existing_path(c(app_slim_path(path), path))
    }, character(1))
}), use.names = FALSE)

resolved_cross_paths <- vapply(names(cross_integration_registry), function(cross_key) {
    cfg <- cross_integration_registry[[cross_key]]
    pick_first_existing_path(c(cfg$slim_path, cfg$path))
}, character(1))

normalized_dependencies <- function(paths) {
    paths <- unique(normalizePath(paths[file.exists(paths)], winslash = "/", mustWork = FALSE))
    paths[nzchar(paths)]
}

dependency_mtime <- function(paths) {
    info <- file.info(paths)
    if (!nrow(info)) {
        return(as.POSIXct(NA))
    }

    max(info$mtime, na.rm = TRUE)
}

target_mtime <- function(paths) {
    info <- file.info(paths)
    if (!nrow(info)) {
        return(as.POSIXct(NA))
    }

    min(info$mtime, na.rm = TRUE)
}

builder_specs <- list(
    list(
        id = "gene_catalogs",
        label = "Gene catalogs and Camex feature lookup",
        script = "scripts/build_gene_catalog_cache.R",
        targets = c(
            "metadata/gene_catalogs/medicago.tsv",
            "metadata/gene_catalogs/glycine.tsv",
            "metadata/gene_catalogs/lotus.tsv",
            "metadata/cross_feature_lookup.tsv"
        ),
        dependencies = c(
            resolved_within_paths,
            resolved_cross_paths[["camex"]],
            annotation_paths,
            "scripts/build_gene_catalog_cache.R",
            "scripts/atlas_dataset_utils.R"
        )
    ),
    list(
        id = "atlas_summary",
        label = "Atlas summary strip cache",
        script = "scripts/build_atlas_summary_cache.R",
        targets = c("metadata/atlas_summary.tsv"),
        dependencies = c(
            resolved_within_paths,
            resolved_cross_paths,
            "scripts/build_atlas_summary_cache.R",
            "scripts/atlas_dataset_utils.R"
        )
    ),
    list(
        id = "startup_ui",
        label = "Startup UI choice and cluster caches",
        script = "scripts/build_startup_ui_cache.R",
        targets = c(
            "metadata/ui_choice_cache.tsv",
            "metadata/ui_cluster_lookup_cache.tsv"
        ),
        dependencies = c(
            resolved_within_paths,
            resolved_cross_paths,
            cluster_annotation_paths,
            "scripts/build_startup_ui_cache.R",
            "scripts/atlas_dataset_utils.R",
            "R/atlas_annotations.R"
        )
    )
)

builder_needs_refresh <- function(spec) {
    dependencies <- normalized_dependencies(spec$dependencies)
    targets_exist <- all(file.exists(spec$targets))

    if (!length(dependencies)) {
        return(list(
            refresh = FALSE,
            reason = "No readable dependencies were found."
        ))
    }

    if (!targets_exist) {
        return(list(
            refresh = TRUE,
            reason = "One or more target files are missing."
        ))
    }

    newest_dependency <- dependency_mtime(dependencies)
    oldest_target <- target_mtime(spec$targets)

    if (is.na(newest_dependency) || is.na(oldest_target)) {
        return(list(
            refresh = TRUE,
            reason = "Could not compare modification times."
        ))
    }

    if (newest_dependency > oldest_target) {
        return(list(
            refresh = TRUE,
            reason = sprintf(
                "A dependency changed after the cache was written (%s > %s).",
                format(newest_dependency, tz = ""),
                format(oldest_target, tz = "")
            )
        ))
    }

    list(
        refresh = FALSE,
        reason = "Targets are newer than every dependency."
    )
}

run_builder <- function(spec) {
    status <- system2("Rscript", spec$script)
    identical(status, 0L)
}

cat("Checking app metadata caches...\n")
if (force_refresh) {
    cat("Force mode enabled: every managed cache will be rebuilt.\n")
}
if (check_only) {
    cat("Check-only mode enabled: no files will be regenerated.\n")
}
if (check_only && fail_on_stale) {
    cat("Fail-on-stale mode enabled: stale caches will return a non-zero exit code.\n")
}

results <- lapply(builder_specs, function(spec) {
    refresh_check <- builder_needs_refresh(spec)
    should_refresh <- isTRUE(force_refresh) || isTRUE(refresh_check$refresh)
    status <- "up_to_date"
    ok <- TRUE

    if (should_refresh && !check_only) {
        cat(sprintf("\n[%s] Refreshing %s...\n", spec$id, spec$label))
        ok <- run_builder(spec)
        status <- if (ok) "rebuilt" else "failed"
    } else if (should_refresh) {
        status <- "stale"
        cat(sprintf("\n[%s] %s\n", spec$id, refresh_check$reason))
    } else {
        cat(sprintf("\n[%s] %s\n", spec$id, refresh_check$reason))
    }

    data.frame(
        builder_id = spec$id,
        label = spec$label,
        status = status,
        reason = refresh_check$reason,
        stringsAsFactors = FALSE
    )
})

result_tbl <- do.call(rbind, results)

cat("\nSummary:\n")
print(result_tbl, row.names = FALSE)

failed <- result_tbl$status == "failed"

if (any(failed)) {
    stop("One or more cache rebuild steps failed.", call. = FALSE)
}

if (check_only && any(result_tbl$status == "stale")) {
    cat("\nSome caches are stale. Re-run without --check-only to rebuild them.\n")
    if (fail_on_stale) {
        quit(status = 1L)
    }
} else if (any(result_tbl$status == "rebuilt")) {
    cat("\nAll stale caches were refreshed.\n")
} else {
    cat("\nAll managed caches are already up to date.\n")
}
