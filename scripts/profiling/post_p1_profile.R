#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(htmlwidgets)
    library(profvis)
    library(scCustomize)
})

profile_dir <- file.path("scripts", "profiling")
dir.create(profile_dir, recursive = TRUE, showWarnings = FALSE)

env <- new.env(parent = globalenv())
sys.source("app.R", envir = env)

genes <- c(
    "05469365-9109-4067-9320-3e2a89459240",
    "06b706c5-0c87-4699-b519-2411dda9000d",
    "0b3aabc5-10e1-47ae-bf81-c1e8eaa77209",
    "1.8-mRNA-MtrunA17Chr2-13608494-13610171",
    "1.8-mRNA-MtrunA17Chr3-31913516-31914185"
)

profile_target <- function() {
    obj <- env$get_within_object("medicago", "ComBat_BBKNN")
    resolution <- env$resolve_target_mapping(
        source_species = "medicago",
        source_genes = genes,
        target_species = "medicago",
        integration_method = "ComBat_BBKNN",
        cross_space = FALSE
    )

    invisible(env$build_expression_heatmap_plot(
        obj = obj,
        feature_ids = resolution$plot_features,
        label_map = resolution$label_map,
        group_by = "Rank_1st",
        colorblind_safe = FALSE
    ))

    invisible(env$build_expression_ridge_plot(
        obj = obj,
        feature_ids = resolution$plot_features,
        label_map = resolution$label_map,
        group_by = "Rank_1st",
        colorblind_safe = FALSE
    ))

    invisible(suppressWarnings(
        scCustomize::DotPlot_scCustom(
            seurat_object = obj,
            features = resolution$plot_features,
            group.by = "Rank_1st"
        )
    ))
}

# Warm the target twice so the recorded timing reflects a steady-state
# warm-cache tab switch rather than first-use Seurat/scCustomize setup.
target_obj <- env$get_within_object("medicago", "ComBat_BBKNN")
invisible(env$get_cached_data_matrix(target_obj))
profile_target()
profile_target()

gc()
warm_elapsed <- unname(system.time(profile_target())[["elapsed"]])

prof <- profvis::profvis(
    expr = profile_target(),
    interval = 0.01
)

html_path <- file.path(profile_dir, "post_p1_profile.html")
htmlwidgets::saveWidget(
    widget = prof,
    file = html_path,
    selfcontained = FALSE,
    libdir = file.path(profile_dir, "post_p1_profile_files")
)

cat(sprintf("WarmElapsedSeconds: %.3f\n", warm_elapsed))
cat(sprintf("ProfileHTML: %s\n", html_path))
