#!/usr/bin/env Rscript

script_flag <- "--file="
script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub(script_flag, "", script_arg[grep(script_flag, script_arg)][1])
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(project_root)

args <- commandArgs(trailingOnly = TRUE)
skip_cache <- "--skip-cache" %in% args
skip_smoke <- "--skip-smoke" %in% args
skip_tests <- "--skip-tests" %in% args
skip_browser <- "--skip-browser" %in% args
skip_annotations <- "--skip-annotations" %in% args
skip_bundle_size <- "--skip-bundle-size" %in% args

run_step <- function(label, command, args = character(), required = TRUE) {
    cat(sprintf("\n==> %s\n", label))
    status <- system2(command, args = shQuote(args))
    ok <- identical(status, 0L)

    if (!ok && required) {
        stop(sprintf("Release check failed: %s", label), call. = FALSE)
    }

    invisible(ok)
}

cat("Checking R source parseability...\n")
invisible(parse("app.R"))
invisible(parse("R/atlas_core.R"))
invisible(parse("R/atlas_ui_tabs.R"))
invisible(parse("R/atlas_annotations.R"))
cat("Parse checks passed.\n")

if (!skip_cache) {
    run_step(
        "metadata cache freshness",
        "Rscript",
        c("scripts/refresh_app_metadata.R", "--check-only", "--fail-on-stale")
    )
}

if (!skip_smoke) {
    run_step("dataset smoke test", "Rscript", "scripts/smoke_test.R")
}

if (!skip_annotations) {
    run_step("cluster annotation validation", "Rscript", "scripts/validate_cluster_annotations.R")
}

if (!skip_bundle_size) {
    run_step("ShinyApps.io bundle size estimate", "Rscript", "scripts/check_bundle_size.R")
}

if (!skip_tests) {
    if (!requireNamespace("testthat", quietly = TRUE)) {
        stop("Package 'testthat' is required for release unit tests. Install it with install.packages('testthat').", call. = FALSE)
    }

    run_step(
        "unit tests",
        "Rscript",
        c("-e", "testthat::test_dir('tests/testthat', reporter = 'summary')")
    )
}

if (!skip_browser) {
    if (!dir.exists("node_modules/playwright")) {
        stop("Playwright is not installed. Run `npm install` before browser release checks.", call. = FALSE)
    }

    run_step("Firefox browser smoke test", "npm", c("run", "test:e2e:firefox"))
}

cat("\nAll requested release checks passed.\n")
