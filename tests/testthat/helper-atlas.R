candidate_roots <- unique(normalizePath(
    c(getwd(), file.path(getwd(), ".."), file.path(getwd(), "..", "..")),
    mustWork = FALSE
))
project_root <- candidate_roots[file.exists(file.path(candidate_roots, "R", "atlas_core.R"))][1]
if (is.na(project_root) || !nzchar(project_root)) {
    stop("Could not locate project root containing R/atlas_core.R.", call. = FALSE)
}

atlas_env <- new.env(parent = globalenv())
old_wd <- setwd(project_root)
on.exit(setwd(old_wd), add = TRUE)
sys.source(file.path(project_root, "R", "atlas_core.R"), envir = atlas_env)
