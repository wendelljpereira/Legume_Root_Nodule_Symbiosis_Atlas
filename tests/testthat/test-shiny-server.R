test_that("Shiny server initializes without eager dataset loading", {
    app_env <- new.env(parent = globalenv())
    old_wd <- setwd(project_root)
    on.exit(setwd(old_wd), add = TRUE)
    sys.source(file.path(project_root, "app.R"), envir = app_env)

    shiny::testServer(app_env$server, {
        expect_equal(current_comparison_source_species(), "medicago")
        expect_equal(staged_source_genes(), character(0))
        expect_equal(selected_source_genes(), character(0))
        expect_true(isTRUE(app_unlocked()))
    })
})
