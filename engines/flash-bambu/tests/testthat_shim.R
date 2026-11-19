# Minimal dependency-free testthat-compatible runner for envs/bambu (no global package edits)
testthat_env <- new.env(parent = globalenv())
testthat_results <- list()
current_ctx <- ""

context <- function(desc) {
    current_ctx <<- desc
}

test_path <- function(...) {
    base_dir <- Sys.getenv("BAMBU_TEST_DIR", unset = getwd())
    file.path(base_dir, ...)
}

test_that <- function(desc, code) {
    cat(sprintf("  [TEST] %s ... ", desc))
    expr <- substitute(code)
    t0 <- proc.time()[["elapsed"]]
    status <- "PASS"
    err_msg <- ""
    env <- new.env(parent = parent.frame())
    tryCatch(
        eval(expr, envir = env),
        error = function(e) {
            status <<- "FAIL"
            err_msg <<- conditionMessage(e)
        }
    )
    dt <- proc.time()[["elapsed"]] - t0
    cat(sprintf("%s (%.2fs)\n", status, dt))
    if (status != "PASS") cat(sprintf("         REASON: %s\n", err_msg))
    testthat_results[[length(testthat_results) + 1L]] <<- data.frame(
        context = current_ctx, test = desc, status = status, elapsed = dt, error = err_msg,
        stringsAsFactors = FALSE
    )
}

expect_true <- function(object, info = NULL) {
    if (!isTRUE(object)) stop("expect_true failed: ", paste(deparse(object), collapse = " "), " ", info)
    invisible(TRUE)
}

expect_false <- function(object) {
    if (!identical(as.logical(object), FALSE)) stop("expect_false failed")
    invisible(TRUE)
}

expect_null <- function(object) {
    if (!is.null(object)) stop("expect_null failed: class=", paste(class(object), collapse = "/"))
    invisible(TRUE)
}

expect_equal <- function(object, expected, tolerance = sqrt(.Machine$double.eps), ...) {
    if (is(object, "SummarizedExperiment") && is(expected, "SummarizedExperiment")) {
        res <- all.equal(object, expected, tolerance = tolerance)
        if (!isTRUE(res)) stop("SE expect_equal failed: ", paste(head(res, 10), collapse = " | "))
        return(invisible(TRUE))
    }
    res <- all.equal(object, expected, tolerance = tolerance)
    if (!isTRUE(res)) stop("expect_equal failed: ", paste(head(res, 10), collapse = " | "))
    invisible(TRUE)
}

expect_identical <- function(object, expected) {
    if (!identical(object, expected)) stop("expect_identical failed")
    invisible(TRUE)
}

expect_s4_class <- function(object, class) {
    if (!is(object, class)) stop("expect_s4_class failed: expected ", class, ", got ", paste(class(object), collapse="/"))
    invisible(TRUE)
}

expect_is <- function(object, class) {
    if (!inherits(object, class) && !is(object, class)) stop("expect_is failed: expected ", class, ", got ", paste(class(object), collapse="/"))
    invisible(TRUE)
}

expect_type <- function(object, type) {
    if (!identical(typeof(object), type)) stop("expect_type failed: expected ", type, ", got ", typeof(object))
    invisible(TRUE)
}

expect_named <- function(object, expected) {
    nms <- names(object)
    if (!identical(nms, expected)) stop("expect_named failed: got [", paste(nms, collapse=","), "] vs [", paste(expected, collapse=","), "]")
    invisible(TRUE)
}

run_test_suite <- function(test_dir, skip_files = c("test_plotBambu.R")) {
    old_td <- Sys.getenv("BAMBU_TEST_DIR", unset = NA)
    Sys.setenv(BAMBU_TEST_DIR = normalizePath(test_dir))
    on.exit({
        if (is.na(old_td)) Sys.unsetenv("BAMBU_TEST_DIR") else Sys.setenv(BAMBU_TEST_DIR = old_td)
    }, add = TRUE)

    files <- sort(list.files(test_dir, pattern = "^test_.*\\.R$", full.names = TRUE))
    for (f in files) {
        bn <- basename(f)
        if (bn %in% skip_files) {
            cat(sprintf("=== SKIPPING %s (Suggests pkg optional) ===\n", bn))
            next
        }
        cat(sprintf("=== RUNNING %s ===\n", bn))
        current_ctx <<- bn
        e <- new.env(parent = asNamespace("bambu"))
        sys.source(f, envir = e)
    }
    df <- do.call(rbind, testthat_results)
    print(df)
    n_fail <- sum(df$status != "PASS")
    cat(sprintf("\nSUMMARY: total=%d passed=%d failed=%d\n", nrow(df), sum(df$status == "PASS"), n_fail))
    invisible(df)
}
