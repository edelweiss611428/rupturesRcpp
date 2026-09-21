#' @noRd
.segEstimate <- function(spec, tsMat, covariates, a, b) {
    Y <- tsMat
    if (is.null(colnames(Y))) colnames(Y) <- paste0("X", seq_len(ncol(Y)))
    n <- b - a
    idx <- (a + 1L):b

    switch(spec$costFunc,
        L1 = {
            Yt <- Y[idx, , drop = FALSE]
            med <- apply(Yt, 2, median)
            list(median = med, scale = colMeans(abs(sweep(Yt, 2, med))))
        },
        L2 = {
            Yt <- Y[idx, , drop = FALSE]
            m <- colMeans(Yt)
            list(mean = m, var = colMeans(sweep(Yt, 2, m)^2))
        },
        SIGMA = {
            Yt <- Y[idx, , drop = FALSE]
            m <- colMeans(Yt)
            list(mean = m, cov = crossprod(sweep(Yt, 2, m)) / n)
        },
        VAR = {
            r <- spec$pVAR
            p <- ncol(Y)
            J <- 1L + r * p
            lagNames <- paste0(rep(colnames(Y), r), ".lag", rep(seq_len(r), each = p))
            if (n < J) {
                return(list(
                    coef = matrix(NA_real_, J, p, dimnames = list(c("(Intercept)", lagNames), colnames(Y))),
                    var = setNames(rep(NA_real_, p), colnames(Y))
                ))
            }
            t <- (max(a, r) + 1L):b # rows whose r lags exist (lags may precede a)
            Z <- cbind(1, do.call(cbind, lapply(seq_len(r), function(j) Y[t - j, , drop = FALSE])))
            colnames(Z) <- c("(Intercept)", lagNames)
            Yt <- Y[t, , drop = FALSE]
            B <- qr.coef(qr(Z), Yt) # NA for aliased columns, never errors
            list(coef = B, var = colMeans((Yt - Z %*% B)^2))
        },
        LinearL2 = {
            Yt <- Y[idx, , drop = FALSE]
            if (is.null(covariates)) { # force-fit path of $fit(): intercept only
                Z <- matrix(1, n, 1L, dimnames = list(NULL, "(Intercept)"))
            } else {
                X <- covariates
                if (is.null(colnames(X))) colnames(X) <- paste0("Z", seq_len(ncol(X)))
                X <- X[idx, , drop = FALSE]
                Z <- if (spec$intercept) cbind("(Intercept)" = 1, X) else X
            }
            J <- ncol(Z)
            if (n < J) {
                return(list(
                    coef = matrix(NA_real_, J, ncol(Y), dimnames = list(colnames(Z), colnames(Y))),
                    var = setNames(rep(NA_real_, ncol(Y)), colnames(Y))
                ))
            }
            B <- qr.coef(qr(Z), Yt)
            list(coef = B, var = colMeans((Yt - Z %*% B)^2))
        }
    )
}

# ========================================================
#          Internal: build a segSummary from endPts
# ========================================================

#' @noRd
.segSummary <- function(evalFn, spec, tsMat, covariates, endPts) {

    nObs <- nrow(tsMat)

    if (!is.numeric(endPts)) {
        stop("`endPts` must be an integer vector specifying endpoints! Could be obtained via `$predict()`.")
    }
    endPts <- as.integer(sort(endPts))
    if (min(endPts) < 1L) stop("`min(endPts)` must not be less than 1!")
    if (max(endPts) != nObs) stop("By construction, `max(endPts)` must be `n`! Can use `$predict()` to obtain `endPts`.")
    if (anyDuplicated(endPts)) stop("`as.integer(endPts)` contains duplicated elements!")

    start <- c(0L, endPts[-length(endPts)])

    structure(
        Map(function(a, b) c(list(start = a, end = b, n = b - a, cost = evalFn(a, b)),
                             .segEstimate(spec, tsMat, covariates, a, b)),
            start, endPts),
        class = "segSummary", costFunc = spec)
}

# ========================================================
#                    segSummary S3 class
# ========================================================

#' Segment summary object
#'
#' @description A list with one element per segment, produced by `$summary()` of `PELT`, `binSeg`, `Window`
#' and `costFactory`. Each element is a named list with `start` (exclusive, 0-based), `end` (inclusive),
#' `n` (number of observations), `cost`, followed by the fitted parameters of the cost function:
#'
#' - `"L1"`: `median` (coordinate-wise median), `scale` (mean absolute deviation from the median).
#' - `"L2"`: `mean`, `var` (per-coordinate residual variance, divisor `n`).
#' - `"SIGMA"`: `mean`, `cov` (empirical covariance, divisor `n`, without `epsilon`).
#' - `"VAR"`: `coef` (rows `(Intercept)`, `X1.lag1`, ..., one column per response), `var`.
#' - `"LinearL2"`: `coef` (rows `(Intercept)` and covariates, one column per response), `var`.
#'
#' Dispersion entries use the maximum-likelihood divisor so that, e.g., `sum(var) * n` equals `cost` for `"L2"`.
#' Coefficients are `NA` where the segment is too short (`n < J`) or the system is rank deficient (there the `C++`
#' cost uses a pseudo-inverse). The `costFunc` configuration is stored as `attr(x, "costFunc")`.
#'
#' @param x A `segSummary` object.
#' @param digits Integer. Significant digits for printing. Default: `max(3L, getOption("digits") - 3L)`.
#' @param row.names,optional Ignored; present for compatibility with the generic.
#' @param ... Ignored.
#'
#' @return `print()` invisibly returns `x`. `as.data.frame()` returns a `data.frame` with one row per segment:
#' `segment`, `start`, `end`, `n`, `cost`, then flattened parameters (`mean.X1`, `coef.(Intercept).X1`, ...).
#'
#' @examples
#' set.seed(1)
#' tsMat = as.matrix(c(rnorm(100, 0), rnorm(100, 5)))
#' obj = binSeg$new(); obj$fit(tsMat); obj$predict(pen = 50)
#' s = obj$summary()
#' s
#' as.data.frame(s)
#' s[[2]]$mean
#' @name segSummary
NULL

#' @rdname segSummary
#' @importFrom utils capture.output
#' @importFrom stats median setNames
#' @export
print.segSummary <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

    spec <- attr(x, "costFunc")
    pars <- if (length(spec) > 1L) {
        paste0(" (", paste(names(spec)[-1L], spec[-1L], sep = " = ", collapse = ", "), ")")
    } else ""

    cat("Segmentation summary\n")
    cat(sprintf("costFunc     : \"%s\"%s\n", spec$costFunc, pars))
    cat(sprintf("segments     : %d (%d change-point%s)\n",
                length(x), length(x) - 1L, if (length(x) == 2L) "" else "s"))

    for (i in seq_along(x)) {
        s <- x[[i]]
        cat(sprintf("\nSegment %d: (%d, %d]    n = %d    cost = %s\n",
                    i, s$start, s$end, s$n, format(s$cost, digits = digits)))
        for (nm in setdiff(names(s), c("start", "end", "n", "cost"))) {
            cat("  ", nm, "\n", sep = "")
            writeLines(paste0("  ", capture.output(print(s[[nm]], digits = digits))))
        }
    }
    invisible(x)
}

#' @rdname segSummary
#' @export
as.data.frame.segSummary <- function(x, row.names = NULL, optional = FALSE, ...) {

    .flat <- function(nm, v) {
        if (is.matrix(v)) {
            setNames(c(v), paste(nm, rownames(v)[row(v)], colnames(v)[col(v)], sep = "."))
        } else {
            setNames(v, paste(nm, names(v), sep = "."))
        }
    }

    rows <- lapply(seq_along(x), function(i) {
        s <- x[[i]]
        pars <- unlist(lapply(setdiff(names(s), c("start", "end", "n", "cost")), function(nm) .flat(nm, s[[nm]])))
        data.frame(segment = i, start = s$start, end = s$end, n = s$n, cost = s$cost, t(pars), check.names = FALSE)
    })
    do.call(rbind, rows)
}
