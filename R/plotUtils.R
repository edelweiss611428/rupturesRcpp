# Internal helpers shared by the `$plot()` methods of `PELT`, `binSeg`, `Window` and `Dynp`.
# None of these are exported.

.plotTypes = c("series", "fitted", "residuals", "diagnostics")

.checkPlotType = function(type) {

  if (is.null(type) || !is.character(type) || length(type) != 1L || is.na(type) ||
      !(type %in% .plotTypes)) {
    stop("`type` must be one of \"series\", \"fitted\", \"residuals\" or \"diagnostics\"!")
  }

  type

}

# Fitted values implied by the parameter estimates of each segment (a, b] given by `endPts`:
# the segment mean for "L2" and "SIGMA", the coordinate-wise median for "L1", and X %*% coef for
# the regression costs "LinearL2", "LinearSIGMA", "LinearL1" and "VAR" (whose regressors are an
# intercept and the lagged series).
#
# `getParams(a, b)` must return the C++ module's `get_params(a, b)`. Returns an n x p matrix with
# NA where a segment has too few observations to estimate its parameters and, for "VAR", on the
# first pVAR observations, which have no lags.

.segmentFitted = function(costPass, tsMat, covariates, endPts, getParams) {

  costName = costPass[["costFunc"]]

  if (costName == "Custom") {
    stop("Fitted values and residuals are only available for the built-in cost functions, not \"Custom\"!")
  }

  n = nrow(tsMat)
  p = ncol(tsMat)

  # Regressors of each observation, with columns in the same order as the rows of `coef`
  design = NULL

  if (costName %in% c("LinearL2", "LinearSIGMA", "LinearL1")) {

    if (is.null(covariates)) {
      # `$fit()` without covariates force-fits an intercept-only model
      design = matrix(1, nrow = n, ncol = 1)
    } else {
      if (nrow(covariates) != n) {
        stop("Numbers of observations in `covariates` and `tsMat` do not match!")
      }
      design = if (isTRUE(costPass[["intercept"]])) cbind(1, covariates) else covariates
    }

  } else if (costName == "VAR") {

    pVAR = costPass[["pVAR"]]
    design = matrix(NA_real_, nrow = n, ncol = 1 + pVAR * p)

    if (n > pVAR) {
      rows = (pVAR + 1):n
      design[rows, 1] = 1
      for (L in seq_len(pVAR)) {
        design[rows, 1 + (L - 1) * p + seq_len(p)] = tsMat[rows - L, , drop = FALSE]
      }
    }

  }

  fitted = matrix(NA_real_, nrow = n, ncol = p)
  starts = c(0L, endPts[-length(endPts)])

  for (k in seq_along(endPts)) {

    rows = (starts[k] + 1):endPts[k]
    params = getParams(starts[k], endPts[k])

    if (costName %in% c("L2", "SIGMA")) {
      fitted[rows, ] = matrix(params[["mean"]], nrow = length(rows), ncol = p, byrow = TRUE)

    } else if (costName == "L1") {
      fitted[rows, ] = matrix(params[["median"]], nrow = length(rows), ncol = p, byrow = TRUE)

    } else {
      coefMat = params[["coef"]]
      if (nrow(coefMat) != ncol(design)) {
        stop("`covariates` do not match the fitted model! Re-run `$fit()` first.")
      }
      fitted[rows, ] = design[rows, , drop = FALSE] %*% coefMat
    }

  }

  fitted

}

# Draws `$plot()`. For `type = "series"` this is the raw series with the segments shaded and the
# change-points marked. The other types need `fitted`, the n x p matrix from `.segmentFitted()`.

.plotSegmentation = function(tsMat, endPts, d, dimNames, main, xlab, tsWidth, tsCol, bgCol,
                             bgAlpha, ncol, type = "series", fitted = NULL, fitCol = "#D55E00",
                             costLabel = NULL) {

  n = nrow(tsMat)
  nSeg = length(endPts)
  segment = rep(seq_len(nSeg), times = diff(c(0L, endPts)))
  segFill = rep(bgCol, length.out = nSeg)

  # Long-format data frame for all selected dimensions
  allTsDf = do.call(rbind, lapply(seq_along(d), function(i) {
    data.frame(
      time = 1:n,
      value = tsMat[, d[i]],
      dimension = dimNames[i],
      segment = segment,
      fill = segFill[segment],
      fitted = if (is.null(fitted)) NA_real_ else fitted[, d[i]]
    )
  }))
  allTsDf$residual = allTsDf$value - allTsDf$fitted

  # Segment info
  segInt = data.frame(
    xmin = c(1, endPts[-nSeg]),
    xmax = endPts,
    fill = segFill
  )

  plotTheme = theme_minimal() +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold")
    )

  # Shaded segments and dashed change-point lines, shared by the plots against time
  segShade = geom_rect(
    data = segInt,
    aes(xmin = .data[["xmin"]], xmax = .data[["xmax"]], ymin = -Inf, ymax = Inf, fill = .data[["fill"]]),
    inherit.aes = FALSE,
    alpha = bgAlpha
  )
  cpLines = geom_vline(xintercept = endPts[-nSeg], linetype = "dashed", color = "black", linewidth = tsWidth)

  if (type == "series") {

    return(
      ggplot(allTsDf, aes(x = .data[["time"]], y = .data[["value"]])) +
        scale_fill_identity() +
        segShade +
        geom_line(color = tsCol, linewidth = tsWidth) +
        cpLines +
        facet_wrap(~ dimension, scales = "free_y", ncol = ncol) +
        plotTheme +
        labs(x = xlab, y = NULL, title = main)
    )

  }

  fitDf = allTsDf[!is.na(allTsDf$fitted), , drop = FALSE]

  if (nrow(fitDf) == 0L) {
    stop("No fitted values to plot: every segment has too few observations to estimate its parameters!")
  }

  noFit = setdiff(seq_len(nSeg), unique(fitDf$segment))
  if (length(noFit) > 0L) {
    message(sprintf("Segment(s) %s have too few observations to estimate their parameters and are left out.",
                    toString(noFit)))
  }

  if (type == "fitted") {

    return(
      ggplot(allTsDf, aes(x = .data[["time"]], y = .data[["value"]])) +
        scale_fill_identity() +
        segShade +
        geom_line(color = tsCol, linewidth = tsWidth) +
        # One line per segment, so the fit is not joined across change-points
        geom_line(data = fitDf, aes(y = .data[["fitted"]], group = .data[["segment"]]),
                  color = fitCol, linewidth = 2 * tsWidth) +
        cpLines +
        facet_wrap(~ dimension, scales = "free_y", ncol = ncol) +
        plotTheme +
        labs(x = xlab, y = NULL, title = main, subtitle = paste0("Fitted values (", costLabel, ")"))
    )

  }

  if (type == "residuals") {

    return(
      ggplot(fitDf, aes(x = .data[["time"]], y = .data[["residual"]], group = .data[["segment"]])) +
        scale_fill_identity() +
        segShade +
        geom_hline(yintercept = 0, color = fitCol, linewidth = tsWidth) +
        geom_line(color = tsCol, linewidth = tsWidth) +
        cpLines +
        facet_wrap(~ dimension, scales = "free_y", ncol = ncol) +
        plotTheme +
        labs(x = xlab, y = NULL, title = main, subtitle = paste0("Residuals (", costLabel, ")"))
    )

  }

  # type == "diagnostics": residuals against fitted values, and the residual spread of each segment
  zeroLine = geom_hline(yintercept = 0, color = fitCol, linewidth = tsWidth)

  pFitted = ggplot(fitDf, aes(x = .data[["fitted"]], y = .data[["residual"]])) +
    scale_fill_identity() +
    zeroLine +
    geom_point(aes(fill = .data[["fill"]]), shape = 21, color = "grey30", stroke = 0.25, size = 1.5) +
    facet_wrap(~ dimension, scales = "free", ncol = ncol) +
    plotTheme +
    labs(x = "Fitted value", y = "Residual")

  fitDf$segment = factor(fitDf$segment, levels = seq_len(nSeg))

  pSpread = ggplot(fitDf, aes(x = .data[["segment"]], y = .data[["residual"]])) +
    scale_fill_identity() +
    zeroLine +
    geom_boxplot(aes(fill = .data[["fill"]]), alpha = bgAlpha, color = "grey30",
                 outlier.size = 0.75, linewidth = tsWidth) +
    # Keep at most about 10 segment labels on the axis
    scale_x_discrete(breaks = function(x) x[unique(round(seq(1, length(x), length.out = min(length(x), 10L))))]) +
    facet_wrap(~ dimension, scales = "free_y", ncol = ncol) +
    plotTheme +
    labs(x = "Segment", y = NULL)

  (pFitted | pSpread) +
    plot_annotation(title = main, subtitle = paste0("Residual diagnostics (", costLabel, ")"),
                    theme = plotTheme)

}

# Name of the cost function shown in the subtitle of the fitted, residual and diagnostic plots
.costLabel = function(costPass) {

  if (costPass[["costFunc"]] == "VAR") {
    return(paste0("VAR(", costPass[["pVAR"]], ")"))
  }

  costPass[["costFunc"]]

}
