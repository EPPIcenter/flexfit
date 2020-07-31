
locatePt <- function(ptx, pty, datx, daty) {
  absdif <- abs(ptx - datx)
  indmin <- which(absdif == min(absdif))       # multiple
  indmin[which.min(abs(pty - daty[indmin]))]   # single
}

#' Plot the fit for standards and the samples
#'
#' Produces a plot that includes points for standards, proposed fit, removed
#' outliers, bounds for "flat" portions of the curve, and values for samples and
#' for the background.
#'
#' @details to be added
#'
#' @param fitpar  values of function parameters.
#' @param FUNmod  model function.
#' @param iout    indices of removed standard points.
#' @param bg      background values.
#' @param smpflag character vector, flags for each sample.
#' @param trimval for final results, the values at which the samples are
#'   trimmed.
#' @param trimext integer vector of length two indicating if the values are
#'   trimmed at the extremum (lower and upper).
#' @inheritParams fitStd
#' @inheritParams processSmp
#'
#' @export

plotFit <- function(std, xvar, yvar, dilvar, fitpar = NULL, FUNmod = NULL,
                    iout = NULL, bg = NULL, vsmp = NULL, smpflag = NULL,
                    trimval = NULL, trimext = NULL,
                    stdcol = c("firebrick3", "darkslategray"),
                    rugcol = c("cadetblue", "purple", "firebrick2"), ...) {
  if (!hasArg(ylim)) {
    ylim <- range(std[, yvar], vsmp, bg, na.rm = TRUE)
    plot(std[, xvar], std[, yvar], col = stdcol[1], xaxt = "n", ylim = ylim,
         lwd = 1.3, ...)
  } else {
    plot(std[, xvar], std[, yvar], col = stdcol[1], xaxt = "n", lwd = 1.3, ...)
  }
  if(!is.null(vsmp)) {
    rug(vsmp, side = 2, col = rugcol[1])
    if (!is.null(smpflag)) {
      rug(vsmp[grep("lower|upper", smpflag)], side = 2, col = rugcol[2])
      rug(vsmp[grep("min|max",     smpflag)], side = 2, col = rugcol[3])
    }
  }

  if (dilvar %in% colnames(std)) {
    labs = parse(text = paste("frac(1, ", std[, dilvar], ")", sep = ""))
  } else {
    labs = round(std[, xvar], 3)
  }
  axis(side = 1, at = std[, xvar], cex.axis = 0.7, tcl = -0.1, labels = labs)
  abline(h = bg, lty = 3)
  if (!is.null(iout)) {
    points(std[iout, xvar], std[iout, yvar], col = 2, pch = 4, cex = 2)
  }
  if (!is.null(trimval)) {
    abline(h = trimval, col = rugcol[2:3][trimext], lty = 6, lwd = 1.2)
    legend("right", inset = 0.03, box.col = "white", bg = "white", cex = 0.9,
           col = rugcol[3:2], lty = 6, lwd = 1.5, seg.len = 2.5,
           title = "trimmed at", legend = c("extrema", "bounds"))
  }
  if (is.null(fitpar)) {
    legend("bottom", inset = 0.03, box.col = "grey", box.lwd = 0.8,
           bg = "white", cex = 0.9, col = c(stdcol[1], 1), lty = c(NA, 3),
           lwd = c(1.5, 1), pch = c(1, NA), seg.len = 2.5,
           legend = c("standards", "background"))
  } else {
    x <- seq(min(std[, xvar]), max(std[, xvar]), length = 200)
    y <- FUNmod(x, fitpar)
    lines(x, y, lty = 5, lwd = 1.8, col = stdcol[2])
    legend("bottom", inset = 0.05, box.col = "grey", box.lwd = 0.8,
           bg = "white", cex = 0.9, col = c(stdcol[1], 1, stdcol[2]),
           lty = c(NA, 3, 5), lwd = c(1.5, 1, 2), pch = c(1, NA, NA),
           seg.len = 2.5, legend = c("standards", "background", "fit"))
  }
}
