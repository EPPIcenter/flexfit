#' Process raw Luminex files
#'
#' Process data for a single antigen: fit a standard curve, establish bounds,
#' normalize samples, and save a plot showing the fit and the samples.
#'
#' @details Note that if \code{fitlog} contains _x_ and thus \code{xvar} in
#'   standards is log-transform before fitting, the results are still returned
#'   on a regular scale (values of inverse function are exponentiated and then
#'   multiplied by dilution).
#'
#' @param smp data frame with sample ID, sample values and other optional
#'   variables (e.g. dilution).
#' @param std matrix or data frame with standards for fitting.
#' @param bg  values for background spots.
#' @param smpdil single value for sample dilutions (if dilutions are not
#'   provided in \code{smp} data frame). Ignored if \code{dilvar} is provided
#'   and the variable is included in \code{smp} data frame.
#' @param fitlog character string indicating if standard values should be
#'   log-transformed for fitting. If the string contains _x_, \code{xvar} will
#'   be transformed, if it contains _y_ - \code{yvar}.
#' @param ismp indices for samples to be plotted (e.g. to exclude standards that
#'   are included in \code{smp} data frame).
#' @param plotdir directory for the plots to be saved.
#' @param pname  character string for fit plot. If \code{NULL}, will be based on
#'   \code{ptitle} with underscore substituted for space.
#' @param ptitle character string for plot title.
#' @param xvar,yvar character strings for the variables used to fit a standard
#'   curve. If \code{NULL}, first two columns are assumed to be \code{x} and
#'   \code{y} variables.
#' @param smpvar character string indicating sample variable.
#' @param addvar named vector for additional variables (e.g. date or plate),
#'   where element names will be used as variable names in \code{smp}.
#' @param dilvar,resvar character strings for dilution variable and results.
#' @param model the model to be fit.
#' @param Alow lower asymptote for the sigmoid model. If \code{NULL}, the lower
#'   asymptote will be estimated, adding an extra parameter to the model. To fix
#'   the asymptote at the level of background, specify \code{"bg"}. Numeric
#'   value of \code{Alow} will force the asymptote to be fixed at the provided
#'   level.
#' @param asym if \code{TRUE}, asymmetry in the fit is allowed, adding an extra
#'   parameter to the model.
#' @param trim.flat logical value determining how the values of \code{yvar} are
#'   trimmed. If \code{TRUE}, they will be trimmed at the bounds where the curve
#'   starts to flatten out (automatically determined as maxima of the third
#'   derivative of the function). If \code{FALSE}, \code{yvar} will be trimmed
#'   at extrema, defined as the range of standards or asymptotes of the fit
#'   (whichever are less extreme).
#' @param interactive logical value. If \code{TRUE}, the user is prompted to
#'   evaluate the standards (and/or the fit) and asked to remove outliers if
#'   needed. \code{TRUE} value takes precedence over \code{rm.before} and
#'   \code{rm.after}: if both are \code{FALSE}, \code{rm.after} is reset to
#'   \code{TRUE}.
#' @param monot.prompt if \code{TRUE}, the user is prompted to evaluate the
#'   standards and possibly remove outliers if the standards are not monotonic
#'   (increasing). \code{FALSE} value is ignored if \code{interactive} is
#'   \code{TRUE}.
#' @param rm.before logical value indicating if potential outliers should be
#'   removed before the model is fitted. Ignored if \code{interactive} is
#'   \code{FALSE}.
#' @param rm.after logical value indicating if potential outliers should be
#'   removed after the model is fitted. Ignored if \code{interactive} is
#'   \code{FALSE}.
#' @param maxrm maximum number of outliers to remove.
#' @param set.bounds if \code{TRUE}, the user is given the option to manually
#'   set the bound that is not set automatically. In that case, the prompt
#'   appears even if \code{interactive} is \code{FALSE}.
#' @param overwrite.bounds logical value indicating the option to overwrite
#'   automatically set bounds. Ignored if \code{interactive} is \code{FALSE}.
#' @param ifix sorted integer vector of length 3 with indices of standards to be
#'   used for getting starting values for optimization.
#'
#' @param extrapolate.low if \code{TRUE}, sample values beyond lower bounds will
#'   be processed by extrapolation of the standard curve (not recommended).
#'   Takes precedence over \code{trim.flat} value.
#' @param extrapolate.up if \code{TRUE}, sample values beyond upper bounds will
#'   be processed by extrapolation of the standard curve (not recommended).
#'   Takes precedence over \code{trim.flat} value.

#' @param optmethod method to be used in optimization.
#' @param maxit maximum number of iterations in optimization.
#' @param stdcol vector of two colors for the standard points and the fit on the
#'   plot.
#' @param rugcol vector of three colors for the rugplot, which indicates sample
#'   values (inside the bounds, between the bounds and extrema, and beyond
#'   extrema).
#' @param width,height optional parameters for the final saved plot.
#' @param ... further graphical parameters.
#'
#' @return A list of length three. The first element is a data frame that
#'   contains the results; the second is a character string with a flag
#'   containing information about removed points, failure to fit the model,
#'   manually set bounds, and/or an optional custom note provided by the user
#'   during an interactive model-fitting procedure. The last element is the
#'   number of sample values for which the results are trimmed.
#'
#' @examples
#'
#' @importFrom grDevices pdf dev.off adjustcolor
#' @importFrom graphics abline axis legend lines locator mtext plot points rug
#' @importFrom methods hasArg
#' @importFrom stats optim
#' @importFrom utils read.csv
#'
#' @export

processSmp <- function(smp, std, bg = NULL, smpdil = 1, fitlog = "xy",
                       ismp = 1:nrow(smp),
                       plotdir = "./", pname = NULL, ptitle = "fit and samples",
                       xvar = NULL, yvar = NULL, smpvar = NULL, addvar = NULL,
                       dilvar = NULL, resvar = "concentration", #NULL,
                       model = "sigmoid", Alow = NULL, asym = TRUE,
                       trim.flat = TRUE, interactive = TRUE,
                       monot.prompt = FALSE,
                       rm.before = FALSE, rm.after = interactive, maxrm = 2,
                       set.bounds = interactive, overwrite.bounds = FALSE,
                       ifix = NULL,
                       extrapolate.low = FALSE, extrapolate.up = FALSE,
                       optmethod = "Nelder-Mead", maxit = 5e3,
                       # nv3 = 10, nv4 = 100, # then include in fitStd() as well
                       stdcol = c("firebrick3", "darkslategray"),
                       rugcol = c("cadetblue", "purple", "firebrick2"),
                       width = 7, height = 6, ...) {
  if (!inherits(smp, "data.frame")) {
    smp <- as.data.frame(smp)
  }
  if (!(inherits(std, "matrix") || inherits(std, "data.frame"))) {
    std <- as.data.frame(std)
  }
  std <- std[order(std[, xvar]), ]       # for iout; also sorted in fitStd()
  if (is.null(xvar) || is.null(yvar)) {  # first two columns assumed to be x, y
    colnames(std) <- c("x", "y")
    xvar <- "x"
    yvar <- "y"
  }
  if (is.null(smpvar)) {                 # 2nd column assumed to be samples
    smpvar           <- "Sample"
    colnames(smp)[2] <- "Sample"
  }
  if (!is.null(addvar)) {                # variables added
    for (i in 1:length(addvar)) {
      smp[, names(addvar[i])] <- addvar[i]
    }
  }
  if (is.null(dilvar)) {                 # can be removed after processing if 1
    dilvar <- "Dilution"
  }
  if (!dilvar %in% colnames(smp)) {
    smp[, dilvar] <- smpdil
  }
  xlab <- xvar
  ylab <- yvar
  if (dilvar %in% colnames(std)) {
    tcklab <- parse(text = paste("frac(1, ", std[, dilvar], ")", sep = ""))
  } else {
#    tcklab <- round(std[, xvar], 3)  #*** low conc shown as 0's
#    tcklab <- std[, xvar]
    tcklab <- sprintf("%.3e", std[, xvar])
  }
  if (grepl("x", fitlog)) {
    std[, xvar] <- log(std[, xvar])
    xlab <- paste(xvar, "(log scale)")
  }
  if (grepl("y", fitlog)) {
    oldvar <- smpvar                     # reassign smpvar
    smpvar <- paste("log", smpvar, sep = "")
    std[, yvar]   <- log(std[, yvar])
    smp[, smpvar] <- log(smp[, oldvar])
    if (!is.null(bg)) bg <- log(bg)
    ylab <- paste("log", yvar)
  }
  if (is.null(pname)) {                  # plot title used for pdf name
    pname <- gsub("[[:space:]]", "_", ptitle)
  }

  finfit <- fitStd(std, xvar, yvar, model, Alow, asym,
                   interactive, monot.prompt, rm.before, rm.after, maxrm,
                   set.bounds, overwrite.bounds, bg, smp[ismp, smpvar],
                   optmethod, maxit,
                   extrapolate.low = extrapolate.low,
                   extrapolate.up  = extrapolate.up,
                   info = ptitle, ifix = ifix,
                   tcklab = tcklab, stdcol = stdcol, rugcol = rugcol,
                   main = ptitle, xlab = xlab, ylab = ylab, ...)
  if (grepl("sig", model)) {
    FUNmod <- fsig
    FUNinv <- fsigInv
  } else {
    FUNmod <- flin
    FUNinv <- flinInv
  }

  dfout <- normalizeSmp(smp, smpvar, resvar, dilvar, FUNinv, finfit$par,
                        finfit$bounds, finfit$flag, fitlog, trim.flat,
                        extrapolate.low, extrapolate.up)

  fplot <- paste(plotdir, pname, ".pdf", sep = "")
  pdf(file = fplot, width = width, height = height)
  if (!is.null(finfit$par)) {
    trimval <- dfout[1, c("trim_lo", "trim_up")]
#    trimext <- (trimval == dfout[1, c("min", "max")]) + 1 # at extrema (min/max)
    a <- trimval == dfout[1, paste0(c("lower", "upper"), "_bound")]
    a[is.na(a)] <- FALSE
    b <- trimval == dfout[1, c("min", "max")]
    # 1: bounds, 2: min/max, 3: asymptote
    trimtype <- (!a & !b) + 2 - a
    plotFit(std, xvar, yvar, finfit$par, FUNmod, FUNinv, finfit$iout, bg,
            dfout[ismp, smpvar], dfout[ismp, "Flag"], trimval, trimtype,
            extrapolate.low = extrapolate.low,
            extrapolate.up  = extrapolate.up,
            tcklab = tcklab, stdcol = stdcol, rugcol = rugcol,
            main = ptitle, xlab = xlab, ylab = ylab, ...)
  } else {
    plotFit(std, xvar, yvar, iout = finfit$iout, bg = bg,
            vsmp = dfout[ismp, smpvar],
            extrapolate.low = extrapolate.low,
            extrapolate.up  = extrapolate.up,
            tcklab = tcklab, stdcol = stdcol,
            rugcol = rugcol, main = ptitle, xlab = xlab, ylab = ylab, ...)
  }
  dev.off()
  #  return(list(smp = dfout, fitflag = finfit$flag))  #*** if not below
  #***================ optional: number of trimmed samples =================***#
  return(list(smp = dfout, fitflag = finfit$flag, ntrim = sum(dfout$trimmed)))
}


