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
#'   provided in \code{smp} data frame).
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
#'   curve.
#' @param smpvar character string indicating sample variable.
#' @param addvar named vector for additional variables (e.g. date or plate),
#'   where element names will be used as variable names in \code{smp}.
#' @param dilvar,resvar character strings for dilution variable and results
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

#' @param optmethod method to be used in optimization.
#' @param maxit maximum number of iterations in optimization.
#' @param rugcol vector of three colors for the rugplot, which indicates sample
#'   values (inside the bounds, between the bounds and extrema, and beyond
#'   extrema).
#' @param stdcol vector of two colors for the standard points and the fit on the
#'   plot.
#' @param width,height optional parameters for the final saved plot.
#' @param ... further graphical parameters.
#'
#' @return A list of length three. The first element is a data frame that
#'   contains the results; the second is a character string with a flag
#'   containing information about removed points, failure to fit the model,
#'   manually set bounds, and/or an optional custom note provided by the user
#'   during an interactive model-fitting procedure. ***optional return of the
#'   number of trimmed samples (then a list of length three)***
#'
#' @examples
#'
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline axis legend lines locator mtext plot points rug
#' @importFrom methods hasArg
#' @importFrom stats optim
#' @importFrom utils read.csv
#'
#' @export

# "coral3", "chartreuse3", "aquamarine3", "mediumvioletred", "olivedrab4" and 3
# "seagreen4" and 3, "springreen4" and 3, "darkslategray4" and "" (!)
processSmp <- function(smp, std, bg = NULL, smpdil = 1, fitlog = "xy",
                       ismp = 1:nrow(smp),
                       plotdir = "./", pname = NULL, ptitle = "fit and samples",
                       xvar = NULL, yvar = NULL, smpvar = NULL, addvar = NULL,
                       dilvar = NULL, resvar = NULL,
                       model = "sigmoid", Alow = NULL, asym = TRUE,
                       trim.flat = TRUE, interactive = TRUE,
                       monot.prompt = FALSE,
                       rm.before = FALSE, rm.after = interactive, maxrm = 2,
                       set.bounds = interactive, overwrite.bounds = FALSE,
                       ifix = NULL, optmethod = "Nelder-Mead", maxit = 5e3,
                       # nv3 = 10, nv4 = 100, # then include in fitStd() as well
                       rugcol = c("cadetblue", "darkorchid3", "firebrick3"),
                       stdcol = c("aquamarine4", "aquamarine3"),
                       width = 7, height = 6, ...) {
  options(warn = 1)    # for interactivity
  if (inherits(smp, "matrix")) {
    smp <- as.data.frame(smp)
  }
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
    smp[, dilvar] <- smpdil
  }
  if (is.null(pname)) {                  # plot title used for pdf name
    pname <- gsub("[[:space:]]", "_", ptitle)
  }
  xlab <- xvar
  ylab <- yvar
  if (grepl("x", fitlog)) {
    std[, xvar] <- log(std[, xvar])
    xlab <- paste("log", xvar)
  }
  if (grepl("y", fitlog)) {
    smpvar <- paste("log", smpvar, sep = "")
    std[, yvar]   <- log(std[, yvar])
    smp[, smpvar] <- log(smp[, smpvar])  # reassigned smpvar
    ylab <- paste("log", yvar)
  }
  std <- std[order(std[, xvar]), ]       # sorted by increasing xvar

  finfit <- fitStd(std, xvar, yvar, dilvar, model, Alow, asym, interactive,
                   monot.prompt, rm.before, rm.after, maxrm, set.bounds,
                   overwrite.bounds, bg, smp[ismp, smpvar], optmethod, maxit,
                   info = ptitle, ifix = ifix, rugcol = rugcol, stdcol = stdcol,
                   main = ptitle, xlab = xlab, ylab = ylab, ...)
  if (grepl("sig", model)) {
    FUNmod <- fsig
    FUNinv <- fsigInv
  } else {
    FUNmod <- flin
    FUNinv <- flinInv
  }

  smps <- normalizeSmp(smp, smpvar, resvar, dilvar, FUNinv, finfit$par,
                       finfit$bounds, finfit$flag, fitlog, trim.flat)

  fplot <- paste(plotdir, pname, ".pdf", sep = "")
  pdf(file = fplot, width = width, height = height)
  if (!is.null(finfit$par)) {
    trimval <- smps[1, c("trim_lo", "trim_up")]
    trimext <- (trimval == smps[1, c("min", "max")]) + 1  # at extrema
    plotFit(std, xvar, yvar, dilvar, finfit$par, FUNmod, finfit$iout, bg,
            smp[, smpvar], smp[, "Flag"], trimval, trimext, rugcol, stdcol,
            main = ptitle, xlab = xlab, ylab = ylab, ...)
  } else {
    plotFit(std, xvar, yvar, dilvar, iout = finfit$out, bg = bg,
            smp = smp[, smpvar], rugcol = rugcol, stdcol = stdcol,
            main = ptitle, xlab = xlab, ylab = ylab, ...)
  }
  dev.off()
  options(warn = 0)
  #  return(list(smps = smps, fitflag = finfit$flag))  #*** if not below
  #***================ optional: number of trimmed samples =================***#
  return(list(smp = smp, fitflag = finfit$flag, ntrim = sum(smp$trimmed)))
}


