#' Fit Standards
#'
#' Fit a specified function to standards (serial dilutions). Optionally an
#' interactive procedure that allows to remove outliers, evaluate resulting
#' fits, perform revisions, and record a message regarding the fit.
#'
#' @details to be added
#'
#' @param dilvar character string to check if there is a dilution variable in
#'   \code{std}. If found, used for x-axis labels only.
#' @param vsmp    sample values.
#' @param info   information about a particular run for warning messages.
#' @inheritParams processSmp
#'
#' @return A list containing parameters of the fit and bounds of the fit (named
#'   vectors), as well as indices of removed points (if any) and flags r.
#'
#' @examples
#'
#' @export

fitStd <- function(std, xvar, yvar, dilvar,
                   model = "sigmoid", Alow = NULL, asym = TRUE,
                   interactive = TRUE, monot.prompt = FALSE,
                   rm.before = FALSE, rm.after = interactive, maxrm = 2,
                   set.bounds = FALSE, overwrite.bounds = FALSE, bg = NULL,
                   vsmp = NULL, optmethod = "Nelder-Mead", maxit = 5e3,
                   info = "", ifix = NULL, stdcol, rugcol, ...) {
  if (!is.null(Alow) && Alow == "bg") Alow <- mean(bg)  # on a log scale
  flag <- ""
  iout <- NULL
  std <- std[order(std[, xvar]), ]         # sort std by increasing xvar
  monot <- all(diff(std[, yvar]) > -1e-9)
  if (!monot) warning(paste(info, yvar, "not strictly increasing in standards"))
  if (interactive) {
    if (!rm.before && !rm.after) {
      rm.after <- TRUE
    }
  } else {
    rm.before <- FALSE
    rm.after  <- FALSE
    overwrite.bounds <- FALSE
    if (monot.prompt && !monot) {
      rm.before <- TRUE  # or rm.after <- TRUE
    }
  }

  if (rm.before) {
    plotFit(std, xvar, yvar, dilvar, bg = bg, vsmp = vsmp,
            stdcol = stdcol, rugcol = rugcol, ...)
    for (i in 1:maxrm) {
      ans1 <- readline("Remove any outliers? (y/n) ")
      if (tolower(ans1) == "y"){
        flag <- "pts_rm"
        mtext("Click on an outlier", col = "red", cex = 1.2)
        out  <- locator(n = 1)
        iout <- c(iout, locatePt(out$x, out$y, std[, xvar], std[, yvar]))
        points(std[iout, xvar], std[iout, yvar], col = 2, pch = 4, cex = 2)
      } else {
        break
      }
    }
  }

  #--------------------------- get starting values ----------------------------#
  std1 <- std  # to account for possible outliers from rm.before
  if (!is.null(iout)) std1 <- std[-iout, ]
  if (grepl("sig", model)) {
    FUNmod <- fsig
    if (!is.null(Alow)) {  # lower asymptote fixed
      startval <- getStart3par(std1[, xvar], std1[, yvar], Alow, ifix = ifix)
      if (is.na(startval[1])) {
        warning(paste(info,
                      "Unable to fit the model with provided value for Alow"))
        ans <- readline("Estimate lower asymptote (recommended)? (y/n) ")
        if (grepl("y", tolower(ans))) {
          Alow <- NULL
        } else {
          startval <- NA
        }
      }
    }
    if (is.null(Alow)) {  # lower asymptote not fixed (originally or later)
      startval <- getStart4par(std1[, xvar], std1[, yvar], ifix = ifix,
                               info = info)
    }
  } else {  # model 5 (6 parameters) *** update when available
    FUNmod <- flin
  }  # model 5 (6 parameters)

  if (is.na(startval[1])) {
    flag <- "no_fit"
    warning(paste(info, "Unable to fit the model, manual check recommended"))
    fit <- NULL
  } else {
    if (grepl("sig", model)) {
      lsopt <- makeLsopt(std1[, xvar], std1[, yvar], Alow, asym)
    } else {}                   # for the other model
    #----------------------------- fit with optim() ---------------------------#
    if (asym) startval <- c(startval, a = 1)
    fit <- optim(startval, lsopt, method = optmethod,
                 control = list(maxit = maxit))
  }
  fitpar <- addParam(fit$par, Alow, asym)
  maxdet <- min(max(std1[, yvar]), fitpar["Aup"])   # extrema
  mindet <- max(min(std1[, yvar]), fitpar["Alow"])

  revise <- rm.after
  if (revise) {
    fit.orig  <- fitpar         # save original fitpar and iout from rm.before
    iout.orig <- iout           # if out at rm.before, that's it
  }

  while (revise) {
    flag <- ""
    fitpar  <- fit.orig         # restore original fit and iout
    iout    <- iout.orig
    #*** INSERTED - is it a useful functionality?
    smpflag <- rep("", length(vsmp))
    mindet <- max(min(std1[, yvar]), fitpar["Alow"])
    maxdet <- min(max(std1[, yvar]), fitpar["Aup"])
    smpflag[!(mindet <= vsmp & vsmp <= maxdet)] <- "min"
    plotFit(std, xvar, yvar, dilvar, fitpar = fitpar, FUNmod = FUNmod, bg = bg,
            vsmp = vsmp, smpflag = smpflag, stdcol = stdcol, rugcol = rugcol,
            ...)
    #*** end INSERTED, uncomment plotFit() below
#    plotFit(std, xvar, yvar, dilvar, fitpar = fitpar, FUNmod = FUNmod, bg = bg,
#            vsmp = vsmp, stdcol = stdcol, rugcol = rugcol, ...)

    for (i in 1:maxrm) {
      ans1 <- readline("Remove any outliers? (y/n) ")
      if (tolower(ans1) == "y"){
        flag <- "pts_rm"
        mtext("Click on an outlier", col = "red", cex = 1.2)
        out  <- locator(n = 1)
        iout <- c(iout, locatePt(out$x, out$y, std[, xvar], std[, yvar]))
        points(std[iout, xvar], std[iout, yvar], col = 2, pch = 4, cex = 2)
        if (grepl("sig", model)) {
          if (!is.null(Alow)) {
            startval <- getStart3par(std[-iout, xvar], std[-iout, yvar], Alow,
                                     ifix = NULL)
          } else {
            startval <- getStart4par(std[-iout, xvar], std[-iout, yvar],
                                     ifix = NULL, info = info)
          }
        } else {}  # model 5 (6 parameters)
        if (is.na(startval[1])) {
          flag <- "no_fit"
          warning("Unable to fit the model")
          fit <- NULL
        } else {  # rerun makeLsopt with points removed
          if (grepl("sig", model)) {
            lsopt <- makeLsopt(std[-iout, xvar], std[-iout, yvar], Alow, asym)
          } else {}
          if (asym) startval <- c(startval, a = 1)
          fit <- optim(startval, lsopt, method = optmethod,
                       control = list(maxit = maxit))
        }
        fitpar <- addParam(fit$par, Alow, asym)
        #*** INSERTED - is it a useful functionality?
        smpflag <- rep("", length(vsmp))
        mindet <- max(min(std1[, yvar]), fitpar["Alow"])
        maxdet <- min(max(std1[, yvar]), fitpar["Aup"])
        smpflag[!(mindet <= vsmp & vsmp <= maxdet)] <- "min"
        plotFit(std, xvar, yvar, dilvar, fitpar = fitpar, FUNmod = FUNmod,
                bg = bg, vsmp = vsmp, smpflag = smpflag,
                stdcol = stdcol, rugcol = rugcol, ...)
        #*** end INSERTED, uncomment plotFit() below
#        plotFit(std, xvar, yvar, dilvar, fitpar = fitpar, FUNmod = FUNmod,
#                iout = iout, bg = bg, vsmp = vsmp,
#                stdcol = stdcol, rugcol = rugcol,...)
      } else {         # answer no
        if (i == 1) {  # answer no for the first time
          revise <- FALSE
        }
        break          # answer no for 2, 3, ... time: option to revise
      }
    }
    if (revise) {
      ans2 <- readline("Revise fit? (y/n) ")
      if (tolower(ans2) == "n") {
        revise <- FALSE
        ans3 <- readline(paste("Update the existing", flag,
                               "flag? (type a new message or 'n' for no) "))
        if (tolower(ans3) != "n") {
          flag <- ans3
        }
      }
  }
  }
  ystd <- std[, yvar]
  if (!is.null(iout)) ystd <- ystd[-iout]
  bounds <- c(mindet = min(ystd), maxdet = max(ystd))
  if (is.na(startval[1])) {
    return(list(par = NULL, bounds = bounds, iout = iout, flag = flag))
  }
  bounds["mindet"] <- max(bounds["mindet"], fitpar["Alow"])
  bounds["maxdet"] <- min(bounds["maxdet"], fitpar["Aup"])

  #---------------------------- calculate bounds ------------------------------#
  if (grepl("sig", model)) {
    bounds <- c(bounds, calcBoundsSig(fitpar, range(std[, xvar])))
    if (length(bounds) == 2) {
      warning(paste(info,
                    "Fit is not good, try estimating lower asymptote if fixed"))
      return(list(par = NULL, bounds = bounds, iout = iout, flag = "no_fit"))
    }
  } else {}
  if (rm.after) {
    abline(h = bounds[c("lowerbound", "upperbound")], col = rugcol[2], lty = 4)
    legend("right", bty = "n", cex = 0.9, col = rugcol[2], lty = 4,
           legend = "bounds of the fit")
  }

  #---------------------- set bounds manually if desired ----------------------#
  if (set.bounds || overwrite.bounds) {
    overlower <- overupper <- FALSE
    if (is.na(bounds["lowerbound"])) {                             # lower bound
      cat("Unable to calculate lower bound \n")
      overlower <- TRUE
    }
    if (overlower || overwrite.bounds) {
      ans1 <- readline("Set lower bound manually? (y/n) ")
      if (tolower(ans1) == "y"){
        revise <- TRUE
        flag <- paste(flag, ", lb_manual", sep = "")
        if (!rm.after) {  # no active current plot
          #*** INSERTED - is it a useful functionality?
          #***** The problem: what if there IS current plot - add purple to rug???
          smpflag <- rep("", length(vsmp))
          mindet <- max(min(std1[, yvar]), fitpar["Alow"])
          maxdet <- min(max(std1[, yvar]), fitpar["Aup"])
          smpflag[!(mindet <= vsmp & vsmp <= maxdet)] <- "min"
          plotFit(std, xvar, yvar, dilvar, fitpar = fitpar, FUNmod = FUNmod,
                  bg = bg, vsmp = vsmp, smpflag = smpflag,
                  stdcol = stdcol, rugcol = rugcol, ...)
          #*** end INSERTED, uncomment plotFit() below
#          plotFit(std, xvar, yvar, dilvar, fitpar = fit$par, FUNmod = FUNmod,
#                  iout = iout, bg = bg, vsmp = vsmp,
#                  stdcol = stdcol, rugcol = rugcol, ...)
          abline(h = bounds[c("lowerbound", "upperbound")], col = rugcol[2],
                 lty = 4)
          legend("right", bty = "n", cex = 0.9, col = rugcol[2], lty = 4,
                 legend = "bounds of the fit")
        }
      }
      while (revise) {
        mtext("Indicate lower bound with a click", col = "red", cex = 1.2)
        bounds["lowerbound"] <- locator(n = 1)$y
        plotFit(std, xvar, yvar, dilvar, fitpar = fit$par, FUNmod = FUNmod,
                iout = iout, bg = bg, vsmp = vsmp,
                stdcol = stdcol, rugcol = rugcol, ...)
        abline(h = bounds[c("lowerbound", "upperbound")], col = rugcol[2],
               lty = 4)
        legend("right", bty = "n", cex = 0.9, col = rugcol[2], lty = 4,
               legend = "bounds of the fit")
        revise <- tolower(readline("Revise? (y/n) ")) == "y"
      }
    }
    if (is.na(bounds["upperbound"])) {                             # upper bound
      cat("Unable to calculate upper bound \n")
      overupper <- TRUE
    }
    if (overupper || overwrite.bounds) {
      ans1 <- readline("Set upper bound manually? (y/n) ")
      if (tolower(ans1) == "y"){
        revise <- TRUE
        flag <- paste(flag, ", ub_manual", sep = "")
        if (!(rm.after || overlower || overwrite.bounds)) {  # no active plot
          plotFit(std, xvar, yvar, dilvar, fitpar = fit$par, FUNmod = FUNmod,
                  iout = iout, bg = bg, vsmp = vsmp,
                  stdcol = stdcol, rugcol = rugcol, ...)
          abline(h = bounds[c("lowerbound", "upperbound")], col = rugcol[2],
                 lty = 4)
          legend("right", bty = "n", cex = 0.9, col = rugcol[2], lty = 4,
                 legend = "bounds of the fit")
        }
      }
      while (revise) {
        mtext("Indicate upper bound with a click", col = "red", cex = 1.2)
        bounds["upperbound"] <- locator(n = 1)$y
        plotFit(std, xvar, yvar, dilvar, fitpar = fit$par, FUNmod = FUNmod,
                iout = iout, bg = bg, vsmp = vsmp,
                stdcol = stdcol, rugcol = rugcol, ...)
        abline(h = bounds[c("lowerbound", "upperbound")], col = rugcol[2],
               lty = 4)
        legend("right", bty = "n", cex = 0.9, col = rugcol[3:2], lty = 4,
               legend = "bounds of the fit")
        revise <- tolower(readline("Revise? (y/n) ")) == "y"
      }
    }
  }
  return(list(par = fitpar, bounds = bounds, iout = iout, flag = flag))
}

fsig <- function(x, par) {
  par["Alow"] + (par["Aup"] - par["Alow"])/
    (1 + exp((par["Xmid"] - x)/par["Scale"]))^par["a"]
}
flin <- function(x, par) {}
fsigInv <- function(y, par) {
  -par["Scale"]*log(((par["Aup"] - par["Alow"])/
                       (y - par["Alow"]))^(1/par["a"]) - 1) + par["Xmid"]
}
flinInv <- function(y, par) {}

addParam <- function(par, Alow = NULL, asym = TRUE) {
  if ( is.null(par))  return(NULL)
  if (!is.null(Alow)) par <- c(par, Alow = Alow)
  if (!asym)          par <- c(par, a = 1)
  return(par)
}

makeLsopt <- function(x, y, Alow = NULL, asym = TRUE) {
  if (!is.null(Alow)) {
    if (!asym) {                            # Alow fixed, a = 1
      lsopt <- function(par) {
        yhat <- Alow + (par["Aup"] - Alow)/
          (1 + exp((par["Xmid"] - x)/par["Scale"]))
        sum((y - yhat)^2)
      }
    } else {
      lsopt <- function(par) {               # Alow fixed, a variable (in par)
        yhat <- Alow + (par["Aup"] - Alow)/
          (1 + exp((par["Xmid"] - x)/par["Scale"]))^par["a"]
        sum((y - yhat)^2)
      }
    }
  } else {
    if (!asym) {                            # Alow variable (in par), a = 1
      lsopt <- function(par) {
        yhat <- par["Alow"] + (par["Aup"] - par["Alow"])/
          (1 + exp((par["Xmid"] - x)/par["Scale"]))
        sum((y - yhat)^2)
      }
    } else {
      lsopt <- function(par) {               # Alow variable, a variable
        yhat <- par["Alow"] + (par["Aup"] - par["Alow"])/
          (1 + exp((par["Xmid"] - x)/par["Scale"]))^par["a"]
        sum((y - yhat)^2)
      }
    }
  }
  return(lsopt)
}

#' Calculate Bounds Numerically
#'
#' Numerically finds roots for 4th derivative of the model function that are
#' used to indicate "flat" portions of the curve.
#'
#' @details to be added
#'
#' @param FUNmod model function
#' @param par    values of function parameters
#' @param xrange range of function domain to be searched for roots
#' @param seql   length of the search grid
#'
#' @return A named list with upper and lower bounds.
#'
#' @examples
#'
#' @export

# General (might be out if not needed for non-sigmoid model)
calcBounds <- function(FUNmod, par, xrange, seql = 200) {
  x <- seq(min(xrange), max(xrange), length = seql)
  y <- FUNmod(x, par)
  deriv4 <- diff(diff(diff(diff(y))))
  ipos <- deriv4 > 0
  if (ipos[1]) {
    lowerbound <- y[which(!ipos)[1] + 2]
  } else {
    lowerbound <- NA
  }
  if (!ipos[seql - 4]) {
    upperbound <- y[rev(which(ipos))[1] + 2]
  } else {
    upperbound <- NA
  }
  c(lowerbound = lowerbound, upperbound = upperbound)
}

# this is for 4th derivative; for the 3rd, both sides positive to negative
# analytical 4th derivative and roots for sigmoid function
f4deriv <- function(x, X0, s, a) {
  w <- exp((X0 - x)/s)
  (a^3*w^3 - 6*a^2*w^2 - 4*a*w^2 + 7*a*w - w^2 + 4*w - 1)*a*w/(1 + w)^(a + 4)
}

#' Calculate Bounds for Logistic function
#'
#' Calculates roots for 4th derivative of logistic function that are used to
#' indicate "flat" portions of the curve.
#'
#' @details to be added
#'
#' @inheritParams calcBounds
#' @return A named list with upper and lower bounds.
#'
#' @examples
#'
#' @export

calcBoundsSig <- function(par, xrange) {
  wrange <- f4deriv(xrange, par["Xmid"], par["Scale"], par["a"])
  if (any(!is.finite(wrange))) return(NULL)
  lbound <- c(wrange[1] > 0, wrange[2] < 0)
  bounds <- c(lowerbound = NA, upperbound = NA)
  if (all(!lbound)) return(bounds)                # bounds beyond range(x)
  # Cardano's formula (solving cubic equation for 4th derivative)
  a <- par["a"]^3
  b <- -(6*par["a"]^2 + 4*par["a"] + 1)
  c <- 7*par["a"] + 4
  d <- -1
  Q <- (3*a*c - b^2)/(9*a^2)
  R <- (9*a*b*c - 27*a^2*d - 2*b^3)/(54*a^3)
  theta <- acos(R/sqrt(-Q^3))
  wroot <- 2*sqrt(-Q)*cos(theta/3 + c(0, 2, 4)*pi/3) - b/(3*a)
  wroot <- sort(wroot)[c(3, 1)]
  yw <- par["Alow"] + (par["Aup"] - par["Alow"])/(1 + wroot)^par["a"]
  bounds[lbound] <- yw[lbound]
  return(bounds)
}
