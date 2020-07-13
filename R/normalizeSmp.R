#' Normalize Samples
#'
#' Calculate sample concentration based on the fit of standard dilutions.
#'
#' @details
#'
#' @param smp     data frame containing samples
#' @param FUNinv  inverse function to infer sample concentration
#' @param par     values of model function parameters
#' @param bounds  named vector with values for extrema and bounds indicating
#'   "flat" regions of the curve
#' @param fitflag flag for the fit as returned by \code{fitStd}
#' @inheritParams processSmp
#'
#' @return A data frame.
#'
#' @export

normalizeSmp <- function(smp, smpvar, resvar, dilvar, FUNinv, par, bounds,
                         fitflag, fitlog, trim.flat = TRUE) {
  if (is.null(resvar)) {
    resvar <- "conc"
  }
  smp[, resvar]   <- NA
  smp$trimmed     <- FALSE                 # for convenience (redundant)
  smp$flagged_run <- fitflag
  smp$min         <- bounds["mindet"]
  smp$max         <- bounds["maxdet"]
  smp$lower_bound <- bounds["lowerbound"]  # if NA, good - no flat portion
  smp$upper_bound <- bounds["upperbound"]
  smp$Flag        <- ""
  if (is.null(par)) {
    smp$trim_lo <- smp$trim_up <- NA
    return(smp)
  }

  imin  <- smp[, smpvar] < bounds["mindet"]
  imax  <- smp[, smpvar] > bounds["maxdet"]
  ilobd <- smp[, smpvar] < bounds["lowerbound"]
  iupbd <- smp[, smpvar] > bounds["upperbound"]
  smp$Flag[ilobd] <- "below_lower_bound"
  smp$Flag[iupbd] <- "above_upper_bound"           # NA's OK
  smp$Flag[imin]  <- "below_min"                   # overwrites
  smp$Flag[imax]  <- "above_max"                   # overwrites

  if (trim.flat && !is.na(bounds["lowerbound"])) {
    ilo <- ilobd
    smp$trim_lo <- bounds["lowerbound"]
  } else {
    ilo <- imin
    smp$trim_lo <- bounds["mindet"]
  }
  if (trim.flat && !is.na(bounds["upperbound"])) {
    iup <- iupbd
    smp$trim_up <- bounds["upperbound"]
  } else {
    iup <- imax
    smp$trim_up <- bounds["maxdet"]
  }
  imid <- which(!(ilo | iup))  # get rid of NA's
  ilo  <- which(ilo)
  iup  <- which(iup)

  # if xvar is log-transformed in std, exponentiate output of FUNinv()
  fx <- ifelse(grepl("x", fitlog), function(x) exp(x), function(x) x)
  smp[ilo,  resvar] <- fx(FUNinv(smp$trim_lo[1],    par))*smp[ilo,  dilvar]
  smp[iup,  resvar] <- fx(FUNinv(smp$trim_up[1],    par))*smp[iup,  dilvar]
  smp[imid, resvar] <- fx(FUNinv(smp[imid, smpvar], par))*smp[imid, dilvar]
  smp$trimmed[c(ilo, iup)] <- TRUE
  return(smp)
}
