% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotFit.R
\name{plotFit}
\alias{plotFit}
\title{Plot the fit for standards and the samples}
\usage{
plotFit(
  std,
  xvar,
  yvar,
  fitpar = NULL,
  FUNmod = NULL,
  FUNinv = NULL,
  iout = NULL,
  bg = NULL,
  vsmp = NULL,
  smpflag = NULL,
  trimval = NULL,
  trimtype = NULL,
  extrapolate.low = FALSE,
  extrapolate.up = FALSE,
  ylim = NULL,
  tcklab = NULL,
  stdcol = c("firebrick3", "darkslategray"),
  rugcol = c("cadetblue", "purple", "firebrick2"),
  ...
)
}
\arguments{
\item{std}{matrix or data frame with standards for fitting.}

\item{xvar, yvar}{character strings for the variables used to fit a standard
curve. If \code{NULL}, first two columns are assumed to be \code{x} and
\code{y} variables.}

\item{fitpar}{values of function parameters.}

\item{FUNmod}{model function.}

\item{FUNinv}{model inverse function.}

\item{iout}{indices of removed standard points.}

\item{bg}{background values.}

\item{vsmp}{sample values.}

\item{smpflag}{character vector, flags for each sample.}

\item{trimval}{for final results, the values at which the samples are
trimmed.}

\item{trimtype}{integer vector of length two indicating if the values are
trimmed at the extremum (lower and upper).}

\item{extrapolate.low}{if \code{TRUE}, sample values beyond lower bounds will
be processed by extrapolation of the standard curve (not recommended).
Takes precedence over \code{trim.flat} value.}

\item{extrapolate.up}{if \code{TRUE}, sample values beyond upper bounds will
be processed by extrapolation of the standard curve (not recommended).
Takes precedence over \code{trim.flat} value.}

\item{ylim}{limits of the y-axis.}

\item{tcklab}{tick labels for x-axis.}

\item{stdcol}{vector of two colors for the standard points and the fit on the
plot.}

\item{rugcol}{vector of three colors for the rugplot, which indicates sample
values (inside the bounds, between the bounds and extrema, and beyond
extrema).}

\item{...}{further graphical parameters.}
}
\description{
Produces a plot that includes points for standards, proposed fit, removed
outliers, bounds for "flat" portions of the curve, and values for samples and
for the background.
}
\details{
to be added
}
