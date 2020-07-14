
#' Reformat raw luminex files for processing
#'
#' Read in raw luminex files, extract standards and background (if provided), and format samples for processing.
#'
#' @details to be added
#'
#' @param antigen character string.
#' @param fname name of the file that contains raw data.
#' @param fdir directory where the file is located (alternatively, full path can
#'   be included in \code{fname}).
#' @param dtype character string for data type in the file.
#' @param stdstr character string indicating standards in the file's "Sample"
#'   column. Not case sensitive. If \code{""} (empty string), standards will be
#'   determined by the pattern "1/" only.
#' @param bgstr character string indicating background in the file's "Sample"
#'   column. Not case sensitive.
#' @param stddil a vector of standard dilutions. If \code{NULL}, dilutions are
#'   inferred from the file's "Sample" column. Non-null value can be used to
#'   exclude some dilutions from model fitting.
#' @param nwells number of wells. If \code{NULL}, inferred from the file.
#' @param nsep number of lines separating different data types in the file.
#' @param ncolmax maximum number of columns in the file.
#' @param dformat date format in the file.
#' @inheritParams processSmp
#'
#' @return A list with standards for fitting, background values, sample values,
#'   indices for samples.
#'
#' @export

prepLum <- function(antigen, fname, fdir = NULL, dtype = "Median",
                    stdstr = "std|stand", bgstr  = "blank|background",
                    stddil = NULL, smpdil = 1000, nwells = NULL, nsep = 2,
                    ncolmax = 105, dformat = "%m/%d/%Y") {
  fpath <- paste(fdir, fname, sep = "")
  MFI <- read_data(fpath, dtype, nwells, nsep, ncolmax)
  if (length(grep("1/", MFI$Sample)) == 0) {
    stop('Please provide standard concentrations in "1/dilution" format,
         e.g "1/200"')
  }
  fdate <- read.csv(fpath, header = FALSE, skip = 2, nrow = 1)[1, 2]
  pdate <- as.Date(as.character(fdate), format = dformat)
  if (is.na(pdate)) pdate <- fdate  #*** optionally - issue a warning?
  if (nchar(stdstr) > 0 && length(grep(stdstr, MFI$Sample)) == 0) {
    warning(paste("File ", fname, ' does not contain "', stdstr, '"; using
                  "1/" to determine standards', sep = ""))
    stdstr <- ""
  }
  return(extractStd(MFI, stdstr, bgstr, stddil, smpdil, antigen, fname, pdate))
}
  #*** figure out: how to keep both columns in std (or not)
  #    and mainly what to keep in smp (since it's also dfout)
  #    do not need logConc for sure, but maybe logMFI???
  #    in plot, add "log" to axis labels (?)
  #*** do NOT transform smp df; either make a separate vector, logged - or
  #    do every time in plotFit() and once in normalizeSmp()
  #*** then maybe better - to exclude non-smp in Luminex

# dilvar = Dilution, smpvar = logMFI (?), resvar = Conc

#' Extract Standards
#'
#' Separates specified serial dilutions for fitting, background, and samples
#'
#' details
#'
#' @param MFI data frame with variable "Sample" that contains information about
#'   standard concentrations in "1/dilution" format (e.g. "1/200").
#' @param stddil standard dilutions to use for standard curve fitting. If
#'   \code{NULL}, all the dilutions indicated in a "Sample" variable are used.
#' @param pdate date of the plate processing.
#' @inheritParams prepLum
#'
#' @return A list with standards for fitting, background values, sample values,
#'   indices for samples.
#'
#' @export

extractStd <- function(MFI, stdstr, bgstr, stddil, smpdil, antigen, fname,
                       pdate) {
  mfi <- suppressWarnings(as.numeric(MFI[, grep(antigen, colnames(MFI))]))
  istd <- grep("1/", MFI$Sample)
  if (length(istd) == 0) {    # handled by processLum(); for standalone use only
    stop('Please provide standard concentrations in "1/dilution" format,
         e.g "1/200"')
  }
  if (stdstr != "") {
    istd <- istd[grep(tolower(stdstr), tolower(MFI$Sample)[istd])]
    if (length(istd) == 0) {
      istd <- grep("1/", MFI$Sample)  # no warning, back to original istd
    }
  }
  sname <- MFI$Sample[istd]
  locdil <- regexpr("1/", sname)
  sdilut <- as.numeric(substr(sname, locdil + 2, nchar(sname)))
  if (!is.null(stddil)) {
    idil <- sdilut %in% stddil
  } else {
    idil <- seq_along(sdilut)
  }
  std <- data.frame(Dilution = sdilut[idil], Conc = 1/sdilut[idil],
                    MFI = mfi[istd[idil]])
  ibg  <- grep(tolower(bgstr),  tolower(MFI$Sample))
  #  ineg <- grep(tolower(negstr), tolower(MFI$Sample))
  ismp <- (1:nrow(MFI))[-c(istd, ibg)]
  bg <- mfi[ibg]
  std <- std[order(std$Conc), ]  # sorted by increasing concentration

  dfout                <- MFI[, 1:2]
  names(dfout)[2]      <- "ID"  # was "Sample"; 1st column is "Location"
  dfout$antigen        <- antigen
  dfout$file_name      <- fname
  dfout$date_plate     <- pdate
  dfout$MFI            <- mfi     # alternatively, logMFI, then fitlog = "x"
  dfout$Dilution       <- smpdil
  dfout$Dilution[ibg]  <- 0
  dfout$Dilution[istd] <- sdilut
  return(list(smp = dfout, std = std, bg = bg, ismp = ismp))
}

#' Read in raw luminex data
#'
#' Extract specified data type from the file.
#'
#' @details to be added
#'
#' @inherit prepLum params
#'
#' @return A data frame.
#'
#' @export

read_data <- function(fname, dtype = "Median", nwells = NULL, nsep = 2,
                      ncolmax = 105) {
  all <- read.csv(fname, header = FALSE, blank.lines.skip = FALSE,
                  as.is = TRUE, col.names = 1:ncolmax)
  itype <- grep("Data", all[[1]])
  types <- all[itype, 2]
  if (is.null(nwells)) {
    nwells <- itype[2] - itype[1] - nsep - 1
  }
  nskip <- itype[tolower(types) %in% tolower(dtype)]
  if (length(nskip) == 0) {                       # what if > 1 (can't happen?)
    stop("Specified data type not found")
  }
  read.csv(fname, as.is = TRUE, skip = nskip, nrows = nwells)
}
