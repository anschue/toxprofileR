#' Import Array Data
#'
#' @description This function builds on the limma function read.maimages and
#' imports raw data from a time and concentration dependent toxicogenomic
#' experiment. It removes outlier arrays based on the IQR distribution and
#' returns an EListRaw.
#'
#' @param targetfile string, name of targetfile .csv, should contain columns
#' named FileName, sample_ID,	substance, concentration_level,	concentration_umol_l,
#' time_hpe,	type, comment, and	scan_ID
#'
#' @param datadir string, path to raw data
#'
#' @param scanID numerical, ID of scan (default: 1)
#'
#' @param output logical, should function return target (default: T)
#'
#' @param removeOutliers logical, should outliers be removed
#'
#'
#' @return EListraw (limma class)
#'
#' @export
#'
#'
importData <-
  function(targetfile,
             datadir,
             scanID = 1,
             output = T,
             removeOutliers = T) {
    # read target file --------------------------------------------------------
    targets <- limma::readTargets(file = targetfile, sep = "\t")
    targets <- targets[targets$scan_ID == scanID, ]

    # create unique names based on concentration and time ---------------------
    targets$names <-
      make.names(
        paste(
          targets$substance,
          targets$concentration_level,
          targets$time_hpe,
          targets$type,
          sep = "_"
        ),
        unique = T
      )

    # read rawdata ------------------------------------------------------------
    rawdata <- limma::read.maimages(
      files = targets,
      source = "agilent",
      path = datadir,
      names = targets$names,
      green.only = T,
      columns = list(
        E = "gProcessedSignal",
        Processederror = "gProcessedSigError",
        Median = "gMedianSignal",
        Medianb = "gBGMedianSignal",
        processedSignal = "gProcessedSignal",
        isNonUniform = "gIsFeatNonUnifOL",
        isNonUniformBG = "gIsBGNonUnifOL",
        isPopOutlier = "gIsFeatPopnOL",
        isPopOutlierBG = "gIsBGPopnOL",
        manualFlag = "IsManualFlag",
        posandsigf = "gIsPosAndSignif",
        aboveBG = "gIsWellAboveBG",
        bgSubSignal = "gBGSubSignal"
      )
    )
    ## not covered: SpotExtentX  gBGMeanSignal

    raw <- rawdata

    # output target table -----------------------------------------------------
    if (output) {
      knitr::kable(targets)
    }

    # remove outlier array(s) -------------------------------------------------
    iqrs <- apply(log2(raw$E), MARGIN = 2, IQR)
    int <-
      findInterval(iqrs, vec = c(mean(iqrs) - (2 * sd(iqrs)), mean(iqrs) +
        (2 * sd(iqrs))))

    exclude <- c(which(int != 1))
    include <- c(1:length(iqrs))
    include <- include[!include %in% exclude]

    group <- rep("include", length(iqrs))
    group[exclude] <- "exclude"

    if(output){
    limma::plotDensities((log2(raw$E)),
      legend = T,
      group = group,
      col = c(2, "grey")
    )}

    rawdata <-
      limma::subsetListOfArrays(
        raw,
        i = c(1:dim(raw)[1]),
        j = include,
        IJ = c(
          "E",
          "Processederror",
          "Median",
          "Medianb",
          "isNonUniform",
          "isNonUniformBG",
          "isPopOutlier",
          "isPopOutlierBG",
          "manualFlag",
          "posandsigf",
          "aboveBG",
          "bgSubSignal"
        ),
        IX = "genes",
        JX = "targets",
        I = NULL
      )
    targets <- targets[include, ]

    return(rawdata)
  }
