#' Import Array Data
#'
#' @description This function builds on the limma function read.maimages and
#' imports raw data from a time and concentration dependent toxicogenomic
#' experiment. It removes outlier arrays based on the KS.test, and outlier test based on sum and 75% quantile. The function
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
#' @param qc_coeff A named character vector giving the coefficient for acceptable intervals
#' @param qc_sum Number of tests leading to array exclusion
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
             removeOutliers = T,
             qc_coeff = c(ks = 3, sum = 3, iqr = 3, q = 3, d = 1),
             qc_sum = 1) {
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
    raw <- limma::read.maimages(
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
      ),
      verbose = F
    )
    ## not covered: SpotExtentX  gBGMeanSignal


    # output target table -----------------------------------------------------
    if (output) {
      targets
    }

    # outlier detection --------------------------------------------------------
    qc_array <- toxprofileR::arrayqc(log2(raw$E), qc_coeff = qc_coeff, qc_sum = qc_sum)

    if (length(qc_array$exclude) > 0) {
      message(paste("detected", targets$names[qc_array$exclude], "as outlier\n"))
    }

    group <- rep("include", nrow(targets))
    group[qc_array$exclude] <- "exclude"

    # plots ---------------------------------------------------------------
    if (output) {
      # density plot
      limma::plotDensities((log2(raw$E)),
        legend = T,
        group = group,
        col = c(2, "grey"),
        main = targets$substance[1]
      )

      # boxplot
      par(mar = c(12, 3, 2, 1))
      boxplot(
        log2(raw$E),
        col = (as.numeric(as.factor(group)) + 2),
        main = targets$substance[1],
        las = 2
      )
      par(mar = c(5.1, 4.1, 4.1, 1))

      # MDS plot
      limma::plotMDS(
        log2(raw$E),
        labels = raw$targets$names,
        col = as.numeric(as.factor(group)) + 2
      )
    }

    if (removeOutliers) {
      raw <-
        limma::subsetListOfArrays(
          raw,
          i = c(1:dim(raw)[1]),
          j = qc_array$include,
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
    }


    return(raw)
  }
