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
            print(targets)
        }

        # outlier detection --------------------------------------------------------
        include <- c(1:nrow(targets))

        outliermetrics1 <-
            arrayQualityMetrics::outliers(log2(raw$E), method = "KS")
        outliermetrics2 <-
            arrayQualityMetrics::outliers(log2(raw$E), method = "sum")
        outliermetrics3 <-
            arrayQualityMetrics::outliers(log2(raw$E), method = "upperquartile")

        exclude <- table(c(
            as.numeric(outliermetrics1@which),
            as.numeric(outliermetrics2@which),
            as.numeric(outliermetrics3@which)
        ))

        # exclude array if 2 of the tests detect it as outlier
        exclude <- as.numeric(names(exclude)[exclude >= 2])
        include <- include[!include %in% exclude]

        message(paste("detected", targets$names[exclude], "as outlier\n"))

        group <- rep("include", nrow(targets))
        group[exclude] <- "exclude"

        # plots ---------------------------------------------------------------
        if (output) {
            limma::plotDensities((log2(raw$E)),
                                 legend = T,
                                 group = group,
                                 col = c(2, "grey"),
                                 main = targets$substance[1]
            )

            boxplot((log2(raw$E)), col = (as.numeric(as.factor(group)) + 2), main = targets$substance[1])


        }

        if (removeOutliers) {
            raw <-
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
        }


        return(raw)
    }
