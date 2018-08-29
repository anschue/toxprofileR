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
        targets <- targets[targets$scan_ID == scanID,]

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
        include <- c(1:nrow(targets))

        # KS statistics
        outliermetrics1 <-
            arrayQualityMetrics::outliers(log2(raw$E), method = "KS")

        # Expression sums
        sums <- colSums(log2(raw$E), na.rm = T)
        th_up <-
            (quantile(sums, 0.75, na.rm = T) + 1.5 * IQR(sums, na.rm = T))
        th_down <-
            (quantile(sums, 0.25, na.rm = T) - 1.5 * IQR(sums, na.rm = T))
        outliermetrics2 <-
            list(
                threshold = c(th_up, th_down),
                which = which(sums > th_up | sums < th_down)
            )

        # Quantiles
        qs <-
            apply(log2(raw$E),
                  2,
                  quantile,
                  na.rm = TRUE,
                  probs = c(0.25, 0.75))
        th_up <-
            (quantile(qs[2, ], 0.75, na.rm = T) + 1.5 * IQR(qs[2, ], na.rm = T))
        th_down <-
            (quantile(qs[1, ], 0.25, na.rm = T) - 1.5 * IQR(qs[1, ], na.rm = T))
        outliermetrics3 <-
            list(
                threshold = c(th_up, th_down),
                which = which(qs[2, ] > th_up | qs[1, ] < th_down)
            )

        # Euclidean distance
        MDS_out <-
            limma::plotMDS(log2(raw$E), top = (nrow(raw$E) / 5), plot = F)
        distances <- apply(MDS_out$distance.matrix, 2, mean)
        outliermetrics4 <-
            arrayQualityMetrics::boxplotOutliers(distances, 1)


        exclude_tab <- table(c(
            as.numeric(outliermetrics1@which),
            as.numeric(outliermetrics2$which),
            as.numeric(outliermetrics3$which),
            as.numeric(outliermetrics4$which)
        ))

        # exclude array if 2 of the tests detect it as outlier
        exclude <- as.numeric(names(exclude_tab)[exclude_tab >= 2])
        include <- include[!include %in% exclude]

        if(length(exclude)>0){
        message(paste("detected", targets$names[exclude], "as outlier\n"))
        }

        group <- rep("include", nrow(targets))
        group[exclude] <- "exclude"

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
