#' Read Microarray data from GEO or ArrayExpress
#'
#' @param datadir data directory
#' @param rawformat format of rawdatat (Agilen or Affymetrix)
#' @param betweenArrayNorm name of normalization method
#' @param metadata metadataframe
#'
#' @return returns a normalized Elist
#'
#' @import limma
#' @export
#'
read_raw_public <- function(datadir,
                            rawformat = c("Agilent", "Affymetrix", "Affymetrix_ST"),
                            betweenArrayNorm = "cyclicloess",
                            metadata) {

    # Agilent ---------------------------------------------------------------------------
    if (rawformat == "Agilent") {

        # loading raw data --------------------------------------------------------------
        metadata$FileName <- metadata$Array.Data.File

        raw <- limma::read.maimages(
            metadata,
            names = metadata$gsm.gsm,
            source = "agilent",
            path = datadir,
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
            sep = "\t",
            quote = ""
        )

        # remove outlier array(s) -------------------------------------------------
        include <- c(1:nrow(metadata))

        outliermetrics1 <- arrayQualityMetrics::outliers(log2(raw$E), method = "KS")
        outliermetrics2 <- arrayQualityMetrics::outliers(log2(raw$E), method = "sum")
        outliermetrics3 <- arrayQualityMetrics::outliers(log2(raw$E), method = "upperquartile")

        exclude <- table(c(as.numeric(outliermetrics1@which),
                           as.numeric(outliermetrics2@which),
                           as.numeric(outliermetrics3@which)))

        exclude <- as.numeric(names(exclude)[exclude>=2])

        include <- include[!include %in% exclude]

        message(paste("detected", metadata$gsm.gsm[exclude], "as outlier\n"))

        group <- rep("include", length(iqrs))
        group[exclude] <- "exclude"

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


        # normalization and log-transformation ---------------------------------
        message("Normalizing")
        data.log.norm <-
            limma::normalizeBetweenArrays(rawdata, method = betweenArrayNorm)

        # getting median of duplicate Probes on Array --------------------------
        message("Summarizing duplicate probes")
        data.log.norm.median <-
            toxprofileR::elist.median(data.log.norm)


        # annotation -----------------------------------------------------------
        message("Updating Annotation")
        data.log.norm.median$genes$ProbeName <-
            make.names(data.log.norm.median$genes$ProbeName)
        rownames(data.log.norm.median$E) <-
            make.names(rownames(data.log.norm.median$E))

        annotation_new <-
            readRDS(
                paste0(
                    "./data/PlatformData/final_annotation/",
                    metadata$gpl_id[1],
                    "annotation.Rds"
                )
            )

        colnames(data.log.norm.median$genes) <-
            paste(colnames(data.log.norm.median$genes), "old", sep = "_")

        data.log.norm.median$genes <-
            merge(
                data.log.norm.median$genes,
                annotation_new,
                by.x = "ProbeName_old",
                by.y = "ProbeID",
                all.x = T,
                sort = F
            )

        data.log.norm.median$genes <-
            data.log.norm.median$genes[match(
                rownames(data.log.norm.median$E),
                data.log.norm.median$genes$ProbeName_old
            ), ]
        colnames(data.log.norm.median$genes)[colnames(data.log.norm.median$genes) ==
                                                 "ProbeName_old"] <- "ProbeName"

        return(data.log.norm.median)
    }

    # Affymetrix ------------------------------------------------------------------------
    if (rawformat == "Affymetrix_ST" | rawformat == "Affymetrix") {

        rownames(metadata) <- as.character(metadata$gsm.gsm)

        message("read data")
        data.raw <-
            oligo::read.celfiles(
                filenames = paste0(datadir, metadata$Array.Data.File),
                sampleNames = metadata$gsm.gsm
            )

        pData(data.raw) <- metadata

        # remove outlier array(s) -------------------------------------------------
        include <- c(1:nrow(metadata))

        outliermetrics1 <- arrayQualityMetrics::outliers(log2(raw$E), method = "KS")
        outliermetrics2 <- arrayQualityMetrics::outliers(log2(raw$E), method = "sum")
        outliermetrics3 <- arrayQualityMetrics::outliers(log2(raw$E), method = "upperquartile")

        exclude <- table(c(as.numeric(outliermetrics1@which),
                           as.numeric(outliermetrics2@which),
                           as.numeric(outliermetrics3@which)))

        exclude <- as.numeric(names(exclude)[exclude>=2])

        include <- include[!include %in% exclude]

        message(paste("detected", metadata$gsm.gsm[exclude], "as outlier\n"))

        data.raw1 <- data.raw[, include]
        metadata <- metadata[include, ]

        # normalization ---------------------------------------------------
        message("normalization (rma)")
        data.norm <- oligo::rma(data.raw1)

        message("creating EList")
        data.norm.E <-
            new("EList", list(
                E = exprs(data.norm),
                targets = pData(data.norm),
                genes = data.frame(ProbeID = rownames(exprs(data.norm)))
            ))


        ####update Annotation
        message("update Annotation")
        annotation_new <-
            readRDS(
                paste0(
                    "./data/PlatformData/final_annotation/",
                    metadata$gpl_id[1],
                    "annotation.Rds"
                )
            )



        data.norm.E$genes <-
            merge(
                data.norm.E$genes,
                annotation_new,
                by = "ProbeID",
                all.x = T,
                sort = F
            )

        rownames(data.norm.E$E) <-
            make.names(rownames(data.norm.E$E))
        data.norm.E$genes$ProbeID <-
            make.names(data.norm.E$genes$ProbeID)
        data.norm.E$genes <-
            data.norm.E$genes[match(rownames(data.norm.E$E), data.norm.E$genes$ProbeID), ]

        colnames(data.norm.E$genes)[colnames(data.norm.E$genes) == "ProbeID"] <-
            "ProbeName"

        return(data.norm.E)
    }

}
