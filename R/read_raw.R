#' Read Microarray data from GEO or ArrayExpress
#'
#' @param datadir data directory
#' @param rawformat format of rawdatat (Agilen or Affymetrix)
#' @param betweenArrayNorm name of normalization method
#' @param metadata metadataframe
#' @param output logical if output plot should be produced
#' @param qc_coeff A named character vector giving the coefficient for acceptable intervals
#' @param qc_sum Number of tests leading to array exclusion
#'
#' @return returns a normalized Elist
#'
#' @import limma
#' @export
#'
read_raw_public <- function(datadir,
                            rawformat = c("Agilent", "Affymetrix", "Affymetrix_ST"),
                            betweenArrayNorm = "cyclicloess",
                            metadata,
                            output = T,
                            qc_coeff = c(ks = 3, sum = 3, iqr = 3, q = 3, d = 1),
                            qc_sum = 1) {

  # Agilent ---------------------------------------------------------------------------
  if (rawformat == "Agilent") {

    # loading raw data --------------------------------------------------------------
    metadata$FileName <- metadata$Array.Data.File

    metadata <- metadata[!duplicated(metadata$FileName), ]

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
    message("outlier detection")
    qc_array <- toxprofileR::arrayqc(exprs = log2(raw$E), qc_coeff = qc_coeff, qc_sum = qc_sum)

    if (length(qc_array$exclude) > 0) {
      message(paste("detected", metadata$gsm.gsm[qc_array$exclude], "as outlier\n"))
    }

    group <- rep("include", nrow(metadata))
    group[qc_array$exclude] <- "exclude"

    # plots ---------------------------------------------------------------
    if (output) {
      # density plot
      limma::plotDensities((log2(raw$E)),
        legend = T,
        group = group,
        col = c(2, "grey"),
        main = metadata$study_id[1]
      )

      # boxplot
      par(mar = c(12, 3, 2, 1))
      boxplot(
        log2(raw$E),
        col = (as.numeric(as.factor(group)) + 2),
        main = metadata$study_id[1],
        las = 2
      )
      par(mar = c(5.1, 4.1, 4.1, 1))

      # MDS plot
      limma::plotMDS(
        log2(raw$E),
        labels = metadata$SampleName,
        col = as.numeric(as.factor(group)) + 2
      )
    }

    rawdata <-
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
    message("outlier detection")
    qc_array <- toxprofileR::arrayqc(log2(exprs(data.raw)), qc_coeff = qc_coeff, qc_sum = qc_sum)

    if (length(qc_array$exclude) > 0) {
      message(paste("detected", metadata$gsm.gsm[qc_array$exclude], "as outlier\n"))
    }

    group <- rep("include", ncol(exprs(data.raw)))
    group[qc_array$exclude] <- "exclude"

    # plots ---------------------------------------------------------------
    if (output) {
      # density plot
      limma::plotDensities(log2(exprs(data.raw)),
        legend = T,
        group = group,
        col = c(2, "grey"),
        main = metadata$study_id[1]
      )

      # boxplot
      par(mar = c(12, 3, 2, 1))
      boxplot(
        log2(exprs(data.raw)),
        col = (as.numeric(as.factor(group)) + 2),
        main = metadata$study_id[1],
        las = 2
      )
      par(mar = c(5.1, 4.1, 4.1, 1))

      # MDS plot
      limma::plotMDS(
        log2(exprs(data.raw)),
        labels = colnames(exprs(data.raw)),
        col = as.numeric(as.factor(group)) + 2
      )
    }

    data.raw1 <- data.raw[, qc_array$include]
    metadata <- metadata[qc_array$include, ]

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


    #### update Annotation
    message("update Annotation")
    annotation_new <-
      readRDS(
        paste0(
          "./data/PlatformData/final_annotation/",
          metadata$gpl_id[1],
          "annotation.Rds"
        )
      )

    data.norm.E$genes$ProbeID <- make.names(data.norm.E$genes$ProbeID)
    annotation_new$ProbeID <- make.names(annotation_new$ProbeID)

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

    data.norm.E$genes <-
      data.norm.E$genes[match(rownames(data.norm.E$E), data.norm.E$genes$ProbeID), ]

    colnames(data.norm.E$genes)[colnames(data.norm.E$genes) == "ProbeID"] <-
      "ProbeName"

    return(data.norm.E)
  }
}
