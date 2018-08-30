#' Array Quality Control
#'
#' @param exprs An expression data.frame
#' @param qc_coeff A named character vector giving the coefficient for acceptable intervals
#' @param qc_sum Number of tests leading to array exclusion
#'
#' @return A list containing indices for columns to include and exclude
#' @export
#'
arrayqc <- function(exprs, qc_coeff = c(ks = 3, sum = 3, iqr = 3, q = 3, d = 1), qc_sum = 1){
    include <- c(1:ncol(exprs))

    # KS statistics (taken from arrayQualityMetrics package)
    fx <- ecdf(as.vector(exprs))
    kss <- suppressWarnings(apply(exprs, 2, function(v){ks.test(v, y = fx, alternative = "two.sided")$statistic}))

    kss_ex <- which(findInterval(kss,vec = c(quantile(kss, 0.25, na.rm = T) - qc_coeff["ks"] * IQR(kss, na.rm = T),
                                              quantile(kss, 0.75, na.rm = T) + qc_coeff["ks"] * IQR(kss, na.rm = T)))>1)


    # Expression sums
    sums <- colSums(exprs, na.rm = T)

    sums_ex <- which(findInterval(sums,vec = c(quantile(sums, 0.25, na.rm = T) - qc_coeff["sum"] * IQR(sums, na.rm = T),
                                               quantile(sums, 0.75, na.rm = T) + qc_coeff["sum"] * IQR(sums, na.rm = T)))!=1)

    # IQR

    iqrs <- apply(exprs, 2, IQR, na.rm = T)

    iqrs_ex <- which(findInterval(iqrs,vec = c(quantile(iqrs, 0.25, na.rm = T) - qc_coeff["iqr"] * IQR(iqrs, na.rm = T),
                                               quantile(iqrs, 0.75, na.rm = T) + qc_coeff["iqr"] * IQR(iqrs, na.rm = T)))!=1)


    # Quantiles
    qs <-
        apply(exprs,
              2,
              quantile,
              na.rm = TRUE,
              probs = c(0.25, 0.75))

    qs_ex <- which(findInterval(qs[1,], vec = c((quantile(qs[1, ], 0.25, na.rm = T) - qc_coeff["q"] * IQR(qs[1, ], na.rm = T)),
                                                (quantile(qs[1, ], 0.75, na.rm = T) + qc_coeff["q"] * IQR(qs[1, ], na.rm = T))))!=1|
                       findInterval(qs[2,], vec = c((quantile(qs[2, ], 0.25, na.rm = T) - qc_coeff["q"] * IQR(qs[2, ], na.rm = T)),
                                                    (quantile(qs[2, ], 0.75, na.rm = T) + qc_coeff["q"] * IQR(qs[2, ], na.rm = T))))!=1)

    # Euclidean distance
    MDS_out <-
        limma::plotMDS(exprs, top = ceiling(nrow(exprs) / 5), plot = F)
    distances <- apply(MDS_out$distance.matrix, 2, mean)
    dist_ex <- which(findInterval(distances,vec = c(quantile(distances, 0.25, na.rm = T) - qc_coeff["d"] * IQR(distances, na.rm = T),
                                                    quantile(distances, 0.75, na.rm = T) + qc_coeff["d"] * IQR(distances, na.rm = T)))!=1)



    exclude_tab <- table(c(
        kss_ex,
        sums_ex,
        iqrs_ex,
        qs_ex,
        dist_ex
    ))

    # exclude array if 2 of the tests detect it as outlier
    exclude <- as.numeric(names(exclude_tab)[exclude_tab >= qc_sum])
    include <- include[!include %in% exclude]

    out <- list(include = include, exclude = exclude)

    out
}
