#' Array Quality Control
#'
#' @param exprs An expression data.frame
#'
#' @return A list containing indices for columns to include and exclude
#' @export
#'
arrayqc <- function(exprs){
    include <- c(1:ncol(exprs))

    # KS statistics
    ks_ex <-
        as.numeric(arrayQualityMetrics::outliers(exprs, method = "KS")@which)

    # Expression sums
    sums <- colSums(exprs, na.rm = T)

    sums_ex <- which(findInterval(sums,vec = c(quantile(sums, 0.25, na.rm = T) - 3 * IQR(sums, na.rm = T),
                                               quantile(sums, 0.75, na.rm = T) + 3 * IQR(sums, na.rm = T)))!=1)

    # Quantiles
    qs <-
        apply(exprs,
              2,
              quantile,
              na.rm = TRUE,
              probs = c(0.25, 0.75))

    qs_ex <- which(findInterval(qs[1,], vec = c((quantile(qs[1, ], 0.25, na.rm = T) - 3 * IQR(qs[1, ], na.rm = T)),
                                                (quantile(qs[1, ], 0.75, na.rm = T) + 3 * IQR(qs[1, ], na.rm = T))))!=1|
                       findInterval(qs[2,], vec = c((quantile(qs[2, ], 0.25, na.rm = T) - 3 * IQR(qs[2, ], na.rm = T)),
                                                    (quantile(qs[2, ], 0.75, na.rm = T) + 3 * IQR(qs[2, ], na.rm = T))))!=1)

    # Euclidean distance
    MDS_out <-
        limma::plotMDS(exprs, top = ceiling(nrow(exprs) / 5), plot = F)
    distances <- apply(MDS_out$distance.matrix, 2, mean)
    dist_ex <- which(findInterval(distances,vec = c(quantile(distances, 0.25, na.rm = T) - 1 * IQR(distances, na.rm = T),
                                                    quantile(distances, 0.75, na.rm = T) + 1 * IQR(distances, na.rm = T)))!=1)



    exclude_tab <- table(c(
        ks_ex,
        sums_ex,
        qs_ex,
        dist_ex
    ))

    # exclude array if 2 of the tests detect it as outlier
    exclude <- as.numeric(names(exclude_tab)[exclude_tab >= 2])
    include <- include[!include %in% exclude]

    out <- list(include = include, exclude = exclude)

    out
}
