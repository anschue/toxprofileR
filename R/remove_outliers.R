#' Remove outliers
#'
#' @param nodedf A data.frame containing a "logFC" column
#' @param plots logical, should plots with outliers be shown?
#'
#' @return The input data.frame with outlying measurements removed
#' @export
#'
remove_outliers <- function(nodedf, plots = F) {
  library("outliers")

  #### detect and remove outliers (according to grubbs.test)###
  outliernew <- as.numeric(outliers::grubbs.test(x = nodedf$logFC)["p.value"])

  if (plots & log10(outliernew) < -3) {
    plot(nodedf$logFC, main = nodedf$nodeID)
  }
  while (log10(outliernew) < -3) {
    if (plots) {
      points(x = which(nodedf$logFC == outliers::outlier(nodedf$logFC)), y = outliers::outlier(nodedf$logFC), col = "red")
    }
    nodedf <- nodedf[-which(nodedf$logFC == outliers::outlier(nodedf$logFC)), ]
    outliernew <- as.numeric(outliers::grubbs.test(x = nodedf$logFC)["p.value"])
  }

  nodedf
}
