#' Remove outliers
#'
#' @param nodeframe A data.frame containing a "logFC" column
#' @param plots logical, should plots with outliers be shown?
#'
#' @return The input data.frame with outlying measurements removed
#' @export
#'
remove_outliers <- function(nodeframe, plots = F) {
  library("outliers")

  #### detect and remove outliers (according to grubbs.test)###
  outliernew <- as.numeric(outliers::grubbs.test(x = nodeframe$logFC)["p.value"])

  if (plots & log10(outliernew) < -3) {
    plot(nodeframe$logFC, main = nodeframe$nodeID)
  }
  while (log10(outliernew) < -3) {
    if (plots) {
      points(x = which(nodeframe$logFC == outliers::outlier(nodeframe$logFC)), y = outliers::outlier(nodeframe$logFC), col = "red")
    }
    nodeframe <- nodeframe[-which(nodeframe$logFC == outliers::outlier(nodeframe$logFC)), ]
    outliernew <- as.numeric(outliers::grubbs.test(x = nodeframe$logFC)["p.value"])
  }

  nodeframe
}
