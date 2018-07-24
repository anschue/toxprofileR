#' Calculate logFC
#'
#' @param elist an EList with normalized expression values
#'
#' @return an EList with logFC values
#' @export
#'
#' @examples
calc_logfc <- function(elist) {
  elist_logfc <- elist
  times <- unique(elist$targets$time_hpe)
  elist_logfc$E <- t(pbapply::pbapply(elist$E, MARGIN = 1, FUN = function(gene) {
    gene_logfc <- gene
    for (time in times) {
      gene_logfc[elist$targets$time_hpe == time] <- gene[elist$targets$time_hpe == time] - mean(gene[elist$targets$time_hpe == time & elist$targets$concentration_umol_l == 0])
    }
    return(gene_logfc)
  }))
  return(elist_logfc)
}
