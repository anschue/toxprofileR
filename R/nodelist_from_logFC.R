#' Create nodelist from logFC dataframe
#'
#' @param logFCframe a dataframe containing the columns 'logFC' and 'ensembl_gene_id'
#' @param nodeframe universe nodeframe
#' @param substance name of substance
#' @param time_hpe time after exposure
#' @param concentration_umol_l concentration
#' @param concentration_level concentration level
#' @param time_hpe_factor time after exposure
#'
#' @return returns a nodelist
#' @export
nodelist_from_logFC <- function(logFCframe, nodeframe, substance = "substance", time_hpe = 1, concentration_umol_l = 1, concentration_level = "C1", time_hpe_factor = 1) {
  nodeframe$ensembl <- as.character(nodeframe$ensembl)
  logFCframe$ensembl_gene_id <- as.character(logFCframe$ensembl_gene_id)

  # create list with logFC and metadata for each node -----------------------
  nodelist <- lapply(X = seq(1, max(nodeframe$toxnode)), FUN = function(nodeID) {
    if (sum(logFCframe$ensembl_gene_id %in% nodeframe$ensembl[nodeframe$toxnode == nodeID]) > 0) {
      logFC <- logFCframe$logFC[logFCframe$ensembl_gene_id %in% nodeframe$ensembl[nodeframe$toxnode == nodeID]]
      ensembl_gene_id <- logFCframe$ensembl_gene_id[logFCframe$ensembl_gene_id %in% nodeframe$ensembl[nodeframe$toxnode == nodeID]]
      probeframe <- data.frame(logFC = logFC, concentration_umol_l = concentration_umol_l, concentration_level = concentration_level, time_hpe_factor = time_hpe_factor, time_hpe = time_hpe, probe_id = NA, ensembl_gene_id = ensembl_gene_id, nodeID = nodeID, substance = substance)
      return(probeframe)
    } else {
      return(NA)
    }
  })

  nodelist
}
