#' Create logFC Nodelist
#'
#' @param elist An limma EList object with normalized logFC values
#'
#' @return A List with logFC values for each node of the toxicogenomic universe
#' @export
create_nodelist <- function(elist) {
  nodeframe$ensembl <- as.character(nodeframe$ensembl)
  nodeframe$ProbeID <- as.character(nodeframe$ProbeID)
  elist$genes$ProbeName <- as.character(elist$genes$ProbeName)

  nodeframe_elist <- nodeframe[!is.na(nodeframe$ProbeID), ]

  if (is.null(elist$targets$type)) {
    elist$targets$type <- "treatment"
  } # for smoothed datasets

  # create list with logFC and metadata for each node -----------------------
  nodelist <- lapply(X = seq(1, max(nodeframe_elist$toxnode)), FUN = function(nodeID) {
    if (sum(nodeframe_elist$toxnode == nodeID) > 0) {
      logFC <- c(t(elist$E[elist$genes$ProbeName %in% nodeframe_elist$ProbeID[nodeframe_elist$toxnode == nodeID], elist$targets$type != "recovery"]))
      concentration_umol_l <- rep(elist$targets$concentration_umol_l[elist$targets$type != "recovery"], times = sum(nodeframe_elist$toxnode == nodeID))
      concentration_level <- rep(elist$targets$concentration_level[elist$targets$type != "recovery"], times = sum(nodeframe_elist$toxnode == nodeID))
      time_hpe_factor <- ordered(rep(elist$targets$time_hpe[elist$targets$type != "recovery"], times = sum(nodeframe_elist$toxnode == nodeID)))
      time_hpe <- rep(elist$targets$time_hpe[elist$targets$type != "recovery"], times = sum(nodeframe_elist$toxnode == nodeID))
      probe_id <- rep(as.character(nodeframe_elist$ProbeID)[nodeframe_elist$toxnode == nodeID], each = nrow(t(elist$E[elist$genes$ProbeName %in% nodeframe_elist$ProbeID[nodeframe_elist$toxnode == nodeID], elist$targets$type != "recovery"])))
      ensembl_gene_id <- rep(as.character(nodeframe_elist$ensembl)[nodeframe_elist$toxnode == nodeID], each = nrow(t(elist$E[elist$genes$ProbeName %in% nodeframe_elist$ProbeID[nodeframe_elist$toxnode == nodeID], elist$targets$type != "recovery"])))
      substance <- elist$targets$substance[1]
      probeframe <- data.frame(logFC = logFC, concentration_umol_l = concentration_umol_l, concentration_level = concentration_level, time_hpe_factor = time_hpe_factor, time_hpe = time_hpe, probe_id = probe_id, ensembl_gene_id = ensembl_gene_id, nodeID = nodeID, substance = substance)
      return(probeframe)
    } else {
      return(NA)
    }
  })

  nodelist
}
