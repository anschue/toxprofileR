#' Preprocess Data
#'
#' @param elist EList for preprocessing
#' @param batchcorrect logical, if batchcorrection should be applied (default: F)
#' @param batch factor, batch information
#'
#' @return a processed EList
#' @export
#'
preprocess <- function(elist, batchcorrect = F, batch) {
  elist.avg <- elist.median(elist)

  map <- new_annotation ## load map created with "unambigousArrayAnnot.R"
  elist.avg$genes <- cbind(elist.avg$genes, map[match(x = elist.avg$genes$ProbeName, table = map$ProbeID), c("ProbeID", "EnsemblID", "wikigene_name", "entrezgene", "ensembl_gene_id", "unique")])

  if (batchcorrect) {
    batch <- batch
    modcombat <- model.matrix(~factor(elist$targets$time_hpe) * factor(elist$targets$concentration_level), data = elist.avg$targets)
    combat_edata <- sva::ComBat(dat = elist.avg$E, batch = batch, mod = modcombat, par.prior = TRUE, prior.plots = F)
    elist.combat <- elist.avg
    elist.combat$E <- combat_edata
    return(elist.combat)
  } else {
    return(elist.avg)
  }
}
