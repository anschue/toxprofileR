#' Create toxicogenomic universe
#'
#' @description This function depends on the package "kohonen" to infer a SOM from toxicogenomic data.
#'
#' @param dslist A list of ELists with smoothed logFC data to include in the learning procedure.
#' @param dimens The dimension of the SOM (default: 60)
#' @param dist.fct The distance function to be used (default: manhattan)
#' @param includearchive logical, should published toxicogenomic data be included into SOM (default: T)
#' @param seed The seed to be used for the intial map
#' @param output logical, should analytical plots be given
#' @param alpha learning parameter for SOM (default: 0.8, 0.005)
#' @param logFC_frame data.frame with additional logFC to include for map learning
#'
#' @return A result list containing the som_model, the final input dataset (dataset.SOM), compiled metadata, and ProbeIDs.
#' @export
#'
create_tox_universe <-
  function(dslist,
             logFC_frame,
             dimens = 60,
             dist.fct = "manhattan",
             alpha = c(0.8, 0.005),
             includearchive = T,
             seed = 2312,
             output = T) {

    # check if probes are in the same order ---------------------------------------

    lapply(1:(length(dslist) - 1), function(id) {
      if (!identical(
        dslist[[id]]$genes$ProbeName,
        dslist[[id + 1]]$genes$ProbeName
      )) {
        stop("Probes not in the same order")
      }
    })

    # create one common dataframe--------------------------------------------------
    E_all <- do.call("cbind", lapply(dslist, function(elist) {
      elist$E
    }))

    # create metadata dataframe-----------------------------------------------------
    metadata <-
      data.frame(
        substance = unlist(lapply(dslist, function(elist) {
          elist$targets$substance
        })),
        time_hpe = unlist(lapply(dslist, function(elist) {
          elist$targets$time_hpe
        })),
        concentration_umol_l = unlist(lapply(dslist, function(elist) {
          elist$targets$concentration_umol_l
        })),
        concentration_level = unlist(lapply(dslist, function(elist) {
          elist$targets$concentration_level
        }))
      )

    metadata$names <- make.names(paste(metadata$substance, metadata$concentration_level, metadata$time_hpe, sep = "_"))

    # take probes with highest IQR for duplicate probe --------------------
    IQRs <- apply(E_all, 1, IQR, na.rm = T)
    data <- E_all[order(dslist[[1]]$genes$unique, IQRs, decreasing = T), ]
    genes <- unlist(as.character(dslist[[1]]$genes$ensembl_gene_id))[order(dslist[[1]]$genes$unique, IQRs, decreasing = T)]

    data <- data[!duplicated(genes), ]
    genes <- genes[!duplicated(genes)]

    data <- data[!is.na(genes), ]
    genes <- genes[!is.na(genes)]

    ProbeIDs_metalogFC <- data.frame(ProbeID = rownames(data), ensembl = genes)

    rownames(data) <- genes

    # merge own and meta-analysis data ------------------------------------
    dataset.SOM <- merge.data.frame(logFC_frame, data, by = 0, all = T)
    rownames(dataset.SOM) <- dataset.SOM$Row.names
    dataset.SOM <- dataset.SOM[, -1]
    dataset.SOM <- as.matrix(dataset.SOM)

    ###### initialize SOM###################################
    som_grid <-
      kohonen::somgrid(
        xdim = dimens,
        ydim = dimens,
        topo = "rectangular",
        neighbourhood.fct = "gaussian"
      )

    ###### get reproducible init##########################################
    nobjects <- nrow(dataset.SOM)
    ncodes <- nrow(som_grid$pts)
    set.seed(seed)
    starters <- sample(1:nobjects, ncodes, replace = FALSE)
    init <-
      lapply(list(dataset.SOM), function(x)
        x[starters, , drop = FALSE])

    ###### create SOM################################
    message("create SOM")
    tictoc::tic()
    som_model <- kohonen::som(
      dataset.SOM,
      grid = som_grid,
      init = init,
      rlen = 1000,
      alpha = alpha,
      radius = as.integer(c(
        quantile(c(0:dimens), 0.67), -quantile(c(0:dimens), 0.67)
      )),
      keep.data = TRUE,
      dist.fcts = "manhattan",
      maxNA.fraction = 0.75
    )
    tictoc::toc()

    nodeframe <- data.frame(ensembl = rownames(som_model$data[[1]]), toxnode = som_model$unit.classif)
    nodeframe <- merge.data.frame(nodeframe, ProbeIDs_metalogFC, by = "ensembl", all.x = T, all.y = F, sort = T)


    resultlist <- list(
      som_model = som_model,
      metadata = metadata,
      nodeframe = nodeframe
    )

    if (output) {
      hist(table(nodeframe$toxnode), main = "Number of genes per toxnode", xlab = "Gene count")
      plot(som_model, type = "changes")
      plot(som_model, type = "count")
      plot(som_model, type = "dist.neighbours")
      plot(som_model, type = "codes")
    }

    return(resultlist)
  }
