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
#'
#' @return A result list containing the som_model, the final input dataset (dataset.SOM), compiled metadata, and ProbeIDs.
#' @export
#'
create_tox_universe <-
  function(dslist,
             dimens = 60,
             dist.fct = "manhattan",
             includearchive = T,
             seed = 2312,
             output = T) {

    # load libraries --------------------------------------------------------------
    library("limma")
    library("kohonen")
    library("tictoc")

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

    # take probes with highest IQR for duplicate probe --------------------
    IQRs <- apply(E_all, 1, IQR, na.rm = T)
    data <- E_all[order(IQRs, decreasing = T), ]
    genes <-
      unlist(as.character(dslist[[1]]$genes$ensembl_gene_id))[order(IQRs,
        decreasing =
          T
      )]
    data <- data[!duplicated(genes), ]
    genes <- genes[!duplicated(genes)]

    data <- data[genes != "NULL", ]
    genes <- genes[genes != "NULL"]

    data <- data[!grepl(pattern = "c", genes), ]
    genes <- genes[!grepl(pattern = "c", genes)]

    ProbeIDs_metalogFC <- data.frame(ProbeID = rownames(data), ensembl = genes)

    rownames(data) <- genes

    # merge own and meta-analysis data ------------------------------------
    dataset.SOM <- merge.data.frame(logFCframe_all, data, by = 0, all = T)
    rownames(dataset.SOM) <- dataset.SOM$Row.names
    dataset.SOM <- dataset.SOM[, -1]
    dataset.SOM <- as.matrix(dataset.SOM)

    ###### initialize SOM###################################
    som_grid <-
      somgrid(
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
      alpha = c(0.5, 0.001),
      radius = as.integer(c(
        quantile(c(0:dimens), 0.67), -quantile(c(0:dimens), 0.67)
      )),
      keep.data = TRUE,
      dist.fcts = "manhattan",
      maxNA.fraction = 0.75
    )
    tictoc::toc()

    resultlist <- list(
      som_model = som_model,
      dataset.SOM = dataset.SOM,
      metadata = metadata,
      comparisons_merge = comparisons_merge,
      ProbeIDs_metalogFC = ProbeIDs_metalogFC
    )

    if (output) {
      plot(som_model, type = "changes")
      plot(som_model, type = "count")
      plot(som_model, type = "dist.neighbours")
      plot(som_model, type = "codes")
    }

    return(resultlist)
  }
