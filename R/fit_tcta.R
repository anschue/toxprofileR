#' Fit regression model on toxicogenomic profile
#'
#' @param elist An EList containing normalized logFC values from one exposure experiment.
#' @param extrema A dataframe containing extreme values for each node
#' @param cluster Logical, if modeling should be parallelized
#'
#' @return A dataframe with fitted model parameters for each node
#' @export
fit_tcta <- function(elist, extrema, cluster = c(TRUE, FALSE), cln = 2){
    library("hydromad")
    library("snow")
    library("outliers")
    library("mgcv")

# create nodelist -------------------------------------------------------------
nodelist <- toxprofileR::create_nodelist(elist)

nodelist_extrema <- lapply(seq(1,length(nodelist)), function(nodeID){
    if(is.data.frame(nodelist[[nodeID]])){
    cbind(nodelist[[nodeID]], do.call("rbind", replicate(n = nrow(nodelist[[nodeID]]), expr = extrema[nodeID,],simplify = F)))
    }else{return(NA)}
})

# set parameter boundaries ----------------------------------------------------
concentration_umol_l<-elist$targets$concentration_umol_l[elist$targets$type!="recovery"]
concentrations<-sort(unique(concentration_umol_l[concentration_umol_l!=0]), decreasing = T)
dilution_factor<-concentrations[1]/concentrations[2]
concrange<-max(concentrations)/min(concentrations)


param_bounds <- list(slope = c(min = -log((1/0.505)-1)/log(concrange^(0.5)),
                               max = -log((1/0.99)-1)/log(dilution_factor^0.5)),
                     sigma = c(min = log(1.5)/((-2*log(0.01))^0.5),
                               max = log(24)/((-2*log(0.99))^0.5)),
                     EC50 = c(min = min(concentrations)/concrange,
                              max = max(concentration_umol_l)*concrange),
                     mS50 = c(min = 1/(max(concentration_umol_l)*concrange),
                              max = 1/(min(concentrations)/concrange)),
                     mu = c(min = 1.5,
                            max = 75),
                     mconc = c(min = min(concentrations)/concrange,
                               max = max(concentration_umol_l)*concrange),
                     sconc = c(min = log(dilution_factor)/((-2*log(0.01))^0.5),
                               max = log(concrange)/((-2*log(0.99))^0.5)),
                     err = c(min = 0,
                             max = 10))


if(cluster){
# check for external cluster --------------------------------------------------
cl<-snow::getMPIcluster()

# if not available set up cluster ---------------------------------------------
if(is.null(cl)){
clustertype <- 1
cl <- snow::makeMPIcluster(cln, type = "SOCK")
} else{ clustertype <- 0}

# load libraries on cluster ---------------------------------------------------
snow::clusterEvalQ(cl = cl, expr = library("limma"))
snow::clusterEvalQ(cl, expr = library("hydromad"))
snow::clusterEvalQ(cl, expr = library("outliers"))
snow::clusterEvalQ(cl, expr = library("mgcv"))
snow::clusterEvalQ(cl, expr = library("toxprofileR"))

# load data on cluster --------------------------------------------------------
snow::clusterExport(cl, c("nodelist_extrema", "param_bounds"))

# apply modeling function -----------------------------------------------------
tictoc::tic()
tcta_paramlist_som <- snow::parLapply(cl = cl, x = nodelist_extrema, fun = toxprofileR::get_tcta_params, param_bounds = param_bounds)
tictoc::toc()
if(clustertype == 1){snow::stopCluster(cl)}
} else {
# apply modeling without paralellization
tictoc::tic()
tcta_paramlist_som <- pbapply::pblapply(X = nodelist_extrema,FUN = toxprofileR::get_tcta_params, param_bounds = param_bounds)
tictoc::toc()
}

tcta_paramframe_som <- as.data.frame(do.call("rbind",tcta_paramlist_som))

tcta_paramframe_som
}






