#' Fit 3-D spline to time and concentration dependent data
#'
#' @param elist an EList with logFC data
#'
#' @return an EList with smoothed logFC data
#' @export
#'
spline_fit <- function(elist) {
    library("mgcv")
    library("outliers")

    #define probe/concentration/time vectors-----------------------------------
    probes <-
        as.data.frame(elist$E[, elist$targets$type !=
                                  "recovery"])
    conc_all <-
        elist$targets$concentration_umol_l[elist$targets$type !=
                                               "recovery"]
    time_all <-
        ordered(elist$targets$time_hpe[elist$targets$type != "recovery"])

    timen_all <-
        elist$targets$time_hpe[elist$targets$type !=
                                   "recovery"]

    concentrations <- sort(unique(conc_all[conc_all!=0]),decreasing = T)
    annotation <- elist$genes
    probes$ProbeID <- rownames(probes)

    #fit spline--------------------------------------------------------------------
    splineprobes <- apply(probes, 1, function(gene, conc_all,time_all,timen_all) {

        #get concentration vector--------------------------------------------------
        conc <- conc_all
        time <- time_all
        timen <- timen_all

        #remove all meta-data from gene vector-------------------------------------
        gene <-
            as.numeric(gene[!names(gene) == "extProbes" &
                                !names(gene) == "ProbeID" &
                                !names(gene) == "maxProbes" &
                                !names(gene) == "minProbes"])

        #perform quality control for probe-----------------------------------------

        #detect and remove outliers (according to grubbs.test)---------------------
        outlier <-
            as.numeric(outliers::grubbs.test(x = gene)["p.value"])
        outliernew <- outlier

        while (log10(outliernew) < -3) {
            time <- time[-which(gene == outliers::outlier(gene))]
            timen <- timen[-which(gene == outliers::outlier(gene))]
            conc <- conc[-which(gene == outliers::outlier(gene))]
            gene <- gene[-which(gene == outliers::outlier(gene))]
            outliernew <-
                as.numeric(outliers::grubbs.test(x = gene)["p.value"])
        }

        #####create dataframe with explaining variables#####
        expldata <- data.frame(dose = conc, time = timen)
        expldata$logFC <- gene
        expldata$ldose <- log(expldata$dose)
        expldata$ltime <- log(expldata$time)

        #####fit thin plate ("tp") spline############
        what = mgcv::gam(formula = logFC ~ te(ldose, ltime, bs = "tp"),
                         data = expldata[expldata$dose != 0,])

        splinedata <-
            expand.grid(ldose = unique(log(conc_all[conc_all != 0])), ltime = unique(log(timen_all[conc_all !=
                                                                                                       0]))) ## do not consider replication because this already influences the spline fit
        splinedata$logFC <- predict(what, newdata = splinedata)
        return(splinedata$logFC)
    })


    splineprobes_list <- new("Elist",
                             list(
                                 E = t(splineprobes),
                                 targets = expand.grid(
                                     concentration_umol_l = unique(conc_all[conc_all != 0]),
                                     time_hpe = unique(timen_all[conc_all != 0])
                                 ),
                                 genes = elist$genes
                             ))

    return(splineprobes_list)

}