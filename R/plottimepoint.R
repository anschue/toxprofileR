#' Plot noderesponse for selected timepoint
#'
#' @param dslist a list of ELists containing logFCs
#' @param tcta_list a list of parameter dataframes
#' @param nodeframe universe nodeframe
#' @param timeselect
#' @param concselect
#' @param mixtureranges
#' @param output
#' @param logx
#' @param plotCIs
#' @param nodeID ID of node to be plotted
#'
#' @return returns plots of measured and fitted noderesponse
#' @export
plot_noderesponse_time <- function(dslist, tcta_list, nodeframe, nodeID, timeselect, concselect, mixtureranges = NULL, output = "plot",logx = F, plotCIs = F) {
    breaksfunction <- function(xlim) {
        df <- (xlim[2] / xlim[1])^(1 / 4)
        breaks <- xlim[2] / df^(c(1, 2, 3))
        return(sort(round(breaks, digits = 1)))
    }

    nodeplotlist_all <- lapply(X = names(dslist), FUN = function(substance) {
        if (sum(nodeframe$toxnode == nodeID) > 0) {

            # aggregate measure data
            elist <- dslist[[substance]]
            logFC <- c(t(elist$E[elist$genes$ProbeName %in% nodeframe$ProbeID[nodeframe$toxnode == nodeID], elist$targets$type != "recovery"]))
            concentration_umol_l <- rep(elist$targets$concentration_umol_l[elist$targets$type != "recovery"], times = sum(nodeframe$toxnode == nodeID&!is.na(nodeframe$ProbeID)))
            concentration_level <- rep(elist$targets$concentration_level[elist$targets$type != "recovery"], times = sum(nodeframe$toxnode == nodeID&!is.na(nodeframe$ProbeID)))
            time_hpe_factor <- ordered(rep(elist$targets$time_hpe[elist$targets$type != "recovery"], times = sum(nodeframe$toxnode == nodeID&!is.na(nodeframe$ProbeID))))
            time_hpe <- rep(elist$targets$time_hpe[elist$targets$type != "recovery"], times = sum(nodeframe$toxnode == nodeID&!is.na(nodeframe$ProbeID)))
            probe_id <- rep(as.character(nodeframe$ProbeID)[nodeframe$toxnode == nodeID&!is.na(nodeframe$ProbeID)], each = nrow(t(elist$E[elist$genes$ProbeName %in% nodeframe$ProbeID[nodeframe$toxnode == nodeID], elist$targets$type != "recovery"])))
            ensembl_gene_id <- rep(as.character(nodeframe$ensembl)[nodeframe$toxnode == nodeID&!is.na(nodeframe$ProbeID)], each = nrow(t(elist$E[elist$genes$ProbeName %in% nodeframe$ProbeID[nodeframe$toxnode == nodeID], elist$targets$type != "recovery"])))
            substance <- elist$targets$substance[1]
            D_measured <- data.frame(logFC = logFC, concentration_umol_l = concentration_umol_l, concentration_level = concentration_level, time_hpe_factor = time_hpe_factor, time_hpe = time_hpe, probe_id = probe_id, ensembl_gene_id = ensembl_gene_id, nodeID = nodeID, substance = substance)

            conc_all <-
                elist$targets$concentration_umol_l[elist$targets$type != "recovery"]
            time_all <-
                ordered(elist$targets$time_hpe[elist$targets$type != "recovery"])
            timen_all <-
                elist$targets$time_hpe[elist$targets$type != "recovery"]

            concentrations <- sort(unique(conc_all[conc_all != 0]))


            if(logx){
            timeconc <-
                expand.grid(
                    time_hpe = timeselect,
                    concentration_umol_l = emdbook::lseq(
                        min(concselect[[substance]]),
                        max(concselect[[substance]]),
                        length.out = 15
                    )
                )
            }else{
                timeconc <-
                    expand.grid(
                        time_hpe = timeselect,
                        concentration_umol_l = seq(
                            min(concselect[[substance]]),
                            max(concselect[[substance]]),
                            length.out = 15
                        )
                    )
            }
            # remove outliers -----------------------------------------------------
            outliernew <-
                as.numeric(outliers::grubbs.test(x = as.numeric(D_measured$logFC))["p.value"])

            while (log10(outliernew) < -3) {
                D_measured <-
                    D_measured[-which(D_measured$logFC == outliers::outlier(D_measured$logFC)), ]
                outliernew <-
                    as.numeric(outliers::grubbs.test(x = as.numeric(D_measured$logFC))["p.value"])
            }

            # aggregate fitted data
            D_fit <- data.frame(
                logFC_hill = NA,
                upr_hill = NA,
                lwr_hill = NA,
                logFC_gauss = NA,
                upr_gauss = NA,
                lwr_gauss = NA,
                time_hpe = timeconc$time_hpe,
                concentration_umol_l = timeconc$concentration_umol_l,
                substance = substance
            )

            # fit for hill-gauss model ----------------------------------------------------

            D_fit$logFC_hill <- toxprofileR::hill_gauss(
                dose = D_fit$concentration_umol_l,
                time = D_fit$time_hpe,
                hillslope = as.numeric(tcta_list[[substance]][nodeID, "hillslope_best_hill"]),
                maxS50 = as.numeric(tcta_list[[substance]][nodeID, "maxS50_best_hill"]),
                mu = as.numeric(tcta_list[[substance]][nodeID, "mu_best_hill"]),
                sigma = as.numeric(tcta_list[[substance]][nodeID, "sigma_best_hill"]),
                maxGene = as.numeric(tcta_list[[substance]][nodeID, "max_best_hill"])
            )

            D_fit$lwr_hill <-
                D_fit$logFC_hill - (qnorm(0.95) * as.numeric(tcta_list[[substance]][nodeID, "err_best_hill"]))
            D_fit$upr_hill <-
                D_fit$logFC_hill + (qnorm(0.95) * as.numeric(tcta_list[[substance]][nodeID, "err_best_hill"]))

            # fit for gauss-gauss model ----------------------------------------------------

            D_fit$logFC_gauss <- toxprofileR::gauss_gauss(
                dose = D_fit$concentration_umol_l,
                time = D_fit$time_hpe,
                mconc = as.numeric(tcta_list[[substance]][nodeID, "mconc_best_gauss"]),
                sconc = as.numeric(tcta_list[[substance]][nodeID, "sconc_best_gauss"]),
                mu = as.numeric(tcta_list[[substance]][nodeID, "mu_best_gauss"]),
                sigma = as.numeric(tcta_list[[substance]][nodeID, "sigma_best_gauss"]),
                maxGene = as.numeric(tcta_list[[substance]][nodeID, "max_best_gauss"])
            )

            D_fit$lwr_gauss <-
                D_fit$logFC_gauss - (qnorm(0.95) * as.numeric(tcta_list[[substance]][nodeID, "err_best_gauss"]))
            D_fit$upr_gauss <-
                D_fit$logFC_gauss + (qnorm(0.95) * as.numeric(tcta_list[[substance]][nodeID, "err_best_gauss"]))


            ControlCIs <-
                aggregate(
                    D_measured$logFC[D_measured$concentration_umol_l == 0],
                    by = list(time_hpe = D_measured$time_hpe[D_measured$concentration_umol_l ==
                                                                 0]),
                    quantile,
                    c(0.025, 0.975)
                )

            ControlCIs <-
                cbind(
                    ControlCIs$time_hpe,
                    as.data.frame(ControlCIs$x)
                )
            colnames(ControlCIs) <- c("time_hpe", "min", "max")
            ControlCIs$substance <- substance


           D_fit <- D_fit[D_fit$time_hpe %in% timeselect,]
           D_measured <- D_measured[D_measured$time_hpe %in% timeselect,]
           ControlCIs <- ControlCIs[ControlCIs$time_hpe %in% timeselect,]


            D_fit$time_name <-
                factor(
                    paste(D_fit$time_hpe, "hpe"),
                    levels = paste(sort(unique(D_fit$time_hpe)), "hpe")
                    )

            D_measured$time_name <-
                factor(
                    paste(D_measured$time_hpe, "hpe"),
                    levels = paste(sort(unique(D_measured$time_hpe)), "hpe")

                    )

            ControlCIs$time_name <-
                factor(
                    paste(ControlCIs$time_hpe, "hpe"),
                    levels = paste(sort(unique(ControlCIs$time_hpe)), "hpe")

                )



            D_fit$nodeID <- nodeID
            D_measured$nodeID <- nodeID
            ControlCIs$nodeID <- nodeID



                return(list(
                    D_fit = D_fit,
                    D_measured = D_measured,
                    ControlCIs = ControlCIs
                ))

        } else {
            return(NA)
        }
    })

    D_fit_all <-
        do.call("rbind", lapply(nodeplotlist_all, function(x) {
            x["D_fit"][[1]]
        }))
    D_fit_all$substance <- as.character(D_fit_all$substance)
    D_fit_all <- D_fit_all[!is.na(D_fit_all$substance), ]

    D_measured_all <-
        do.call("rbind", lapply(nodeplotlist_all, function(x) {
            x["D_measured"][[1]]
        }))
    D_measured_all$substance <- as.character(D_measured_all$substance)
    D_measured_all <- D_measured_all[!is.na(D_measured_all$substance), ]

    Control_CIs_all <-
        do.call("rbind", lapply(nodeplotlist_all, function(x) {
            x["ControlCIs"][[1]]
        }))
    Control_CIs_all$substance <- as.character(Control_CIs_all$substance)
    Control_CIs_all <-
        Control_CIs_all[!is.na(Control_CIs_all$substance), ]

    ### plot

    poiplot_hill <-
        ggplot2::ggplot(data = D_measured_all, aes(x = concentration_umol_l, y = logFC)) +
        geom_point(
            size = 1,
            lwd = 0
        ) +
        facet_wrap(~substance + factor(time_hpe), nrow = 3, scales = "free_x") +
        geom_line(
            data = D_fit_all,
            aes(x = concentration_umol_l, y = logFC_hill),
            lwd = 1.5
        ) +
        geom_ribbon(
            data = D_fit_all,
            aes(
                x = concentration_umol_l,
                ymin = lwr_hill,
                ymax = upr_hill,
                fill = substance
            ),
            alpha = 0.3,
            inherit.aes = F
        ) +

        theme(
            legend.position = "right",
            axis.title.y = element_text(size = 16, face = "bold"),
            axis.title.x = element_text(size = 16, face = "bold"),
            axis.text.x = element_text(
                angle = 70,
                vjust = 1,
                hjust = 1,
                size = 14
            ),
            axis.text.y = element_text(size = 14),
            strip.text = element_text(size = 11, face = "bold")
        ) +

        ylab("logFC") +
        xlab(expression("exposure concentration [" * mu * "M]")) +
        # xlab("\n\nexposure concentration [µM]") +
        geom_hline(aes(yintercept = 0))

    if(plotCIs){
        poiplot_hill <- poiplot_hill +
        geom_hline(
            data = Control_CIs_all,
            aes(yintercept = min),
            lty = 2,
            lwd = 0.75
        ) +
        geom_hline(
            data = Control_CIs_all,
            aes(yintercept = max),
            lty = 2,
            lwd = 0.75
        )
    }

    if(logx){
        poiplot_hill <- poiplot_hill +scale_x_log10(breaks = breaksfunction)
    }

    if(!is.null(mixtureranges)){



        poiplot_hill <- poiplot_hill +
            geom_rect(data = mixtureranges, mapping = aes(xmin = x1, ymin = y1, xmax = x2, ymax = y2), alpha = 0.2, size = 0.2, color = "black", fill = "grey", inherit.aes =F)



    }

    if (output == "plot") {
        print(poiplot_hill)
        return(NULL)
    }

    if (output == "data") {
        return(poiplot_hill)
    }
}
