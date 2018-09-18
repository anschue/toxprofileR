#' Create Shiny Data
#'
#' @param dslist List of logFC data (one EList per substance)
#' @param tcta_list fitted parameters for same substances
#' @param ci_list list of CI differences
#'
#' @return a list of dataframes needed by the shiny app
#' @export
#'
#' @import plotly
create_shiny_data <- function(dslist, tcta_list, ci_list) {
    # 1. Create metadataframe
    message("create metadataframe")

    targets_all <-
        unique(do.call("rbind", lapply(dslist, function(ds) {
            ds$targets
        }))[, c("substance",
                "concentration_level",
                "concentration_umol_l",
                "time_hpe",
                "type")], margin = 1)
    targets_all <- targets_all[targets_all$type != "recovery", ]
    targets_all$name <-
        paste(
            targets_all$substance,
            targets_all$concentration_level,
            targets_all$time_hpe,
            sep = "_"
        )

    if (sum(!names(dslist) %in% targets_all$substance) > 0) {
        stop("substance names of dslist do not align with names in targetfile")
    }

    # 2. create dataframe with all plot data
    message("create map data")

    map_data_all_list <-
        pbapply::pblapply(targets_all$name, function(mapname) {
            # retrieve treatment information
            substance <-
                targets_all$substance[targets_all$name == mapname]
            concentration <-
                targets_all$concentration_level[targets_all$name == mapname]
            time <-
                targets_all$time_hpe[targets_all$name == mapname]

            maplist <- lapply(
                c(1:3600),
                FUN = function(node) {
                    probes <-
                        which(rownames(dslist[[substance]][["E"]]) %in% as.character(nodeframe$ProbeID[nodeframe$toxnode == node]))
                    samples <-
                        which(
                            dslist[[substance]][["targets"]]$time_hpe == time &
                                dslist[[substance]][["targets"]]$concentration_level == concentration
                        )

                    data.frame(
                        mapID = mapname,
                        node = node,
                        x = as.numeric(grid$pts[node, "x"]),
                        y = as.numeric(grid$pts[node, "y"]),
                        logFCmedian = median(dslist[[substance]][["E"]][probes, samples], na.rm = T),
                        logFChill = if (is.data.frame(ci_list[[substance]][[node]])) {
                            ci_list[[substance]][[node]][ci_list[[substance]][[node]][, "concentration_level"] == concentration &
                                                             ci_list[[substance]][[node]][, "time_hpe"] == time, "logFC_hill"]
                        } else{
                            0
                        },
                        cidiff = if (is.data.frame(ci_list[[substance]][[node]])) {
                            ci_list[[substance]][[node]][ci_list[[substance]][[node]][, "concentration_level"] == concentration &
                                                             ci_list[[substance]][[node]][, "time_hpe"] == time, "diff_hill"]
                        } else{
                            0
                        }
                    )

                }
            )

            mapframe <- as.data.frame(do.call("rbind", maplist))

            mapframe
        })

    map_data_all_frame <-
        as.data.frame(do.call("rbind", map_data_all_list))
    map_data_all_frame$mapID <-
        as.character(map_data_all_frame$mapID)


    # create list of mapplots -----------------------------------------------------
    colvec_special <-
        c(0, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 1)
    scale_range <- c(0.2, 4.5)

    map_data_all_frame$dotsize <- abs(map_data_all_frame$cidiff)
    map_data_all_frame$dotsize[is.na(map_data_all_frame$dotsize)] <-
        0


    plotlist <-
        pbapply::pblapply(unique(map_data_all_frame$mapID), function(mapID) {
            toxmap <-
                map_data_all_frame[map_data_all_frame$mapID == mapID,] %>%
                ggplot(aes(x = x, y = y)) +
                geom_point(aes(size = dotsize, colour = logFCmedian)) +
                labs(x = "", y = "") +
                scale_colour_distiller(
                    palette = "RdBu",
                    direction = -1,
                    limits = c(-5, 5),
                    values = colvec_special
                ) +
                scale_size(range = scale_range, limits = c(0, 4)) +
                theme_bw() +
                theme(
                    plot.background = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    legend.position = "none",
                    plot.margin = unit(c(0, 0, 0, 0), "cm")
                )
            return(toxmap)
        })

    names(plotlist) <- unique(map_data_all_frame$mapID)


    # 3. create nodeplotlist -------------------------------------------------------------------
    message("gather data for nodeplots")

    nodeplotlist_all <- lapply(names(dslist), function(substance) {
        message(paste("Creating nodeplotlist for substance", substance, "\n"))

        data_logFC <- dslist[[substance]]
        nodelist <- toxprofileR::create_nodelist(data_logFC)

        conc_all <-
            data_logFC$targets$concentration_umol_l[data_logFC$targets$type != "recovery"]
        time_all <-
            ordered(data_logFC$targets$time_hpe[data_logFC$targets$type != "recovery"])
        timen_all <-
            data_logFC$targets$time_hpe[data_logFC$targets$type != "recovery"]

        concentrations <- sort(unique(conc_all[conc_all != 0]))
        timeconc <-
            expand.grid(
                time_hpe = sort(unique(timen_all)),
                concentration_umol_l = seq(
                    min(concentrations),
                    max(concentrations),
                    length.out = 15
                )
            )

        nodeplotlist <-  pbapply::pblapply(c(1:3600),
                                           function(node) {
                                               if (is.data.frame(nodelist[[node]])) {
                                                   D_measured <- nodelist[[node]]

                                                   # remove outliers -----------------------------------------------------
                                                   outliernew <-
                                                       as.numeric(outliers::grubbs.test(x = as.numeric(D_measured$logFC))["p.value"])

                                                   while (log10(outliernew) < -3) {
                                                       D_measured <-
                                                           D_measured[-which(D_measured$logFC == outliers::outlier(D_measured$logFC)), ]
                                                       outliernew <-
                                                           as.numeric(outliers::grubbs.test(x = as.numeric(D_measured$logFC))["p.value"])
                                                   }


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
                                                       hillslope = as.numeric(tcta_list[[substance]][node, "hillslope_best_hill"]),
                                                       maxS50 = as.numeric(tcta_list[[substance]][node, "maxS50_best_hill"]),
                                                       mu = as.numeric(tcta_list[[substance]][node, "mu_best_hill"]),
                                                       sigma = as.numeric(tcta_list[[substance]][node, "sigma_best_hill"]),
                                                       maxGene = as.numeric(tcta_list[[substance]][node, "max_best_hill"])
                                                   )

                                                   D_fit$lwr_hill <-
                                                       D_fit$logFC_hill - (qnorm(0.95) * as.numeric(tcta_list[[substance]][node, "err_best_hill"]))
                                                   D_fit$upr_hill <-
                                                       D_fit$logFC_hill + (qnorm(0.95) * as.numeric(tcta_list[[substance]][node, "err_best_hill"]))

                                                   # fit for gauss-gauss model ----------------------------------------------------

                                                   D_fit$logFC_gauss <- toxprofileR::gauss_gauss(
                                                       dose = D_fit$concentration_umol_l,
                                                       time = D_fit$time_hpe,
                                                       mconc = as.numeric(tcta_list[[substance]][node, "mconc_best_gauss"]),
                                                       sconc = as.numeric(tcta_list[[substance]][node, "sconc_best_gauss"]),
                                                       mu = as.numeric(tcta_list[[substance]][node, "mu_best_gauss"]),
                                                       sigma = as.numeric(tcta_list[[substance]][node, "sigma_best_gauss"]),
                                                       maxGene = as.numeric(tcta_list[[substance]][node, "max_best_gauss"])
                                                   )

                                                   D_fit$lwr_gauss <-
                                                       D_fit$logFC_gauss - (qnorm(0.95) * as.numeric(tcta_list[[substance]][node, "err_best_gauss"]))
                                                   D_fit$upr_gauss <-
                                                       D_fit$logFC_gauss + (qnorm(0.95) * as.numeric(tcta_list[[substance]][node, "err_best_gauss"]))


                                                   ControlCIs <-
                                                       aggregate(
                                                           D_measured$logFC[D_measured$concentration_umol_l == 0],
                                                           by = list(time_hpe = D_measured$time_hpe[D_measured$concentration_umol_l ==
                                                                                                        0]),
                                                           quantile,
                                                           c(0.025, 0.975)
                                                       )

                                                   ControlCIs <-
                                                       cbind(ControlCIs$time_hpe,
                                                             as.data.frame(ControlCIs$x))
                                                   colnames(ControlCIs) <- c("time_hpe", "min", "max")
                                                   ControlCIs$substance <- substance

                                                   # ###3D
                                                   # D_fit_3D <- data.frame(logFC=NA,
                                                   #                        time_hpe=expand.grid(seq(3,72,length.out = 50),
                                                   #                                             seq(min(concentrations), max(concentrations),length.out = 50))[,1],
                                                   #                        concentration_umol_l = expand.grid(seq(3,72,length.out = 50),
                                                   #                                                           seq(min(concentrations), max(concentrations), length.out = 50))[,2],
                                                   #                        substance = substance)
                                                   #
                                                   # D_fit_3D$logFC <- toxprofileR::hill_gauss(dose = D_fit_3D$concentration_umol_l,
                                                   #                                           time = D_fit_3D$time_hpe,
                                                   #                                           hillslope = as.numeric(tcta_list[[substance]][node, "hillslope_best_hill"]),
                                                   #                                           maxS50 = as.numeric(tcta_list[[substance]][node,"maxS50_best_hill"]),
                                                   #                                           mu = as.numeric(tcta_list[[substance]][node,"mu_best_hill"]),
                                                   #                                           sigma = as.numeric(tcta_list[[substance]][node,"sigma_best_hill"]),
                                                   #                                           maxGene = as.numeric(tcta_list[[substance]][node,"max_best_hill"]))


                                                   D_fit$time_name <-
                                                       factor(
                                                           paste(D_fit$time_hpe, "hpe"),
                                                           levels = c(
                                                               "3 hpe",
                                                               "6 hpe",
                                                               "12 hpe",
                                                               "24 hpe",
                                                               "48 hpe",
                                                               "72 hpe"
                                                           )
                                                       )
                                                   D_measured$time_name <-
                                                       factor(
                                                           paste(D_measured$time_hpe, "hpe"),
                                                           levels = c(
                                                               "3 hpe",
                                                               "6 hpe",
                                                               "12 hpe",
                                                               "24 hpe",
                                                               "48 hpe",
                                                               "72 hpe"
                                                           )
                                                       )
                                                   ControlCIs$time_name <-
                                                       factor(
                                                           paste(ControlCIs$time_hpe, "hpe"),
                                                           levels = c(
                                                               "3 hpe",
                                                               "6 hpe",
                                                               "12 hpe",
                                                               "24 hpe",
                                                               "48 hpe",
                                                               "72 hpe"
                                                           )
                                                       )

                                                   # D_fit_3D$node<-node
                                                   D_fit$node <- node
                                                   D_measured$node <- node
                                                   ControlCIs$node <- node



                                                   out <-
                                                       list(
                                                           D_fit = D_fit,
                                                           D_measured = D_measured,
                                                           ControlCIs = ControlCIs,
                                                           D_fit_3D = D_fit_3D
                                                       )
                                                   return(out)

                                               } else{
                                                   return(NA)
                                               }
                                           })

        return(nodeplotlist)

    })

    D_fit_all <-
        do.call("rbind", lapply(nodeplotlist_all, function(x) {
            do.call("rbind", lapply(x, function(node) {
                node["D_fit"][[1]]
            }))
        }))
    D_fit_all$substance <- as.character(D_fit_all$substance)
    D_fit_all <- D_fit_all[!is.na(D_fit_all$substance), ]

    D_measured_all <-
        do.call("rbind", lapply(nodeplotlist_all, function(x) {
            do.call("rbind", lapply(x, function(node) {
                node["D_measured"][[1]]
            }))
        }))
    D_measured_all$substance <- as.character(D_measured_all$substance)
    D_measured_all <- D_measured_all[!is.na(D_measured_all$substance), ]

    Control_CIs_all <-
        do.call("rbind", lapply(nodeplotlist_all, function(x) {
            do.call("rbind", lapply(x, function(node) {
                node["ControlCIs"][[1]]
            }))
        }))
    Control_CIs_all$substance <- as.character(Control_CIs_all$substance)
    Control_CIs_all <-
        Control_CIs_all[!is.na(Control_CIs_all$substance), ]

    # D_fit_3D_all<-do.call("rbind",lapply(nodeplotlist_all, function(x){do.call("rbind", lapply(x, function(node){node["D_fit_3D"][[1]]}))}))
    # D_fit_3D_all$substance <- as.character(D_fit_3D_all$substance)
    # D_fit_3D_all <- D_fit_3D_all[!is.na(D_fit_3D_all$substance),]

    shiny_data <- list(
        targets_all = targets_all,
        plotlist = plotlist,
        D_measured_all = D_measured_all,
        D_fit_all = D_fit_all,
        #D_fit_3D_all = D_fit_3D_all,
        Control_CIs_all = Control_CIs_all
    )

    return(shiny_data)

}


