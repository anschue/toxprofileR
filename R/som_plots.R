#' Plot Nodecode
#'
#' @param nodelist A nodelist (created by toxprofileR::create_nodelist())
#' @param tox_universe A tox_universe (created by toxprofileR::create_universe())
#' @param nodeIDs Numerical, a vector of nodeIDs to be plotted
#'
#' @return Plots nodecodes and data
#' @export
plot_nodecode <- function(nodelist, tox_universe, nodeIDs) {
  library("cowplot")

  lapply(nodeIDs, function(nodeID) {
    node_data_raw <- nodelist[[nodeID]]

    node_data <- aggregate.data.frame(node_data_raw$logFC,
      by = list(
        concentration_umol_l = node_data_raw$concentration_umol_l,
        concentration_level = node_data_raw$concentration_level,
        time_hpe = node_data_raw$time_hpe,
        probe_id = node_data_raw$probe_id,
        ensembl_gene_id = node_data_raw$ensembl_gene_id,
        nodeID = node_data_raw$nodeID,
        substance = node_data_raw$substance
      ),
      FUN = median
    )

    colnames(node_data)[colnames(node_data) == "x"] <- "logFC"

    node_data$logFC_code <- as.numeric(t(as.data.frame(tox_universe$som_model$codes[[1]][, match(paste(node_data$substance, node_data$concentration_level, node_data$time_hpe, sep = "_"), colnames(tox_universe$som_model$codes[[1]]))])[node_data$nodeID[1], ]))

    nodeplot <- ggplot2::ggplot(data = node_data, aes(x = concentration_umol_l, y = logFC, group = interaction(probe_id, time_hpe), colour = factor(time_hpe))) +
      geom_line(lwd = 0.2, lty = "dashed") +
      geom_point() +
      geom_line(aes(x = concentration_umol_l, y = logFC_code, group = interaction(probe_id, time_hpe), colour = factor(time_hpe)), lwd = 1) +
      xlab("Concentration [µmol/L]") +
      ylab("logFC") +
      theme_bw() +
      scale_colour_discrete(name = "Time point [hpe]") +
      ggtitle(paste0("toxnode #", nodeID))

    print(nodeplot)
  })


  return(NULL)
}



#' Plot noderesponse
#'
#' @param dslist a list of ELists containing logFCs
#' @param tcta_list a list of parameter dataframes
#' @param nodeframe universe nodeframe
#' @param nodeID ID of node to be plotted
#' @param plot3D logical, should 3D plots be given
#'
#' @return returns plots of measured and fitted noderesponse
#' @export
plot_noderesponse <- function(dslist, tcta_list, nodeframe, nodeID, plot3D = TRUE, output = "plot"){

    breaksfunction<-function(xlim){
        df<-(xlim[2]/xlim[1])^(1/4)
        breaks<-xlim[2]/df^(c(1,2,3))
        return(sort(round(breaks,digits = 1)))
    }

    nodeplotlist_all <- lapply(X = names(dslist), FUN = function(substance) {
        if (sum(nodeframe$toxnode == nodeID) > 0) {

            # aggregate measure data
            elist <- dslist[[substance]]
            logFC <- c(t(elist$E[elist$genes$ProbeName %in% nodeframe$ProbeID[nodeframe$toxnode == nodeID], elist$targets$type != "recovery"]))
            concentration_umol_l <- rep(elist$targets$concentration_umol_l[elist$targets$type != "recovery"], times = sum(nodeframe$toxnode == nodeID))
            concentration_level <- rep(elist$targets$concentration_level[elist$targets$type != "recovery"], times = sum(nodeframe$toxnode == nodeID))
            time_hpe_factor <- ordered(rep(elist$targets$time_hpe[elist$targets$type != "recovery"], times = sum(nodeframe$toxnode == nodeID)))
            time_hpe <- rep(elist$targets$time_hpe[elist$targets$type != "recovery"], times = sum(nodeframe$toxnode == nodeID))
            probe_id <- rep(as.character(nodeframe$ProbeID)[nodeframe$toxnode == nodeID], each = nrow(t(elist$E[elist$genes$ProbeName %in% nodeframe$ProbeID[nodeframe$toxnode == nodeID], elist$targets$type != "recovery"])))
            ensembl_gene_id <- rep(as.character(nodeframe$ensembl)[nodeframe$toxnode == nodeID], each = nrow(t(elist$E[elist$genes$ProbeName %in% nodeframe$ProbeID[nodeframe$toxnode == nodeID], elist$targets$type != "recovery"])))
            substance <- elist$targets$substance[1]
            D_measured <- data.frame(logFC = logFC, concentration_umol_l = concentration_umol_l, concentration_level = concentration_level, time_hpe_factor = time_hpe_factor, time_hpe = time_hpe, probe_id = probe_id, ensembl_gene_id = ensembl_gene_id, nodeID = nodeID, substance = substance)

            conc_all <-
                elist$targets$concentration_umol_l[elist$targets$type != "recovery"]
            time_all <-
                ordered(elist$targets$time_hpe[elist$targets$type != "recovery"])
            timen_all <-
                elist$targets$time_hpe[elist$targets$type != "recovery"]

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
                cbind(ControlCIs$time_hpe,
                      as.data.frame(ControlCIs$x))
            colnames(ControlCIs) <- c("time_hpe", "min", "max")
            ControlCIs$substance <- substance

            if(plot3D == TRUE){
            ###3D
            D_fit_3D <- data.frame(logFC=NA,
                                    time_hpe=expand.grid(seq(3,72,length.out = 50),
                                                         seq(min(concentrations), max(concentrations),length.out = 50))[,1],
                                    concentration_umol_l = expand.grid(seq(3,72,length.out = 50),
                                                                       seq(min(concentrations), max(concentrations), length.out = 50))[,2],
                                    substance = substance)

             D_fit_3D$logFC <- toxprofileR::hill_gauss(dose = D_fit_3D$concentration_umol_l,
                                                       time = D_fit_3D$time_hpe,
                                                       hillslope = as.numeric(tcta_list[[substance]][nodeID, "hillslope_best_hill"]),
                                                       maxS50 = as.numeric(tcta_list[[substance]][nodeID,"maxS50_best_hill"]),
                                                       mu = as.numeric(tcta_list[[substance]][nodeID,"mu_best_hill"]),
                                                       sigma = as.numeric(tcta_list[[substance]][nodeID,"sigma_best_hill"]),
                                                       maxGene = as.numeric(tcta_list[[substance]][nodeID,"max_best_hill"]))
            }

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

            if(plot3D){D_fit_3D$nodeID<-nodeID}
            D_fit$nodeID <- nodeID
            D_measured$nodeID <- nodeID
            ControlCIs$nodeID <- nodeID


            if(plot3D){
                return(list(D_fit = D_fit,
                            D_fit_3D = D_fit_3D,
                            D_measured = D_measured,
                            ControlCIs = ControlCIs))
            } else {
                return(list(D_fit = D_fit,
                            D_measured = D_measured,
                            ControlCIs = ControlCIs))
            }


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

    if(plot3D){
    D_fit_3D_all <-
        do.call("rbind", lapply(nodeplotlist_all, function(x) {
            x["D_fit_3D"][[1]]
        }))
    D_fit_3D_all$substance <- as.character(D_fit_3D_all$substance)
    D_fit_3D_all <- D_fit_3D_all[!is.na(D_fit_3D_all$substance), ]
    }

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

    ###plot

    poiplot_hill <-
        ggplot2::ggplot(data = D_measured_all, aes(x = concentration_umol_l, y = logFC)) +
        geom_point(
                   size = 1,
                   lwd = 0) +
        facet_wrap( ~substance+factor(time_hpe), nrow = 3, scales = "free_x") +
        geom_line(data = D_fit_all,
                  aes(x = concentration_umol_l, y = logFC_hill),
                  lwd = 1.5) +
        geom_ribbon(
            data = D_fit_all,
            aes(
                x = concentration_umol_l,
                ymin = lwr_hill,
                ymax = upr_hill
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

        #theme_bw()+
        ylab("logFC") +
        xlab(expression("exposure concentration [" * mu * "M]")) +
        #xlab("\n\nexposure concentration [µM]") +
        geom_hline(aes(yintercept = 0)) +
        scale_x_log10(breaks = breaksfunction) +
        #ggtitle(paste(POI,"hill"))+
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

    poiplot_gauss <-
        ggplot2::ggplot(data = D_measured_all, aes(x = concentration_umol_l, y = logFC)) +
        geom_point(
            size = 1,
            lwd = 0) +
        facet_wrap( ~substance+factor(time_hpe), nrow = 3, scales = "free_x") +
        geom_line(data = D_fit_all,
                  aes(x = concentration_umol_l, y = logFC_gauss),
                  lwd = 1.5) +
        geom_ribbon(
            data = D_fit_all,
            aes(
                x = concentration_umol_l,
                ymin = lwr_gauss,
                ymax = upr_gauss
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

        #theme_bw()+
        ylab("logFC") +
        xlab(expression("exposure concentration [" * mu * "M]")) +
        #xlab("\n\nexposure concentration [µM]") +
        geom_hline(aes(yintercept = 0)) +
        scale_x_log10(breaks = breaksfunction) +
        #ggtitle(paste(POI,"hill"))+
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



if(plot3D){
list3D <- list()
for(subst in unique(D_fit_3D_all$substance)){
matrix_3D <- xtabs(logFC~concentration_umol_l+time_hpe, data=D_fit_3D_all[D_fit_3D_all$substance == subst,])
plot3D::persp3D(z = matrix_3D ,x= as.numeric(rownames(matrix_3D)), y=as.numeric(colnames(matrix_3D)), col=NULL, colvar = NULL, facets = F, curtain = F,phi = 30,theta = 50,xlab="concentration",ylab="time",zlab="logFC",main=subst, zlim = c(min(D_fit_3D_all$logFC, na.rm = T), max(D_fit_3D_all$logFC, na.rm = T)), cex.lab = 2, resfac = 1, zmin = 0)
list3D[[subst]] <- recordPlot()
dev.off()
}
}

    if(output == "plot"){
    print(poiplot_hill)
    print(poiplot_gauss)
    list3D
    return(NULL)
    }

    if(output == "data"){

        return(list(hill = poiplot_hill, gauss = poiplot_gauss, list3D = list3D))

    }


}



#' Plot Tox Portrait
#'
#' @param nodelist A nodelist created with toxprofileR::create_nodelist()
#' @param tox_universe Toxicogenomic universe
#' @param grid grid of toxicogenomic universe
#' @param tcta_paramframe Optional, dataframe containing fitted parameter values for each node
#' @param substance Optional, substance name
#' @param time_hpe Time point to be plotted
#' @param concentration_level Concentration level to be plotted
#' @param type character string, type of response to be plotted; Possible values are  "code", "median", "modeled", "parameter"
#' @param parameter Optional, name of parameter to be plotted
#' @param onlysig logical, if reponse should be filtered to signficant values
#' @param output character string, if the output should be plotted ("plot") or plot data should be given as output ("data")
#' @param legend logical, should a legend be given out (default: F)
#' @param siglevel vector of significant effect levels
#' @param clip value, above which parameter values should be clipped
#'
#' @return either a ggplot printed or ggplot data, depending on the parameter "output"
#' @export
#'
plot_portrait <- function(nodelist, tox_universe = NULL, grid = NULL, tcta_paramframe = NULL, substance = NULL, time_hpe, concentration_level, type = c("code","median","modeled","parameter"), parameter = NULL, onlysig = c(TRUE, FALSE), siglevel= NULL, output = c("plot","data"), legend = FALSE, clip = NULL, colvec = NULL){
    library("cowplot")

    hill_gauss<-function(dose,time,hillslope,maxS50,mu,sigma,maxGene){
        maxGene/(1+exp(-hillslope*(log(dose)-log(1/((maxS50)*exp(-0.5*(((log(time)-log(mu))/sigma)^2)))))))}

    # colvec_cubic<- scales::rescale(seq(-1,1,length.out = 100)^3)
    # colvec_logistic <- scales::rescale(1/(1+exp(-10*(seq(0,1,length.out = 100)-0.1))))
    # colvec_exp <- scales::rescale((seq(0,1,length.out = 100)^5))

    if(type == "distances"){
        plotdata <- data.frame(distance = unlist(lapply(seq_len(max(tox_universe$som_model$unit.classif)),function(x){
            if(sum(tox_universe$som_model$unit.classif==x)>0){
                median(dist(tox_universe$som_model$data[[1]][tox_universe$som_model$unit.classif==x,],method="manhattan"),na.rm = T)}else{NA}
        })),
        x = tox_universe$som_model$grid$pts[,1],
        y = tox_universe$som_model$grid$pts[,2]
        )


        p1 <- ggplot(plotdata, aes(x,y)) +
            geom_point(aes(size=distance,colour=distance))+
            labs(x="", y="")+
            scale_colour_distiller(palette = "Blues", direction = 1)+
            scale_size(range=c(0.2,2))+
        theme_bw()

        if(output == "plot"){
            print(p1)
            return(NULL)
        } else{
            return(p1)
        }
    }

    if(type == "n"){

        plotdata <- data.frame(n = unlist(lapply(seq_len(max(tox_universe$som_model$unit.classif)),function(x){
            sum(tox_universe$som_model$unit.classif==x)
        })),
        x = tox_universe$som_model$grid$pts[,1],
        y = tox_universe$som_model$grid$pts[,2]
        )


        p1 <- ggplot(plotdata, aes(x,y)) +
            geom_point(aes(size=n,colour=n))+
            labs(x="", y="")+
            scale_colour_distiller(palette = "Blues",direction = 1)+
            scale_size(range=c(0.2,2))+
        theme_bw()

        p2 <- ggplot(plotdata)+geom_histogram(aes(x = n))

        if(output == "plot"){
            print(p1)
            return(NULL)
        } else{
            return(p1)
        }
    }

    if(type == "parameter"){
        plotdata <- data.frame(value = unlist(tcta_paramframe[, parameter]),
                               x = grid$pts[,1],
                               y = grid$pts[,2]
        )

        colvec <- seq(0,1,length.out = 1000)
        if(onlysig){plotdata$value[siglevel==0]<-0}
        # if(!is.null(clip)){
        #     roundvalue <- round(clip/(max(plotdata$value)-min(plotdata$value))*1000, digits = 0)
        #     colvec[1:(1000-roundvalue)]<-seq(0,0.05,length.out = (1000-roundvalue))
        #     colvec[(1000-roundvalue):1000] <- seq(0.05, 1, length.out = (roundvalue+1))
        # }


        p1 <- ggplot(plotdata, aes(x,y)) +
            geom_point(aes(size=abs(siglevel),colour=value))+
            labs(x="", y="")+
            scale_colour_gradientn(colours = c("white","goldenrod","darkorange","firebrick","sienna","brown","coral4"),
                                   values = sort(c(0,clip/4,clip/2,clip,clip*2,clip*4, max(plotdata$value, na.rm = T)))/max(plotdata$value, na.rm = T))+
            scale_size("sum(CI)",range=c(1,3),limits = c(0, 40))+
            theme_bw()


        if(output == "plot"){
            print(p1)
            return(NULL)
        } else{
            return(p1)
        }
    }


    if(type=="code"){
        if(paste(substance, concentration_level, time_hpe ,sep="_")%in%colnames(tox_universe$som_model$codes[[1]])){
        plotdata <- data.frame(logFC = as.numeric(tox_universe$som_model$codes[[1]][,match(paste(substance, concentration_level, time_hpe ,sep="_"),colnames(tox_universe$som_model$codes[[1]]))]),
                               x=tox_universe$som_model$grid$pts[,1],
                               y=tox_universe$som_model$grid$pts[,2])
        } else {
            message("Given conditions not available in toxicogenomic universe, please choose other conditions or other plot type")
        }
    }

    # if(type=="modeled"){
    #     plotdata <- data.frame(logFC=abund.funct(dose = as.numeric(datalist$targets$Concentration..micromol.L.[datalist$targets$Concentration.level==conc][1]),time = timepoint,hillslope = as.numeric(paramlist[,"hillslope_final"]),maxS50 = as.numeric(paramlist[,"maxS50_final"]),mu = as.numeric(paramlist[,"mu_final"]),sigma = as.numeric(paramlist[,"sigma_final"]),maxGene = as.numeric(paramlist[,"maxGene_final_hill"])),x=som_model$grid$pts[,1],y=som_model$grid$pts[,2])
    # }
    #




    if(type == "median"){

        plotdata <- data.frame(logFC = unlist(lapply(nodelist, function(nodedf){
            if(is.data.frame(nodedf)){
            nodedf.agg <- aggregate.data.frame(nodedf$logFC, by = list(concentration_level=nodedf$concentration_level, time_hpe = nodedf$time_hpe),FUN = median)
            nodedf.agg$x[nodedf.agg$concentration_level == concentration_level&nodedf.agg$time_hpe == time_hpe]
            } else {NA}
            })),
            x = grid$pts[,1],
            y = grid$pts[,2])
    }

    # clip logFC greater 5/smaller -5
    if(sum(abs(plotdata$logFC)>5, na.rm = T)>0){
        message("clipping logFC greater 5/smaller -5")
        plotdata$logFC[plotdata$logFC>5&!is.na(plotdata$logFC)] <- 5
        plotdata$logFC[plotdata$logFC< -5&!is.na(plotdata$logFC)] <- -5
    }

    p1 <- ggplot(plotdata, aes(x,y)) +
    geom_point(aes(size=abs(logFC),colour=logFC))+
        labs(x="", y="")+
        scale_colour_distiller(palette = "RdBu",direction = -1, limits=c(-5,5), values = colvec)+
        scale_size(range=c(0.1, 3), limits = c(0,5))
        theme_bw()

    if(legend){
    p <- p1+ theme(
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.margin = unit(c(0,0,0,0), "cm"))
      }else{
       p <- p1 + theme(
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none",
            plot.margin = unit(c(0,0,0,0), "cm"))
      }

    if(output == "plot"){
        print(p)
        return(NULL)
    } else{
        return(p)
    }
}


#' Plot portrait grid
#'
#' @param nodelist A nodelist created with toxprofileR::create_nodelist()
#' @param tox_universe Toxicogenomic universe
#' @param grid grid of toxicogenomic universe
#' @param tcta_paramframe Optional, dataframe containing fitted parameter values for each node
#' @param substance Optional, substance name
#' @param type character string, type of response to be plotted; Possible values are  "code", "median", "modeled", "parameter"
#' @param parameter Optional, name of parameter to be plotted
#' @param filename filename to save plot
#' @param onlysig logical, if reponse should be filtered to signficant values
#'
#' @return Saves a fingerprint grid in the given filename
#' @export
#'
plot_portrait_grid <- function(nodelist, filename, tox_universe = NULL, grid = NULL, tcta_paramframe = NULL, substance = NULL,  type = c("code","median","modeled"), parameter = NULL, onlysig = c(TRUE, FALSE), colvec = NULL){
library("cowplot")

    treatment_frame <- expand.grid(list(concentration_umol_l = sort(unique(nodelist[[1]]$concentration_umol_l[nodelist[[1]]$concentration_umol_l!=0]),decreasing = F), time_hpe = sort(unique(nodelist[[1]]$time_hpe), decreasing = T)))

    plotlist <- lapply(seq(1,nrow(treatment_frame)), function(treatID){
        concentration_umol_l <- treatment_frame$concentration_umol_l[treatID]
        concentration_level <- nodelist[[1]]$concentration_level[nodelist[[1]]$concentration_umol_l == concentration_umol_l][1]
        time_hpe <- treatment_frame$time_hpe[treatID]


        p <- toxprofileR::plot_portrait(nodelist, tox_universe, grid, tcta_paramframe = tcta_paramframe, substance = substance, time_hpe = time_hpe, concentration_level = concentration_level, type = type, parameter = parameter, onlysig = onlysig, output = "data", colvec = colvec)

        return(p)
    })

    concentration_umol_l <- treatment_frame$concentration_umol_l[1]
    concentration_level <- nodelist[[1]]$concentration_level[nodelist[[1]]$concentration_umol_l == concentration_umol_l][1]
    time_hpe <- treatment_frame$time_hpe[1]

    legendplot <- toxprofileR::plot_portrait(nodelist, tox_universe = tox_universe, grid = grid, tcta_paramframe = tcta_paramframe, substance = substance, time_hpe = time_hpe, concentration_level = concentration_level, type = type, parameter = parameter, onlysig = onlysig, output = "data", legend = TRUE)

    som_legend<-cowplot::get_legend(legendplot)

    somplot<-plot_grid(plot_grid(plotlist = plotlist,
                                 nrow = length(unique(nodelist[[1]]$time_hpe)),
                                 ncol= length(unique(nodelist[[1]]$concentration_umol_l[nodelist[[1]]$concentration_umol_l!=0]))),
                       plot_grid(NULL, som_legend, ncol=1),
                       rel_widths=c(length(unique(nodelist[[1]]$concentration_umol_l)), 1))


    save_plot(filename, somplot, nrow = length(unique(nodelist[[1]]$time_hpe)), ncol = (length(unique(nodelist[[1]]$concentration_umol_l[nodelist[[1]]$concentration_umol_l!=0]))+1))



}
