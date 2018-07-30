#' Plot Nodecode
#'
#' @param nodelist A nodelist (created by toxprofileR::create_nodelist())
#' @param tox_universe A tox_universe (created by toxprofileR::create_universe())
#' @param nodeIDs Numerical, a vector of nodeIDs to be plotted
#'
#' @return Plots nodecodes and data
#' @export
plot_nodecode <- function(nodelist, tox_universe, nodeIDs){
    library("cowplot")

    lapply(nodeIDs, function(nodeID){

    node_data_raw <- nodelist[[nodeID]]

    node_data <- aggregate.data.frame(node_data_raw$logFC, by = list(concentration_umol_l = node_data_raw$concentration_umol_l,
                                                                     concentration_level = node_data_raw$concentration_level,
                                                                     time_hpe = node_data_raw$time_hpe,
                                                                     probe_id = node_data_raw$probe_id,
                                                                     ensembl_gene_id = node_data_raw$ensembl_gene_id,
                                                                     nodeID = node_data_raw$nodeID,
                                                                     substance = node_data_raw$substance),
                                      FUN = median)

    colnames(node_data)[colnames(node_data)=="x"] <- "logFC"

    node_data$logFC_code <-  as.numeric(t(as.data.frame(tox_universe$som_model$codes[[1]][,match(paste(node_data$substance, node_data$concentration_level, node_data$time_hpe,sep="_"),colnames(tox_universe$som_model$codes[[1]]))])[node_data$nodeID[1],]))

    nodeplot <- ggplot2::ggplot(data = node_data, aes(x=concentration_umol_l,y=logFC,group=interaction(probe_id,time_hpe),colour=factor(time_hpe)))+
        geom_line(lwd=0.2, lty="dashed")+
        geom_point()+
        geom_line(aes(x=concentration_umol_l,y=logFC_code,group=interaction(probe_id,time_hpe),colour=factor(time_hpe)),lwd=1)+
        xlab("Concentration [µmol/L]")+
        ylab ("logFC")+
        theme_bw()+
        scale_colour_discrete(name="Time point [hpe]")+
        ggtitle(paste0("toxnode #",nodeID))

    print(nodeplot)
    })


    return(NULL)
}

plot_nodecode(nodelist = nodelist_diuron_smooth, tox_universe = tox_universe, nodeIDs = c(1,2,3))

#' Plot Tox Portrait
#'
#' @param nodelist A nodelist created with toxprofileR::create_nodelist()
#' @param tox_universe Toxicogenomic universe
#' @param tcta_paramframe Optional, dataframe containing fitted parameter values for each node
#' @param substance Optional, substance name
#' @param time_hpe Time point to be plotted
#' @param concentration_level Concentration level to be plotted
#' @param type character string, type of response to be plotted; Possible values are  "code", "median", "modeled", "parameter"
#' @param parameter Optional, name of parameter to be plotted
#' @param onlysig logical, if reponse should be filtered to signficant values
#' @param output character string, if the output should be plotted ("plot") or plot data should be given as output ("data")
#' @param legend logical, should a legend be given out (default: F)
#'
#' @return either a ggplot printed or ggplot data, depending on the parameter "output"
#' @export
#'
plot_portrait <- function(nodelist, tox_universe, tcta_paramframe = NULL, substance = NULL, time_hpe, concentration_level, type = c("code","median","modeled","parameter"), parameter = NULL, onlysig = c(TRUE, FALSE), output = c("plot","data"), legend = FALSE){
    library("cowplot")

    hill_gauss<-function(dose,time,hillslope,maxS50,mu,sigma,maxGene){
        maxGene/(1+exp(-hillslope*(log(dose)-log(1/((maxS50)*exp(-0.5*(((log(time)-log(mu))/sigma)^2)))))))}

    colvec_special<-c(0,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,1)


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
            x = tox_universe$som_model$grid$pts[,1],
            y = tox_universe$som_model$grid$pts[,2])
    }

    p1 <- ggplot(plotdata, aes(x,y)) +
    geom_point(aes(size=abs(logFC),colour=logFC))+
        labs(x="", y="")+
        scale_colour_distiller(palette = "RdBu",direction = -1,limits=c(-4,4),values=colvec_special)+
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
#' @param tcta_paramframe Optional, dataframe containing fitted parameter values for each node
#' @param substance Optional, substance name
#' @param type character string, type of response to be plotted; Possible values are  "code", "median", "modeled", "parameter"
#' @param parameter Optional, name of parameter to be plotted
#' @param onlysig logical, if reponse should be filtered to signficant values
#' @param output character string, if the output should be plotted ("plot") or plot data should be given as output ("data")
#'
#' @return
#' @export
#'
plot_portrait_grid <- function(nodelist, tox_universe, tcta_paramframe = NULL, substance = NULL,  type = c("code","median","modeled","parameter"), parameter = NULL, onlysig = c(TRUE, FALSE), output = c("plot", "data")){
library("cowplot")

    treatment_frame <- expand.grid(list(concentration_umol_l = sort(unique(nodelist[[1]]$concentration_umol_l)), time_hpe = sort(unique(nodelist[[1]]$time_hpe))))

    plotlist <- lapply(seq(1,nrow(treatment_frame)), function(treatID){

        concentration_umol_l <- treatment_frame$concentration_umol_l[treatID]
        concentration_level <- nodelist[[1]]$concentration_level[nodelist[[1]]$concentration_umol_l == concentration_umol_l][1]
        time_hpe <- treatment_frame$time_hpe[treatID]

        p <- toxprofileR::plot_portrait(nodelist, tox_universe, tcta_paramframe = tcta_paramframe, substance = substance, time_hpe = time_hpe, concentration_level = concentration_level, type = type, parameter = parameter, onlysig = onlysig, output = "data")

        return(p)
    })

    concentration_umol_l <- treatment_frame$concentration_umol_l[1]
    concentration_level <- nodelist[[1]]$concentration_level[nodelist[[1]]$concentration_umol_l == concentration_umol_l][1]
    time_hpe <- treatment_frame$time_hpe[1]

    legendplot <- toxprofileR::plot_portrait(nodelist, tox_universe, tcta_paramframe = tcta_paramframe, substance = substance, time_hpe = time_hpe, concentration_level = concentration_level, type = type, parameter = parameter, onlysig = onlysig, output = "data", legend = TRUE)

    som_legend<-cowplot::get_legend(legendplot)

    somplot<-plot_grid(plot_grid(plotlist = plotlist,
                                 nrow = length(unique(nodelist[[1]]$time_hpe))),
                                 ncol= length(unique(nodelist[[1]]$concentration_umol_l)),
                       plot_grid(NULL, som_legend, ncol=1),
                       rel_widths=c(length(unique(nodelist[[1]]$concentration_umol_l)), 1))

    if(output == "plot"){
        print(somplot)
        return(NULL)
    } else {
        return(somplot)
    }




}




