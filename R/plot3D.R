#' Plot noderesponse as surface
#'
#' @param nodeID
#' @param zlim
#' @param concentrations
#' @param tcta_paramframe
#' @param timeselect
#' @param surfacecolor
#' @param concselect
#' @param logconc
#'
#' @return
#' @export
#'
plot3D_node <- function(nodeID,
                        zlim,
                        concentrations,
                        tcta_paramframe,
                        timeselect = NULL,
                        concselect = NULL,
                        surfacecolor = "grey",
                        logconc = FALSE
){

D_fit_3D <- data.frame(
    logFC = NA,
    time_hpe = expand.grid(
        seq(3, 72, length.out = 50),
        seq(min(concentrations), max(concentrations), length.out = 50)
    )[, 1],
    concentration_umol_l = expand.grid(
        seq(3, 72, length.out = 50),
        seq(min(concentrations), max(concentrations), length.out = 50)
    )[, 2]

)

D_fit_3D$logFC <- toxprofileR::hill_gauss(
    dose = D_fit_3D$concentration_umol_l,
    time = D_fit_3D$time_hpe,
    hillslope = as.numeric(tcta_paramframe[nodeID, "hillslope_best_hill"]),
    maxS50 = as.numeric(tcta_paramframe[nodeID, "maxS50_best_hill"]),
    mu = as.numeric(tcta_paramframe[nodeID, "mu_best_hill"]),
    sigma = as.numeric(tcta_paramframe[nodeID, "sigma_best_hill"]),
    maxGene = as.numeric(tcta_paramframe[nodeID, "max_best_hill"])
)

D_fit_3D$lconc <- log(D_fit_3D$concentration_umol_l)



if(logconc){
    concentrations <- log(concentrations)
    if(!is.null(concselect)){  concselect <- log(concselect)}
    matrix_3D <- xtabs(logFC ~ lconc + time_hpe, data = D_fit_3D)
    plot3D::persp3D(z = matrix_3D, x = as.numeric(rownames(matrix_3D)), y = as.numeric(colnames(matrix_3D)), col = surfacecolor, colvar = NULL, facets = F, curtain = F, phi = 30, theta = 50, xlab = "log(concentration)", ylab = "time", zlab = "logFC", main = "", zlim = zlim, cex.lab = 2, resfac = 1, zmin = 0, xlim = c(min(concentrations,concselect), max(concentrations,concselect)))
}else{
matrix_3D <- xtabs(logFC ~ concentration_umol_l + time_hpe, data = D_fit_3D)
plot3D::persp3D(z = matrix_3D, x = as.numeric(rownames(matrix_3D)), y = as.numeric(colnames(matrix_3D)), col = surfacecolor, colvar = NULL, facets = F, curtain = F, phi = 30, theta = 50, xlab = "concentration", ylab = "time", zlab = "logFC", main = "", zlim = zlim, cex.lab = 2, resfac = 1, zmin = 0)
}

if(!is.null(timeselect)){

    if(!is.null(concselect)){
        plot3D::rect3D(x0 = concselect[1], y0 = timeselect, z0 = zlim[1], x1 = concselect[2], z1 = zlim[2],
                       bty = "g", facets = TRUE,
                       border = "black", col ="grey", alpha=0.8,
                       lwd = 2, add = T)
    }

plot3D::rect3D(x0 = min(concentrations), y0 = timeselect, z0 = zlim[1], x1 = max(concentrations), z1 = zlim[2],
       bty = "g", facets = TRUE,
       border = "black", col ="grey", alpha=0.5,
       lwd = 2, add = T)



}
plot_3D <- recordPlot()
return(plot_3D)
}
