#' Get regression model parameters for one node
#'
#' @param nodedf A dataframe containing logFC, time, concentration of a node
#' @param param_bounds A list containing parameter boundaries for hill-gauss and gauss-gauss model
#'
#' @return Returns a named vector with fitted model parameters
#' @export
#'
get_tcta_params <- function(nodedf, param_bounds) {
  if (is.data.frame(nodedf)) {
    nodedf <- toxprofileR::remove_outliers(nodedf)
    nodedf <- nodedf[!is.na(nodedf$logFC), ]

    # define regression models --------------------------------------------

    ## Hill-Gauss ---------------------------
    abund.funct <- function(dose, time, hillslope, maxS50, mu, sigma, maxGene) {
      maxGene / (1 + exp(-hillslope * (log(dose) - log(1 / ((maxS50) * exp(-0.5 * (((log(time) - log(mu)) / sigma)^2)))))))
    }

    ## Gauss-Gaus ---------------------------
    abund.funct_gaugau <- function(dose, time, mconc, sconc, mu, sigma, maxGene) {
      maxGene * exp(-((((log(dose) - log(mconc))^2) / (2 * (sconc^2))) + (((log(time) - log(mu))^2) / (2 * (sigma^2)))))
    }


    # define Likelihood functions -----------------------------------------

    ## Hill-Gauss ---------------------------

    ### upregulation --------------
    Likelihood_theta_up <- function(theta) {
      sum(dnorm(x = nodedf$logFC - abund.funct(dose = nodedf$concentration_umol_l, time = nodedf$time_hpe, hillslope = theta[1], maxS50 = theta[2], mu = theta[3], sigma = theta[4], maxGene = nodedf$max[1]), mean = 0, sd = theta[5], log = T))
    }

    ### downregulation ------------
    Likelihood_theta_down <- function(theta) {
      sum(dnorm(x = nodedf$logFC - abund.funct(dose = nodedf$concentration_umol_l, time = nodedf$time_hpe, hillslope = theta[1], maxS50 = theta[2], mu = theta[3], sigma = theta[4], maxGene = nodedf$min[1]), mean = 0, sd = theta[5], log = T))
    }

    ## Gauss-Gauss --------------------------

    ### upregulation --------------
    Likelihood_theta_up_gau <- function(theta) {
      sum(dnorm(x = nodedf$logFC - abund.funct_gaugau(dose = nodedf$concentration_umol_l, time = nodedf$time_hpe, mconc = theta[1], sconc = theta[2], mu = theta[3], sigma = theta[4], maxGene = nodedf$max[1]), mean = 0, sd = theta[5], log = T))
    }

    ### downregulation ------------
    Likelihood_theta_down_gau <- function(theta) {
      sum(dnorm(x = nodedf$logFC - abund.funct_gaugau(dose = nodedf$concentration_umol_l, time = nodedf$time_hpe, mconc = theta[1], sconc = theta[2], mu = theta[3], sigma = theta[4], maxGene = nodedf$min[1]), mean = 0, sd = theta[5], log = T))
    }


    # define starting parameters ------------------------------------------

    starthillslope <- 1
    startmaxS50 <- 1 / median(unique(nodedf$concentration_umol_l[nodedf$concentration_umol_l != 0]))
    startmu <- median(unique(nodedf$time_hpe))
    startsigma <- 0.3
    starterr <- sd(nodedf$logFC[nodedf$concentration_umol_l == min(nodedf$concentration_umol_l) & nodedf$time_hpe_factor == 3])

    startmconc <- nodedf$concentration_umol_l[which.max(abs(nodedf$logFC))]
    startmu_gau <- nodedf$time_hpe[which.max(abs(nodedf$logFC))]
    startsconc <- 0.3

    initialg <- as.numeric(c(hillslope = starthillslope, maxS50 = startmaxS50, mu = startmu, sigma = startsigma, err = starterr))
    names(initialg) <- c("hillslope", "maxS50", "mu", "sigma", "err")

    initialg_gau <- as.numeric(c(mconc = startmconc, sconc = startsconc, mu = startmu_gau, sigma = startsigma, err = starterr))
    names(initialg_gau) <- c("mconc", "sconc", "mu", "sigma", "err")

    # Parameter estimation ------------------------------------------------

    seeds <- c(1, 12, 123) # perfrom parameter estimation with different random seeds and afterwards select best fit

    ## Hill-Gauss ----------------------------
    ### up --------------
    results_up_hill <- lapply(seq(1, length(seeds)), function(seedID) {
      set.seed(seeds[seedID])
      hydromad::SCEoptim(FUN = Likelihood_theta_up, par = initialg, control = list(fnscale = -1, returnpop = F, ncomplex = 10), lower = as.numeric(c(hillslope = param_bounds$slope["min"], maxS50 = param_bounds$mS50["min"], mu = param_bounds$mu["min"], sigma = param_bounds$sigma["min"], err = param_bounds$err["min"])), upper = as.numeric(c(hillslope = param_bounds$slope["max"], maxS50 = param_bounds$mS50["max"], mu = param_bounds$mu["max"], sigma = param_bounds$sigma["max"], err = param_bounds$err["max"])))
    })

    fit_up_hill <- results_up_hill[[which.min(unlist(lapply(results_up_hill, function(results) {
      results$value
    })))]]

    ### down ------------
    results_down_hill <- lapply(seq(1, length(seeds)), function(seedID) {
      set.seed(seeds[seedID])
      hydromad::SCEoptim(FUN = Likelihood_theta_down, par = initialg, control = list(fnscale = -1, returnpop = F, ncomplex = 10), lower = as.numeric(c(hillslope = param_bounds$slope["min"], maxS50 = param_bounds$mS50["min"], mu = param_bounds$mu["min"], sigma = param_bounds$sigma["min"], err = param_bounds$err["min"])), upper = as.numeric(c(hillslope = param_bounds$slope["max"], maxS50 = param_bounds$mS50["max"], mu = param_bounds$mu["max"], sigma = param_bounds$sigma["max"], err = param_bounds$err["max"])))
    })

    fit_down_hill <- results_down_hill[[which.min(unlist(lapply(results_down_hill, function(results) {
      results$value
    })))]]

    ## Gauss-Gauss --------------------------
    ### up --------------
    results_up_gauss <- lapply(seq(1, length(seeds)), function(seedID) {
      set.seed(seeds[seedID])
      hydromad::SCEoptim(FUN = Likelihood_theta_up_gau, par = initialg_gau, control = list(fnscale = -1, returnpop = F, ncomplex = 10), lower = as.numeric(c(mconc = param_bounds$mconc["min"], sconc = param_bounds$sconc["min"], mu = param_bounds$mu["min"], sigma = param_bounds$sigma["min"], err = param_bounds$err["min"])), upper = as.numeric(c(mconc = param_bounds$mconc["max"], sconc = param_bounds$sconc["max"], mu = param_bounds$mu["max"], sigma = param_bounds$sigma["max"], err = param_bounds$err["max"])))
    })

    fit_up_gauss <- results_up_gauss[[which.min(unlist(lapply(results_up_gauss, function(results) {
      results$value
    })))]]

    ### down ------------
    results_down_gauss <- lapply(seq(1, length(seeds)), function(seedID) {
      set.seed(seeds[seedID])
      hydromad::SCEoptim(FUN = Likelihood_theta_down_gau, par = initialg_gau, control = list(fnscale = -1, returnpop = F, ncomplex = 10), lower = as.numeric(c(mconc = param_bounds$mconc["min"], sconc = param_bounds$sconc["min"], mu = param_bounds$mu["min"], sigma = param_bounds$sigma["min"], err = param_bounds$err["min"])), upper = as.numeric(c(mconc = param_bounds$mconc["max"], sconc = param_bounds$sconc["max"], mu = param_bounds$mu["max"], sigma = param_bounds$sigma["max"], err = param_bounds$err["max"])))
    })

    fit_down_gauss <- results_down_gauss[[which.min(unlist(lapply(results_down_gauss, function(results) {
      results$value
    })))]]


    # retrieve AICs for model ----------------------------------------------

    ## fit spline as "positive reference" ----------------------------------
    nodedf$ldose <- log(nodedf$concentration_umol_l)
    nodedf$ltime <- log(nodedf$time_hpe)
    what <- mgcv::gam(formula = logFC ~ te(ldose, ltime, bs = "tp"), data = nodedf[nodedf$concentration_umol_l != 0, ])

    ## retrieve AICcs for Spline/Nullmodel as reference -------------------
    n <- nrow(nodedf)
    knull <- 2
    kfull <- 5

    AIC_spline <- -2 * as.numeric(logLik(what)) + (2 * (attributes(logLik(what))$df + 1) * n / (n - attributes(logLik(what))$df - 2))
    AIC_null <- -2 * as.numeric(logLik(lm(nodedf$logFC ~ 1))) + (2 * (knull + 1) * n / (n - knull - 2))

    ## retrieve AICs for model fits ---------------------------------------
    AIC_up_hill <- -2 * as.numeric(-fit_up_hill$value) + (2 * (kfull + 1) * n / (n - kfull - 2))
    AIC_up_gauss <- -2 * as.numeric(-fit_up_gauss$value) + (2 * (kfull + 1) * n / (n - kfull - 2))
    AIC_down_hill <- -2 * as.numeric(-fit_down_hill$value) + (2 * (kfull + 1) * n / (n - kfull - 2))
    AIC_down_gauss <- -2 * as.numeric(-fit_down_gauss$value) + (2 * (kfull + 1) * n / (n - kfull - 2))

    # compare AICs and calculate AICweights -------------------------------
    sortframe <- data.frame(models = c(AIC_up_hill, AIC_null), index = c(1:2))
    delta <- c(sortframe$models - min(sortframe$models, na.rm = T))
    AICw_up_hill <- (exp(-delta / 2) / sum(exp(-delta / 2), na.rm = T))[1]
    rm(list = c("sortframe", "delta"))

    sortframe <- data.frame(models = c(AIC_up_gauss, AIC_null), index = c(1:2))
    delta <- c(sortframe$models - min(sortframe$models, na.rm = T))
    AICw_up_gauss <- (exp(-delta / 2) / sum(exp(-delta / 2), na.rm = T))[1]
    rm(list = c("sortframe", "delta"))

    sortframe <- data.frame(models = c(AIC_down_hill, AIC_null), index = c(1:2))
    delta <- c(sortframe$models - min(sortframe$models, na.rm = T))
    AICw_down_hill <- (exp(-delta / 2) / sum(exp(-delta / 2), na.rm = T))[1]
    rm(list = c("sortframe", "delta"))

    sortframe <- data.frame(models = c(AIC_down_gauss, AIC_null), index = c(1:2))
    delta <- c(sortframe$models - min(sortframe$models, na.rm = T))
    AICw_down_gauss <- (exp(-delta / 2) / sum(exp(-delta / 2), na.rm = T))[1]
    rm(list = c("sortframe", "delta"))

    sortframe <- data.frame(models = c(AIC_spline, min(AIC_up_hill, AIC_down_hill)), index = c(1:2))
    delta <- c(sortframe$models - min(sortframe$models, na.rm = T))
    AICw_vspline_hill <- (exp(-delta / 2) / sum(exp(-delta / 2), na.rm = T))[2]
    rm(list = c("sortframe", "delta"))

    sortframe <- data.frame(models = c(AIC_spline, min(AIC_up_gauss, AIC_down_gauss)), index = c(1:2))
    delta <- c(sortframe$models - min(sortframe$models, na.rm = T))
    AICw_vspline_gauss <- (exp(-delta / 2) / sum(exp(-delta / 2), na.rm = T))[2]
    rm(list = c("sortframe", "delta"))

    # summarize modeling results ------------------------------------------

    finalmodelcoefs_up_hill <- fit_up_hill$par
    finalmodelcoefs_up_gauss <- fit_up_gauss$par
    finalmodelcoefs_down_hill <- fit_down_hill$par
    finalmodelcoefs_down_gauss <- fit_down_gauss$par

    finalmodelcoefs_best_hill <- unlist(list(c(fit_up_hill$par, max = nodedf$max[1]), c(fit_down_hill$par, max = nodedf$min[1]))[[which.max(c(AICw_up_hill, AICw_down_hill))]])
    finalmodelcoefs_best_gauss <- unlist(list(c(fit_up_gauss$par, max = nodedf$max[1]), c(fit_down_gauss$par, max = nodedf$min[1]))[[which.max(c(AICw_up_gauss, AICw_down_gauss))]])

    names(finalmodelcoefs_best_hill) <- paste0(names(finalmodelcoefs_best_hill), "_best_hill")
    names(finalmodelcoefs_best_gauss) <- paste0(names(finalmodelcoefs_best_gauss), "_best_gauss")

    names(finalmodelcoefs_up_hill) <- paste0(names(finalmodelcoefs_up_hill), "_up_hill")
    names(finalmodelcoefs_down_hill) <- paste0(names(finalmodelcoefs_down_hill), "_down_hill")

    names(finalmodelcoefs_up_gauss) <- paste0(names(finalmodelcoefs_up_gauss), "_up_gauss")
    names(finalmodelcoefs_down_gauss) <- paste0(names(finalmodelcoefs_down_gauss), "_down_gauss")

    finalmodelcoefs <- c(finalmodelcoefs_up_hill,
      finalmodelcoefs_down_hill,
      finalmodelcoefs_best_hill,
      finalmodelcoefs_up_gauss,
      finalmodelcoefs_down_gauss,
      finalmodelcoefs_best_gauss,
      AIC_spline = AIC_spline,
      AIC_null = AIC_null,
      AIC_up_hill = AIC_up_hill,
      AIC_down_hill = AIC_down_hill,
      AIC_up_gauss = AIC_up_gauss,
      AIC_down_gauss = AIC_down_gauss,
      AICw_up_hill = AICw_up_hill,
      AICw_up_gauss = AICw_up_gauss,
      AICw_down_hill = AICw_down_hill,
      AICw_down_gauss = AICw_down_gauss,
      AICw_vspline_hill = AICw_vspline_hill,
      AICw_vspline_gauss = AICw_vspline_gauss,
      Convergence_up_hill = fit_up_hill$convergence,
      Convergence_down_hill = fit_down_hill$convergence,
      Convergence_up_gauss = fit_up_gauss$convergence,
      Convergence_down_gauss = fit_down_gauss$convergence,
      maxGene_up = nodedf$max[1],
      maxGene_down = nodedf$min[1],
      maxGene_up_sub = max(aggregate(nodedf$logFC, by = list(nodedf$concentration_umol_l, nodedf$time_hpe_factor), median)[, "x"], na.rm = T),
      maxGene_down_sub = min(aggregate(nodedf$logFC, by = list(nodedf$concentration_umol_l, nodedf$time_hpe_factor), median)[, "x"], na.rm = T),
      direction_hill = c(1, -1)[which.max(c(AICw_up_hill, AICw_down_hill))],
      direction_gauss = c(1, -1)[which.max(c(AICw_up_gauss, AICw_down_gauss))],
      best_model = c(1, 1, 2, 2)[which.max(c(AICw_up_hill, AICw_down_hill, AICw_up_gauss, AICw_down_gauss))],
      best_model_dir = c(1, -1, 1, -1)[which.max(c(AICw_up_hill, AICw_down_hill, AICw_up_gauss, AICw_down_gauss))]
    )

    if (max(c(AICw_up_hill, AICw_down_hill)) == max(c(AICw_up_gauss, AICw_down_gauss))) {
      finalmodelcoefs$best_model <- 0
    }

    return(finalmodelcoefs)
  } else {
    NA
  }
}
