#' Get Difference of Confidence Interval
#'
#' @param tcta_frame dataframe with fitted parameters
#' @param nodelist nodelist of the same experiment
#'
#' @return a list with dataframes for each node containing the CI-diffs for each treatment
#' @export
#'
get_CI_diffs <-
  function(tcta_frame,
             nodelist) {
    overlap_tcta <- lapply(
      X = seq_len(length(nodelist)),
      FUN = function(nodeID) {
        if (is.data.frame(nodelist[[nodeID]])) {
          nodedf <- nodelist[[nodeID]]

          treatmentns <-
            nrow(nodedf[nodedf$concentration_umol_l != 0, ])

          controls <-
            nodedf[nodedf$concentration_umol_l == 0, ]
          controlns <-
            aggregate(controls$logFC,
              by = list(time_hpe = controls$time_hpe),
              length
            )

          # take quantiles --------------------------------------------------------
          ControlCIs <-
            aggregate(
              controls$logFC,
              by = list(time_hpe = controls$time_hpe),
              quantile,
              c(0.025, 0.975)
            )
          ControlCIs <-
            cbind(
              ControlCIs$time_hpe,
              as.data.frame(ControlCIs$x)
            )
          colnames(ControlCIs) <-
            c("time_hpe", "min_hill", "max_hill")
          ControlCIs$n <- controlns$x
          ControlCIs$max_gauss <- ControlCIs$max_hill
          ControlCIs$min_gauss <- ControlCIs$min_hill

          # calculate confidence interval for treatment ---------------------------
          modeldf <-
            nodedf[
              !duplicated(nodedf[, c("concentration_umol_l", "time_hpe")]),
              c(
                "concentration_umol_l",
                "concentration_level",
                "time_hpe"
              )
            ]
          modeldf <-
            modeldf[modeldf$concentration_umol_l != 0, ]

          modeldf$logFC_hill <- NA
          modeldf$logFC_gauss <- NA
          modeldf$n <- treatmentns

          # hill-gauss ------------------------------------------------------------
          modeldf$logFC_hill <-
            toxprofileR::hill_gauss(
              dose = modeldf$concentration_umol_l,
              time = modeldf$time_hpe,
              hillslope = tcta_frame[[nodeID, "hillslope_best_hill"]],
              maxS50 = tcta_frame[[nodeID, "maxS50_best_hill"]],
              mu = tcta_frame[[nodeID, "mu_best_hill"]],
              sigma = tcta_frame[[nodeID, "sigma_best_hill"]],
              maxGene = tcta_frame[[nodeID, "max_best_hill"]]
            )

          modeldf$se_hill <-
            tcta_frame[[nodeID, "err_best_hill"]]
          modeldf$max_hill <-
            modeldf$logFC_hill + (qnorm(0.975) * modeldf$se_hill)
          modeldf$min_hill <-
            modeldf$logFC_hill - (qnorm(0.975) * modeldf$se_hill)



          # gauss-gauss ------------------------------------------------------------
          modeldf$logFC_gauss <-
            toxprofileR::gauss_gauss(
              dose = modeldf$concentration_umol_l,
              time = modeldf$time_hpe,
              mconc = tcta_frame[[nodeID, "mconc_best_gauss"]],
              sconc = tcta_frame[[nodeID, "sconc_best_gauss"]],
              mu = tcta_frame[[nodeID, "mu_best_gauss"]],
              sigma = tcta_frame[[nodeID, "sigma_best_gauss"]],
              maxGene = tcta_frame[[nodeID, "max_best_gauss"]]
            )

          modeldf$se_gauss <-
            tcta_frame[[nodeID, "err_best_gauss"]]
          modeldf$max_gauss <-
            modeldf$logFC_gauss + (qnorm(0.975) * modeldf$se_gauss)
          modeldf$min_gauss <-
            modeldf$logFC_gauss - (qnorm(0.975) * modeldf$se_gauss)


          # calculate difference for each measured treatment -----------------------

          # hill-gauss -----
          modeldf$diff_hill <-
            apply(
              modeldf,
              MARGIN = 1,
              FUN = function(treatment) {
                if (as.numeric(treatment["min_hill"]) > 0) {
                  dif <-
                    as.numeric(treatment["min_hill"]) - ControlCIs$max_hill[ControlCIs$time_hpe ==
                      as.numeric(treatment["time_hpe"])]
                  if (length(dif) == 0) {
                    dif <- 0
                  } # if controls were removed as outliers
                  dif[dif < 0] <- 0
                  return(dif)
                } else {
                  if (as.numeric(treatment["max_hill"]) < 0) {
                    dif <-
                      as.numeric(treatment["max_hill"]) - ControlCIs$min_hill[ControlCIs$time ==
                        as.numeric(treatment["time_hpe"])]
                    if (length(dif) == 0) {
                      dif <- 0
                    } # if controls were removed as outliers
                    dif[dif > 0] <- 0
                    return(dif)
                  } else {
                    return(0)
                  }
                }
              }
            )

          # gauss-gauss ------
          modeldf$diff_gauss <-
            apply(
              modeldf,
              MARGIN = 1,
              FUN = function(treatment) {
                if (as.numeric(treatment["min_gauss"]) > 0) {
                  dif <-
                    as.numeric(treatment["min_gauss"]) - ControlCIs$max_gauss[ControlCIs$time_hpe ==
                      as.numeric(treatment["time_hpe"])]
                  if (length(dif) == 0) {
                    dif <- 0
                  } # if controls were removed as outliers
                  dif[dif < 0] <- 0
                  return(dif)
                } else {
                  if (as.numeric(treatment["max_gauss"]) < 0) {
                    dif <-
                      as.numeric(treatment["max_gauss"]) - ControlCIs$min_gauss[ControlCIs$time ==
                        as.numeric(treatment["time_hpe"])]
                    if (length(dif) == 0) {
                      dif <- 0
                    } # if controls were removed as outliers
                    dif[dif > 0] <- 0
                    return(dif)
                  } else {
                    return(0)
                  }
                }
              }
            )

          # add controls -----
          ControlCIs$concentration_umol_l <- 0
          ControlCIs$concentration_level <- "Control"
          ControlCIs$logFC_hill <- 0
          ControlCIs$diff_hill <- NA
          ControlCIs$se_hill <- NA
          ControlCIs$logFC_gauss <- 0
          ControlCIs$diff_gauss <- NA
          ControlCIs$se_gauss <- NA

          modeldf <- rbind(modeldf, ControlCIs)

          modeldf$node <- nodeID

          return(modeldf)
        } else {
          return(NA)
        }
      }
    )

    return(overlap_tcta)
  }
