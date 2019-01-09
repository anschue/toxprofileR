#' Fit Apical
#'
#' @param responsedf a dataframe containing effect counts
#' @param plot logical, if plot should be given
#' @param designtime which timepoint should be taken to calculate EC values and designconcentrations
#'
#' @return a list of model fits and designconcentrations
#' @export
fit_apical <-
  function(responsedf,
             plot = T,
             designtime = 96) {
    # define alternative models ----------------------------------

    ## full model
    model_logit <-
      function(t1, t2, concentration) {
        1 / (1 + exp(-t1 - (t2 * log10(
          concentration
        ))))
      }

    model_wb <-
      function(t1, t2, concentration) {
        1 * (1 - exp(-exp(t1 + (
          t2 * log10(concentration)
        ))))
      }

    model_gl <-
      function(t1, t2, t3, concentration) {
        1 / ((1 + exp(-t1 - (
          t2 * log10(concentration)
        )))^t3)
      }

    ## for mle estimation
    model_logit_mle <-
      function(t1, t2) {
        1 / (1 + exp(-t1 - (t2 * log10(
          concentration
        ))))
      }

    model_wb_mle <-
      function(t1, t2) {
        1 * (1 - exp(-exp(t1 + (
          t2 * log10(concentration)
        ))))
      }

    model_gl_mle <-
      function(t1, t2, t3) {
        1 / ((1 + exp(-t1 - (
          t2 * log10(concentration)
        )))^t3)
      }

    ## likelihood function ------------------------------------------
    LL_logit <- function(t1, t2, y = effect_count, size = n_testorganism) {
      mu <- model_logit_mle(t1 = t1, t2 = t2)
      loglike <- -sum(dbinom(
        x = y,
        prob = mu,
        size = size,
        log = T
      )) # get the total negative log likelihood
      return(loglike)
    }

    LL_wb <- function(t1, t2, y = effect_count, size = n_testorganism) {
      mu <- model_wb_mle(t1 = t1, t2 = t2)
      loglike <- -sum(dbinom(
        x = y,
        prob = mu,
        size = size,
        log = T
      )) # get the total negative log likelihood
      return(loglike)
    }

    LL_gl <- function(t1,
                          t2,
                          t3,
                          y = effect_count,
                          size = n_testorganism) {
      mu <- model_gl_mle(
        t1 = t1,
        t2 = t2,
        t3 = t3
      )
      loglike <- -sum(dbinom(
        x = y,
        prob = mu,
        size = size,
        log = T
      )) # get the total negative log likelihood
      return(loglike)
    }

    ## reverse functions --------------------------------------------------------------------
    reverse_logit <-
      function(t1, t2, effect) {
        10^((-t1 - log((1 / (
          effect
        )) - 1)) / t2)
      }
    reverse_wb <-
      function(t1, t2, effect) {
        10^((log(-log(1 - ((effect) / (1)
        ))) - t1) / t2)
      }
    reverse_gl <-
      function(t1, t2, t3, effect) {
        10^((-t1 - log((((1) / (effect))^(1 / t3)
        ) - 1)) / t2)
      }


    ECx_ml <-
      as.data.frame(matrix(nrow = length(
        unique(responsedf$exposure_end_hpf)
      ), ncol = length(seq(0, 1, 0.001))), row.names = make.names(unique(responsedf$exposure_end_hpf)))
    colnames(ECx_ml) <- seq(0, 1, 0.001)

    par(mar = c(5.1, 5.1, 1, 1))
    plot((
      responsedf$effect_lethal_count / responsedf$n_testorganism
    ) * 100 ~ responsedf$concentration_umol_L,
    log = "x",
    ylim = c(0, 100),
    type = "n",
    main = "",
    xlab = "Concentration [µmol/l]",
    ylab = "Percent Effect",
    cex.lab = 2,
    cex.axis = 1.5
    )


    bestmodels_ml <- list()
    for (t in c(1:length(unique(responsedf$exposure_end_hpf)))) {
      df_t <-
        responsedf[responsedf$exposure_end_hpf == sort(unique(responsedf$exposure_end_hpf))[t], ]

      points((df_t$effect_lethal_count / df_t$n_testorganism) * 100 ~ df_t$concentration_umol_L,
        col = t,
        pch = t,
        cex = 1.5,
        lwd = 2
      )

      df_t$percent_effect <-
        (df_t$effect_lethal_count / df_t$n_testorganism)

      ### Maximum likelihood############

      concentration <- df_t$concentration_umol_L
      df_t$effect_count <- df_t$effect_lethal_count

      fit_logit_ml <- tryCatch({
        bbmle::mle2(
          minuslogl = LL_logit,
          start = list(t1 = -1, t2 = -1),
          control = list(maxit = 5000),
          data = df_t
        )
      }, error = function(cond) {
        return(NA)
      })

      fit_wb_ml <- tryCatch({
        bbmle::mle2(
          minuslogl = LL_wb,
          start = list(t1 = 10, t2 = 10),
          control = list(maxit = 5000),
          data = df_t
        )
      }, error = function(cond) {
        return(NA)
      })

      fit_gl_ml <- tryCatch({
        bbmle::mle2(
          minuslogl = LL_gl,
          start = list(
            t1 = 10,
            t2 = 10,
            t3 = 0.01
          ),
          control = list(maxit = 5000),
          data = df_t
        )
      }, error = function(cond) {
        return(NA)
      })

      modellist_ml <- list(fit_logit_ml, fit_wb_ml, fit_gl_ml)

      AIC_logit_ml <-
        tryCatch({
          bbmle::AICc(fit_logit_ml, nobs = length(df_t$percent_effect[df_t$concentration_umol_L >
            0]))
        }, error = function(cond) {
          return(NA)
        })
      AIC_wb_ml <-
        tryCatch({
          bbmle::AICc(fit_wb_ml, nobs = length(df_t$percent_effect[df_t$concentration_umol_L >
            0]))
        }, error = function(cond) {
          return(NA)
        })
      AIC_gl_ml <-
        tryCatch({
          bbmle::AICc(fit_gl_ml, nobs = length(df_t$percent_effect[df_t$concentration_umol_L >
            0]))
        }, error = function(cond) {
          return(NA)
        })

      bestmodel_ml <-
        tryCatch({
          modellist_ml[[which.min(c(AIC_logit_ml, AIC_wb_ml, AIC_gl_ml))]]
        }, error = function(cond) {
          return(NA)
        })


      if (isS4(bestmodel_ml)) {
        if (max(df_t$percent_effect) > 0.5) {
          bestmodels_ml[[t]] <- bestmodel_ml
          modelname_ml <-
            c("Logit", "WB", "GL")[which.min(c(AIC_logit_ml, AIC_wb_ml, AIC_gl_ml))]

          if (modelname_ml == "Logit") {
            lines(
              model_logit(
                t1 = bbmle::coef(fit_logit_ml)["t1"],
                t2 = bbmle::coef(fit_logit_ml)["t2"],
                concentration = seq(
                  min(df_t$concentration_umol_L),
                  max(df_t$concentration_umol_L),
                  length.out = 2000
                )
              ) * 100 ~ seq(
                min(df_t$concentration_umol_L),
                max(df_t$concentration_umol_L),
                length.out = 2000
              ),
              col = t,
              lwd = 3,
              lty = 1
            )
            ECx_ml[t, ] <- reverse_logit(
              t1 = bbmle::coef(fit_logit_ml)["t1"],
              t2 = bbmle::coef(fit_logit_ml)["t2"],
              effect = seq(0, 1, 0.001)
            )
          }

          if (modelname_ml == "WB") {
            lines(
              model_wb(
                t1 = bbmle::coef(fit_wb_ml)["t1"],
                t2 = bbmle::coef(fit_wb_ml)["t2"],
                concentration = seq(
                  min(df_t$concentration_umol_L),
                  max(df_t$concentration_umol_L),
                  length.out = 2000
                )
              ) * 100 ~ seq(
                min(df_t$concentration_umol_L),
                max(df_t$concentration_umol_L),
                length.out = 2000
              ),
              col = t,
              lwd = 3,
              lty = 1
            )
            ECx_ml[t, ] <- reverse_wb(
              t1 = bbmle::coef(fit_wb_ml)["t1"],
              t2 = bbmle::coef(fit_wb_ml)["t2"],
              effect = seq(0, 1, 0.001)
            )
          }

          if (modelname_ml == "GL") {
            lines(
              model_gl(
                t1 = bbmle::coef(fit_gl_ml)["t1"],
                t2 = bbmle::coef(fit_gl_ml)["t2"],
                t3 = bbmle::coef(fit_gl_ml)["t3"],
                concentration = seq(
                  min(df_t$concentration_umol_L),
                  max(df_t$concentration_umol_L),
                  length.out = 2000
                )
              ) * 100 ~ seq(
                min(df_t$concentration_umol_L),
                max(df_t$concentration_umol_L),
                length.out = 2000
              ),
              col = t,
              lty = 1,
              lwd = 3
            )
            ECx_ml[t, ] <- reverse_gl(
              t1 = bbmle::coef(fit_gl_ml)["t1"],
              t2 = bbmle::coef(fit_gl_ml)["t2"],
              t3 = bbmle::coef(fit_gl_ml)["t3"],
              effect = seq(0, 1, 0.001)
            )
          }
        }
      }
    }


    minEC50_lethal <- min(ECx_ml["0.5"], na.rm = T)
    EC25_lethal <- ECx_ml["0.25"][paste0("X", designtime), 1]
    EC05_lethal <- ECx_ml["0.005"][paste0("X", designtime), 1]




    VF1 <- (EC25_lethal / EC05_lethal)^(1 / 6)


    designconcentrations <-
      c(
        EC25_lethal,
        EC25_lethal / (VF1^1),
        EC25_lethal / (VF1^2),
        EC25_lethal / (VF1^4),
        EC25_lethal / (VF1^6)
      )
    abline(v = designconcentrations, lty = 2)
    legend(
      "topleft",
      title = "Exposure window [hpf]",
      legend = c(paste0(
        unique(responsedf$exposure_start_hpf), "-",
        sort(unique(responsedf$exposure_end_hpf)),
        " lethal"
      ), paste0(
        unique(responsedf$exposure_start_hpf), "-",
        sort(unique(responsedf$exposure_end_hpf)),
        " sublethal"
      )),
      col = c(1:length(
        unique(responsedf$exposure_end_hpf)
      )),
      pch = c(1:(2 * length(
        unique(responsedf$exposure_end_hpf)
      ))),
      lwd = 2
    )


    ## sublethal effects ------------------------------------------------------------------



    ECx_ml_sub <-
      as.data.frame(matrix(nrow = length(
        unique(responsedf$exposure_end_hpf)
      ), ncol = length(seq(0, 1, 0.001))), row.names = make.names(unique(responsedf$exposure_end_hpf)))
    colnames(ECx_ml_sub) <- seq(0, 1, 0.001)

    bestmodels_ml_sub <- list()

    for (t in c(1:length(unique(responsedf$exposure_end_hpf)))) {
      df_t <-
        responsedf[responsedf$exposure_end_hpf == sort(unique(responsedf$exposure_end_hpf))[t], ]

      df_t <- df_t[!is.na(df_t$effect_sublethal_count), ]

      points((df_t$effect_sublethal_count / df_t$n_testorganism) * 100 ~ df_t$concentration_umol_L,
        col = t,
        pch = (t + length(unique(responsedf$exposure_end_hpf))),
        cex = 1.5,
        lwd = 2
      )

      df_t$percent_effect <-
        (df_t$effect_sublethal_count / df_t$n_testorganism)

      ### Maximum likelihood############

      concentration <- df_t$concentration_umol_L
      df_t$effect_count <- df_t$effect_sublethal_count

      fit_logit_ml <- tryCatch({
        bbmle::mle2(
          minuslogl = LL_logit,
          start = list(t1 = -1, t2 = -1),
          control = list(maxit = 5000),
          data = df_t
        )
      }, error = function(cond) {
        return(NA)
      })
      fit_wb_ml <- tryCatch({
        bbmle::mle2(
          minuslogl = LL_wb,
          start = list(t1 = 10, t2 = 10),
          control = list(maxit = 5000),
          data = df_t
        )
      }, error = function(cond) {
        return(NA)
      })
      fit_gl_ml <- tryCatch({
        bbmle::mle2(
          minuslogl = LL_gl,
          start = list(
            t1 = 10,
            t2 = 10,
            t3 = 0.01
          ),
          control = list(maxit = 5000),
          data = df_t
        )
      }, error = function(cond) {
        return(NA)
      })

      modellist_ml <- list(fit_logit_ml, fit_wb_ml, fit_gl_ml)

      AIC_logit_ml <-
        tryCatch({
          bbmle::AICc(fit_logit_ml, nobs = length(df_t$percent_effect[df_t$concentration_umol_L >
            0]))
        }, error = function(cond) {
          return(NA)
        })
      AIC_wb_ml <-
        tryCatch({
          bbmle::AICc(fit_wb_ml, nobs = length(df_t$percent_effect[df_t$concentration_umol_L >
            0]))
        }, error = function(cond) {
          return(NA)
        })
      AIC_gl_ml <-
        tryCatch({
          bbmle::AICc(fit_gl_ml, nobs = length(df_t$percent_effect[df_t$concentration_umol_L >
            0]))
        }, error = function(cond) {
          return(NA)
        })

      bestmodel_ml <-
        tryCatch({
          modellist_ml[[which.min(c(AIC_logit_ml, AIC_wb_ml, AIC_gl_ml))]]
        }, error = function(cond) {
          return(NA)
        })


      if (isS4(bestmodel_ml)) {
        if (max(df_t$percent_effect) > 0.5) {
          bestmodels_ml[[t]] <- bestmodel_ml
          modelname_ml <-
            c("Logit", "WB", "GL")[which.min(c(AIC_logit_ml, AIC_wb_ml, AIC_gl_ml))]

          if (modelname_ml == "Logit") {
            lines(
              model_logit(
                t1 = bbmle::coef(fit_logit_ml)["t1"],
                t2 = bbmle::coef(fit_logit_ml)["t2"],
                concentration = seq(
                  min(df_t$concentration_umol_L),
                  max(df_t$concentration_umol_L),
                  length.out = 2000
                )
              ) * 100 ~ seq(
                min(df_t$concentration_umol_L),
                max(df_t$concentration_umol_L),
                length.out = 2000
              ),
              col = t,
              lwd = 3,
              lty = 2
            )
            ECx_ml_sub[t, ] <- reverse_logit(
              t1 = bbmle::coef(fit_logit_ml)["t1"],
              t2 = bbmle::coef(fit_logit_ml)["t2"],
              effect = seq(0, 1, 0.001)
            )
          }

          if (modelname_ml == "WB") {
            lines(
              model_wb(
                t1 = bbmle::coef(fit_wb_ml)["t1"],
                t2 = bbmle::coef(fit_wb_ml)["t2"],
                concentration = seq(
                  min(df_t$concentration_umol_L),
                  max(df_t$concentration_umol_L),
                  length.out = 2000
                )
              ) * 100 ~ seq(
                min(df_t$concentration_umol_L),
                max(df_t$concentration_umol_L),
                length.out = 2000
              ),
              col = t,
              lwd = 3,
              lty = 2
            )
            ECx_ml_sub[t, ] <- reverse_wb(
              t1 = bbmle::coef(fit_wb_ml)["t1"],
              t2 = bbmle::coef(fit_wb_ml)["t2"],
              effect = seq(0, 1, 0.001)
            )
          }

          if (modelname_ml == "GL") {
            lines(
              model_gl(
                t1 = bbmle::coef(fit_gl_ml)["t1"],
                t2 = bbmle::coef(fit_gl_ml)["t2"],
                t3 = bbmle::coef(fit_gl_ml)["t3"],
                concentration = seq(
                  min(df_t$concentration_umol_L),
                  max(df_t$concentration_umol_L),
                  length.out = 2000
                )
              ) * 100 ~ seq(
                min(df_t$concentration_umol_L),
                max(df_t$concentration_umol_L),
                length.out = 2000
              ),
              col = t,
              lty = 2,
              lwd = 3
            )
            ECx_ml_sub[t, ] <- reverse_gl(
              t1 = bbmle::coef(fit_gl_ml)["t1"],
              t2 = bbmle::coef(fit_gl_ml)["t2"],
              t3 = bbmle::coef(fit_gl_ml)["t3"],
              effect = seq(0, 1, 0.001)
            )
          }
        }
      }
    }


    minEC50_sub <- tryCatch({
      min(ECx_ml_sub["0.5"], na.rm = T)
    }, error = function(cond) {
      return(minEC50_lethal)
    })
    EC25_sub <- ECx_ml_sub["0.25"][paste0("X", designtime), 1]
    EC05_sub <- ECx_ml_sub["0.005"][paste0("X", designtime), 1]



    return(
      list(
        bestmodels = bestmodels_ml,
        ECx_lethal = ECx_ml,
        ECx_sublethal = ECx_ml_sub,
        minEC50_lethal = minEC50_lethal,
        minEC50_sublethal = minEC50_sub,
        EC25_lethal = EC25_lethal,
        EC05 = EC05_lethal,
        designconcentrations = designconcentrations
      )
    )
  }
