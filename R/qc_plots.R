#' QC plots
#'
#' @param elist an EList
#'
#' @return Different QC Plots
#' @export
qc_plots <- function(elist) {

  # boxplot----------------------------------------------------------------------
  boxplot(elist$E, ylim = c(0, 20))

  # density plot-----------------------------------------------------------------
  limma::plotDensities(elist$E, legend = F)

  # control probes---------------------------------------------------------------
  boxplot(elist$E[elist$genes$ProbeName == "GE_BrightCorner", ], ylim = c(0, 20), xlab = "Samples", ylab = "log2 Expression Bright Corners")
  boxplot(elist$E[elist$genes$ProbeName == "DarkCorner", ], ylim = c(0, 20), xlab = "Samples", ylab = "log2 Expression Dark Corners")

  # Spike ins--------------------------------------------------------------------
  data <- elist
  x <- c(1:10)
  plot(x, x, ylim = c(0, 20), type = "n", main = "Spike Ins")
  spikes <- c("(+)E1A_r60_3", "(+)E1A_r60_a104", "(+)E1A_r60_a107", "(+)E1A_r60_a135", "(+)E1A_r60_a20", "(+)E1A_r60_a22", "(+)E1A_r60_a97", "(+)E1A_r60_n11", "(+)E1A_r60_n9", "(+)E1A_r60_1")

  for (i in c(1:dim(data$targets)[1])) {
    meds <- rep(NA, 10)
    for (j in c(1:10)) {
      meds[j] <- median(data$E[data$genes$ProbeName == spikes[j], i])
    }

    lines(x, meds, col = i)
  }
  abline(h = seq(0, 20, 0.5), lty = 3, lwd = 0.5)


  # MDS plots--------------------------------------------------------------------

  plotMDS(elist$E, labels = elist$targets$names, col = as.numeric(as.factor(elist$targets$type)), main = "Cloess norm")
}
