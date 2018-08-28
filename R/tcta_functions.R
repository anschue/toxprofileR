#' Hill-Gauss Model
#'
#' @param dose concentration
#' @param time time after exposure start
#' @param hillslope hillslope
#' @param maxS50 maximum Sensitivity
#' @param mu aka t_max
#' @param sigma sigma
#' @param maxGene maximum logFC
#'
#' @return logFC at given time and concentration with given parameters
#' @export
#'
hill_gauss <- function(dose, time, hillslope, maxS50, mu, sigma, maxGene) {
  maxGene / (1 + exp(-hillslope * (log(dose) - log(1 / ((maxS50) * exp(-0.5 * (((log(time) - log(mu)) / sigma)^2)))))))
}


#' Gauss-Gauss Model
#'
#' @param dose concentration
#' @param time time after exposure start
#' @param mconc mconc
#' @param sconc sconc
#' @param mu t_max
#' @param sigma simga
#' @param maxGene maximum logFC
#'
#' @return logFC at given time and concentration with given parameters
#' @export
#'
gauss_gauss <- function(dose, time, mconc, sconc, mu, sigma, maxGene) {
  maxGene * exp(-((((log(dose) - log(mconc))^2) / (2 * (sconc^2))) + (((log(time) - log(mu))^2) / (2 * (sigma^2)))))
}
