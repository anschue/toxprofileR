#' Hill-Gauss Model
#'
#' @param dose
#' @param time
#' @param hillslope
#' @param maxS50
#' @param mu
#' @param sigma
#' @param maxGene
#'
#' @return
#' @export
#'
hill_gauss<-function(dose,time,hillslope,maxS50,mu,sigma,maxGene){
    maxGene/(1+exp(-hillslope*(log(dose)-log(1/((maxS50)*exp(-0.5*(((log(time)-log(mu))/sigma)^2)))))))}


#' Gauss-Gauss Model
#'
#' @param dose
#' @param time
#' @param mconc
#' @param sconc
#' @param mu
#' @param sigma
#' @param maxGene
#'
#' @return
#' @export
#'
gauss_gauss<-function(dose,time,mconc,sconc,mu,sigma,maxGene){
    maxGene*exp(-((((log(dose)-log(mconc))^2)/(2*(sconc^2)))+(((log(time)-log(mu))^2)/(2*(sigma^2)))))}
