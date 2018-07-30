#' Get significance level
#'
#' @param CIdiffs List of CIdiffs
#' @param model "hill-gauss" or "gauss-gauss"
#' @param siglevel level of significant (default: 0)
#'
#' @return a vector with significant levels for the map
#' @export
get_sig_level <- function(CIdiffs, model = c("hill-gauss", "gauss-gauss"), siglevel = 0){

    if(model == "hill-gauss"){
    sig_level <- unlist(lapply(CIdiffs, function(CIframe){
        if(is.data.frame(CIframe)){
            sum(CIframe$diff_hill,na.rm = T)
        } else {0}
    }))}


    if(model == "gauss-gauss"){
    sig_level <- unlist(lapply(CIdiffs, function(CIframe){
            if(is.data.frame(CIframe)){
                sum(CIframe$diff_gauss,na.rm = T)
            } else { 0 }
        }))}

    sig_level
}
