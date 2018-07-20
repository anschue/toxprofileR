#' Title
#'
#' @param dslist
#' @param method
#' @param output
#'
#' @return
#' @export
#'
#' @examples
normalizeBatch <- function(dslist, method = "cloess", output = T){

  # check if probe names are in the same order---------------------------------
   lapply(1:(length(dslist)-1),function(id){
    if(!identical(dslist[[id]]$genes$ProbeName,dslist[[id+1]]$genes$ProbeName)){stop("Probes not in the same order")}
   })

  common.dataframe<-do.call("cbind", lapply(dslist,function(ds){ds[["E"]]}))

  if(output){boxplot(log2(common.dataframe))}

  ##Normalize common dataframe
  common.cloess<-limma::normalizeBetweenArrays(log2(common.dataframe),method="cyclicloess")

  ##Write normalized data into single dataframes
  dslist_norm <- lapply(1:length(dslist), function(id){
    start_id <- if(id == 1){1}else{1+sum(unlist(lapply(1:(id-1),function(x){ncol(dslist[[x]][["E"]])})))}
    end_id <- sum(unlist(lapply(1:id,function(x){ncol(dslist[[x]][["E"]])})))
    ds <- dslist[[id]]
    ds[["E"]] <- common.cloess[,start_id:end_id]
    return(ds)
  })

  return(dslist_norm)
}
