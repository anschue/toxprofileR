#' Normalize Batch
#'
#' @param dslist list of EListRaw objects which should be normalized together
#' @param method string, normalization method (default: "cyclicloess")
#' @param output logical, if output plots should be given
#'
#' @return a list of EList objects containing the normalized data
#' @export
#'
normalizeBatch <- function(dslist, method = "cyclicloess", output = T){

  # check if probe names are in the same order---------------------------------
   lapply(1:(length(dslist)-1),function(id){
    if(!identical(dslist[[id]]$genes$ProbeName,dslist[[id+1]]$genes$ProbeName)){stop("Probes not in the same order")}
   })

  common.dataframe<-do.call("cbind", lapply(dslist,function(ds){ds[["E"]]}))

  if(output){boxplot(log2(common.dataframe))}

  ##Normalize common dataframe
  common.cloess<-limma::normalizeBetweenArrays(log2(common.dataframe),method=method)

  if(output){boxplot(common.cloess)}

  ##Write normalized data into single dataframes
  dslist_norm <- lapply(1:length(dslist), function(id){
    start_id <- if(id == 1){1}else{1 + sum(unlist(lapply(1:(id-1), function(x){ncol(dslist[[x]][["E"]])})))}
    end_id <- sum(unlist(lapply(1:id,function(x){ncol(dslist[[x]][["E"]])})))
    ds <- new("EList", lapply(dslist[[id]],function(x){return(x)}))
    ds[["E"]] <- common.cloess[,start_id:end_id]
    return(ds)
  })

  return(dslist_norm)
}

#test_norm<-normalizeBatch(list(test1,test2))
