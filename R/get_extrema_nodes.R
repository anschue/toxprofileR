#' Get extreme values across several experiments per toxnode
#'
#' @param dslist A list of ELists with normalized logFC values
#'
#' @return A dataframe with minimum, maximum and absolute extreme value for each node
#' @export
#'
get_extrema_nodes <- function(dslist){

    library("pbapply")
    library("outliers")

    dslist_extrema <- pbapply::pblapply(X = seq(1,length(dslist)), FUN = function(substanceID){

    nodelist <- toxprofileR::create_nodelist(dslist[[substanceID]])

    # determine maximum for each gene (median of condition) -------------------

    logFCagg<-lapply(X = nodelist,FUN = function(node){
        if(is.data.frame(node)){

            # remove outliers -------------------------------------------------
            conc<-node$concentration_umol_l
            time<-node$time_hpe_factor
            timen<-node$time_hpe
            logFC<-node$logFC

            outliernew<-as.numeric(outliers::grubbs.test(x = as.numeric(logFC))["p.value"])

            while(log10(outliernew)< -3){
                time<-time[-which(logFC==outliers::outlier(logFC))]
                timen<-timen[-which(logFC==outliers::outlier(logFC))]
                conc<-conc[-which(logFC==outliers::outlier(logFC))]
                logFC<-logFC[-which(logFC==outliers::outlier(logFC))]
                outliernew<-as.numeric(outliers::grubbs.test(x = logFC)["p.value"])
            }

            # calculate median of each condition and get extrema --------------
            df_agg<-aggregate.data.frame(x=logFC,by = list(conc=conc,time=time),FUN = median)
            maxnode<-max(df_agg$x)
            minnode<-min(df_agg$x)
            extnode<-maxnode
            extnode[abs(minnode)>maxnode]<-minnode[abs(minnode)>maxnode]
            result<-c(minnode,maxnode,extnode)
            return(result)}else{return(c(NA,NA,NA))}
    })

    extnodes<-as.data.frame(do.call(rbind, logFCagg))
    colnames(extnodes) <- c("min", "max", "extremum")
    return(extnodes)
})

minnodes <- do.call("cbind",lapply(dslist_extrema, function(nodelist){nodelist$min}))
maxnodes <- do.call("cbind",lapply(dslist_extrema, function(nodelist){nodelist$max}))
extnodes <- do.call("cbind",lapply(dslist_extrema, function(nodelist){nodelist$ext}))

minnodes_all<-apply(minnodes,MARGIN = 1,FUN = min)
maxnodes_all<-apply(maxnodes,MARGIN = 1,FUN = max)

extnodes_all<-maxnodes_all
extnodes_all[abs(minnodes_all)>maxnodes_all&!is.na(maxnodes_all)]<-minnodes_all[abs(minnodes_all)>maxnodes_all&!is.na(maxnodes_all)]

extrema <- data.frame(min = minnodes_all, max = maxnodes_all, ext = extnodes_all)
return(extrema)
}
