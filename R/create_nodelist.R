#' Create logFC Nodelist
#'
#' @param elist An limma EList object with normalized logFC values
#'
#' @return A List with logFC values for each node of the toxicogenomic universe
#' @export
create_nodelist <- function(elist){

nodeframe_elist <- nodeframe[!is.na(nodeframe$ProbeID),]


    # create list with logFC and metadata for each node -----------------------
    nodelist <- lapply(X = seq(1,max(nodeframe_elist$toxnode)),FUN = function(nodeID){
        if(sum(nodeframe_elist$toxnode==nodeID)>0){
            logFC<-c(t(elist$E[elist$genes$ProbeName%in%nodeframe_elist$ProbeID[nodeframe_elist$toxnode==nodeID],elist$targets$type!="recovery"]))
            conc_all<-rep(elist$targets$concentration_umol_l[elist$targets$type!="recovery"],times=sum(nodeframe_elist$toxnode==nodeID))
            time_all<-ordered(rep(elist$targets$time_hpe[elist$targets$type!="recovery"],times=sum(nodeframe_elist$toxnode==nodeID)))
            timen_all<-rep(elist$targets$time_hpe[elist$targets$type!="recovery"],times=sum(nodeframe_elist$toxnode==nodeID))
            probe_id <- rep(as.character(nodeframe_elist$ProbeID)[nodeframe_elist$toxnode==nodeID], each = nrow(t(elist$E[elist$genes$ProbeName%in%nodeframe_elist$ProbeID[nodeframe_elist$toxnode==nodeID],elist$targets$type!="recovery"])))
            ensembl_gene_id <- rep(as.character(nodeframe_elist$ensembl)[nodeframe_elist$toxnode==nodeID], each = nrow(t(elist$E[elist$genes$ProbeName%in%nodeframe_elist$ProbeID[nodeframe_elist$toxnode==nodeID],elist$targets$type!="recovery"])))
            probeframe<-data.frame(logFC=logFC,conc_all=conc_all,time_all=time_all,timen_all=timen_all, probe_id = probe_id, ensembl_gene_id = ensembl_gene_id, nodeID = nodeID)
            return(probeframe)
        }else{return(NA)}
    })

    nodelist
}
