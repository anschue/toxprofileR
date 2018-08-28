

#' Calculate logFC of public microarray data
#'
#' @param data EList
#' @param comparison data.frame describing the comparison for logFC
#'
#' @return dataframe of logFCs
#' @export
calc_logfc_public <- function(data, comparison){


    treatment <- as.character(data$targets$Comp_SubstanceName_Trivial) ==  as.character(comparison$Comp_SubstanceName_Trivial.treatment) &
      (data$targets$Exp_Conc                                 ==  comparison$Exp_Conc.treatment|is.na(comparison$Exp_Conc.treatment)) &
      (data$targets$Exp_Messzeit_hpf                         ==  comparison$Exp_Messzeit_hpf|is.na(comparison$Exp_Messzeit_hpf)) &
      (data$targets$Exp_Expositionsstart_hpf                 ==  comparison$Exp_Expositionsstart_hpf|is.na(comparison$Exp_Expositionsstart_hpf)) &
      (data$targets$Exp_Expositionsstop_hpf                  ==  comparison$Exp_Expositionsstop_hpf|is.na(data$targets$Exp_Expositionsstop_hpf))


    control <- as.character(data$targets$Comp_SubstanceName_Trivial) ==  as.character(comparison$Comp_SubstanceName_Trivial.control) &
      (data$targets$Exp_Conc                                 ==  comparison$Exp_Conc.control|is.na(comparison$Exp_Conc.control)) &
      (data$targets$Exp_Messzeit_hpf                         ==  comparison$Exp_Messzeit_hpf|is.na(comparison$Exp_Messzeit_hpf))

  if(sum(treatment)>1&sum(control)>1){

  ##take only probes for each gene with highest IQR
  IQRs <- apply(data$E, 1, IQR, na.rm=T)
  data <- data[order(IQRs, decreasing=T),]
  data <- data[!duplicated(data$genes$ensembl_gene_id),]
  data <- data[!is.na(data$genes$ensembl_gene_id),]

  mean_logFC <- apply(data$E, MARGIN = 1, function(probe){mean(probe[treatment])-mean(probe[control])})
  logFCframe <- data.frame(row.names = data$genes$ensembl_gene_id, logFC = mean_logFC)

  colnames(logFCframe) <- comparison$SampleName.treatment
  return(logFCframe)

  }else{
    logFCframe<-NA
    message("no replication - no FC")
    return(logFCframe)
  }

}
