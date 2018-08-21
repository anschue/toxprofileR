#' get biomaRt annotaiton
#'
#' @param values IDs
#' @param filter Type of IDs (default: "ensembl_gene_id")
#' @param mart a biomaRt object
#'
#' @return a dataframe with selected biomaRt annotations
#' @export
#'
getBM_annotation <-
    function(values = mapFinal_reduced$ensembl_gene_id,
             filter = "ensembl_gene_id",
             mart) {

        library("biomaRt")

            if (filter != "ensembl_gene_id"){
        BM_ensembl <-
            biomaRt::getBM(
                attributes = unique(
                    c(
                        filter,
                        "ensembl_gene_id"
                     )
                ),
                filters = filter,
                values = values,
                mart = mart,
                uniqueRows = T
            )
           
           
        values <- BM_ensembl$ensembl_gene_id
        }
        
        
        BM_IDS <-
            biomaRt::getBM(
                attributes =
                    c(
                        "ensembl_gene_id",
                        "external_gene_name",
                        "description",
                        "hsapiens_homolog_associated_gene_name"
                    ),
                filters = "ensembl_gene_id",
                values = values,
                mart = mart,
                uniqueRows = T
            )
        
        BM_IDS <-
            aggregate.data.frame(
                BM_IDS,
                by = list(BM_IDS[, "ensembl_gene_id"]),
                FUN = function(x) {
                    paste(unique(x))
                }
            )[, -1]
        
        BM_GO <-
            biomaRt::getBM(
                attributes = c("ensembl_gene_id", "name_1006"),
                filters = "ensembl_gene_id",
                values = values,
                mart = mart,
                uniqueRows = T
            )
        
        BM_GO <-
            aggregate.data.frame(
                BM_GO,
                by = list(BM_GO[, "ensembl_gene_id"]),
                FUN = function(x) {
                    paste(unique(x))
                }
            )[, -1]
        
        
        BM_interpro <-
            biomaRt::getBM(
                attributes = c("ensembl_gene_id", "interpro_description"),
                filters = "ensembl_gene_id",
                values = values,
                mart = mart,
                uniqueRows = T
            )
        
        BM_interpro <-
            aggregate.data.frame(
                BM_interpro,
                by = list(BM_interpro[, "ensembl_gene_id"]),
                FUN = function(x) {
                    paste(unique(x))
                }
            )[, -1]
        
        BM_phenotype <-
            biomaRt::getBM(
                attributes = c("ensembl_gene_id", "phenotype_description"),
                filters = "ensembl_gene_id",
                values = values,
                mart = mart,
                uniqueRows = T
            )
        
        BM_phenotype <-
            aggregate.data.frame(
                BM_phenotype,
                by = list(BM_phenotype[, "ensembl_gene_id"]),
                FUN = function(x) {
                    paste(unique(x))
                }
            )[, -1]
        
        BM_genetype <-
            biomaRt::getBM(
                attributes = c("ensembl_gene_id", "gene_biotype"),
                filters = "ensembl_gene_id",
                values = values,
                mart = mart,
                uniqueRows = T
            )
        
        BM_genetype <-
            aggregate.data.frame(
                BM_genetype,
                by = list(BM_genetype[, "ensembl_gene_id"]),
                FUN = function(x) {
                    paste(unique(x))
                }
            )[, -1]
        
        
        BM_all <-
            Reduce(
                function(...)
                    merge(..., by = "ensembl_gene_id", all = TRUE),
                list(BM_IDS, BM_GO, BM_interpro, BM_genetype, BM_phenotype)
            )
            
         if (filter != "ensembl_gene_id"){
			 BM_all <- merge(BM_ensembl, BM_all, by = "ensembl_gene_id", all = TRUE)
			 BM_all <- aggregate.data.frame(
			 BM_all,
			 by = list(BM_all[,filter]),
			  FUN = function(x) {
                    paste(unique(x))
                }
                )[, -1]
			 }
        
        return(BM_all)
    }
