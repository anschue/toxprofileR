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

        BM_IDS <-
            biomaRt::getBM(
                attributes = unique(
                    c(
                        filter,
                        "ensembl_gene_id",
                        "external_gene_name",
                        "description",
                        "hsapiens_homolog_associated_gene_name"
                    )
                ),
                filters = filter,
                values = values,
                mart = mart,
                uniqueRows = T
            )

        BM_IDS <-
            aggregate.data.frame(
                BM_IDS,
                by = list(BM_IDS[, filter]),
                FUN = function(x) {
                    paste(unique(x))
                }
            )[, -1]

        BM_GO <-
            biomaRt::getBM(
                attributes = c("ensembl_gene_id", "name_1006"),
                filters = filter,
                values = values,
                mart = mart,
                uniqueRows = T
            )

        BM_GO <-
            aggregate.data.frame(
                BM_GO,
                by = list(BM_GO[, filter]),
                FUN = function(x) {
                    paste(unique(x))
                }
            )[, -1]


        BM_interpro <-
            biomaRt::getBM(
                attributes = c("ensembl_gene_id", "interpro_description"),
                filters = filter,
                values = values,
                mart = mart,
                uniqueRows = T
            )

        BM_interpro <-
            aggregate.data.frame(
                BM_interpro,
                by = list(BM_interpro[, filter]),
                FUN = function(x) {
                    paste(unique(x))
                }
            )[, -1]

        BM_phenotype <-
            biomaRt::getBM(
                attributes = c("ensembl_gene_id", "phenotype_description"),
                filters = filter,
                values = values,
                mart = mart,
                uniqueRows = T
            )

        BM_phenotype <-
            aggregate.data.frame(
                BM_phenotype,
                by = list(BM_phenotype[, filter]),
                FUN = function(x) {
                    paste(unique(x))
                }
            )[, -1]

        BM_genetype <-
            biomaRt::getBM(
                attributes = c("ensembl_gene_id", "gene_biotype"),
                filters = filter,
                values = values,
                mart = mart,
                uniqueRows = T
            )

        BM_genetype <-
            aggregate.data.frame(
                BM_genetype,
                by = list(BM_genetype[, filter]),
                FUN = function(x) {
                    paste(unique(x))
                }
            )[, -1]


        BM_all <-
            Reduce(
                function(...)
                    merge(..., by = filter, all = TRUE),
                list(BM_IDS, BM_GO, BM_interpro, BM_genetype, BM_phenotype)
            )

        return(BM_all)
    }
