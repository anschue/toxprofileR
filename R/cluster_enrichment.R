#' Functional enrichment of universe clusters
#'
#' @param cluster_table a dataframe with the columns 'ensembl' and 'clustername'
#'
#' @return returns a list with functional enrichments for ZFIN, interpro, GO, and reactome annotations
#' @export
#'
enrich_clusters <- function(cluster_table) {
  if (requireNamespace("clusterProfiler", quietly = TRUE)) {
    ensembl_universe <- as.character(unique(cluster_table$ensembl))

    cluster_enrichments <- pbapply::pblapply(unique(cluster_table$clustername), function(cluster) {
      genes_ensembl <- as.character(unique(cluster_table$ensembl[cluster_table$clustername == cluster]))

      ezfin <- clusterProfiler::enricher(
        gene = genes_ensembl,
        universe = ensembl_universe,
        TERM2GENE = enrichment_terms$ZFIN$term2gene,
        TERM2NAME = enrichment_terms$ZFIN$term2name,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        minGSSize = 10,
        maxGSSize = 500
      )

      einterpro <- clusterProfiler::enricher(
        gene = genes_ensembl,
        universe = ensembl_universe,
        TERM2GENE = enrichment_terms$interpro$term2gene,
        TERM2NAME = enrichment_terms$interpro$term2name,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        minGSSize = 10,
        maxGSSize = 500
      )

      ereactome <- clusterProfiler::enricher(
        gene = genes_ensembl,
        universe = ensembl_universe,
        TERM2GENE = enrichment_terms$reactome$term2gene,
        TERM2NAME = enrichment_terms$reactome$term2name,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        minGSSize = 10,
        maxGSSize = 500
      )

      ego <- clusterProfiler::enrichGO(
        gene = genes_ensembl,
        universe = ensembl_universe,
        OrgDb = "org.Dr.eg.db",
        keyType = "ENSEMBL",
        ont = "all",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.01,
        qvalueCutoff = 0.05,
        readable = T,
        pool = T
      )

      results <- list(clustername = cluster, zfin = ezfin, interpro = einterpro, reactome = ereactome, GO = ego)
      return(results)
    })

    ## ZFIN
    zfin_enrichments <- lapply(cluster_enrichments, function(cluster_enrichment) {
      if (!is.null(cluster_enrichment$zfin)) {
        results <- cluster_enrichment$zfin@result
        if (dim(results)[1] > 0) {
          results$genenr_annot <- unlist(lapply(strsplit(results$GeneRatio, split = "/", fixed = T), function(x) {
            as.numeric(x[2])
          }))
          results$genenr_total <- length(cluster_enrichment$zfin@gene)
          results$clustername <- cluster_enrichment$clustername
          return(results)
        } else {
          return(NULL)
        }
      } else {
        return(NULL)
      }
    })

    zfin_enrichments <- do.call("rbind", zfin_enrichments)

    zfin_enrichments$ratio_annot <- zfin_enrichments$Count / zfin_enrichments$genenr_annot
    zfin_enrichments$ratio_total <- zfin_enrichments$Count / zfin_enrichments$genenr_total

    ## GO
    GO_enrichments <- lapply(cluster_enrichments, function(cluster_enrichment) {
      if (!is.null(cluster_enrichment$GO)) {
        results <- cluster_enrichment$GO@result
        if (dim(results)[1] > 0) {
          results$genenr_annot <- unlist(lapply(strsplit(results$GeneRatio, split = "/", fixed = T), function(x) {
            as.numeric(x[2])
          }))
          results$genenr_total <- length(cluster_enrichment$GO@gene)
          results$clustername <- cluster_enrichment$clustername
          return(results)
        } else {
          return(NULL)
        }
      } else {
        return(NULL)
      }
    })

    GO_enrichments <- do.call("rbind", GO_enrichments)
    GO_enrichments$ratio_annot <- GO_enrichments$Count / GO_enrichments$genenr_annot
    GO_enrichments$ratio_total <- GO_enrichments$Count / GO_enrichments$genenr_total

    ## Interpro
    interpro_enrichments <- lapply(cluster_enrichments, function(cluster_enrichment) {
      if (!is.null(cluster_enrichment$interpro)) {
        results <- cluster_enrichment$interpro@result
        if (dim(results)[1] > 0) {
          results$genenr_annot <- unlist(lapply(strsplit(results$GeneRatio, split = "/", fixed = T), function(x) {
            as.numeric(x[2])
          }))
          results$genenr_total <- length(cluster_enrichment$interpro@gene)
          results$clustername <- cluster_enrichment$clustername
          return(results)
        } else {
          return(NULL)
        }
      } else {
        return(NULL)
      }
    })

    interpro_enrichments <- do.call("rbind", interpro_enrichments)

    interpro_enrichments$ratio_annot <- interpro_enrichments$Count / interpro_enrichments$genenr_annot
    interpro_enrichments$ratio_total <- interpro_enrichments$Count / interpro_enrichments$genenr_total


    ## reactome
    reactome_enrichments <- lapply(cluster_enrichments, function(cluster_enrichment) {
      if (!is.null(cluster_enrichment$reactome)) {
        results <- cluster_enrichment$reactome@result
        if (dim(results)[1] > 0) {
          results$genenr_annot <- unlist(lapply(strsplit(results$GeneRatio, split = "/", fixed = T), function(x) {
            as.numeric(x[2])
          }))
          results$genenr_total <- length(cluster_enrichment$reactome@gene)
          results$clustername <- cluster_enrichment$clustername
          return(results)
        } else {
          return(NULL)
        }
      } else {
        return(NULL)
      }
    })

    reactome_enrichments <- do.call("rbind", reactome_enrichments)

    reactome_enrichments$ratio_annot <- reactome_enrichments$Count / reactome_enrichments$genenr_annot
    reactome_enrichments$ratio_total <- reactome_enrichments$Count / reactome_enrichments$genenr_total

    return(list(ZFIN = zfin_enrichments, interpro = interpro_enrichments, reactome = reactome_enrichments, GO = GO_enrichments))
  } else {
    message("For this function you need to install the package clusterProfiler")
  }
}
