rm(list = ls())

# toxicogenomic universe assignment -------------------------------------------
load("./data-raw/toxuniverse/universe_nodeframe.Rd")

# toxicogenomic universe grid -------------------------------------------------
load("./data-raw/toxuniverse/universe_grid.Rd")


# new array annotation --------------------------------------------------------
new_annotation <- readRDS("./data-raw/Annotation/GPL20686annotation.Rds")

## Annotation for Enrichment
mart <- biomaRt::useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "drerio_gene_ensembl",
  host = "jul2018.archive.ensembl.org"
)

# ZFIN annotation (downloaded from zfin.org on 2018/11/22)
ZFIN_terms <- read.table(file = "./data-raw/Annotation/wildtype-expression_fish_2018.11.20.txt", skip = 1, header = T, sep = "\t", quote = "")
ZFIN_IDS <- biomaRt::getBM(attributes = c("zfin_id_id", "ensembl_gene_id"), values = as.character(ZFIN_terms$Gene.ID), filters = "zfin_id_id", mart = mart, uniqueRows = T)
ZFIN_IDS <- unique(ZFIN_IDS, margin = 1)
ZFIN_terms <- merge.data.frame(ZFIN_terms, ZFIN_IDS, by.x = "Gene.ID", by.y = "zfin_id_id", all = T)

zfin_term2gene <- ZFIN_terms[, c("Super.Structure.ID", "ensembl_gene_id")]
zfin_term2gene <- zfin_term2gene[!is.na(zfin_term2gene$Super.Structure.ID) & !is.na(zfin_term2gene$ensembl_gene_id), ]
zfin_term2name <- ZFIN_terms[, c("Super.Structure.ID", "Super.Structure.Name")]
zfin_term2name <- unique(zfin_term2name[!is.na(zfin_term2name$Super.Structure.ID) & !is.na(zfin_term2name$Super.Structure.Name), ], margin = 1)

# Interpro annotation
interpro_ids <- biomaRt::getBM(attributes = c("interpro", "interpro_description", "ensembl_gene_id"), filters = "ensembl_gene_id", values = nodeframe$ensembl, mart = mart)
interpro_term2gene <- unique(interpro_ids[, c("interpro", "ensembl_gene_id")], margin = 1)
interpro_term2name <- unique(interpro_ids[, c("interpro", "interpro_description")], margin = 1)

# Reactome annotation
library("reactome.db")
reactome_ids <- biomaRt::getBM(attributes = c("reactome", "ensembl_gene_id"), filters = "ensembl_gene_id", values = nodeframe$ensembl, mart = mart)
reactome_term2gene <- unique(reactome_ids[, c("reactome", "ensembl_gene_id")], margin = 1)
reactome_term2name <- data.frame(reactome = unique(reactome_term2gene$reactome), name = NA, stringsAsFactors = F)
reactome_term2name$name <- unlist(lapply(reactome_term2name$reactome, function(r_id) {
  name <- reactomePATHID2NAME[[r_id]]
  if (!is.null(name)) {
    return(name)
  } else {
    return(NA)
  }
}))
reactome_term2name <- reactome_term2name[!is.na(reactome_term2name$name), ]

enrichment_terms <- list(
  ZFIN = list(term2gene = zfin_term2gene, term2name = zfin_term2name),
  interpro = list(term2gene = interpro_term2gene, term2name = interpro_term2name),
  reactome = list(term2gene = reactome_term2gene, term2name = reactome_term2name)
)


# add to internal data --------------------------------------------------------
devtools::use_data(
  new_annotation,
  enrichment_terms,
  internal = T,
  overwrite = T
)

# add to external data
devtools::use_data(
  grid,
  nodeframe,
  internal = F,
  overwrite = T
)
