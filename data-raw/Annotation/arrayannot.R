library("biomaRt")

# load mapping from OakLabs----------------------------------------------------
map <- read.table("Danio_rerio.GRCz10.cdna.ncrna.map")
colnames(map) <- c("EnsemblID", "ProbeID", "allEnsemblIDs")

# remove version nr. from transcript ID----------------------------------------
map$EnsemblID <- unlist(lapply(strsplit(x = as.character(map$EnsemblID), split = ".", fixed = T), function(x) {
  x[[1]]
}))
map$allEnsemblIDs <- lapply(strsplit(as.character(map$allEnsemblIDs), split = ";", fixed = T), function(x) {
  unlist(lapply(strsplit(x, split = ".", fixed = T), function(y) {
    y[[1]]
  }))
})

map$countTranscriptID <- unlist(lapply(map$allEnsemblIDs, function(x) {
  length(x)
}))

# sum(map$countTranscriptID==1)
# length(unique(map$EnsemblID[map$countTranscriptID==1]))
# length(unique(map$EnsemblID))

# retrieve GeneIDs from Transcript IDs-----------------------------------------
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", version = "Ensembl Genes 83", host = "dec2015.archive.ensembl.org", dataset = "drerio_gene_ensembl")

map$gene_names <- pbapply::pblapply(X = map$allEnsemblIDs, function(x) {
  biomaRt::getBM(attributes = c("ensembl_transcript_id", "wikigene_name", "entrezgene", "wikigene_description", "ensembl_gene_id"), filters = "ensembl_transcript_id", values = x, mart = mart, uniqueRows = T)
})

save(map, file = "map.Rd")

load("map.Rd")

map$wikigene_name <- lapply(map$gene_names, function(x) {
  unique(x[, "wikigene_name"])
})
map$count_wikigene <- unlist(lapply(map$gene_names, function(x) {
  length(unique(x[, "wikigene_name"]))
}))

map$entrezgene <- lapply(map$gene_names, function(x) {
  unique(x[, "entrezgene"])
})
map$count_entrezgene <- unlist(lapply(map$gene_names, function(x) {
  length(unique(x[, "entrezgene"]))
}))

map$ensembl_gene_id <- lapply(map$gene_names, function(x) {
  unique(x[, "ensembl_gene_id"])
})
map$count_ensemblgene <- unlist(lapply(map$gene_names, function(x) {
  length(unique(x[, "ensembl_gene_id"]))
}))

map$wikigene_description <- lapply(map$gene_names, function(x) {
  unique(x[, "wikigene_description"])
})

map$unique <- apply(map, MARGIN = 1, FUN = function(x) {
  (x["count_wikigene"] == 1 & !is.na(x["wikigene_name"])) |
    (x["count_entrezgene"] == 1 & !is.na(x["entrezgene"])) |
    (x["count_ensemblgene"] == 1 & !is.na(x["ensembl_gene_id"]))
})

sum(!map$unique & is.na(map$wikigene_name) & is.na(map$entrezgene) & is.na(map$ensembl_gene_id))

head(map)

x <- map[map$unique & map$count_ensemblgene != 1, ]

map$wikigene_name[!map$unique][1:500]

length(unique(map$ensembl_gene_id[map$unique]))
length(unique(map$EnsemblID))

map[map$EnsemblID == "ENSDART00000164144", ]
library("Biostrings")
ncrna <- Biostrings::readDNAStringSet(filepath = "../RawData/ArrayAnnotation/Danio_rerio.GRCz10.ncrna.fa", format = "fasta")

ncEnsembl <- unlist(lapply(strsplit(names(ncrna), split = ".", fixed = T), function(x) {
  x[[1]]
}))

length(unique(map$ensembl_gene_id[map$EnsemblID %in% ncEnsembl & map$unique]))
length(unique(map$EnsemblID[map$EnsemblID %in% ncEnsembl & map$countTranscriptID == 1]))

new_annotation <- map



save(map, file = "./data-raw/Annotation/mapComplete.Rd")

# firstannotation<-read.table("../RawData/ArrayAnnotation/Design_69507/AllAnnotations/069507_D_AA_20140902.txt",header=T,sep="\t",quote="",fill=T)
# length(unique(firstannotation$PrimaryAccession))

# load(file="./rawdata/Annotation/mapComplete.Rd")
