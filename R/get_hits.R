#' Get Exon Hits
#'
#'
#' @param exon_bed path to .bed file
#' @param mart biomart object
#' @param maxmismatch maximum number of mismatches
#' @param probelengths length of probes
#' @param name database name
#'
#' @return a data.frame containing probe-exon associations
#'
#' @export
#'
get_hits_exons <-
    function(exon_bed,
             mart,
             maxmismatch,
             probelengths,
             name = "exon",
             arraytype = c("regular","Affy_ST")) {

        # read bed --------------------------------------------------------------------
        message("read bed")
        drexon <- read.table(exon_bed, header = F)
        colnames(drexon) <-
            c(
                "seqnames",
                "start",
                "end",
                "ProbeID",
                "score",
                "strand",
                "thickStart",
                "thickEnd",
                "itemRgb",
                "blockCount",
                "blockSizes",
                "blockStarts"
            )

        # remove version number from transcript id --------------------------------
        if (sum(grepl(pattern = ".", drexon$seqnames, fixed = T)) > 0) {
            message("remove version number from transcript id")
            drexon$seqnames <-
                unlist(lapply(strsplit(
                    x = as.character(drexon$seqnames),
                    split = ".",
                    fixed = T
                ), function(x) {
                    x[1]
                }))
        }

        # remove duplicate rows ---------------------------------------------------
        drexon <- unique(drexon, MARGIN = 1)

        # remove hits with score lower than set minimum score ---------------------
        message("remove hits with score lower than set minimum score")
        drexon$length <-
            probelengths[match(drexon$ProbeID, names(probelengths))]
        drexon$nmismatch <- drexon$length - drexon$score
        drexon <- drexon[drexon$nmismatch <= maxmismatch,]

        # ###Transform ProbeIDs to valid Name (!!!Need to be done for Array-Data the same!!!)
        # drexon$ProbeID<-make.names(drexon$ProbeID)
        #

        # take only hits on "forward" strand (since there is no reverse strand for transcripts) --------
        if(arraytype!="Affy_ST"){
            drexon <- drexon[drexon$strand == "+",]
        }

        # for Affy-ST take hits from "reverse" strand###
        if(arraytype=="Affy_ST"){
            drexon <- drexon[drexon$strand == "-",]
        }

        # get GeneIDs for exon alignments -----------------------------------------
        message("retreive Ensembl GeneIDs for exon alignments")
        BM <-
            biomaRt::getBM(
                attributes = c("ensembl_gene_id", "ensembl_transcript_id"),
                filters = "ensembl_transcript_id",
                values = unique(drexon$seqnames),
                mart = mart,
                uniqueRows = T
            )
        #BM <- BM[!duplicated(BM[, "ensembl_transcript_id"]), ]

        annot_table <-
            merge(
                drexon,
                BM,
                by.x = "seqnames",
                by.y = "ensembl_transcript_id",
                all = T,
                sort = F
            )



        # aggregate table by ProbeID ----------------------------------------------
        message("aggregate table by ProbeID")
        aggr_table_exon <-
            aggregate(
                annot_table[, c("ensembl_gene_id", "score")],
                by = list(ProbeID = annot_table$ProbeID),
                FUN = function(x) {
                    list(paste(x))
                }
            )

        # determine number of hits per probe
        aggr_table_exon$n_exon <-
            sapply(
                aggr_table_exon$ensembl_gene_id,
                FUN = function(x) {
                    length(unique(x))
                }
            )

        aggr_table_exon$ensembl_gene_id_exon <-
            apply(aggr_table_exon, MARGIN = 1, function(probe) {
                probe$ensembl_gene_id_exon_all[which.max(probe$score)]
            })

        colnames(aggr_table_exon)[colnames(aggr_table_exon) == "ensembl_gene_id"] <-
            paste0("ensembl_gene_id_", name, "_all")
        colnames(aggr_table_exon)[colnames(aggr_table_exon) == "n_exon"] <-
            paste0("n_", name)
        colnames(aggr_table_exon)[colnames(aggr_table_exon) == "ensembl_gene_id_exon"] <-
            paste0("ensembl_gene_id_", name)

        # clean up Table --------------------------------------
        aggr_table_exon$ProbeID <- as.character(aggr_table_exon$ProbeID)

        return(aggr_table_exon)
    }

#' Get Genome hits
#'
#' @param genome_bed path to .bed file
#' @param mart mart object
#' @param maxmismatch maximum number of mismatches
#' @param probelengths length of probes
#'
#' @return a data.frame containing probe-genome associations
#'
#' @import data.table
#' @import GenomicRanges
#' @export
#'
get_hits_genome <-
    function(genome_bed,
             mart,
             maxmismatch,
             probelengths,
             arraytype = c("regular","Affy_ST")) {

        # read bed ----------------------------------------------------------------
        message("read .bed file")
        drGenome <-
            data.table::fread(genome_bed, header = F, sep = "\t")[, 1:12]
        #drGenome <- read.table(genome_bed, header = F,sep = "\t")[,1:12]
        colnames(drGenome) <-
            c(
                "chr",
                "start",
                "end",
                "name",
                "score",
                "strand",
                "thickStart",
                "thickEnd",
                "itemRgb",
                "blockCount",
                "blockSizes",
                "blockStarts"
            )

        # remove duplicate rows -------------------------------------------------------
        drGenome <- unique(drGenome, MARGIN = 1)

        # remove rows with scores lower than minimal score ----------------------------
        message("remove rows with scores lower than minimal score")
        drGenome$length <-
            probelengths[match(drGenome$name, names(probelengths))]
        drGenome$nmismatch <- drGenome$length - drGenome$score
        drGenome <- drGenome[drGenome$nmismatch <= maxmismatch,]

        # if there are several hits for one probe, only take the ones from localized scaffolds (if available)
        message(
            "if there are several hits for one probe, only take the ones from localized scaffolds (if available)"
        )

        drGenome$hitID <- c(1:nrow(drGenome))

        drGenome_byProbe <-
            drGenome[, .(chr = list(chr),
                         hitID = list(hitID),
                         n = .N), by = name][, hits_remove := mapply(function(chrs, ns, hitIDs) {
                             if (ns > 1 &
                                 sum(!grepl(pattern = "K", unlist(chrs))) > 0) {
                                 list(unlist(hitIDs)[grepl(pattern = "K", unlist(chrs))])
                             }
                         },
                         chrs = chr,
                         ns = n,
                         hitIDs = hitID)]

        hits_to_remove <- unlist(drGenome_byProbe$hits_remove)

        drGenome <- drGenome[!hitID %in% hits_to_remove, ]


        # # make valid Probenames -------------------------------------------------------
        # drGenome$name <- make.names(drGenome$name)
        #

        nHitsGenome <- drGenome[, .(n = .N), by = name]

        # replacing chromosome names
        message("replacing chromosome names")
        drGenome <-
            drGenome[, chr := gsub("chrX", "chr23", x = chr)][, chr := gsub("chrY", "chr24", x = chr)][, chr := gsub("chrM", "MT", x = chr)][, chr := gsub("chr", "", chr)][, chr := gsub("Un_KN", "KN", x = chr, fixed = T)][, chr := gsub("v", ".", x = chr, fixed = T)]

        # save information about number of probes not mapping -------------------------
        # annotInfo$nomap <-
        #     annotInfo$NrProbes - length(unique(c(drGenome$name, drGenome$ProbeID)))
        #

        # annotInfo$nonuniquemap_genome <-
        #     sum(nHitsGenome$n > 1)
        #
        # calculate new block starts to fit with chromosome counts ----------------
        message("calculate new block starts to fit with chromosome count")

        drGenome <-
            drGenome[, blockStartsNew := mapply(function(bl_starts, starts) {
                paste0(paste(
                    as.integer(strsplit(as.character(bl_starts), split = ",")[[1]]) + as.integer(starts),
                    collapse = ","
                ), ",")
            },
            bl_starts = blockStarts,
            starts = start)]


        # change strand-orientation for Affy-ST
        if(arraytype=="Affy_ST"){
            strandnew <- drGenome$strand
            strandnew[drGenome$strand=="+"] <- "-"
            strandnew[drGenome$strand=="-"] <- "+"
            drGenome$strand <- strandnew
        }

        # transform to GRangesList object -----------------------------------------
        message("transform to GRangesList object")

        drGenome_rangeslist <-
            with(
                drGenome,
                GenomicRanges::makeGRangesListFromFeatureFragments(
                    seqnames = chr,
                    fragmentStarts = blockStartsNew,
                    fragmentWidths = blockSizes,
                    strand = strand
                )
            )
        names(drGenome_rangeslist) <- drGenome$name

        # find overlaps between Alignment and Annotation ------------------------------
        message("find Overlaps")

        olap <-
            GenomicRanges::findOverlaps(drGenome_rangeslist, exbygene)

        olap_df <-
            data.table::data.table(
                ProbeID = names(drGenome_rangeslist)[queryHits(olap)],
                ensembl_gene_id_genome_all = names(exbygene)[subjectHits(olap)]
            )

        # determine overlap score
        x <- drGenome_rangeslist[queryHits(olap)]
        y <- exbygene[subjectHits(olap)]
        olap_df$overlap_length <-
            unlist(lapply(width(intersect(x, y)), sum))

        # create mapping data frame ----------------------------------------------------------------
        message("Create mapping data frame")


        # annotInfo$not_anotated <-
        #     length(unique(as.character(drGenome$name[!drGenome$name %in% unique(olap_df$ProbeID)])))
        #
        aggr_table_genome <-
            olap_df[, .(
                ensembl_gene_id_genome_all = list(ensembl_gene_id_genome_all),
                overlap_length = list(overlap_length),
                n_genome = uniqueN(ensembl_gene_id_genome_all)
            ), by = ProbeID][, ensembl_gene_id_genome := mapply(function(all_ids, olap_lengths) {
                unlist(all_ids)[which.max(unlist(olap_lengths))]
            }, all_ids = ensembl_gene_id_genome_all, olap_lengths = overlap_length)]


        # annotInfo$nonuniquemap_genome <-
        #     sum(aggr_table_genome$n_genome > 1)
        #
        aggr_table_genome <-
            merge(
                aggr_table_genome,
                nHitsGenome,
                by.x = "ProbeID",
                by.y = "name",
                all = T
            )

        return(aggr_table_genome)
    }



