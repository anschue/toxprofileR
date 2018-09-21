rm(list = ls())
# new array annotation --------------------------------------------------------
new_annotation <- readRDS("./data-raw/Annotation/GPL20686annotation.Rds")

# # public toxicogenomic data ---------------------------------------------------
# x <- load("./data-raw/publicdata/logFC_frame.Rd")

# # toxicogenomic universe assignment -------------------------------------------
# load("./data-raw/toxuniverse/universe_nodeframe.Rd")
#
# # toxicogenomic universe grid -------------------------------------------------
# load("./data-raw/toxuniverse/universe_grid.Rd")

# add to internal data --------------------------------------------------------
devtools::use_data(
    #grid,
    #logFC_frame,
    new_annotation,
    #nodeframe,
    internal = T,
    overwrite = T
)
