rm(list = ls())
# new array annotation --------------------------------------------------------
load("./data-raw/Annotation/mapComplete.Rd")
new_annotation <- map

# public toxicogenomic data ---------------------------------------------------
x <- load("./data-raw/publicdata/logFC_frame.Rd")

# toxicogenomic universe assignment -------------------------------------------
load("./data-raw/toxuniverse/universe_nodeframe.Rd")

# toxicogenomic universe grid -------------------------------------------------
load("./data-raw/toxuniverse/universe_grid.Rd")

# add to internal data --------------------------------------------------------
devtools::use_data(new_annotation, logFC_frame, nodeframe, comparisons_merge, grid, internal = T, overwrite = T)
