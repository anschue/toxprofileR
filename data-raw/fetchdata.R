# new array annotation --------------------------------------------------------
load("./data-raw/Annotation/mapComplete.Rd")
new_annotation <- map

# public toxicogenomic data ---------------------------------------------------
x <- load("./data-raw/publicdata/logFCframe_all.Rd")

# toxicogenomic universe assignment -------------------------------------------
load("./data-raw/toxuniverse/universe_nodeframe.Rd")

# add to internal data --------------------------------------------------------
devtools::use_data(new_annotation, logFCframe_all, nodeframe, comparisons_merge, internal = T, overwrite = T)
