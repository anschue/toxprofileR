load("./data-raw/Annotation/mapComplete.Rd")
new_annotation <- map

x <- load("./data-raw/publicdata/logFCframe_all.Rd")

load("./data-raw/toxuniverse/universe_nodeframe.Rd")

devtools::use_data(new_annotation, logFCframe_all, nodeframe, comparisons_merge, internal = T, overwrite = T)
