load("./data-raw/Annotation/mapComplete.Rd")
new_annotation <- map
x <- load("./data-raw/publicdata/logFCframe_all.Rd")

devtools::use_data(new_annotation, logFCframe_all, comparisons_merge, internal = T, overwrite = T)
