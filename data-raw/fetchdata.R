rm(list = ls())
# new array annotation --------------------------------------------------------
new_annotation <- readRDS("./data-raw/Annotation/GPL20686annotation.Rds")

# toxicogenomic universe assignment -------------------------------------------
load("./data-raw/toxuniverse/universe_nodeframe.Rd")

# toxicogenomic universe grid -------------------------------------------------
load("./data-raw/toxuniverse/universe_grid.Rd")

# add to internal data --------------------------------------------------------
devtools::use_data(
    new_annotation,
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
