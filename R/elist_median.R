

#' Calculate Probe Medians of ExpressionSet
#'
#' @description This function takes an EList, calculates a score from the quality flags
#' and removes measurements with lower scores than the minimum score for that probe. If
#' there are multiple measurements for the same probe the median is calculated.
#'
#'
#' @param data an EList
#'
#' @return an Elist
#'
#' @import limma
#' @export
#'
#' @examples No example yet
#'
#' @author This function is based on a function of Kristin Reiche

elist.median <- function(data) {
  # Initialize dataframe, each column must be of type numeric
  exprs.df <- as.data.frame(data$E)
  for (i in seq(1, nrow(data$targets))) {
    exprs.df[, i] <- as.numeric(as.vector(exprs.df[, i]))
  }


  ## calculate scores for quality flags
  scores <-
    data$isNonUniform + data$isNonUniform + data$isPopOutlier + data$isPopOutlier +
    data$isNonUniformBG

  minscores <-
    aggregate(scores, by = list(data$genes$ProbeName), min) ## for each probename calculate minimum score

  minscores_all <-
    minscores[match(x = data$genes$ProbeName, table = minscores$Group.1), -1]

  scores <- as.data.frame(scores)
  flag <- minscores_all < scores ### flag those probes that have a higher score than the minimum of that probe-type

  exprs.df[flag] <- NA
  apply(
    exprs.df,
    MARGIN = 2,
    FUN = function(x) {
      sum(is.na(x))
    }
  )

  # Add spot IDs as additional column to dataframe
  exprs.df <- cbind(exprs.df, data$genes$ProbeName)

  # Set column names of dataframe
  colnames(exprs.df)[length(colnames(exprs.df))] <- "Name"

  # Calculate median for all probes that are spotted more than once on array
  exprs.df.ag <-
    aggregate(exprs.df[, 1:length(colnames(exprs.df)) - 1],
      by = list(exprs.df$Name),
      median,
      na.rm = T
    )

  apply(
    exprs.df.ag[, -1],
    MARGIN = 2,
    FUN = function(x) {
      sum(is.na(x))
    }
  )

  exprs <- as.matrix(exprs.df.ag[, -1])
  rownames(exprs) <- exprs.df.ag[, 1]

  # Construct Elist object
  data.new <- methods::new(
    Class = "EList",
    list(
      E = exprs,
      minscores = minscores[, -1],
      targets = data$targets,
      genes = data$genes[match(rownames(exprs), data$genes$ProbeName), ]
    )
  )


  return(data.new)
}
