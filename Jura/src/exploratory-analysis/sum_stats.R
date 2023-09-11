#' Calculate summary statistics.
#'
#' @param x data frame of numeric values.
#'
#' @return \code{x} data frame of summary statistics for each variable.
#'
#' @keywords sum_stats
#' @export
sum_stats <- function(x, method=c("full","simple")) {
  is_char <- is.character(x)
  if (is_char) {
    bd <- sum(grepl("[<]", x))
    x <- readr::parse_number(x)
  } else {
    bd <- NA
  }
  y <- na.omit(x)
  n <- max(digits(x), 2)
  if (any(!is.na(x))) {
    v <- c(Number = length(x),
           Missing = sum(is.na(x)),
           `Below Detection` = bd,
           Completeness = round(length(y)/length(x)*100, 1),
           Unique = length(unique(y)),
           Min = round(min(y), n),
           Max = round(max(y), n),
           Range = round(diff(range(y)), n),
           Mean = round(mean(y), n),
           Median = round(median(y), n),
           Var = round(var(y), 3),
           CV = signif(sd(y) / mean(y), 3))
  } else {
    v <- c(Number = 0,
           Missing = sum(is.na(x)),
           `Below Detection` = bd,
           Completeness = 0,
           Unique = NA,
           Min = NA,
           Max = NA,
           Range = NA,
           Mean = NA,
           Median = NA,
           Var = NA,
           CV = NA)
  }

  method <- match.arg(method, c("full","simple"))
  if (method == "simple") {
    v <- v[-c(3,4,5,12)]
  } else {
    if (!is_char) {
      v <- v[-3]
    }
  }

  return(v)
}
