digits <- function(x) {
  d <- if (is.numeric(x)) as.character(x) else x
  d <- nchar(sub('.+\\.', '', x))
  d <- ifelse(grepl("\\.", x), d, 0)
  d <- round(mean(d, na.rm=TRUE), 0)

  return(d)
}
