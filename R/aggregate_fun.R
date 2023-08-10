#' @export
aggregate_fun <- function(index) {
  if (is.numeric(index) == TRUE) {
    sum(index)
  }
  else {
    length(unique(index))
  }
}
