#' @section Methods:
#' \preformatted{
#' getnext(d = 1L, type = "r", drop = d == 1L)
#' collect(type = "r")
#' reset()
#' }
#' \describe{
#' \item{\code{d}}{number of fetched arrangements}
#' \item{\code{type}}{if "r", "c" or "l" is specified, the returned value would be a
#'  "row-major" matrix, a "column-major" matrix or a list respectively}
#' \item{\code{drop}}{vectorize a matrix or unlist a list}
#' }
