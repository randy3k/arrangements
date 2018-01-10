#' @details
#' \preformatted{
#' getnext(d = 1L, type = NULL, drop = d == 1L && is.null(type))
#' collect(type = "r")
#' reset()
#' }
#' \describe{
#' \item{d}{number of fetched arrangements}
#'
#' \item{type}{if "r", "c" or "l" is specified, the returned value would be a
#'  "row-major" matrix, a "column-major" matrix or a list respectively}
#'
#' \item{drop}{vectorize a matrix or unlist a list}
#' }
