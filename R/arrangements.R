#' @name arrangements-package
#' @docType package
#' @importFrom methods new
#' @importFrom R6 R6Class
#' @import gmp
#' @useDynLib arrangements, .registration = TRUE, .fixes = "C_"
"_PACKAGE"

Arrangements <- R6::R6Class(
    c("Arrangements", "iter", "abstractiter"),
    private = list(
        state = NULL
    ),
    public = list(
        nextElem = function() {
            out <- self$getnext()
            is.null(out) && stop("StopIteration", call. = FALSE)
            out
        }
    )
)


validate_n_value <- function(x = NULL, k = n, n = NULL, v = NULL, freq = NULL, replace = FALSE) {
    .Call(C_validate_n_value, x, k, n, v, freq, replace)
}
