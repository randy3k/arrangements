#' @details
#' The `Partitions` class can be initialized by using the convenient wrapper `ipartitions` or
#' \preformatted{
#' Partitions$new(n, k = NULL, descending = FALSE)
#' }
#' @template iterator_methods
#' @rdname ipartitions
#' @export
Partitions <- R6::R6Class(
    "Partitions",
    inherit = Arrangements,
    private = list(
        state = NULL,
        null_pending = FALSE
    ),
    public = list(
        n = NULL,
        k = NULL,
        descending = NULL,
        initialize = function(n, k = NULL, descending = FALSE) {
            (n %% 1 == 0  && n >= 0) || stop("expect non-negative integer")
            self$n <- as.integer(n)
            if (!is.null(k)) {
                (k %% 1 == 0 && k >= 0) || stop("expect non-negative integer")
                self$k <- as.integer(k)
            }
            self$descending <- descending
            self$reset()
        },
        reset = function() {
            private$state <- new.env()
            private$null_pending <- FALSE
        },
        collect = function(type = "r") {
            out <- self$getnext(-1L, type, drop = FALSE)
            self$reset()
            out
        },
        getnext = function(d = 1L, type = NULL, drop = d == 1L && is.null(type)) {
            if (private$null_pending) {
                out <- NULL
                self$reset()
            } else {
                out <- next_partitions(
                    self$n, self$k, d, private$state, self$descending, type)
                if (type == "r" || is.null(type)){
                    if (nrow(out) == 0) {
                        out <- NULL
                        self$reset()
                    } else if (nrow(out) < d || ncol(out) == 0) {
                        private$null_pending <- TRUE
                    }
                    if (!is.null(out) && drop) {
                        dim(out) <- NULL
                    }
                } else if (type == "c"){
                    if (ncol(out) == 0) {
                        out <- NULL
                        self$reset()
                    } else if (ncol(out) < d || nrow(out) == 0) {
                        private$null_pending <- TRUE
                    }
                    if (!is.null(out) && drop) {
                        dim(out) <- NULL
                    }
                } else if (type == "l"){
                    if (length(out) == 0) {
                        out <- list()
                        self$reset()
                    } else if (length(out) < d) {
                        private$null_pending <- TRUE
                    }
                    if (length(out) > 0 && drop) {
                        out <- unlist(out)
                    }
                }
            }
            out
        },
        print = function(...) {
            if (is.null(self$k)) {
                cat("Partitions of", self$n, "\n")
            } else {
                cat("Partitions of", self$n, "into", self$k, "parts\n")
            }
            invisible(self)
        }
    )
)

next_partitions <- function(n, k, d, state, descending, type) {
    if (is.null(k)) {
        if (n == 0) {
            if (type == "r" || is.null(type)) {
                out <- integer(0)
                dim(out) <- c(1, 0)
            } else if (type == "c") {
                out <- integer(0)
                dim(out) <- c(0, 1)
            } else {
                out <- list(integer(0))
            }
        } else if (descending) {
            out <- .Call(
                "next_desc_partitions",
                PACKAGE = "arrangements",
                n,
                d,
                state,
                type)
        } else {
            out <- .Call(
                "next_asc_partitions",
                PACKAGE = "arrangements",
                n,
                d,
                state,
                type)
        }
    } else {
        if (n < k) {
            if (type == "r" || is.null(type)) {
                out <- integer(0)
                dim(out) <- c(0, k)
            } else if (type == "c") {
                out <- integer(0)
                dim(out) <- c(k, 0)
            } else {
                out <- list()
            }
        } else if (n == 0 && k == 0) {
            if (type == "r" || is.null(type)) {
                out <- integer(0)
                dim(out) <- c(1, 0)
            } else if (type == "c") {
                out <- integer(0)
                dim(out) <- c(0, 1)
            } else {
                out <- list(integer(0))
            }
        } else if (k == 0) {
            if (type == "r" || is.null(type)) {
                out <- integer(0)
                dim(out) <- c(0, 0)
            } else if (type == "c") {
                out <- integer(0)
                dim(out) <- c(0, 0)
            } else {
                out <- list()
            }
        } else if (descending) {
            out <- .Call(
                "next_desc_k_partitions",
                PACKAGE = "arrangements",
                n,
                k,
                d,
                state,
                type)
        } else {
            out <- .Call(
                "next_asc_k_partitions",
                PACKAGE = "arrangements",
                n,
                k,
                d,
                state,
                type)
        }
    }
    out
}

#' Partitions generator
#'
#' This function partitions an non-negative interger into `k` parts or any part size.
#' The results are in lexicographical or reversed lexicographical order.
#'
#' @param n an non-negative integer to be partitioned
#' @param k number of parts
#' @param descending logical to use reversed lexicographical order
#' @template param_type
#' @seealso [ipartitions] for iterating partitions and [npartitions] to calculate number of partitions
#' @examples
#' # all partitions of 6
#' partitions(6)
#' # reversed lexicographical order
#' partitions(6, descending = TRUE)
#'
#' # fixed number of parts
#' partitions(10, 5)
#' # reversed lexicographical order
#' partitions(10, 5, descending = TRUE)
#'
#' # column major
#' partitions(6, type = "c")
#' partitions(6, 3, type = "c")
#'
#' # list output
#' partitions(6, type = "l")
#' partitions(6, 3, type = "l")
#'
#' # zero sized partitions
#' dim(partitions(0))
#' dim(partitions(5, 0))
#' dim(partitions(5, 6))
#' dim(partitions(0, 0))
#' dim(partitions(0, 1))
#'
#' @export
partitions <- function(n, k = NULL, descending = FALSE, type = "r") {
    next_partitions(n, k, -1L, NULL, descending, type)
}

#' Partitions iterator
#'
#' This function returns a [Partitions](Partitions-class.html) iterator which
#' allows users to fetch the next partition(s) via the `getnext()` method. All remaing
#' partitions of the iterator can be fetched via the `collect()` method.
#'
#' @param n an non-negative integer to be partitioned
#' @param k number of parts
#' @param descending logical to use reversed lexicographical order
#' @seealso [partitions] for generating all partitions and [npartitions] to calculate number of partitions
#' @examples
#' ipart <- ipartitions(10)
#' ipart$getnext()
#' ipart$getnext(2)
#' ipart$getnext(type = "c", drop = FALSE)
#' # collect remaining partitions
#' ipart$collect()
#'
#' library(foreach)
#' foreach(x = ipartitions(6, 2), .combine=c) %do% {
#'   prod(x)
#' }
#' @export
ipartitions <- function(n, k = NULL, descending = FALSE) {
    Partitions$new(n, k, descending)
}

#' Number of partitions
#'
#' @param n an non-negative integer to be partitioned
#' @param k number of parts
#' @param bigz an logical to indicate using [gmp::bigz]
#' @seealso [partitions] for generating all partitions and [ipartitions] for iterating partitions
#' @examples
#' # number of partitions of 10
#' npartitions(10)
#' # number of partitions of 10 into 5 parts
#' npartitions(10, 5)
#'
#' # integer overflow
#' \dontrun{npartitions(160)}
#' npartitions(160, bigz = TRUE)
#'
#' # zero sized partitions
#' npartitions(0)
#' npartitions(5, 0)
#' npartitions(5, 6)
#' npartitions(0, 0)
#' npartitions(0, 1)
#' @export
npartitions <- function(n, k = NULL, bigz = FALSE) {
    (n %% 1 == 0  && n >= 0) || stop("expect non-negative integer")
    if (is.null(k)) {
        if (bigz) {
            out <- .Call("npart_bigz", PACKAGE = "arrangements", n)
        } else {
            out <- .Call("npart", PACKAGE = "arrangements", n)
        }
    } else {
        (k %% 1 == 0  && k >= 0) || stop("expect non-negative integer")
        if (bigz) {
            out <- .Call("npart_k_bigz", PACKAGE = "arrangements", n, k)
        } else {
            out <- .Call("npart_k", PACKAGE = "arrangements", n, k)
        }
    }
    convertz(out, bigz)
}
