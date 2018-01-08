#' Partitions class
#'
#' \preformatted{
#' Partitions$new(n, k=NULL, descending = FALSE)
#' }
#' @param n integer to be partitioned
#' @param k integer: number of partitions
#' @param descending logical: lexicographical or reverse lexicographical order
#' @name Partitions-class
NULL

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
        initialize = function(n, k=NULL, descending = FALSE) {
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
        getnext = function(d = 1L, type = "r", drop = d == 1L) {
            if (private$null_pending) {
                out <- NULL
                self$reset()
            } else {
                out <- next_partitions(
                    self$n, self$k, d, private$state, self$descending, type)
                if (type == "r"){
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
        if (descending) {
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
            if (type == "r") {
                out <- integer(0)
                dim(out) <- c(0, k)
            } else if (type == "c") {
                out <- integer(0)
                dim(out) <- c(k, 0)
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
#' @export
partitions <- function(n, k=NULL, descending = FALSE, type = "r") {
    next_partitions(n, k, -1L, NULL, descending, type)
}

#' Partitions iterator
#' @export
ipartitions <- function(n, k=NULL, descending = FALSE) {
    Partitions$new(n, k, descending)
}

#' Number of partitions
#' @export
npartitions <- function(n, k=NULL, bigz=FALSE) {
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
