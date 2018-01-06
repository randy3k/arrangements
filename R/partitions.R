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
        m = NULL,
        descending = NULL,
        initialize = function(n, m=NULL, descending = FALSE) {
            self$n <- n
            self$m <- m
            self$descending <- descending
            self$reset()
        },
        reset = function() {
            private$state <- new.env()
            private$null_pending <- FALSE
        },
        collect = function(type = 'r') {
            P <- tryCatch(npartitions(self$n, self$m), error = function(e) {
                if (startsWith(e$message, "integer overflow")) {
                    stop("too many results")
                } else {
                    stop(e)
                }
            })
            out <- self$getnext(P, type, drop = FALSE)
            self$reset()
            out
        },
        getnext = function(d = 1L, type = 'r', drop = d == 1L) {
            if (private$null_pending) {
                out <- NULL
            } else {
                out <- next_partitions(
                    self$n, self$m, d, private$state, self$descending, type)
            }
            if (is.null(out) || length(out) == 0) {
                self$reset()
            } else if (d > 1) {
                if (type == 'r' && nrow(out) < d){
                    private$null_pending <- TRUE
                } else if (type == 'c' && ncol(out) < d){
                    private$null_pending <- TRUE
                } else if (type == 'l' && length(out) < d){
                    private$null_pending <- TRUE
                }
            }
            if (!is.null(out) && drop) {
                if (type == 'l') {
                    out <- out[[1]]
                } else {
                    dim(out) <- NULL
                }
            }
            out
        },
        print = function(...) {
            if (is.null(self$m)) {
                cat("Partitions of", self$n, "\n")
            } else {
                cat("Partitions of", self$n, "into", self$m, "parts\n")
            }
            invisible(self)
        }
    )
)

next_partitions <- function(n, m, d, state, descending, type) {
    if (d > 1) {
        if ((type != 'l' && is.null(m) && d * n > .Machine$integer.max) ||
                (type != 'l' && !is.null(m) && d * m > .Machine$integer.max) ||
                (type == 'l' && d > .Machine$integer.max)) {
            stop("too many results")
        }
    }

    if (is.null(m)) {
        if (descending) {
            out <- .Call(
                "next_desc_partitions",
                PACKAGE = "arrangements",
                as.integer(n),
                as.integer(d),
                state,
                type)
        } else {
            out <- .Call(
                "next_asc_partitions",
                PACKAGE = "arrangements",
                as.integer(n),
                as.integer(d),
                state,
                type)
        }
    } else {
        if (descending) {
            out <- .Call(
                "next_desc_k_partitions",
                PACKAGE = "arrangements",
                as.integer(n),
                as.integer(m),
                as.integer(d),
                state,
                type)
        } else {
            out <- .Call(
                "next_asc_k_partitions",
                PACKAGE = "arrangements",
                as.integer(n),
                as.integer(m),
                as.integer(d),
                state,
                type)
        }
    }

    if (!is.null(out)) {
        if (type == 'r') {
            if (is.null(m)) {
                dim(out) <- c(length(out) / n, n)
            } else {
                dim(out) <- c(length(out) / m, m)
            }
        } else if (type == 'c') {
            if (is.null(m)) {
                dim(out) <- c(n, length(out) / n)
            } else {
                dim(out) <- c(m, length(out) / m)
            }
        }
    }
    out
}

#' @export
partitions <- function(n, m=NULL, descending = FALSE, type = 'r') {
    P <- tryCatch(npartitions(n, m), error = function(e) {
        if (startsWith(e$message, "integer overflow")) {
            stop("too many results")
        } else {
            stop(e)
        }
    })
    next_partitions(n, m, P, NULL, descending, type)
}


#' @export
ipartitions <- function(n, m=NULL, descending = FALSE) {
    Partitions$new(n, m, descending)
}


#' @export
npartitions <- function(n, m=NULL, bigz=FALSE) {
    (n > 0 && n %% 1 == 0) || stop("n should be a positive integer")
    if (is.null(m)) {
        if (n > 120L) {
            out <- gmp::as.bigz(.Call("npart_bigz", PACKAGE = "arrangements", as.integer(n)))
        } else {
            out <- .Call("npart", PACKAGE = "arrangements", as.integer(n))
        }
    } else {
        (m > 0 && m %% 1 == 0) || stop("m should be a positive integer")
        if (n < m) {
            out <- 0
        } else if (n > 158L) {
            out <- gmp::as.bigz(
                .Call("nfixedpart_bigz", PACKAGE = "arrangements", as.integer(n), as.integer(m)))
        } else {
            out <- .Call("nfixedpart", PACKAGE = "arrangements", as.integer(n), as.integer(m))
        }
    }
    convertz(out, bigz)
}
