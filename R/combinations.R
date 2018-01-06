#' @export
Combinations <- R6::R6Class(
    "Combinations",
    inherit = Arrangements,
    private = list(
        state = NULL,
        null_pending = FALSE
    ),
    public = list(
        n = NULL,
        r = NULL,
        x = NULL,
        f = NULL,
        replace = NULL,
        initialize = function(n, r, x=NULL, f=NULL, replace = FALSE) {
            self$n <- n
            self$r <- r
            self$x <- x
            self$f <- f
            self$replace <- replace
            self$reset()
        },
        reset = function() {
            private$state <- new.env()
            private$null_pending <- FALSE
        },
        collect = function(type = 'r') {
            P <- try(ncombinations(self$n, self$r, self$f, self$replace), silent = TRUE)
            if (inherits(P, "try-error")) stop("too many results, use `icombinations`")
            out <- self$getnext(P, type, drop = FALSE)
            self$reset()
            out
        },
        getnext = function(d = 1L, type = 'r', drop = d == 1L) {
            if (private$null_pending) {
                out <- NULL
            } else {
                out <- next_combinations(
                    self$n, self$r, d, private$state, self$x, self$f, self$replace, type)
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
            cat("Combinations of", self$r, " items from", self$n, "items\n")
            invisible(self)
        }
    )
)

next_combinations <- function(n, r, d, state, x, f, replace, type) {
    if (d > 1) {
        if ((type != 'l' && is.null(r) && d * r > .Machine$integer.max) ||
                (type == 'l' && d > .Machine$integer.max)) {
            stop("too many results, use `icombinations`")
        }
    }
    if (!is.null(f)) {
        f <- as.integer(f)
    }

    if (replace) {
        out <- .Call(
            "next_replace_combinations",
            PACKAGE = "arrangements",
            as.integer(n),
            as.integer(r),
            as.integer(d),
            state,
            x,
            type)
    } else if (is.null(f)) {
        out <- .Call(
            "next_combinations",
            PACKAGE = "arrangements",
            as.integer(n),
            as.integer(r),
            as.integer(d),
            state,
            x,
            type)
    } else {
        out <- .Call(
            "next_multiset_combinations",
            PACKAGE = "arrangements",
            as.integer(n),
            as.integer(r),
            as.integer(d),
            state,
            x,
            f,
            type)
    }

    if (!is.null(out)) {
        if (type == 'r') {
            dim(out) <- c(length(out) / r, r)
        } else if (type == 'c') {
            dim(out) <- c(r, length(out) / r)
        }
    }
    out
}

#' @export
combinations <- function(n, r, x=NULL, f=NULL, replace=FALSE, type = 'r') {
    n <- check_nrxf(n, r, x, f, replace)
    P <- try(ncombinations(n, r, f, replace), silent = TRUE)
    if (inherits(P, "try-error")) stop("too many results, use `icombinations`")
    next_combinations(n, r, P, NULL, x, f, replace, type)
}


#' @export
icombinations <- function(n, r, x=NULL, f=NULL, replace = FALSE) {
    n <- check_nrxf(n, r, x, f, replace)
    Combinations$new(n, r, x, f, replace)
}

#' @export
ncombinations <- function(n, r, f=NULL, replace=FALSE, bigz=FALSE) {
    if (missing(n) && !is.null(f)) {
        n <- sum(f)
    }
    if (bigz) {
        if (replace) {
            out <- gmp::chooseZ(n + r - 1 , r)
        } else if (is.null(f)) {
            out <- gmp::chooseZ(n, r)
        } else {
            out <- .Call("ncomb_f_bigz", PACKAGE = "arrangements", as.integer(f), as.integer(r))
        }

    } else {
        if (replace) {
            out <- choose(n + r - 1, r)
        } else if (is.null(f)) {
            out <- choose(n, r)
        } else {
            out <- .Call("ncomb_f", PACKAGE = "arrangements", as.integer(f), as.integer(r))
        }
    }
    convertz(out, bigz)
}
