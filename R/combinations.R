#" @export
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
            (n %% 1 == 0  && n >= 0) || stop("expect non-negative integer")
            (r %% 1 == 0  && r >= 0) || stop("expect non-negative integer")
            self$n <- as.integer(n)
            self$r <- as.integer(r)
            self$x <- x
            self$f <- as_uint_array(f)
            self$replace <- replace
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
                out <- next_combinations(
                    self$n, self$r, d, private$state, self$x, self$f, self$replace, type)
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
                        out <- NULL
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
            cat("Combinations of", self$r, " items from", self$n, "items\n")
            invisible(self)
        }
    )
)

next_combinations <- function(n, r, d, state, x, f, replace, type) {
    if (r == 0) {
        if (type == "r") {
            out <- integer(0)
            dim(out) <- c(1, 0)
        } else if (type == "c") {
            out <- integer(0)
            dim(out) <- c(0, 1)
        } else {
            out <- list()
        }
    } else if (replace) {
        out <- .Call(
            "next_replace_combinations",
            PACKAGE = "arrangements",
            n,
            r,
            d,
            state,
            x,
            type)
    } else if (n < r) {
        if (type == "r") {
            out <- integer(0)
            dim(out) <- c(0, r)
        } else if (type == "c") {
            out <- integer(0)
            dim(out) <- c(r, 0)
        } else {
            out <- list()
        }
    } else if (is.null(f)) {
        out <- .Call(
            "next_combinations",
            PACKAGE = "arrangements",
            n,
            r,
            d,
            state,
            x,
            type)
    } else {
        out <- .Call(
            "next_multiset_combinations",
            PACKAGE = "arrangements",
            n,
            r,
            d,
            state,
            x,
            as_uint_array(f),
            type)
    }
    out
}

#" @export
combinations <- function(n, r, x=NULL, f=NULL, replace=FALSE, type = "r") {
    if (missing(n)) {
        if (is.null(f) && !is.null(x)) {
            n <- length(x)
        } else if (!is.null(f)) {
            n <- sum(f)
        }
    }
    next_combinations(n, r, -1L, NULL, x, f, replace, type)
}


#" @export
icombinations <- function(n, r, x=NULL, f=NULL, replace = FALSE) {
    if (missing(n)) {
        if (is.null(f) && !is.null(x)) {
            n <- length(x)
        } else if (!is.null(f)) {
            n <- sum(f)
        }
    }
    Combinations$new(n, r, x, f, replace)
}

#" @export
ncombinations <- function(n, r, x=NULL, f=NULL, replace=FALSE, bigz=FALSE) {
    if (missing(n)) {
        if (is.null(f) && !is.null(x)) {
            n <- length(x)
        } else if (!is.null(f)) {
            n <- sum(f)
        }
    }
    if (!is.null(f)) {
        f <- as_uint_array(f)
    }
    if (bigz) {
        if (replace) {
            out <- gmp::chooseZ(n + r - 1, r)
        } else if (n < r) {
            out <- 0
        } else if (is.null(f)) {
            out <- gmp::chooseZ(n, r)
        } else {
            out <- .Call("ncomb_f_bigz", PACKAGE = "arrangements", as_uint_array(f), r)
        }

    } else {
        if (replace) {
            out <- choose(n + r - 1, r)
        } else if (n < r) {
            out <- 0
        } else if (is.null(f)) {
            out <- choose(n, r)
        } else {
            out <- .Call("ncomb_f", PACKAGE = "arrangements", as_uint_array(f), r)
        }
    }
    convertz(out, bigz)
}
