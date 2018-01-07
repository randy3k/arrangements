convertz <- function(z, bigz){
    if (bigz){
        return(gmp::as.bigz(z))
    } else if (abs(z) < .Machine$integer.max) {
        return(as.integer(round(z)))
    } else {
        stop("integer overflow, consider using big integer")
    }
}

validate_n <- function(n, x, f){
    if (missing(n)) {
        if (is.null(f) && !is.null(x)) {
            n <- length(x)
        } else if (!is.null(f)) {
            n <- sum(f)
        }
    } else {
        if (is.null(f) && !is.null(x) && n != length(x)) {
            stop("n does not equal to length(x)")
        } else if (!is.null(f) && n != sum(f)) {
            stop("n does not equal to sum(f)")
        }
    }
    (n > 0 && n %% 1 == 0) || stop("n should be a positive integer")
    n
}

validate_r <- function(r, n, replace=TRUE) {
    if (missing(r)) {
        r <- n
    }
    (r > 0 && r %% 1 == 0) || stop("m should be a positive integer")
    replace || r <= n || stop("r is greater than n")
    r
}
