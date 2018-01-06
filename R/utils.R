convertz <- function(z, bigz){
    if (bigz){
        return(gmp::as.bigz(z))
    } else if (abs(z) < .Machine$integer.max) {
        return(as.integer(z))
    } else {
        stop("integer overflow, consider using big interger")
    }
}

check_nrxf <- function(n, r, x, f, replace){
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
    if (!replace && r > n) {
        stop("r is greater than n")
    } else if (r < 0) {
        stop("r should be non-negative")
    }
    n
}
