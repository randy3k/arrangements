convertz <- function(z, bigz){
    if (bigz){
        return(gmp::as.bigz(z))
    } else if (abs(z) < .Machine$integer.max) {
        return(as.integer(round(z)))
    } else {
        stop("integer overflow, consider using big integer")
    }
}

as_uint_array <- function(a) {
    .Call("as_uint_array", PACKAGE = "arrangements", a)
}
