convertz <- function(z, bigz){
    if (bigz){
        return(gmp::as.bigz(z))
    } else if (abs(z) < .Machine$integer.max) {
        return(as.integer(z))
    } else {
        stop("integer overflow, consider using big interger")
    }
}
