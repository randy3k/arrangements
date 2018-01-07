context("Permutations")

test_that("Permutations - npermutations", {
    expect_equal(npermutations(7), 5040)
    expect_equal(npermutations(x = LETTERS[1:7]), 5040)
    expect_error(npermutations(13), "integer overflow")
    expect_error(npermutations(x = LETTERS[1:13]), "integer overflow")
    expect_equal(npermutations(13, bigz = TRUE), gmp::as.bigz("6227020800"))

})
