context("tests pFranksCopula()")

test_that("pFranksCopula works ok", {

    ## "diagonal" probabilities
    m0 <- matrix(rep(1:15/16, 4), ncol = 4)
    rho <- 1
    d <- 2; expect_equal_to_reference(pFranksCopula(m0[ , 1:d], rho), paste0("pFranksCopula_diag_", d))
    d <- 3; expect_equal_to_reference(pFranksCopula(m0[ , 1:d], rho), paste0("pFranksCopula_diag_", d))
    d <- 4; expect_equal_to_reference(pFranksCopula(m0[ , 1:d], rho), paste0("pFranksCopula_diag_", d))

})



