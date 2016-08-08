context("Testing scanonevar")

test_that(
  desc = 'testing scanonevar',
  code = {
    set.seed(27599)
    test.cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(20, 5)))
    test.cross <- qtl::calc.genoprob(cross = test.cross, step = 2)
    test.cross[['pheno']][['sex']] <- sample(x = c(0, 1), size = qtl::nind(test.cross), replace = TRUE)

    sov_result1 <- scanonevar(cross = test.cross)
    expect_equal_to_reference(object = sov_result1,
                              file = 'z_scanonevar_test_result1.RDS')


    sov_result2 <- scanonevar(cross = test.cross,
                              mean.formula = phenotype ~ sex + D3M3_add + mean.QTL.add + mean.QTL.dom,
                              var.formula = ~ sex + D2M7_dom + var.QTL.add + var.QTL.dom)
    expect_equal_to_reference(object = sov_result2,
                              file = 'z_scanonevar_test_result2.RDS')


    sov_result3 <- scanonevar(cross = test.cross,
                              mean.formula = phenotype ~ sex + D3M3_dom,
                              var.formula = ~ sex + D2M7_add + var.QTL.add + var.QTL.dom)
    expect_equal_to_reference(object = sov_result3,
                              file = 'z_scanonevar_test_result3.RDS')


    sov_result4 <- scanonevar(cross = test.cross,
                              mean.formula = phenotype ~ sex + D3M3_add + mean.QTL.add + mean.QTL.dom,
                              var.formula = ~ sex + D2M7_add + D2M7_dom)
    expect_equal_to_reference(object = sov_result4,
                              file = 'z_scanonevar_test_result4.RDS')

  }
)
