context("Testing scanonevar.perm")

test_that(
  desc = 'testing scanonevar.perm',
  code = {
    set.seed(27599)
    test.cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(20, 5)))
    test.cross[['pheno']][['sex']] <- sample(x = c(0, 1),
                                             size = qtl::nind(test.cross),
                                             replace = TRUE)


    scanonevar(cross = test.cross,
               mean.formula = phenotype ~ sex + D1M2 + mean.QTL.add + mean.QTL.dom,
               var.formula = ~ sex + D2M3 + var.QTL.add + var.QTL.dom)

  }
)