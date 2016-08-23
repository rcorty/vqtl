context("Testing scanonevar.perm")

test_that(
  desc = 'testing scanonevar.perm',
  code = {
    set.seed(27599)
    test.cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(20, 5), n.mar = 10, eq.spacing = TRUE))
    test.cross[['pheno']][['sex']] <- sample(x = c(0, 1),
                                             size = qtl::nind(test.cross),
                                             replace = TRUE)
    test.cross[['pheno']][['phenotype']] <- rnorm(n = qtl::nind(test.cross),
                                                  mean = test.cross$geno$`2`$data[,3],
                                                  sd = exp(test.cross$geno$`3`$data[,3]))


    x <- scanonevar(cross = test.cross,
                    mean.formula = phenotype ~ sex + D1M2 + mean.QTL.add + mean.QTL.dom,
                    var.formula = ~ sex + D2M3 + var.QTL.add + var.QTL.dom)

    y <- scanonevar.perm(sov = x,
                         n.perms = 50)

  }
)