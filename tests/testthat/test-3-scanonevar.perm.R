context("Testing scanonevar.perm")

test_that(
  desc = 'testing scanonevar.perm',
  code = {
    set.seed(27599)
    test.cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(40, 5), n.mar = 10, eq.spacing = TRUE, include.x = FALSE),
                                 n.ind = 300,
                                 type = 'f2')
    test.cross[['pheno']][['sex']] <- sample(x = c(0, 1),
                                             size = qtl::nind(test.cross),
                                             replace = TRUE)
    test.cross[['pheno']][['phenotype1']] <- rnorm(n = qtl::nind(test.cross))
    test.cross[['pheno']][['phenotype2']] <- rnorm(n = qtl::nind(test.cross), mean = 0.3*test.cross$geno$`2`$data[,3])
    test.cross[['pheno']][['phenotype4']] <- rnorm(n = qtl::nind(test.cross),
                                                   mean = test.cross$pheno$sex + test.cross$geno$`2`$data[,3],
                                                   sd = exp(test.cross$pheno$sex + test.cross$geno$`3`$data[,3]))


    x <- scanonevar(cross = test.cross,
                    mean.formula = phenotype2 ~ sex + D1M2 + mean.QTL.add + mean.QTL.dom,
                    var.formula = ~ sex + D2M3 + var.QTL.add + var.QTL.dom)

    y1 <- scanonevar.perm(sov = x, n.perms = 10)

    plot(y1)

    expect_true(is.scanonevar(y1))


  }
)