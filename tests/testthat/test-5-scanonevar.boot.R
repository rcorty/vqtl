context("Testing scanonevar.boot")

set.seed(27599)
test.cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(40, 5), n.mar = 10, eq.spacing = FALSE, include.x = FALSE),
                             n.ind = 500)
test.cross[['pheno']][['sex']] <- sample(x = c(0, 1),
                                         size = qtl::nind(test.cross),
                                         replace = TRUE)
test.cross[['pheno']][['sire']] <- factor(x = sample(x = 1:20,
                                                     size = qtl::nind(test.cross),
                                                     replace = TRUE))

test.cross <- qtl::calc.genoprob(cross = test.cross, step = 3)

test.cross[['pheno']][['phenotype']] <-rnorm(n = qtl::nind(test.cross), mean = 0.5*test.cross$geno$`1`$data[,4])

test_that(
  desc = 'testing scanonevar.boot with dglm with gaussian model',
  code = {

    sov <- vqtl::scanonevar(cross = test.cross,
                            mean.formula = phenotype ~ sex + mean.QTL.add + mean.QTL.dom,
                            var.formula = ~ sex + var.QTL.add + var.QTL.dom)

    bs <- vqtl::scanonevar.boot(sov = sov, n.resamples = 20, chr = 1, qtl_type = 'mQTL')

    expect_equal(object = length(bs$bootstrap_maxes), expected = 20)
    expect_true(object = all(bs$bootstrap_maxes >= 0))
    expect_true(object = all(bs$bootstrap_maxes <= 40))


    bs <- vqtl::scanonevar.boot(sov = sov, n.resamples = 20, chr = 1, qtl_type = 'vQTL')

    expect_equal(object = length(bs$bootstrap_maxes), expected = 20)
    expect_true(object = all(bs$bootstrap_maxes >= 0))
    expect_true(object = all(bs$bootstrap_maxes <= 40))
  }
)