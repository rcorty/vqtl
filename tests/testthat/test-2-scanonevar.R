context("Testing scanonevar")

test_that(
  desc = 'testing scanonevar',
  code = {

    set.seed(27599)
    test.cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(20, 5), eq.spacing = FALSE))
    test.cross[['pheno']][['sex']] <- sample(x = c(0, 1),
                                             size = qtl::nind(test.cross),
                                             replace = TRUE)
    test.cross[['pheno']][['sire']] <- factor(x = sample(x = 1:5,
                                                         size = qtl::nind(test.cross),
                                                         replace = TRUE))


    test.cross <- qtl::calc.genoprob(cross = test.cross, step = 5)

    sov <- scanonevar(cross = test.cross,
                      mean.formula = phenotype ~ sex + D1M2 + mean.QTL.add + mean.QTL.dom,
                      var.formula = ~ sire + D2M3 + var.QTL.add + var.QTL.dom)

    # should be a scanonevar object
    expect_true(is.scanonevar(sov))

    # joint lods should be pointwise higher than mean and variance lods
    expect_true(all(sov$result %>% .[['joint.lod']] >= sov$result %>% .[['mean.lod']], na.rm = TRUE))
    expect_true(all(sov$result %>% .[['joint.lod']] >= sov$result %>% .[['var.lod']], na.rm = TRUE))


    x <- summary(object = sov, units = 'lod')

    # components of summary
    expect_true(all(c('mQTL', 'vQTL', 'mvQTL') %in% names(x)))

    # absolute maximum size of each component of summary is size of result
    expect_true(nrow(x$mQTL) <= nrow(sov$result))
    expect_true(nrow(x$vQTL) <= nrow(sov$result))
    expect_true(nrow(x$mvQTL) <= nrow(sov$result))

    # should be no randomness in summary -- recalculate it and it's the same
    expect_identical(object = x, expected = summary(object = sov, units = 'lod'))


    # the factor, 'sire', should be handled correctly
    sov1 <- scanonevar(cross = test.cross,
                       mean.formula = phenotype ~ sire + mean.QTL.add + mean.QTL.dom,
                       var.formula = ~ sire + D2M3 + var.QTL.add + var.QTL.dom,
                       return.covar.effects = TRUE)

    expect_true(is.scanonevar(sov))
  }
)


