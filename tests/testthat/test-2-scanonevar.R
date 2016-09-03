context("Testing scanonevar")

test_that(
  desc = 'testing scanonevar',
  code = {
    set.seed(27599)
    test.cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(20, 5)))
    test.cross[['pheno']][['sex']] <- sample(x = c(0, 1),
                                             size = qtl::nind(test.cross),
                                             replace = TRUE)


    sov <- scanonevar(cross = test.cross,
                      mean.formula = phenotype ~ sex + D1M2 + mean.QTL.add + mean.QTL.dom,
                      var.formula = ~ sex + D2M3 + var.QTL.add + var.QTL.dom)

    # should be a scanonevar object
    expect_true(is.scanonevar(sov))

    # joint lods should be pointwise higher than mean and variance lods
    expect_true(all(sov$result %>% .[['joint.lod']] >= sov$result %>% .[['mean.lod']]))
    expect_true(all(sov$result %>% .[['joint.lod']] >= sov$result %>% .[['var.lod']]))

  }
)



# context("Testing validate.scanonevar.input_")
#
# test_that(desc = "Testing is.cross, the first check in validate.scanonevar.input_",
#           code = {
#             expect_error(object = validate.scanonevar.input_(cross = 27599,
#                                                              mean.formula = 27599,
#                                                              var.formula = 27599,
#                                                              chrs = 'applesauce'))
#           })
#
# test_that(
#   desc = 'Testing checks that mean.formula and var.formula are of the appropriate shape',
#   code = {
#     test.cross <- qtl::sim.cross(map = qtl::sim.map())
#     test.cross[['pheno']][['sex']] <- sample(x = c(0, 1),
#                                              size = qtl::nind(test.cross),
#                                              replace = TRUE)
#
#
#     # mean.formula has no RHS
#     expect_error(object = validate.scanonevar.input_(cross = test.cross,
#                                                      mean.formula =  ~ mean.QTL.add + mean.QTL.dom,
#                                                      var.formula = ~ var.QTL.add + var.QTL.dom,
#                                                      chrs = names(test.cross[['geno']])))
#
#     # var.formula has a RHS
#     expect_error(object = validate.scanonevar.input_(cross = test.cross,
#                                                      mean.formula =  phenotype ~ mean.QTL.add + mean.QTL.dom,
#                                                      var.formula = phenotype ~ var.QTL.add + var.QTL.dom,
#                                                      chrs = names(test.cross[['geno']])))
#   }
# )
#
#
#
# test_that(
#   desc = paste('Testing the checks that all variables excpet mean.QTL.add,',
#                'mean.QTL.dom, var.QTL.add, and var.QTL.dom are either the name',
#                'of a genetic marker followed by _add or _dom or the name of a phenotype'),
#   code = {
#     test.cross <- qtl::sim.cross(map = qtl::sim.map())
#     test.cross[['pheno']][['sex']] <- sample(x = c(0, 1),
#                                              size = qtl::nind(test.cross),
#                                              replace = TRUE)
#
#     # mean.formula RHS isn't a phenotype in cross
#     expect_error(object = validate.scanonevar.input_(cross = test.cross,
#                                                      mean.formula = applesauce ~ mean.QTL.add + mean.QTL.dom,
#                                                      var.formula = ~ var.QTL.add + var.QTL.dom,
#                                                      chrs = names(test.cross[['geno']])))
#
#     # mean.formula has a covariate thats not a keyword, phenotype, nor marker
#     expect_error(object = validate.scanonevar.input_(cross = test.cross,
#                                                      mean.formula = phenotype ~ applesauce + mean.QTL.add + mean.QTL.dom,
#                                                      var.formula = ~ var.QTL.add + var.QTL.dom,
#                                                      chrs = names(test.cross[['geno']])))
#
#     # var.formula has a covariate thats not a keyword, phenotype, nor marker
#     expect_error(object = validate.scanonevar.input_(cross = test.cross,
#                                                      mean.formula = phenotype ~ mean.QTL.add + mean.QTL.dom,
#                                                      var.formula = ~ applesauce + var.QTL.add + var.QTL.dom,
#                                                      chrs = names(test.cross[['geno']])))
#
#     # mean.formula uses a variance keyword
#     expect_error(object = validate.scanonevar.input_(cross = test.cross,
#                                                      mean.formula = phenotype ~ var.QTL.add + mean.QTL.dom,
#                                                      var.formula = ~ var.QTL.add + var.QTL.dom,
#                                                      chrs = names(test.cross[['geno']])))
#
#     # var.formula uses a mean keyword
#     expect_error(object = validate.scanonevar.input_(cross = test.cross,
#                                                      mean.formula = applesauce ~  mean.QTL.add + mean.QTL.dom,
#                                                      var.formula = ~ mean.QTL.add + var.QTL.dom,
#                                                      chrs = names(test.cross[['geno']])))
#
#     # same mean keyword appears more than once in mean.formula
#     expect_error(object = validate.scanonevar.input_(cross = test.cross,
#                                                      mean.formula = applesauce ~ mean.QTL.add + mean.QTL.add,
#                                                      var.formula = ~ mean.QTL.add + var.QTL.dom,
#                                                      chrs = names(test.cross[['geno']])))
#
#
#     # same var keyword appears more than once in var.formula
#     expect_error(object = validate.scanonevar.input_(cross = test.cross,
#                                                      mean.formula = applesauce ~ mean.QTL.add + mean.QTL.dom,
#                                                      var.formula = ~ mean.QTL.add + var.QTL.add,
#                                                      chrs = names(test.cross[['geno']])))
#
#     # no mean keywords nor var keywords
#     expect_error(object = validate.scanonevar.input_(cross = test.cross,
#                                                      mean.formula = applesauce ~ 1,
#                                                      var.formula = ~ 1,
#                                                      chrs = names(test.cross[['geno']])))
#
#
#
#
#     # simplest valid input (defaults of scanonevar)
#     expect_true(object = validate.scanonevar.input_(cross = test.cross,
#                                                     mean.formula = phenotype ~ mean.QTL.add + mean.QTL.dom,
#                                                     var.formula = ~ var.QTL.add + var.QTL.dom,
#                                                     chrs = names(test.cross[['geno']])))
#
#     # add a phenotype covariate to mean.formula
#     expect_true(object = validate.scanonevar.input_(cross = test.cross,
#                                                     mean.formula = phenotype ~ sex + mean.QTL.add + mean.QTL.dom,
#                                                     var.formula = ~ var.QTL.add + var.QTL.dom,
#                                                     chrs = names(test.cross[['geno']])))
#
#     # ada a phenotype covariate to var.formula
#     expect_true(object = validate.scanonevar.input_(cross = test.cross,
#                                                     mean.formula = phenotype ~ mean.QTL.add + mean.QTL.dom,
#                                                     var.formula = ~ sex + var.QTL.add + var.QTL.dom,
#                                                     chrs = names(test.cross[['geno']])))
#
#     # add a marker covariate to mean.formula
#     expect_true(object = validate.scanonevar.input_(cross = test.cross,
#                                                     mean.formula = phenotype ~ D1M3_add + mean.QTL.add + mean.QTL.dom,
#                                                     var.formula = ~ var.QTL.add + var.QTL.dom,
#                                                     chrs = names(test.cross[['geno']])))
#
#     # add a marker covariate to var.formula
#     expect_true(object = validate.scanonevar.input_(cross = test.cross,
#                                                     mean.formula = phenotype ~ mean.QTL.add + mean.QTL.dom,
#                                                     var.formula = ~ D1M3_add + var.QTL.add + var.QTL.dom,
#                                                     chrs = names(test.cross[['geno']])))
#
#     # add both types of covariates to both formulae
#     expect_true(object = validate.scanonevar.input_(cross = test.cross,
#                                                     mean.formula = phenotype ~ sex + D1M3_dom + mean.QTL.add + mean.QTL.dom,
#                                                     var.formula = ~ sex + D1M3_add + var.QTL.add + var.QTL.dom,
#                                                     chrs = names(test.cross[['geno']])))
#
#     # mean keywords but no var keywords
#     expect_true(object = validate.scanonevar.input_(cross = test.cross,
#                                                     mean.formula = phenotype ~ sex + D1M3_add + mean.QTL.add + mean.QTL.dom,
#                                                     var.formula = ~ sex + D1M3_dom,
#                                                     chrs = names(test.cross[['geno']])))
#
#     # var keywords but no mean keywords
#     expect_true(object = validate.scanonevar.input_(cross = test.cross,
#                                                     mean.formula = phenotype ~ sex + D1M3_add,
#                                                     var.formula = ~ sex + D1M3_dom + var.QTL.add + var.QTL.dom,
#                                                     chrs = names(test.cross[['geno']])))
#
#     # mean.formula has a marker with neither '_add' nor '_dom' appended
#     expect_true(object = validate.scanonevar.input_(cross = test.cross,
#                                                     mean.formula = phenotype ~ D1M2 + mean.QTL.add + mean.QTL.dom,
#                                                     var.formula = ~ var.QTL.add + var.QTL.dom,
#                                                     chrs = names(test.cross[['geno']])))
#
#
#     # var.formula has a marker with neither '_add' nor '_dom' appended
#     expect_true(object = validate.scanonevar.input_(cross = test.cross,
#                                                     mean.formula = phenotype ~ mean.QTL.add + mean.QTL.dom,
#                                                     var.formula = ~ D1M3 + var.QTL.add + var.QTL.dom,
#                                                     chrs = names(test.cross[['geno']])))
#
#
#
#   }
# )
#
#
#
#
#
# test_that(
#   desc = 'Testing the checks that all requested chrs are in the cross',
#   code = {
#     test.cross <- qtl::sim.cross(map = qtl::sim.map())
#     test.cross[['pheno']][['sex']] <- sample(x = c(0, 1),
#                                              size = qtl::nind(test.cross),
#                                              replace = TRUE)
#
#     expect_error(object = validate.scanonevar.input_(cross = test.cross,
#                                                      mean.formula = phenotype ~ mean.QTL.add + mean.QTL.dom,
#                                                      var.formula = ~ var.QTL.add + var.QTL.dom,
#                                                      chrs = 'z'))
#   }
# )
