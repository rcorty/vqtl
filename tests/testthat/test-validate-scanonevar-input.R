context("Testing validate.scanonevar.input_")

test_that(desc = "Testing is.cross, the first check in validate.scanonevar.input_",
          code = {
            expect_error(object = validate.scanonevar.input_(cross = 27599,
                                                             mean.formula = 27599,
                                                             var.formula = 27599))
          })

test_that(desc = 'Testing the checks that all variables excpet mean.QTL.add, mean.QTL.dom, var.QTL.add, and var.QTL.dom are either the name of a phenotype or the name of a genetic marker',
          code = {
            test.cross <- qtl::sim.cross(map = qtl::sim.map())
            test.cross[['pheno']][['sex']] <- sample(x = c(0, 1),
                                                     size = qtl::nind(test.cross),
                                                     replace = TRUE)
            expect_error(object = validate.scanonevar.input_(cross = test.cross,
                                                             mean.formula = formula(applesauce ~ sex + mean.QTL.add + mean.QTL.dom),
                                                             var.formula = ~ var.QTL.add + var.QTL.dom))

            expect_error(object = validate.scanonevar.input_(cross = test.cross,
                                                             mean.formula = formula(phenotype ~ applesauce + mean.QTL.add + mean.QTL.dom),
                                                             var.formula = ~ var.QTL.add + var.QTL.dom))

            expect_error(object = validate.scanonevar.input_(cross = test.cross,
                                                             mean.formula = formula(applesauce ~ sex + mean.QTL.add + mean.QTL.dom),
                                                             var.formula = ~ applesauce + var.QTL.add + var.QTL.dom))

            expect_true(object = validate.scanonevar.input_(cross = test.cross,
                                                            mean.formula = formula(phenotype ~ sex + mean.QTL.add + mean.QTL.dom),
                                                            var.formula = ~ var.QTL.add + var.QTL.dom))

            expect_true(object = validate.scanonevar.input_(cross = test.cross,
                                                            mean.formula = formula(phenotype ~ sex + D1M3 + mean.QTL.add + mean.QTL.dom),
                                                            var.formula = ~ var.QTL.add + var.QTL.dom))

            expect_true(object = validate.scanonevar.input_(cross = test.cross,
                                                            mean.formula = formula(phenotype ~ sex + mean.QTL.add + mean.QTL.dom),
                                                            var.formula = ~ D2M1 + var.QTL.add + var.QTL.dom))
          })
