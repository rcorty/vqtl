context("Testing wrangle.scanonevar.input_")


test_that(desc = 'testing wrangle.scanonevar.input_',
          code = {
            test.cross <- qtl::sim.cross(map = qtl::sim.map())
            test.cross[['pheno']][['sex']] <- sample(x = c(0, 1),
                                                     size = qtl::nind(test.cross),
                                                     replace = TRUE)


          })
