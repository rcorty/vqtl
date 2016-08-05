context("Testing scanonevar")


test_that(desc = 'testing scanonevar',
          code = {
            test.cross <- qtl::sim.cross(map = qtl::sim.map())
            test.cross[['pheno']][['sex']] <- sample(x = c(0, 1),
                                                     size = qtl::nind(test.cross),
                                                     replace = TRUE)
            scanonevar(cross = test.cross)
          })
