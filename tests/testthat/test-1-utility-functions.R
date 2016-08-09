context("Testing utility functions")

test_that(desc = "Functions inhreited from package 'qtl'",
          code = {
            test.cross <- qtl::sim.cross(map = qtl::sim.map(), n.ind = 75)
            expect_s3_class(object = test.cross,
                            class = c('f2', 'cross'))
            expect_s3_class(object = qtl::calc.genoprob(cross = test.cross),
                            class = c('f2', 'cross'))
            expect_equal(object = qtl::nind(test.cross),
                         expected = 75)
          })

test_that(desc = 'is.cross()',
          code = {
            x <- 27599

            expect_false(object = is.cross(x = x))

            class(x) <- 'cross'
            expect_false(object = is.cross(x = x))

            names(x) <- 'pheno'
            expect_false(object = is.cross(x = x))

            x <- list(pheno = 1:5, geno = 1:5)
            class(x) <- 'cross'
            expect_false(object = is.cross(x = x))

            y <- qtl::sim.cross(map = qtl::sim.map())
            y[['pheno']] <- y[['pheno']][-17,]
            expect_false(object = is.cross(x = y))

            y <- qtl::sim.cross(map = qtl::sim.map())
            y[['geno']][[1]][['map']] <- y[['geno']][[1]][['map']][-1]
            expect_false(object = is.cross(x = y))

            z <- qtl::sim.cross(map = qtl::sim.map())
            expect_true(object = is.cross(x = z))
          })


test_that(desc = 'is.f2.cross()',
          code = {
            x <- qtl::sim.cross(map = qtl::sim.map(), type = 'bc')
            expect_false(object = is.f2.cross(x = x))

            y <- qtl::sim.cross(map = qtl::sim.map(sex.sp = TRUE), type = '4way')
            expect_false(object = is.f2.cross(x = y))

            z <- qtl::sim.cross(map = qtl::sim.map(), type = 'f2')
            expect_true(object = is.f2.cross(x = z))
          })


test_that(
  desc = 'testing scanonevar',
  code = {
    set.seed(27599)
    test.cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(20, 4), n.mar = 5))
    test.cross <- qtl::calc.genoprob(cross = test.cross, step = 2)

    # the simplest model
    sov_result1 <- scanonevar(cross = test.cross)

    expect_false(object = is.scanonevar(x = 3))

    expect_true(object = is.scanonevar(x = sov_result1))
  }
)