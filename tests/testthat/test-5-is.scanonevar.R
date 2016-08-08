context("Testing is.scanonevar")

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