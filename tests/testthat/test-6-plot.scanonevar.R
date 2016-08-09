context('Testing plot.scanonevar')

test_that(
  desc = 'plot.scanonevar should run without error',
  code = {
    set.seed(27599)
    test.cross <- qtl::sim.cross(map = qtl::sim.map(len = sort(round(runif(n = 15, min = 20, max = 50)), decreasing = TRUE)),
                                 n.ind = 50)
    test.cross <- qtl::calc.genoprob(test.cross, step = 2)
    test.sov <- scanonevar(cross = test.cross)
    test.so <- qtl::scanone(cross = test.cross)

    p1 <- plot(x = test.sov, plot.title = 'My First Great Plot')
    expect_s3_class(object = p1, class = 'ggplot')

    p2 <- plot(x = test.sov, y = test.so, plot.title = 'My Second Great Plot')
    expect_s3_class(object = p2, class = 'ggplot')



    p3 <- plot(x = test.sov, plotting.units = 'asymp.p', plot.title = 'My Third Great Plot')
    expect_s3_class(object = p3, class = 'ggplot')

    p4 <- plot(x = test.sov, y = test.so, plotting.units = 'asymp.p', plot.title = 'My Fourth Great Plot')
    expect_s3_class(object = p4, class = 'ggplot')

    # do more test if you can think of how to do them
  }
)