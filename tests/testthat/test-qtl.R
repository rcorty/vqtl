context('Functions needed from package qtl')

test_that(desc = 'Simulate cross with qtl::cross.',
          code = {
            test.cross <- qtl::sim.cross(map = qtl::sim.map())
          })

test_that(desc = 'Calculate genotype probabilities with qtl::calc.genoprob.',
          code = {
            test.cross <- qtl::sim.cross(map = qtl::sim.map())
            test.cross <- qtl::calc.genoprob(cross = test.cross)
          })

test_that(desc = 'Get number of individuals with qtl::nind.',
          code = {
            test.cross <- qtl::sim.cross(map = qtl::sim.map())
            n <- qtl::nind(object = test.cross)
          })