context('Testing plots')

set.seed(27599)
test.cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(40, 5), n.mar = 10, eq.spacing = TRUE, include.x = FALSE),
                             n.ind = 400,
                             type = 'f2')
test.cross[['pheno']][['sex']] <- factor(sample(x = c('female', 'male'),
                                                size = qtl::nind(test.cross),
                                                replace = TRUE))
test.cross[['pheno']][['phenotype1']] <- rnorm(n = qtl::nind(test.cross))
test.cross[['pheno']][['phenotype2']] <- rnorm(n = qtl::nind(test.cross),
                                               mean = 0.2*as.numeric(test.cross$pheno$sex) + 0.5*test.cross$geno$`2`$data[,3],
                                               sd = exp(0.5*as.numeric(test.cross$pheno$sex) + 0.8*test.cross$geno$`3`$data[,3] - 3))

test.cross <- qtl::calc.genoprob(cross = test.cross, step = 2)

sov <- vqtl::scanonevar(cross = test.cross,
                  mean.formula = phenotype2 ~ sex + mean.QTL.add + mean.QTL.dom,
                  var.formula = ~ sex + var.QTL.add + var.QTL.dom,
                  return.covar.effects = TRUE)

sovp <-vqtl::scanonevar.perm(sov = sov,
                       n.perms = 100,
                       n.cores = 3)

so <- qtl::scanone(cross = test.cross,
                   pheno.col = 'phenotype2',
                   addcovar = as.numeric(test.cross$pheno$sex))

sop <- qtl::scanone(cross = test.cross,
                    pheno.col = 'phenotype2',
                    addcovar = as.numeric(test.cross$pheno$sex),
                    n.perm = 100)

the.evd <- evd::fgev(x = sop)

so$empir.p <- evd::pgev(q = so$lod,
                        loc = fitted(the.evd)[1],
                        scale = fitted(the.evd)[2],
                        shape = fitted(the.evd)[3],
                        lower.tail = FALSE)


test_that(
  desc = 'plot.scanonevar',
  code = {
    plot(x = sovp, y = so)
  }
)

test_that(
  desc = 'mean_var_sample_plot',
  code = {
    mean_var_plot_model_free(cross = test.cross,
                             phenotype.name = 'phenotype1',
                             grouping.factor.names = c('sex', 'D3M3'))
  }
)


test_that(
  desc = 'mean_var_predictive_plot',
  code = {
    mean_var_plot_model_based(cross = test.cross,
                              phenotype.name = 'phenotype1',
                              focal.groups = c('D3M3'),
                              nuisance.groups = 'sex')
  }
)



test_that(
  desc = 'effects_plot',
  code = {
    effects_plot(sov = sov, effect.names = 'sex')
    effects_plot(sov = sov, effect.names = 'mean.QTL')
    effects_plot(sov = sov, effect.names = 'var.QTL')
    effects_plot(sov = sov, effect.names = 'QTL.add')
    effects_plot(sov = sov, effect.names = 'QTL.dom')
    effects_plot(sov = sov, effect.names = 'QTL')
    effects_plot(sov = sov)

    effects_plot(sov = sov, effect.names = 'sex', mean.or.var = 'mean')
    effects_plot(sov = sov, effect.names = 'QTL.add', mean.or.var = 'var')
    effects_plot(sov = sov, effect.names = 'QTL.dom', mean.or.var = 'mean')
    effects_plot(sov = sov, effect.names = 'QTL', mean.or.var = 'var')
    effects_plot(sov = sov, mean.or.var = 'mean')
  }
)

# test_that(
#   desc = 'effects_plot',
#   code = {
#
#     x <- scanonevar(cross = test.cross,
#                     mean.formula = phenotype2 ~ sex + mean.QTL.add + mean.QTL.dom,
#                     var.formula = ~ sex + var.QTL.add + var.QTL.dom,
#                     return.covar.effects = TRUE)
#
#     effects_plot
#   }
# )


#
# test_that(
#   desc = 'plot.scanonevar should run without error',
#   code = {
#     set.seed(27599)
#     test.cross <- qtl::sim.cross(map = qtl::sim.map(len = sort(round(runif(n = 15, min = 20, max = 50)), decreasing = TRUE)),
#                                  n.ind = 50)
#     test.cross <- qtl::calc.genoprob(test.cross, step = 2)
#     test.sov <- scanonevar(cross = test.cross)
#     test.so <- qtl::scanone(cross = test.cross)
#
#     p1 <- plot(x = test.sov, plot.title = 'My First Great Plot')
#     expect_s3_class(object = p1, class = 'ggplot')
#
#     p2 <- plot(x = test.sov, y = test.so, plot.title = 'My Second Great Plot')
#     expect_s3_class(object = p2, class = 'ggplot')
#
#
#
#     p3 <- plot(x = test.sov, plotting.units = 'asymp.p', plot.title = 'My Third Great Plot')
#     expect_s3_class(object = p3, class = 'ggplot')
#
#     p4 <- plot(x = test.sov, y = test.so, plotting.units = 'asymp.p', plot.title = 'My Fourth Great Plot')
#     expect_s3_class(object = p4, class = 'ggplot')
#
#     # do more test if you can think of how to do them
#   }
# )
#
#
#
# test_that(
#   desc = 'margin.plot should run without error',
#   code = {
#     set.seed(27599)
#     test.cross <- qtl::sim.cross(map = qtl::sim.map(len = sort(round(runif(n = 15, min = 20, max = 50)), decreasing = TRUE)))
#     test.cross[['pheno']][['sex']] <- sample(x = c('red', 'blue'), size = qtl::nind(test.cross), replace = TRUE)
#
#     margin.plot(cross = test.cross,
#                 focal.phenotype.name = 'phenotype',
#                 marginal.marker.names = list('D1M1'),
#                 col = test.cross[['pheno']][['sex']])
#
#     # do more tests
#   }
# )
#
#
#
#
# test_that(
#   desc = 'predictive.plot should run without error',
#   code = {
#     set.seed(27599)
#     test.cross <- qtl::sim.cross(map = qtl::sim.map(len = sort(round(runif(n = 15, min = 20, max = 50)), decreasing = TRUE)))
#     test.cross[['pheno']][['sex']] <- factor(sample(x = c('red', 'blue'), size = qtl::nind(test.cross), replace = TRUE))
#     test.cross[['pheno']][['handler']] <- factor(sample(x = c('mark', 'greg', 'amy', 'susie'), size = qtl::nind(test.cross), replace = TRUE))
#
#     prediction.ci.plot(cross = test.cross,
#                        mean.formula = phenotype ~ sex + D1M2,
#                        var.formula = ~ sex + D3M1,
#                        group.by = list('sex', 'D2M2'))
#     # do more tests
#   }
# )