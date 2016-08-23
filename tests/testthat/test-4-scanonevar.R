#
# context("Testing scanonevar")
#
# test_that(
#   desc = 'testing scanonevar',
#   code = {
#     set.seed(27599)
#     test.cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(20, 4), n.mar = 5))
#     test.cross <- qtl::calc.genoprob(cross = test.cross, step = 2)
#     test.cross[['pheno']][['sex']] <- sample(x = c(0, 1), size = qtl::nind(test.cross), replace = TRUE)
#
#     # the simplest model
#     sov_result1 <- scanonevar(cross = test.cross)
#     expect_identical(object = names(sov_result1),
#                      expected = c('loc.name', 'chr', 'pos', 'mean.lod', 'mean.asymp.p', 'var.lod', 'var.asymp.p', 'joint.lod', 'joint.asymp.p'))
#
#     # phenotype and marker covars in mean and var submodels
#     sov_result2 <- scanonevar(cross = test.cross,
#                               mean.formula = phenotype ~ sex + D3M3_add + mean.QTL.add + mean.QTL.dom,
#                               var.formula = ~ sex + D2M1_dom + var.QTL.add + var.QTL.dom)
#     expect_identical(object = names(sov_result2),
#                      expected = c('loc.name', 'chr', 'pos', 'mean.lod', 'mean.asymp.p', 'var.lod', 'var.asymp.p', 'joint.lod', 'joint.asymp.p'))
#
#     # phenotype and marker covars in both submodels, only doing var testing
#     sov_result3 <- scanonevar(cross = test.cross,
#                               mean.formula = phenotype ~ sex + D3M3_dom,
#                               var.formula = ~ sex + D2M3_add + var.QTL.add + var.QTL.dom)
#     expect_identical(object = names(sov_result3),
#                      expected = c('loc.name', 'chr', 'pos', 'var.lod', 'var.asymp.p'))
#
#     # phenotype and marker covars in both submodels, only doing mean testing
#     sov_result4 <- scanonevar(cross = test.cross,
#                               mean.formula = phenotype ~ sex + D3M3_add + mean.QTL.add + mean.QTL.dom,
#                               var.formula = ~ sex + D2M2_add + D2M3_dom)
#     expect_identical(object = names(sov_result4),
#                      expected = c('loc.name', 'chr', 'pos', 'mean.lod', 'mean.asymp.p'))
#
#     # filtering to a subset of chromosomes
#     sov_result5 <- scanonevar(cross = test.cross,
#                               mean.formula = phenotype ~ sex + D3M3_add + mean.QTL.add + mean.QTL.dom,
#                               var.formula = ~ sex + D2M2_add + D2M1_dom + var.QTL.add,
#                               chrs = c('1', 'X'))
#     expect_identical(object = unique(sov_result5[['chr']]),
#                      c('1', 'X'))
#
#
#     # returning covariate effect estimates
#     sov_result6 <- scanonevar(cross = test.cross,
#                               mean.formula = phenotype ~ sex + D3M3_add + mean.QTL.add + mean.QTL.dom,
#                               var.formula = ~ sex + D2M2_add + D2M1_dom + var.QTL.add,
#                               return.covar.effects = TRUE)
#     expect_identical(object = names(sov_result6),
#                      expected = c(c('loc.name', 'chr', 'pos', 'mean.lod', 'mean.asymp.p', 'var.lod',
#                                     'var.asymp.p', 'joint.lod', 'joint.asymp.p', '(Intercept)_mef',
#                                     'sex_mef', 'D3M3_add_mef', 'mean.QTL.add_mef', 'mean.QTL.dom_mef',
#                                     '(Intercept)_vef', 'sex_vef', 'D2M2_add_vef', 'D2M1_dom_vef',
#                                     'var.QTL.add_vef')))
#
#   }
# )
