# context('Testing plots')
#
# set.seed(27599)
# test.cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(40, 5), n.mar = 10, eq.spacing = TRUE, include.x = FALSE),
#                              n.ind = 400,
#                              type = 'f2')
# test.cross[['pheno']][['sex']] <- factor(sample(x = c('female', 'male'),
#                                                 size = qtl::nind(test.cross),
#                                                 replace = TRUE))
# test.cross[['pheno']][['phenotype1']] <- rnorm(n = qtl::nind(test.cross))
# test.cross[['pheno']][['phenotype2']] <- rnorm(n = qtl::nind(test.cross),
#                                                mean = 0.2*as.numeric(test.cross$pheno$sex) + 0.5*test.cross$geno$`2`$data[,3],
#                                                sd = exp(0.5*as.numeric(test.cross$pheno$sex) + 0.4*test.cross$geno$`3`$data[,3] - 2))
#
#
#
# test_that(
#   desc = 'phenotype_at_marker_plot should run without error',
#   code = {
#     phenotype_at_marker_plot(cross = test.cross,
#                              phenotype_name = 'phenotype1',
#                              marker_name = 'D1M1',
#                              Ibars = FALSE)
#
#     phenotype_at_marker_plot(cross = test.cross,
#                              phenotype_name = 'phenotype1',
#                              marker_name = 'D1M1',
#                              color_by = 'sex',
#                              point_alpha = 0.5)
#
#     phenotype_at_marker_plot(cross = test.cross,
#                              phenotype_name = 'phenotype1',
#                              marker_name = 'D1M1',
#                              color_by = 'sex')
#
#     phenotype_at_marker_plot(cross = test.cross,
#                              phenotype_name = 'phenotype1',
#                              marker_name = 'D1M1',
#                              color_by = 'sex',
#                              genotype_labels = c('AA', 'AB', 'BB'))
#   }
# )
#
#
#
#
# test.cross <- qtl::calc.genoprob(cross = test.cross, step = 5)
#
# sov <- vqtl::scanonevar(cross = test.cross,
#                         mean.formula = phenotype2 ~ sex + mean.QTL.add + mean.QTL.dom,
#                         var.formula = ~ sex + var.QTL.add + var.QTL.dom,
#                         return.covar.effects = TRUE)
#
# so <- qtl::scanone(cross = test.cross,
#                    pheno.col = 'phenotype2',
#                    addcovar = as.numeric(test.cross$pheno$sex))
#
# # note: there's no testing of plotting with p-values bc it takes too long to do the permutation scans
# # todo: include a sov and so with perms in data/ to test plotting
#
# test_that(
#   desc = 'plot.scanonevar',
#   code = {
#     plot(x = sov, y = so)
#   }
# )
#
# test_that(
#   desc = 'mean_var_sample_plot',
#   code = {
#     mean_var_plot_model_free(cross = test.cross,
#                              phenotype.name = 'phenotype1',
#                              grouping.factor.names = c('sex', 'D3M3'))
#   }
# )
#
#
# test_that(
#   desc = 'mean_var_predictive_plot',
#   code = {
#     mean_var_plot_model_based(cross = test.cross,
#                               phenotype.name = 'phenotype1',
#                               focal.groups = 'D3M3',
#                               se_line_size = 2,
#                               point_size = 7)
#
#     mean_var_plot_model_based(cross = test.cross,
#                               phenotype.name = 'phenotype1',
#                               focal.groups = c('sex', 'D3M3'))
#
#     mean_var_plot_model_based(cross = test.cross,
#                               phenotype.name = 'phenotype1',
#                               focal.groups = 'D3M3',
#                               nuisance.groups = 'D2M2')
#
#     mean_var_plot_model_based(cross = test.cross,
#                               phenotype.name = 'phenotype1',
#                               focal.groups = c('D3M3'),
#                               nuisance.groups = 'sex')
#   }
# )
#
#
#
# test_that(
#   desc = 'effects_plot',
#   code = {
#     effects_over_genome_plot(sov = sov, covar_name_regex = 'sex')
#     effects_over_genome_plot(sov = sov, covar_name_regex = 'mean.QTL')
#     effects_over_genome_plot(sov = sov, covar_name_regex = 'var.QTL')
#     effects_over_genome_plot(sov = sov, covar_name_regex = 'QTL.add')
#     effects_over_genome_plot(sov = sov, covar_name_regex = 'QTL.dom')
#     effects_over_genome_plot(sov = sov, covar_name_regex = 'QTL')
#     effects_over_genome_plot(sov = sov)
#
#     effects_over_genome_plot(sov = sov, covar_name_regex = 'sex', effect_type_regex = 'mean')
#     effects_over_genome_plot(sov = sov, covar_name_regex = 'QTL.add', effect_type_regex = 'var')
#     effects_over_genome_plot(sov = sov, covar_name_regex = 'QTL.dom', effect_type_regex = 'mean')
#     effects_over_genome_plot(sov = sov, covar_name_regex = 'QTL', effect_type_regex = 'var')
#     effects_over_genome_plot(sov = sov, effect_type_regex = 'mean')
#   }
# )
#
