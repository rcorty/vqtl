context("Testing wrangle.scanonevar.input_")

test_that(
  desc = 'testing wrangle.scan.types_',
  code = {

    # mean only testing
    expect_identical(object = wrangle.scan.types_(mean.formula = applesauce ~ mean.QTL.add,
                                                  var.formula = ~ 1),
                     expected = 'mean')

    expect_identical(object = wrangle.scan.types_(mean.formula = applesauce ~ mean.QTL.dom,
                                                  var.formula = ~ 1),
                     expected = 'mean')

    expect_identical(object = wrangle.scan.types_(mean.formula = applesauce ~ mean.QTL.add + mean.QTL.dom,
                                                  var.formula = ~ 1),
                     expected = 'mean')

    expect_identical(object = wrangle.scan.types_(mean.formula = applesauce ~ mean.QTL.add + mean.QTL.dom,
                                                  var.formula = ~ banana),
                     expected = 'mean')

    # var only testing
    expect_identical(object = wrangle.scan.types_(mean.formula = applesauce ~ 1,
                                                  var.formula = ~ var.QTL.add),
                     expected = 'var')

    expect_identical(object = wrangle.scan.types_(mean.formula = applesauce ~ 1,
                                                  var.formula = ~ var.QTL.dom),
                     expected = 'var')

    expect_identical(object = wrangle.scan.types_(mean.formula = applesauce ~ 1,
                                                  var.formula = ~ var.QTL.add + var.QTL.dom),
                     expected = 'var')

    # mean, var, and joint testing
    expect_identical(object = wrangle.scan.types_(mean.formula = applesauce ~ mean.QTL.add + mean.QTL.dom,
                                                  var.formula = ~ var.QTL.add + var.QTL.dom),
                     expected = c('mean', 'var', 'joint'))

    expect_identical(object = wrangle.scan.types_(mean.formula = applesauce ~ mean.QTL.add,
                                                  var.formula = ~ var.QTL.add),
                     expected = c('mean', 'var', 'joint'))

    expect_identical(object = wrangle.scan.types_(mean.formula = applesauce ~ mean.QTL.add + mean.QTL.dom,
                                                  var.formula = ~ var.QTL.add),
                     expected = c('mean', 'var', 'joint'))

    expect_identical(object = wrangle.scan.types_(mean.formula = applesauce ~ mean.QTL.add,
                                                  var.formula = ~ var.QTL.add + var.QTL.dom),
                     expected = c('mean', 'var', 'joint'))

    expect_identical(object = wrangle.scan.types_(mean.formula = applesauce ~ mean.QTL.add + mean.QTL.dom,
                                                  var.formula = ~ var.QTL.add + var.QTL.dom),
                     expected = c('mean', 'var', 'joint'))


    expect_identical(object = wrangle.scan.types_(mean.formula = applesauce ~ banana + mean.QTL.add + mean.QTL.dom,
                                                  var.formula = ~ banana + var.QTL.add + var.QTL.dom),
                     expected = c('mean', 'var', 'joint'))
  }
)


test_that(
  desc = 'testing wrangle.scan.formulae_',
  code = {
    expect_identical(object = wrangle.scan.formulae_(mean.formula = a ~ b + mean.QTL.add,
                                                     var.formula = ~ c + var.QTL.add),
                     expected = list(mean.alt.formula = a ~ b + mean.QTL.add,
                                     mean.null.formula = a ~ b,
                                     var.alt.formula = ~ c + var.QTL.add,
                                     var.null.formula = ~ c))

    expect_identical(object = wrangle.scan.formulae_(mean.formula = a ~ b,
                                                     var.formula = ~ c + var.QTL.add),
                     expected = list(mean.alt.formula = a ~ b,
                                     var.alt.formula = ~ c + var.QTL.add,
                                     var.null.formula = ~ c))

    expect_identical(object = wrangle.scan.formulae_(mean.formula = a ~ b + mean.QTL.add,
                                                     var.formula = ~ c),
                     expected = list(mean.alt.formula = a ~ b + mean.QTL.add,
                                     mean.null.formula = a ~ b,
                                     var.alt.formula = ~ c))
  }
)


test_that(
  desc = 'testing wrangle.loc.info.df_',
  code = {
    test.cross <- qtl::sim.cross(map = qtl::sim.map())
    test.cross <- qtl::calc.genoprob(cross = test.cross, step = 2)
    info.df <- wrangle.loc.info.df_(cross = test.cross,
                                    chrs = qtl::chrnames(test.cross))
    genoprob.df <- wrangle.genoprob.df_(cross = test.cross)

    expect_identical(object = names(x = info.df),
                     expected = c('loc.name', 'chr', 'pos'))

    expect_identical(object = info.df[['loc.name']],
                     expected = unique(genoprob.df[['loc.name']]))

    expect_identical(object = unique(info.df[['chr']]),
                     expected = names(test.cross[['geno']]))

    info.df2 <- wrangle.loc.info.df_(cross = test.cross,
                                     chrs = c('7', '11', 'X'))

    expect_identical(object = unique(info.df2[['chr']]),
                     expected = c('7', '11', 'X'))

    # do more tests if you can think of how to do them
  }
)


test_that(
  desc = 'testing wrangle.genoprob.df_',
  code = {
    test.cross <- qtl::sim.cross(map = qtl::sim.map())
    test.cross <- qtl::calc.genoprob(cross = test.cross, step = 2)
    genoprob.df <- wrangle.genoprob.df_(cross = test.cross)
    info.df <- wrangle.loc.info.df_(cross = test.cross,
                                    chrs = qtl::chrnames(test.cross))

    expect_identical(object = names(x = genoprob.df),
                     expected = c('iid', 'loc.name', 'allele', 'genoprob'))

    expect_identical(object = unique(genoprob.df[['loc.name']]),
                     expected = info.df[['loc.name']])

    expect_true(object = all(genoprob.df[['genoprob']] < 1))
    expect_true(object = all(genoprob.df[['genoprob']] > 0))

    # do more tests if you can think of how to do them
  }
)


test_that(
  desc = 'testing wrangle.modeling.df_',
  code = {
    test.cross <- qtl::sim.cross(map = qtl::sim.map())
    test.cross <- qtl::calc.genoprob(cross = test.cross, step = 2)
    test.cross[['pheno']][['sex']] <- sample(x = c(0, 1), size = qtl::nind(object = test.cross), replace = TRUE)
    test.cross[['pheno']][['bodyweight']] <- rnorm(n = qtl::nind(object = test.cross), mean = 10, sd = 1)
    genoprob.df <- wrangle.genoprob.df_(cross = test.cross)


    # covariates of all three types (keyword, phenotype, and marker) in
    # mean and variance submodels, but no overlap between submodels
    formulae1 <- wrangle.scan.formulae_(mean.formula = phenotype ~ bodyweight + D1M3_add + D1M3_dom + mean.QTL.add,
                                        var.formula = ~ D5M4_add + sex + var.QTL.add)
    modeling.df1 <- wrangle.modeling.df_(cross = test.cross,
                                         genoprobs = genoprob.df,
                                         scan.formulae = formulae1)
    needed.names1 <- c('phenotype', 'mean.QTL.add', 'var.QTL.add', 'bodyweight', 'sex', 'D1M3_add', 'D1M3_dom', 'D5M4_add')
    expect_true(object = all(names(modeling.df1) %in% needed.names1))
    expect_true(object = all(!is.na(dplyr::select(modeling.df1, -dplyr::matches('QTL')))))


    # covariates of all three types (keyword, phenotype, and marker) in
    # mean and variance submodels, with overlap between submodels
    formulae2 <- wrangle.scan.formulae_(mean.formula = phenotype ~ bodyweight + D1M3_add + mean.QTL.add,
                                        var.formula = ~ bodyweight + D1M3_add + var.QTL.add)
    modeling.df2 <- wrangle.modeling.df_(cross = test.cross,
                                         genoprobs = genoprob.df,
                                         scan.formulae = formulae2)
    needed.names2 <- c('phenotype', 'mean.QTL.add', 'var.QTL.add', 'bodyweight', 'D1M3_add')
    expect_true(object = all(names(modeling.df2) %in% needed.names2))
    expect_true(object = all(!is.na(dplyr::select(modeling.df2, -dplyr::matches('QTL')))))

    # need a test of ability to wrangle df when only mean or only var testing is being done
    # e.g. when mean or var models are NULL
    # do more tests if you can think of how to do them
  }
)