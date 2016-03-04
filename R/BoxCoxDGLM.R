BoxCoxDGLM <- function(mformula,
                       vformula,
                       data,
                       lambdas = seq(from = -2, to = 2, by = 0.1),
                       plotit = TRUE) {

  response <- as.character(mformula[[2]])
  orig.response <- data[[response]]

  lls <- rep(NA, length(lambdas))
  for (lambda.idx in 1L:length(lambdas)) {

    data[[response]] <- BoxCoxTransformY(y = orig.response,
                                         lambda = lambdas[lambda.idx])

    fit <- tryCatch(expr = dglm(formula = mformula, dformula = vformula, data = data, return.null.deviance = FALSE),
                    finally = NA)
    lls[lambda.idx] <- -0.5*fit$m2loglik
  }

  if (plotit) {
    plot(x = lambdas, y = lls, type = 'l')
  }

  return(list(lambdas, lls))
}