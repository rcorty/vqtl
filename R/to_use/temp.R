my.formula <- fruit ~ apple * banana + carrot + apple
my.string <- deparse(expr = my.formula)
var.to.replace <- 'apple'
new.terms <- paste0('(',
                    paste0(var.to.replace,
                           c('_part1', '_part2'),
                           collapse = '+'),
                    ')')

new.string <- gsub(pattern = var.to.replace, replacement = new.terms, x = my.string)

formula(new.string)

unknown.function <- function(formula,
                             to.replace,
                             replace.with) {
  replacer <- list(replace.with)
  replacer <- setNames(object = replacer, nm = to.replace)
  result <- do.call(what = 'substitute',
                    args = list(formula, replacer))
  return(result)
}


res <- unknown.function(formula = my.formula,
                        to.replace = var.to.replace,
                        replace.with = new.terms)

res