library(vqtl)
library(MASS)

n <- 1000
x <- runif(n = n)
y <- exp(rnorm(n = n))    # right answer is 0
y2 <- 10 + rnorm(n = n)   # right answer is 1
y3 <-  exp(rnorm(n = n, sd = 2*x))   # right answer is 0
y4 <- log(rnorm(n = n, mean = 10))

a <- boxcox(object = y ~ x)

b <- BoxCoxDGLM(mformula = y ~ x, vformula = ~1, data = data.frame(y, x))

c <- BoxCoxDGLM(mformula = y ~ x, vformula = ~x, data = data.frame(y, x))


d <- boxcox(object = y2 ~ x)

e <- BoxCoxDGLM(mformula = y2 ~ x, vformula = ~1, data = data.frame(y2, x))

f <- BoxCoxDGLM(mformula = y2 ~ x, vformula = ~x, data = data.frame(y2, x))



g <- boxcox(object = y3 ~ x)

h <- BoxCoxDGLM(mformula = y3 ~ x, vformula = ~1, data = data.frame(y3, x))

i <- BoxCoxDGLM(mformula = y3 ~ x, vformula = ~x, data = data.frame(y3, x))


j <- boxcox(object = y4 ~ x)

k <- BoxCoxDGLM(mformula = y4 ~ x, vformula = ~1, data = data.frame(y4, x))

l <- BoxCoxDGLM(mformula = y4 ~ x, vformula = ~x, data = data.frame(y4, x))
