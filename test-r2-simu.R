library(bigstatsr)
G <- big_attachExtdata()
r <- 0.6
a <- sqrt(1 / r^2 - 1)
g <- G[, 1]
x <- g + rnorm(length(g), sd = a * sd(g))
plot(x, g)
cor(x, g)

r2 <- sqrt(runif(ncol(G), 0.1, 1))
r2_simu <- sapply(cols_along(G), function(j) {
  a <- sqrt(1 / r2[j] - 1)
  g <- G[, j]
  x <- g + rnorm(length(g), sd = a * sd(g))
  x2 <- pmax(0, pmin(x, 2))
  cor(x2, g)^2
})
plot(r2, r2_simu, pch = 20); abline(0, 1, col = "red", lwd = 2)
