sapply(setNames(nm = c(0.05, 0.1, 0.2)), function(K) {
  sapply(setNames(nm = c(0.3, 0.5, 0.6, 0.8)), function(h2) {
    T0 <- qnorm(1 - K)
    z <- dnorm(T0)
    i <- z / K
    v <- -i * K / (1 - K)
    num <- (i - v) * h2
    deno.part <- 1 - h2 * i * (i - T0) + (1 - h2 * v * (v - T0))
    round(100 * pnorm(num / sqrt(h2 * deno.part)), 2)
  })
})
