
## ----simple-efficiency
n <- 100
eff <- function(n1) 4 * n1 * (n - n1) / n^2
curve(eff, from = 0, to = n, ylab = "Eff", xlab = expression(n[1]))


