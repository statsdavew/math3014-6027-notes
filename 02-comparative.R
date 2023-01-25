## ----pulp-expt-data
pulp <- data.frame(operator = rep(factor(1:4), 5),
                   repetition = rep(1:5, rep(4, 5)),
                   reflectance = c(59.8, 59.8, 60.7, 61.0, 60.0, 60.2, 60.7, 60.8,
                                    60.8, 60.4, 60.5, 60.6, 60.8, 59.9, 60.9, 60.5, 59.8, 60.0, 60.3, 60.5)
                     )
knitr::kable(
 tidyr::pivot_wider(pulp, names_from = operator, values_from = reflectance)[, -1],
 col.names = paste("Operator", 1:4),
 caption = "Pulp experiment: reflectance values (unitless) from four different operators."
)


## ----pulp-boxplot
boxplot(reflectance ~ operator, data = pulp)


## ----one-way-analysis
knitr::kable(
 tidyr::pivot_wider(pulp, names_from = operator, values_from = reflectance)[, -1],
 col.names = paste("Operator", 1:4),
 caption = "Pulp experiment: reflectance values (unitless) from four different operators."
)


## ----one-way-t-test
with(pulp,
 pairwise.t.test(reflectance, operator, p.adjust.method = 'none'))


## ----one-way-emmeans
pulp.lm <- lm(reflectance ~ operator, data = pulp)
anova(pulp.lm)
pulp.emm <- emmeans::emmeans(pulp.lm, ~ operator)
pairs(pulp.emm, adjust = 'none')


## ----one-way-contrasts
# same as `pairs` above
emmeans::contrast(pulp.emm, "pairwise", adjust = "none")
# estimating single contrast c = (1, -.5, -.5)
# comparing operator 1 with operators 2 and 3
contrast1v23.emmc <- function(levs)
  data.frame('t1 v avg t2 t3' = c(1, -.5, -.5, 0))
emmeans::contrast(pulp.emm, 'contrast1v23')


## ----type-I
alpha <- 0.05
1 - (1 - alpha)^6


## ----type-I-jb
alpha <- 0.05
1 - (1 - alpha)^20


## ----bonferroni
pairs(pulp.emm, adjust = 'bonferroni')


## ----tukey
pairs(pulp.emm)


## ----bonf-tukey
pairs.u <- pairs(pulp.emm, adjust = 'none')
pairs.b <- pairs(pulp.emm, adjust = 'bonferroni')
pairs.t <- pairs(pulp.emm)
data.frame(transform(pairs.b)[, 1:5], Bonf.p.value = transform(pairs.b)[, 6], Tukey.p.value = transform(pairs.t)[, 6], unadjust.p.value = transform(pairs.u)[, 6])


## ----opt-ni
opt_ni <- function(C, n) {
  CtC <- t(C) %*% C
  n * sqrt(diag(CtC)) / sum(sqrt(diag(CtC)))
}


## ----opt-pairwise
C <- matrix(
  c(
-1, 1, 0, 0,
-1, 0, 1, 0,
-1, 0, 0, 1,
0, -1, 1, 0,
0, -1, 0, 1,
0, 0, -1, 1),
  nrow = 6, byrow = T
)
opt_ni(C, 20)


## ----opt-op4
C <- matrix(
  c(1/3, 1/3, 1/3, -1),
  nrow = 1
)
opt_ni(C, 20)


## ----compare-allocations
crd_var <- function(C, n) {
  CtC <- t(C) %*% C
  sum(diag(CtC) / n)
}
n_equal <- rep(5, 4)
n_opt <- c(4, 3, 3, 10)
crd_var(C, n_opt) / crd_var(C, n_equal)


## ----exp-size
opt_n <- function(cv, prop, snr, target) target ^ 2 * c(t(cv) %*% diag( 1 / prop) %*% cv) / snr ^ 2
cv <- c(-1, 1, 0, 0)
w <- rep(1/4, 4)
snr <- c(0.5, 1, 1.5, 2, 2.5, 3)
cbind('Signal-to-noise' = snr, 'n' = opt_n(cv, w, snr, 3))


