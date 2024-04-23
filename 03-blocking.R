## ----bar-expt-data
bar <- data.frame(coating = rep(factor(1:4), 8),
                   block = rep(factor(1:8), rep(4, 8)),
                   strength = c(136, 147, 138, 149, 136, 143, 122, 153, 150, 142, 131, 136,
                                   155, 148, 130, 129, 145, 149, 136, 139, 150, 149, 147, 144,
                                   147, 150, 125, 140, 148, 149, 118, 145)
                     )
knitr::kable(
 tidyr::pivot_wider(bar, names_from = coating, values_from = strength),
 col.names = c("Block", paste("Coating", 1:4)),
 caption = "Steel bar experiment: tensile strength values (kliograms per square inch, ksi) from steel bars with four different coatings."
)


## ----bar-expt-boxplots
boxplot(strength ~ block, data = bar)
boxplot(strength ~ coating, data = bar)


## ----tyre-expt-data
tyre <- data.frame(compound = as.factor(c(1, 2, 3, 1, 2, 4, 1, 3, 4, 2, 3, 4)),
                   block = rep(factor(1:4), rep(3, 4)),
                   wear = c(238, 238, 279, 196, 213, 308, 254, 334, 367, 312, 421, 412)
                     )
options(knitr.kable.NA = '')
knitr::kable(
 tidyr::pivot_wider(tyre, names_from = compound, values_from = wear),
 col.names = c("Block", paste("Compound", 1:4)),
 caption = "Tyre experiment: relative wear measurements (unitless) from tires made with four different rubber compounds."
)


## ----tyre-expt-boxplots
boxplot(wear ~ block, data = tyre)
boxplot(wear ~ compound, data = tyre)



## ----block-bars-anova
bar.lm <- lm(strength ~ block + coating, data = bar)
anova(bar.lm)


## ----block-tyres-anova
tyre.lm <- lm(wear ~ block + compound, data = tyre)
anova(tyre.lm)


## ----blocks-s2
bar.s2 <- summary(bar.lm)$sigma^2
tyre.s2 <- summary(tyre.lm)$sigma^2


## ----blocks-bar-contrasts
bar.emm <- emmeans::emmeans(bar.lm, ~ coating)
contrastv1.emmc <- function(levs, ...)
  data.frame('t1 v t2' = c(1, -1, 0, 0), 't1 v t3' = c(1, 0, -1, 0),
  't1 v t4' = c(1, 0, 0, -1))
emmeans::contrast(bar.emm, 'contrastv1')


## ----blocks-bar-bonf
temp <- pmin(3 * transform(emmeans::contrast(bar.emm, 'contrastv1'))[, 6], 1)


## ----blocks-crd
bar_crd.lm <- lm(strength ~ coating, data = bar)
bar_crd.emm <- emmeans::emmeans(bar_crd.lm, ~ coating)
emmeans::contrast(bar_crd.emm, 'contrastv1')


## ----blocks-crd-s2
crd.s2 <- summary(bar_crd.lm)$sigma^2
rcbd.s2 <- summary(bar.lm)$sigma^2


## ----block-tyre-concurrence
N <- matrix(
  c(1, 1, 1, 0,
    1, 1, 0, 1,
    1, 0, 1, 1,
    0, 1, 1, 1),
  nrow = 4, byrow = T
)
N %*% t(N)



## ----idb-tyre
tyre.bibd <- ibd::bibd(v = 4, b = 4, r = 3, k = 3, lambda = 2) # note, v is the notation for the number of treatments
tyre.bibd$N # incidence matrix


## ----idb-larger-ex
larger.bibd <- ibd::bibd(v = 8, b = 14, r = 7, k = 4, lambda = 3)
larger.bibd$N


## ----bibd-tyre-q
trtsum <- aggregate(wear ~ compound, data = tyre, FUN = sum)[, 2]
blocksum <- aggregate(wear ~ block, data = tyre, FUN = sum)[, 2]
q <- trtsum - N %*% blocksum / 3
C <- matrix(
  c(1, -1, 0, 0,
    1, 0, -1, 0,
    1, 0, 0, -1,
    0, 1, -1, 0,
    0, 1, 0, -1,
    0, 0, 1, -1),
  ncol = 4, byrow = T
)
k <- 3; lambda <- 2; t <- 4
pe <- k * C %*% q / (lambda * t) # point estimates
se <- sqrt(2 * k * tyre.s2 / (lambda * t)) # st error (same for each contrast)
t.ratio <- pe / se
p.value <- 1 - ptukey(abs(t.ratio) * sqrt(2), 4, 5)
data.frame(Pair = c('1v2', '1v3', '1v4', '2v3', '2v4', '3v4'),
  Estimate = pe, St.err = se, t.ratio = t.ratio,
           p.value = p.value, reject = p.value < 0.05)


## ----bibd-tyre-contrasts
tyre.emm <- emmeans::emmeans(tyre.lm, ~ compound)
pairs(tyre.emm)


