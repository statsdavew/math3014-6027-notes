## ----factorial-desilylation-data
desilylation <- FrF2::FrF2(nruns = 16, nfactors = 4, randomize = F,
                           factor.names = list(temp = c(10, 20), time = c(19, 25),
                                               solvent = c(5, 7), reagent = c(1, 1.33)))
yield <- c(82.93, 94.04, 88.07, 93.97, 77.21, 92.99, 83.60, 94.38,
           88.68, 94.30, 93.00, 93.42, 84.86, 94.26, 88.71, 94.66)
desilylation <- data.frame(desilylation, yield = yield)
rownames(desilylation) <- paste("Trt", 1:16)
knitr::kable(desilylation,
             col.names = c("Temp (degrees C)", "Time (hours)", "Solvent (vol.)",
                           "Reagent (equiv.)", "Yield (%)"),
             caption = "Desilylation experiment: 16 treatments defined
             by settings of four factors, with response (yield).")


## ----desilylation-anova
desilylation <- data.frame(desilylation, trt = factor(1:16))
desilylation.lm <- lm(yield ~ trt, data = desilylation)
anova(desilylation.lm)


## ----desilylation-lm
desilylation.emm <- emmeans::emmeans(desilylation.lm, ~ trt)
reagent_me.emmc <- function(levs, ...) data.frame('reagent m.e.' = rep(c(-1, 1), rep(8, 2)) / 8)
emmeans::contrast(desilylation.emm, 'reagent_me')


## ----desilylation-all-me-contrasts
contrast.mat <- FrF2::FrF2(nruns = 16, nfactors = 4, randomize = F,
                           factor.names = c("temp", "time", "solvent", "reagent"))
fac.contrasts.emmc <- function(levs, ...)
  dplyr::mutate_all(data.frame(contrast.mat), function(x) scale(as.numeric(as.character(x)), scale = 8))
main_effect_contrasts <- fac.contrasts.emmc()
rownames(main_effect_contrasts) <- paste("Trt", 1:16)
knitr::kable(main_effect_contrasts, caption = 'Desilylation experiment: main effect contrast coefficients', col.names = c("Temperature", "Time", "Solvent", "Reagent"))


## ----desilylation-all-me-estimates
t(as.matrix(main_effect_contrasts)) %*% yield
emmeans::contrast(desilylation.emm, 'fac.contrasts')


## ----desilylation-me-plots
## calculate the means
temp_bar <- aggregate(yield ~ temp, data = desilylation, FUN = mean)
time_bar <- aggregate(yield ~ time, data = desilylation, FUN = mean)
solvent_bar <- aggregate(yield ~ solvent, data = desilylation, FUN = mean)
reagent_bar <- aggregate(yield ~ reagent, data = desilylation, FUN = mean)

## convert factors to numeric
fac_to_num <- function(x) as.numeric(as.character(x))
temp_bar$temp <- fac_to_num(temp_bar$temp)
time_bar$time <- fac_to_num(time_bar$time)
solvent_bar$solvent <- fac_to_num(solvent_bar$solvent)
reagent_bar$reagent <- fac_to_num(reagent_bar$reagent)

## main effect plots
plotmin <- min(temp_bar$yield, time_bar$yield, solvent_bar$yield, reagent_bar$yield)
plotmax <- max(temp_bar$yield, time_bar$yield, solvent_bar$yield, reagent_bar$yield)
par(cex = 2)
layout(matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE), respect = T)
plot(temp_bar, pch = 16, type = "b", ylim = c(plotmin, plotmax))
plot(time_bar, pch = 16, type = "b", ylim = c(plotmin, plotmax))
plot(solvent_bar, pch = 16, type = "b", ylim = c(plotmin, plotmax))
plot(reagent_bar, pch = 16, type = "b", ylim = c(plotmin, plotmax))


## ----desilylation-solvent-reagent-interaction
sol_reg_int.emmc <- function(levels, ...)
  data.frame('reagent x solvent' = .125 * c(rep(1, 4), rep(-1, 8), rep(1, 4)))
emmeans::contrast(desilylation.emm, 'sol_reg_int')


## ----desilylation-2fi-interactions
fac.contrasts.int.emmc <- function(levs, ...) {
  with(sqrt(8) * main_effect_contrasts, {
    data.frame('tem_x_tim' = temp * time,
               'tem_x_sol' = temp * solvent,
               'tem_x_rea' = temp * reagent,
               'tim_x_sol' = time * solvent,
               'tim_x_rea' = time * reagent,
               'sol_x_rea' = solvent * reagent)
  })
}
two_fi_contrasts <- fac.contrasts.int.emmc()
rownames(two_fi_contrasts) <- paste("Trt", 1:16)
knitr::kable(two_fi_contrasts, caption = 'Desilylation experiment: two-factor interaction contrast coefficients')


## ----desilylation-all-2fi-estimates
t(as.matrix(two_fi_contrasts)) %*% yield
emmeans::contrast(desilylation.emm, 'fac.contrasts.int')


## ----desilylation-2fi-plots
plotmin <- min(desilylation$yield)
plotmax <- max(desilylation$yield)
par(cex = 2)
layout(matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2, byrow = TRUE), respect = T)

with(desilylation, {
  interaction.plot(temp, time, yield, type = "b", pch = 16, legend = F,
                   ylim = c(plotmin, plotmax))
  legend("bottomright", legend = c("Time = 19", "Time = 25"), lty = 2:1, cex = .75)
  interaction.plot(temp, solvent, yield, type = "b", pch = 16, legend = F,
                   ylim = c(plotmin, plotmax))
  legend("bottomright", legend = c("Solvent = 5", "Solvent = 7"), lty = 2:1, cex = .75)
  interaction.plot(temp, reagent, yield, type = "b", pch = 16, legend = F,
                   ylim = c(plotmin, plotmax))
  legend("bottomright", legend = c("Reagent = 1", "Reagent = 1.33"), lty = 2:1, cex = .75)
  interaction.plot(time, solvent, yield, type = "b", pch = 16, legend = F,
                   ylim = c(plotmin, plotmax))
  legend("bottomright", legend = c("Solvent = 5", "Solvent = 7"), lty = 2:1, cex = .75)
  interaction.plot(time, reagent, yield, type = "b", pch = 16, legend = F,
                   ylim = c(plotmin, plotmax))
  legend("bottomright", legend = c("Reagent = 1", "Reagent = 1.33"), lty = 2:1, cex = .75)
  interaction.plot(solvent, reagent, yield, type = "b", pch = 16, legend = F,
                   ylim = c(plotmin, plotmax))
  legend("bottomright", legend = c("Reagent = 1", "Reagent = 1.33"), lty = 2:1, cex = .75)
})


## ----desilylation-fe
## hadamard products
unscaled_me_contrasts <- 8 * main_effect_contrasts
factorial_contrasts <- model.matrix(~.^4, unscaled_me_contrasts)[, -1] / 8

## all factorial effects - directly, as there is no treatment replication
t(factorial_contrasts) %*% yield

## using emmeans
factorial_contrasts.emmc <- function(levs, ...) data.frame(factorial_contrasts)
desilylation.effs <- emmeans::contrast(desilylation.emm, 'factorial_contrasts')
desilylation.effs


## ----desilylation-effects-plot1
effs <- dplyr::arrange(transform(desilylation.effs)[,1:2], dplyr::desc(estimate))
knitr::kable(effs, caption = "Desilylation experiment: sorted estimated factorial effects")
qqnorm(effs[ ,2], ylab = "Factorial effects", main = "") # note that qqnorm/qqline/qqplot don't require sorted data
qqline(effs[ ,2])


## ----desilylation-half-normal-plot
p <- .5 + .5 * (1:16 - .5) /16 # probabilities we will plot against
qqplot(x = qnorm(p), y = abs(effs[,2]), ylab = "Absolute factorial effects",
       xlab = "Half-normal quantiles")


## ----desilylation-pse
s0 <- 1.5 * median(abs(effs[, 2]))
trimmed <- abs(effs[, 2]) < 2.5 * s0
pse <- 1.5 * median(abs(effs[trimmed, 2]))
pse


## ----desilylation-lenth
eff_est <- effs[, 2]
names(eff_est) <- effs[, 1]
lenth_tests <- unrepx::eff.test(eff_est, method = "Lenth")
knitr::kable(lenth_tests, caption = "Desilylation experiment: hypothesis tests using Lenth's method.")


## ----desilylation-lenth-plot
unrepx::hnplot(eff_est, method = "Lenth", horiz = F, ID = 2.7, alpha = 0.05)


## ----factorial-regression-c-matrix
y <- kronecker(desilylation$yield, rep(1, 3)) # hypothetical response vector
C <- factorial_contrasts
Ctilde <- kronecker(C, rep(1, 3))
t(Ctilde) %*% y / 3 # to check


## ----desilylation-X-C
X <- 8 * C
Xt <- t(X)
XtX <- Xt %*% X
2 * solve(XtX) %*% t(X) %*% desilylation$yield


## ----desilylation-reg
desilylation.df <- dplyr::mutate(desilylation,
                                      across(.cols = temp:reagent,
                                             ~ as.numeric(as.character(.x))))
desilylation.df <- dplyr::select(desilylation.df, -c(trt))
desilylation.df <- dplyr::mutate(desilylation.df, across(.cols = temp:reagent,
                                               ~ scales::rescale(.x, to = c(-1, 1))))
desilylation_reg.lm <- lm(yield ~ (.) ^ 4, data = desilylation.df)
knitr::kable(2 * coef(desilylation_reg.lm)[-1], caption = "Desilylation example: factorial effects calculated using a regression model.")


## ----desilylation-effects
temp_x_time <- effects::Effect(c("temp", "time"), desilylation_reg.lm, xlevels = list(time = c(-1, 1)), se = F)
plot(temp_x_time, main = "", rug = F, x.var = "temp", ylim = c(80, 100))


## ----desilylation-regss
ss <- nrow(desilylation) * coef(desilylation_reg.lm)^2
regss <- sum(ss) - nrow(desilylation) * mean(desilylation$yield)^2
names(regss) <- "Regression"
knitr::kable(c(regss, ss[-1]), col.names = "Sum Sq.", caption = "Desilylation experiment: regression sums of squares for each factorial effect calculated directly.")

## ----desilylation-anova-ss
knitr::kable(anova(desilylation_reg.lm)[, 1:2], caption = "Desilylation experiment: ANOVA table from `anova` function.")


