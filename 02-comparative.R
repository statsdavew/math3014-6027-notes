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

boxplot(reflectance ~ operator, data = pulp)

knitr::kable(
 tidyr::pivot_wider(pulp, names_from = operator, values_from = reflectance)[, -1],
 col.names = paste("Operator", 1:4),
 caption = "Pulp experiment: reflectance values (unitless) from four different operators."
)

with(pulp,
 pairwise.t.test(reflectance, operator, p.adjust.method = 'none'))

pulp.lm <- lm(reflectance ~ operator, data = pulp)
anova(pulp.lm)
pulp.emm <- emmeans::emmeans(pulp.lm, ~ operator)
pairs(pulp.emm, adjust = 'none')

# same as `pairs` above
emmeans::contrast(pulp.emm, "pairwise", adjust = "none")
# estimating single contrast c = (1, -.5, -.5)
# comparing operator 1 with operators 2 and 3
contrast1v23.emmc <- function(levs)
  data.frame('t1 v avg t2 t3' = c(1, -.5, -.5, 0))
emmeans::contrast(pulp.emm, 'contrast1v23')



pairs(pulp.emm, adjust = 'bonferroni')

pairs(pulp.emm)

pairs.b <- pairs(pulp.emm, adjust = 'bonferroni')
pairs.t <- pairs(pulp.emm)
pairs.u <- pairs(pulp.emm, adjust = "none")

data.frame(transform(pairs.b)[, 1:5], Bonf.p.value = transform(pairs.b)[, 6],
           Tukey.p.value = transform(pairs.t)[, 6], unadjust.p.value = transform(pairs.u)[, 6])


