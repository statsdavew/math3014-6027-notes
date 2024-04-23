## ----size-of-factorial
size <- data.frame(1:15, 2^(1:15))
knitr::kable(size, col.names = c("No. factors", "No. of trts"), caption = "Number of treatments in a $2^f$ factorial designs for different numbers, $f$, of factors.")


## ----spring-factors
factor.name <- c("Quench temperature (F)", "Heat temperature (F)", "Heating time (s)",
                 "Transfer time (s)", "Hold-down time (s)")
low.level <- c("130-150", 1840, 23, 10, 2)
high.level <- c("150-170", 1880, 25, 12, 3)
spring.factors <- data.frame(factor = factor.name, low = low.level, high = high.level)
row.names(spring.factors) <- LETTERS[1:5]
knitr::kable(spring.factors, col.names = c("Factor", "Low level", "High level"),
             align = c("l", "r", "r"), caption = "Spring experiment: factors and levels")


## ----spring-data
spring <- FrF2::FrF2(nruns = 16, nfactors = 5, generators = "BCD", randomize = F)
spring$height <- c(7.54, 7.20, 7.69, 7.63, 7.94, 7.40, 7.95, 7.62, 7.52, 7.52,
                   7.63, 7.65, 7.79, 7.29, 8.07, 7.73)
knitr::kable(spring, caption = "Spring experiment: 16 run design.", align = rep("r", 6))


## ----spring-def-rel
fac_to_numeric <- function(x) as.numeric(as.character(x))
BCDE <- fac_to_numeric(spring$B) * fac_to_numeric(spring$C) *
  fac_to_numeric(spring$D) * fac_to_numeric(spring$E)
BCDE


## ----spring-alias1
BC <- fac_to_numeric(spring$B) * fac_to_numeric(spring$C)
DE <- fac_to_numeric(spring$D) * fac_to_numeric(spring$E)
BC
DE
all.equal(BC, DE)

## ----spring-alias2
B <- fac_to_numeric(spring$B)
CDE <- fac_to_numeric(spring$C) * fac_to_numeric(spring$D) * fac_to_numeric(spring$E)
B
CDE
all.equal(B, CDE)


## ----frfr-aliases
spring.lm <- lm(height ~ (.)^5, data = spring)
FrF2::aliases(spring.lm)


## ----spring-alias
t(alias(spring.lm)$Complete)


## ----FrF2-example
ff.2.6.2 <- FrF2::FrF2(nruns = 16, nfactors = 6, generators = c("ABC", "BCD"),
                       randomize = F, alias.info = 3)


## ----FrFr-example-di
library(FrF2)
design.info(ff.2.6.2)$aliased


## ----FrF2-example-design
knitr::kable(ff.2.6.2, align = "r",
             caption = "Treatments from a $2^{6-2}$ fractional factorial design.")


## ----FrF2-res
ff.2.6.2.r <- FrF2::FrF2(nfactors = 6, resolution = 4,
                         randomize = F, alias.info = 3)
design.info(ff.2.6.2.r)$aliased


## ----FrF2-alone
ff.2.6.2.a <- FrF2::FrF2(nfactors = 6, nruns = 16,
                         randomize = F, alias.info = 3)
design.info(ff.2.6.2.a)$aliased


## ----spring-contrast
spring$treatment <- factor(1:16)
spring.ut <- lm(height ~ treatment, data = spring)
fac.contrasts.emmc <- function(nlevs, ...) {
  spring.num <- apply(spring[, c("A", "B", "C", "D", "E")], 2, fac_to_numeric)
  data.frame(model.matrix(~ . + A:B + A:C + A:D + A:E + B:C
                          + B:D + B:E, data.frame(spring.num))[, -1] / 8)
}
spring.emm <- emmeans::emmeans(spring.ut, ~ treatment)
emmeans::contrast(spring.emm, 'fac.contrasts')


## ----spring-reg
spring.lm <- lm(height ~ (A + B + C + D + E) ^ 5, data = spring)
c(na.omit(2 * coef(spring.lm)[-1]))


## ----spring-sigma2
spring.lm <- lm(height ~ (A + B + C + D + E) ^ 2, data = spring)
anova(spring.lm)


## ----frac-block
ff.2.6.2.b.4 <- FrF2::FrF2(nruns = 16, nfactors = 6, generators = c("ABC", "ABD"), blocks = c("ACD", "BCD"),
                           randomize = F, alias.block.2fis = T, alias.info = 3)
design.info(ff.2.6.2.b.4)$aliased
design.info(ff.2.6.2.b.4)$aliased.with.blocks


## ----frac-block-design
block12 <- ff.2.6.2.b.4[1:8, ]
block34 <- ff.2.6.2.b.4[9:16, ]
knitr::kable(cbind(block12, block34),
      caption = "Fractional factorial
      $2^{6-2}$ design in $b=4$ blocks of size $k=4$.",
      align = "r")


## ----spring-blocks
spring.blocks <- spring
spring.blocks$blocks <- with(data.frame(spring), fac_to_numeric(A) * fac_to_numeric(B) * fac_to_numeric(C))
springb.lm <- lm(height ~ blocks + (A + B + C + D + E) ^ 3, data = spring.blocks)
anova(springb.lm)

