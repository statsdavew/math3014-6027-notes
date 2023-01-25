## ----example-factorial
example.design <- FrF2::FrF2(nruns = 8, nfactors = 3, randomize = F)
knitr::kable(example.design, caption = "Treatments from a $2^3$ factorial design", align = rep("r", 3))


## ----block-frac-1
block1 <- c(1, 2, 2, 1, 2, 1, 1, 2)
example.design.a <- cbind(example.design, Block = block1)
knitr::kable(example.design.a, caption = "Treatments from a $2^3$ factorial design split into two blocks of size four.", align = rep("r", 4))


## ----block-frac-1-mm
X <- model.matrix( ~ Block + (A + B + C)^3, data = example.design.a)
Xdf <- data.frame(X[, -1])
colnames(Xdf) <- c("Block", "A", "B", "C", "A:B", "A:C", "B:C", "A:B:C")
Xdf <- dplyr::mutate(Xdf, Treatment = 1:8, .before = Block)
knitr::kable(Xdf, caption = "Unscaled factorial effect contrasts for a $2^3$ design with one possible assignment of treatments to blocks")


## ----block-assignment
no.assign <- choose(8, 4)
assignments <- combinat::combn(8, 4)
yfake <- rnorm(8)
Xadf <- cbind(Xdf[, c(-1, -9)], y = yfake)
avgvar <- NULL
for(i in 1:no.assign) {
  B <- rep(1, 8)
  B[assignments[, i]] <- -1
  Xadf$Block <- B
  temp.lm <- lm(y ~ Block + (A + B + C)^2, data = Xadf)
  temp.lm$residuals <- yfake
  temp.lm$df.residual <- 8
  vmat <- vcov(temp.lm) / (summary(temp.lm)$sigma^2)
  vars <- (diag(vmat[-c(1:2), -c(1:2)]))
  tidyr::replace_na(vars, Inf)
  avgvar[i] <- sum(tidyr::replace_na(vars, Inf)) / 6
}
knitr::kable(table(avgvar), col.names = c("Avg. variance", "Freq."))


## ----partial-confounding
Xa <- as.matrix(Xdf[, -c(1)])
Xb <- Xa
B <- rep(1, 8)
B[assignments[, 2]] <- -1
Xb[,1] <- B


## ----cor-optimal
knitr::kable(cor(Xa), caption = "Scaled inner-products between contrast vectors for $2^3$ with treatments assigned to blocks so $\\mathrm{Blocks} = ABC$.")


## ----cor-partial
knitr::kable(cor(Xb), caption = "Scaled inner-products between contrast vectors for $2^3$ with treatments assigned to blocks arbitrarily.")


## ----block-frac-2
block2 <- c(4, 3, 2, 1, 1, 2, 3, 4)
example.design.b <- cbind(example.design, Block = block2)
knitr::kable(example.design.b, caption = "Treatments from a $2^3$ factorial design split into four blocks of size two.", align = rep("r", 4))


## ----block-frac-2-mm
X <- model.matrix( ~ Block + (A + B + C)^3, data = example.design.b)
Xdf <- data.frame(X[, -1])
colnames(Xdf) <- c("Block", "A", "B", "C", "A:B", "A:C", "B:C", "A:B:C")
Xdf <- dplyr::mutate(Xdf, Treatment = 1:8, .before = Block)
knitr::kable(Xdf, caption = "Unscaled factorial effect contrasts for a $2^3$ design with one possible assignment of treatments to four blocks of size two.")


## ----frf2-block-1
block1.frf2 <- FrF2::FrF2(nruns = 8, nfactors = 3, blocks = 2,
                          alias.info = 3, randomize = F)
block1 <- data.frame(model.matrix(~ Blocks + (A + B + C)^3, block1.frf2))
block1 <- dplyr::mutate(block1, Treatment = 1:8, .before = Blocks1)
knitr::kable(block1[, -1], col.names = c("Treatment", "Block", "A", "B", "C",
                                         "A:B", "A:C", "B:C", "A:B:C"),
             caption = "$2^3$ factorial design in two blocks of size four")


## ----frf2-block-2
block2.frf2 <- FrF2::FrF2(nruns = 8, nfactors = 3, blocks = 4,
                          alias.info = 3, randomize = F, alias.block.2fis = T)
block2 <- data.frame(model.matrix(~ Blocks + (A + B + C)^3, block2.frf2))
block2 <- dplyr::mutate(block2, Treatment = 1:8, .before = Blocks1)
knitr::kable(block2[, -1],
             col.names = c("Treatment", "Block1", "Block2", "Block3",
                           "A", "B", "C",
                           "A:B", "A:C", "B:C", "A:B:C"),
             caption = "$2^3$ factorial design in four blocks of size two")



## ----FrF2-info
library(FrF2)
design.info(block1.frf2)$aliased.with.blocks
design.info(block2.frf2)$aliased.with.blocks


## ----FrF2-choose-confounding
block3.frf2 <- FrF2::FrF2(nruns = 2^8, nfactors = 8,
                     alias.info = 3, randomize = F, blocks = c("ACEGH", "BCFGH", "BDEGH"))


## ----confounded-analysis
ex2.df <- dplyr::mutate_at(Xdf,c("Block"), factor)
X <- model.matrix(~ Block + (A + B + C)^3, ex2.df)
coef.values <- c(0, 2, 6, 12, 4, 0, 3, 0, 2, 0, 0)
y <- X %*% coef.values + rnorm(8, 0, 1)
betahat <- t(X[, c(5, 6, 7, 11)]) %*% y / nrow(X)
betahat # coef estimates, obtained directly
ex2.df <- data.frame(ex2.df, y = y)
ex2.lm <- lm(y ~  Block + A + B + C + A:B:C, data = ex2.df)
anova(ex2.lm)
coef(ex2.lm)[5:8] # coef estimates
2 * coef(ex2.lm)[5:8] # factorial effects


## ----fact-lm2
ex2.lm2 <- lm(y ~  (A + B + C)^3, data = ex2.df)
anova(ex2.lm2)
coef(ex2.lm2)


