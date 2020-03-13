remove(list = ls())
mortgage <- read.csv("mortgage.csv", header = T)
library("MASS")

# **************************************************************************
# ASSET CORRELATIONS AND WORST CASE DEFAULT RATE
# **************************************************************************

pd        <- seq(0, .1, by = 0.001)
AssetCorr <- function(x, a, b, c)
{
  a * (1-exp(-b* x))/ (1-exp(b)) + c *(1-(1-exp(-b*pd)/(1-exp(-b))))
}
ac_corporate <- AssetCorr(pd,0.12,50,.24) 
ac_retail    <- AssetCorr(pd,0.13,35,.16) 
ac_mortgages <- rep(.15, length(pd)) 
ac_revolving <- rep(.04, length(pd))

par(xpd = T)
par(mar = c(6,4,4,2) +0.1)
par(mgp = c(1.8,0.8,0))
plot(0, xlim = c(0,0.1), ylim = c(0, 0.25), xlab = "PD", ylab = "asset correlation")
lines(smooth.spline(pd, ac_corporate), lty=1)
lines(smooth.spline(pd, ac_mortgages), lty=2)
lines(smooth.spline(pd, ac_revolving), lty=3)
lines(smooth.spline(pd, ac_retail),    lty=4)
legend("bottom", inset = -.4, legend = 
         c("ac_corporate","ac_mortgages","ac_revolving","ac_retail" ),
       lty = 1:4, ncol = 2, bty = "n",horiz = T, cex = 1)

wcdr_corporate <- pnorm((qnorm(pd) + sqrt(ac_corporate))/ sqrt(1-ac_corporate))
wcdr_retail    <- pnorm((qnorm(pd) + sqrt(ac_retail))   / sqrt(1-ac_retail))
wcdr_mortgages <- pnorm((qnorm(pd) + sqrt(ac_mortgages))/ sqrt(1-ac_mortgages))
wcdr_revolving <- pnorm((qnorm(pd) + sqrt(ac_revolving))/ sqrt(1-ac_revolving ))

par(xpd = T)
par(mar = c(6,4,4,2) + 0.1)
par(mgp = c(1.8,0.8,0))
plot(0, xlim = c(0,0.1), ylim = c(0, 0.20), xlab = "PD", ylab = "WCDR")
lines(smooth.spline(pd, wcdr_corporate), lty = 1)
lines(smooth.spline(pd, wcdr_mortgages), lty = 2)
lines(smooth.spline(pd, wcdr_revolving), lty = 3)
lines(smooth.spline(pd, wcdr_retail),    lty = 4)
legend("bottom", inset = -.4,
       legend = c("wcdr_corporate","wcdr_mortgages","wcdr_revolving","wcdr_retail"),
       lty = 1:4, ncol = 4, bty = "n", cex = 1)

DEAD <- 1
DLGD <- 1
MA   <- 1

capital_corporate <- (wcdr_corporate - pd)* DEAD * DLGD * MA
capital_mortgages <- (wcdr_mortgages - pd)* DEAD * DLGD * MA
capital_revolving <- (wcdr_revolving - pd)* DEAD * DLGD * MA
capital_retail    <- (wcdr_retail -    pd)* DEAD * DLGD * MA

par(xpd = TRUE)
par(mar = c(6,4,4,2) + 0.1)
par(mgp = c(1.8,0.8,0))
plot(0, xlim = c(0,0.1), ylim = c(0, 0.08), xlab = "PD", ylab = "Capital")
lines(smooth.spline(pd, capital_corporate), lty = 1)
lines(smooth.spline(pd, capital_mortgages), lty = 2)
lines(smooth.spline(pd, capital_revolving), lty = 3)
lines(smooth.spline(pd, capital_retail),    lty = 4)
legend("bottom", inset = -.4,
       legend = c("capital_corporate","capital_mortgages","capital_revolving","capital_retail"),
       lty = 1:4, ncol = 2, bty = "n", cex = 0.75)


# **************************************************************************
# STRESS TESTING APPLICATIONS IN R
# **************************************************************************

# *****************************************************
# Scenario Based Stress Testing
# *****************************************************
model_probit <- glm(default_time ~ FICO_orig_time + 
                    LTV_time + gdp_time + uer_time, 
                    family = binomial(link = "probit"),
                    data = mortgage)
summary(model_probit)

# use ony those loans with a time = 60
mortgage_base <-mortgage[mortgage$time == 60,
                         c("FICO_orig_time" , "LTV_time","gdp_time", "uer_time")]
predicted    <- predict(model_probit, mortgage_base)

PD_time_base <- pnorm(predicted)
PD_time_wcdr <- pnorm((qnorm(PD_time_base)+sqrt(0.15)) / sqrt(1-.15))
PD_time_base_pcs <- ecdf(PD_time_base)
PD_time_wcdr_pcs <- ecdf(PD_time_wcdr)

# GDP is = 2.836358
min(mortgage_base$gdp_time)
max(mortgage_base$gdp_time)
# UERate is = 5.7
min(mortgage_base$uer_time)
max(mortgage_base$uer_time)

gdp_time_stress <- rep(-6.1, dim(mortgage_base)[1])
uer_time_stress <- rep(8, dim(mortgage_base)[1])
mortgage_stress <- mortgage_base
mortgage_stress[, c("gdp_time", "uer_time")] <- cbind(gdp_time_stress, uer_time_stress)

PD_time_stress     <- predict(model_probit, mortgage_stress, type = "response")
PD_time_stress_pcs <- ecdf(PD_time_stress) 
plot(0, xlim = c(0,.1), ylim = c(0,1), xlab = "PD", ylab = "Frequency")
lines(PD_time_base_pcs, col = "black")
lines(PD_time_wcdr_pcs, col = "lightgrey")
lines(PD_time_stress_pcs, col = "darkgrey")
legend("bottomright", 
       legend = c("PD_time_base_pcs","PD_time_wcdr_pcs","PD_time_stress_pcs"),
       lty = rep(1,3), col = c("black","lightgrey","darkgrey"), 
       bty = "n", cex = 0.75, horiz = F)


# **************************************************************************
# STRESS TESTING AND PARAMETER UNCERTAINTY
# **************************************************************************

# *****************************************************
# Basic Stress Testing of Parameter Uncertainty 
# *****************************************************

model_probit2 <- glm(default_time ~ FICO_orig_time + LTV_time
                    + gdp_time, 
                    family = binomial(link = "probit"),
                    data = mortgage)
summary(model_probit2)

covb <- vcov(model_probit2)
print(covb)

mortgage_base2 <- mortgage_base[, -4]
mortgage_base2 <- mortgage_base2[complete.cases(mortgage_base2),]
PD_time_base2  <- predict(model_probit2, mortgage_base2, type = "response")
PD_time_base2_pcs <- ecdf(PD_time_base2)
estimate <- model_probit2$coefficients
std_err  <- summary(model_probit2)$coefficients[,2]

estimate_stress    <- estimate    + std_err*.99
estimate_stress[4] <- estimate[4] - std_err[4]*.99


PD_time_stress2 <- as.numeric(pnorm(
  estimate_stress[1] + as.matrix(mortgage_base2)%*% estimate_stress[2:4]))
PD_time_stress2_pcs <- ecdf(PD_time_stress2)

plot(0, xlim = c(0,.1), ylim = c(0,1), xlab = "PD", ylab = "Frequency")
lines(PD_time_base2_pcs, col = "black")
lines(PD_time_stress2_pcs, col = "lightgrey")
legend("bottomright", 
       legend = c("PD_time_base2_pcs","PD_time_stress2_pcs"),
       lty = rep(1,2), col = c("black","lightgrey"), 
       bty = "n", cex = 0.75, horiz = F)


# Multivariate Stress Testing of Parameter Uncertainty

library("MASS")
param_sim <- mvrnorm(1000, estimate, covb)
PD_time_stress3 <- pnorm(param_sim[, 1] + as.matrix(mortgage_base2)%*%t(param_sim[,2:4]))

PD_time_stress4 <- apply(PD_time_stress3, 1, quantile, probs = 0.99)
PD_time_stress4_pcs <- ecdf(PD_time_stress4)

plot(0, xlim = c(0,.1), ylim = c(0,1), xlab = "PD", ylab = "Frequency")
lines(PD_time_base2_pcs, col = "black")
lines(PD_time_stress2_pcs, col = "lightgrey")
lines(PD_time_stress4_pcs, col = "darkgrey")
legend("bottomright", 
       legend = c("PD_time_base2_pcs","PD_time_stress2_pcs","PD_time_stress4_pcs"),
       lty = rep(1,3), col = c("black","lightgrey", "darkgrey"), 
       bty = "n", cex = 0.75, horiz = F)

