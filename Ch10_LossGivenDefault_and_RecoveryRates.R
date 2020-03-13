remove(list = ls())
install.packages("moments")
install.packages("betareg")
install.packages("censReg")

library("moments")
library("betareg")
library("censReg")

lgd <- read.csv("lgd.csv")

# **************************************************************************
# MARGINAL LGD MODELS
# **************************************************************************

# *****************************************************
# Descriptive Statistics
# *****************************************************

lgd_time   <- lgd$lgd_time
y_logistic <- lgd$y_logistic
y_probit   <- lgd$Y_probit
lnrr       <- lgd$lnrr
LTV        <- lgd$LTV
purpose1   <- lgd$purpose1
event      <- lgd$event

# Check the transformation and ln(RR)
lgd_time[1]

y_logistic[1]
-log(1/lgd_time[1]-1)

y_probit[1]
qnorm(lgd_time[1])

log(lgd$Recovery_rate[1])
lnrr[1]

options(scipen = 999)

fun_statistics <- function(x){
  rbind(
    N = length(x),
    Sum_Obs = sum(x),
    Mean = mean(x),
    Median = median(x),
    Variance = var(x),
    StandDev = sd(x),
    Range = diff(range(x)),
    Skewness = skewness(x),
    Kurtosis = kurtosis(x) - 3,
    Quant_100_perc = quantile(x,1),
    Quant_99.9_perc = quantile(x,.999),
    Quant_99.5_perc = quantile(x,.995),
    Quant_99_perc = quantile(x,.99),
    Quant_97.50_perc = quantile(x,.975),
    Quant_95_perc = quantile(x,.95),
    Quant_75_perc = quantile(x,.75),
    Quant_50_perc = quantile(x,.5),
    Quant_25_perc = quantile(x,.25),
    Quant_10_perc = quantile(x,.1),
    Quant_5_perc = quantile(x,.05),
    Quant_1_perc = quantile(x,.01),
    Quant_0_perc = quantile(x,0)
  )
}
statistics_lgd_time <- fun_statistics(lgd_time)
print(statistics_lgd_time)

# For the Histogram with the Density, you cannot use the following
hist(lgd_time, breaks = 20, probability = F)
hist(lgd_time, breaks = 20, probability = T)

# But do this way
histdata_lgd_time <- hist(lgd_time, breaks = 20,plot = F)
histdata_lgd_time$density = histdata_lgd_time$counts / sum(histdata_lgd_time$counts)*100
plot(histdata_lgd_time, 
     ylim = c(0, max(histdata_lgd_time$density) + 10),
     ylab = "Percent", xlab = "Lgd",freq = F)
box()

statistics_y_logistic <- fun_statistics(y_logistic)
print(statistics_y_logistic)

histdata_y_logistic <- hist(y_logistic, breaks = 20,plot = F)
histdata_y_logistic$density = histdata_y_logistic$counts / sum(histdata_y_logistic$counts)*100
plot(histdata_y_logistic, ylim = c(0, max(histdata_y_logistic$density)+2),
     ylab = "Percent", xlab = "y_logistic",freq = F)
box()

statistics_y_probit <- fun_statistics(y_probit)
print(statistics_y_probit)

histdata_y_probit <- hist(y_probit, breaks = 20,plot = F)
histdata_y_probit$density = histdata_y_probit$counts / sum(histdata_y_probit$counts)*100
plot(histdata_y_probit, ylim = c(0, max(histdata_y_probit$density)+2),
     ylab = "Percent", xlab = "y_probit",freq = F)
box()


statistics_lnrr <- fun_statistics(lnrr)
print(statistics_lnrr)

histdata_lnrr <- hist(lnrr, breaks = 20,plot = F)
histdata_lnrr$density = histdata_lnrr$counts / sum(histdata_lnrr$counts)*100
plot(histdata_lnrr, ylim = c(0, max(histdata_lnrr$density)+2),
     ylab = "Percent", xlab = "lnrr",freq = F)
box()


statistics_LTV <- fun_statistics(LTV)
print(statistics_LTV)

histdata_LTV <- hist(LTV, breaks = 20,plot = F)
histdata_LTV$density = histdata_LTV$counts / sum(histdata_LTV$counts)*100
plot(histdata_LTV, ylim = c(0, max(histdata_LTV$density)+2),
     ylab = "Percent", xlab = "LTV",freq = F)
box()

# *****************************************************
# Linear Regression
# *****************************************************
linear_model <- lm(lgd_time ~ LTV + purpose1)
summary(linear_model)
# plot(linear_model, which = 1:6)
hist(residuals(linear_model), breaks = 20, freq = F, xaxs = "i", yaxs = "i",
     ylab = "Density", xlab = "Residuals", main = "", las = 1)
curve(dnorm(x, mean = 0, sd = sd(residuals(linear_model))), add = T, col = "red")

# *****************************************************
# Transformed Linear Regression
# *****************************************************
logistic_linear_model <- lm(y_logistic ~ LTV + purpose1)
summary(logistic_linear_model)
#plot(logistic_linear_model, which = 1:6)
hist(residuals(logistic_linear_model), breaks = 20, freq = F, xaxs = "i", yaxs = "i",
     ylab = "Density", xlab = "Residuals", main = "", las = 1)
curve(dnorm(x, mean = 0, sd = sd(residuals(logistic_linear_model))), add = T, col = "red")

# *****************************************************
# probit Linear Regression
# *****************************************************
probit_linear_model <- lm(y_probit ~ LTV + purpose1)
summary(probit_linear_model)
#plot(probit_linear_model, which = 1:6)
hist(residuals(probit_linear_model), breaks = 20, freq = F, xaxs = "i", yaxs = "i",
     ylab = "Density", xlab = "Residuals", main = "", las = 1)
curve(dnorm(x, mean = 0, sd = sd(residuals(probit_linear_model))), add = T, col = "red")

# *****************************************************
# Non Linear Regression
# *****************************************************
ll <- function(parvector){
  b0 <- parvector[1]
  b1 <- parvector[2]
  b2 <- parvector[3]
  sigma <- parvector[4]
  xb <- b0 + b1 * LTV + b2 * purpose1
  mu = 1/(1+exp(-xb))
  sum(log(dnorm(lgd_time, mu, sigma)))
}
nonlinear_model <-optim(par = c(0,0,0,1), 
                        fn = ll, method = "BFGS", 
                        control = list(fnscale = -1))
estimates <- as.data.frame(nonlinear_model$par)
row.names(estimates) <- c("b0", "b1", "b2", "b3")
colnames(estimates)  <- c("Parameter Estimates")
print(estimates)
ress <- 
  ((estimates[[1]][1] + estimates[[1]][2]*LTV + estimates[[1]][3] * purpose1) - 
  lgd_time)

mean(ress)
hist(ress - mean(ress), breaks = 20, freq = F, xaxs = "i", yaxs = "i",
     ylab = "Density", xlab = "Residuals", main = "", las = 1)
curve(dnorm(x, mean = 0, sd = sd(ress)), add = T, col = "red")

# **************************************************************************
# FRACTIONAL LOGIT REGRESSION
# **************************************************************************

# *****************************************************
# Approach 1
# *****************************************************

ll_Fractional <- function(parvector){
  b0 <- parvector[1]
  b1 <- parvector[2]
  b2 <- parvector[3]
  xb <- b0 + b1 * LTV + b2 * purpose1
  mu = 1/(1+exp(-xb))
  sum(log((mu ^ lgd_time) * ((1-mu) ^ (1- lgd_time))))
}

fractional_logit_model_1 <-optim(par = c(0,0,0), 
                        fn = ll_Fractional, method = "BFGS", 
                        control = list(fnscale = -1))

estimates_fractional_1            <- as.data.frame(fractional_logit_model_1$par)
row.names(estimates_fractional_1) <- c("b0", "b1", "b2")
colnames(estimates_fractional_1)  <- c("Parameter Estimates")
print(estimates_fractional_1)

ress_2 <- 
  ((estimates_fractional_1[[1]][1] + estimates_fractional_1[[1]][2]*LTV + 
      estimates_fractional_1[[1]][3] * purpose1) - 
     lgd_time)

mean(ress_2)
hist(ress_2 - mean(ress_2), breaks = 20, freq = F, xaxs = "i", yaxs = "i",
     ylab = "Density", xlab = "Residuals", main = "", las = 1)
curve(dnorm(x, mean = 0, sd = sd(ress_2)), add = T, col = "red")

# *****************************************************
# Approach 2 - GLM Regression 
# *****************************************************

fractional_logit_model_2 <-glm(lgd_time ~ LTV + purpose1,
                               family = binomial(link = "logit"), 
                               start = c(0,0,0))
summary(fractional_logit_model_2)

estimates_fractional_2 <- as.data.frame(fractional_logit_model_2$coefficients)
row.names(estimates_fractional_2) <- c("b0 Intercept", "b1 - LTV", "b2 purpose 1")
colnames(estimates_fractional_2) <- c("Parameter Estimates")
print(estimates_fractional_2)

ress_3 <- 
  ((estimates_fractional_2[[1]][1] + estimates_fractional_2[[1]][2]*LTV + 
      estimates_fractional_2[[1]][3] * purpose1) - 
     lgd_time)

mean(ress_3)
hist(ress_3 - mean(ress_3), breaks = 20, freq = F, xaxs = "i", yaxs = "i",
     ylab = "Density", xlab = "Residuals", main = "", las = 1)
curve(dnorm(x, mean = 0, sd = sd(ress_2)), add = T, col = "red")

# **************************************************************************
# BETA REGRESSION
# **************************************************************************

# *****************************************************
# Approach 1
# *****************************************************

ll_Beta_1 <- function(parvector){
  b0 <- parvector[1]
  b1 <- parvector[2]
  b2 <- parvector[3]
  c0 <- parvector[4]
  c1 <- parvector[5]
  c2 <- parvector[6]
  xb <- b0 + b1 * LTV + b2 * purpose1
  wc <- c0 + c1 * LTV + c2 * purpose1
  mu = 1/(1+exp(-xb))
  delta <- exp(wc)
  alpha <- mu * delta
  beta <- (1-mu)*delta
  sum(log(
    gamma(alpha + beta) / ( gamma(alpha) * gamma(beta)) *
      (lgd_time ^ (alpha -  1)) * (1- lgd_time)^(beta -1)
    )
  )
}

beta_model_1 <-optim(par = c(0,0,0,0,0,0), 
                                 fn = ll_Beta_1, method = "BFGS", 
                                 control = list(fnscale = -1))

estimates_beta_model_1            <- as.data.frame(beta_model_1$par)
row.names(estimates_beta_model_1) <- c("b0", "b1", "b2", "c0", "c1", "c2")
colnames(estimates_beta_model_1)  <- c("Parameter Estimates")
print(estimates_beta_model_1)

b0 <- estimates_beta_model_1[1,1]
b1 <- estimates_beta_model_1[2,1]
b2 <- estimates_beta_model_1[3,1]
c0 <- estimates_beta_model_1[4,1]
c1 <- estimates_beta_model_1[5,1]
c2 <- estimates_beta_model_1[6,1]
xb <- b0 + b1 * LTV + b2 * purpose1
wc <- c0 + c1 * LTV + c2 * purpose1
mu <- 1/(1+exp(-xb))
delta        <- exp(wc)
pred_mu_1    <- mu
pred_delta_1 <- delta

plot(pred_mu_1, lgd_time, main = "Real-fit plot of beta regrassion",
     xlab = "Predicted values", ylab = "lgd_time", pch = 16, cex = 0.5, las =1)

summary(lm(lgd_time ~ pred_mu_1))

# *****************************************************
# Beta Regression - approach 2
# *****************************************************
Beta_model_2 <- betareg(lgd_time ~ LTV + purpose1 |  
                        LTV + purpose1,
                        link ="logit", link.phi = "log",
                        type = "ML")

summary(Beta_model_2)

estimates_beta_model_2 <- as.data.frame(Beta_model_2$coefficients)
print(estimates_beta_model_2)

pred_mu_2    <- predict(Beta_model_2, type = "response")
pred_delta_2 <- predict(Beta_model_2, type = "precision")

plot(pred_mu_2, lgd_time, main = "Real-fit plot of beta regression with betareg()",
     xlab = "Predicted values", ylab = "lgd_time", pch = 16, cex = 0.5, las =1)

summary(lm(lgd_time ~ pred_mu_2))

# **************************************************************************
# TOBIT RERESSION 
# **************************************************************************

# *****************************************************
# Approach 1
# *****************************************************
ll_Tobit_1<- function(parvector){
  b0 <- parvector[1]
  b1 <- parvector[2]
  b2 <- parvector[3]
  sigma <- parvector[4]
  xb <- b0 + b1 * LTV + b2 * purpose1
  mu <- 1/(1+exp(-xb))
  sum(
    ifelse(event == 1, log(dnorm(lgd_time, xb, sigma)),log(pnorm(0, xb, sigma)))
  )
}

Tobit_model_1 <-optim(par = c(0,0,0,1), 
                      fn = ll_Tobit_1, method = "BFGS", 
                      control = list(fnscale = -1))

estimates_Tobit_model_1            <- as.data.frame(Tobit_model_1$par)
row.names(estimates_Tobit_model_1) <- c("b0", "b1", "b2", "sigma")
colnames(estimates_Tobit_model_1)  <- c("Parameter Estimates")
print(estimates_Tobit_model_1)

b0 <- estimates_Tobit_model_1[[1]][1]
b1 <- estimates_Tobit_model_1[[1]][2]
b2 <- estimates_Tobit_model_1[[1]][3]
b3 <- estimates_Tobit_model_1[[1]][4]

xbeta  <- b0 + b1*LTV + b2* purpose1

lambda <- dnorm((0.00001 - xbeta) / b3) / (1-pnorm((0.00001 - xbeta) / b3))

cexpct_lgd_time <- xbeta + lambda*b3
plot(cexpct_lgd_time, lgd_time, main = "Real-fit plot of Tobit regression-Approach 1",
     xlab = "Predicted values", ylab = "lgd_time", pch = 16, cex = 0.5, las =1)

summary(lm(lgd_time ~ cexpct_lgd_time))

# *****************************************************
# Approach 2
# *****************************************************
Tobit_model_2 <- censReg(lgd_time ~ LTV + purpose1, 
          left = 0.00001,start = c(0,0,0,log(1)))

print(summary(Tobit_model_2), logSigma = F)

b0 <- Tobit_model_2$estimate[1]
b1 <- Tobit_model_2$estimate[2]
b2 <- Tobit_model_2$estimate[3]
b3 <- Tobit_model_2$estimate[4]

xbeta  <- b0 + b1*LTV + b2 * purpose1
lambda <- dnorm((0.00001 - xbeta) / exp(b3)) / (1-pnorm((0.00001 - xbeta) / exp(b3)))

cexpct_lgd_time <- xbeta + lambda*exp(b3)
plot(cexpct_lgd_time, lgd_time, main = "Real-fit plot of Tobit regression-Approach 2",
     xlab = "Predicted values", ylab = "lgd_time", pch = 16, cex = 0.5, las =1)

summary(lm(lgd_time ~ cexpct_lgd_time))

# **************************************************************************
# Heckman Sample Selection Model
# **************************************************************************

ll_Heckman <- function(parvector){
  b0 <- parvector[1]
  b1 <- parvector[2]
  b2 <- parvector[3]
  sigma <- parvector[4]
  a0 <- parvector[5]
  rho <- parvector[6]
  y_o <- lgd_time
  xb_s <- a0
  xb_o <- b0 + b1 * LTV + b2 * purpose1
  
  sum(
    ifelse(event == 0, log(pnorm(-xb_s)),
        (log(pnorm((xb_s + rho / sigma *(y_o - xb_o)) /
        sqrt(1-rho^2))) - 0.5 * log(2*pi)-log(sigma) - 0.5 *(y_o -xb_o)^2/sigma ^2)
    )
  )
}

Heckman_model <-optim(par = c(0,0,0,1,0,0), 
                     fn = ll_Heckman, method = "Nelder-Mead", 
                     control = list(fnscale = -1))

estimates_Heckman_model <- as.data.frame(Heckman_model$par)
row.names(estimates_Heckman_model) <- 
  c("lgdtime.Interc", "lgdtime.LTV", "lgdtime.purpose",
  "Sigmalgdtime", "event.Intercept", "_Rho")
colnames(estimates_Heckman_model) <- c("Parameter Estimates")
print(estimates_Heckman_model)

# **************************************************************************
# Censored Beta Regression
# **************************************************************************

ll_Beta_2 <- function(parvector){
  a0 <- parvector[1]
  b0 <- parvector[2]
  b1 <- parvector[3]
  b2 <- parvector[4]
  c0 <- parvector[5]
  c1 <- parvector[6]
  c2 <- parvector[7]
  
  y <- lgd_time
  xa <- a0
  xb <- b0 + b1 * LTV + b2 * purpose1
  xc <- c0 + c1 * LTV + c2 * purpose1
  
  mu_xa = 1/(1+exp(-xa))
  mu_xb = 1/(1+exp(-xb))
  
  delta <- exp(xc)
  
  alpha <- mu_xb * delta
  beta <- (1-mu_xb)*delta
  
  sum(
    ifelse(event ==0, log(1- mu_xa), 
           log(mu_xa * gamma(alpha + beta)/
           (gamma(alpha)*gamma(beta)) * (y ^ (alpha -1)) * (1-y)^(beta -1))
    )
  )
}

Beta_censored_model <-optim(par = c(0,0,0,0,0,0,0), 
                      fn = ll_Beta_2, method = "BFGS", 
                      control = list(fnscale = -1))

estimates_Beta_censored_model            <- as.data.frame(Beta_censored_model$par)
row.names(estimates_Beta_censored_model) <- c("a0","b0", "b1", "b2", "c0", "c1", "c2")
colnames(estimates_Beta_censored_model)  <- c("Parameter Estimates")
print(estimates_Beta_censored_model)