remove(list = ls())

install.packages("mnormt")
install.packages("TSA")
install.packages("moments")
install.packages("lmtest")

library("mnormt")
library("TSA")
library("moments")
library("lmtest")

# **************************************************************************
# MODELLING LOSS DISTRIBUTION WITH CORRELATED DEFAULTS
# **************************************************************************

# *****************************************************
# Analytic Solution
# *****************************************************
PD1 <- 0.01
rho1<- 0.1
PD2 <- 0.05
rho2<- 0.2

p_i <- seq(0.0001, 0.2, by = 0.0001)
length(p_i)
g_1 <- numeric(length(p_i))
g_2 <- numeric(length(p_i))

g_1 <- sqrt(1-rho1)/ sqrt(rho1)*exp(0.5*(qnorm(p_i)**2)-0.5/
                rho1 *(qnorm(PD1)-sqrt(1-rho1)*qnorm(p_i))**2)
g_2 <- sqrt(1-rho2)/ sqrt(rho2)*exp(0.5*(qnorm(p_i)**2)-0.5/
                rho2 *(qnorm(PD2)-sqrt(1-rho2)*qnorm(p_i))**2)

plot(p_i, g_1 , type = "l")
lines(p_i, g_2, lty  = 2)
legend(x = "topright", legend = c("PD1", "PD2"), horiz = F,lty = c(1,2), bty = "n")

# *****************************************************
# Numerical Solution
# *****************************************************
rho <- 0.2
n   <- 100
p   <- 0.05

func_v <- function(x){
  cpd <- pnorm((1 / sqrt(1-rho))*(qnorm(p) - sqrt(rho)*x))
  v   <- pbinom(k, n, cpd) * (exp(-(x*x)/2)) / (sqrt(2*pi))
  return (v)
}

k <- 0
z <- 0
out <- data.frame(quantile = numeric(), default_rate = numeric(), z = numeric())

while(z< 0.9999){
  z <- integrate(f = func_v, lower = -Inf, upper=Inf, abs.tol = 1.34E-15)$value
  quantile <- k
  default_rate <- quantile / n
  out[k+1, ] <- c(quantile, default_rate, z)
  k <- k+1
}

prob    <- vector(mode = "numeric", length = nrow(out))
prob[1] <- out$z[1]
for (i in 2:nrow(out)){
  prob[i] <- out$z[i] - out$z[i-1]
}

numeric <- data.frame(out, prob)
plot(numeric$default_rate, numeric$prob, type = "h", 
     main = "Numerically Computed Loss Distribution", xlab = "Loss", ylab = "Probability")

# Expected Default Rate
(mean.numeric     <- round(sum(numeric$prob *numeric$default_rate),4))
(variance.numeric <- round(sum((numeric$default_rate - mean)^2*numeric$prob),4))
numeric$prob%*%numeric$default_rate

# *****************************************************
# Monte Carlo Simulation
# *****************************************************
n_sim <- 100000
n     <- 10000
p     <- 0.05
w_2   <- 0.2
w     <- rep(sqrt(w_2), times = n)
ksi   <- rep(sqrt(1-w_2), times = n)

pd <- rep(p, n)
sum_defaults <- numeric(n_sim)
default_rate <- numeric(n_sim)

for (i in 1:n_sim){
  z   <- rnorm(1)
  eps <- rnorm(n)
  lnv <- w*z + ksi*eps
  thresold <- qnorm(pd)
  defaults <- lnv < thresold
  sum_defaults[i] <- sum(defaults)
  default_rate[i] <- sum_defaults[i]/n
}

sim <- data.frame(n, sum_defaults, default_rate)
colnames(sim) <- c("N", "Defaults", "Default Rate")

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
    Kurtosis = kurtosis(x),
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
fun_statistics.mean.variance <- function(x){
  rbind(
    N       = length(x),
    Sum_Obs = sum(x),
    Mean    = mean(x),
    Median  = median(x),
    Variance= round(var(x),3),
    StandDev= sd(x)
  )
}

options(scipen = 999)
statistics_default_rate_MC           <- fun_statistics(default_rate)
colnames(statistics_default_rate_MC) <- "Loss"
print(statistics_default_rate_MC)
hist(default_rate, breaks = 50, freq = F, xlab = "Loss", ylab = "Percent",
     main = "Distribution of Loss")

round(fun_statistics.mean.variance(default_rate),4)

# ************************************************************************
# ESTIMATING CORRELATIONS

# with three methods:
# 1 - Methods of Moment
# 2 - Maximum Likelihood
# 3 - Probit-Logit Regression
# ************************************************************************

# *****************************************************
# 1 - Methods of Moment
# *****************************************************
dat   <- read.csv("mortgage.csv", header = T)
tmp   <- dat[order(dat$time),]
print(unique(dat$time))
print(max(unique(dat$time)))
print(min(unique(dat$time)))

means <- merge(aggregate(tmp[, c("default_time", "gdp_time")],
                         by = list(time = tmp[,"time"]), FUN = length),
               aggregate(tmp[, c("default_time", "gdp_time")],
                         by = list(time = tmp[,"time"]), FUN = mean),
               by = "time", suffixes = c(".Freq",".Mean"))

fun_lag <- function(x, k){
  c(rep(NA, k), x)[1:length(x)]
}

time          <- means$time
freq          <- means$default_time.Freq

# This is the expected default time given time = i, with i = {1,2,..,60}
default_time  <- means$default_time.Mean 
# This is the expected gdp time given time = i, with i = {1,2,..,60}
gdp_time      <- means$gdp_time.Mean
n_default     <- default_time * freq
sum(n_default)/sum(freq)

default_time_1<- fun_lag(default_time,1) # default_time(t - 1)
probit_dr     <- qnorm(default_time) 
tmp2          <- data.frame(time, freq, default_time, gdp_time, 
                            n_default, default_time_1, probit_dr)
default_rate  <- default_time
lambda_hat    <- mean(default_rate)
lambda_hat2   <- lambda_hat**2
n_inverse     <- 1/freq
threshold     <- qnorm(lambda_hat)
mse           <- as.numeric(t(default_rate - lambda_hat) %*% 
                    (default_rate - lambda_hat)/nrow(tmp2))
var_lambda_hat<- as.numeric(
      (mse - t(n_inverse) %*% rep(1,60) * lambda_hat * (1-lambda_hat) / nrow(tmp2))
      / (1-t(n_inverse) %*% rep(1,60) / nrow(tmp2))
)
f <- function(x){
    abs(
      pmnorm(
        c(threshold, threshold), 
        mean = c(0,0),
        varcov = matrix(c(1,x,x,1), nrow = 2, ncol = 2)
        ) -
    lambda_hat **2 - var_lambda_hat
    )
}
opt                <- optimize(f = f, interval = c(0,1))
optmization_result <- data.frame(opt[1], opt[2])
colnames(optmization_result) <- c("Estimate", "Objective Function")
print(optmization_result)

# Variance - Covariance Matrix
matrix(c(1,optmization_result$Estimate,
         optmization_result$Estimate,1), nrow = 2)
cat("the Asset correlation is", optmization_result$Estimate)

# *****************************************************
# 2 - Maximum Likelihood
# *****************************************************
c_t    <- qnorm(tmp2$default_time)
c_bar  <- mean(c_t)
var_c  <- as.numeric(t(c_t) %*% c_t  / nrow(tmp2) - c_bar**2)
pd_ML  <- pnorm(c_bar / sqrt(1+var_c))
rho_ML <- var_c /sqrt(1+var_c)
print(data.frame(c_bar,var_c,pd_ML,rho_ML))      
cat("the Estimated Asset Correlation is", rho_ML)
cat("the Estimated Default Probability is", pd_ML)
# From the data:
mean(default_rate)
# Default rate weighted by the frequency
sum(default_time * freq)/sum(freq)

pnorm(0,0,1)
qnorm(pnorm(0,0,1),0,1)

# *****************************************************
# Probit-Logit Regression without Covariates
# *****************************************************
Probit_dr <- tmp2$probit_dr
east1     <- lm(Probit_dr ~ 1)
summary(east1)
plot(rstudent(east1), type =  "h", xlab = "Observation", 
     ylab = "t-Student Residuals",las = 1)
abline(h = mean(rstudent(east1)))
plot(Probit_dr, type = "p", xlab = "Observation",
     ylab = "Probit dr", pch = 16)
plot( cooks.distance(east1), type = "h", xlab = "Observation",
     ylab = "Cook's D", pch = 16)
abline(h = 0)
qqnorm(residuals(east1),  xlab = "Quantile",ylab = "Residuals", pch = 16)
qqline(residuals(east1))

# *****************************************************
# Probit-Logit Regression with Covariates
# *****************************************************
default_time_1 <- tmp2$default_time_1
east2          <- lm(Probit_dr ~ qnorm(default_time_1))
summary(east2)

plot(rstudent(east2), type =  "h", xlab = "Observation", 
     ylab = "t-Student Residuals",las = 1)
abline(h = mean(rstudent(east2)))
plot(Probit_dr, type = "p", xlab = "Observation",
     ylab = "Probit dr", pch = 16)
plot( cooks.distance(east2), type = "h", xlab = "Observation",
      ylab = "Cook's D", pch = 16)
abline(h = 0)
qqnorm(residuals(east2),  xlab = "Quantile",ylab = "Residuals", pch = 16)
qqline(residuals(east2))

# additional the gdp covariates
gdp_time <- tmp2$gdp_time
east3    <- lm(Probit_dr ~ qnorm(default_time_1) + gdp_time)
summary(east3)

plot(rstudent(east3), type =  "h", xlab = "Observation", 
     ylab = "t-Student Residuals",las = 1)
abline(h = mean(rstudent(east3)))
plot(Probit_dr, type = "p", xlab = "Observation",
     ylab = "Probit dr", pch = 16)
plot( cooks.distance(east3), type = "h", xlab = "Observation",
      ylab = "Cook's D", pch = 16)
abline(h = 0)
qqnorm(residuals(east3),  xlab = "Quantile",ylab = "Residuals", pch = 16)
qqline(residuals(east3))

# *****************************************************
# AR model
# *****************************************************
east4 <- arima(Probit_dr , order = c(3,0,0), xreg  = gdp_time, 
               include.mean = TRUE, method = "ML")
summary4 <- data.frame(
  SSE          = sum(east4$residuals^2),
  MSE          = east4$sigma2,
  Root_MSE     = sqrt(east4$sigma2),
  Log_Likehood = east4$loglik,
  AIC          = east4$aic,
  DFE          = east4$nobs - length(east4$coef),
  Observations = east4$nobs)
print (summary4)

coeftest(east4, vcov. = east4$var.coef, df = east4$nobs - length(east4$coef))

plot(rstandard(east4), type = "h", xlab = "Observations", ylab="Standard Res")
abline(h = mean(rstandard(east4)))
plot(Probit_dr, type = "o", xlab = "Observation", ylab = "Probit dr", las = 1, pch = 16)
hist(east4$residuals, breaks = 10, freq = FALSE, xaxs = "i", yaxs = "i",
     xlab = "Residuals", ylab = "Density", main = "", las = 1)
box()
curve(dnorm(x, mean = 0, sd = sd(east4$residuals)), add = T)
qqnorm(residuals(east4),  xlab = "Quantile",ylab = "Residuals", pch = 16)
qqline(residuals(east4))

acf_east4  <- acf(east4$residuals, lag.max = 25, 
                 type = "correlation", plot = F, drop.lag.0 = F)
pacf_east4 <- pacf(east4$residuals, lag.max = 25, 
                   type = "partial", plot = F, drop.lag.0 = F)

plot(acf_east4, type = "h", main = "", ylim = c(-1,1), 
     las = 1, ci.col = "black")
plot(pacf_east4, type = "h", main = "", ylim = c(-1,1), 
     las = 1, ci.col = "black")

# *****************************************************
# Extension
# *****************************************************
# Model Specifications other than Gaussian

# Set Parameters
nsim <- 100000
n    <- 1000

# Parameters for the Gaussian
pd_N    <- 0.0212
c       <- qnorm(pd_N)
rho_N   <- 0.055
ksi_N   <- sqrt(1-rho_N)
var_cpd <- pmnorm(c(c,c), mean = c(0,0), 
            varcov = matrix(c(1, rho_N, rho_N, 1), nrow = 2, ncol = 2)) - 
            pd_N**2
# Parameters for the CR+ model
pd_C    <- pd_N
alpha_C <- pd_C**2/var_cpd
beta_C  <- 1/alpha_C

# Parameters for the Clayton Copula
pd_T    <- pd_N
theta   <- 0.020401

# Genarating Default Distributin for Gaussian Model
z_N   <- rnorm(nsim)
cpd_N <- pnorm((qnorm(pd_N) - sqrt(rho_N) * z_N) / ksi_N)
out_N <- rbinom(n = nsim, size = n, prob = cpd_N)/n

# Genarating Default Distributin for CR+ Model
z_C   <- rgamma(n = nsim, shape = alpha_C, scale = beta_C)
cpd_C <- n * pd_C * z_C
out_C <- rpois(n = nsim, lambda = cpd_C)/n

# Genarating Default Distributin for Clayton Model
z_T   <- rgamma(n = nsim, shape = 1/theta, scale = 1)
cpd_T <- exp(-(pd_T**(-theta) - 1) * z_T)
out_T <- rbinom(n = nsim, size =n, prob = cpd_T)/n

# Genarating Output data
svsim_N <- data.frame(out_N)
svsim_C <- data.frame(out_C)
svsim_T <- data.frame(out_T)
dim(svsim_N)
dim(svsim_C)
dim(svsim_T)

statistcs_svsim_T <- fun_statistics(svsim_T$out_T)
statistcs_svsim_N <- fun_statistics(svsim_N$out_N)
statistcs_svsim_C <- fun_statistics(svsim_C$out_C)
colnames(statistcs_svsim_T) <- "svsim_T"
colnames(statistcs_svsim_N) <- "svsim_N"
colnames(statistcs_svsim_C) <- "svsim_C"

print(statistcs_svsim_T)
print(statistcs_svsim_C)
print(statistcs_svsim_N)

cbind(statistcs_svsim_T, statistcs_svsim_C, statistcs_svsim_N)

temp_t        <- as.data.frame(table(svsim_T)/nrow(svsim_T))
names(temp_t) <- c("Default_Rate", "Clayton_Copula")
dim(temp_t)
temp_n        <- as.data.frame(table(svsim_N)/nrow(svsim_N))
names(temp_n) <- c("Default_Rate", "Gaussian_Copula")
dim(temp_n)
temp_c        <- as.data.frame(table(svsim_C)/nrow(svsim_C))
names(temp_c) <- c("Default_Rate", "CR+")
dim(temp_c)

temp_t1       <- data.frame(cbind(as.numeric(as.character(temp_t$Default_Rate)),
                            temp_t$Clayton_Copula))
names(temp_t1)<- c("Default_Rate", "Clayton_Copula")
dim(temp_t1)

temp_g1       <- data.frame(cbind(as.numeric(as.character(temp_n$Default_Rate)),
                            temp_n$Gaussian_Copula))
names(temp_g1)<- c("Default_Rate", "Gaussian_Copula")
dim(temp_g1)

temp_c1       <- data.frame(cbind(as.numeric(as.character(temp_c$Default_Rate)),
                            temp_c$CR))
names(temp_c1)<- c("Default_Rate", "CR")
dim(temp_c1)

temp <- merge(merge(temp_t1, temp_g1, by = "Default_Rate", all = T),
              temp_c1, by = "Default_Rate", all = T)
dim(temp)

plot(temp$Default_Rate,  temp$Gaussian_Copula, 
     xlab = "Default Rates", ylab = "Probabilities",xlim = c(0,0.08),type = "l", 
                                              lty=1, col = "Red")
lines(temp$Default_Rate, temp$Clayton_Copula, lty = 2, col = "Blue")
lines(temp$Default_Rate, temp$CR,             lty = 3, col = "magenta")
legend(x = "topright", legend = c("Gaussian","Clayton","CR+"),
       lty = c(1,2,3), col = c("Red","Blue", "magenta"), bty = "n")