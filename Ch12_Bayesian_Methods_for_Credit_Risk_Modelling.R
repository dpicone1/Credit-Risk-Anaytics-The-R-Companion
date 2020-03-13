remove(list = ls())
mortgage <- read.csv("mortgage.csv", header = T)

install.packages("MCMCpack")
install.packages("coda")
install.packages("Bolstad2")
install.packages("rjags")

library("MCMCpack")
library("coda")
library("Bolstad2")
library("rjags")

set.seed(12345)
nrows <- nrow(mortgage)
ref <- sample(nrows,nrows*0.01,replace = F )
mortgage_sample <- mortgage[ref, ]
model_probit <- glm(default_time ~ FICO_orig_time + LTV_orig_time + gdp_time,
                    family = binomial(link = "probit"),data = mortgage_sample)
summary(model_probit)

# MCMC method to estimate Bayesian Probit Model
# ************************************************************************
# Model 1
# ************************************************************************
model_probit_MCMC <- MCMCprobit(
  default_time ~ FICO_orig_time + LTV_orig_time + gdp_time,
  data = mortgage_sample,
  seed = 12345,
  thin = 2,
  b0 = 0,
  B0 = 0.001,
  burnin = 1000,
  mcmc = 10000)


summary(model_probit_MCMC)

crosscorr(model_probit_MCMC)
autocorr.diag(model_probit_MCMC)
geweke.diag(model_probit_MCMC)
raftery.diag(model_probit_MCMC)
heidel.diag(model_probit_MCMC)
effectiveSize(model_probit_MCMC)


par(mfrow = c(2,4))
plot(model_probit_MCMC,auto.layout =  F, ask = F)
par(mfrow = c(1,1))

par(mfrow = c(2,2))
autocorr.plot(model_probit_MCMC,auto.layout =  F, ask = F)
par(mfrow = c(1,1))

# ************************************************************************
# Model 2
# ************************************************************************
P <- diag(x = rep(1/1000, 4))
P[3,3] <- 1/.5
model_probit_MCMC2 <- MCMCprobit(
  default_time ~ FICO_orig_time + LTV_orig_time + gdp_time,
  data = mortgage_sample,
  seed = 12345,
  thin = 2,
  b0 = c(0,0,3,0),
  B0 = P,
  burnin = 1000,
  mcmc = 10000)


summary(model_probit_MCMC2)

crosscorr(model_probit_MCMC2)
autocorr.diag(model_probit_MCMC2)
geweke.diag(model_probit_MCMC2)
raftery.diag(model_probit_MCMC2)
heidel.diag(model_probit_MCMC2)
effectiveSize(model_probit_MCMC2)


par(mfrow = c(2,4))
plot(model_probit_MCMC2,auto.layout =  F, ask = F)
par(mfrow = c(1,1))

par(mfrow = c(2,2))
autocorr.plot(model_probit_MCMC2,auto.layout =  F, ask = F)
par(mfrow = c(1,1))


# MCMC method to estimate Bayesian Survival Model
# ************************************************************************
est <- BayesCPH(
  y = mortgage_sample$default_time,
  t = mortgage_sample$time,
  x = cbind(
    mortgage_sample$FICO_orig_time,
    mortgage_sample$LTV_orig_time,
    mortgage_sample$gdp_time),
  plot = T)

# **********************************************************************
# Introduction to JAGS
# **********************************************************************
default_mean <- tapply(mortgage$default_time, mortgage$time, FUN = mean)
probit_dr    <- qnorm(default_mean)

# correlation Estimation with Bayesian Statistics
# Model 1
x = as.double(probit_dr)

model1.string <- "
model {
for (i in 1:60) {
x[i] ~ dnorm(beta0, 1/sigma2)
}
beta0 ~ dnorm(0,1/10000)
sigma2 ~ dunif(0,1)
}
"
model1.spec <- textConnection(model1.string)
jags <- jags.model(
  model1.spec,
  data = list('x'=x),
  n.chains = 1,
  n.adapt = 1000
)
correlation1 <- coda.samples(jags, c('beta0', 'sigma2'), 10000)
summary(correlation1)
crosscorr(correlation1)
autocorr.diag(correlation1)
geweke.diag(correlation1)
raftery.diag(correlation1)
heidel.diag(correlation1)
effectiveSize(correlation1)


par(mfrow = c(2,2))
plot(correlation1,auto.layout =  F, ask = F)

par(mfrow = c(1,2))
autocorr.plot(correlation1,auto.layout =  F, ask = F)
par(mfrow = c(1,1))

# Model 2
x = as.double(probit_dr)

model2.string <- "
model {
for (i in 1:60) {
x[i] ~ dnorm(beta0, 1/sigma2)
}
beta0 ~ dnorm(0,1/10000)
sigma2 ~ dunif(0,10)
}
"
model2.spec <- textConnection(model2.string)
jags2 <- jags.model(
  model2.spec,
  data = list('x'=x),
  n.chains = 1,
  n.adapt = 1000
)
correlation2 <- coda.samples(jags2, c('beta0', 'sigma2'), 10000)
summary(correlation2)
crosscorr(correlation2)
autocorr.diag(correlation2)
geweke.diag(correlation2)
raftery.diag(correlation2)
heidel.diag(correlation2)
effectiveSize(correlation2)


par(mfrow = c(2,2))
plot(correlation2,auto.layout =  F, ask = F)

par(mfrow = c(1,2))
autocorr.plot(correlation2,auto.layout =  F, ask = F)
par(mfrow = c(1,1))

# ***********************************************************
# PD estimation of Low Default Portfolios
# ***********************************************************
# Model 1
tmp <- subset(mortgage, FICO_orig_time > 810)
default_freq <- tapply(tmp$default_time, tmp$time, sum) 
Total_freq   <- tapply(tmp$default_time, tmp$time, length) 

x <- as.double(default_freq)
N <- as.double(Total_freq)

modelLowPD.string <- "
model {
for (i in 1:47) {
x[i] ~ dbin(pd, N[i])
}
pd ~ dbeta(1,1)
}
"
modelLowPD.spec <- textConnection(modelLowPD.string)
jagsPD <- jags.model(
  modelLowPD.spec,
  data = list('x' = x, 'N' = N),
  n.chains = 1,
  n.adapt = 1000
)
PD <- coda.samples(jagsPD, c('pd'), 10000)
summary(PD)
crosscorr(PD)
autocorr.diag(PD)
geweke.diag(PD)
raftery.diag(PD)
heidel.diag(PD)
effectiveSize(PD)


par(mfrow = c(1,2))
plot(PD,auto.layout =  F, ask = F)

par(mfrow = c(1,1))
autocorr.plot(PD,auto.layout =  F, ask = F)

 
# Model 2 - run a probit model and then re fit it with the Bayesian approach
PD_model_probit <- glm(default_time ~ LTV_orig_time + gdp_time,
                       data = tmp, family = binomial(link = 'probit'))
summary(PD_model_probit)

x  <- as.double(tmp$default_time)
y1 <- as.double(tmp$LTV_orig_time)
y2 <- as.double(tmp$gdp_time)

modelLowPD_probit.string <- "
model {
for (i in 1:1780) {
pd[i] = pnorm(beta0 + beta1 * y1[i] + beta2*y2[i], 0 , 1)
x[i] ~ dbern(pd[i])
}
beta0 ~ dnorm(0,1/10000)
beta1 ~ dnorm(0,1/10000)
beta2 ~ dnorm(0,1/10000)
}
"
modelLowPD_probit.spec <- textConnection(modelLowPD_probit.string)
jagsPD_probit <- jags.model(
  modelLowPD_probit.spec,
  data = list('x' = x, 'y1' = y1, 'y2' = y2),
  n.chains = 1,
  n.adapt = 1000
)
PD_probit <- coda.samples(jagsPD_probit, c('beta0', 'beta1', 'beta2'), 10000)
summary(PD_probit)
crosscorr(PD_probit)
autocorr.diag(PD_probit)
geweke.diag(PD_probit)
raftery.diag(PD_probit)
heidel.diag(PD_probit)
effectiveSize(PD_probit)


par(mfrow = c(3,2))
plot(PD_probit,auto.layout =  F, ask = F)

par(mfrow = c(1,3))
autocorr.plot(PD_probit,auto.layout =  F, ask = F)
par(mfrow = c(1,1))
