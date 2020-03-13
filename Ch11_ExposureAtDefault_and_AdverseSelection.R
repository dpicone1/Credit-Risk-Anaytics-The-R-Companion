remove(list = ls())

#install.packages("data.table")
#install.packages("pastecs")
#install.packages("betareg")
#install.packages("nnet")

library("data.table")
library("pastecs")
library("betareg")
library("nnet")
library("survival")

# **************************************************************************
# DATA PREPROCESSING - VARIABLE TIME HORIZON METHOD
# **************************************************************************
mortgage <- read.csv("mortgage.csv", header = T)      

ncol(mortgage)
nrow(mortgage)
ncol(mortgage)*nrow(mortgage)
dim(mortgage)

mortgage <- mortgage[complete.cases(mortgage),]
mortgage <- mortgage[order(mortgage$id),]
nrow(mortgage)

temp_mortgage <- mortgage[,c("id", "balance_time")]
# For the variable Balance creating 4 lagging variables [t - 1, t - 2, t - 3, t - 4]
lag_mortgage  <- setDT(temp_mortgage)[order(id), 
    c("lag1", "lag2", "lag3", "lag4") := shift(balance_time, 1:4), id][order(id)]
lag_mortgage <- as.data.frame(lag_mortgage)
mortgage2    <- cbind(mortgage, lag_mortgage[,3:6])
remove(lag_mortgage, temp_mortgage)

# For id 46 and 47, we can read how the balance changes during a period of 4 months
subset(mortgage2, id == 46|id == 47, select = 
         c("id", "time", "default_time","balance_orig_time",
           "balance_time", "lag1", "lag2", "lag3", "lag4"))

# Cearly the Exposure at default changes!!!!!!!

# **************************************************************************
# DEFINITON
# **************************************************************************

mortgage2$drawn   <- mortgage2$lag4
mortgage2$limit   <- mortgage2$balance_orig_time
mortgage2$exposure<- mortgage2$balance_time

# We remove loans where drawn is NA and limit is zero
nrow(mortgage2)
mortgage2<- subset(mortgage2, limit != 0 & drawn != "NA")
nrow(mortgage2)
subset(mortgage2, id == 46|id == 47, select = 
         c("id", "time", "default_time","balance_orig_time",
           "balance_time", "lag1", "lag2", "lag3", "lag4"))
# Look atloan id 47, Now it has gone as it does not have 
# a balance t - 4

# Create caps for exposure and draw.
# So Further Advances are NOT modelled
mortgage2$exposure <- ifelse(mortgage2$exposure> mortgage2$limit,
                             mortgage2$limit,
                             mortgage2$exposure)

mortgage2$drawn    <- ifelse(mortgage2$drawn> mortgage2$limit,
                             mortgage2$limit,
                             mortgage2$drawn)

# Original Balance = 0 are removed 
mortgage2          <- subset(mortgage2, exposure != 0)
nrow(mortgage2)

# Compute Conversion Measures
mortgage2$CCF <- ifelse(mortgage2$drawn == mortgage2$limit,0,
          (mortgage2$exposure - mortgage2$drawn)/
          (mortgage2$limit - mortgage2$drawn))
mortgage2$CEQ <- (mortgage2$exposure - mortgage2$drawn)/ (mortgage2$limit)
mortgage2$LCF <- (mortgage2$exposure / mortgage2$limit)
mortgage2$UACF<- (mortgage2$exposure / mortgage2$drawn)

conversionMeasures <- cbind(mortgage2$CCF, 
                            mortgage2$CEQ,
                            mortgage2$LCF,
                            mortgage2$UACF)

conversionMeasuresPercentiles <- apply(
  conversionMeasures, 2, quantile, probs = c(0.01, 0.99), na.rm = T)
colnames(conversionMeasuresPercentiles) <- c("CCF", 
                                             "CEQ",
                                             "LCF",
                                             "UACF")
print(conversionMeasuresPercentiles)

# Apply Floors
mortgage2$CCF <- ifelse(mortgage2$CCF <= -18.0502849,-18.0502849, mortgage2$CCF)
mortgage2$CEQ <- ifelse(mortgage2$CEQ <= -0.1297378,-0.1297378,   mortgage2$CEQ)
mortgage2$LCF <- ifelse(mortgage2$LCF <=  0.3724237,0.3724237,    mortgage2$LCF)
mortgage2$UACF<- ifelse(mortgage2$UACF<=  0.7492269,0.7492269,    mortgage2$UACF)


# Apply Caps
mortgage2$CCF  <- ifelse(mortgage2$CCF > 0.9999999,0.9999999, mortgage2$CCF)
mortgage2$CEQ  <- ifelse(mortgage2$CEQ > 0.0102912,0.0102912, mortgage2$CEQ)
mortgage2$LCF  <- ifelse(mortgage2$LCF >  .9999999,.9999999,  mortgage2$LCF)
mortgage2$UACF <- ifelse(mortgage2$UACF > 1.0105358,1.0105358,mortgage2$UACF)

# Compute log tranformation
mortgage2$CCF_t <- -log(1 - mortgage2$CCF)
mortgage2$CEQ_t <- log((1 + mortgage2$CEQ)/(1 - mortgage2$CEQ))
mortgage2$LCF_t <- log(mortgage2$LCF / (1-mortgage2$LCF))
mortgage2$UACF_t<- log(mortgage2$UACF)

# Keep the data where deafault occured , i.e. default_time == 1
mortgage2 <- subset(mortgage2, default_time == 1)
descriptive_Stats <- stat.desc(mortgage2[ , 
  c("CCF", "CEQ","LCF", "UACF", "CCF_t","CEQ_t", "LCF_t","UACF_t")])                       
print(descriptive_Stats, digit = 5)


hist(mortgage2$CCF, freq = F, main = "Distr of CCF", xlab = "CCF", ylab="Perc",
     breaks = 40); box()

hist(mortgage2$CEQ, freq = F, main = "Distr of CEQ", xlab = "CEQ", ylab="Perc",
     breaks = 40); box()

hist(mortgage2$LCF, freq = F, main = "Distr of LCF", xlab = "LCF", ylab="Perc",
     breaks = 40); box()

hist(mortgage2$UACF, freq = F, main = "Distr of UACF", xlab = "UACF", ylab="Perc",
     breaks = 40); box()

hist(mortgage2$CCF_t, freq = F, main = "Distr of CCF_t", xlab = "CCF_t", 
     ylab="Perc",
     breaks = 40); box()

hist(mortgage2$CEQ_t, freq = F, main = "Distr of CEQ_t", xlab = "CEQ_t", 
     ylab="Perc",
     breaks = 40); box()

hist(mortgage2$LCF_t, freq = F, main = "Distr of LCF_t", xlab = "LCF_t", 
     ylab="Perc",
     breaks = 40); box()

hist(mortgage2$UACF_t, freq = F, main = "Distr of UACF_t", xlab = "UACF_t", 
     ylab="Perc",
     breaks = 40); box()

# **************************************************************************
# LOANS WITH FLEXIBLE PAYMENTS
# **************************************************************************

# *****************************************************
# Linear regression
# *****************************************************

lm_CCF <- lm(CCF ~  LTV_time, data = mortgage2)
lm_CEQ <- lm(CEQ ~  LTV_time, data = mortgage2)
lm_LCF <- lm(LCF ~  LTV_time, data = mortgage2)
lm_UACF<- lm(UACF ~ LTV_time, data = mortgage2)

summary(lm_CCF)
summary(lm_CEQ)
summary(lm_LCF)
summary(lm_UACF)


par(mfrow = c(2,3))
plot(lm_LCF, which = 1:6)
par(mfrow = c(1,1))


lm_CCF_t <- lm(CCF_t  ~ LTV_time, data = mortgage2)
lm_CEQ_t <- lm(CEQ_t  ~ LTV_time, data = mortgage2)
lm_LCF_t <- lm(LCF_t  ~ LTV_time, data = mortgage2)
lm_UACF_t<- lm(UACF_t ~ LTV_time, data = mortgage2)

summary(lm_CCF_t)
summary(lm_CEQ_t)
summary(lm_LCF_t)
summary(lm_UACF_t)


par(mfrow = c(2,3))
plot(lm_LCF_t, which = 1:6)
par(mfrow = c(1,1))

# *****************************************************
# Beta Regression
# *****************************************************
# just for the variable LCF
betareg_LCF <- betareg(data = mortgage2, LCF ~ LTV_time | LTV_time)
summary(betareg_LCF)
plot(betareg_LCF$fitted.values, mortgage2$LCF)

# Estimate a linear regression model linking observed and
# fitted values
summary(lm(mortgage2$LCF ~ betareg_LCF$fitted.values))

# **************************************************************************
# CONTROLLING ADVERSE SELECTION IN PD MODELS
# **************************************************************************

# *****************************************************
# Discrete time Hazard Models
# *****************************************************
mortgage <- read.csv("mortgage.csv", header = T)      

dim(mortgage)


mortgage3             <- mortgage
mortgage3$status_time <- as.factor(mortgage3$status_time)
mortgage3$status_time <- relevel(mortgage3$status_time, "0")

# multinomial model from library nnet()
multi <- multinom(status_time ~ FICO_orig_time + LTV_time + gdp_time , 
                  data = mortgage3)
summary(multi)

# *****************************************************
# Estimation of default probabilities
# *****************************************************
# compute the proportion of observations where status_time = 1, which implies default 
default_mean_obs <- length(mortgage$status_time[mortgage$status_time == 1]) /
  length(mortgage$status_time)
# and compare it with teh average estimated default probabilities
# from the multinomial model
default_mean_pre <- mean(fitted(multi)[,2])
cbind(default_mean_obs, default_mean_pre)

# compute the proportion of observations where status_time = 1, which implies default
# but now with a function
default_obs_bytime <- by(
  mortgage,
  list(mortgage$time),
  FUN = function(x) length(x$status_time[x$status_time == 1]) / length(x$status_time)
)

# compute the average estimated probability of default from the multinomial model
mortgage$default_fit <- fitted(multi)[,2]

default_fit_bytime  <- by(
  mortgage,
  list(mortgage$time),
  FUN = function(x) mean(x$default_fit)
)

yaxes <- c(
  min(default_obs_bytime, default_fit_bytime), 
  max(default_obs_bytime, default_fit_bytime)
  )
plot(default_obs_bytime, type = "l", lty = 1, ylim = yaxes,
     xlab = "time", ylab = "DR and PD")
lines(default_fit_bytime, lty = 2)
legend(x = "topleft", cex = 0.8, legend = c("default time", "prob of default"),
       lty = c(1,2), bty = "n", horiz = F)


# *****************************************************
# Estimation of payoff probabilities
# *****************************************************
# Similar to teh default probabilities we do the same with
# the estimated payoff probabilities from the multinomial model
payoff_mean_obs <- length(mortgage$status_time[mortgage$status_time == 2]) /
  length(mortgage$status_time)
payoff_mean_fit <- mean(fitted(multi)[,3])

cbind(payoff_mean_obs,payoff_mean_fit)

payoff_obs_bytime <- by(
  mortgage,
  list(mortgage$time),
  FUN = function(x) length(x$status_time[x$status_time == 2]) / 
    length(x$status_time)
)

mortgage$payoff_fit <- fitted(multi)[,3]

payoff_fit_bytime  <- by(
  mortgage,
  list(mortgage$time),
  FUN = function(x) mean(x$payoff_fit)
)
yaxes2 <- c(
  min(payoff_obs_bytime, payoff_fit_bytime), 
  max(payoff_obs_bytime, payoff_fit_bytime)
)
plot(payoff_obs_bytime, type = "l", lty = 1, ylim = yaxes2,
     xlab = "time", ylab = "Payoff Rate and Estimated Prob Payoff")
lines(payoff_fit_bytime, lty = 2)
legend(x = "topleft", cex = 0.8, legend = c("default time", "prob of default"),
       lty = c(1,2), bty = "n", horiz = F)

# **************************************************************************
# CONTINUOUS-TIME HAZARD HAZARD MODEL
# **************************************************************************

# *****************************************************
# reshaping the data
# *****************************************************
# Prior to fitting a continuous haward model the data is reshaped  
id_last <- rep(0, length(mortgage$id))
for (i in 1:(length(mortgage$id) - 1)){
  if(mortgage$id[i+1] != mortgage$id[i]) id_last[i] = 1
}
mortgage$id_last < id_last
mortgage$time2   <- mortgage$time - mortgage$first_time + 1
mortgage3 <- subset(mortgage, id_last == 1 | status_time ==1 | status_time ==2)
mortgage3 <- mortgage3[!duplicated(mortgage3$id), ]

cross_sectional_data <- subset(mortgage3, id == 46 | id == 47 | id == 56,
                                select = c("id", "first_time", "time2", 
                                           "default_time", "payoff_time",
                                           "status_time")
)

print(cross_sectional_data)
# The table above tells us that:
# loan id 46 defaulted after 5 periods,
# loan id 47 prepaid after 3 periods,
# loan id 56 is still oustanding after 36 periods

# Now we fit a CPH model using coxph() from the library survival
# Before fitting, the function Surv() is applied to the response variable
# to clarify the condition of censoring
# The condition in this case is evaluated by the default_time
hazard <- coxph(Surv(time2, default_time) ~ FICO_orig_time + LTV_orig_time,
                data = mortgage3)
summary(hazard)