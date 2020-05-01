remove(list=ls())

install.packages("aod")
install.packages("Hmisc")
install.packages("gmodels")
install.packages("ordinal")
library("aod")
library("Hmisc")
library("gmodels")
library("ordinal")

mortgage <-read.csv("mortgage.csv")
attach(mortgage)

# **************************************************************************
# LINEAR MODEL
# **************************************************************************
linearmodel <- lm(default_time ~ FICO_orig_time + LTV_orig_time + gdp_time, 
                  data = mortgage)
summary(linearmodel)
anova(linearmodel)

# **************************************************************************
# NON LINEAR MODEL - Probit, Logit, Cloglog
# **************************************************************************

x <- seq(-5,5, 0.1)
y <- pnorm(x, mean = 0, sd = 1)
ylogit   <- function(x) exp(x)/(1+exp(x))
ycloglog <- function(x) 1 - exp(-exp(x))

# Three charts together
plot(x,y,         type = "l", lwd = 2, lty = 1, xlab = "linear predictor")
curve(ylogit(x),  type = "l", lwd = 2, lty = 2, add = TRUE)
curve(ycloglog(x),type = "l", lwd = 2, lty = 3, add = TRUE)
legend("topleft", cex = 1, y.intersp = 0.8, bty= "n", 
       legend = c("probit","logit","cloglog"), lty=c(1,2,3), lwd = 2)

# **************************************************************************
# PROBIT MODEL - With Maximum Likelihood Estimation
# **************************************************************************

probit_model <- glm(default_time ~ FICO_orig_time + LTV_time + gdp_time,
                    family = binomial(link = "probit"), data = mortgage)

summary(probit_model)

# glm model automatically removes data lines (rows) where data contains na
# so to match the rows used we remove them from the mortgage data set
mortgage2 <- mortgage[complete.cases(mortgage),]
dim(mortgage2)
dim(mortgage)
probit_model_interceptonly <- glm(default_time ~ 1,
                    family = binomial(link = "probit"), data = mortgage2)

AIC_probit_model <-probit_model$aic
AIC_probit_model_interceptonly <-probit_model_interceptonly$aic

BIC_probit_model <-BIC(probit_model)
BIC_probit_model_interceptonly <-BIC(probit_model_interceptonly)

dev_probit_model <-probit_model$deviance
dev_probit_model_interceptonly <-probit_model_interceptonly$deviance

table_mfs <- as.data.frame(matrix(c(AIC_probit_model,
                                    AIC_probit_model_interceptonly,
                                    BIC_probit_model,
                                    BIC_probit_model_interceptonly,
                                    dev_probit_model,
                                    dev_probit_model_interceptonly),
                                  nrow = 3, byrow = T),
                           row.names = c("AIC","BIC","-2 LnL"))

colnames(table_mfs) <- c("Intercept and Covariates","Intercept Only")
print(table_mfs)

# For an explanation of 2 LnL, AIC and BIC tests go to the following note
# https://www.jmp.com/support/help/14-2/likelihood-aicc-and-bic.shtml
# In full, 
# the Maximum Likelihood is a technique which seeks to estimate the parameters 
# of a model, which we denote generically by β, by maximizing the likelihood function. 
# The likelihood function, denoted L(β), is the product of the probability 
# density functions (or probability mass functions for discrete distributions) 
# evaluated at the observed data values. 
# Given the observed data, maximum likelihood estimation seeks to find values 
# for the parameters, β, that maximize L(β).
# Rather than maximize the likelihood function L(β), it is more convenient 
# to work with the negative of the natural logarithm of the likelihood function, 
# -Log L(β). 
# The problem of maximizing L(β) is reformulated as a minimization problem 
# where you seek to minimize the negative log-likelihood (-LogLikelihood = -Log L(β)). 
# Therefore, smaller values of the negative log-likelihood or 
# twice the negative log-likelihood (-2LogLikelihood) indicate better model fits.
# You can use the value of negative log-likelihood to choose between models 
# and to conduct custom hypothesis tests that compare models fit using 
# different platforms in estimations 
# This is done through the use of likelihood ratio tests. 
# One reason that -2LogLikelihood is reported in many estimation platforms is that 
# the distribution of the difference between the full and reduced model 
# -2LogLikelihood values is asymptotically Chi-square. 

# McFadden's Pseudo R-Squared
# This is a link to McFadden's ratio, which is the R2 ration in logistic regrassion 
# https://thestatsgeek.com/2014/02/08/r-squared-in-logistic-regression/
1 - logLik(probit_model) / logLik(probit_model_interceptonly)

# Conduct Null Hypothesis with
# 1 - Likelihood-Ratio and its p-value
# 2 - WaldTest and its p-value

# 1 - likelihood-Ratio and p-value
LR_statistics <- probit_model$null.deviance - probit_model$deviance
pval_LR <- dchisq(LR_statistics, df = 3)

# 2 - Wald and pvalue
summary_WaldTest<-wald.test(coef(probit_model),
                             Sigma = vcov(probit_model),
                             Terms = 2:4)
pval_waldTest <- 1 - pchisq(summary_WaldTest$result$chi2[1], df = 3)
table_tgn <- as.data.frame(
  matrix(c(
    LR_statistics,
    3,
    pval_LR,
    summary_WaldTest$result$chi2[1], 
    3, pval_waldTest),
    nrow = 2,byrow = T
  ),
  row.names = c("Likelihood-Ratio","Wald")
  )
  colnames(table_tgn)<- c("Chi-Square","DF","PR > ChiSq")
print(table_tgn)

# Calibration of Probit Model
default_time_mean <- mean(mortgage$default_time, na.rm = T)
PD_time_mean      <- mean(probit_model$fitted.values)

means <- as.data.frame(c(default_time_mean ,PD_time_mean),
                       row.names = c("default_time","PD_time"))
colnames(means)<- "Means"
print(means)

# **************************************************************************
# LOGIT MODEL - With Maximum Likelihood Estimation
# **************************************************************************

logit_model <- glm(default_time ~ FICO_orig_time + LTV_time + gdp_time,
                    family = binomial(link = "logit"), data = mortgage)
summary(logit_model)

# Null Hypothesis test on each of the 4 variables of the Logit model,
# to test whether each variable is different from zero 
WaldTest0 <- (wald.test(coef(logit_model),Sigma = vcov(logit_model),
                        Terms = 1))$result$chi2
WaldTest1 <- (wald.test(coef(logit_model),Sigma = vcov(logit_model),
                        Terms = 2))$result$chi2
WaldTest2 <- (wald.test(coef(logit_model),Sigma = vcov(logit_model),
                        Terms = 3))$result$chi2
WaldTest3 <- (wald.test(coef(logit_model),Sigma = vcov(logit_model),
                        Terms = 4))$result$chi2
WaldTest  <- as.data.frame(rbind(WaldTest0, 
                                WaldTest1, 
                                WaldTest2,
                                WaldTest3),
                           row.names = 
                            c("Intercept", "FICO_orig_time","LTV_time", "gdp_time"))
print(WaldTest)

# function to prepare the Concordance Table
CONCORDANCE <- function(modelname){
  fitted           <- data.frame(cbind(modelname$y, modelname$fitted.values))
  colnames(fitted) <- c("respvar", "score")
  
  obs   <- as.numeric(length(modelname$y))
  ones  <- fitted[fitted[,1]==1, ]
  zeros <- fitted[fitted[,1]==0, ]

    pairs <- as.numeric(nrow(ones))*as.numeric(nrow(zeros))
    conc <- 0
    disc <- 0
    tied <- 0
    
    for (i in (1:nrow(ones)))
    {
      conc <- conc + sum(ones[i, "score"] >  zeros[, "score"])
      disc <- disc + sum(ones[i, "score"] <  zeros[, "score"])
      tied <- tied + sum(ones[i, "score"] == zeros[, "score"])
    }
    concordance <- (conc/pairs) *100
    discordance <- (disc/pairs) *100
    tied_perc   <- (tied/pairs) *100
    
    SOMERS.D <- round((conc - disc)/pairs,3)
    GAMMA    <- round((conc - disc)/(conc + disc),3)
    TAU.A    <- round((conc - disc)/(0.5*obs*(obs -1)),3)
    C        <- round((conc + 0.5*(pairs - conc - disc))/pairs,3)
    
    return (list("Concordance"= concordance, 
                 "Discordance"= discordance,
                 "Tied"       = tied_perc,
                 "Pairs"      = pairs,
                 "Somers.d"   = SOMERS.D,
                 "Gamma"      = GAMMA,
                 "Tau-a"      = TAU.A,
                 "c"          = C)
    )
}

CONCORD.logit_model <- CONCORDANCE(logit_model)
unlist(CONCORD.logit_model)

# **************************************************************************
# CLOGLOG MODEL - With Maximum Likelihood Estimation
# **************************************************************************

Cloglog_model <- glm(default_time ~ FICO_orig_time + LTV_time + gdp_time,
                   family = binomial(link = "cloglog"), data = mortgage)

# the glm model automatically removes na data from the dataset being analised
summary(Cloglog_model)
length(Cloglog_model$fitted.values)
length(mortgage$FICO_orig_time)

WaldTest0 <- (wald.test(coef(Cloglog_model),Sigma = vcov(Cloglog_model),
                        Terms = 1))$result$chi2
WaldTest1 <- (wald.test(coef(Cloglog_model),Sigma = vcov(Cloglog_model),
                        Terms = 2))$result$chi2
WaldTest2 <- (wald.test(coef(Cloglog_model),Sigma = vcov(Cloglog_model),
                        Terms = 3))$result$chi2
WaldTest3 <- (wald.test(coef(Cloglog_model),Sigma = vcov(Cloglog_model),
                        Terms = 4))$result$chi2
WaldTest <- as.data.frame(
  rbind(
    WaldTest0, 
    WaldTest1, 
    WaldTest2,
    WaldTest3),
    row.names = 
    c("Intercept", "FICO_orig_time","LTV_time", "gdp_time"))
print(WaldTest)

CONCORD.Cloglog_model <- CONCORDANCE(Cloglog_model)
unlist(CONCORD.Cloglog_model)

# **************************************************************************
# Qualitative Information
# **************************************************************************
# create the variable orig_time2_, which takes the value of orig_time if this is
# bewteen 20 to 25 and zero otherwise, and then we fransform it into a factor variable
# As a resut, when it used later in a probit regression, we do use its values [0, 20:25]
# but we treat it as 6 independent "factor" variables!!
(unique(mortgage$orig_time))
orig_time2  <- mortgage$orig_time
orig_time2_ <- ifelse(mortgage$orig_time %in% 20:25, mortgage$orig_time, 0)
orig_time2  <- factor(orig_time2_, levels = c(25,0,20:24))

# you could have use the two following methods:
# 1 orig_time3 <- factor(orig_time2_, levels = c(0,20:25))
# 2 orig_time4 <- factor(orig_time2_)

unique(mortgage$orig_time)%in% 20:25

probit_model_cat <- glm(default_time ~ 
      FICO_orig_time + 
      LTV_time + 
      gdp_time +
      orig_time2,
      family = binomial(link = "probit"), 
      data = mortgage)

summary(probit_model_cat)

# wald test of the model variables to "test" whether they are "statistically" different 
# from 0
library("aod") # to use the function wald.test() 
waldtest <- matrix(0, nrow = 10, ncol = 3)
for (i in 1:10){
  waldtest[i,] <- (wald.test(
    coef(probit_model_cat),
    Sigma = vcov(probit_model_cat),
    Terms = i))$result$chi2
}
waldtest <- round(waldtest,3)
print(waldtest)

waldtest <- as.data.frame(
  waldtest,row.names = names(probit_model_cat$coefficients)
  )

colnames(waldtest) <- c("Wald Chi-Square", "df","Pr>ChiSq")
print(waldtest)

# **************************************************************************
# THROUGH-THE-CYCLE (TTC) VS POINT-IN-TIME (PIT)
# **************************************************************************
# Remove all those rows where at least one column has an na value 
mortgage_omit <- na.omit(mortgage)
# Use this second method if you want to remove those rows with na ONLY in 
# columns c(7,9,10,11,17,21) 
mortgage2     <- mortgage[complete.cases(mortgage[,c(7,9,10,11,17,21)]),]

# mortgage2 is the dataset we ar using to run the next two models:
# 1 - THROUGH-THE-CYCLE (TTC), time-constant variables 
# 2 - POINT-IN-TIME (PIT)    , time-varying  variables
nrow(mortgage2)
nrow(mortgage_omit)
nrow(mortgage)

# 1
logit_TTC_with_Time <- glm(default_time ~ 
                FICO_orig_time + 
                LTV_time + 
                gdp_time + 
                time,
            family = binomial(link = "logit"), data = mortgage2)
# 2
probit_TTC_with_Time <- glm(
                default_time ~ 
                FICO_orig_time + 
                LTV_time + 
                gdp_time + 
                time,
            family = binomial(link = "probit"), data = mortgage2)
# 3
probit_TTC <- glm(
                default_time ~ 
                FICO_orig_time + 
                LTV_time + 
                gdp_time,
                family = binomial(link = "probit"), data = mortgage2)
# 4
probit_Flat <- glm(default_time ~ 
               FICO_orig_time + 
               LTV_time,
            family = binomial(link = "probit"), data = mortgage2)

# 5
probit_PIT <- glm(default_time ~ 
                FICO_orig_time + 
                LTV_time + 
                gdp_time +
                uer_time + 
                hpi_time,
                family = binomial(link = "probit"), data = mortgage2)

summary(probit_TTC_with_Time)
summary(probit_Flat)
summary(probit_TTC)
summary(probit_PIT)

PD_probit_Flat                 <- cbind(mortgage2$time, probit_Flat$fitted.values)
PD_probit_Flat_means           <- numeric()
PD_TTC_time_with_Time          <- cbind(mortgage2$time, probit_TTC_with_Time$fitted.values)
PD_TTC_time_with_Time_means    <- numeric()
PD_TTC_time                    <- cbind(mortgage2$time, probit_TTC$fitted.values)
PD_TTC_time_means              <- numeric()
PD_PIT_time                    <- cbind(mortgage2$time, probit_PIT$fitted.values)
PD_PIT_time_means              <- numeric()
Log_PD_TTC_time_with_Time      <- cbind(mortgage2$time, logit_TTC_with_Time$fitted.values)
Log_PD_TTC_time_with_Time_means<- numeric()
DR_data                        <- cbind(mortgage2$time, mortgage2$default_time)
DR                             <- numeric()

data_time <- sort(as.numeric(unique(mortgage2$time)))

for (i in 1:length(unique(mortgage2$time))){
  temp <- PD_probit_Flat[PD_probit_Flat[,1] == data_time[i],]
  PD_probit_Flat_means[i] <- mean(temp[,2], na.rm = T)
  
  temp <- Log_PD_TTC_time_with_Time[Log_PD_TTC_time_with_Time[,1] == data_time[i],]
  Log_PD_TTC_time_with_Time_means[i] <- mean(temp[,2], na.rm = T)
  
  temp <- PD_TTC_time_with_Time[PD_TTC_time_with_Time[,1] == data_time[i],]
  PD_TTC_time_with_Time_means[i] <- mean(temp[,2], na.rm = T)
  
  temp <- PD_TTC_time[PD_TTC_time[,1] == data_time[i],]
  PD_TTC_time_means[i] <- mean(temp[,2], na.rm = T)
  
  temp <- PD_PIT_time[PD_PIT_time[,1] == data_time[i],]
  PD_PIT_time_means[i] <- mean(temp[,2], na.rm = T)  

  temp <- DR_data[DR_data[,1] == data_time[i],]
  DR[i] <- sum(temp[,2], na.rm = T)/as.numeric(nrow(temp))
}

PD_TTC_time_means              <- ts(PD_TTC_time_means, frequency = 1)
plot(PD_TTC_time_means, ylim = c(0, 0.06), lty = 3, ylab = "DR and PD", col = "magenta")
PD_PIT_time_means              <- ts(PD_PIT_time_means, frequency = 1)
lines(PD_PIT_time_means, lty = 2, col = "red")
DR                             <- ts(DR, frequency = 1)
lines(DR, lty = 1, col = "black")
legend("topleft", cex = 0.8, y.intersp = 1, bty = "n", 
       legend = c("PD_TTC_time","PD_PIT_time", "Default_time"),
       lty = c(3,2,1), lwd = 2, 
       col = c("magenta", "red","black"))

PD_TTC_time_with_Time_means    <- ts(PD_TTC_time_with_Time_means, frequency = 1)
plot(PD_TTC_time_with_Time_means, ylim = c(0, 0.06), lty = 4, ylab = "DR and PD", 
     col = "red")
Log_PD_TTC_time_with_Time_means<- ts(Log_PD_TTC_time_with_Time_means, frequency = 1)
lines(log_PD_TTC_time_with_Time_means, lty = 5, col = "blue")
PD_probit_Flat_means           <- ts(PD_probit_Flat_means, frequency = 1)
lines(PD_probit_Flat_means, lty = 6, col = "green")
DR                             <- ts(DR, frequency = 1)
lines(DR, lty = 1, col = "black")
legend("topleft", cex = 0.6, y.intersp = 1, bty = "n", 
       legend = c("PD_TTC_time_with_Time", "log_PD_TTC_time_with_Time",
                  "PD_probit_Flat", "Default_time"),lty = c(4:6, 1), lwd = 2, 
                  col = c("red","blue", "green", "black"))

# **************************************************************************
# ESTIMATION OF RATING MIGRATION PROBABILITIES
# **************************************************************************
# create rating classes according to the FICO_orig_time variable
mortgage2 <- mortgage[order(mortgage$id,mortgage$time), ]
rating    <- rep(1, length(mortgage2$FICO_orig_time))
rating[mortgage2$FICO_orig_time > 350 ]<- 1
rating[mortgage2$FICO_orig_time > 500 ]<- 2
rating[mortgage2$FICO_orig_time > 650 ]<- 3

library("Hmisc")
# shift the vector mortgage2$id and rating by one element 
lagid     <- Lag(mortgage2$id,1)
lagrating <- Lag(rating,1)

rating[mortgage2$default_time == 1]    <- 4

# Now I have created an history of rating migration for the same loan id.
# when the loan id is not the same, I set it to NA
lagrating[(mortgage2$id != lagid)] = NA

library("gmodels")
CrossTable(lagrating, rating, digits = 4, prop.c = F, prop.t = F, prop.chisq = F)

# Preparing for a probit model
rating    <- factor(rating, c(as.character(1:4)))
lagrating <- factor(lagrating)
contrasts(lagrating) <- contr.sum(3) # the book says contr.sum(3), I think is incorrect 
contrasts(lagrating)

# fitting a cumulative probit model
ordered_probit <- clm(rating ~ lagrating, link = "probit")
summary(ordered_probit)

# fitting a cumulative probit model
contrasts(lagrating) <- contr.sum(3)
contrasts(lagrating)
ordered_probit <- clm(rating ~ lagrating + mortgage2$gdp_time, link = "probit")
summary(ordered_probit)

# **************************************************************************
# Variable Selection
# **************************************************************************
# performing t-test on fico_orig_time with 2 methods:
# 1 Pooled and 2 Satterthwaite
tt_pooled <- t.test(mortgage$FICO_orig_time ~ mortgage$default_time, var.equal = T)
tt_sat    <- t.test(mortgage$FICO_orig_time ~ mortgage$default_time, var.equal = F)

# ***************************************************************
# SUPER IMPORTANT tests on the FICO score assciated to default/no default
print(tt_pooled)
print(tt_sat)

values <- c(
  round(tt_pooled$parameter,1), round(tt_pooled$statistic,2), tt_pooled$p.value,
  round(tt_sat$parameter,1),    round(tt_sat$statistic,2), tt_sat$p.value )
results <- matrix(round(values,4), nrow = 2, byrow = T)
results <- as.data.frame(cbind(c("Equal","Unequal"), results), 
                         row.names = c("Pooled", "Satterthwaite"))
colnames(results) <- c("Variances","DF","t value", "Pr > |t|")
print(results)
# ***************************************************************

bounds           <- cbind(c(NA, 488, 576, 664, 752),c(488,576,664,752,NA))
colnames(bounds) <- c("lower", "upper")
bounds
non_event_count <- function(lower, upper){
  if(missing(lower)){
    return (length(mortgage$FICO_orig_time[mortgage$FICO_orig_time < upper &
            mortgage$default_time == 0]))
  }
  else if(missing(upper)){
      return (length(mortgage$FICO_orig_time[mortgage$FICO_orig_time > lower &
            mortgage$default_time == 0]))
  }else{
      return (length(mortgage$FICO_orig_time[mortgage$FICO_orig_time >= lower &
          mortgage$FICO_orig_time < upper & mortgage$default_time == 0]))
  }
}

event_count <- function(lower, upper){
  if(missing(lower)){
    return (length(mortgage$FICO_orig_time[mortgage$FICO_orig_time < upper &
            mortgage$default_time == 1]))
  }
  else if(missing(upper)){
      return (length(mortgage$FICO_orig_time[mortgage$FICO_orig_time > lower &
            mortgage$default_time == 1]))
  }else{
      return (length(mortgage$FICO_orig_time[mortgage$FICO_orig_time >= lower &
          mortgage$FICO_orig_time < upper & mortgage$default_time == 1]))
  }
}
noneventcount <- numeric()
for (i in 1:length(bounds[,1])){
  noneventcount[i] <- non_event_count(lower = bounds[i,1], upper = bounds[i,2])
}
eventcount    <- numeric()
for (i in 1:length(bounds[,1])){
  eventcount[i] <- event_count(lower = bounds[i,1], upper = bounds[i,2])
}

total        <- apply(cbind(noneventcount, eventcount),1,sum)
cbind(noneventcount,eventcount,total)

noneventrate <- noneventcount/total
eventrate    <- eventcount/total

cbind(noneventrate,eventrate,noneventrate+eventrate)

Total_nonevent <- sum(noneventcount, na.rm = T)
Total_event    <- sum(eventcount, na.rm = T)

# compute the weight of events for each set of bounds 
woe <- numeric()
for (i in 1:length(bounds[,1])){
  woe[i]<- log((noneventcount[i]/Total_nonevent) / (eventcount[i]/ Total_event))
}
# prepare the information value for each set of bounds
iv <- numeric()
for (i in 1:length(bounds[,1])){
  iv[i]<- sum((noneventcount[i]/Total_nonevent) - (eventcount[i]/ Total_event))*woe[i] 
}

table_hpbin <- as.data.frame(
    matrix(
      c(
        bounds[,1], bounds[,2], noneventcount, round(noneventrate,3),
        eventcount, round(eventrate,3), round(woe,3), round(iv,3)
        ),
      ncol = 8)
    )
colnames(table_hpbin) <- c("Lower Bound", "Upper Bound", "Non_Event Count",
                           "Non_Event Rate", "Event Count","Evnt Rate", 
                           "Weight of Evidence", "Information Value")
print(table_hpbin)

# the FICO score information value is the sum of the column "information value"  
information_value <- sum(iv, na.rm = T)
print(information_value)

# **************************************************************************
# FITTING AND FORECASTING
# **************************************************************************
# remove missing values from the dataset, 
# set the random seed to a value, so that the experiment cen be replicated,
# create two stratified samples: one with 80% an another with 20% of original
# data
mortgage2 <- mortgage[complete.cases(mortgage),]
set.seed(12345)
ref               <- sample(nrow(mortgage2), dim(mortgage2)[1]*.8, replace=F)
mortgage2_sample1 <- mortgage2[ref,]
mortgage2_sample0 <- mortgage2[-ref,]

model_logistic <- glm(default_time ~ 
          FICO_orig_time + LTV_time + gdp_time,
          family = binomial("logit"), data= mortgage2_sample1)

model_logistic_Flat <- glm(default_time ~ 
          FICO_orig_time + LTV_orig_time,
          family = binomial("logit"), data= mortgage2_sample1)

PD1     <- cbind(mortgage2_sample1$time, model_logistic$fitted.values)
PDFlat1 <- cbind(mortgage2_sample1$time, model_logistic_Flat$fitted.values)

# prepare three tapply tables, where the analysis is run by the variable time,
# and we include the defaut time and the fitted values of the two probit models
default_rate1 <- tapply(mortgage2_sample1$default_time, mortgage2_sample1$time,
                        function(x) sum(x, na.rm = T)/length(x))
PD_time1      <- tapply(PD1[,2], PD1[,1], 
                        function(x) sum(x, na.rm = T)/ length(x))
PDFlat_time1  <- tapply(PDFlat1[,2], PDFlat1[,1], 
                        function(x) sum(x, na.rm = T)/ length(x))

default_rate1 <- as.ts(default_rate1)
PD_time1      <- as.ts(PD_time1)
PDFlat_time1  <- as.ts(PDFlat_time1)

plot(default_rate1, ylim = c(0, max(PD_time1)), lty = 1, ylab = "DR and PD")
lines(PD_time1,     lty = 2)
lines(PDFlat_time1, lty = 3)
legend("topright", cex = 0.8, y.intersp = 1, bty = "n",
     legend = c("Default Time","PD time", "PD time Flat"), lty = 1:3, lwd = 1)

# Outsample analysis on "mortgage2_sample0"
# we are now testing the model of the data excluded from the fitting, that 20%

PD0 <- predict.glm(model_logistic, newdata = mortgage2_sample0, type = "response")
PD0 <- cbind(mortgage2_sample0$time, PD0)

default_rate0 <- tapply(mortgage2_sample0$default_time, mortgage2_sample0$time,
                        function(x) sum(x, na.rm = T)/length(x))
PD_time0      <- tapply(PD0[,2], PD0[,1], 
                        function(x) sum(x, na.rm = T)/ length(x))

plot(default_rate0, ylim = c(0, max(PD_time0)), lty = 1, ylab = "DR and PD")
lines(PD_time0, lty = 2)
legend("topright", cex = 0.8, y.intersp = 1, bty = "n",
     legend = c("Default Time OutSample","PD Time"), lty = 1:2, lwd = 1)

# **************************************************************************
# FORMATION OF RATING CLASSES
# **************************************************************************
# Rating Approach 1
# numbers from [0,9] according to FICO rating Bands
FICO_orig_time_rank1 <- rep(0, nrow(mortgage2))
FICO_orig_time_rank1[mortgage2$FICO_orig_time >= 300 & 
                       mortgage2$FICO_orig_time <400]<- 0
FICO_orig_time_rank1[mortgage2$FICO_orig_time >= 400 & 
                       mortgage2$FICO_orig_time <450]<- 1
FICO_orig_time_rank1[mortgage2$FICO_orig_time >= 450 & 
                       mortgage2$FICO_orig_time <500]<- 2
FICO_orig_time_rank1[mortgage2$FICO_orig_time >= 500 & 
                       mortgage2$FICO_orig_time <550]<- 3
FICO_orig_time_rank1[mortgage2$FICO_orig_time >= 550 & 
                       mortgage2$FICO_orig_time <600]<- 4
FICO_orig_time_rank1[mortgage2$FICO_orig_time >= 600 & 
                       mortgage2$FICO_orig_time <650]<- 5
FICO_orig_time_rank1[mortgage2$FICO_orig_time >= 650 & 
                       mortgage2$FICO_orig_time <700]<- 6
FICO_orig_time_rank1[mortgage2$FICO_orig_time >= 700 & 
                       mortgage2$FICO_orig_time <750]<- 7
FICO_orig_time_rank1[mortgage2$FICO_orig_time >= 750 & 
                       mortgage2$FICO_orig_time <800]<- 8
FICO_orig_time_rank1[mortgage2$FICO_orig_time >= 800 & 
                       mortgage2$FICO_orig_time <850]<- 9
min(mortgage2$FICO_orig_time)
max(mortgage2$FICO_orig_time)

# Rating Approach 2 : on all data set
q <- seq(0,1,.1)
quantiles <- numeric()
for (i in 1:length(q)){
  quantiles[i] <- quantile(mortgage2$FICO_orig_time, q[i])
}
quantile(mortgage2$FICO_orig_time, q[1])
quantile(mortgage2$FICO_orig_time, q[2])
quantile(mortgage2$FICO_orig_time, q[10])
quantile(mortgage2$FICO_orig_time, q[11])

# numbers from [0,9] Based on Quantile of FICO 
FICO_orig_time_rank2 <- FICO_orig_time_rank1 
FICO_orig_time_rank2 <- rep(0, nrow(mortgage2))
for (i in 0:9){
  FICO_orig_time_rank2[mortgage2$FICO_orig_time %in% quantiles[i+1]:quantiles[i+2]] <- i
}

# create a copy of mortgage2 dataset into rank3
rank3     <- mortgage2
# in rank3 add the two previous ranks, 1 and 2 
rank3$FICO_orig_time_rank1<- FICO_orig_time_rank1
rank3$FICO_orig_time_rank2<- FICO_orig_time_rank2
# Order rank3 according to default time
# This will move the data with default time = 0 to the top, 
# and the data with default time = 1 to the bottom
rank3     <- rank3[order(rank3$default_time),]
rank3_dt0 <- rank3[rank3$default_time == 0,]
rank3_dt1 <- rank3[rank3$default_time == 1,]

# Rating Approach 3: 
# unlike method 2, the are now two sets of quantile being prepared:
# ONLY on the data that did not default
q <- seq(0,1,.1)
quantiles <- numeric()
for (i in 1:length(q)){
  quantiles[i] <- quantile(rank3_dt0$FICO_orig_time, q[i])
}

FICO_orig_time_rank3_0 <- rep(0,nrow( rank3_dt0) )
for (i in 0:9){
  FICO_orig_time_rank3_0[rank3_dt0$FICO_orig_time %in% quantiles[i+1]:quantiles[i+2]]<-i
}

# ONLY on the data that did default
q <- seq(0,1,.1)
quantiles <- numeric()
for (i in 1:length(q)){
  quantiles[i] <- quantile(rank3_dt1$FICO_orig_time, q[i])
}

FICO_orig_time_rank3_1 <- rep(0, dim( rank3_dt1)[1] )
for (i in 0:9){
  FICO_orig_time_rank3_1[rank3_dt1$FICO_orig_time %in% quantiles[i+1]:quantiles[i+2]]<-i
}

dim(rank3)
rank3$FICO_orig_time_rank3 <- c(FICO_orig_time_rank3_0, FICO_orig_time_rank3_1)
dim(rank3)

#order rank3 in ascending order by FICO_orig_time and descending order by default_time
rank3 <- rank3[order(rank3$FICO_orig_time, rank3$default_time, decreasing = c(F,T)), ]
# find now the row where the defult occurs first
first_default_ref <- which(rank3[-1,]$default_time > rank3[-nrow(rank3),]$default_time)+1 

rank3$new_FICO_score <- rep(0,nrow(rank3))
rank3$new_FICO_score[first_default_ref] <- rank3$FICO_orig_time_rank3[first_default_ref]

for(i in 1:length(first_default_ref)){
  if (i == length(first_default_ref)){
    temp_ref <- (first_default_ref[i] + 1):nrow(rank3)
  }
  else{
    temp_ref <- (first_default_ref[i] + 1):first_default_ref[i + 1]
  }
  temp <- rank3$new_FICO_score[first_default_ref[i]]
  rank3$FICO_orig_time_rank3[temp_ref][rank3$default_time[temp_ref] == 0]<- temp
}

rank3         <- rank3[order(rank3$FICO_orig_time_rank1), ]
# dataframe for plotting purposes
orank1        <- data.frame(matrix(0,nrow = 9))
orank1$N_1    <- tapply(rank3$default_time,  rank3$FICO_orig_time_rank1, length)
orank1$DR_1   <- tapply(rank3$default_time,  rank3$FICO_orig_time_rank1, mean)
orank1$FICO_1 <- tapply(rank3$FICO_orig_time,rank3$FICO_orig_time_rank1, mean)
orank1 <- orank1[,2:4]

rank3         <- rank3[order(rank3$FICO_orig_time_rank2), ]
# dataframe for plotting purposes
orank2        <- data.frame(matrix(0,nrow = 10))
orank2$N_2    <- tapply(rank3$default_time,  rank3$FICO_orig_time_rank2, length)
orank2$DR_2   <- tapply(rank3$default_time,  rank3$FICO_orig_time_rank2, mean)
orank2$FICO_2 <- tapply(rank3$FICO_orig_time,rank3$FICO_orig_time_rank2, mean)
orank2 <- orank2[,2:4]

rank3         <- rank3[order(rank3$FICO_orig_time_rank3), ]
# dataframe for plotting purposes
orank3        <- data.frame(matrix(0,nrow = 10))
orank3$N_3    <- tapply(rank3$default_time,  rank3$FICO_orig_time_rank3, length)
orank3$DR_3   <- tapply(rank3$default_time,  rank3$FICO_orig_time_rank3, mean)
orank3$FICO_3 <- tapply(rank3$FICO_orig_time,rank3$FICO_orig_time_rank3, mean)
orank3 <- orank3[,2:4]

orank1$N_1 <- orank1$N_1/nrow(rank3) 
orank2$N_2 <- orank2$N_2/nrow(rank3) 
orank3$N_3 <- orank3$N_3/nrow(rank3) 

plot(orank1$FICO_1, orank1$N_1, xlab = "Average Fico Score", main = "Frequency of loans per FICO bands",
     ylab = "Relative Freq", xlim = c(300, 1000), ylim = c(0, 0.3),type = "l", lty = 1)
lines(orank2$FICO_2, orank2$N_2, lty = 2)
lines(orank3$FICO_3, orank3$N_3, lty = 3)
legend("topright", cex = 0.8, y.intersp = 1, bty = "n", 
      legend = c("approach1", "approach2", "approach3"), lty = 1:3, lwd = 1)

plot(orank1$FICO_1, orank1$DR_1, xlab = "Average Fico Score", main = "Def rate per FICO bands ",
     ylab = "Relative Freq", xlim = c(300, 1000), ylim = c(0, max(orank1$DR_1))
     ,type = "l", lty = 1)
lines(orank2$FICO_2, orank2$DR_2, lty = 2)
lines(orank3$FICO_3, orank3$DR_3, lty = 3)
legend("topright", cex = 0.8, y.intersp = 1, bty = "n",
      legend = c("approach1", "approach2", "approach3"), lty = 1:3, lwd = 1)
