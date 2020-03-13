remove(list = ls())
mortgage <- read.csv("mortgage.csv")

#install.packages("KMsurv")#
library("KMsurv")
library("survival")

# **************************************************************************
# ACTUARIAL METHOD
# **************************************************************************

# *****************************************************
# Reshaping the data
# *****************************************************
lifetest_temp1 <- mortgage
# Sort the data by the loan id
lifetest_temp1 <- lifetest_temp1[order(lifetest_temp1$id), ]
lifetest_temp1$id[1:26]

# WE are going to record the time between the latest observation and 
# the first time the loan was recorded 
lifetest_temp1$time2     <- lifetest_temp1$time - lifetest_temp1$first_time + 1
lifetest_temp1$indicator <- rep(0, nrow(lifetest_temp1))
# return only the last record where the id is duplicated 
temp_ref                 <- which(!duplicated(lifetest_temp1$id, fromLast = T))
# and set ths to 1
lifetest_temp1$indicator[temp_ref] <- 1

# Create a new data set where we keep only the last record for those defaulted loans
lifetest_temp2           <- lifetest_temp1[lifetest_temp1$indicator == 1 |
                                   lifetest_temp1$default_time == 1,  ]

# lifetest_temp2 dataframe ONLY contains the last time when default occurs
# as an example look at the id 46
lifetest_temp2[lifetest_temp2$id == 46,]

# To develop a Survival Model, we need the time history, from first record,
# until the earlier between the last current record and the default time 
# return only the first record where the id is duplicated.
# From this only record we can read the history 
temp_ref <- which(!duplicated(lifetest_temp2$id, fromLast = F))
lifetest <- lifetest_temp2[temp_ref,]
variables <- c("id", "first_time", "time2", "default_time", 
              "payoff_time","FICO_orig_time","LTV_orig_time")

variables <- c("id", "time", "first_time", "time2", "default_time", 
               "payoff_time","FICO_orig_time","LTV_orig_time")
# As an exampe we can read the hostry for loans with id : 46, 47, 56
print(lifetest[,variables][which(lifetest$id %in% c(46, 47, 56)), ])

# I can also read the history for loan id 46, from the original dataset lifetest_temp1
lifetest_temp1[lifetest_temp1$id == 46,variables]

# *****************************************************
# Model Estimtion using lifetab
# *****************************************************
tis   <- c(seq(0,50,10), Inf)
ninit <- nrow(lifetest)
nlost <- numeric()
for (i in 1:(length(tis) - 1)){
  nlost[i] <- as.numeric(
        nrow(lifetest[lifetest$time2 >= tis[i] & 
                      lifetest$time2 <  tis[i + 1] & 
                      lifetest$default_time == 0, ]))
}

nevent <- numeric()
for (i in 1:(length(tis) - 1)){
  nevent[i] <- as.numeric(
    nrow(lifetest[lifetest$time2 >= tis[i] & 
                  lifetest$time2 < tis[i + 1] & 
                  lifetest$default_time == 1, ]))
}
LifeTable <- lifetab(tis, ninit, nlost, nevent)
print(LifeTable)

midpoints <- c(5, 15, 25,35,45)
plot(LifeTable$pdf[1:5] ~ midpoints, type = "o", ylim = c(0, 0.025),
     main = "Estimated Prob Density Function", xlab = "time2", 
     ylab = "Prob Density")

plot(LifeTable$surv[1:5] ~ midpoints, type = "o", ylim = c(0, 1),
     main = "Life Table Surv Curves", xlab = "time2", 
     ylab = "Surv Prob")

plot(LifeTable$hazard[1:5] ~ midpoints, type = "o", ylim = c(0, 0.03),
     main = "Life Table Hazard Curves", xlab = "time2", 
     ylab = "Hazard Rate")

base.model <- LifeTable$surv
# *****************************************************
# Controlling for Information in NON-PARAMETRIC Models
# *****************************************************
lifetest2 <- lifetest
lifetest2 <- lifetest2[order(lifetest2$FICO_orig_time), ]

q         <- seq(0, 1, .20)
quantiles <- numeric()
for (i in 1:6){
  quantiles[i] <- quantile(lifetest2$FICO_orig_time, q[i])
}

lifetest2$FICO_orig_time_rank <- rep (0, nrow(lifetest2))
for (i in 1:5){
  temp_ref <- which(lifetest2$FICO_orig_time %in% quantiles[i]:quantiles[i + 1])
  lifetest2$FICO_orig_time_rank[temp_ref] <- i - 1
}

lifetest2 <- lifetest2[order(lifetest2$id), ]

ranks   <- 0:4
sp_rank <- matrix(rep(0, 6*5), nrow = 6, ncol = 5)

tis   <- c(seq(0,50,10), Inf)
ninit <- as.numeric(nrow(lifetest2[lifetest2$FICO_orig_time_rank == ranks[1],]))  

nlost <- numeric()
for (i in 1:(length(tis) - 1)){
  nlost[i] <- as.numeric(
    nrow(lifetest[lifetest$time2 >= tis[i] & 
                    lifetest2$time2 < tis[i + 1] & 
                    lifetest2$default_time == 0 &
                    lifetest2$FICO_orig_time_rank == ranks[1], ]))
}
nevent <- numeric()
for (i in 1:(length(tis) - 1)){
  nevent[i] <- as.numeric(
    nrow(lifetest[lifetest$time2 >= tis[i] & 
                    lifetest2$time2 < tis[i + 1] & 
                    lifetest2$default_time == 1 &
                    lifetest2$FICO_orig_time_rank == ranks[1], ]))
}

LifeTable   <- lifetab(tis, ninit, nlost, nevent)
print(LifeTable)
midpoints   <- c(5, 15, 25,35,45)
sp_rank[,1] <- LifeTable$surv

plot(sp_rank[,1] ~ seq(0,50,10),  ylim = c(0, 1),pch = 1, typ = "o",
     main = "Life Table Surv Curves", xlab = "time2", 
     ylab = "Surv Prob")

for (j in 2:length(ranks)){
  ninit <- as.numeric(nrow(lifetest2[lifetest2$FICO_orig_time_rank == ranks[j],]))
  nlost <- numeric()
  for (i in 1:(length(tis) - 1)){
    nlost[i] <- as.numeric(
      nrow(lifetest[lifetest$time2 >= tis[i] & 
                      lifetest2$time2 < tis[i + 1] & 
                      lifetest2$default_time == 0 &
                      lifetest2$FICO_orig_time_rank == ranks[j], ]))
  }
  nevent <- numeric()
  for (i in 1:(length(tis) - 1)){
    nevent[i] <- as.numeric(
      nrow(lifetest[lifetest$time2 >= tis[i] & 
                      lifetest2$time2 < tis[i + 1] & 
                      lifetest2$default_time == 1 &
                      lifetest2$FICO_orig_time_rank == ranks[j], ]))
  }
  
  LifeTable   <- lifetab(tis, ninit, nlost, nevent)
  sp_rank[,j] <- LifeTable$surv
  
  lines(sp_rank[,j] ~ seq(0,50,10), pch = j, typ = "o")
}

# I am adding now the base model developped for all FICO_orig_time 
lines(base.model ~ seq(0,50,10), pch = 6, typ = "o", col = "red")

legend(x ="bottomleft", y.intersp = 0.9, bty = "n", cex = 0.65,
       title = "Rank of Variable FICO_orig_time", 
       legend = as.character(0:5), lwd = 1, pch = 1:6, 
       lty = 1:6, col= c(rep("black",5), "red"),  horiz = F)

# let's recheck unique values of FICO_orig_time_rank and quantile
print (unique(as.numeric(lifetest2$FICO_orig_time_rank)))
print (quantiles)

# *****************************************************
# Test Equality over Groups
# *****************************************************
# We can use the function survdiff() from the library
# survival to test  the survival curve differences
# We can run tests with and without Gehan-Wicoxon test 
survdiff (formula = Surv(time2,default_time) ~ 
  lifetest2$FICO_orig_time_rank, data = lifetest2, 
  rho = 0) # without the Gehan-Wicoxon test

survdiff (formula = Surv(time2,default_time) ~ 
  lifetest2$FICO_orig_time_rank, data = lifetest2, 
  rho = 1) # with the Gehan-Wicoxon test

# **************************************************************************
# COX PROPORTIONAL HAZARD MODELS
# **************************************************************************

# *****************************************************
# Cox Proportional Hazard Models using coxph

cph_model <-coxph(Surv(time2,default_time) ~ 
            FICO_orig_time + LTV_orig_time, 
            data = lifetest2, ties = "efron" )
summary(cph_model)

min(lifetest2$FICO_orig_time)
max(lifetest2$FICO_orig_time)
min(lifetest2$LTV_orig_time)
max(lifetest2$LTV_orig_time)

curve_high <- survfit(cph_model, newdata = 
  data.frame(FICO_orig_time = 600, LTV_orig_time = 90),
  conf.type = "none")

curve_low <- survfit(cph_model, newdata = 
  data.frame(FICO_orig_time = 800, LTV_orig_time = 60), 
  conf.type = "none")

plot(curve_high,xlim = c(0,45), xlab = "time", 
     main = "Survivor Functions",
     ylab = "Survival Probability",lty = 1)
lines(curve_low, lty = 2, add = T)
legend(x ="bottomleft", bty = "n", cex = 0.8,
       legend = c("High","Low"), lwd = 1,lty = 1:2, horiz = T)

# *****************************************************
# Time-Varying Covariates: 
# Aggregation of Time-Varying Information  
# We sort the data by id in descending order and remove
# rows with duplicated id
moment <- mortgage
moment <- moment[order(moment$id, decreasing = T),]
moment <- moment[!duplicated(moment$id),]
moment <- subset(moment, select = c("id", "LTV_time","gdp_time"))
colnames(moment) <- c("id","LTV", "gdp")
lifetest  <- lifetest[order(lifetest$id),]
lifetest2 <- merge(lifetest, moment)

# Firstly, model is run with no-time varying variables
cph_model <- coxph(Surv(time2, default_time) ~ 
                  FICO_orig_time + LTV + gdp,
                  data = lifetest2, ties = "efron")
summary(cph_model)

# now with Time-varying Covariates: Counting Process Data 
phreg       <- mortgage
phreg$time1 <- phreg$time - phreg$first_time
phreg$time2 <- phreg$time - phreg$first_time + 1
variables   <- c("id", "first_time", "time", "time1","time2", "default_time", 
               "payoff_time", "FICO_orig_time", "LTV_orig_time", "LTV_time")
id46 <- phreg[which(phreg$id == 46 & phreg$time %in% c(27:29)), variables]
id47 <- phreg[which(phreg$id == 47 & phreg$time %in% c(25:27)), variables]
id56 <- phreg[which(phreg$id == 56 & phreg$time %in% c(58:60)), variables]
print(rbind(id46, id47, id56))

cph_model <- coxph(Surv(time1, time2, default_time) ~ 
              FICO_orig_time + LTV_time + gdp_time, 
              data = phreg, ties = "efron")

summary(cph_model)
curve_downturn <- survfit(cph_model, 
          newdata = data.frame(FICO_orig_time= 800, LTV_time= 60,gdp_time =-3), 
          conf.type = "none")

curve_upturn <- survfit(cph_model, 
          newdata = data.frame(FICO_orig_time= 800, LTV_time= 60,gdp_time = 3), 
          conf.type = "none")

plot(curve_downturn, xlim = c(0,45), xlab = "time2",
     main = "Survivor Function", ylab = "Survival Probability", lty=1)
lines(curve_upturn, lty = 2, add = T)
legend(x ="bottomleft", bty = "n", cex = 1,
       legend = c("DownTurn","UpTurn"), lwd = 1, lty = 1:2, horiz = T)

# **************************************************************************
# ACCELERATED FAILURE TIME MODELS
# **************************************************************************

# *****************************************************
# Graphical Procedures

tis   <- c(seq(0,102,1), Inf)
ninit <- nrow(lifetest)

nlost <- numeric()
for (i in 1:length(tis) -1){
  nlost[i] <- as.numeric(nrow(
    lifetest[lifetest$time2 >= tis[i] & lifetest$time2 < tis[i + 1] & 
    lifetest$default_time == 0, ])) 
}

nevent <- numeric()
for (i in 1:length(tis) -1){
  nevent[i] <- as.numeric(nrow(
    lifetest[lifetest$time2 >= tis[i] & lifetest$time2 < tis[i + 1] & 
    lifetest$default_time == 1, ])) 
}

LifeTable <- lifetab(tis, ninit, nlost, nevent)
SURVIVAL  <- LifeTable$surv

plot(SURVIVAL, type = "o", main = "Estimated Survivals",
     xlim = c(0,60), xlab = "time2", ylab = "Survival probs")

plot(-log(SURVIVAL), type = "o", main = "Negative log of estimated Survivals",
     xlim = c(0,60), xlab = "time2", ylab = "-log(Survival probs)")

plot(log(-log(SURVIVAL)) ~ log(1:as.numeric(length(SURVIVAL))), type = "o", 
     main = "Negative log of estimated Survivals", xlab = "time2", 
     ylab = "Log of negative log of estimated Survival")

# *****************************************************
# Accelerated Failure Time Models with survreg

# Fitting the survival regression using the function survreg
lifereg <- survreg(Surv(time2, default_time) ~ 
                     FICO_orig_time + LTV_orig_time, 
                   data = lifetest, dist = "exponential")

summary(lifereg)
# Calibration of AFT Models:  Comparison of Default Indicators 
# and Estimated Default Probabilities
xbeta   <- predict(lifereg, type = "lp")
lambda  <- exp(-xbeta)
S1      <- 1 - pexp(phreg$time1, rate = lambda)
S2      <- 1 - pexp(phreg$time2, rate = lambda)
PD_time <- (S1 - S2)/S1
mean(phreg$default_time, na.rm = T)
mean(PD_time, na.rm = T)