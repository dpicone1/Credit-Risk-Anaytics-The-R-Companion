remove(list=ls())
mortgage <- read.csv("mortgage.csv", header = T)
lgd      <- read.csv("lgd.csv", header = T) 

install.packages("pROC")
install.packages("Hmisc")
install.packages("ResourceSelection")

library("pROC")
library("Hmisc")
library("ResourceSelection")

# **************************************************************************
# BACKTESTING PD MODELS
# **************************************************************************

# *****************************************************
# Backtesting PD at level 0
# *****************************************************

score_expected_actual     <- matrix(nrow = 10, ncol = 2)
score_expected_actual[,1] <- c(0.06,.1,.09,.12,.12,.08, .07,.08, .12, .16)
score_expected_actual[,2] <- c(.07,.08,.07,.09,.11,.11,.1,.12,.11,.15)
print(score_expected_actual)

PSI_row <- (score_expected_actual[,2] - score_expected_actual[,1]) * 
  (log(score_expected_actual[,2]) - log(score_expected_actual[,1]))
PSI <- sum(PSI_row)
print(PSI_row)
print(PSI)

score_expected_actual <- cbind(score_expected_actual,
                               c(.06,.07,.1,.11,.1,.09,.11,.11,.1,.15))
PSI_row_e0 <- PSI_row
PSI_e0 <- PSI
print(PSI_row_e0)
print(PSI_e0)

PSI_row_e1 <- (score_expected_actual[,3] - score_expected_actual[,1]) * 
  (log(score_expected_actual[,3]) - log(score_expected_actual[,1]))
PSI_e1 <- sum(PSI_row_e1)
print(PSI_row_e1)
print(PSI_e1)


PSI_row_t <- (score_expected_actual[,3] - score_expected_actual[,2]) * 
  (log(score_expected_actual[,3]) - log(score_expected_actual[,2]))
PSI_t <- sum(PSI_row_t)
print(PSI_row_t)
print(PSI_t)

time_dummy <- rep(0, length(mortgage$time))
len <- length(mortgage$time)
for (i in 1:len){
  if(mortgage$time[i] > 59) time_dummy[i] = 1
}

model_probit <- glm(default_time ~ FICO_orig_time + LTV_orig_time + gdp_time +
                      time_dummy + 
                      FICO_orig_time * time_dummy +
                      LTV_orig_time  * time_dummy,
                    family = binomial(link = "probit"),data = mortgage)

summary(model_probit)


model_probit_interceptonly <- glm(default_time ~ 1,
                    family = binomial(link = "probit"),data = mortgage)

summary(model_probit_interceptonly)

AIC_model_probit               <- model_probit$aic
AIC_model_probit_interceptonly <- model_probit_interceptonly$aic

BIC_model_probit               <- BIC(model_probit)
BIC_model_probit_interceptonly <- BIC(model_probit_interceptonly)

dev_model_probit               <- model_probit$deviance
dev_model_probit_interceptonly <- model_probit_interceptonly$deviance

table.mfs <- as.data.frame(matrix(c(
  AIC_model_probit, 
  AIC_model_probit_interceptonly,
  BIC_model_probit,
  BIC_model_probit_interceptonly, 
  dev_model_probit,
  dev_model_probit_interceptonly
  ),
  nrow = 3, byrow = T),
  row.names = c("AIC","BIC","-2 ln L")
  )
colnames(table.mfs) <- c("Intercept Only", "Intercept and Covariates")
print(table.mfs)

CONCORDANCE <- function(modelname){
    fitted <- data.frame(cbind(modelname$y,modelname$fitted.values))
    colnames(fitted) <- c("respvar", "score")
    obs   <- as.numeric(length(modelname$y))
    ones  <- fitted[fitted[,1]==1, ]
    zeros <- fitted[fitted[,1]==0, ]
    pairs <- as.numeric(nrow(ones))*as.numeric(nrow(zeros))
    conc <- 0
    disc <- 0
    tied <- 0
    
    for (i in 1: nrow(ones))
    {
      conc <- conc + sum(ones[i, "score"] >  zeros[, "score"])
      disc <- disc + sum(ones[i, "score"] <  zeros[, "score"])
      tied <- tied + sum(ones[i, "score"] == zeros[, "score"])
    }
    concordance <- (conc/pairs) *100
    discordance <- (disc/pairs) *100
    tied_perc   <- (tied/pairs) *100
    
    SOMERS.D <- round((conc - disc)/pairs,3)
    GAMMA <- round((conc - disc)/(conc + disc),3)
    TAU.A <- round((conc - disc)/(0.5*obs*(obs -1)),3)
    C     <- round((conc + 0.5*(pairs - conc - disc))/pairs,3)
    
    return (list("Concordance"= concordance, 
                 "Discordance"= discordance,
                 "Tied"  = tied_perc,
                 "Pairs" = pairs,
                 "Somers.d" = SOMERS.D,
                 "Gamma" = GAMMA,
                 "Tau-a" = TAU.A,
                 "c" = C)
    )
  }
unlist(CONCORDANCE(model_probit))

# ***********************************************************************
# BACKTESTING PD AT LEVEL 1
# ***********************************************************************
# Two datasets: with time <60 and with time >59
tmp_pdvali1 <- mortgage[mortgage$time < 60, ]
tmp_pdvali2 <- mortgage[mortgage$time > 59, ]

# Running three glm - probit models on the dataset with time < 60
model_probit1 <- glm(default_time ~ LTV_orig_time ,
                     family = binomial(link = "probit"),data = tmp_pdvali1)
model_probit2 <- glm(default_time ~ FICO_orig_time + LTV_orig_time ,
                     family = binomial(link = "probit"),data = tmp_pdvali1)
model_probit3 <- glm(default_time ~ FICO_orig_time + LTV_orig_time + gdp_time,
                     family = binomial(link = "probit"),data = tmp_pdvali1)


model_probit1_predict <- predict.glm(model_probit1,newdata =tmp_pdvali2, 
                                     type = "response")
model_probit2_predict <- predict.glm(model_probit2,newdata =tmp_pdvali2, 
                                     type = "response")
model_probit3_predict <- predict.glm(model_probit3,newdata =tmp_pdvali2, 
                                     type = "response")
obs <- tmp_pdvali2$default_time
# calculate the ROC for three models above vs actual default with time > 59
roc1 <- roc(obs,model_probit1_predict)
roc2 <- roc(obs,model_probit2_predict)
roc3 <- roc(obs,model_probit3_predict)

area    <- c(auc(roc1), auc(roc2), auc(roc3))
auc_ci  <- rbind(ci(roc1)[c(1,3)], ci(roc2)[c(1,3)],ci(roc3)[c(1,3)])
auc_se  <- rbind(sqrt(var(roc1)), sqrt(var(roc2)),sqrt(var(roc3)))
somersD <- 2*area -1

table_roc           <- as.data.frame(cbind(area, auc_se, auc_ci, somersD))
colnames(table_roc) <- c("Area","SE","95% CI L","95% CI R", "Somers'D")
print(table_roc)

roc.test(roc1, roc2, method = "delong")
roc.test(roc1, roc3, method = "delong")
roc.test(roc2, roc3, method = "delong")

plot(roc1, lty = 1, xlim = c(1.0,0.0), ylim = c(0,1), col = "red")
lines(roc2, lty = 2, col = "blue");lines(roc3, lty = 3, col = "magenta")
legend(x = "topleft", lty= 1:3, col = c("red", "blue", "magenta"), 
       legend=paste0("roc",1:3), horiz = F, cex=1, bty = "n")

# ***********************************************************************
# BACKTESTING PD AT LEVEL 2
# ***********************************************************************
# Brier scores
brier_row_1 <- (model_probit1_predict - obs) ^2
brier_row_2 <- (model_probit2_predict - obs) ^2
brier_row_3 <- (model_probit3_predict - obs) ^2
brier_1 <- sum(brier_row_1) / length(obs)
brier_2 <- sum(brier_row_2) / length(obs)
brier_3 <- sum(brier_row_3) / length(obs)
print(c(brier_1, brier_2, brier_3))

len         <- length(tmp_pdvali1$FICO_orig_time)
fico_class1 <- rep(0, len)

# on the dataset with time < 60, create three fico classes as per below
# *********************************************************************
for (i in (1:len)){
  if(tmp_pdvali1$FICO_orig_time[i] >= 713){
    fico_class1[i] = 1
  }else if(tmp_pdvali1$FICO_orig_time[i] >= 648){
    fico_class1[i] = 2
  }else{
    fico_class1[i] = 3
  }
}
# alternatively, instead of running a loop which takes long time
# one can use the following:
# 1 - ifelse(,,) 
fico_class12 <- ifelse(tmp_pdvali1$FICO_orig_time >= 713,
                      1,
                      ifelse(tmp_pdvali1$FICO_orig_time >= 648,
                      2,
                      3))
# 2 - a function together with sapply()
f <- function(x)
  {
    res <- 0
    if(x >= 713){
      res = 1
    }else if(x >= 648){
      res = 2
    }else{
      res = 3
    }
    return (res)
}
sum(fico_class1)
sum(fico_class12)
sum(sapply(tmp_pdvali1$FICO_orig_time, f))
# ********************************************************************* 
fico_class2 <- 
  ifelse(tmp_pdvali2$FICO_orig_time >= 713,1,
         ifelse(tmp_pdvali2$FICO_orig_time >= 648,2,3))

table_insample       <- table(fico_class1, tmp_pdvali1$default_time)
table_insample_prop  <- prop.table(table_insample, 1)
print (table_insample)
print (table_insample_prop)

table_outsample      <- table(fico_class2, tmp_pdvali2$default_time)
table_outsample_prop <- prop.table(table_outsample, 1)
print (table_outsample)
print (table_outsample_prop)

print (table_outsample)
# bin test on class 1 
binconf(
  x = table_outsample[1,2],
  n = sum(table_outsample[1,]),
  alpha = 0.05,
  method = "all"
)

binom.test(
  x = table_outsample[1,2],
  n = sum(table_outsample[1,]),
  p = 0.01)

# bin test on class 2
binconf(
  x = table_outsample[2,2],
  n = sum(table_outsample[2,]),
  alpha = 0.05,
  method = "all"
)

binom.test(
  x = table_outsample[2,2],
  n = sum(table_outsample[2,]),
  p = 0.01)

# bin test on class 3
binconf(
  x = table_outsample[3,2],
  n = sum(table_outsample[3,]),
  alpha = 0.05,
  method = "all"
)

binom.test(
  x = table_outsample[3,2],
  n = sum(table_outsample[3,]),
  p = 0.01)

# Hosmer-Lemeshow test on the prediction from model_probit3
hl <- hoslem.test(x = obs, y = model_probit3_predict, g = 10)
cbind(hl$observed, hl$expected)

# Running model_probit3 on the dataset with time > 59
xbeta3    <- predict.glm(model_probit3, newdata = tmp_pdvali2, type = "link")
model_HLS <- glm(default_time ~ xbeta3,
                 family = binomial(link = "logit"),
                 data = tmp_pdvali2)

model_HLS_predit <- predict.glm(model_HLS, type = "response")

# Hosmer-Lemeshow test on the prediction from model_probit3
# with object hl()
hl2 <- hoslem.test(x = obs, y = model_HLS_predit, g = 10)
cbind(hl2$observed, hl2$expected)

obs_freq <- hl$observed[,2] / apply(hl$observed,1,sum)
exp_freq <- hl$expected[,2] / apply(hl$observed,1,sum)
plot(obs_freq, exp_freq, xlim = c(0, max(obs_freq)), ylim = c(0, max(exp_freq)))
abline(0,1,lty = 2)

# Hosmer-Lemeshow test on the prediction from model_probit3
# with object hl2()
obs_freq2 <- hl2$observed[,2] / apply(hl2$observed,1,sum)
exp_freq2 <- hl2$expected[,2] / apply(hl2$observed,1,sum)
plot(obs_freq2, exp_freq2, xlim = c(0, max(obs_freq2)), ylim = c(0, max(exp_freq2)))
abline(0,1,lty = 2)

# Generate sequences
N   <- 100
rho <- seq(0, 0.3,   by = 0.05)
p   <- seq(.005, .1, by = 0.005)

Myfn <- function(k,p,N,rho){
  cpd <- function(x,p,rho){
    pnorm((1/sqrt(1-rho))*(qnorm(p) - sqrt(rho)*x))
  }
  v <- function(x){
    pbinom(k,N,cpd(x,p,rho))*dnorm(x)
  }
  return (v)
}

default_rate <- matrix(nrow = length(p), ncol = length(rho))
for (j in (1:length(rho))){
  for (i in (1:length(p))){
    z <- 0
    k <- 0
    while (z <.99){
      z <- integrate(Myfn(k = k, p = p[i], N = N, rho = rho[j]),
                     lower = -Inf, upper = Inf)$value
      k <- k + 1
    }
    default_rate[i,j] <- k/N
  }
}
plot(0, xlab = "PD", ylab = "Crit Values", xlim =c(0, .1),ylim = c(0, .6))
for (j in (1:length(rho))){
  lines(p, default_rate[,j], lty = j)
}
legend("topleft", legend = rho, lty = seq(1:length(rho)), bty = "n", cex = .7, horiz = F)


p_crit <- matrix(nrow = length(p), ncol = length(rho))
plot(0, xlab = "PD", ylab = "Crit Values", xlim =c(0, .1),ylim = c(0, .6))
for (i in (1:length(rho))){
  p_crit[ ,i] <- pnorm((qnorm(p) + sqrt(rho[i])*qnorm(0.99))/sqrt(1-rho[i]))
  lines(p, p_crit[,i], lty = i)
}
legend("topleft", legend = rho, lty = seq(1:length(rho)), bty = "n", cex = .7, horiz = F)


alpha<-0.99
n    <- 1000
pd0  <- seq(0.005, 0.1, by = 0.005)
pd1  <- pd0 + 0.005
pd2  <- pd0 + 0.01
pd3  <- pd0 + 0.025
pd4  <- pd0 + 0.05
p_beta_005 <- pnorm((pd0 - pd1 + qnorm(alpha) * sqrt(pd0 * (1-pd0)/n))/
                      (sqrt(pd1*(1-pd1)/n)))
p_beta_01  <- pnorm((pd0 - pd2 + qnorm(alpha) * sqrt(pd0 * (1-pd0)/n))/
                      (sqrt(pd2*(1-pd2)/n)))
p_beta_025 <- pnorm((pd0 - pd3 + qnorm(alpha) * sqrt(pd0 * (1-pd0)/n))/
                      (sqrt(pd3*(1-pd3)/n)))
p_beta_05  <- pnorm((pd0 - pd4 + qnorm(alpha) * sqrt(pd0 * (1-pd0)/n))/
                      (sqrt(pd4*(1-pd4)/n)))

plot(0, xlab = "PD0", xlim = c(0,.1), ylab = "P_Beta", ylim = c(0,1))
lines(pd0, p_beta_005, lty = 1)
lines(pd0, p_beta_01, lty = 2)
lines(pd0, p_beta_025, lty = 3)
lines(pd0, p_beta_05, lty = 4)
legend("bottomright", 
       legend = c("p_beta_005","p_beta_01","p_beta_025","p_beta_05"), 
       lty = 1:4, bty = "n", cex = .7, horiz = F)

alpha <-0.99
pd0   <- seq(0.005, 0.1, by = 0.005)
pd1   <- pd0 + 0.05
rho1  <- 0.01
rho2  <- 0.05
rho3  <- 0.1
rho4  <- 0.2
rho5  <- 0.3

p_crit_rho01 <- pnorm((qnorm(pd0) + sqrt(rho1) * qnorm(alpha)) / sqrt(1- rho1))
p_beta_rho01 <- pnorm((sqrt(1-rho1) * qnorm(p_crit_rho01) - qnorm(pd1))/sqrt(rho1))

p_crit_rho05 <- pnorm((qnorm(pd0) + sqrt(rho2) * qnorm(alpha)) / sqrt(1- rho2))
p_beta_rho05 <- pnorm((sqrt(1-rho2) * qnorm(p_crit_rho05) - qnorm(pd1))/sqrt(rho2))

p_crit_rho10 <- pnorm((qnorm(pd0) + sqrt(rho3) * qnorm(alpha)) / sqrt(1- rho3))
p_beta_rho10 <- pnorm((sqrt(1-rho3) * qnorm(p_crit_rho10) - qnorm(pd1))/sqrt(rho3))

p_crit_rho20 <- pnorm((qnorm(pd0) + sqrt(rho4) * qnorm(alpha)) / sqrt(1- rho4))
p_beta_rho20 <- pnorm((sqrt(1-rho4) * qnorm(p_crit_rho20) - qnorm(pd1))/sqrt(rho4))

p_crit_rho30 <- pnorm((qnorm(pd0) + sqrt(rho5) * qnorm(alpha)) / sqrt(1- rho5))
p_beta_rho30 <- pnorm((sqrt(1-rho5) * qnorm(p_crit_rho30) - qnorm(pd1))/sqrt(rho5))

plot(0, xlab = "PD0", xlim = c(0,.1), ylab = "P_Beta_Rho", ylim = c(0,1))
lines(pd0, p_beta_rho01, lty = 1)
lines(pd0, p_beta_rho05, lty = 2)
lines(pd0, p_beta_rho10, lty = 3)
lines(pd0, p_beta_rho20, lty = 4)
lines(pd0, p_beta_rho30, lty = 5)
legend("bottomright", 
       legend = c("p_beta_rho01","p_beta_rho05","p_beta_rho10","p_beta_rho20","p_beta_rho30"), 
       lty = 1:5, bty = "n", cex = .7, horiz = F)

# ***********************************************************************
# BACKTESTING LGD AND EAD
# ***********************************************************************

# *****************************************************
# BACKTESTING AT LEVEL 1
# *****************************************************

# Create two LGD dataset. With one the model is estimated
# On the second the model is tested/audited
set.seed(12345)
len    <- nrow(lgd)
is_ref <- sample(len, 0.7*len, replace = F)
lgd_is <- lgd[is_ref,]  # in_sample dataset
lgd_os <- lgd[-is_ref,] # out_sample dataset
model1 <- lm(lgd_time ~ LTV + purpose1, data = lgd_is)
os_predict <- predict(model1, newdata = lgd_os)

# Histogram of Actual and Predicted
min(lgd_os$lgd_time)
max(lgd_os$lgd_time)
hist(lgd_os$lgd_time, axes = F, xlim = c(-.36,1.24), ylim = c(0,7), prob = T)
axis(side = 1, at = seq(-.36, 1.24, length = 11))
axis(side = 2, at = seq(0, 7, length = 8))
box()
lines(density(lgd$lgd_time, adjust = 0.5), lty = 1)
lines(density(lgd$lgd_time, adjust = 1), lty = 2)
lines(density(lgd$lgd_time, adjust = 1.5), lty = 3)
lines(density(lgd$lgd_time, adjust = 2), lty = 4)
legend("topright", 
       legend = c("adjust=0.5","adjust=1","adjust=1.5","adjust=2"), 
       lty = 1:4, bty = "n", cex = .7, horiz = F)


hist(os_predict, axes = F, xlim = c(-.2,.8), ylim = c(0,3), prob = T)
axis(side = 1, at = seq(-.2, .8, length = 11))
axis(side = 2, at = seq(0, 3, length = 7))
box()
lines(density(os_predict, adjust = 0.5), lty = 1)
lines(density(os_predict, adjust = 1), lty = 2)
lines(density(os_predict, adjust = 1.5), lty = 3)
lines(density(os_predict, adjust = 2), lty = 4)
legend("topright", 
       legend = c("adjust=0.5","adjust=1","adjust=1.5","adjust=2"), 
       lty = 1:4, bty = "n", cex = .7, horiz = F)

# Boxplots of Actual and Predicted
boxplot(lgd_os$lgd_time)
points(mean(lgd_os$lgd_time))

boxplot(os_predict)
points(mean(os_predict))

# ROC Test is missing 
# Prepare it !!

# Correlation Test
cor.test(lgd_os$lgd_time, os_predict, method = "pearson")
cor.test(lgd_os$lgd_time, os_predict, method = "spearman")
cor.test(lgd_os$lgd_time, os_predict, method = "kendall")

# *****************************************************
# BACKTESTING AT LEVEL 2
# *****************************************************
model2 <- lm(lgd_time ~ os_predict, data = lgd_os)
summary(model2)

par(mfrow = c(2,3))
plot(model2, which = 1:6,auto.layout =  F, ask = F)
par(mfrow = c(1,1))

plot(model2$fitted.values, lgd_os$lgd_time, main = "Observed vs Fitted")

# Histgram of the Residuals
hist(model2$residuals, prob = T);box()
res_mean  <- mean(model2$residuals)
res_sd    <- sd(model2$residuals)
res_xaxis <- seq(min(model2$residuals), max(model2$residuals), 
                 length.out = length(model2$residuals))
res_ndensity <- dnorm(res_xaxis, mean = res_mean, sd = res_sd)
lines(res_xaxis, res_ndensity)

# Distribution of LGD time, actual vs predicted
# There is a HUGE difference in the LGD time distribution, between
# Actual LGD vs Predicted LGD
# Actual
lgd_os$lgd_rank <- 
  cut(lgd_os$lgd_time,c(0,quantile(lgd_os$lgd_time, c(.2, .4, .6, .8,1))))
boxplot(lgd_os$lgd_time ~ lgd_os$lgd_rank, main = "Distribution of lgd_time by lgd_rank")
lgd_os_gmean <- as.double(lapply(split( lgd_os$lgd_time,lgd_os$lgd_rank), mean))
points(lgd_os_gmean)

# Predicted
boxplot(model2$fitted.values ~ lgd_os$lgd_rank, 
        main = "Distribution of Predicted by lgd_rank")
os_predict_gmean <- as.double(lapply(
  split( model2$fitted.values ,lgd_os$lgd_rank), mean))
points(os_predict_gmean)

# and the plot of actual predicted means show large errors  
plot(os_predict_gmean,lgd_os_gmean , xlab = "Predicted Values", ylab = "lgd_time",
     ylim = c(0,1), xlim = c(0,1))