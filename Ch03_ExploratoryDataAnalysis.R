# *********************************************************************
# CHAPTER 3 EXPLORATORY DATA ANALYSIS
# *********************************************************************

remove(list=ls())

install.packages("moments")
install.packages("gmodels")
install.packages("vcd")
install.packages("mixtools")
library("moments")
library("gmodels")
library("vcd")
library("mixtools")

# for skewness and kurtosis
install.packages("e1071")
library(e1071)

mortgage <-read.csv("mortgage.csv")
attach(mortgage)

# *********************************************************************
# ONE DIMENSIONAL ANALYSIS
# *********************************************************************

# Observed Frequencies and Empirical Distributions
Frequency <- numeric()
Percent   <- numeric()
Cum_Frequency <- numeric()
Cum_Percent   <- numeric()

defaut_indicator <- unique(default_time)
for(i in 1:2){
  temp <- mortgage[mortgage$default_time == defaut_indicator[i],]
  Frequency[i] <- length(temp$default_time)
  Percent[i]   <- round(Frequency[i]/nrow(mortgage),4)  * 100               
  if (i == 1){
    Cum_Frequency[i] <- Frequency[i]
    Cum_Percent[i]   <- Percent[i]
  }else{
    Cum_Frequency[i] <- Frequency[i-1] + Frequency[i]
    Cum_Percent[i]   <- Percent[i-1] + Percent[i]
  }
}

results <- cbind.data.frame(defaut_indicator,Frequency,Percent,Cum_Frequency,Cum_Percent)
colnames(results) <- c("defaut_time","Frequency","Percent","Cum_Frequency","Cum_Percent")
print(results)

hist(FICO_orig_time, freq = FALSE, breaks = 100, 
     main = "Distr. of Fico at Orig", xlab = "Fico at Orig")

plot.ecdf(FICO_orig_time, main = "Cum. Distr. of Fico at Orig", 
          xlab = "Fico at Orig", ylab = "Cum Prob",pch = ".")

hist(LTV_orig_time, freq = FALSE, breaks = 100, 
     main = "Distr. of LTV at Orig", xlab = "LTV at Orig")
max(LTV_orig_time)
hist(LTV_orig_time, freq = FALSE, breaks = 100, xlim = c(50,110),
     main = "Distr. of LTV at Orig", xlab = "LTV at Orig")

plot.ecdf(LTV_orig_time, main = "Cum. Distr. of LTV at Orig", 
          xlab = "LTV at Orig", ylab = "Cum Prob", verticals= TRUE, pch = ".")

# Location Measures
# Mode Function as R does not have one 
get_mode <- function(x){
  unique_x <- unique(x)
  unique_x[which.max(tabulate(match(x, unique_x)))]}

proc_means<- function(x){
  N      <- length(x)
  Mean   <- mean(x,   na.rm = TRUE)
  Median <- median(x, na.rm = TRUE)
  Mode   <- get_mode(x)
  Pct001 <- quantile(x, 0.01)
  Pct099 <- quantile(x, 0.99)
  proc_means.results <- as.vector(round(cbind(N,Mean,Median,Mode,Pct001,Pct099),4))
  }


var_names    <- c("default_time", "FICO_orig_time", "LTV_orig_time")
loc_measures <- as.data.frame(matrix(NA, nrow = 6, ncol = 3))
for (i in 1:3){
  loc_measures[,i] <-(lapply(mortgage[var_names[i]], proc_means))
}
print(loc_measures)

loc_measures           <- as.data.frame(t(loc_measures),row.names = var_names)
print(loc_measures)
colnames(loc_measures) <- c("N","Mean","Median","Mode","Pct001","Pct099")
print(loc_measures)
loc_measures[,c(2,3,4,5,6)] <- round(loc_measures[,c(2,3,4,5,6)],2)
print(loc_measures)

# Generate a Q-Q plot
qqnorm(mortgage$FICO_orig_time, 
       xlim = c(-6,6), ylim = c(200,1200), 
       main = "QQ plot for FICO at Orig", 
       xlab = "Normal Quantile", ylab = "FICO at Orig")
qqline(mortgage$FICO_orig_time) 

qqnorm(mortgage$LTV_orig_time, 
       xlim = c(-6,6), ylim = c(0,250), 
       main = "QQ plot for LTV at Orig", 
       xlab = "Normal Quantile", ylab = "LTV at Orig")
qqline(mortgage$LTV_orig_time) 


# Dispersion measures
proc_means_ext<- function(x){
  N         <- length(x)
  Minimum   <- min(x,   na.rm = TRUE)
  Maximum   <- max(x, na.rm = TRUE)
  Range     <- range(x)[2]-range(x)[1]
  QuantileRange <- quantile(x, 0.75)-quantile(x, 0.25)
  Var       <- var(x, na.rm = TRUE)
  SD        <- sqrt(Var)
  CoeffVar  <- SD/mean(x,   na.rm = TRUE) *100
  return(as.vector(round(cbind( N,Minimum,Maximum,Range,QuantileRange,Var,SD,CoeffVar),4)))
}


var_names     <- c("default_time", "FICO_orig_time", "LTV_orig_time")
disp_measures <- as.data.frame(matrix(NA, nrow = 8, ncol = 3))
for (i in 1:3){
  disp_measures[,i] <-(lapply(mortgage[var_names[i]], proc_means_ext))
}

disp_measures <- as.data.frame(t(disp_measures),row.names = var_names)
colnames(disp_measures) <- 
  c("N","Min","Max","Range","QuantRange","Var","SD","CoeffVar")
print(round(disp_measures,2))

# Skewness and Kurtosis
proc_means_SkewKurt<- function(x){
  N          <- length(x)
  Skewness   <- skewness(x,   na.rm = TRUE)
  Kurtosis   <- kurtosis(x, na.rm = TRUE) - 3
  return(as.vector(round(cbind( N,Skewness, Kurtosis),4)))
}

var_names <- c("default_time", "FICO_orig_time", "LTV_orig_time")
SkewKurt_measures <- as.data.frame(matrix(NA, nrow = 3, ncol = 3))
for (i in 1:3){
  SkewKurt_measures[,i] <-(lapply(mortgage[var_names[i]], proc_means_SkewKurt))
}

SkewKurt_measures           <- as.data.frame(t(SkewKurt_measures),row.names = var_names)
print(SkewKurt_measures)
colnames(SkewKurt_measures) <- c("N","Skew", "Kurtosis")
print(SkewKurt_measures)

# *********************************************************************
# TWO DIMENSIONAL DATA ANALYSIS
# *********************************************************************

# Joint Empirical Distributions
# We create a two-dimensional ferquency table
FICO_orig_time_factor   <- mortgage$FICO_orig_time
FICO_orig_time_factor[
  FICO_orig_time_factor < quantile(mortgage$FICO_orig_time,.2)
  ] <- 0
FICO_orig_time_factor[
  (FICO_orig_time_factor>= quantile(mortgage$FICO_orig_time,.2))
& 
  (FICO_orig_time_factor < quantile(mortgage$FICO_orig_time,.4))
  ] <- 1
FICO_orig_time_factor[
  (FICO_orig_time_factor >= quantile(mortgage$FICO_orig_time,.4))
  & 
    (FICO_orig_time_factor < quantile(mortgage$FICO_orig_time,.6))
  ] <- 2
FICO_orig_time_factor[
  (FICO_orig_time_factor>= quantile(mortgage$FICO_orig_time,.6))
  & 
    (FICO_orig_time_factor < quantile(mortgage$FICO_orig_time,.8))
  ] <- 3
FICO_orig_time_factor[
  FICO_orig_time_factor  >= quantile(mortgage$FICO_orig_time,.8)
  ] <- 4

# crosstable object is in the library "gmodels"
CrossTable(mortgage$default_time, FICO_orig_time_factor, 
           prop.t = TRUE,prop.r = TRUE, prop.c = TRUE)


boxplot(FICO_orig_time ~ default_time, 
        data = mortgage, 
        range = 0, 
        xlab = "defaut time", ylab = "FICO at Orig Time", 
        main = "Distribution of Fico_orig_time and by default_time")
means <- tapply(mortgage$FICO_orig_time, mortgage$default_time, mean)
points(means, pch = 18)


boxplot(LTV_orig_time ~ default_time, 
        data = mortgage, 
        range = 0, 
        xlab = "defaut time", ylab = "LTV at Orig Time", 
        main = "Distribution of LTV_orig_time and by default_time")
means <- tapply(mortgage$LTV_orig_time, mortgage$default_time, mean)
points(means, pch = 18)

# Correlation Measures
# Chi-sqaure, phi, contingency coeff's and Creamer's V
tab <- xtabs(~ FICO_orig_time_factor + mortgage$default_time)

# assocstats object is in the library("vcd")
assocstats(tab)

# Create a stratified sample equal to 1% of the initial dataset
set.seed(12345)
sample_Fico <- sample(mortgage$FICO_orig_time, size = 0.01*nrow(mortgage), replace = FALSE)
sample_LTV  <- sample(mortgage$LTV_orig_time,  size = 0.01*nrow(mortgage), replace = F)
length(sample_Fico)

# and run various correlation tests.
#cor.test(mortgage$FICO_orig_time,mortgage$LTV_orig_time , method = "pearson")
cor.test(sample_Fico, sample_LTV, method = "pearson")
cor.test(sample_Fico, sample_LTV, method = "spearman", exact = F)
cor.test(sample_Fico, sample_LTV, method = "kendall")

# plot the sample data and overlay a ellipse chart to check the relationship
# ellipse object is in the library("mixtools")
smpl_data <- cbind(sample_Fico, sample_LTV)
plot(smpl_data, xlab = "Fico", ylab = "LTV", main = "Scatter PLot")
ellipse(mu = colMeans(smpl_data), sigma = cov(smpl_data), alpha = 0.1,
        npoints = 250, lwd = 2, col = "blue")
ellipse(mu = colMeans(smpl_data), sigma = cov(smpl_data), alpha = 0.2,
        npoints = 250, lwd = 2, col = "green")
ellipse(mu = colMeans(smpl_data), sigma = cov(smpl_data), alpha = 0.3,
        npoints = 250, lwd = 2, col = "red")
legend(x = "topleft", y.intersp = 0.5, 
       cex = 0.9, title ="prediction elipses",
       legend = c("90%", "80%", "70%"), bty="n", lty = c(1,2,3), lwd=1, horiz = T)

# *********************************************************************
# HIGHLIGHTS OF INDUCTIVE STATISTICS
# *********************************************************************
# Once parameters are estimated, they will not match the true value,
# so now we compute confidence intervalsand Hypothesis testing

# Confidence Intervals
proc_Univariate<- function(x){
  N      <- length(x)
  Mean   <- mean(x, na.rm = TRUE)
  Var    <- var(x,  na.rm = TRUE)
  SD     <- sqrt(Var)
  Lower_Conf<- c(Mean - qnorm(0.99)*SD/sqrt(N))
  Upper_Conf<- c(Mean + qnorm(0.99)*SD/sqrt(N))
  result <-matrix(round(c(Mean,Var,SD,Lower_Conf,Upper_Conf),4),ncol = 5, nrow = 1 )
  colnames(result)  <- c("Mean","Var","SD","Lower_Conf","Upper_Conf")
  print(result)
}

proc_Univariate(mortgage$LTV_orig_time)

# Hypothesis Testing
t.test(mortgage$LTV_orig_time, mu = 60,       alternative = "two.sided")
# We reject the Null Hypothesis H0 : mean = 60, as p-value < 2.2e-16
t.test(mortgage$LTV_orig_time, mu = 78.97546, alternative = "two.sided")
# We cannot reject the Null Hypothesis H0 : mean = 78.97546 as p-value = 1

# *********************************************************************
# Not in the Book, but I personally like to show in a plot
# defaulted vs non defaulted loans. 
# If the pool is too large, I create a subset of it
# *********************************************************************
def  <- mortgage[mortgage$default_time == 1,]
surv <- mortgage[mortgage$default_time == 0,]
nrow(def)
nrow(surv)

plot(def$FICO_orig_time, def$LTV_orig_time, xlab = "Fico", ylab = "LTV", 
     main = "Scatter PLot", col = "red", pch = 18, xlim = c(300, 900))

points(surv$FICO_orig_time, surv$LTV_orig_time, 
       xlab = "Fico", ylab = "LTV", main = "Scatter PLot", col="blue",pch = 1)

legend(x = "topleft", y.intersp = 0.8, 
       cex = 0.8, title ="SubSet (50%) of Mortgage Data",
       legend = c("Defaulted", "No Defauted"), bty="n", col = c("red", "blue"),
       pch = c(18,1),lty)
