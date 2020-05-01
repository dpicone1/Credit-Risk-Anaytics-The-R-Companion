# *********************************************************************
# CHAPTER 2 INTRODUCTION TO R
# *********************************************************************

mortgage <-read.csv("mortgage.csv")
attach(mortgage)

# Select the data where the loan has defaulted
mortgage.temp2 <- mortgage[mortgage$default_time == 1,]

len <- length(mortgage$FICO_orig_time)
FICO_cat  <- vector(length = len)
FICO_cat2 <- rep(0, len)

#in case you wanted to remove FICO_cat2 from the memory, use the command remove()  
#remove(FICO_cat2)

for (i in 1:len){
  FICO_cat[i] <- if((FICO_orig_time[i] > 500) & (FICO_orig_time[i] <= 700)) 1
  else if (FICO_orig_time[i] > 700) 2
  else 0
} 

for (i in 1:len){
  if((FICO_orig_time[i] > 500) & (FICO_orig_time[i] <= 700)) {
    FICO_cat2[i] = 1
  }else if (FICO_orig_time[i] > 700) {
      FICO_cat2[i] = 2
  }else{
    FICO_cat2[i] = 0
  } 
} 
f1 <- function(x){ifelse( (x > 500) & (x <= 700), 1, ifelse(x>700, 2,0))}
f2 <- function(x){
  res <- 0
  if( (x > 500) & (x <= 700)){
    res <- 1
  }else if (x > 700){
    res <- 2
  }else{
    res <- 0
  }
  return (res)
}
f2v <- Vectorize(f2)
sum(FICO_cat)
sum(FICO_cat2)
sum(f1(FICO_orig_time))
sum(sapply(FICO_orig_time, f2))
sum(f2v(FICO_orig_time))


mortgage$status_time <- NULL

lista       <- c("default_time", "FICO_orig_time", "LTV_orig_time", "gdp_time")
selec.cols1 <- subset(mortgage, select = c(default_time, FICO_orig_time, LTV_orig_time, gdp_time))
selec.cols2 <- subset(mortgage, select = lista)
selec.cols3 <- mortgage[ ,c("default_time", "FICO_orig_time", "LTV_orig_time", "gdp_time")]
selec.cols4 <- mortgage[ ,lista]

# this does not work!!
lista      <- c(mortgage$default_time, mortgage$FICO_orig_time)
selec.cols <- mortgage[ ,lista]


n.mortgage     <- apply(selec.cols2, 2, length)
mean.mortgage  <- apply(selec.cols2, 2, mean)
st.dev.mortgage<- apply(selec.cols2, 2, sd)
min.mortgage   <- apply(selec.cols2, 2, min)
max.mortgage   <- apply(selec.cols2, 2, max)

select.summary <- cbind(n.mortgage,
                        mean.mortgage,
                        st.dev.mortgage,
                        min.mortgage,
                        max.mortgage)
colnames(select.summary) <- c("N", "Mean", "SD", "Min", "Max")
print(round(select.summary,3))

mortgage.lm <- lm(default_time ~ FICO_orig_time+ LTV_orig_time+ gdp_time,data = mortgage)
summary(mortgage.lm)

example <- function(lhs, rhs){
  d          <- as.data.frame(cbind(lhs, rhs))
  example.lm <- lm(lhs ~ ., data = d)
  return (summary(example.lm))
}

example(default_time, FICO_orig_time)
example(default_time, with(mortgage, cbind(FICO_orig_time , LTV_orig_time, gdp_time)))

# This creates a dataframe as
# either with and cbind
RXS1 <- with(mortgage, cbind(FICO_orig_time , LTV_orig_time, gdp_time))
# or with c() and as.data.frame()
lista <- c("FICO_orig_time", "LTV_orig_time", "gdp_time")
RXS2  <- as.data.frame(mortgage[, lista])

example(default_time, RXS1)
example(default_time, RXS2)