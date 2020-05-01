remove(list=ls())

install.packages("pastecs")
library("pastecs")

hmeq <-read.csv("hmeq.csv")
attach(hmeq)

# SAMPLING
# we create a stratified sample of 1000 obs without replacement
set.seed(12345)
ref <- sample(nrow(hmeq), replace = F, size = 1000)

# the row of the new dataframe mySample will be randomly chosen to form a new dataframe
mySample <- hmeq[ref,] 
nrow(hmeq)

Freq_Proc_population               <- data.frame(unique(hmeq[,"BAD"]))
colnames(Freq_Proc_population)     <- "BAD"
Freq_Proc_population[,"Frequency"] <- c(nrow(hmeq[hmeq$BAD==0,]),nrow(hmeq[hmeq$BAD==1,]))
Freq_Proc_population[,"Perc"]      <- Freq_Proc_population[,"Frequency"]/nrow(hmeq)
print (Freq_Proc_population)

Freq_Proc_sample               <- data.frame(unique(mySample[,"BAD"]))
colnames(Freq_Proc_sample)     <- "BAD"
Freq_Proc_sample[,"Frequency"] <- c(nrow(mySample[mySample$BAD==0,]),
                                    nrow(mySample[mySample$BAD==1,]))
Freq_Proc_sample[,"Perc"]      <- Freq_Proc_sample[,"Frequency"]/nrow(mySample)
print (Freq_Proc_sample)

# DESCRIPTIVE STATISTICS
round(stat.desc(mySample$LOAN),2)
round(stat.desc(hmeq$LOAN),2)

# MISSING VALUES
# replacing missing values (for eaxch columns) with their respective column means
mySample_cnt  <- mySample[, -c(1,5,6)]
for(i in 1:ncol(mySample_cnt)){
  mySample_cnt[is.na(mySample_cnt[,i]), i] <- mean(mySample_cnt[,i], na.rm = T)
}

# OUTLIER DETECTION AND TREATMENT 
mySample_mean  <- apply(mySample_cnt,2, mean) # per column
mySample_sd    <- apply(mySample_cnt,2, sd)   # per column
zscores        <- matrix(nrow = nrow(mySample_cnt), ncol = ncol(mySample_cnt))
# compute a zscore for each cell of the matrix mySample_cnt
for(i in 1:ncol(mySample_cnt)){
  zscores[, i] <- (mySample_cnt[,i]- mySample_mean[i])/mySample_sd[i]
}
# calculate the max zscore in each row
max_abs <- apply(zscores, 1, function(x) max(abs(x)))

# extract those rows where the max_abs is below 3, and write them into a new matrix
# In summary, for each column, we have calculated the zscore in each cell of the matrix
# and now we remove the full row where at least one cell of the row is greater or equal to
# the zscore. the new matrix is called "filtered_example"  
filtered_example <- zscores[max_abs< 3, ] 
nrow(zscores)
nrow(filtered_example)

# alternatively, you can also use the function subset. The result is the same.
filtered_example2 <- subset(zscores, max_abs< 3)
# Testing that it is the same
nrow(zscores)
nrow(filtered_example)
nrow(filtered_example2)
sum(filtered_example)
sum(filtered_example2)

# CATEGORIZATION
# Example of Contingency Table. 
residence          <- as.data.frame(matrix(nrow = 12, ncol = 3))
colnames(residence)<- c("default", "resstatus", "count")
residence[,1]      <- c(rep("good",6), rep("bad",6))
residence[,2]      <- rep(c("owner", "rentunf", "rentfurn","withpar","other","noanswer"),2)
residence[,3]      <- as.double(c(6000,1600,350,950,90,10,300,400,140,100,50,10))


coarse1           <- as.data.frame(matrix(nrow = 6, ncol = 3))
colnames(coarse1) <- c("default", "resstatus", "count")
coarse1[,1]       <- c(rep("good",3), rep("bad",3))
coarse1[,2]       <- rep(c("owner", "withpar","other"),2)
coarse1[,3]       <- as.double(c(6000,950,1050,300,100,600))

coarse2           <- as.data.frame(matrix(nrow = 6, ncol = 3))
colnames(coarse2) <- c("default", "resstatus", "count")
coarse2[,1]       <- c(rep("good",3), rep("bad",3))
coarse2[,2]       <- rep(c("owner", "withpar","other"),2)
coarse2[,3]       <- as.double(c(6000,950,2050,300,100,600))


coarse1_tbl       <-table(coarse1[,1],coarse1[,2])
coarse1_tbl[1,]   <-coarse1[c(6,4,5),3]
coarse1_tbl[2,]   <-coarse1[c(3,1,2),3]
coarse1_tbl

coarse2_tbl       <-table(coarse2[,1],coarse2[,2])
coarse2_tbl[1,]   <-coarse2[c(6,4,5),3]
coarse2_tbl[2,]   <-coarse2[c(3,1,2),3]
coarse2_tbl

prop.table(coarse1_tbl,1)
prop.table(coarse1_tbl,2)
chisq.test(coarse1_tbl)

prop.table(coarse2_tbl,1)
prop.table(coarse2_tbl,2)
chisq.test(coarse2_tbl)
