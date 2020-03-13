# **************************************************************************
# 1 - UNDERSAMPLING AND OVERSAMPLING
# 2 - ADJUSTING POSTERIOUS PROBABILITIES
# 3 - MAPPING TO AN EXTERNAL RATING AGENCY
# 4 - CONFIDENCE LEVEL BASED APPROACH
# **************************************************************************

remove(list = ls())

mortgage     <- read.csv("mortgage.csv")
data.ratings <- read.csv("ratings.csv")
library("MASS")

# **************************************************************************
# UNDERSAMPLING AND OVERSAMPLING
# **************************************************************************
mortgage0 <- mortgage[mortgage$default_time == 0,]
mortgage1 <- mortgage[mortgage$default_time == 1,]
unique(mortgage0[,"default_time"])
unique(mortgage1[,"default_time"])

set.seed(12345)
ref0 <- sample(nrow(mortgage0), size = 1000, replace = F)
ref1 <- sample(nrow(mortgage1), size = 1000, replace = F)

mortgage_sample0 <- mortgage0[ref0,]
mortgage_sample1 <- mortgage1[ref1,]
unique(mortgage_sample0[,"default_time"])
unique(mortgage_sample1[,"default_time"])

# We now join the two samples
mortgage_sample  <- rbind(mortgage_sample0, mortgage_sample1)

FREQ_Procedure <- data.frame(unique(mortgage_sample[,"default_time"]))
FREQ_Procedure[, "Frequency"]<- c(nrow(mortgage_sample0), nrow(mortgage_sample1))

# **************************************************************************
# ADJUSTING POSTERIOUS PROBABILITIES
# **************************************************************************
probdef       <- c(.1,.3,.5,.6,.85,.9)
probnondef    <- 1- probdef
oldprior      <- 0.20
newprior      <- 0.01
temp1         <- probdef/oldprior * newprior
temp2         <- (1-probdef)/(1-oldprior) * (1-newprior)
newprobdef    <- temp1/(temp1+temp2)
newprobnondef <- temp2/(temp1+temp2)
adjust_posterior <- cbind(probdef, probnondef, newprobdef, newprobnondef)
print(adjust_posterior)

# **************************************************************************
# MAPPING TO AN EXTERNAL RATING AGENCY
# **************************************************************************
data.ratings$rating <- factor(data.ratings$rating, ordered = T)
model_logit         <- polr(data.ratings$rating ~ 
                          data.ratings$COMMEQTA + data.ratings$LLPLOANS +
                          data.ratings$COSTTOINCOME + data.ratings$ROE + 
                          data.ratings$LIQASSTA + data.ratings$SIZE, 
                          method = "logistic")
summary(model_logit)

prediction <- predict(model_logit, data.ratings[,-1])
notchdiff  <- abs(as.double(data.ratings$rating) - as.double(prediction))
notchdiff  <- notchdiff[order(notchdiff)]
FREQ_Procedure2 <- data.frame(unique(notchdiff))
dim(FREQ_Procedure2)
Freq <- rep(0,7)
for (i in 1:7){
  Freq[i] <- length(subset(notchdiff, notchdiff == FREQ_Procedure2[i,1])) 
}

FREQ_Procedure2[,"Frequency"] <- Freq
FREQ_Procedure2[,"Percent"]   <- round(FREQ_Procedure2[,"Frequency"]/length(notchdiff)*100,3)
FREQ_Procedure2[,"Cum Freq"]  <- cumsum(FREQ_Procedure2[,"Frequency"])
FREQ_Procedure2[,"Cum Perc"]  <- cumsum(FREQ_Procedure2[,"Percent"])
colnames(FREQ_Procedure2)[colnames(FREQ_Procedure2) == "unique.notchdiff."] <- "Errors Bands"
print (FREQ_Procedure2)
barplot(height = FREQ_Procedure2[,5], names.arg = FREQ_Procedure2[,1],
        main="Cum Distr of notchdiff", 
        ylab = "Percent", xlab = "notchdiff")

# **************************************************************************
# CONFIDENCE LEVEL BASED APPROACH
# **************************************************************************
# this is an example for confidence levels
nA  <- 100
nB  <- 150
nC  <- 80
sig <- 0.99
PDA <- 1 - (1-sig)^(1/(nA + nB + nC))
PDB <- 1 - (1-sig)^(1/(nB + nC))
PDC <- 1 - (1-sig)^(1/(nC))

print (cbind(PDA, PDB, PDC))