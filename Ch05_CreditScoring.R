remove(list = ls())

install.packages("pROC")
library("pROC")

hmeq <-read.csv("hmeq.csv")
attach(hmeq)

# BUILDING LOGISTIC REGRESSION MODELS
# We first remove all na from the dataset 
hmeq_omit        <- na.omit(hmeq)

# JOB and REASON variables are categorical variables. 
# This can be seen as both variables show "levels:"
unique(hmeq_omit$JOB)
unique(hmeq_omit$REASON)

is.factor(hmeq_omit$JOB)
is.factor(hmeq_omit$REASON)

class(hmeq_omit$REASON)
class(hmeq_omit$JOB)

# They do NOT need to be transformed in factor type, as they are ALREADY
#hmeq_omit$JOB    <- as.factor(hmeq_omit$JOB)
#hmeq_omit$REASON <- as.factor(hmeq_omit$REASON)


hmeq_null <- glm(BAD ~ 1, data = hmeq_omit, family = binomial(link = "logit"))
hmeq_full <- glm(BAD ~ LOAN+MORTDUE+VALUE+
                 REASON+JOB+YOJ+DEROG+DELINQ+
                 CLAGE+NINQ+CLNO+DEBTINC, 
                 data = hmeq_omit, family = binomial(link = "logit"))
summary(hmeq_null)
summary(hmeq_full)

# step function is used to select variables
# the number k is the multiple of the number of degrees of freedom 
# used for the penalty and it is tested at 5,4,3 and 2
# the default set up if with k = 5
step(hmeq_null, scope = list(upper = hmeq_full), direction = "both",
     test = "Chisq", data = hmeq_omit, k = 5, trace = F)

step(hmeq_null, scope = list(upper = hmeq_full), direction = "both",
     test = "Chisq", data = hmeq_omit, k = 4, trace = F)

step(hmeq_null, scope = list(upper = hmeq_full), direction = "both",
     test = "Chisq", data = hmeq_omit, k = 3, trace = F)

step(hmeq_null, scope = list(upper = hmeq_full), direction = "both",
     test = "Chisq", data = hmeq_omit, k = 2, trace = F)

hmeq_final <- glm(BAD ~ DEBTINC+DELINQ+DEROG+CLAGE+NINQ+CLNO+JOB,
                  data = hmeq_omit, family = binomial(link = "logit"))
summary(hmeq_final)  
auc(hmeq_omit$BAD, hmeq_final$fitted.values)
round(exp(cbind(coef(hmeq_final), confint(hmeq_final))),4)

# step without indicating the k value, sets by default k = 5
step(hmeq_null, scope = list(upper = hmeq_full), direction = "both",
     test = "Chisq", data = hmeq_omit, trace = F)

hmeq_final2 <- glm(BAD ~ DEBTINC + DELINQ + DEROG + CLAGE + JOB + NINQ + CLNO + LOAN,
                  data = hmeq_omit, family = binomial(link = "logit"))
summary(hmeq_final2)
auc(hmeq_omit$BAD, hmeq_final2$fitted.values)
round(exp(cbind(coef(hmeq_final2), confint(hmeq_final2))),4)
