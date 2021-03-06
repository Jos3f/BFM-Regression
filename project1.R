library(MASS) # needed for studres
library(car) # needed for avPlots
library(leaps) # needed for regsubsets
library(glmnet)
library(simpleboot)
source("resplot.R")
source("functions_sf2930.R") 


# Read data 
bodyfatmen=read.csv("bodyfatmen.csv")
bodyfatwomen=read.csv("bodyfatwomen.csv")
data_to_use = bodyfatmen

#Simple linear regression without transformation
value.column <- "density"
fit.lm <- lm(bodyfatmen$density ~ ., data=bodyfatmen)
# fit.lm <- lm(bodyfatwomen$DEXfat ~ ., data=bodyfatwomen[,-2])
summary(fit.lm)

# Normal probability plot of the residuals
par(mfrow = c(2, 2))
qqnorm(rstandard(fit.lm))
qqnorm(rstudent(fit.lm))
qqnorm(fit.lm$residuals)
# QQ-plot seems ok. just a little hint of a heavy tailed distribution. 

# Plot of residuals versus the predicted response
plot(fit.lm$fitted.values, fit.lm$residuals)
plot(fit.lm$fitted.values, rstandard(fit.lm), xlab="Fitted values", ylab="standardized residuals")
plot(fit.lm$fitted.values, rstudent(fit.lm), xlab="Fitted values", ylab="studentized residuals")
# Uneven distribution of the values but otherwise it looks good. 


# Print residuals versus every attribute (Not really useful?)
par(mfrow = c(2, 2))
for (column_index in colnames(data_to_use)[-1]){
  print((data_to_use[[column_index]])[1]) 
  plot(data_to_use[[column_index]], rstudent(fit.lm), xlab=column_index,ylab="studentized res")
}
# Seems like residuals depend on density. The Box - Cox Method might be useful


# Same as exercise_4_2.R, might be enough for residual testing
p <- fit.lm$rank - 1
n <- nrow(data_to_use)
leverage.cutoff <- 2*p/n # Montgomery p. 213
cooks.cutoff <- qf(0.5, p, n - p, lower.tail = FALSE) # Montgomery p. 215
studres.cutoff <- qt(0.05/2, n - p, lower.tail = FALSE) # Montgomery p. 135
leverage.cutoff
cooks.cutoff
studres.cutoff

#par(mfrow = c(2, 3))
par(mfrow = c(1, 1))
plot(fit.lm, which = 1:5, pch = 20)
infl.lm <- influence.measures(fit.lm)
summary(infl.lm)
# Some of the residuals are large (200, 203 and 220) but from residuals vs.
# fitted it seems ok. Note that these values are outside the 5% quantile in the
# t(n-p) distrubution, but since we have 248 observations we expect
# about 12 observations to be this "extreme" (possible outliers).
#
# Many data points (5,31,36,39,41,52,83...) have hat value over the leverage cutoff
# of 0.1048 (Montgomery p. 213) which makes them interesting (considered leverage points) to look closer at in case they are influential.
# 
# We use cooks distance to determine if we should remove any points that are influential, which we should 
# consider analysing since we had some leverage points. Doing so we can se that point 39, 83 and 217 stick
# out from the rest with large cooks distances. However, these distances are smaller than the cutoff of 
# cutoff of 0.9519228 (Montgomery p.215)
# 
# 


# Partial regression plots
avPlots(fit.lm)
# No need for transformations?

# Does the correlation matrix give any indication of multicollinearity?
X1 <- model.matrix(fit.lm) # get the X matrix
X <- X1[, -1] # remove the intercept column (column 1)
WtW <- cor(X)
WtW
# Many cases of large co-variance values. Check weight and chest for example. We need to investigate further. 

# Calculate the variance inflation factor and the condition number of
# X'X. Is there any evidence for multicollinearity?

# VIF
vif(fit.lm)
# compare to the diabonal elements of
solve(WtW)
# Condition number
eval <- eigen(WtW)
eval
condition_number <- sqrt(max(eval$values) / min(eval$values))
condition_number # ( see p. 298 Montgomery)
condition_number_no_sqrt <- (max(eval$values) / min(eval$values))
condition_number_no_sqrt # ( see p. 298 Montgomery)


# Variable selection:

predict.regsubsets = function (object ,newdata ,id ,...){
  form=as.formula(object$call [[2]])
  mat=model.matrix(form ,newdata )
  coefi =coef(object ,id=id)
  xvars =names(coefi )
  mat[,xvars ]%*% coefi
}


methods =  list("backward","forward", "seqrep","exhaustive")

for(m in methods){
  
  k=10
  variables = 13
  
  set.seed (1)
  folds=sample (1:k,nrow(data_to_use ),replace =TRUE)
  cv.errors =matrix (NA ,k,variables, dimnames =list(NULL , paste (1:variables) ))
  
  r.sqr = matrix (NA ,k,variables, dimnames =list(NULL , paste (1:variables) ))
  r.sqr.adj = matrix (NA ,k,variables, dimnames =list(NULL , paste (1:variables) ))
  Cp = matrix (NA ,k,variables, dimnames =list(NULL , paste (1:variables) ))
  
  
  for(j in 1:k){
    best.fit = regsubsets (density~.,data=data_to_use[folds !=j,],method = m,
                           nvmax =variables)
    
    for(i in 1:variables) {
      pred=predict (best.fit ,data_to_use[folds ==j,], id=i)
      cv.errors [j,i]=mean( (data_to_use$density[folds ==j]-pred)^2)
      r.sqr [j,i]= (1 - ( sum((data_to_use$density[folds ==j]-pred)^2) / sum((data_to_use$density[folds ==j]-mean(data_to_use$density[folds ==j] ))^2 )) )
      n = nrow(data_to_use[folds !=j,])
      r.sqr.adj [j,i]= 1 - (1-r.sqr[j,i])*(n - 1)/(n-i-1)
      n_test = nrow(data_to_use[folds ==j,])
      Cp[j,i] = (sum((data_to_use$density[folds ==j]-pred)^2)/(summary(fit.lm)$sigma)^2) - n_test + 2*(i-1)
      
      
      
    }
  }
  
  mean.cv.errors =apply(cv.errors ,2, mean)
  par(mfrow =c(2,3))
  plot(mean.cv.errors ,type='l', main=m)
  points(which.min(mean.cv.errors),
         mean.cv.errors[which.min(mean.cv.errors)],
         col = "red", cex = 2, pch = 20)
  
  mean.r.sqr = apply(r.sqr ,2, mean)
  plot(mean.r.sqr ,type='l', main=m, ylab = "R-squared", xlab = "Number of variables")
  points(which.max(mean.r.sqr),
         mean.r.sqr[which.max(mean.r.sqr)],
         col = "red", cex = 2, pch = 20)
  
  mean.r.sqr.adj = apply(r.sqr.adj ,2, mean)
  plot(mean.r.sqr.adj ,type='l', main=m, ylab = "adj-R-squared", xlab = "Number of variables")
  points(which.max(mean.r.sqr.adj),
         mean.r.sqr.adj[which.max(mean.r.sqr.adj)],
         col = "red", cex = 2, pch = 20)
  
  mean.Cp = apply(Cp ,2, mean)
  #plot(mean.Cp ,type='p', main=m, ylab = "Mallows's Cp", xlab = "Number of variables", ylim = c(0,20))
  # lines(x = c(0,100), y = c(0,100),  lty=2 , col = "blue")


  
  reg.best=regsubsets (density~.,data=data_to_use , nvmax =variables, method = m)
  reg.best.summary <- summary(reg.best)
  plot(reg.best.summary$bic,
       xlab = "Number of Variables ", ylab = "BIC", type = "l")
  points(which.min(reg.best.summary$bic),
         reg.best.summary$bic[which.min(reg.best.summary$bic)],
         col = "red", cex = 2, pch = 20)
  
  plot(reg.best.summary$cp,
       xlab = "Number of Variables ", ylab = "Mallows's Cp", type = "l")
  points(which.min(reg.best.summary$cp),
         reg.best.summary$cp[which.min(reg.best.summary$cp)],
         col = "red", cex = 2, pch = 20)
  
  
  # plot(reg.best)
  print(reg.best.summary)
  coef(reg.best ,8)
}




# Total search for finding best subset:
regfit.all <- regsubsets(bodyfatmen$density ~ ., data = data_to_use, nvmax=13)
summary(regfit.all)
plot(regfit.all)
coef(regfit.all, 4) # BIC suggests 4 variables
coef(regfit.all, 8) # Cp suggests 8 variables
coef(regfit.all, 9) # Adjusted RSq suggests 9 variables

# Forward
regfit.forward <- regsubsets(bodyfatmen$density ~ ., data = data_to_use, method = "forward", nvmax=13)
summary(regfit.forward)
plot(regfit.forward)
coef(regfit.all, 4) # BIC suggests 4 variables
coef(regfit.all, 9) # Cp suggests 9 variables
coef(regfit.all, 9) # Adjusted RSq suggests 9 variables

# Backward
regfit.backward <- regsubsets(bodyfatmen$density ~ ., data = data_to_use, method = "backward", nvmax=13)
summary(regfit.backward)
plot(regfit.backward)
coef(regfit.all, 4) # BIC suggests 4 variables
coef(regfit.all, 8) # Cp suggests 8 variables
coef(regfit.all, 9) # Adjusted RSq suggests 9 variables

# Stepwise selection (or sequential replacement)
regfit.step <- regsubsets(bodyfatmen$density ~ ., data = data_to_use, method = "seqrep", nvmax=13)
summary(regfit.step)
plot(regfit.step)
coef(regfit.all, 4) # BIC suggests 4 variables
coef(regfit.all, 8) # Cp suggests 8 variables
coef(regfit.all, 9) # Adjusted RSq suggests 9 variables

# p 122 in Modern multivariate statistical techniques suggests 10 or 5 fold CV to balance bias and variance. 


# Bootstrapping confidence intervalls:
#lboot <- lm.boot(fit.lm, R = 1000)
#summary(lboot)
#hist(lboot)
# Other better method
fit.lm.boot <- Boot(fit.lm, method="residual", R=1000) # Smaller R gives larger confidence intervalls
summary(fit.lm.boot)
confint(fit.lm.boot)
hist(fit.lm.boot)
