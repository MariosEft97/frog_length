rm(list=ls())
#load libraries
library(MASS)
library(moments)
library(robustbase)
library(ellipse)
# Data are imported into RStudio.
frogs = read.csv(file.choose(), header=TRUE, sep = " ")
head(frogs)

# Dataset is divided into training and validation sets.
set.seed(0199514)
d_test = sample(1:dim(frogs)[1], round(dim(frogs)[1]/2))
data_test = frogs[d_test, ]
data_training = frogs[-d_test, ]

attach(data_training)

# Histograms of continuous variables.
par(mfrow = c(2,3))
hist(Length)
hist(Canopy)
hist(Shrub)
hist(Effort)
hist(Size)
# The histogram of Length seems to be right-skewed. A log transformation might be useful.

# Box plots of continuous variables.
par(mfrow = c(2,4))
boxplot(Length, main = "Length Boxplot")
boxplot(Canopy, main = "Canopy Boxplot")
boxplot(Shrub, main = "Shrub Boxplot")
boxplot(Effort, main = "Effort Boxplot")
boxplot(Size, main = "Size Boxplot")
# No outliers are found in the plots of continuous variables (needs to be formally checked).
boxplot(Length~Sex, main = "Sex Boxplot")
boxplot(Length~Natural, main = "Natural Boxplot")
boxplot(Length~Forest, main = "Forest Boxplot")
# There are outliers in Sex and Natural categorical variables.
#means look the same with natural maybe non significant differences

# Scatter plots of continuous variables.
pairs(data_training[2:6], panel = function(x,y) {points(x,y); lines(lowess(x,y), col = "red")})
# Relationship of Length with continuous predictor variables is more or less linear.

# First we fit a model with all the available predictor variables.
fit_all = lm(Length~Sex+Canopy+Shrub+Effort+Size+Natural+Forest, data=data_training)
summary(fit_all)
par(mfrow = c(2,2))
plot(fit_all)

# fit_all_scaled = lm(Length~Sex+scale(Canopy)+scale(Shrub)+scale(Effort)+scale(Size)+Natural+Forest, data=data_training)
# summary(fit_all)
# par(mfrow = c(2,2))
# plot(fit_all)

# Model Selection approaches
# Forward selection based on F-statistic/t-statistic
fit_null = lm(Length ~ 1, data = data_training)
addterm(fit_null, ~ . +Sex+Canopy+Shrub+Effort+Size+Natural+Forest, test = "F")
# add Sex
fit1 = update(fit_null, ~ . + Sex)
addterm(fit1, ~ . +Sex+Canopy+Shrub+Effort+Size+Natural+Forest, test = "F")
# add Shrub
fit2 = update(fit1, ~ . + Shrub)
addterm(fit2, ~ . +Sex+Canopy+Shrub+Effort+Size+Natural+Forest, test = "F")
# add Canopy
fit3 = update(fit2, ~ . + Canopy)
addterm(fit3, ~ . +Sex+Canopy+Shrub+Effort+Size+Natural+Forest, test = "F")
# add Effort
fit4 = update(fit3, ~ . + Effort)
addterm(fit4, ~ . +Sex+Canopy+Shrub+Effort+Size+Natural+Forest, test = "F")
# add Natural
fit5 = update(fit4, ~ . + Natural)
addterm(fit5, ~ . +Sex+Canopy+Shrub+Effort+Size+Natural+Forest, test = "F")
# Length ~ Sex + Shrub + Canopy + Effort + Natural

# Backward elimination based on F-statistic/t-statistic
dropterm(fit_all, test = "F")
# Remove Size
fit1 = update(fit_all, ~ . - Size)
dropterm(fit1, test = "F")
# Remove Forest
fit2 = update(fit1, ~ . - Forest)
dropterm(fit2, test = "F")
# Length ~ Sex + Canopy + Shrub + Effort + Natural

# Backward elimination based on AIC
stepAIC(fit_all, scope = list(upper = ~ Sex+Canopy+Shrub+Effort+Size+Natural+Forest, lower = ~ 1), direction = "backward")
# AIC=-798.56
# Length ~ Sex + Canopy + Shrub + Effort + Natural

# Forward selection based on AIC
fit_null = lm(Length ~ 1, data = data_training)
stepAIC(fit_null, scope = list(upper = ~ Sex+Canopy+Shrub+Effort+Size+Natural+Forest, lower = ~ 1), direction = "forward")
# AIC=-798.56
# Length ~ Sex + Canopy + Shrub + Effort + Natural

# Stepwise selection based on AIC (started at full model)
stepAIC(fit_all, scope = list(upper = ~ Sex+Canopy+Shrub+Effort+Size+Natural+Forest, lower = ~ 1), direction = "both")
# AIC=-798.56
# Length ~ Sex + Canopy + Shrub + Effort + Natural

# Stepwise selection based on AIC (started at null model)
stepAIC(fit_null, scope = list(upper = ~ Sex+Canopy+Shrub+Effort+Size+Natural+Forest, lower = ~ 1), direction = "both")
#  AIC=-798.56
# Length ~ Sex + Shrub + Canopy + Effort + Natural

# Model returned from variable selection
model1 = lm(Length~Sex + Shrub + Canopy + Effort + Natural, data = data_training)
summary(model1)
par(mfrow = c(2,2))
plot(model1)

# Assumptions check
# Standardized residuals
model1_stdres = stdres(model1)

# Studentized residuals
model1_studres = studres(model1)

# Normality of error assumption check
ks.test(model1_stdres, "pnorm")
par(mfrow=c(1,1))
plot(model1, which=2)
# There might be a slight problem with normality.

# Homoscedasticity assumption check
par(mfrow=c(1,1))
plot(residuals(model1), ylab = "Residuals")
# Homoscedasticity is relatively good.

# Independence assumption check for Y
par(mfrow = c(1,1))
plot(model1, which=1)
# There is a slight curved band indicating that the linear model might not be a good fit.
# => Transformation of Y variable might be a solution

# Independence assumption check for X
e = residuals(model1)
par(mfrow = c(1,3))
a = lm(e~Canopy)
plot(Canopy,e,ylab="Residuals")
abline(h=0,lty=3)
plot(Shrub,e,ylab="Residuals")
abline(h=0,lty=3)
plot(Effort,e,ylab="Residuals")
abline(h=0,lty=3)
# plot(factor(Sex),e,ylab="Residuals")
# plot(factor(Natural),e,ylab="Residuals")

#to further test the curvature problem we will use partial residuals
#it is a good indication for adding higher order terms
#partial residuals:
e = residuals(model1)
coefficient = coefficients(model1)
partial_e_canopy = e + coefficient[4] * Canopy
partial_e_shrub = e + coefficient[3]*Shrub
partial_e_effort = e + coefficient[5]* Effort
par(mfrow=c(3,1))
plot(Canopy, partial_e_canopy, ylab = "Partial residual (Canopy)")
abline(lm(unname(partial_e_canopy) ~ Canopy),col='blue')
lines(lowess(Canopy, partial_e_canopy), col = "red") ##linear relationship looks good
plot(Shrub,partial_e_shrub, ylab = "Partial residual (Shrub)")
abline(lm(partial_e_shrub ~ Shrub),col='blue')
lines(lowess(Shrub, partial_e_shrub), col = "red")##linear relationship looks good
plot(Effort, partial_e_effort, ylab = "Partial residual (Effort)")
abline(lm(unname(partial_e_effort) ~ Effort),col='blue')
lines(lowess(Effort, partial_e_effort), col = "red")# Xij are discrete but again the red and black lines are almost identical 
##there is no indication of curvature. 

# To improve curvature of residuals vs fitted values we tried the following:
# 1) adding higher order terms to the exploratory variables.
# 2) Transformation of X variables
# 3) Transformation of Y variables

# Higher order terms
fit_effort_squared  = lm(Length~Sex+Canopy+Shrub+Effort+I(Effort^2)+Natural)
summary(fit_effort_squared)
par(mfrow = c(2,2))
plot(fit_effort_squared)

# Box-Cox transformation
par(mfrow = c(1,1))
data_boxcox = boxcox(model1, lambda = seq(0,4,by=.1))
lambda_max = max(data_boxcox$y)
data_boxcox$x[data_boxcox$y==lambda_max]
# The results of the Box-Cox tranformation indicate that the optimal value of
# lambda is 0 which in turn indicates a log transformation of the dependent variable.

model2 = lm(log(Length)~Sex+Canopy+Shrub+Effort+Natural, data=data_training)
summary(model2)
par(mfrow = c(2,2))
plot(model2)

# We also checked Sex and Natural for a possible interaction
# Interaction plot of Sex-Natural
par(mfrow = c(1,1))
with(data_training, interaction.plot(Sex, Natural, log(Length), fun = mean, type="b", legend = TRUE, main = "Interaction Plot of Sex and Natural", pch=c(1,19)))
# The lines are not completely parallel (interaction indication)

# Interaction Effects
fit_interaction = lm(log(Length)~Sex*Natural+Canopy+Shrub+Effort, data=data_training)
summary(fit_interaction)
par(mfrow = c(2,2))
plot(fit_interaction)

fit_interaction = lm(log(Length)~Sex+Natural+Canopy*Shrub+Effort, data=data_training)
summary(fit_interaction)
par(mfrow = c(2,2))
plot(fit_interaction)

fit_interaction = lm(log(Length)~Sex+Natural+Canopy+Shrub*Effort, data=data_training)
summary(fit_interaction)
par(mfrow = c(2,2))
plot(fit_interaction)

fit_interaction = lm(log(Length)~Sex+Natural+Canopy*Effort+Shrub, data=data_training)
summary(fit_interaction)
par(mfrow = c(2,2))
plot(fit_interaction)

# Standardized residuals model2
model2_stdres = stdres(model2)

skewness(model2_stdres) #slightly higher than zero (normal distribution has skewness=0) smallest than before..much more satisfying 

kurtosis(model2_stdres)# close to 3 just like the normal distribution Kurtosis not different than before
ks.test(model2_stdres,'pnorm') #normality looks good not violated

par(mfrow = c(1,2))
hist(Length)
hist(log(Length),main="Histogram of Transformed Length")

# Transformation of X variables
model3 = lm(log(Length) ~ Sex + log(Shrub/(1-Shrub)) + log(Canopy/(1-Canopy)) + Effort + Natural, data=data_training)
summary(model3)
par(mfrow = c(2,2))
plot(model3)

# Multicollinearity diagnostics
# Correlation matrix
cor_matrix = cor(data_training[2:5])
cor_matrix
# Correlation among the predictor variables are small.

# VIF
vif = diag(solve(cor_matrix))
vif
# VIF values are more or less equal to 1 => No multicollinearity indication.

# Eigenvalues
eigenvalues = eigen(cor_matrix)$values
eigenvalues
sqrt(max(eigenvalues)/eigenvalues)
# All condition numbers are less than 30 => No multicollinearity indication.

# Then we fitted all X variables one at a time with Length
fit_canopy = lm(log(Length)~Canopy, data = data_training)
summary(fit_canopy)
fit_shrub = lm(log(Length)~Shrub, data = data_training)
summary(fit_shrub)
fit_effort = lm(log(Length)~Effort, data = data_training)
summary(fit_effort)
fit_size = lm(log(Length)~Size, data = data_training)
summary(fit_size)
fit_sex = lm(log(Length)~Sex, data = data_training)
summary(fit_sex)
fit_natural = lm(log(Length)~Natural, data = data_training)
summary(fit_natural)
fit_forest = lm(log(Length)~Forest, data = data_training)
summary(fit_forest)
# Although no sign of multicollinearity was found in the tests above,
# fitting each X variables showed that Effort and Natural had
# non-significant coefficients

model4 = lm(log(Length)~Sex+Canopy+Shrub)
summary(model4)
par(mfrow = c(2,2))
plot(model4)

model5 = lm(log(Length)~Sex+log(Shrub/(1-Shrub)) + log(Canopy/(1-Canopy)))
summary(model5)
par(mfrow = c(2,2))
plot(model5)

model6 = lm(log(Length)~Sex+Canopy*Shrub)
summary(model6)
par(mfrow = c(2,2))
plot(model6)
# We have antagonistic interaction between Canopy and Shrub

ks.test(stdres(model6), "pnorm")
par(mfrow=c(1,1))
plot(model6, which=2)

# Outliers-Influencial observations model4
# Standardized residuals
model4_stdres = stdres(model4)
# Studentized residuals
model4_studres = studres(model4)

par(mfrow = c(1,1))
plot(model4_stdres, ylim = c(-4,4), ylab = "Standardized residuals")
abline(h = c(-2.5,2.5), col = "red")

plot(model4_studres, ylim = c(-4,4), ylab = "Studentized residuals")
abline(h = c(-2.5,2.5), col = "red")

par(mfrow = c(2,2))
# Diagonal elements of hat matrix
model4_influence = influence(model4)
plot(model4_influence$hat, ylab = "Diagonal elements of hat matrix")
n = dim(data_training[1:3])[1]
p = dim(data_training[1:3])[2]+1
abline(h = 2*p/n, col = "red")

# DFFITS
model4_dffits = dffits(model4)
plot(model4_dffits, ylab = "DFFITS")
abline(h = 2*sqrt(p/n), col = "red")
# points above the lines are influential observations

# Cook's distance
model4_cd = cooks.distance(model4)
plot(model4_cd, ylab = "Cook's distance")
abline(h = 1, col = "red")

# DFBETAS
model4_dfbetas = dfbetas(model4) 
model4_dfbetas
2/sqrt(n)
plot(model4_dfbetas, ylab = "DFBETAS")
abline(h = 2/sqrt(n), col = "red")

# # RLTS (20% breakdown value)
# model4_robust = ltsReg(log(Length)~Sex+Canopy+Shrub, data=data_training, alpha=0.8)
# summary(model4_robust)
# # Detection of outliers
# plot(model4_robust, which = "rindex")
# plot(model4_robust, which = "rdiag")
# 
# # RLTS (30% breakdown value)
# model4_robust = ltsReg(log(Length)~Sex+Canopy+Shrub, data=data_training, alpha=0.7)
# summary(model4_robust)
# # Detection of outliers
# plot(model4_robust, which = "rindex")
# plot(model4_robust, which = "rdiag")
# 
# RLTS (40% breakdown value)
model4_robust = ltsReg(log(Length)~Sex+Canopy+Shrub, data=data_training, alpha=0.6)
summary(model4_robust)
# Detection of outliers
plot(model4_robust, which = "rindex")
plot(model4_robust, which = "rdiag")

# Robust Standardized Residuals
model4_robust_stdres = model4_robust$residuals/model4_robust$scale
plot(model4_robust_stdres, ylim = c(-4,4), ylab = "Standardized residuals")
abline(h = c(-2.5,2.5), col = "red")

# Diagnostic plot
n = dim(data_training[1:3])[1]
p = dim(data_training[1:3])[2]+1

plot(model4_robust$RD, model4_robust_stdres, ylim = c(-4,4), xlab = "Robust distances", ylab = "Standardized residuals")
abline(v = sqrt(qchisq(0.975, p - 1)), col = "red")
abline(h = c(-2.5,2.5), col = "red")

# Outliers-Influencial observations model6
# Standardized residuals
model6_stdres = stdres(model6)
# Studentized residuals
model6_studres = studres(model6)

par(mfrow = c(1,1))
plot(model6_stdres, ylim = c(-4,4), ylab = "Standardized residuals")
abline(h = c(-2.5,2.5), col = "red")

plot(model6_studres, ylim = c(-4,4), ylab = "Studentized residuals")
abline(h = c(-2.5,2.5), col = "red")

par(mfrow = c(2,2))
# Diagonal elements of hat matrix
model6_influence = influence(model6)
plot(model6_influence$hat, ylab = "Diagonal elements of hat matrix")
n = dim(data_training[1:3])[1]
p = dim(data_training[1:3])[2] + 2
abline(h = 2*p/n, col = "red")

# DFFITS
model6_dffits = dffits(model6)
plot(model6_dffits, ylab = "DFFITS")
abline(h = 2*sqrt(p/n), col = "red")
# points above the lines are influential observations

# Cook's distance
model6_cd = cooks.distance(model6)
plot(model6_cd, ylab = "Cook's distance")
abline(h = 1, col = "red")

# DFBETAS
model6_dfbetas = dfbetas(model6) 
model6_dfbetas
2/sqrt(n)
plot(model6_dfbetas, ylab = "DFBETAS")
abline(h = 2/sqrt(n), col = "red")

#in two of the graphs we observe some influential points 


# RLTS (40% breakdown value)
model6_robust = ltsReg(log(Length)~Sex+Canopy*Shrub, data=data_training, alpha=0.6)
summary(model6_robust)
# Detection of outliers
plot(model6_robust, which = "rindex")
plot(model6_robust, which = "rdiag")

# Robust Standardized Residuals
model6_robust_stdres = model6_robust$residuals/model6_robust$scale
plot(model6_robust_stdres, ylim = c(-4,4), ylab = "Standardized residuals")
abline(h = c(-2.5,2.5), col = "red")

# Diagnostic plot
n = dim(data_training[1:3])[1]
p = dim(data_training[1:3])[2]+2

plot(model6_robust$RD, model6_robust_stdres, ylim = c(-4,4), xlab = "Robust distances", ylab = "Standardized residuals")
abline(v = sqrt(qchisq(0.975, p - 1)), col = "red")
abline(h = c(-2.5,2.5), col = "red")

detach(data_training)

# After thoughful investigation of the models constructed we proceeded with the following 3 models for validation:
# 1) model4
# 2) model5
# 3) model6
# In addition we fitted robust models of the three models above
attach(data_test)
n = dim(data_test)[1]

# Compare estimated coefficients and standard errors
model5_validation = lm(log(Length)~Sex+log(Shrub/(1-Shrub)) + log(Canopy/(1-Canopy)),data = data_test)
summary(model5)
sum5 = summary(model5_validation)
sum5
# The coefficient of Natural is not significant as we suspected and thus we reject model 2.

model4_validation = lm(log(Length)~Sex+Canopy+Shrub,data = data_test)
summary(model4)
sum4 = summary(model4_validation)
sum4

model6_validation = lm(log(Length)~Sex+Canopy*Shrub,data = data_test)
summary(model6)
sum6 = summary(model6_validation)
sum6

model6_robust_validation = ltsReg(log(Length)~Sex+Canopy*Shrub, data=data_test, alpha=0.6)
summary(model6_robust)
sum6_robust = summary(model6_validation)
sum6_robust

# PRESS
PRESS1 = sum((residuals(model4_validation) / (1 - lm.influence(model4_validation)$hat))^2)
PRESS2 = sum((residuals(model5_validation) / (1 - lm.influence(model5_validation)$hat))^2)
PRESS3 = sum((residuals(model6_validation) / (1 - lm.influence(model6_validation)$hat))^2)
PRESS = c(PRESS1, PRESS2, PRESS3)
names(PRESS) = c("model4", "model5", "model6")
sort(round(PRESS,3))

# MSE
MSE1 = summary(model4_validation)$sigma^2
MSE2 = summary(model5_validation)$sigma^2
MSE3 = summary(model6_validation)$sigma^2
MSE = c(MSE1, MSE2, MSE3)
names(MSE) <- c("model4", "model5", "model6")
sort(round(MSE,4))

# MSEP
MSEP1 = mean((predict(model4_validation, newdata = data_test) - Length)^2)
MSEP2 = mean((predict(model5_validation, newdata = data_test) - Length)^2)
MSEP3 = mean((predict(model6_validation, newdata = data_test) - Length)^2)
MSEP = c(MSEP1, MSEP2, MSEP3)
names(MSEP) <- c("model4", "model5", "model6")
sort(round(MSEP,4))

# Results
validation_results = data.frame(rbind(PRESS/n, MSE, MSEP), row.names = c("PRESS/n", "MSE", "MSEP"))
names(validation_results) = c("model4", "model5", "model6")
validation_results
# Model 2 is the optimal model with respect to PRESS/n, MSE and MSEP

# validation R^2
rsum4 = sum4$r.squared
rsum5 = sum5$r.squared
rsum6 = sum6$r.squared
R_val = c(rsum4,rsum5,rsum6)
names(R_val) = c("model4_val R^2", "model5_val R^2", "model6_val R^2")
round(R_val,4)
#training R^2
R_tr=c(summary(model4)$r.squared,summary(model5)$r.squared,summary(model6)$r.squared)
names(R_tr)=c("model4_tr R^2", "model5_tr R^2", "model6_tr R^2")
round(R_tr,4)
#very close


# AR^2
arsum4 = sum4$adj.r.squared
arsum5 = sum5$adj.r.squared
arsum6 = sum6$adj.r.squared
AR = c(arsum4,arsum5,arsum6)
names(AR) = c("model4 AR^2", "model5 AR^2", "model6 AR^2")
round(AR,4)

# AIC
AIC(model4, model4_validation)
AIC(model5, model5_validation)
AIC(model6, model6_validation)

# Outliers-Influencial observations
# Standardized residuals
model4_validation_stdres = stdres(model4_validation)
model5_validation_stdres = stdres(model5_validation)
model6_validation_stdres = stdres(model6_validation)


par(mfrow = c(1,3))

plot(model4_validation_stdres, ylim = c(-4,4), ylab = "Standardized residuals",main='Model 4 Outliers Detection')
abline(h = c(-2.5,2.5), col = "red")
plot(model5_validation_stdres, ylim = c(-4,4), ylab = "Standardized residuals",main='Model 5 Outliers Detection')
abline(h = c(-2.5,2.5), col = "red")
plot(model6_validation_stdres, ylim = c(-4,4), ylab = "Standardized residuals",main='Model 6 Outlier Detection')
abline(h = c(-2.5,2.5), col = "red")


# Studentized residuals
model4_validation_studres = studres(model4_validation)
model5_validation_studres = studres(model5_validation)
model6_validation_studres = studres(model6_validation)


plot(model4_validation_studres, ylim = c(-4,4), ylab = "Studentized residuals",main='Model 4 Outlier Detection')
abline(h = c(-2.5,2.5), col = "red")
plot(model5_validation_studres, ylim = c(-4,4), ylab = "Studentized residuals",main='Model 5 Outlier Detection')
abline(h = c(-2.5,2.5), col = "red")
plot(model6_validation_studres, ylim = c(-4,4), ylab = "Studentized residuals",main='Model 6 Outlier Detection')
abline(h = c(-2.5,2.5), col = "red")

# Model 4 has potential outliers-influential points.
# Model 5 has potential outliers-influential points.
# Model 6 does not have outliers-influential points.

##coefficients comparison
#model 4
training_4<-model4$coefficients
test_4<-model4_validation$coefficients
cbind(training_4,test_4)

#model 5
training_5<-model5$coefficients
test_5<-model5_validation$coefficients
cbind(training_5,test_5)

#model6
training_6<-model6$coefficients
test_6<-model6_validation$coefficients
cbind(training_6,test_6)
target_values<-log(data_test$Length)
#plot the fitted values againts the actual values
par(mfrow = c(1,1))
plot(model6_validation$fitted.values,target_values,ylab='Target',xlab='Fitted',main='Fitted Values vs Target Values')
abline(a=0,b=1,col='red')


#statistical inference
summary(model6_validation)
#Partial F-tests
mod1<-lm(log(Length)~Canopy+Shrub,data=data_test)
mod2<-lm(log(Length)~Sex+Canopy+Shrub,data=data_test)
mod3<-lm(log(Length)~Canopy:Shrub,data=data_test)
anova(model6_validation,mod2)
anova(model6_validation,mod1)
anova(model6_validation,mod3)
#we confirm that the model can't be simplified more
# Individual confidence intervals

alpha = 0.05
confint(model6_validation, level=1-alpha)

# Simultaneous confidence intervals with Bonferroni correction
confint(model6_validation, level=1-alpha/2)


# Simultaneous confidence intervals with Bonferroni correction

con_int<-confint(model6_validation, parm = c(2,5), level = 1 - alpha / 2) 
# wider still no 0 in the confidence intervals

# Joint confidence region

alpha <- 0.05
plot(ellipse(model6_validation, which = c(2,5), level = 1 - alpha), type = "l",main='Joint Confidence Interval')
points(model6_validation$coefficients[2],model6_validation$coefficients[5])
abline(v = con_int[1,1], lty = 2,col='red')
abline(v = con_int[1,2], lty = 2,col='red')
abline(h = con_int[2,1], lty = 2,col='red')
abline(h = con_int[2,2], lty = 2,col='red')
detach(data_test)
