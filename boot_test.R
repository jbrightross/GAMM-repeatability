# Shinichi's attmpt to get 95% boostrapping CI
library(lme4)
library(gamm4)
library(dplyr)
library(rptR)
library(purrr)

rm(list = ls())
#setwd("C:/Users/jbrig/Desktop/Badgers/Badger condition/Processed datasets")
load("Condition modelling output.RData")

# getting data aned formula
formula <- summerobj$mainFormula
data <- summerobj$data

# running the model
mod <- gamm4(formula = formula, 
             data = data, 
             random = ~ (1|Sett) + (1|Individual))

# looking at results
#mod.rr <- ranef(mod$mer,condVar=TRUE)
summary(mod$mer)
summary(mod$gam)

# fixed effect bit of the model
fixed <- predict(mod$mer, re.form = NA) #re.form = NA; not to include any random effects
v.fixed <- var(fixed) # varinace explained by the fixed effect and looks about right = 0.1713035
v.sett <- unname(sqrt(unlist(VarCorr(mod$mer)))["Sett"])
v.individual <- unname(sqrt(unlist(VarCorr(mod$mer)))["Individual"])
v.residual <- getME(mod$mer, "sigma")^2
icc1 <- v.individual/(v.residual + v.sett +  v.individual) # ajudsted repeatablity 1
icc2 <- (v.sett + v.individual)/(v.residual + v.sett +  v.individual) # adjusted repeatablity 2
icc3 <-  v.individual/(v.residual + v.sett +  v.individual + v.fixed)# enhanced repeatablity 1
icc4 <- (v.sett + v.individual)/(v.residual + v.sett +  v.individual + v.fixed) # enhanced repeatablity 2

## helper function to get ICC and variance components

get_v_icc <- function(model) { 
  #beta <- getME(model, "beta")
  sigma <- getME(model, "sigma")
  v.fixed <- var(predict(mod$mer, re.form = NA)) # v for fixed effects
  v.residual <- sigma^2
  v.individual <-  unname(sqrt(unlist(VarCorr(model)))["Individual"])
  v.sett <- unname(sqrt(unlist(VarCorr(model)))["Sett"])
  icc1 <- v.individual/(v.residual + v.sett +  v.individual) # ajudsted repeatablity 1
  icc2 <- (v.sett + v.individual)/(v.residual + v.sett +  v.individual) # adjusted repeatablity 2
  icc3 <-  v.individual/(v.residual + v.sett +  v.individual + v.fixed)# enhanced repeatablity 1
  icc4 <- (v.sett + v.individual)/(v.residual + v.sett +  v.individual + v.fixed) # enhanced repeatablity 2
  est<-c(v.residual,  v.sett, v.individual, v.fixed, icc1, icc2, icc3, icc4) 
  print(est)
}

# simulations

Nboot <- 1000
outcomes <- c("v.residual", "v.sett", "v.individual", "v.fixed", "icc1", "icc2", "icc3", "icc4")
Store <- data.frame(matrix(NA, nrow = Nboot, ncol = length(outcomes)))
colnames(Store) <- outcomes

# warning - will take forever (maybe you want to make Nboot ~10 or 100 to test)
for(i in 1:Nboot){
set.seed(777) # to get the same results every time
# creating new BCI based on the model - this is why it is called "parametric" boostrapping
BCI <- fixed + 
  model.matrix(~factor(data$Individual)-1)%*%rnorm(length(levels(factor(data$Individual))), 0, sqrt(v.individual)) + # Individual effects
  model.matrix(~factor(data$Sett)-1)%*%rnorm(length(levels(factor(data$Sett))), 0, sqrt(v.sett)) + # Sett effeccts
  rnorm(nrow(data), 0, sqrt(v.residual)) # residuals

#cor(fixed, data$BCI)  
data$BCI <- BCI

mod <- gamm4(formula = formula, 
             data = data, 
             random = ~ (1|Sett) + (1|Individual))

Store[i,]<- get_v_icc(mod$mer)
}

# 95% CI for all the varaibles
# it seems pretty tight CI but I guess you have a lot of N so probably correct?
CI <- map(Store, ~quantile(., c(0.025, 0.5, 0.975)))
CI

# to test to see wehther this is really working - you can reduce the N and see whether CI increases - use the half of the data???

# plotting - run out of time - please do!

