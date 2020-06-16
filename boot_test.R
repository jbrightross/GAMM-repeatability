# Shinichi's attmpt to get 95% boostrapping CI
library(lme4)
library(gamm4)
library(dplyr)
library(rptR)
library(purrr)

rm(list = ls())
#setwd("C:/Users/jbrig/Desktop/Badgers/Badger condition/Processed datasets")
load("Condition modelling output.RData")

# getting data and formula
formula <- summerobj$mainFormula
data    <- summerobj$data

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
v.fixed <- var(fixed) # variance explained by the fixed effect and looks about right = 0.1713035
v.sett <- unname(sqrt(unlist(VarCorr(mod$mer)))["Sett"]) #variance explained by sett
v.individual <- unname(sqrt(unlist(VarCorr(mod$mer)))["Individual"]) #variance explained by individual
v.residual <- getME(mod$mer, "sigma")^2 #residual variance
icc1 <- v.individual/(v.residual + v.sett +  v.individual) # adjusted repeatability 1
icc2 <- (v.sett + v.individual)/(v.residual + v.sett +  v.individual) # adjusted repeatability 2
icc3 <- v.individual/(v.residual + v.sett + v.individual + v.fixed) # enhanced repeatability 1
icc4 <- (v.sett + v.individual)/(v.residual + v.sett + v.individual + v.fixed) # enhanced repeatablity 2

## helper function to get ICC and variance components

get_v_icc <- function(model) { 
  #beta <- getME(model, "beta")
  sigma <- getME(model, "sigma")
  v.fixed <- var(predict(mod$mer, re.form = NA)) # v for fixed effects
  v.residual <- sigma^2
  v.individual <-  unname(sqrt(unlist(VarCorr(model)))["Individual"])
  v.sett <- unname(sqrt(unlist(VarCorr(model)))["Sett"])
  icc1 <- v.individual/(v.residual + v.sett +  v.individual) # adjusted repeatability 1
  icc2 <- (v.sett + v.individual)/(v.residual + v.sett +  v.individual) # adjusted repeatability 2
  icc3 <- v.individual/(v.residual + v.sett + v.individual + v.fixed) # enhanced repeatability 1
  icc4 <- (v.sett + v.individual)/(v.residual + v.sett + v.individual + v.fixed) # enhanced repeatablity 2
  est <- c(v.residual,  v.sett, v.individual, v.fixed, icc1, icc2, icc3, icc4) 
  
  return(est)
}

# simulations

Nboot    <- 100
outcomes <- c("v.residual", "v.sett", "v.individual", "v.fixed", "icc1", "icc2", "icc3", "icc4")
Store    <- data.frame(matrix(NA, nrow = Nboot, ncol = length(outcomes)))
colnames(Store) <- outcomes

# warning - will take forever (maybe you want to make Nboot ~10 or 100 to test)

set.seed(777) # to get the same results every time

for (i in 1:Nboot) {

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
  
  if (i %% 100 == 0) {
    print(paste(100*i/Nboot, "percent done."))
  }
}

# 95% CI for all the variables
# it seems pretty tight CI but I guess you have a lot of N so probably correct? 
# ^You had put the set.seed() command inside the loop so it was generating identical repeatabilities each time
CI <- map(Store, ~quantile(., c(0.025, 0.5, 0.975)))

# Testing if this is really working: expanding stuff out to include different amounts of data #

Nboot    <- 50
outcomes <- c("v.residual", "v.sett", "v.individual", "v.fixed", "icc1", "icc2", "icc3", "icc4")
Store    <- data.frame(matrix(NA, nrow = Nboot*3, ncol = length(outcomes)))
colnames(Store) <- outcomes

data_fraction <- c(rep(1/2, Nboot), rep(1/4, Nboot), rep(1/8, Nboot))

dataOriginal <- summerobj$data
data         <- dataOriginal

set.seed(777) # to get the same results every time

# I rewrote this to sample the data down and remodel, but I must have screwed something up for the 
# model.matrix command because I'm getting warning messages about the "fixed + model.matrix..." line 
# including different lengths of objects

for (i in 1:(Nboot*3)) {
  
  # refitting model with different amounts of data when relevant
  
  if (i != 1) {
    if (i %% Nboot == 1) { # So, for first iteration of each reduction of data
      data <- dataOriginal[sample(nrow(dataOriginal), nrow(dataOriginal)*data_fraction[i], replace = FALSE),]
      
      mod <- gamm4(formula = formula, 
                   data = data, 
                   random = ~ (1|Sett) + (1|Individual))
      
      fixed <- predict(mod$mer, re.form = NA) #re.form = NA; not to include any random effects 
      v.fixed <- var(fixed) # variance explained by the fixed effect and looks about right = 0.1713035
      v.sett <- unname(sqrt(unlist(VarCorr(mod$mer)))["Sett"]) #variance explained by sett
      v.individual <- unname(sqrt(unlist(VarCorr(mod$mer)))["Individual"]) #variance explained by individual
      v.residual <- getME(mod$mer, "sigma")^2 #residual variance
      icc1 <- v.individual/(v.residual + v.sett +  v.individual) # adjusted repeatability 1
      icc2 <- (v.sett + v.individual)/(v.residual + v.sett +  v.individual) # adjusted repeatability 2
      icc3 <- v.individual/(v.residual + v.sett + v.individual + v.fixed) # enhanced repeatability 1
      icc4 <- (v.sett + v.individual)/(v.residual + v.sett + v.individual + v.fixed) # enhanced repeatablity 2
    }
  }
  
  # creating new BCI based on the model - this is why it is called "parametric" boostrapping
  
  BCI <- fixed + 
    model.matrix(~factor(data$Individual)-1)%*%rnorm(length(levels(factor(data$Individual))), 0, sqrt(v.individual)) + # Individual effects
    model.matrix(~factor(data$Sett)-1)%*%rnorm(length(levels(factor(data$Sett))), 0, sqrt(v.sett)) + # Sett effects
    rnorm(nrow(data), 0, sqrt(v.residual)) # residuals
  
  #cor(fixed, data$BCI)  
  
  data$BCI <- BCI
  
  modTemp <- gamm4(formula = formula, 
                   data = data, 
                   random = ~ (1|Sett) + (1|Individual))
  
  Store[i,] <- get_v_icc(modTemp$mer)
  
  if (i %% 50 == 0) {
    print(paste(100*i/(Nboot*3), "percent done."))
  }
}

Store$datafraction <- data_fraction

# plotting - run out of time - please do!

library(ggplot2)

ggplot(Store, aes(x = datafraction, group = datafraction, y = icc2)) + geom_boxplot()

# Looks alright to me--spread gets wider with more data and the average moves up because we're rarefying 
# the data and ending up with more singletons (more repeatability if more individuals have R = 1)