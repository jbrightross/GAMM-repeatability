########################################################
# Title: Bootstrapping repeatability and BLUPs for
# paper
# Name: Julius G. Bright Ross
# Description: Using bootstrapping procedure to 
# get CIs for repeatability estimates and to 
# incorporate individual data spread into the 
# relationship between individual random effects and 
# lifetime life history metrics
########################################################

## Clearing workspace and loading packages/data ##

rm(list = ls())

library(lme4)
library(gamm4)
library(dplyr)
library(rptR)
library(purrr)
library(ggplot2)

load("Condition modelling output.RData")

## Setting up objects ##

springData <- springobj$data
summerData <- summerobj$data
autumnData <- autumnobj$data

springMod <- gamm4(formula = springobj$mainFormula,
                   data = springData,
                   random = ~ (1|Sett) + (1|Individual))
summerMod <- gamm4(formula = summerobj$mainFormula,
                   data = summerData,
                   random = ~ (1|Sett) + (1|Individual))
autumnMod <- gamm4(formula = autumnobj$mainFormula,
                   data = autumnData,
                   random = ~ (1|Sett) + (1|Individual))

spring.random <- ranef(springMod$mer, condVar = TRUE)
summer.random <- ranef(summerMod$mer, condVar = TRUE)
autumn.random <- ranef(autumnMod$mer, condVar = TRUE)

## Making function to get ICC and variance components ##

get_v_icc <- function(model, season) { 
  
  #beta <- getME(model, "beta")
  sigma <- getME(model$mer, "sigma")
  
  v.fixed <- var(predict(model$mer, re.form = NA)) # v for fixed effects
  v.residual <- sigma^2
  v.individual <-  unname(sqrt(unlist(VarCorr(model$mer)))["Individual"])
  v.sett <- unname(sqrt(unlist(VarCorr(model$mer)))["Sett"])
  
  icc1 <- v.individual/(v.residual + v.sett +  v.individual) # adjusted repeatability 1
  icc2 <- (v.sett + v.individual)/(v.residual + v.sett +  v.individual) # adjusted repeatability 2
  icc3 <- v.individual/(v.residual + v.sett + v.individual + v.fixed) # enhanced repeatability 1
  icc4 <- (v.sett + v.individual)/(v.residual + v.sett + v.individual + v.fixed) # enhanced repeatablity 2
  est <- c(v.residual,  v.sett, v.individual, v.fixed, icc1, icc2, icc3, icc4, season) 
  
  return(est)
}

## Getting parameters to bootstrap on ##

spring <- list(fixed = predict(springMod$mer, re.form = NA))

spring$v.fixed <- var(spring$fixed) #Variance explained by the fixed effects
spring$v.sett  <- unname(sqrt(unlist(VarCorr(springMod$mer)))["Sett"]) #Variance explained by sett
spring$v.individual <- unname(sqrt(unlist(VarCorr(springMod$mer)))["Individual"]) #Variance explained by individual
spring$v.residual   <- getME(springMod$mer, "sigma")^2 #Residual variance

summer <- list(fixed = predict(summerMod$mer, re.form = NA))

summer$v.fixed <- var(summer$fixed) #Variance explained by the fixed effects
summer$v.sett  <- unname(sqrt(unlist(VarCorr(summerMod$mer)))["Sett"]) #Variance explained by sett
summer$v.individual <- unname(sqrt(unlist(VarCorr(summerMod$mer)))["Individual"]) #Variance explained by individual
summer$v.residual   <- getME(summerMod$mer, "sigma")^2 #Residual variance

autumn <- list(fixed = predict(autumnMod$mer, re.form = NA))

autumn$v.fixed <- var(autumn$fixed) #Variance explained by the fixed effects
autumn$v.sett  <- unname(sqrt(unlist(VarCorr(autumnMod$mer)))["Sett"]) #Variance explained by sett
autumn$v.individual <- unname(sqrt(unlist(VarCorr(autumnMod$mer)))["Individual"]) #Variance explained by individual
autumn$v.residual   <- getME(autumnMod$mer, "sigma")^2 #Residual variance

## Bootstrapping ##

Nboot    <- 1000

outcomes <- c("v.residual", "v.sett", "v.individual", "v.fixed", "icc1", "icc2", "icc3", "icc4", "season")
iccStore <- data.frame(matrix(NA, nrow = 0, ncol = length(outcomes)))
colnames(iccStore) <- outcomes

outcomes <- spring.random[[1]] %>% row.names()
sprStore <- data.frame(matrix(NA, nrow = Nboot, ncol = length(outcomes)))
colnames(sprStore) <- outcomes

outcomes <- summer.random[[1]] %>% row.names()
sumStore <- data.frame(matrix(NA, nrow = Nboot, ncol = length(outcomes)))
colnames(sumStore) <- outcomes

outcomes <- autumn.random[[1]] %>% row.names()
autStore <- data.frame(matrix(NA, nrow = Nboot, ncol = length(outcomes)))
colnames(autStore) <- outcomes

set.seed(9453278) #"WildCRU" in numpad

for (i in 1:Nboot) {
  
  # Creating new response variable for ICC bootstrapping: varying all random levels in the model for re-estimation #
  
  springBCI <- spring$fixed + 
    model.matrix(~factor(springData$Individual)-1)%*%rnorm(length(levels(factor(springData$Individual))), 0, sqrt(spring$v.individual)) + #Individual effects
    model.matrix(~factor(springData$Sett)-1)%*%rnorm(length(levels(factor(springData$Sett))), 0, sqrt(spring$v.sett)) + #Sett effeccts
    rnorm(nrow(springData), 0, sqrt(spring$v.residual)) #Residuals
  
  summerBCI <- summer$fixed + 
    model.matrix(~factor(summerData$Individual)-1)%*%rnorm(length(levels(factor(summerData$Individual))), 0, sqrt(summer$v.individual)) + #Individual effects
    model.matrix(~factor(summerData$Sett)-1)%*%rnorm(length(levels(factor(summerData$Sett))), 0, sqrt(summer$v.sett)) + #Sett effeccts
    rnorm(nrow(summerData), 0, sqrt(summer$v.residual)) #Residuals
  
  autumnBCI <- autumn$fixed + 
    model.matrix(~factor(autumnData$Individual)-1)%*%rnorm(length(levels(factor(autumnData$Individual))), 0, sqrt(autumn$v.individual)) + #Individual effects
    model.matrix(~factor(autumnData$Sett)-1)%*%rnorm(length(levels(factor(autumnData$Sett))), 0, sqrt(autumn$v.sett)) + #Sett effeccts
    rnorm(nrow(autumnData), 0, sqrt(autumn$v.residual)) #Residuals
  
  springData$BCI <- springBCI
  summerData$BCI <- summerBCI
  autumnData$BCI <- autumnBCI
  
  # Refitting models for each season with the new response variables #
  
  springBootMod <- gamm4(formula = springobj$mainFormula, 
                         data = springData, 
                         random = ~ (1|Sett) + (1|Individual))
  summerBootMod <- gamm4(formula = summerobj$mainFormula, 
                         data = summerData, 
                         random = ~ (1|Sett) + (1|Individual))
  autumnBootMod <- gamm4(formula = autumnobj$mainFormula, 
                         data = autumnData, 
                         random = ~ (1|Sett) + (1|Individual))
  
  # Storing ICC values for each season's bootstrapped model #
  
  iccStore <- rbind(iccStore, 
                    get_v_icc(model = springBootMod, season = "Spring"),
                    get_v_icc(model = summerBootMod, season = "Summer"),
                    get_v_icc(model = autumnBootMod, season = "Autumn"))
  
  # Now, bootstrapping while fixing random effects (only the residual is allowed to vary) to produce bootstrapped estimates of 
  # individual random effect sizes
  
  
  springBCI <- spring$fixed + 
    model.matrix(~factor(springData$Individual)-1)%*%spring.random[[1]][,1] + #Individual effects - keeping real rather than random
    model.matrix(~factor(springData$Sett)-1)%*%spring.random[[2]][,1] + #Sett effects - keeping real rather than random
    rnorm(nrow(springData), 0, sqrt(spring$v.residual)) #Residuals - this is where variation can come in
  
  summerBCI <- summer$fixed + 
    model.matrix(~factor(summerData$Individual)-1)%*%summer.random[[1]][,1] + #Individual effects - keeping real rather than random
    model.matrix(~factor(summerData$Sett)-1)%*%summer.random[[2]][,1] + #Sett effects - keeping real rather than random
    rnorm(nrow(summerData), 0, sqrt(summer$v.residual)) #Residuals - this is where variation can come in
  
  autumnBCI <- autumn$fixed + 
    model.matrix(~factor(autumnData$Individual)-1)%*%autumn.random[[1]][,1] + #Individual effects - keeping real rather than random
    model.matrix(~factor(autumnData$Sett)-1)%*%autumn.random[[2]][,1] + #Sett effects - keeping real rather than random
    rnorm(nrow(autumnData), 0, sqrt(autumn$v.residual)) #Residuals - this is where variation can come in
  
  springData$BCI <- springBCI
  summerData$BCI <- summerBCI
  autumnData$BCI <- autumnBCI
  
  springBootMod2 <- gamm4(formula = springobj$mainFormula, 
                          data = springData, 
                          random = ~ (1|Sett) + (1|Individual))
  summerBootMod2 <- gamm4(formula = summerobj$mainFormula, 
                          data = summerData, 
                          random = ~ (1|Sett) + (1|Individual))
  autumnBootMod2 <- gamm4(formula = autumnobj$mainFormula, 
                          data = autumnData, 
                          random = ~ (1|Sett) + (1|Individual))
  
  sprStore[i,] <- ranef(springBootMod2$mer, condVar = TRUE)[[1]][,1]
  sumStore[i,] <- ranef(summerBootMod2$mer, condVar = TRUE)[[1]][,1]
  autStore[i,] <- ranef(autumnBootMod2$mer, condVar = TRUE)[[1]][,1]
  
  if (i %% 25 == 0) {
    print(paste(100*i/Nboot %>% round(2), "% done.", sep = ""))
  }
}
