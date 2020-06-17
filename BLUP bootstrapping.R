########################################################
# Title: Bootstrapping repeatability and BLUPs
# Name: Julius G. Bright Ross
# Description: Setting up a bootstrapping procedure for 
# my models that will calculate both repeatability 
# (with CIs) and individual intercepts (with spread)
# that can be used to predict lifetime metrics
########################################################

## Clearing workspace and loading packages/data ##

library(lme4)
library(gamm4)
library(dplyr)
library(rptR)
library(purrr)

load("Condition modelling output.RData")

## Setting up objects ##

formula <- summerobj$mainFormula
data    <- summerobj$data

mod     <- gamm4(formula = formula, 
                 data = data, 
                 random = ~ (1|Sett) + (1|Individual))

mod.rr <- ranef(mod$mer, condVar=TRUE)
summary(mod$mer)
summary(mod$gam)

## Making function to get ICC and variance components ##

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

## Getting parameters to bootstrap on ##

fixed <- predict(mod$mer, re.form = NA) #re.form = NA; not to include any random effects 

v.fixed      <- var(fixed) # variance explained by the fixed effects
v.sett       <- unname(sqrt(unlist(VarCorr(mod$mer)))["Sett"]) #variance explained by sett
v.individual <- unname(sqrt(unlist(VarCorr(mod$mer)))["Individual"]) #variance explained by individual
v.residual   <- getME(mod$mer, "sigma")^2 #residual variance

## Actual bootstrapping procedure ##

Nboot    <- 50
outcomes <- mod.rr[[1]] %>% row.names()
Store    <- data.frame(matrix(NA, nrow = Nboot, ncol = length(outcomes)))
colnames(Store) <- outcomes

set.seed(777) # to get the same results every time

for (i in 1:Nboot) {
  
  # Creating new BCI based on the model
  
  BCI <- fixed + 
    model.matrix(~factor(data$Individual)-1)%*%mod.rr[[1]][,1] + # Individual effects - keeping real rather than random
    model.matrix(~factor(data$Sett)-1)%*%mod.rr[[2]][,1] + # Sett effects - keeping real rather than random
    rnorm(nrow(data), 0, sqrt(v.residual)) # residuals - this is where variation can come in
  
  #cor(fixed, data$BCI)  
  
  data$BCI <- BCI
  
  mod <- gamm4(formula = formula, 
               data = data, 
               random = ~ (1|Sett) + (1|Individual))
  
  Store[i,] <- ranef(mod$mer, condVar = TRUE)[[1]][,1]
  
  if (i %% 10 == 0) {
    print(paste(100*i/Nboot %>% round(2), "% done.", sep = ""))
  }
}

avgBootstrapped <- apply(Store, 2, mean)
original        <- mod.rr[[1]][,1]

## What I've done seems to have constrained the estimates of individual random effects towards 0:

plot(avgBootstrapped ~ original)
