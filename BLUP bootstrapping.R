########################################################
# Title: Bootstrapping BLUPs
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
library(ggplot2)

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

## Bootstrapping while fixing random effects ##

Nboot    <- 1000
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
  
  if (i %% 25 == 0) {
    print(paste(100*i/Nboot %>% round(2), "% done.", sep = ""))
  }
}

avgBootstrapped <- apply(Store, 2, mean)
sdBootstrapped  <- apply(Store, 2, sd)
original        <- mod.rr[[1]][,1]

numObs <- c()

for (i in 1:length(original)) {
  ind <- names(avgBootstrapped)[[i]]
  
  numObs <- c(numObs, length(which(data$Individual == ind)))
}

## What I've done seems to have constrained the estimates of individual random effects towards 0:

plot(avgBootstrapped ~ original)

## Interestingly, the extra line that's visible close to slope = 0 is individuals with only
## one observation:
plot(avgBootstrapped[which(numObs > 1)] ~ original[which(numObs > 1)])

## Although, it also seems to affect estimates less when we have more observations for
## a given individual: 

plot(I(abs(original - avgBootstrapped)) ~ numObs, 
     ylab = "Bootstrapped difference",
     xlab = "Number of observations")

## And seems to provide a tendency towards 0 for individuals with few observations but towards 
## an actual estimate with more observations

plot(avgBootstrapped ~ numObs, 
     ylab = "Bootstrapped average",
     xlab = "Number of observations")

## Finally, there's an oddly small variance for individuals with only one capture, but after that
## seems like a decrease in variance the more captures an individual has:

plot(sdBootstrapped ~ numObs, 
     ylab = "Bootstrapped variance",
     xlab = "Number of observations")

## Let's try modelling a bit. First, load data and calculate some lifetime metrics ##

load("Prepared badger data_9-4-20.RData")

maxAge         <- c()
totalOffspring <- c()
sexHolder      <- c()

for (i in 1:length(original)) {
  ind <- names(avgBootstrapped)[i]
  indPos <- which(badgerIndices == ind)
  
  maxAge         <- c(maxAge, max(ageYearly[indPos,]))
  totalOffspring <- c(totalOffspring, sum(parentalMatrix[indPos,]))
  sexHolder      <- c(sexHolder, badgerSex[indPos])
}

## Some models (done separately by sex) ##

sex <- c()
ageBeta <- c()
offBeta <- c()

for (i in 1:nrow(Store)) {
  ageModM  <- lm(maxAge[which(sexHolder == "Male")] ~ Store[i, which(sexHolder == "Male")] %>% 
                   unlist() %>% unname())
  ageModF  <- lm(maxAge[which(sexHolder == "Female")] ~ Store[i, which(sexHolder == "Female")] %>% 
                   unlist() %>% unname())
  
  ageBeta  <- c(ageBeta, 
                ageModM$coefficients[2] %>% unname(), 
                ageModF$coefficients[2] %>% unname())
  
  offModM  <- glm(totalOffspring[which(sexHolder == "Male")] ~ Store[i, which(sexHolder == "Male")] %>% 
                   unlist() %>% unname(),
                  family = poisson)
  offModF  <- glm(totalOffspring[which(sexHolder == "Female")] ~ Store[i, which(sexHolder == "Female")] %>% 
                   unlist() %>% unname(),
                  family = poisson)
  
  offBeta  <- c(offBeta, 
                offModM$coefficients[2] %>% unname(),
                offModF$coefficients[2] %>% unname())
  
  sex <- c(sex, "Male", "Female")
  
}

modelFrame <- data.frame(Coefficient = c(ageBeta, offBeta), 
                         Sex = c(sex, sex), 
                         Variable = c(rep("Max age", length(ageBeta)),
                                      rep("Total offspring", length(offBeta))))

ggplot(modelFrame, aes(x = Variable, y = Coefficient, colour = Sex)) + 
  geom_boxplot()

quantile(offBeta[sex == "Male"], probs = c(0.025, 0.975))
pnorm(0, mean = mean(offBeta[sex == "Male"]), sd = sd(offBeta[sex == "Male"]))

## Comparing to results of modelling original random effects ##

summary(glm(totalOffspring[which(sexHolder == "Male")] ~ original[which(sexHolder == "Male")], family = poisson))
mean(offBeta[which(sex == "Male")])
summary(glm(totalOffspring[which(sexHolder == "Female")] ~ original[which(sexHolder == "Female")], family = poisson))
mean(offBeta[which(sex == "Female")])

numObsBin <- rep(NA, length(numObs))
numObsBin[which(numObs == 1)] <- "1"
numObsBin[which(1 < numObs & numObs < 4)] <- "2-3"
numObsBin[which(3 < numObs & numObs < 6)] <- "4-5"
numObsBin[which(numObs > 5)] <- "6+"

originalFrame <- data.frame(TotalOffspring = totalOffspring, 
                            Sex = sexHolder, 
                            RandEffect = original, 
                            NumObservations = numObsBin)

ggplot(originalFrame, aes(x = RandEffect, y = TotalOffspring, colour = Sex)) + 
  geom_point() + facet_wrap(.~NumObservations, nrow = 2)

## Alternative bootstrapping without fixing random effects ##

Nboot    <- 50
outcomes <- mod.rr[[1]] %>% row.names()
Store    <- data.frame(matrix(NA, nrow = Nboot, ncol = length(outcomes)))
colnames(Store) <- outcomes

set.seed(777) # to get the same results every time

for (i in 1:Nboot) {
  
  # Creating new BCI based on the model
  
  BCI <- fixed + 
    model.matrix(~factor(data$Individual)-1)%*%rnorm(length(levels(factor(data$Individual))), 0, sqrt(v.individual)) + # Individual effects
    model.matrix(~factor(data$Sett)-1)%*%rnorm(length(levels(factor(data$Sett))), 0, sqrt(v.sett)) + # Sett effeccts
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

## This seems to have completely randomized it, so there's no longer any link between the original 
## random effect value and the bootstrapped average

plot(avgBootstrapped ~ original)

