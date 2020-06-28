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
  est <- data.frame(v.residual,  v.sett, v.individual, v.fixed, icc1, icc2, icc3, icc4, season) 
  
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

## Bootstrapping--will take forever; it's already run and saved--if you skip to below, it ##
## has code for reading the objects in                                                    ##

Nboot    <- 1000

outcomes <- c("v.residual", "v.sett", "v.individual", "v.fixed", "icc1", "icc2", "icc3", "icc4", "season")
iccStore <- data.frame(matrix(NA, nrow = 0, ncol = length(outcomes)))
names(iccStore) <- outcomes

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
  
  stopifnot(sum(is.na(springData$BCI)) == 0 & sum(is.na(summerData$BCI) == 0) & sum(is.na(autumnData$BCI) == 0))
  
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
  
  stopifnot(sum(is.na(springData$BCI)) == 0 & sum(is.na(summerData$BCI) == 0) & sum(is.na(autumnData$BCI) == 0))
  
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

#save(iccStore, sprStore, sumStore, autStore, file = "Bootstrapping output.RData")

## Now, getting summary statistics for repeatability ##

load("Bootstrapping output.RData")

ggplot(iccStore, aes(x = season, y = icc1)) + geom_boxplot()

bootDict <- list()

sprCI <- iccStore %>% subset(season == "Spring") %>%
  select(icc1) %>% pull() %>% quantile(probs = c(0.025, 0.975))
sumCI <- iccStore %>% subset(season == "Summer") %>%
  select(icc1) %>% pull() %>% quantile(probs = c(0.025, 0.975))
autCI <- iccStore %>% subset(season == "Autumn") %>%
  select(icc1) %>% pull() %>% quantile(probs = c(0.025, 0.975))

bootDict$spring.mean.icc <- iccStore %>% 
  subset(season == "Spring") %>% 
  select(icc1) %>% pull() %>% mean()
bootDict$summer.mean.icc <- iccStore %>% 
  subset(season == "Summer") %>% 
  select(icc1) %>% pull() %>% mean()
bootDict$autumn.mean.icc <- iccStore %>% 
  subset(season == "Autumn") %>% 
  select(icc1) %>% pull() %>% mean()

bootDict$spring.org.icc <- get_v_icc(model = springMod, season = "Spring") %>% select(icc1)
bootDict$summer.org.icc <- get_v_icc(model = summerMod, season = "Summer") %>% select(icc1)
bootDict$autumn.org.icc <- get_v_icc(model = autumnMod, season = "Autumn") %>% select(icc1)
#^These are interestingly, they're all lower than the 25th quantile for the respective 
#bootstrapped season results

### Building models of the relationship between individual random intercepts and lifetime metrics ###

## Setting up lifetime metric vectors ##

sprIndividuals <- row.names(spring.random$Individual)
sumIndividuals <- row.names(summer.random$Individual)
autIndividuals <- row.names(autumn.random$Individual)

load("Prepared badger data_9-4-20.RData")

sprMaxAge         <- c()
sprTotalOffspring <- c()
sprSexHolder      <- c()

for (i in 1:length(sprIndividuals)) {
  ind <- sprIndividuals[i]
  indPos <- which(badgerIndices == ind)
  
  sprMaxAge         <- c(sprMaxAge, max(ageYearly[indPos,]))
  sprTotalOffspring <- c(sprTotalOffspring, sum(parentalMatrix[indPos,]))
  sprSexHolder      <- c(sprSexHolder, badgerSex[indPos])
}

sumMaxAge         <- c()
sumTotalOffspring <- c()
sumSexHolder      <- c()

for (i in 1:length(sumIndividuals)) {
  ind <- sumIndividuals[i]
  indPos <- which(badgerIndices == ind)
  
  sumMaxAge         <- c(sumMaxAge, max(ageYearly[indPos,]))
  sumTotalOffspring <- c(sumTotalOffspring, sum(parentalMatrix[indPos,]))
  sumSexHolder      <- c(sumSexHolder, badgerSex[indPos])
}

autMaxAge         <- c()
autTotalOffspring <- c()
autSexHolder      <- c()

for (i in 1:length(autIndividuals)) {
  ind <- autIndividuals[i]
  indPos <- which(badgerIndices == ind)
  
  autMaxAge         <- c(autMaxAge, max(ageYearly[indPos,]))
  autTotalOffspring <- c(autTotalOffspring, sum(parentalMatrix[indPos,]))
  autSexHolder      <- c(autSexHolder, badgerSex[indPos])
}

## Building models (done separately by sex) ##

sex <- c()
ageBeta <- c()
offBeta <- c()
ageInt  <- c()
offInt  <- c()
season  <- c()

for (i in 1:Nboot) {
  
  # First, spring:
  
  ageModM  <- glm(sprMaxAge[which(sprSexHolder == "Male")] ~ sprStore[i, which(sprSexHolder == "Male")] %>% 
                    unlist() %>% unname(),
                  family = poisson)
  ageModF  <- glm(sprMaxAge[which(sprSexHolder == "Female")] ~ sprStore[i, which(sprSexHolder == "Female")] %>% 
                    unlist() %>% unname(),
                  family = poisson)
  
  ageBeta  <- c(ageBeta, 
                ageModM$coefficients[2] %>% unname(), 
                ageModF$coefficients[2] %>% unname())
  
  ageInt   <- c(ageInt,
                ageModM$coefficients[1] %>% unname(),
                ageModF$coefficients[1] %>% unname())
  
  offModM  <- glm(sprTotalOffspring[which(sprSexHolder == "Male")] ~ sprStore[i, which(sprSexHolder == "Male")] %>% 
                    unlist() %>% unname(),
                  family = poisson)
  offModF  <- glm(sprTotalOffspring[which(sprSexHolder == "Female")] ~ sprStore[i, which(sprSexHolder == "Female")] %>% 
                    unlist() %>% unname(),
                  family = poisson)
  
  offBeta  <- c(offBeta, 
                offModM$coefficients[2] %>% unname(),
                offModF$coefficients[2] %>% unname())
  offInt   <- c(offInt,
                offModM$coefficients[1] %>% unname(),
                offModF$coefficients[1] %>% unname())
  
  sex    <- c(sex, "Male", "Female")
  season <- c(season, "Spring", "Spring")
  
  # Next, summer:
  
  ageModM  <- glm(sumMaxAge[which(sumSexHolder == "Male")] ~ sumStore[i, which(sumSexHolder == "Male")] %>% 
                    unlist() %>% unname(),
                  family = poisson)
  ageModF  <- glm(sumMaxAge[which(sumSexHolder == "Female")] ~ sumStore[i, which(sumSexHolder == "Female")] %>% 
                    unlist() %>% unname(),
                  family = poisson)
  
  ageBeta  <- c(ageBeta, 
                ageModM$coefficients[2] %>% unname(), 
                ageModF$coefficients[2] %>% unname())
  ageInt   <- c(ageInt,
                ageModM$coefficients[1] %>% unname(),
                ageModF$coefficients[1] %>% unname())
  
  offModM  <- glm(sumTotalOffspring[which(sumSexHolder == "Male")] ~ sumStore[i, which(sumSexHolder == "Male")] %>% 
                    unlist() %>% unname(),
                  family = poisson)
  offModF  <- glm(sumTotalOffspring[which(sumSexHolder == "Female")] ~ sumStore[i, which(sumSexHolder == "Female")] %>% 
                    unlist() %>% unname(),
                  family = poisson)
  
  offBeta  <- c(offBeta, 
                offModM$coefficients[2] %>% unname(),
                offModF$coefficients[2] %>% unname())
  offInt   <- c(offInt,
                offModM$coefficients[1] %>% unname(),
                offModF$coefficients[1] %>% unname())
  
  sex    <- c(sex, "Male", "Female")
  season <- c(season, "Summer", "Summer")
  
  # Finally, autumn:
  
  ageModM  <- glm(autMaxAge[which(autSexHolder == "Male")] ~ autStore[i, which(autSexHolder == "Male")] %>% 
                    unlist() %>% unname(),
                  family = poisson)
  ageModF  <- glm(autMaxAge[which(autSexHolder == "Female")] ~ autStore[i, which(autSexHolder == "Female")] %>% 
                    unlist() %>% unname(),
                  family = poisson)
  
  ageBeta  <- c(ageBeta, 
                ageModM$coefficients[2] %>% unname(), 
                ageModF$coefficients[2] %>% unname())
  ageInt   <- c(ageInt,
                ageModM$coefficients[1] %>% unname(),
                ageModF$coefficients[1] %>% unname())
  
  offModM  <- glm(autTotalOffspring[which(autSexHolder == "Male")] ~ autStore[i, which(autSexHolder == "Male")] %>% 
                    unlist() %>% unname(),
                  family = poisson)
  offModF  <- glm(autTotalOffspring[which(autSexHolder == "Female")] ~ autStore[i, which(autSexHolder == "Female")] %>% 
                    unlist() %>% unname(),
                  family = poisson)
  
  offBeta  <- c(offBeta, 
                offModM$coefficients[2] %>% unname(),
                offModF$coefficients[2] %>% unname())
  offInt   <- c(offInt,
                offModM$coefficients[1] %>% unname(),
                offModF$coefficients[1] %>% unname())
  
  sex    <- c(sex, "Male", "Female")
  season <- c(season, "Autumn", "Autumn")
  
}

modelFrame <- data.frame(Coefficient = c(ageBeta, offBeta), 
                         Intercept   = c(ageInt, offInt),
                         Sex = c(sex, sex), 
                         Variable = c(rep("Max age", length(ageBeta)),
                                      rep("Total offspring", length(offBeta))),
                         Season = c(season, season))

## Summarizing and making plots ##

minEff <- min(c(spring.random[[1]][,1], summer.random[[1]][,1], autumn.random[[1]][,1]))
maxEff <- max(c(spring.random[[1]][,1], summer.random[[1]][,1], autumn.random[[1]][,1]))

pointsPerPrediction <- 10
predictionData <- seq(minEff, maxEff, length.out = pointsPerPrediction)

predictionFrame <- data.frame(matrix(nrow = pointsPerPrediction*nrow(modelFrame), 
                                     ncol = ncol(modelFrame) + 4))
names(predictionFrame) <- c("prediction", "group", "average", "bound", names(modelFrame))

# First, we'll just get the prediction line for each set of bootstrapped values #

for (i in 1:nrow(modelFrame)) {
  response <- exp(modelFrame$Intercept[i] + predictionData*modelFrame$Coefficient[i])
  
  temp <- data.frame(prediction  = response,
                     group       = i,
                     average     = "No",
                     bound = "Neither")
  
  temp <- cbind(temp, 
                do.call("rbind", replicate(pointsPerPrediction, 
                                           modelFrame[i,] %>% mutate_if(is.factor, as.character), 
                                           simplify = FALSE)))
  
  posToInsert <- min(which(is.na(predictionFrame$prediction)))
  
  predictionFrame[posToInsert:(posToInsert + pointsPerPrediction - 1),] <- temp
  
}

# Now, we'll add in some special lines: the mean relationship and the 95% quantile relationships #

coefSummary <- modelFrame %>% group_by(Season, Sex, Variable) %>% 
  mutate(LowerInt = quantile(Intercept, probs = c(0.025)),
         MeanInt = mean(Intercept),
         HigherInt = quantile(Intercept, probs = c(0.975)),
         LowerCoef = quantile(Coefficient, probs = c(0.025)),
         MeanCoef = mean(Coefficient),
         HigherCoef = quantile(Coefficient, probs = c(0.975))) %>%
  distinct(Season, Sex, Variable, .keep_all = TRUE)

for (i in 1:nrow(coefSummary)) {
  # First, averages:
  
  response <- exp(coefSummary$MeanInt[i] + predictionData*coefSummary$MeanCoef[i])
  
  temp <- data.frame(prediction  = response,
                     group       = i,
                     average     = "Yes",
                     bound       = "Neither")
  
  temp <- cbind(temp, 
                data.frame(do.call("rbind", replicate(pointsPerPrediction, 
                                           coefSummary[i,] %>% mutate_if(is.factor, as.character) %>% 
                                             select(Coefficient, Intercept, Sex, Variable, Season), 
                                           simplify = FALSE))))

  predictionFrame <- rbind(predictionFrame, temp)
  
  # Now, lower bound:
  
  response <- exp(coefSummary$LowerInt[i] + predictionData*coefSummary$LowerCoef[i])
  
  temp <- data.frame(prediction  = response,
                     group       = i,
                     average     = "No",
                     bound       = "Lower")
  
  temp <- cbind(temp, 
                data.frame(do.call("rbind", replicate(pointsPerPrediction, 
                                                      coefSummary[i,] %>% 
                                                        mutate_if(is.factor, as.character) %>% 
                                                        select(Coefficient, Intercept, Sex, Variable, Season), 
                                                      simplify = FALSE))))
  
  predictionFrame <- rbind(predictionFrame, temp)
  
  # Finally, upper bound:
  
  response <- exp(coefSummary$HigherInt[i] + predictionData*coefSummary$HigherCoef[i])
  
  temp <- data.frame(prediction  = response,
                     group       = i,
                     average     = "No",
                     bound       = "Upper")
  
  temp <- cbind(temp, 
                data.frame(do.call("rbind", replicate(pointsPerPrediction, 
                                           coefSummary[i,] %>% 
                                             mutate_if(is.factor, as.character) %>% 
                                             select(Coefficient, Intercept, Sex, Variable, Season), 
                                           simplify = FALSE))))
  
  predictionFrame <- rbind(predictionFrame, temp)
  
}

predictionFrame$average[which(predictionFrame$average == "1")] <- "No"
predictionFrame$bound[which(predictionFrame$bound == "1")] <- "Neither"

predictionFrame$x <- rep(predictionData, nrow(predictionFrame)/pointsPerPrediction)

## Pulling out metrics for reporting ##

bootDict$autFemaleMeanMaxOff <- predictionFrame %>% 
  subset(Sex == "Female" & Season == "Autumn" & Variable == "Total offspring" & average == "Yes") %>%
  select(prediction) %>% pull() %>% max() %>% round(2)

bootDict$autFemaleMeanMinOff <- predictionFrame %>% 
  subset(Sex == "Female" & Season == "Autumn" & Variable == "Total offspring" & average == "Yes") %>%
  select(prediction) %>% pull() %>% min() %>% round(2)

bootDict$autFemaleLowerMaxOff <- predictionFrame %>% 
  subset(Sex == "Female" & Season == "Autumn" & Variable == "Total offspring" & bound == "Lower") %>%
  select(prediction) %>% pull() %>% max() %>% round(2)

bootDict$autFemaleLowerMinOff <- predictionFrame %>% 
  subset(Sex == "Female" & Season == "Autumn" & Variable == "Total offspring" & bound == "Lower") %>%
  select(prediction) %>% pull() %>% min() %>% round(2)

bootDict$autFemaleUpperMaxOff <- predictionFrame %>% 
  subset(Sex == "Female" & Season == "Autumn" & Variable == "Total offspring" & bound == "Upper") %>%
  select(prediction) %>% pull() %>% max() %>% round(2)

bootDict$autFemaleUpperMinOff <- predictionFrame %>% 
  subset(Sex == "Female" & Season == "Autumn" & Variable == "Total offspring" & bound == "Upper") %>%
  select(prediction) %>% pull() %>% min() %>% round(2)

bootDict$sprMaleMeanMaxOff <- predictionFrame %>% 
  subset(Sex == "Male" & Season == "Spring" & Variable == "Total offspring" & average == "Yes") %>%
  select(prediction) %>% pull() %>% max() %>% round(2)

bootDict$sprMaleMeanMinOff <- predictionFrame %>% 
  subset(Sex == "Male" & Season == "Spring" & Variable == "Total offspring" & average == "Yes") %>%
  select(prediction) %>% pull() %>% min() %>% round(2)

bootDict$sprMaleLowerMaxOff <- predictionFrame %>% 
  subset(Sex == "Male" & Season == "Spring" & Variable == "Total offspring" & bound == "Lower") %>%
  select(prediction) %>% pull() %>% max() %>% round(2)

bootDict$sprMaleLowerMinOff <- predictionFrame %>% 
  subset(Sex == "Male" & Season == "Spring" & Variable == "Total offspring" & bound == "Lower") %>%
  select(prediction) %>% pull() %>% min() %>% round(2)

bootDict$sprMaleUpperMaxOff <- predictionFrame %>% 
  subset(Sex == "Male" & Season == "Spring" & Variable == "Total offspring" & bound == "Upper") %>%
  select(prediction) %>% pull() %>% max() %>% round(2)

bootDict$sprMaleUpperMinOff <- predictionFrame %>% 
  subset(Sex == "Male" & Season == "Spring" & Variable == "Total offspring" & bound == "Upper") %>%
  select(prediction) %>% pull() %>% min() %>% round(2)

## Exporting stuff ##

save(predictionFrame, coefSummary, sprStore, sumStore, autStore, bootDict, file = "Bootstrapping output.RData")
