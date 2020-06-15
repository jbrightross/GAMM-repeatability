library(lme4)
library(gamm4)
library(dplyr)
library(rptR)

rm(list = ls())
setwd("C:/Users/jbrig/Desktop/Badgers/Badger condition/Processed datasets")
load("Condition modelling output.RData")

## Just a quick illustration of the issues with the simulate.merMod() function ##

set.seed(0) 
dat <- gamSim(1,n=400,scale=2) ## simulate 4 term additive truth
# Now add 20 level random effect `fac'...
dat$fac <- fac <- as.factor(sample(1:20,400,replace=TRUE))
dat$y <- dat$y + model.matrix(~fac-1)%*%rnorm(20)*.5
br <- gamm4(y~s(x0)+x1+s(x2),data=dat,random=~(1|fac))

lme4:::simulate.merMod(br$mer)

## Setting all our function values for rptR, using the summer model 
## (without splines) as an example 

formula <- summerobj$mainFormula
grname <- c("Sett", "Individual", "Fixed", "Residual") #Terms to estimate repeatability for
nboot <- 100 #Number of bootstrapping iterations
data <- summerobj$data
ratio <- TRUE
adjusted <- TRUE #Adjusts for the variance taken by fixed effects
CI <- 0.95
npermut <- 100 #Number of permutations to get a p-value (not necessary, but needed for code)
update <- FALSE #Just has to do with updating existing objects
rptObj <- NULL #Same as above
ncores <- NULL #We're not doing parallel stuff--no extra cores on my computer
parallel <- FALSE

## Checks to make sure no NA values in the data (not a problem) ##

no_NA_vals <- stats::complete.cases(data[all.vars(formula)])

if (sum(!no_NA_vals) > 0) {
  warning(paste0(sum(!no_NA_vals), " rows containing missing values were removed"))
  data <- data[no_NA_vals, ]
}

## Setting up the primary model ##

mod <- gamm(formula = formula, 
             data = data, 
             random = list(Sett = ~1, Individual = ~1))

## Now we set up some holders for other kinds of variance ##

grname_org      <- grname
output_resid    <- FALSE
output_overdisp <- FALSE
output_fixed    <- FALSE

for (component in c("Residual", "Overdispersion", "Fixed")) {
  # Pulling out those preset types
  if (any(grname == component)) {
    grname <- grname[-which(grname == component)]
    if (component == "Residual") 
      output_resid <- TRUE
    if (component == "Overdispersion") 
      output_overdisp <- TRUE
    if (component == "Fixed") 
      output_fixed <- TRUE
  }
}

## Getting which terms are random, to subset on ##

terms <- attr(terms(formula), "term.labels")

randterms <- terms[which(regexpr(" | ", terms, perl = TRUE) > 
                           0)]

## Making the main function for making a model with bootstrapped data and extracting 
## variance types

R_pe <- function(formula, data, grname, mod = NULL, resp = NULL) {
  
  # If we're specifying "mod" (for non-default), we refit it with the new data
  if (!is.null(mod)) {
    # With new response variable
    newdata <- data
    newdata$BCI <- resp
    mod <- gamm(formula = formula, 
                data = newdata, 
                random = list(Sett = ~1, Individual = ~1))
  }
  else {
    # Otherwise, setting up our default
    mod <- gamm(formula = formula, 
                data = data, 
                random = list(Sett = ~1, Individual = ~1))
  }
  
  # Here's our variance object 
  
  VarComps <- lme4::VarCorr(mod$lme)
  
  # This is the residual variance
  var_e <- as.numeric(VarComps[which(row.names(VarComps) == "Residual"), 1])
  names(var_e) <- "Residual"
  
  # This is the overdispersion variance? Apparently same as the residual
  var_o <- var_e
  names(var_o) <- "Overdispersion"
  
  # This calculates the variance of predictions without the random effects (variance taken up
  # by the fixed effects)
  var_f <- stats::var(stats::predict(mod$lme, re.form = NA)) 
  names(var_f) <- "Fixed"
  
  # This is the variance of each random effect
  
  var_a <- c(VarComps[which(grepl(pattern = grname[1], row.names(VarComps))) + 1, 1],
             VarComps[which(grepl(pattern = grname[2], row.names(VarComps))) + 1, 1]) %>% as.numeric()
  names(var_a) <- grname
  
  # Sum them and the residual variance--this is our denominator ##
  var_p <- sum(var_a) + var_e %>% unname()
  
  # Okay, so now, our default is that we adjust for the variance taken up by fixed effects
  if (!adjusted) 
    # But if we don't want to adjust, we add the variance taken by fixed effects to the denominator
    var_p <- var_p + var_f
  
  # If, for some reason, we just want variances, not repeatabilities:
  if (ratio == FALSE) {
    # Then we'll just get a variance dataframe back
    R <- as.data.frame(t(var_a))
    names(R) <- grname
    
    # And add in a thing for the residual's information
    if (output_resid) {
      R$Residual <- var_e
    }
    
    # And add in info for overdispersion parameter (?)
    if (output_overdisp) {
      R$Overdispersion <- var_o
    }
    
    # And info for the fixed variance
    if (output_fixed) {
      R$Fixed <- var_f
    }
    
    return(R)
  }
  
  # Now, if we're normal and do want the repeatability:
  if (ratio == TRUE) {
    # Repeatability is each variance term divided by the total variance
    R <- var_a/var_p
    R <- as.data.frame(t(R))
    names(R) <- grname
    
    # And we're adding info for each of the extra terms
    if (output_resid) {
      R$Residual <- var_e/var_p
    }
    if (output_overdisp) {
      R$Overdispersion <- var_e/var_p
    }
    if (output_fixed) {
      R$Fixed <- var_f/var_p
    }
    return(R)
  }
}

# R is our basic repeatability dataframe first, for the full data (not bootstrapped yet)
R <- R_pe(formula, data, grname)

## Now the bootstrapping. We start out by simulating response variables with the model structure ##
## from our fitted mixed effects model--that means

if (nboot > 0) {
  # Although we're specifying library stats, the function knows it's a merMod object, so it should 
  # take on that model structure for the prediction
  Ysim <- as.matrix(lme4:::simulate.merMod(mod$lme, nsim = nboot))
}
