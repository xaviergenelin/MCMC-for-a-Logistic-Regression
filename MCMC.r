# Setup
library(MASS)
fgData <- read.csv("kickData.csv")
burnin <- 10000
kickDistance <- fgData$kick_distance
kickResult <- fgData$kick_result

# Functions
# Expit function to create a probability of making a kick given X distance
expit = function(B0, B1, X){
  return(1 / (1 + exp(-B0 - B1 * X)))
}

# Posterior function
posterior = function(B0, B1, X, Y){
  # Calculate the probability of making a kick with X distances 
  probs = expit(B0, B1, X = X)
  # Calculate the log likelihood of making a kick given the probability calculated and the Y value of the kick result
  logLikelihood = sum(log(dbinom(x = Y, size = 1, prob = probs)))
  # Calculate the log prior values for both Beta0 and Beta1 used in the probability and logLikelihood calculations
  logPriors = sum(log(dnorm(B0, mean = 0, sd = 15)) + log(dnorm(B1, mean = 0, sd = 15)))
  # Calculate the posterior value from the logPriors and logLikelihood
  logLikelihood + logPriors
}

# Sampler Function
betaSampler = function(X, Y, niter, B0StartVal, B1StartVal, proposalsd){
  
  # Create an empty data frame to store the results from the MCMC with niter rows and 2 columns
  betas = data.frame(matrix(nrow = niter, ncol = 2))
  # Create an empty vector to store the selected Beta0 and Beta1 values
  B0 = rep(0, niter)
  B1 = rep(0, niter)
  
  #use the starting value as the first Beta0 and Beta1
  B0[1] = B0StartVal
  B1[1] = B1StartVal
  
  for(i in 2:niter){
    # Temporary value for current beta0 and current beta1 is previous Beta0 and Beta1 respectively
    currentB0 = B0[i-1]
    currentB1 = B1[i-1]
    
    # Get the next random draw add in a random value from the jumping distr. to our Beta0 value
    newB0 = currentB0 + rnorm(1,0,proposalsd)
    
    # Calculate the candidate posterior for the new Beta0 value with the current Beta1 value
    candidatePosterior = posterior(B0 = newB0, B1 = currentB1, X = X, Y = Y)
    # Calculate the current posterior for the current Beta0 value with the current Beta1 value
    currentPosterior = posterior(B0 = currentB0, B1 = currentB1, X = X, Y = Y)
    
    #Find r ratio for beta0: 
    r = candidatePosterior - currentPosterior
    
    #accept this new value with prob min(1,r), leave value the same with prob 1-r
    # Store the selected value
    if(log(runif(1))<r){
      B0[i] = newB0       # accept move with probability min(1,r)
    } else {
      B0[i] = currentB0        # otherwise "reject" move, and stay where we are
    } #end if statement for Beta0
    
    # Save the selected beta0 from above as selectedB0 to use in selecting Beta1
    selectedB0 = B0[i]
    
    # Repeat the selection process for Beta1
    #get the next random draw add in a random value from the jumping distr. to our Beta1 value
    newB1 = currentB1 + rnorm(1,0,proposalsd)
    
    # Calculate the candidate posterior for the new Beta1 value using the Beta0 value we selected earlier
    candidatePosterior = posterior(B0 = selectedB0, B1 = newB1, X = X, Y = Y)
    # Calculate the current posterior for the current Beta1 value using the Beta0 value we selected earlier
    currentPosterior = posterior(B0 = selectedB0, B1 = currentB1, X = X, Y = Y)
    
    #Find r ratio for Beta1: 
    r = candidatePosterior - currentPosterior
    
    #accept this new value with prob min(1,r), leave value the same with prob 1-r
    if(log(runif(1))<r){
      B1[i] = newB1       # accept move with probability min(1,r)
    } else {
      B1[i] = currentB1        # otherwise "reject" move, and stay where we are
    } #end if statement for Beta1
  } #end loop 
  
  # Store the Beta0 values in the first column of our betas dataframe
  betas[,1] = B0
  # Store the Beta1 values in the second column of our betas dataframe
  betas[,2] = B1
  return(betas)
}

# Run Sampler
betaValues = betaSampler(X = kickDistance, Y = kickResult, niter = 100000, B0StartVal = 0, B1StartVal = 0, proposalsd = 0.2)

# Calculations
# Calculate the mean, standard deviation, and a 95% credible interval for Beta0 after using the burn-in
mean(betaValues[,1][-c(1:burnin)])
# 6.690513
sd(betaValues[,1][-c(1:burnin)])
# 0.3646577
quantile(x = betaValues[,1][-c(1:burnin)], c(0.025, 0.975))
# (6.011510, 7.300007)

# Calculate the mean, standard deviation, and a 95% credible interval for Beta1 after using the burn-in
mean(betaValues[,2][-c(1:burnin)])
# -0.1214335
sd(betaValues[,2][-c(1:burnin)])
# 0.007995718
quantile(x = betaValues[,2][-c(1:burnin)], c(0.025, 0.975))
# (-0.1349399, -0.1068982)

# Maximum likelihood model
fit <- glm(kick_result ~ kick_distance, family = binomial(link = 'logit'), data = fgData)
summary(fit)
confint(fit)

# Results from maximum likelihood model

# Beta0
# Estimate = 6.692524
# Standard Error = 0.377172
# 95% Confidence Interval (5.9737159, 7.453366)

# Beta1
# Estimate = -0.121511
# Standard Error = 0.008276
# 95% Confidence Interval = (-0.1381169, -0.105649)