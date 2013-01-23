model {
  #Make a categorical probability vector for nu.  
  pi[1] <- 1
  pi[2] <- 1
  pi[3] <- 1
  pi[4] <- 1
  pi[5] <- 1
  pi[6] <- 1
  pi[7] <- 1
  pi[8] <- 1
  nu ~ dcat(pi)
  
  for (i in 1:r){
    for(j in 1:c){
      y[i,j] ~ dt(y.hat[i,j], tau[i],nu)
      y.hat[i,j] <- beta[i]+phi[i]*x[i,j]
    }
    tau[i] ~ dgamma(a, b)
    beta[i] <- mu[i]*(1-phi[i])
    mu[i] ~ dnorm(lambda.mu, tau.mu)
    phi[i] ~ dunif(0,1)
  }
  
  b ~ dunif(0,100*(a-1))
  a ~ dunif(1,10)
  lambda.mu ~ dunif(-100,100)
  tau.mu <- pow(sigma.mu,-2)
  sigma.mu ~ dunif(0,100)
}