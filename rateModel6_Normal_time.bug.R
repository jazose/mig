model {  
  for (i in 1:1){
    for(j in 1:c){
      y[i,j] ~ dnorm(y.hat[i,j], tau[i])
      y.hat[i,j] <- beta[i]+phi[i]*x[i,j]
    }
    tau[i] ~ dgamma(a, b)
    beta[i] <- mu[i]*(1-phi[i])
    mu[i] ~ dnorm(lambda.mu, tau.mu)
    phi[i] ~ dunif(0,1)
    logDummy[i] <- 0
    dummy[i] <- exp(logDummy[i])
  }
  for (i in 2:r){
    for(j in 1:c){
      y[i,j] ~ dnorm(dummy[i]*y.hat[i,j], tau[i]/dummy[i]^2) #For all times other than the first, introduce a dummy
      y.hat[i,j] <- beta[i] + phi[i]*x[i,j]
    }
    tau[i] ~ dgamma(a, b)
    beta[i] <- mu[i]*(1-phi[i])
    mu[i] ~ dnorm(lambda.mu, tau.mu)
    phi[i] ~ dunif(0,1)
    logDummy[i] ~ dnorm(0,0.0001)
    dummy[i] <- exp(logDummy[i])
  }
  
  b ~ dunif(0,100*(a-1))
  a ~ dunif(1,10)
  lambda.mu ~ dunif(-100,100)
  tau.mu <- pow(sigma.mu,-2)
  sigma.mu ~ dunif(0,100)
}