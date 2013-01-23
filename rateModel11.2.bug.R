model {
  for (i in 1:r){
    for(j in 1:c){
      y[i,j] ~ dt(y.hat[i,j], tau[i],nu[i]) T(-80,150)
      y.hat[i,j] <- beta[i]+phi[i]*x[i,j]
    }
    nu[i] ~ dexp(0.05) T(1,)
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