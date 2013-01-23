model {
  for (i in 1:r){
    for(j in 1:c){
      y[i,j] ~ dt(y.hat[i,j], tauPrime[i],nu[i])
      y.hat[i,j] <- beta[i]+phi[i]*x[i,j]
    }
    nu[i] ~ dexp(0.05) I(0.6689,)#Prior dist. only allows ratio of 95% to 80% width of 8:1
    ciWidth99[i] <- 2*qt(0.995,0,1,nu[i])
    sigmaMax[i] <- 200/ciWidth99[i]
    tauMin[i] <- sigmaMax[i]^(-2)
    tauPrime[i] <- tau[i] + tauMin[i]
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