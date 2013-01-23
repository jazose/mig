#Not a good idea. Leads to unnecessarily wide confidence intervals

model {
  k1 <- -2.95551
  k2 <- 27.23162
  asymptoteValue <- 3.919928
  for (i in 1:r){
    for(j in 1:c){
      y[i,j] ~ dt(y.hat[i,j], tau[i],nu[i])
      y.hat[i,j] <- beta[i]+phi[i]*x[i,j]
    }
    nu[i] ~ dexp(0.05) I(0.6689,)#Prior dist. only allows ratio of 95% to 80% width of 8:1
    tau[i] <- pow(sigmaPrime[i],-2)
    ciWidth99[i] <- 2*qt(0.995,0,1,nu[i])
    sigmaMax[i] <- 200/ciWidth99[i]
    sigmaPrime[i] <- sigmaMax[i]*sigma[i]    
    sigma[i] ~ dbeta(A,B)
    beta[i] <- mu[i]*(1-phi[i])
    mu[i] ~ dnorm(lambda.mu, tau.mu)
    phi[i] ~ dunif(0,1)
  }
  
  A ~ dunif(0,10)
  B ~ dunif(0,10)
  lambda.mu ~ dunif(-100,100)
  tau.mu <- pow(sigma.mu,-2)
  sigma.mu ~ dunif(0,100)
}