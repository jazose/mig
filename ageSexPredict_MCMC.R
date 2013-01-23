#No required libraries

main=function(){
  test()
}


test=function(){
  nCountry=3;
  nSim=100;
  nAge=4;
  
  #mu=c(-5,0,10);
  #nu=1;
  #phi=c(0.1,0.5,0.9);
  #sigma=c(1,1,1);
  mu=matrix(rnorm(nSim*nCountry),nrow=nSim,ncol=nCountry)
  nu=matrix(1,nrow=nSim,ncol=nCountry)
  phi=matrix(runif(nSim*nCountry),nrow=nSim,ncol=nCountry)
  sigma=matrix(1,nrow=nSim,ncol=nCountry)

  
  #Pretend we only have four age categories, all migrating at equal rates
  maleMatrix=matrix(1/(2*nAge),nrow=nAge,ncol=nCountry);
  femaleMatrix=maleMatrix;
  pop=c(1,1,1);
  rates=matrix(1,nrow=nSim,ncol=nCountry)
  z=list(rates=rates,
         mu=mu,
         phi=phi,
         sigma=sigma,
         nu=nu,
         maleMatrix=maleMatrix,
         femaleMatrix=femaleMatrix,
         pop=pop);
  return(predict_MCMC(z));
}




#USE THIS VERSION OF PREDICT FOR MCMC OUTPUT WITH MANY PARAMETER ESTIMATES
#FOR EACH COUNTRY.
#
#The input z is assumed to be a list of all the things needed to advance prediction
#one time step into the future.
#
#z should have the following components:
#z$rates: a matrix of projected or actual migration rates for all countries at time t
#  Each row in the matrix should represent a single simulated trajectory
#z$mu: a vector of fitted mu_c values for all countries (NO! A MATRIX!!!)
#z$phi: a vector of fitted phi_c values for all countries
#z$sigma: a vector of fitted sigma_c values for all countries
#z$nu: a vector of fitted nu_c values for all countries
#z$maleMatrix: a matrix of age-specific migration rates.
#  The entry in row i, column j of maleMatrix represents the percentage
#  of net migration in country j expected to be coming from males in age group i
#  at time t+1.
#z$femaleMatrix: the equivalent matrix for females.
#  Note that the kth column of maleMatrix and the kth column of femaleMatrix should sum to 1.
#z$pop: a vector of projected populations for all countries at time t+1
#z$truncation: Tells us whether the predictive distribution is truncated.
#  If truncation is NULL, don't truncate.
#  truncation=c(-80,150) would truncate the rate distribution to the given range
#
#
#OUTPUT:
#The output should be a list o with the following components:
#o$rates: a matrix of simulated trajectories of migration rates for all countries at time t+1
#  Each row in the matrix should represent a single simulated trajectory (for all countries)
#o$maleArray: an array of simulated age-specific migration rates.
#  The entry in row i, column j, layer k of maleMatrix represents the simulated
#  percentage of net migration in country j coming from males in age group i
#  in simulation k at time t+1
#o$femaleArray: likewise, but for females
predict_MCMC=function(z){
  
  pop=z$pop;#vector of length nCountries
  rates=z$rates;#Matrix of dim: nSimulations x nCountries
  phi=z$phi;#Matrix of dim: nSimulations x nCountries
  mu=z$mu;#Matrix of dim: nSimulations x nCountries
  sigma=z$sigma;#Matrix of dim: nSimulations x nCountries
  nu=z$nu;#Matrix of dim: nSimulations x nCountries
  maleMatrix=z$maleMatrix;#Matrix of dim: nAgeGroups x nCountries
  femaleMatrix=z$femaleMatrix;#Matrix of dim: nAgeGroups x nCountries
  truncation=z$truncation;#Either NULL or a truncation range
  
  nAgeGroups=dim(maleMatrix)[1];
  nCountries=length(pop);
  nSimulations=dim(rates)[1];
  
  if(is.null(truncation)){
    x=simulateXUntruncated(pop,rates,phi,mu,sigma,nu);
  }
  else{
    x=simulateXTruncated(pop,rates,phi,mu,sigma,nu,truncation);
  }
    
  #The x_{c,t}'s are simulated net migration counts for each country.
  #There's no guarantee that they sum to zero, or that each individual age or sex group
  #will sum to zero.
  #Split up the x's into age and sex, then enforce the restriction.
  
  #Initialize result arrays
  #These are not quite the results, because they'll be counts that must later
  #be converted to rates.
  maleArrayCounts=array(0,dim=c(nAgeGroups, nCountries, nSimulations));
  femaleArrayCounts=array(0,dim=c(nAgeGroups, nCountries, nSimulations));
  
  for(k in 1:nSimulations){
    maleArrayCounts[,,k]=rescaleAll(x[k,],maleMatrix,pop);
    femaleArrayCounts[,,k]=rescaleAll(x[k,],femaleMatrix,pop);
  }
  
  rates=makeNetRates(maleArrayCounts,femaleArrayCounts,pop);
  maleArray=convertToRates(maleArrayCounts,pop);
  femaleArray=convertToRates(femaleArrayCounts,pop);
  
  
  
  return(list(rates=rates,
              maleArray=maleArray,
              femaleArray=femaleArray))
}

######################
#simulateXUntruncated#
######################
simulateXUntruncated=function(pop,rates,phi,mu,sigma,nu){
  nSimulations=dim(rates)[1]
  nCountries=dim(rates)[2]
  kappa=t(pop*(t(rates)*t(phi)+t(mu)*(1-t(phi)))); #Wow. Matrices are a mess.
  #kappa has dimensions nSimulations x nCountries
  theta=t(pop*t(sigma))^2;#matrix of dim nSimulations x nCountries
  
  omega0=theta;
  
  #Check on how we did with those omegas.
  #score=ksScore(psi,omegaQP,theta,nu,nKSSamples);
  #cat("omegaQP score=",score,"\n");
  #score=ksScore(psi,omega0,theta,nu,nKSSamples);
  #cat("omega0 score=",score,"\n");
  
  #Simulate x_{c,t}'s.
  randomTObs=matrix(0,
                    nrow=nSimulations,
                    ncol=nCountries);
  for(i in 1:nSimulations){
    for(j in 1:nCountries){
      randomTObs[i,j]=rt(1,df=nu[i,j]);      
    }
  }
  x=kappa+sqrt(omega0)*(randomTObs)#Matrix of dim: nSimulations x nCountries
  return(x);
}


####################
#simulateXTruncated#
####################
simulateXTruncated=function(pop,rates,phi,mu,sigma,nu,range){
  nSimulations=dim(rates)[1]
  nCountries=dim(rates)[2]
  #Simulate x_{c,t}'s.
  randomTObs=matrix(0,
                    nrow=nSimulations,
                    ncol=nCountries);
  for(i in 1:nSimulations){
    for(j in 1:nCountries){
      thisMean=mu[i,j]*(1-phi[i,j])+phi[i,j]*rates[i,j];
      thisSigma=sigma[i,j];
      thisDF=nu[i,j];
      randomTObs[i,j]=rttrunc(n=1,
                            mean=thisMean,
                            sigma=thisSigma,
                            df=thisDF,
                            lower=range[1],
                            upper=range[2]);
    }
  }
  x=matrix(0,
           nrow=nSimulations,
           ncol=nCountries);
  for(i in 1:nSimulations){
    for(j in 1:nCountries){
      x[i,j]=pop[j]*randomTObs[i,j];
    }
  }
  return(x);
}

#########
#rttrunc#
#########
#Rejection sampling from a truncated t distribution
#
#Inputs:
# n: number of samples
# mean: mean of distribution
# sigma: scale parameter of distribution
# df: degrees of freedom of distribution
# lower, upper: truncation range
rttrunc=function(n,mean,sigma,df,lower,upper){
  result=rep(0,n);
  for(i in 1:n){
    counter=0;
    repeat{
      x=rt(n=1,df=df);
      y=sigma*x+mean;
      if(y>=lower && y<=upper){
        result[i]=y;
        break;
      }
      counter=counter+1
      if(counter==1000){
        cat("mean =",mean,"\n")
        cat("sigma =",sigma,"\n")
        cat("df =",df,"\n")
      }
    }
  }
  return(result);
}

##############
#makeNetRates#
##############
#Inputs:
#maleArrayCounts & femaleArrayCounts:
#  these should have dimensions nAgeGroups x nCountries x nSimulations
#pop: a vector of each country's projected population, of length nCountries
#
#Output:
#A matrix of dimensions nSimulations x nCountries
#  where each row gives net migration rates for all countries in a single simulation
makeNetRates=function(maleArrayCounts,femaleArrayCounts,pop){
  nCountries=dim(maleArrayCounts)[2];
  nSimulations=dim(femaleArrayCounts)[3];
  
  result=array(0,dim=c(nSimulations,nCountries));
  
  for(k in 1:nSimulations){
    maleTotals=colSums(maleArrayCounts[,,k]);
    femaleTotals=colSums(femaleArrayCounts[,,k]);
    totalTotals=maleTotals+femaleTotals;
    result[k,]=totalTotals/pop;
  }
  
  return(result);
}






################
#convertToRates#
################
#Inputs:
#A: an array of net migration counts, either maleArrayCounts or femaleArrayCounts
#  A should have dimensions nAgeGroups x nCountries x nSimulations
#pop: a vector of each country's projected population
#
#output: those counts rescaled to rates by dividing each count by the appropriate country's population.
convertToRates=function(A,pop){
  result=array(0,dim=dim(A))
  nCountries=length(pop);
  for(j in 1:nCountries){
    result[,j,]=A[,j,]/pop[j];
  }
  return(result);
}






############
#rescaleAll#
############
#Inputs:
#x: a vector of x's, i.e. a single set of simulated net migration
#  counts for all countries that may not sum to zero.
#ageMatrix: Assumed to be either maleMatrix or femaleMatrix, a matrix of age-specific
#  migration rates for one of the sexes
#  ageMatrix should have dimensions nAgeGroups x nCountries
#pop: a vector of the projected population of all countries
#
#Output:
#a matrix of dimensions nAgeGroups x nCountries, giving the appropriately rescaled
#  net migration counts for each country and age group.
#
#Sanity check: Each row should sum to zero.
rescaleAll=function(x,ageMatrix,pop){
  nAgeGroups=dim(ageMatrix)[1];
  nCountries=dim(ageMatrix)[2];
  
  #Initialize result storage
  result=array(0,dim=c(nAgeGroups,nCountries));
  
  #Reconstruct country weights (psi)
  totalPop=sum(pop);#scalar
  psi=pop/totalPop;#vector of length nCountries  
  
  for(i in 1:nAgeGroups){
    ageProportions=ageMatrix[i,];
    ageSpecificX=ageProportions*x;
    result[i,]=rescaleVector(ageSpecificX,psi);
  }
  
  return(result);
}



#rescaleVector takes in a vector of predictions that have not yet been fixed to sum to zero
#and a vector of country weights.
#
#It returns a rescaled version of the predictions where the overflow is proportionally
#distributed to all countries.
rescaleVector=function(predictions, weights){
  s=sum(predictions);
  return(predictions-s*weights);
}