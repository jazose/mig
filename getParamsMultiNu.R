#This file is set up so that you can just execute all the code once to read in
#functions and then make a call to main.


#Required libraries
library(quadprog)
library(gplots)
library(hett)
library(arm)
library(car)
library(plotrix)
library(rjags)
library(coda)



main=function(){
  #setwd("C:/Users/jonazose/Dropbox/RA/Code/stateOfTheArt")
  setwd("C:/Users/jon-everyday/Dropbox/RA/Code/faloorum")
  
  setup();
  
  hierModelMultiNu("rateModel11.1.bug.R",rateData,ppc=8,"Hier111_3Exc");
}


hierModelMultiNu=function(bugFile,data,ppc=10,outFile){
  #Sample call: hierModelMultiNu("rateModel7.3.bug.R",rateData,ppc=10,"Hier73_1Exc")
  xMatrix=as.matrix(data[,2:(2+ppc-1)]);
  rownames(xMatrix)=NULL;
  colnames(xMatrix)=NULL;
  
  yMatrix=as.matrix(data[,3:(3+ppc-1)]);
  rownames(yMatrix)=NULL;
  colnames(yMatrix)=NULL;  
  
  r = nrow(xMatrix);
  c = ncol(xMatrix);
  
  jags <- jags.model(bugFile,
                     data = list('x' = xMatrix,
                                 'y' = yMatrix,
                                 'r' = r,
                                 'c' = c),
                     n.chains = 4,
                     n.adapt = 1000)
  
  update(jags, 100)
  
  samples=coda.samples(jags,
                       c('mu','nu','phi','tau'),
                       10000,
                       thin=25)
  
  #raftery.diag: Need 3746 samples, thinning interval ~25, burn in of ~100
  
  #Goal: Write out all the param vectors one after another.
  #First, say how many vectors to expect and how many params in each.
  nParams=r*4;
  nEstimates=dim(samples[[1]])[1];
  for(i in 1:length(samples)){
    thisOutFile=paste("./PE/PE_",outFile,"_Chain",i,".txt",sep="");
    write(c(nParams,nEstimates),thisOutFile)
    write(t(samples[[i]]),thisOutFile,append=TRUE)
  }
  
  #return(samples)

}



##########################
#LESS USED FUNCTIONS BELOW THIS POINT
##########################


nameLookup=function(countryNames,countryCode){
  return(toString(countryNames[which(countryNames$V2==countryCode),1]))
}



setup=function(){
  
  countData<<-read.csv("migrationData.csv")
  countData<<-countData[which(countData$Country.code<900),] #Get rid of regions
  
  #Important indices: 54 (Sierra Leone), 59 (North Korea), 61 (Mongolia), 176 (Ecuador)
  
  r=dim(countData)[1];
  c=dim(countData)[2];
  
  rateData<<-read.csv("migrationRateData.csv")
  rateData<<-rateData[which(rateData$Country.code<900),] #Get rid of regions
  
  #Apply some jitter to the rate data.
  rateData<<-rateData+runif(r*c,min=-0.0005,max=0.0005);
  
  #Allow for easy lookup of country names from codes
  countryNames<<-read.csv("countryNames.csv",header=FALSE, sep=";")
  countryNames<<-countryNames[which(countryNames$V2<900),] #Get rid of regions  
}

rewrite=function(fileName){
  paramEstimates=scan(fileName);
  n=(length(paramEstimates)-1)/3;
  muVector=paramEstimates[1:n];
  nu=paramEstimates[(n+1)];
  phiVector=paramEstimates[(n+2):(2*n+1)];
  tauVector=paramEstimates[(2*n+2):(3*n+1)];
  
  out=rep(0,4*n);
  out[1:n]=muVector;
  out[(n+1):(2*n)]=rep(nu,n);
  out[(2*n+1):(3*n)]=phiVector;
  out[(3*n+1):(4*n)]=tauVector;
  write(out,fileName)
}

#Same model 4 as in 07122012 report.
#Countries share a single value for phi, t-distributed residuals
model4=function(data,nu,ppc=8){
  #ppc is the number of usable data points per country.
  #The first data point doesn't count as usable since we condition on it.
  #ppc should be 10 if we're excluding only the most recent data point.
  #ppc should be 8 to exclude the 3 most recent data points.
  
  n=nrow(data);
  
  logLikOfPhi=function(phi){
    maxLogLik=0;
    logLikVector=rep(0,n)
    for(i in 1:n){
      x=as.numeric(data[i,2:(2+ppc-1)]);
      y=as.numeric(data[i,3:(3+ppc-1)]);
      
      logLikOfMuSigma=function(params){
        #params should be a vector of [mu, sigma]
        mu=params[1];
        sigma=params[2];
        r=(y-mu)-phi*(x-mu);
        k=length(r);
        result=k*lgamma((nu+1)/2);
        result=result-k*lgamma(nu/2);
        result=result-k/2*log(pi*nu*sigma^2);
        result=result-((nu+1)/2)*sum(log(1+r^2/(nu*sigma^2)));
        return(result);
      }
      
      thisLogLik=optim(par=c(mean(x),sd(x)),
                       fn=logLikOfMuSigma,
                       control=list(fnscale=-1))$value;
      logLikVector[i]=thisLogLik;
      maxLogLik=maxLogLik+thisLogLik;
      
    }
    return(maxLogLik);
  }
  
  result=optimize(logLikOfPhi,
                  interval=c(-1,1),
                  maximum=TRUE);
  
  bestPhi=result$maximum;
  logLikVector=rep(0,n)
  muVector=rep(0,n);
  sigma2Vector=rep(0,n);
  for(i in 1:n){
    x=as.numeric(data[i,2:(2+ppc-1)]);
    y=as.numeric(data[i,3:(3+ppc-1)]);
    
    logLikOfMuSigma=function(params){
      #params should be a vector of [mu, sigma]
      mu=params[1];
      sigma=params[2];
      r=(y-mu)-bestPhi*(x-mu);
      k=length(r);
      result=k*lgamma((nu+1)/2);
      result=result-k*lgamma(nu/2);
      result=result-k/2*log(pi*nu*sigma^2);
      result=result-((nu+1)/2)*sum(log(1+r^2/(nu*sigma^2)));
      return(result);
    }
    
    model=optim(par=c(mean(x),sd(x)),
                fn=logLikOfMuSigma,
                control=list(fnscale=-1));
    logLikVector[i]=model$value;
    muVector[i]=model$par[1];
    sigma2Vector[i]=(model$par[2])^2
  }
  
  #Write parameters to file
  paramEstimates=rep(0,3*n+1)
  paramEstimates[1:n]=muVector;
  paramEstimates[n+1]=nu;
  paramEstimates[(n+2):(2*n+1)]=rep(bestPhi,n);
  paramEstimates[(2*n+2):(3*n+1)]=sigma2Vector^(-1);
  write(paramEstimates,"junk.txt")#Order of parameters: mus, nu, phis, and taus
  
  
  return(list(nu=nu,
              phi=bestPhi,
              muVector=muVector,
              sigma2Vector=sigma2Vector))
}

mod4SensitivityAnalysis=function(n){
  rateData=read.csv("migrationRateData.csv")
  rateData=rateData[which(rateData$Country.code<900),] #Get rid of regions
  r=dim(rateData)[1];
  c=dim(rateData)[2];      
  
  #Allocate storage space for variables we're interested in tracking across runs.
  #The entry in row i, column j gives the results for country i's parameter in trial j.
  muEstimates=matrix(0,nrow=r,ncol=n);
  sigma2Estimates=matrix(0,nrow=r,ncol=n);
  phiEstimates=rep(0,n);
  
  for(i in 1:n){
    cat(i,"\n");
    #Apply some jitter to the rate data.
    jitteredRateData=rateData+runif(r*c,min=-0.0005,max=0.0005);
    
    #Fit model 4 and hold onto the variables.
    mod=model4(jitteredRateData,2);
    muEstimates[,i]=mod$mu;
    sigma2Estimates[,i]=mod$sigma2;
    phiEstimates[i]=mod$phi;
  }
  return(list(muEstimates=muEstimates,
              sigma2Estimates=sigma2Estimates,
              phiEstimates=phiEstimates));
}







#Same model 5 as in 07122012 report.
#Random walk, t-distributed residuals
model5=function(data,nu,ppc=10){
  #ppc is the number of usable data points per country.
  #The first data point doesn't count as usable since we condition on it.
  #ppc should be 10 if we're excluding only the most recent data point.
  #ppc should be 8 to exclude the 3 most recent data points.
  
  n=nrow(data);
  logLikVector=rep(0,n);
  sigma2Vector=rep(0,n);
  
  for(i in 1:n){
    x=as.numeric(data[i,2:(2+ppc-1)]);
    y=as.numeric(data[i,3:(3+ppc-1)]);
    r=y-x;
    
    logLikOfSigma2=function(sigma2){
      result=ppc*lgamma((nu+1)/2);
      result=result-ppc*lgamma(nu/2);
      result=result-ppc/2*log(pi*nu*sigma2);
      result=result-((nu+1)/2)*sum(log(1+r^2/(nu*sigma2)));
      return(result);
    }
    model=optimize(logLikOfSigma2,
                   interval=c(1e-20,10^50),
                   maximum=TRUE,
                   tol=.Machine$double.eps^0.75);
    
    logLikVector[i]=model$objective;
    sigma2Vector[i]=model$maximum;
  }
  
  
  
  #Write parameters to file
  paramEstimates=rep(0,3*n+1)
  paramEstimates[1:n]=rep(0,n);
  paramEstimates[n+1]=nu;
  paramEstimates[(n+2):(2*n+1)]=rep(1,n);
  paramEstimates[(2*n+2):(3*n+1)]=sigma2Vector^(-1);
  write(paramEstimates,"junk.txt")#Order of parameters: mus, nu, phis, and taus
  
  
  return(list(nu=nu,
              sigma2Vector=sigma2Vector));
}

mod5SensitivityAnalysis=function(n){
  rateData=read.csv("migrationRateData.csv")
  rateData=rateData[which(rateData$Country.code<900),] #Get rid of regions
  r=dim(rateData)[1];
  c=dim(rateData)[2];      
  
  #Allocate storage space for variables we're interested in tracking across runs.
  #The entry in row i, column j gives the results for country i's parameter in trial j.
  sigma2Estimates=matrix(0,nrow=r,ncol=n);
  
  for(i in 1:n){
    cat(i,"\n");
    #Apply some jitter to the rate data.
    jitteredRateData=rateData+runif(r*c,min=-0.0005,max=0.0005);
    
    #Fit model 5 and hold onto the variables.
    mod=model5(jitteredRateData,2);
    sigma2Estimates[,i]=mod$sigma2;
  }
  return(sigma2Estimates);
}