#This file is set up so that you can just execute all the code once to read in
#functions and then make a call to main.
#Example call to main:
#CImain(paramFile="./PE/PE_Hier73_1Exc",nChains=4,nExcluded=1,CIWidths=c(0.8,0.95),outFile="./CI/CI_Hier73_1Exc.txt")

#Define variables for the call to main:
paramFile="./PE/PE_Hier6_1Exc";nChains=4;nExcluded=1;CIWidths=c(0.8,0.95);
outFile="./CI/CI_Hier6_1Exc.txt"

#Set working directory
#setwd("C:/Users/jonazose/Dropbox/RA/Code/faloorum")
setwd("C:/Users/jon-everyday/Dropbox/RA/Code/faloorum")



#Load required libraries
if(!exists("makeNetRates", mode="function")){source("ageSexPredict_MCMC.R")}


CImain=function(paramFile,nChains,nExcluded,CIWidths,outFile){
  setup();
  CIArrays=getCIs(paramFile,nChains,nExcluded,CIWidths,outFile);
}

getCIs=function(paramFile,nChains,nExcluded,CIWidths=0.8,outFile){
  #Sample call: getCIs(paramFile="./PE/PE_Hier73_1Exc",nChains=4,nExcluded=1,CIWidths=c(0.8,0.95),outFile="./CI/CI_Hier73_1Exc.txt")
  #Most inputs are straightforward.
  #CIWidths may be a vector so that we can get multiple confidence intervals at once
  #For example, to get 0.8 and 0.95 CI's, take CIWidths=c(0.8,0.95).
  
  ppc=11-nExcluded #ppc gives the number of data points use per country.
  
  #Figure out how many parameter files are being read in.
  nFiles=length(paramFile)*nChains
  
  #Initialize result list.
  result=list();
  
  #Figure out dimensions of the parameter matrix
  f=paste(paramFile,"_Chain1.txt",sep="");

  #Read parameter values in from file.
  paramEstimates=scan(f);
  nParams=paramEstimates[1];
  nSimPerChain=paramEstimates[2];
  nSimulations=nSimPerChain*nChains;
  paramMatrix=matrix(0,nrow=nParams,ncol=nSimulations);
  paramMatrix[,1:nSimPerChain]=paramEstimates[3:(2+nSimPerChain*nParams)]
  for(i in 2:nChains){
    f=paste(paramFile,"_Chain",i,".txt",sep="");
    paramEstimates=scan(f);    
    paramMatrix[,(nSimPerChain*(i-1)+1):(nSimPerChain*i)]=paramEstimates[3:(2+nSimPerChain*nParams)]
  }
  
  n=(nParams-1)/3;
  mu=t(paramMatrix[1:n,]);
  nu=t(paramMatrix[(n+1),]);
  nuMatrix=matrix(0,nrow=nrow(mu),ncol=ncol(mu))
  for(i in 1:ncol(nuMatrix)){
    nuMatrix[,i]=nu
  }
  phi=t(paramMatrix[(n+2):(2*n+1),]);
  tau=t(paramMatrix[(2*n+2):(3*n+1),]);
  sigma=tau^(-1/2);
  
  #Make a vector of quantiles we'll want corresponding to the CI Widths.
  nCIs=length(CIWidths)
  quantileVec=rep(0,2*nCIs+1);
  CIWidths=sort(CIWidths,decreasing=TRUE);
  quantileVec[nCIs+1]=0.5;
  for(i in 1:length(CIWidths)){
    thisCIWidth=CIWidths[i];
    quantileVec[i]=(0.5+thisCIWidth/2);
    quantileVec[2*nCIs+2-i]=(0.5-thisCIWidth/2);
  }
  
  #Initialize storage for result CIs
  CIArray=array(0,dim=c(2*nCIs+1,nExcluded+18,n))
  #Each of the n layers in an array represents CI data for a single country.
  #If we're only looking for one confidence interval,
  #the 3 rows give the upper, median, and lower values for the CI.
  #For two CIs, we've got upper1, upper2, median, lower2, lower1.
  #The nExcluded+18 columns represent chronologically increasing time points.    
  
  #Start constructing an input to the predict function.
  z=list();
  
  z$rates=matrix(0,nrow=nSimulations,ncol=n);
  for(i in 1:nSimulations){
    z$rates[i,]=rateData[,(13-nExcluded)];
  }
  z$mu=mu;
  z$phi=phi;
  z$sigma=sigma;
  z$nu=nuMatrix;
  z$maleMatrix=maleMatrix;
  z$femaleMatrix=femaleMatrix;
  z$pop=futurePopData[,(5-nExcluded)];
  
  newRates=predict_MCMC(z)$rates;
  
  #Recover and store CIs from the prediction.
  for(i in 1:n){
    CIArray[,1,i]=as.vector(quantile(newRates[,i],probs=quantileVec))
  }
  
  cat("Prediction 1 Complete.\n");
  
  #Repeat into the future until population projections run out.
  for(t in 2:(18+nExcluded)){#18 points if 1 excluded, 20 points if 3 excluded
    z$rates=newRates;
    z$pop=futurePopData[,(4-nExcluded+t)];
    newRates=predict_MCMC(z)$rates;
    
    #Recover and store CIs from the prediction.
    for(i in 1:n){
      CIArray[,t,i]=as.vector(quantile(newRates[,i],probs=quantileVec))
    }
    
    cat("Prediction",t,"Complete.\n");
  }
  
  write(c(dim(CIArray),CIArray),outFile);
  #Write the array of confidence intervals into storage.
  #The first three elements should be the dimensions of the array so that
  #we can recover it precisely.
  
  result=c(result,list(CIArray));
  
  cat("\n",paramFile,"\n");
  if(nExcluded>0){
    evaluate(result[[1]],nExcluded);      
  }
  
  return(result);
}


####################################
#LESS USED FUNCTIONS BELOW THIS POINT
####################################

evaluate=function(ci,nExcluded){
  #Take in the set of confidence intervals,
  #scale parameters for all countries, and number of time points excluded

  #Recover the number of confidence intervals.
  #The first dimension of any element of CIArrays must be 2*nCIs+1
  nCIs=(dim(ci)[1]-1)/2;
  
  #Grab the median predictions for the points for which
  #the true value is known.
  predictions=as.vector(ci[nCIs+1,1,])
  if(nExcluded>1){
    for(i in 2:nExcluded){
      predictions=c(predictions,
                    as.vector(ci[nCIs+1,i,]))
    }
  }
  
  #Get actual values for all those points.
  truePoints=as.vector(rateData[,(13-nExcluded+1)]);#13 for nExcluded=1, 11 for nExcluded=3
  if(nExcluded>1){
    for(i in 2:nExcluded){
      truePoints=c(truePoints,
                   as.vector(rateData[,(13-nExcluded+i)]))
    }
  }
  
  #Compute some metrics
  errorVector=predictions-truePoints;
  
  squaredErrorVector=errorVector^2;
  rmsError=sqrt(mean(squaredErrorVector));
  cat("rmsError:", rmsError,"\n");
  
  meanAbsError=mean(abs(errorVector));
  cat("MAE:", meanAbsError,"\n");  
  
  #Compute CI inclusions (possibly more than one).
  #Also compute interval scores (assuming 80% and 95% CIs)
  for(j in 1:nCIs){
    #Similarly, get upper and lower bounds for the confidence intervals
    upperCI=as.vector(ci[nCIs+1-j,1,])
    if(nExcluded>1){
      for(i in 2:nExcluded){
        upperCI=c(upperCI,
                  as.vector(ci[nCIs+1-j,i,]))
      }
    }
    
    lowerCI=as.vector(ci[nCIs+1+j,1,])
    if(nExcluded>1){
      for(i in 2:nExcluded){
        lowerCI=c(lowerCI,
                  as.vector(ci[nCIs+1+j,i,]))
      }
    }
    ciInclusion=mean(truePoints>lowerCI & truePoints<upperCI);
    cat("CI inclusion",j,":",ciInclusion,"\n\n");
    
    #Assume we had 80% and 95% confidence intervals and calculate interval scores
    alpha=ifelse(j==1,0.2,0.05);
    score=sum(upperCI-lowerCI);
    for(k in 1:length(truePoints)){
      if(truePoints[k]>upperCI[k]){score=score+2/alpha*(truePoints[k]-upperCI[k]);}
      if(truePoints[k]<lowerCI[k]){score=score+2/alpha*(lowerCI[k]-truePoints[k]);}
    }
    #Convert to a mean interval score.
    score=score/length(truePoints);
    cat("Interval score",j,":",score,"\n\n");
    
  }

}

nameLookup=function(countryNames,countryCode){
  return(toString(countryNames[which(countryNames$V2==countryCode),1]))
}

setup=function(){
  
  if(!exists("countData")){  
    countData<<-read.csv("migrationData.csv")
    countData<<-countData[which(countData$Country.code<900),] #Get rid of regions
    
    #Important indices: 54 (Sierra Leone), 59 (North Korea), 61 (Mongolia), 176 (Ecuador)
    
    r=dim(countData)[1];
    c=dim(countData)[2];
    
    rateData<<-read.csv("migrationRateData.csv")
    rateData<<-rateData[which(rateData$Country.code<900),] #Get rid of regions
    
    #Allow for easy lookup of country names from codes
    countryNames<<-read.csv("countryNames.csv",header=FALSE, sep=";")
    countryNames<<-countryNames[which(countryNames$V2<900),] #Get rid of regions
    
    rawPopData=read.csv("popData.csv",sep=";")
    rawPopData=rawPopData[which(rawPopData$Country.code<900),] #Get rid of regions
    
    #Convert population data to five-year averages
    popData=rateData
    n=nrow(popData)
    nCountries=n;
    
    #Extend popData by adding needed columns
    popData$"2010.2015"=rep(0,n)
    popData$"2015.2020"=rep(0,n)
    popData$"2020.2025"=rep(0,n)
    popData$"2025.2030"=rep(0,n)
    popData$"2030.2035"=rep(0,n)
    popData$"2035.2040"=rep(0,n)
    popData$"2040.2045"=rep(0,n)
    popData$"2045.2050"=rep(0,n)
    popData$"2050.2055"=rep(0,n)
    popData$"2055.2060"=rep(0,n)
    popData$"2060.2065"=rep(0,n)
    popData$"2065.2070"=rep(0,n)
    popData$"2070.2075"=rep(0,n)
    popData$"2075.2080"=rep(0,n)
    popData$"2080.2085"=rep(0,n)
    popData$"2085.2090"=rep(0,n)
    popData$"2090.2095"=rep(0,n)
    popData$"2095.2100"=rep(0,n)
    
    for(i in 1:nCountries){
      index=which(rawPopData$Country.code==popData$Country.code[i])
      thisCountryData=rawPopData[index,];
      for(j in 1:(ncol(popData)-1)){
        popData[i,(j+1)]=mean(as.matrix(thisCountryData[(5*j-3):(5*j+1)]))
      }
    }
    
    #Grab the only part of the population data we're likely to need
    futurePopData<<-popData[,-(2:10)]#Includes 1995 onward
  }
  
  
  #Construct age/sex schedules on the basis of MigByAgeAndSex.csv, if we haven't already
  if(!exists("maleMatrix")){
    nAgeGroups=21;
    #Initialize male and female matrices.
    maleMatrix<<-array(0,dim=c(nAgeGroups,nCountries));
    femaleMatrix<<-array(0,dim=c(nAgeGroups,nCountries));
    
    dat=read.csv("MigByAgeAndSex.csv");
    dat=dat[dat$Revision==2010,];#Grab only the most recent revision
    dat=dat[dat$LocID<900,]#Remove regions
    dat$AgeSortOrder=ifelse(dat$AgeSortOrder>20,21,dat$AgeSortOrder);#Rename age
    #category 1000 to 21, where it really belongs.
    dat=dat[,-c(1,2,3,5,8)]#Toss out columns we don't need.
    for(j in 1:nCountries){
      countryCode=rateData[j,1];
      relevantData=dat[dat$LocID==countryCode&dat$Year=="2005-2010",]
      netMig=sum(relevantData$Total)
      if(netMig==0){#This is an annoyance. Some countries have 0 net migration in 2005-2010,
        #So we get around by picking a different time period for them.
        relevantData=dat[dat$LocID==countryCode&dat$Year=="2000-2005",];
        netMig=sum(relevantData$Total);
      }
      if(netMig==0){#This is an annoyance. Some countries have 0 net migration in 2005-2010,
        #So we get around by picking a different time period for them.
        relevantData=dat[dat$LocID==countryCode&dat$Year=="1995-2000",];
        netMig=sum(relevantData$Total);
      }
      if(netMig==0){#This is an annoyance. Some countries have 0 net migration in 2005-2010,
        #So we get around by picking a different time period for them.
        relevantData=dat[dat$LocID==countryCode&dat$Year=="1990-1995",];
        netMig=sum(relevantData$Total);
      }
      if(netMig==0){#This is an annoyance. Some countries have 0 net migration in 2005-2010,
        #So we get around by picking a different time period for them.
        relevantData=dat[dat$LocID==countryCode&dat$Year=="1985-1990",];
        netMig=sum(relevantData$Total);
      }
      if(netMig==0){#This is an annoyance. Some countries have 0 net migration in 2005-2010,
        #So we get around by picking a different time period for them.
        relevantData=dat[dat$LocID==countryCode&dat$Year=="1980-1985",];
        netMig=sum(relevantData$Total);
      }
      if(netMig==0){#This is an annoyance. Some countries have 0 net migration in 2005-2010,
        #So we get around by picking a different time period for them.
        relevantData=dat[dat$LocID==countryCode&dat$Year=="1975-1980",];
        netMig=sum(relevantData$Total);
      }
      if(netMig==0){cat("ISSUE WITH COUNTRY CODE:",countryCode,"\n")}
      for(i in 1:nAgeGroups){
        maleMatrix[i,j]<<-relevantData$Male[relevantData$AgeSortOrder==i]/netMig;
        femaleMatrix[i,j]<<-relevantData$Female[relevantData$AgeSortOrder==i]/netMig;        
      }
    }
  }
}

getCIFromFile=function(fileName,nExcluded){
  #Files are assumed to be located in the ./CIs/ directory
  #Names all start with "CI_", so that should be excluded.
  #For example, to read in ./CIs/CI_Hier6_1Exc
  #Call getCIFromFile("Hier6_1Exc",1)
  
  #The first three elements of each CI file give array dimensions. The rest are the
  #array contents.
  
  fullFileName=paste("./CIs/CI_",fileName,".txt",sep="");
  ciData=scan(fullFileName);
  dimensions=ciData[1:3];
  ciArray=array(ciData[4:length(ciData)],dim=dimensions);
  return(ciArray);
}

evaluateCIFromFile=function(fileName,nExcluded){
  #Call evaluateCIFromFile("Hier8_1Exc",1)
  ciArray=getCIFromFile(fileName,nExcluded);
  evaluate(ciArray,nExcluded);  
}

CImain(paramFile,nChains,nExcluded,CIWidths,outFile)