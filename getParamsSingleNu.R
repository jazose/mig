#This file is set up so that you can just execute all the code once to read in
#functions and then make a call to main.


#Required libraries
#library(quadprog)
#library(gplots)
#library(hett)
#library(arm)
#library(car)
#library(plotrix)
library(rjags)
library(coda)



main=function(){
  #setwd("C:/Users/jonazose/Dropbox/RA/Code/stateOfTheArt")
  setwd("C:/Users/jon-everyday/Dropbox/RA/Code/faloorum")
  
  setup();
  
  hierModelSingleNu("rateModel6.bug.R",rateData,ppc=10,"Hier6_1Exc");
}


hierModelSingleNu=function(bugFile,data,ppc=10,outFile){
  #Sample call: hierModelSingleNu("rateModel6.bug.R",rateData,ppc=10,"Hier6_1Exc")
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
  nParams=r*3+1;
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
