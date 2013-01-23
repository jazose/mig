setwd("C:/Users/jonazose/Dropbox/RA/Code/faloorum")
library(plotrix)

getCIFromFile=function(fileName,nExcluded){
  #Files are assumed to be located in the ./CI/ directory
  #Names all start with "CI_", so that should be excluded.
  #For example, to read in ./CI/CI_Hier6_1Exc
  #Call getCIFromFile("Hier6_1Exc",1)
  
  #The first three elements of each CI file give array dimensions. The rest are the
  #array contents.
  
  fullFileName=paste("./CI/CI_",fileName,".txt",sep="");
  ciData=scan(fullFileName);
  dimensions=ciData[1:3];
  ciArray=array(ciData[4:length(ciData)],dim=dimensions);
  return(ciArray);
}


plotCIs=function(CIArrays,nExcluded,countryCode,col="rainbow",ylim="default",plotMain="default"){
  nModels=length(CIArrays);
  #Recover the number of confidence intervals.
  #The first dimension of any element of CIArrays must be 2*nCIs+1
  nCIs=(dim(CIArrays[[1]])[1]-1)/2;
  
  colors=rainbow(nModels);#Set up graph colors
  countryIndex=which(countryNames[,2]==countryCode);
  
  #Pull data for just the one country we're looking at.
  countryData=rateData[countryIndex,];
  countryCIs=array(0,dim=c(dim(CIArrays[[1]][,,countryIndex]),nModels));
  for(i in 1:nModels){
    countryCIs[,,i]=CIArrays[[i]][,,countryIndex];
  }
  
  dataVector=as.vector(as.matrix(countryData[-1]));
  
  #If we weren't passed plot limits, figure them out now.
  if(length(ylim)!=2){
    plotMax=max(c(countryCIs,dataVector));
    plotMin=min(c(countryCIs,dataVector));
    ylim=c(plotMin,plotMax);    
  }
  
  #If we weren't passed a title, figure it out now.
  if(plotMain=="default"){
    plotMain=paste("Simulated",nameLookup(countryNames,countryCode),"Data")
  }
  
  plot(seq(1950,2005,5),
       dataVector,
       xlim=c(1950,2100),
       ylim=ylim,
       main=plotMain,
       xlab="Time",
       ylab=paste("Net Migration Rate"));
  abline(h=0,lty=3)
  
  for(i in 1:(18+nExcluded)){#19 for nExcluded=1, 21 for nExcluded=3
    for(j in 1:nModels){
      for(k in 1:nCIs){
        lineType=ifelse(k==1,1,2)
        plotCI(x=((2005-5*nExcluded)+5*i),#start at 2005 for nExcluded=1, 1995 for nExcluded=3
               y=countryCIs[nCIs+1,i,j],
               ui=countryCIs[nCIs+1-k,i,j],
               li=countryCIs[nCIs+1+k,i,j],
               col=ifelse(col=="rainbow",colors[j],col),
               pch=4,
               add=TRUE,
               slty=lineType)
      }
    }
  }
  
}


#############################
#Plot all CIs from a CI file#
#############################
visualizeCIs=function(fileName,nExcluded){
  #Call: visualizeCIs("Hier111_3Exc",3)
  #Open up a file for writing
  fullOutFile=paste("junk.pdf");
  pdf(file=fullOutFile, height=8, width=12)
  
  colorVector=c("red")
  ylim="default"
  
  par(mfrow=c(3,2));
  
  ciArray=getCIFromFile(fileName,nExcluded);
  for(i in 1:dim(rateData)[1]){
    plotCIs(list(ciArray),nExcluded,rateData[i,1],col=colorVector,ylim=ylim)
  }
  
  #Push output to file
  dev.off()
}

###############
#Make UN Plots#
###############
#Need to have all the functions from simulateNormal_MCMC.R for this to run
makeUNPlots=function(){
  fileName="Hier6_Normal_0Exc";
  nExcluded=0;
  nTrajectories=3;
  nCIs=2;
  colorVector=c("red")
  ylim="default"
  countryIndices=c(185, 56, 145, 18)#USA, China, Netherlands, Zimbabwe
  
  par(mfrow=c(2,2));
  
  ciArray=getCIFromFile(fileName,nExcluded);
  for(i in countryIndices){
    plotCIs(list(ciArray),nExcluded,rateData[i,1],col=colorVector,ylim=ylim);
    trajectoryPoints=getSomeTrajectories(paramFile=paste("./PE/PE_",fileName,sep=""),nExcluded=nExcluded,countryIndex=i,nTrajectories=nTrajectories)
    for(j in 1:nTrajectories){
      lines(x=seq(2010,2095,5),y=trajectoryPoints[j,],col="gray")
    }
    for(j in 1:(18+nExcluded)){#19 for nExcluded=1, 21 for nExcluded=3
        for(k in 1:nCIs){
          lineType=ifelse(k==1,1,2)
          plotCI(x=((2005-5*nExcluded)+5*j),#start at 2005 for nExcluded=1, 1995 for nExcluded=3
                 y=ciArray[nCIs+1,j,i],
                 ui=ciArray[nCIs+1-k,j,i],
                 li=ciArray[nCIs+1+k,j,i],
                 col="red",
                 pch=4,
                 add=TRUE,
                 slty=lineType)
        }
      }
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