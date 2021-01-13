library(raster);library(reshape2);library(ggplot2)
library(mblm);library(zoo);library(rNOMADS)
library(RNetCDF);library(ncdf4);library(yhat)
library(grid);library(gridExtra);library(ppcor)
library(Kendall);library(pals);library(magic)
library(car);library(ks);library(Hmisc)
library(rgeos);library(sf);library(foreign)
library(tidyselect);library(R.utils);library(data.table)
library(randomForest);library(caret);library(pdp)
library(rfUtilities);library("devtools");library(basicPlotteR)
library(ggplot2);library(reshape2);library(ggridges)

setwd("F:/ARCTIC_GREENING/")

##############################################################################
# Chapter II. Time heterogeneity of greening trends                          #
# Question : Did the greening of the arctic has varied dynamically in time ? #
# If so, how are arctic greening trends distributed over time ?              #
##############################################################################
# a. Load data table ----
Data = read.csv("BIS/DATA/SamplingSites/Data_table.csv",header=T,sep=",",dec=".")
Data.ndvi = dplyr::select(as.data.frame(Data),contains("NDVI_"))
Data.swi = dplyr::select(as.data.frame(Data),contains("SWI_"))

# b. Prepare median year of window slope ----
DFNUM = data.frame(YEAR = 1984:2019,NUM = 1:length(1984:2019))
x=1
y=11
TODO = c()
while(y != c(DFNUM$NUM[length(DFNUM$NUM)] + 1)){
  TODO[x] = mean(DFNUM$YEAR[x:y])
  x = x+1
  y = y+1
}

# d. Compute moving windows for NDVI ----
Data.ndvi.slope = as.data.frame(matrix(nrow = nrow(Data.ndvi), ncol=length(TODO)+1));colnames(Data.ndvi.slope) = c("ID",TODO);Data.ndvi.slope$ID = c(1:length(Data.ndvi.slope$ID))
Data.ndvi.mean = as.data.frame(matrix(nrow = nrow(Data.ndvi), ncol=length(TODO)+1));colnames(Data.ndvi.mean) = c("ID",TODO);Data.ndvi.mean$ID = c(1:length(Data.ndvi.mean$ID))
Data.ndvi.pval = as.data.frame(matrix(nrow = nrow(Data.ndvi), ncol=length(TODO)+1));colnames(Data.ndvi.pval) = c("ID",TODO);Data.ndvi.pval$ID = c(1:length(Data.ndvi.pval$ID))
Data.ndvi.num = as.data.frame(matrix(nrow = nrow(Data.ndvi), ncol=length(TODO)+1));colnames(Data.ndvi.num) = c("ID",TODO);Data.ndvi.num$ID = c(1:length(Data.ndvi.num$ID))

for(i in 1:nrow(Data.ndvi)){
  
  x = 1
  y = 11
  z = 2
  
  while(y != 37){
    
    TODO = DFNUM$YEAR[x:y] 
    
    dat.ndvi = Data.ndvi[i,]
    dat.ndvi = data.frame(YEAR = as.numeric(substr(colnames(dat.ndvi),6,9)), NDVI = as.numeric(dat.ndvi))
    dat.ndvi = as.data.frame(na.omit(dat.ndvi))
    dat.ndvi = dat.ndvi[which(dat.ndvi$YEAR %in% TODO),]
    
    if(length(dat.ndvi$YEAR) > 8){
      MOD = mblm(NDVI ~ YEAR,dataframe = dat.ndvi)
      
      Data.ndvi.slope[i,z] = MOD$coefficients[[2]]
      Data.ndvi.mean[i,z] = mean(dat.ndvi$NDVI)
      Data.ndvi.num[i,z] = mean(dat.ndvi$YEAR)
      
      dat.ts = as.ts(read.zoo(dat.ndvi))
      Data.ndvi.pval[i,z] = MannKendall(dat.ts)[[2]][1]
      
    } else {Data.ndvi.slope[i,z] = NA; Data.ndvi.mean[i,z] = NA; Data.ndvi.num[i,z] = NA; Data.ndvi.pval[i,z] = NA}
    
    x = x + 1
    y = y + 1
    z = z + 1
    
  }
  print(paste0("ID = ",Data.ndvi.slope$ID[i]))
}

write.table(Data.ndvi.slope,"BIS/DATA/Resultats/Greening/Data_ndvi_slope.csv",row.names = T,col.names = T,sep=";")
write.table(Data.ndvi.pval,"BIS/DATA/Resultats/Greening/Data_ndvi_pval.csv",row.names = T,col.names = T,sep=";")
write.table(Data.ndvi.num,"BIS/DATA/Resultats/Greening/Data_ndvi_num.csv",row.names = T,col.names = T,sep=";")
write.table(Data.ndvi.mean,"BIS/DATA/Resultats/Greening/Data_ndvi_mean.csv",row.names = T,col.names = T,sep=";")

# e. Compute greening trends for entire time series ----
Data.ndvi.ts.all = data.frame(ID = 1:nrow(Data.ndvi), SLP = NA, PVAL = NA, MEAN = NA)

for(i in 1:nrow(Data.ndvi)){
  
  print(i)
  
  dat.ndvi = Data.ndvi[i,]
  dat.ndvi = data.frame(YEAR = as.numeric(substr(colnames(dat.ndvi),6,9)), NDVI = as.numeric(dat.ndvi))
  dat.ndvi = as.data.frame(na.omit(dat.ndvi))
  
  if(length(dat.ndvi$YEAR) > 20){
    MOD = mblm(NDVI ~ YEAR,dataframe = dat.ndvi)
    
    Data.ndvi.ts.all[i,2] = MOD$coefficients[[2]]
    Data.ndvi.ts.all[i,4] = mean(dat.ndvi$NDVI)
    
    dat.ts = as.ts(read.zoo(dat.ndvi))
    Data.ndvi.ts.all[i,3] = MannKendall(dat.ts)[[2]][1]
    
  } else {Data.ndvi.ts.all[i,2] = NA; Data.ndvi.ts.all[i,3] = NA; Data.ndvi.ts.all[i,4] = NA}
  
}

write.table(Data.ndvi.ts.all,"BIS/DATA/Resultats/Greening/Data_ndvi_ts_all.csv",row.names = T,col.names = T,sep=";")


# e. Figure 2 ----
COLNAMES = c("ID"   ,  "1989", "1990", "1991", "1992" ,"1993" ,"1994", "1995", "1996", "1997", "1998" ,"1999" ,"2000", "2001" ,
             "2002", "2003", "2004" ,"2005","2006", "2007" ,"2008", "2009" ,"2010" ,"2011" ,"2012", "2013","2014")
Data.ndvi.slope = read.csv("BIS/DATA/Resultats/Greening/Data_ndvi_slope.csv",header = T,sep=" "); colnames(Data.ndvi.slope) = COLNAMES
Data.ndvi.pval = read.csv("BIS/DATA/Resultats/Greening/Data_ndvi_pval.csv",header = T,sep=" "); colnames(Data.ndvi.pval) = COLNAMES

# Plot all pixels and mean
n = sample(1:length(Data.ndvi.slope$ID),length(Data.ndvi.slope$ID)/10 )
windows(width = 2000,height = 2000)
par(mfrow=c(2,1))
par(mar=c(0,5,5,5))

# Base plot with points and line
plot(1, ylim=c(-0.025,0.045), xlim=c(1989.5,2013.5),ylab="Absolute annual increases in NDVI",xlab="",xaxt="n")
for(i in 1:length(n)){
  
  dat.ndvi = Data.ndvi.slope[ n[i] , -1]
  dat.pval = Data.ndvi.pval[ n[i], -1]
  dat.ndvi = data.frame(YEAR = as.numeric(colnames(dat.ndvi)), NDVI = as.numeric(dat.ndvi))
  dat.pval = data.frame(YEAR = as.numeric(colnames(dat.pval)), PVAL = as.numeric(dat.pval))
  
  dat.ndvi.greening005 = dat.ndvi; dat.ndvi.greening005$NDVI[which(dat.pval$PVAL > 0.05 | dat.ndvi$NDVI < 0)] = NA
  dat.ndvi.browning005 = dat.ndvi; dat.ndvi.browning005$NDVI[which(dat.pval$PVAL > 0.05 | dat.ndvi$NDVI > 0)] = NA
  dat.ndvi.unsign = dat.ndvi; dat.ndvi.unsign$NDVI[which(dat.pval$PVAL < 0.05)] = NA
  
  points(dat.ndvi.greening005 ,col = adjustcolor("darkgreen",alpha.f = 0.05),type="o")
  points(dat.ndvi.browning005 ,col = adjustcolor("darkred",alpha.f = 0.05),type="o")
  points(dat.ndvi.unsign ,col = adjustcolor("grey",alpha.f = 0.05),type="o")
  
}
axis(3, at = seq(1985,2015,3), labels = seq(1985,2015,3))
abline(h=0,lty=1)
dat.ndvi = apply(Data.ndvi.slope[,-1],MARGIN = 2,FUN=function(x){mean(x,na.rm=T)})
dat.ndvi = data.frame(YEAR = as.numeric(names(dat.ndvi)), NDVI = as.numeric(dat.ndvi))
dat.error = apply(Data.ndvi.slope[,-1],MARGIN = 2,FUN=function(x){sd(x,na.rm=T)})
dat.error = data.frame(YEAR = as.numeric(names(dat.error)), ERROR = as.numeric(dat.error))
lines(dat.ndvi$YEAR,dat.ndvi$NDVI,col="black",lwd=2,lty=2)
lines(dat.ndvi$YEAR,dat.ndvi$NDVI+dat.error$ERROR,col="black",lwd=1,lty=3)
lines(dat.ndvi$YEAR,dat.ndvi$NDVI+2*dat.error$ERROR,col="black",lwd=1,lty=3)
lines(dat.ndvi$YEAR,dat.ndvi$NDVI-dat.error$ERROR,col="black",lwd=1,lty=3)
lines(dat.ndvi$YEAR,dat.ndvi$NDVI-2*dat.error$ERROR,col="black",lwd=1,lty=3)

# Add boxplot of significant slopes
Data.ndvi.slopePpos005 =  Data.ndvi.slope
Data.ndvi.slopePpos001 =  Data.ndvi.slope
Data.ndvi.slopePneg001 =  Data.ndvi.slope
Data.ndvi.slopePneg005 =  Data.ndvi.slope
Data.ndvi.slopePno =  Data.ndvi.slope
for(i in 2:ncol(Data.ndvi.slope)){
  Data.ndvi.slopePpos005[which(Data.ndvi.pval[,i] > 0.05 | Data.ndvi.pval[,i] < 0.01 | Data.ndvi.slope[,i] < 0),i] = NA
  Data.ndvi.slopePpos001[which(Data.ndvi.pval[,i] > 0.01 | Data.ndvi.slope[,i] < 0),i] = NA
  Data.ndvi.slopePneg005[which(Data.ndvi.pval[,i] > 0.05 | Data.ndvi.pval[,i] < 0.01 | Data.ndvi.slope[,i] > 0),i] = NA
  Data.ndvi.slopePneg001[which(Data.ndvi.pval[,i] > 0.01 | Data.ndvi.slope[,i] > 0),i] = NA
  Data.ndvi.slopePno[which(Data.ndvi.pval[,i] < 0.05),i] = NA
}
par(mar=c(10,5,0,5))
DDD = t(data.frame(
  as.numeric(colSums(!is.na(Data.ndvi.slopePneg001[,-1]))),
  as.numeric(colSums(!is.na(Data.ndvi.slopePneg005[,-1]))),
  as.numeric(colSums(!is.na(Data.ndvi.slopePno[,-1]))),
  as.numeric(colSums(!is.na(Data.ndvi.slopePpos005[,-1]))),
  as.numeric(colSums(!is.na(Data.ndvi.slopePpos001[,-1])))))
rownames(DDD) = c("Neg001","Neg005","Nosign","Pos005","Pos001")
colnames(DDD) = COLNAMES[-1]
DDD = apply(DDD,MARGIN = 2,FUN=function(x){x/length(Data.ndvi.slope$ID)}); DDD = DDD * 100
barplot(DDD,col=c("darkred","lightcoral","grey","palegreen","darkgreen"),ylim=c(0,100),
        density=c(500,500,10,500,500) , angle=c(11,11,11,11,11),xaxs="i",
        space=0,width=2,ylab="% of pixels",las=2)
box(which="plot", bty="]")
