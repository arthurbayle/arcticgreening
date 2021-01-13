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

#################################################
# Chapter I. Sampling sites and NDVI extraction #
#################################################
# a. Land cover data downloaded from https://www.foretouverte.gouv.qc.ca/ ----
# Données écoforestières -> Couches écoforestières du Nord québecois -> Végétation du Nord québecois
# Tiles : 13L, 13M, 14D, 14E, 14L, 23I, 23J, 23O, 23P, 24A, 24B, 24G, 24H, 24I, 24J, 24P, 25A

# Centroids of the Land cover product were computed on ArcGIS Pro 2.5.1. (Veg_centroids.shp; n = 109447)
# Samples with distance from the nearest samples < 300 m were removed to reduce auto-correlation (Veg_centroids_dist.shp; n = 86339)
cent.veg = shapefile("BIS/DATA/Landcover/Centroids/Veg_centroids_dist.shp")

# Recode column name
names(cent.veg) = "LC"

# b. Land cover initial coding to new code ----
# Prostrate tundra shrubs with rocks : 3 -> PSTR
# Prostrate tundra shrubs : 4 -> PST
# Erected tundra shrubs with rocks : 2 -> ESTR
# Erected tundra shrubs : 1 -> EST
# Subarctic vegetation : 5 -> SAV
# Subarctic vegetation dominated by shrubs : 6 -> SAV
# Coniferous forest with low deciduous density : 7 -> BF
# Coniferous forest with low density mosses and lichens : 8 -> BF
# Coniferous forest with middle density mosses : 9 -> BF
# Coniferous forest with low density mosses : 10 -> BF
cent.veg$LC[which(cent.veg$LC == 5 | cent.veg$LC == 6)] = "SAV"
cent.veg$LC[which(cent.veg$LC == 7 | cent.veg$LC == 8 | cent.veg$LC == 9 | cent.veg$LC == 10)] = "BF"
cent.veg$LC[which(cent.veg$LC == 1)] = "EST"
cent.veg$LC[which(cent.veg$LC == 2)] = "ESTR"
cent.veg$LC[which(cent.veg$LC == 3)] = "PSTR"
cent.veg$LC[which(cent.veg$LC == 4)] = "PST"

# Number of samples per land cover
table(cent.veg$LC);barplot(table(cent.veg$LC))

# c. Load distance table before removing nearest sample (computed on ArcGIS Pro 2.5.1 : "Generate proximity table") ----
cent.dist.before = read.csv(file = "BIS/DATA/Landcover/Dist/Veg_centroids_table.csv") # Mean distance : 485 m
cent.dist.after = read.csv(file = "BIS/DATA/Landcover/Dist/Veg_centroids_dist_table.csv") # Mean distance : 574 m

# Plot change in distance from nearest sample
hist(cent.dist.before$NEAR_DIST,breaks=200,xlim=c(0,2000),main="",xlab="Distance from nearest sample",col=adjustcolor("black",alpha.f = 0.3))
par(new=T)
hist(cent.dist.after$NEAR_DIST,breaks=200,xlim=c(0,2000),col=adjustcolor("royalblue",alpha.f = 0.4),xaxt="n",yaxt="n",xlab="",ylab="",main="")
abline(v=300,lty=2,lwd=2)

# d. Load 37 years of NDVImax ----
ndvi.sta = stack(list.files("BIS/DATA/NDVImax/",full.names=T))

# Reproject to long lat WGS
cent.veg = spTransform(cent.veg,CRSobj = crs(ndvi.sta))

# Extract NDVI on samples
for(i in 1:nlayers(ndvi.sta)){print(i);cent.veg$temp = extract(ndvi.sta[[i]],cent.veg)
names(cent.veg)[length(names(cent.veg))] = paste0("NDVI_",substr(names(ndvi.sta[[i]]),9,12))}
shapefile(cent.veg,"BIS/DATA/Landcover/Centroids/Veg_centroids_dist_ndvi.shp")

# e. Compute number of clear observation for each sample ----
cent.veg = shapefile("BIS/DATA/Landcover/Centroids/Veg_centroids_dist_ndvi.shp")
temp = data.frame(clean = 1:length(cent.veg))
for(i in 1:length(cent.veg)){print(i);temp$clean[i] = length(na.omit(as.numeric(as.data.frame(cent.veg[i,])[-c(1,39,40)])))}

# Remove samples with less than 20 years of observations (n = 85830)
cent.veg = cent.veg[-which(temp$clean < 20),]

# Compute mean NDVI for each sample
temp = data.frame(NDVImean = 1:length(cent.veg))
for(i in 1:length(cent.veg)){print(i);temp$NDVImean[i] = mean(na.omit(as.numeric(as.data.frame(cent.veg[i,])[-c(1,39,40)])))}
hist(temp$NDVImean)

# Remove samples with a NDVI mean inferior to 0.05 (n = 85703)
cent.veg = cent.veg[-which(temp$NDVImean < 0.1),]

# f. Write final samples ----
shapefile(cent.veg,"BIS/DATA/Landcover/Centroids/All_sampling_sites.shp")

# g. Compute climate variables from https://climate.northwestknowledge.net/TERRACLIMATE-DATA/ ----
cent.veg = shapefile("BIS/DATA/Landcover/Centroids/All_sampling_sites.shp")

# Summer Warmth Index (SWI)
for(i in YEAR){
  print(i)
  
  TCtmax = stack(paste0("F:/ARCTIC_GREENING/DATA/CLIMATE/MOD4_TerraClimate/TMAX/TerraClimate_tmax_",i,".nc"),varname="tmax")
  TCtmin = stack(paste0("F:/ARCTIC_GREENING/DATA/CLIMATE/MOD4_TerraClimate/TMIN/TerraClimate_tmin_",i,".nc"),varname="tmin")
  
  cent.veg = spTransform(cent.veg,crs(TC))
  TCtmax = crop(TCtmax,cent.veg)
  TCtmin = crop(TCtmin,cent.veg)
  
  TCtmax_tmp = TCtmax[[6:8]]
  TCtmin_tmp = TCtmax[[6:8]]
  
  TCtmean = list()
  for(z in 1:nlayers(TCtmax_tmp)){ 
    TCtmeanT = mean(TCtmax_tmp[[z]],TCtmin_tmp[[z]]) 
    TCtmeanT[which(TCtmeanT[] < 0 )] = NA
    TCtmean[[z]] = TCtmeanT
  }
  
  TCtmean_tmp = stack(TCtmean)
  
  SWI = calc(TCtmean_tmp,fun = sum)
  plot(SWI,main=i)
  
  writeRaster(SWI,paste0("F:/ARCTIC_GREENING/DATA/CLIMATE/MOD4_TerraClimate/INDICATEURS/SWI/SWI_",i,".tif"),overwrite=T)
}

# Summer Water Deficit (SWD)
STOCK = data.frame(YEAR = 1984:2018,SWD=NA)
x=1
for(i in YEAR){
  print(i)
  
  TCswd = stack(paste0("F:/ARCTIC_GREENING/DATA/CLIMATE/MOD4_TerraClimate/DEF/TerraClimate_def_",i,".nc"),varname="def")
  
  cent.veg = spTransform(cent.veg,crs(TCswd))
  TCswd = crop(TCswd,cent.veg)
  
  TCswd_tmp = TCswd[[6:8]]
  
  SWD = calc(TCswd_tmp,fun = sum)
  STOCK$SWD[x] = mean(SWD[],na.rm=T)
  plot(SWD,main=i)
  x=x+1
  #writeRaster(SWD,paste0("F:/ARCTIC_GREENING/DATA/CLIMATE/MOD4_TerraClimate/INDICATEURS/SWD/SWD_",i,".tif"),overwrite=T)
}

# Minimal Summer Soil Moisture (MSSM)

STOCK = data.frame(YEAR = 1984:2018,MSSM=NA)
x=1

for(i in YEAR){
  print(i)
  
  TCsm = stack(paste0("F:/ARCTIC_GREENING/DATA/CLIMATE/MOD4_TerraClimate/SM/TerraClimate_soil_",i,".nc"),varname="soil")
  
  CENT = spTransform(CENT,crs(TCsm))
  TCsm = crop(TCsm,CENT)
  
  TCsm_tmp = TCsm[[6:8]]
  
  for(z in 1:nlayers(TCsm_tmp)){
    MSSM = min(TCsm_tmp,na.rm=T)
    plot(MSSM,main=i)
  }
  STOCK$MSSM[x] = mean(MSSM[],na.rm=T)
  x=x+1
  writeRaster(MSSM,paste0("F:/ARCTIC_GREENING/DATA/CLIMATE/MOD4_TerraClimate/INDICATEURS/MSSM/MSSM_",i,".tif"),overwrite=T)
  
}
plot(STOCK,type="o")

# h. Load variables on samples ----
# Load shapefile
cent.veg = shapefile("BIS/DATA/Landcover/Centroids/All_sampling_sites.shp")

# Extract variables

# DEM
DEM = raster("BIS/DATA/Environmental/DEM_30m.tif")
cent.veg$DEM = raster::extract(DEM,cent.veg)

# NDVImean
temp = data.frame(NDVImean = 1:length(cent.veg))
for(i in 1:length(cent.veg)){print(i);temp$NDVImean[i] = mean(na.omit(as.numeric(as.data.frame(dplyr::select(as.data.frame(cent.veg),contains("NDVI"))[i,]))))}
cent.veg$NDVImean = temp$NDVImean

# SWI
SWI = stack(list.files("BIS/DATA/Climate/Indicateurs/SWI/",full.names = T))
for(i in 1:nlayers(SWI)){print(i);cent.veg$temp = raster::extract(SWI[[i]],cent.veg)
names(cent.veg)[length(names(cent.veg))] = paste0("SWI_",substr(names(SWI[[i]]),5,9))}

# SWD
SWD = stack(list.files("BIS/DATA/Climate/Indicateurs/SWD/",full.names = T))
for(i in 1:nlayers(SWD)){print(i);cent.veg$temp = raster::extract(SWD[[i]],cent.veg)
names(cent.veg)[length(names(cent.veg))] = paste0("SWD_",substr(names(SWD[[i]]),5,9))}

# MSSM
MSSM = stack(list.files("BIS/DATA/Climate/Indicateurs/MSSM/",full.names = T))
for(i in 1:nlayers(MSSM)){print(i);cent.veg$temp = raster::extract(MSSM[[i]],cent.veg)
names(cent.veg)[length(names(cent.veg))] = paste0("MSSM_",substr(names(MSSM[[i]]),6,10))}

# Order variables
SampleSite = cbind(coordinates(cent.veg)[,1],coordinates(cent.veg)[,2],cent.veg$LC,cent.veg$DEM,cent.veg$NDVImean,
                   dplyr::select(as.data.frame(cent.veg),contains("NDVI_")),
                   dplyr::select(as.data.frame(cent.veg),contains("SWI")),
                   dplyr::select(as.data.frame(cent.veg),contains("SWD")),
                   dplyr::select(as.data.frame(cent.veg),contains("MSSM")))
colnames(SampleSite) = c("Longitude","Latitude","LC","DEM","NDVImean",
                         colnames(dplyr::select(as.data.frame(cent.veg),contains("NDVI_"))),
                         colnames(dplyr::select(as.data.frame(cent.veg),contains("SWI_"))),
                         colnames(dplyr::select(as.data.frame(cent.veg),contains("SWD_"))),
                         colnames(dplyr::select(as.data.frame(cent.veg),contains("MSSM_"))))

# Remove NDVI in 2020 as TerraClimate data ends in 2019
SampleSite = SampleSite[,-42]

# i. Write final table ----
write.csv(SampleSite,"BIS/DATA/SamplingSites/Data_table.csv",row.names = F,quote = F)
