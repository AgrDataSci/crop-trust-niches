library(raster)
library(rgdal)
library(dismo)
library(parallel)
library(maptools)
library(grDevices)
library(reshape2)
library(pscl)
library(geodata)

wcpath = "data/wc2.1-global"

############################################################################################
###Preparation: Define functions############################################################
############################################################################################

makeClimate <- function(path = ".", res = 5)
{
  
  rain = worldclim_global("prec", res = res, path)
  names(rain) = paste0("prec", 1:length(names(rain)))
  
  tmin = worldclim_global("tmin", res = res, path)
  names(tmin) = paste0("tmin", 1:length(names(tmin)))
  
  tmax = worldclim_global("tmax", res = res, path)
  names(tmax) = paste0("tmax", 1:length(names(tmax)))
  
  tavg = mosaic(tmin, tmax, fun = "mean")
  names(tavg) = paste0("tmean", 1:length(names(tavg)))
  
  climate1 <- raster::stack(tmin)
  climate2 <- raster::stack(tavg)
  climate3 <- raster::stack(rain)
  
  climate = raster::stack(climate1, climate2, climate3)
  
  return(climate)
  
}

makeClimate2 <- function(path=".")
{
  
  rain = worldclim_global("prec", res = res, path)
  names(rain) = paste0("prec", 1:length(names(rain)))
  
  tmin = worldclim_global("tmin", res = res, path)
  names(tmin) = paste0("tmin", 1:length(names(tmin)))
  
  tmax = worldclim_global("tmax", res = res, path)
  names(tmax) = paste0("tmax", 1:length(names(tmax)))
  
  tavg = mosaic(tmin, tmax, fun = "mean")
  names(tavg) = paste0("tmean", 1:length(names(tavg)))
  
  climate <- stack(c(tmin, tavg, rain))
  return(climate)
  
}

makeClimate3 <- function(path=".")
{
  
  listFiles <- c(paste(path, "/tmin_", 1:12, "/w001001.adf", sep=""),
                 paste(path, "/tmean_", 1:12, "/w001001.adf", sep=""),
                 paste(path, "/prec_", 1:12, "/w001001.adf", sep=""))
  stack(listFiles)
  
}

#Define ecocrop spatial function
ecs <-function(s)
{
  if(any(is.na(s))) return(0) #this speeds things up a little bit
  ecocrop(crop, s[1:12]/10, s[13:24]/10, s[25:36], rainfed=TRUE)@maxsuit
}


############################################################################################
###Preparation: get data needed in next steps###############################################
############################################################################################

# Climate data
currentClimate <- makeClimate(wcpath)

#Land area occupied by crops in 2000
cropland_5min <- raster("data/cropland/global-crop-land.tif", native=TRUE)
extent(cropland_5min) <- extent(-180, 180,-90, 90)

#Cropping the crop land raster to the same extent as the climate data
cropland_5min <- crop(cropland_5min, currentClimate)
projection(cropland_5min) <- projection(currentClimate)

#Regions
regions <- readShapeSpatial("data/world-map/world_borders_adm0")
regions <- regions[regions$REGION_WB != "Antarctica", ]
regions <- regions[regions$REGION_UN != "Seven seas (open ocean)", ]


regionnames <- as.vector(as.data.frame(regions)$REGION_UN)

# regions_20min <- raster("Data/Rasters/PlantWorldRegions_20min.asc")
# Regions_5min<-disaggregate(regions_20min, fact=4)
# uniqueRegions <- sort(unique(regions_20min))

# r <- raster(ncol=180, nrow=180)
# extent(r) <- extent(regions)
# regions_20min <- rasterize(regions, r, 'REGION_WB')
# plot(regions_20min)
# Regions_5min<-disaggregate(regions_20min, fact=4)
uniqueRegions <- sort(unique(regionnames))

#Cropping to make the rasters coincide
projection(regions) <- projection(currentClimate)
regions <- crop(regions, currentClimate)
currentClimate <- crop(currentClimate, regions)
cropland_5min <- crop(cropland_5min, regions)
cropland <- crop(cropland_5min, regions)
cropland_5min <- cropland_5min * area(cropland)
cropland <- aggregate(cropland_5min, fact=4, fun=sum)

# Current cropping area
cropCells <- which(getValues(cropland) > 0)
regions <- extract(regions, xyFromCell(cropland,cropCells))
cropCells <- cropCells[which(!is.na(regions))]
regions <- regions[which(!is.na(regions))]
cropArea <- cropland[cropCells]
suitCrop <- raster(cropland)

#Removing non existent regiong under 10 min resolution
# b<-unique(regions)
# b<-sort(b)
# regionnames<-regionnames[-c(14, 18, 22, 29, 41, 43, 44)]

# Get folder names where climate data is located
ld <- dir("/Data/Climate Data/RCP 4.5", pattern="rcp4_5_2050s", 
          full.names=TRUE, 
          include.dirs=TRUE)
ldout <- paste("/Results/Future", substr(ld, 16, 100), "/", sep="")
ldc <- c("/Results/Current", ldout)

# EcoCrop data
CropData <- read.csv("data/EcocropParametersCropType5.csv")
CropData <- CropData[CropData$Selected=="1", ]
EcoIDi <- which(!duplicated(CropData$EcoID))
EcoID <- as.integer(CropData$EcoID)
Crops <- as.character(CropData$Crop)
MonfredaFileNames <- CropData$MonfredaFileName
CropsUnique<-unique(Crops)

############################################################################################
###Step 1: Suitability under current climate conditions#####################################
############################################################################################

# climate in RAM and aggregated
currentClimate <- brick(currentClimate)
currentClimate <- readAll(currentClimate)
currentClimate <- aggregate(currentClimate, fact=2)
currentClimate <- currentClimate[cropCells]

# Make cluster
cl <- makeCluster(8)
clusterExport(cl, c("ecs", "ecocrop", "currentClimate", "movingFun"))
setwd("/Results/Current")

# Make EcoCrop crop object
crop <- getCrop("maize")

# Loop through crops
for (i in EcoIDi){
 
  EC <- CropData[i,]
  ag <- sqrt(EC$GMIN * EC$GMAX)
  crop@GMIN <- ag
  crop@GMAX <- ag
  crop@KTMP <- EC$KTMP
  crop@TMIN <- EC$TMIN
  crop@TOPMN <- EC$TOPMN
  crop@TOPMX <- EC$TOPMX
  crop@TMAX <- EC$TMAX
  crop@RMIN <- EC$RMIN
  crop@ROPMN <- EC$ROPMN
  crop@ROPMX <- EC$ROPMX
  crop@RMAX <- EC$RMAX
  clusterExport(cl, "crop")
  suit <- parApply(cl, currentClimate, 1, ecs)
  suitCrop[cropCells] <- suit
  writeRaster(suitCrop, paste("Eco", EcoID[i], ".grd", sep=""), overwrite=T)
  cat("-")
  
}

############################################################################################
###Step 2: Suitability under future climate conditions######################################
############################################################################################
setwd("/Results/Future")
# Following line needs to be run only once.

for(i in 1:length(ldout)) dir.create(ldout[i], recursive=TRUE)

# Make EcoCrop crop object
crop <- getCrop("maize")

# Loop through GCMs 
for(l in 1:length(ld))
{
  # Get CGM data
  options(tolerance=0.5) #There are small differences in the extent
  futureClimate <- makeClimate2(ld[l])
  options(tolerance=0.1) #Back to original value for safety
  futureClimate <- brick(futureClimate)
  futureClimate <- readAll(futureClimate)
  futureClimate <- aggregate(futureClimate, 2)
  futureClimate <- futureClimate[cropCells]
  
  # Make cluster and export objects
  cl <- makeCluster(8)
  clusterExport(cl, c("ecs", "ecocrop", "futureClimate", "movingFun"))
  
for (i in EcoIDi)
  {
    
    EC <- CropData[i,]
    ag <- sqrt(EC$GMIN * EC$GMAX)
    crop@GMIN <- ag
    crop@GMAX <- ag
    crop@KTMP <- EC$KTMP
    crop@TMIN <- EC$TMIN
    crop@TOPMN <- EC$TOPMN
    crop@TOPMX <- EC$TOPMX
    crop@TMAX <- EC$TMAX
    crop@RMIN <- EC$RMIN
    crop@ROPMN <- EC$ROPMN
    crop@ROPMX <- EC$ROPMX
    crop@RMAX <- EC$RMAX
    clusterExport(cl, "crop")
    suit <- parApply(cl, futureClimate, 1, ecs)
    suitCrop[cropCells] <- suit
    writeRaster(suitCrop, paste(ldout[l], "Eco", EcoID[i], ".grd", sep=""), overwrite=TRUE)
    cat(ldout[l], "Eco", EcoID[i], ".grd\n", sep="")
    
  }
  
  stopCluster(cl)
  cat(ld[l], "\n")
  
}

############################################################################################
###Step 3: Combine crop species ############################################################
############################################################################################

for(i in 1:length(CropsUnique))

{
  Cropsi <- EcoID[which(Crops %in% CropsUnique[i])]
  
 for(j in 1:length(ldc))

  {
      CropSpecies <- stack(paste0(ldc[j],"/", "Eco", Cropsi, ".grd", sep=""))
      if(length(Cropsi)>1){CropSpecies <- calc(CropSpecies, function(x) max(x, na.rm=TRUE))}
      writeRaster(CropSpecies, paste0(ldc[j], CropsUnique[i], ".grd", sep=""), overwrite=TRUE)
  }
}


############################################################################################
###Step 4: Rasters current and future suitability###########################################
############################################################################################

#Creating matrix for crop-specific area outputs to be placed into 

Suitable <- array(NA, c(length(CropsUnique), length(regionnames), length(ldc)))

for(i in 1:length(CropsUnique))

{
 
for(j in 1:length(ldc))
  {
    
    SuitCropi <- raster(paste(ldc[j], CropsUnique[i], ".grd", sep=""))
    SuitSurfCrop <- SuitCropi[cropCells] * cropArea
    SuitabilityRegion <- tapply(SuitSurfCrop, regions, function(x) sum(x, na.rm=TRUE))
    Suitable[i,,j] <- SuitabilityRegion
    
  }
}

NModelsPositiveSum <- apply(Suitable, c(1,2), function(x) sum(x[1]<x[-1]))
NModelsNegativeSum <- apply(Suitable, c(1,2), function(x) sum(x[1]>x[-1]))
NModels0916 <- (NModelsPositiveSum >= 16) - (NModelsNegativeSum >= 16)
colnames(NModels0916) <- regionnames
rownames(NModels0916) <- CropsUnique

write.csv(NModels, "/Results/Future/NModels2017.csv")



############################################################################################
###Step 5: Get citation data ###############################################################
############################################################################################

setwd("/Results/Citation/")
library(RCurl)
useragents <- c("Mozilla/5.0 (X11; U; SunOS sun4u; en-US; rv:1.7.7) Gecko/20050421",
                "Mozilla/5.0 (X11; Linux x86_64; rv:7.0.1) Gecko/20100101 Firefox/7.0.1",
                "Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.7.7) Gecko/20050427 Red Hat/1.7.7-1.1.3.4",
                "Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.7.7) Gecko/20050420 Debian/1.7.7-2",
                "Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.7.7) Gecko/20050414",
                "Mozilla/5.0 (X11; U; Linux i686; de-AT; rv:1.7.7) Gecko/20050415",
                "Mozilla/5.0 (Windows; U; Windows NT 5.1; fr-FR; rv:1.7.7) Gecko/20050414",
                "Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.7.7) Gecko/20050414",
                "Mozilla/5.0 (Windows; U; Windows NT 5.1; de-AT; rv:1.7.7) Gecko/20050414",
                "Mozilla/5.0 (Windows; U; Windows NT 5.0; en-US; rv:1.7.7) Gecko/20050414",
                "Mozilla/5.0 (Windows; U; Windows NT 5.0; de-AT; rv:1.7.7) Gecko/20050414")

#Get plant distribution regions
region <- read.csv("regiontable.csv")
code <- sort(unique(region$Region))

#Make table
ScTable <- matrix(NA, nrow=length(CropsUnique), ncol=length(code))
rownames(ScTable) <- CropsUnique
colnames(ScTable) <- code

for(i in 1:length(CropsUnique))
{
  for(j in 1:length(code))
  {
    
    useragent <- useragents[round(runif(1)*11+0.5)]
    CropsSelected <- unique(as.character(CropData$Search[which(CropData$Crop == CropsUnique[i])]))
    CropNames <- gsub(" ", "+", CropsSelected)
    CropNames <- paste("title%3A%22" , CropNames, "%22", sep="")
    CropNames <- paste(CropNames, collapse="+OR+", sep="")
    
    regionj <- region[region$Region == code[j],]
    countries <- as.character(unique(regionj$countries))
    countries <- gsub(" ", "+", countries)
    countries1 <- paste("affiliation%3A%22" , countries, "%22", sep="")
    countries2 <- paste("title%3A%22" , countries, "%22", sep="")
    countries <- paste(c(countries1, countries2), collapse="+OR+", sep="")
        
    u <- paste("http://scirus.com/srsapp/search?sort=0&t=all&q=%28", CropNames, "%29+%28", countries, "%29&cn=all&co=AND&t=all&q=&cn=all&g=a&fdt=2002&tdt=2013&dt=abs&dt=fta&dt=aip&dt=bok&dt=con&dt=pre&dt=rev&dt=uhp&ff=all&ds=jnl&ds=nom&ds=web&sa=agr", sep="")
    
    Sc <- try(getURI(u, useragent=useragent))
    if(inherits(Sc,"try-error")) 
    {
      Sc <- -999
      cat("failed\n")
    }
    else{
      Sc <- strsplit(Sc[1], "hits")[[1]]
      Sc <- strsplit(Sc[1], " of ")[[1]]
      Sc <- as.numeric(gsub(",", "", Sc[2]))
      cat(Sc, " results for ", CropNames, countries, "\n")
    }
    ScTable[i,j] <- Sc
    

  }
}

write.csv(ScTable, "ScTable.csv")


############################################################################################
############################################################################################
###Step 6: Data Preparation ################################################################
############################################################################################
############################################################################################

############################################################################################
###Step 6.1: Citation, climate, production data#############################################
############################################################################################
#Values based upon global average yields

ScTable <- read.csv("Datasets/ScTable.csv",row.names=1)
NModels2016 <- read.csv("NModels2016.csv", row.names=1)
Pro.Tot.Val<-read.csv("Regional_Production_Quantity.csv", row.names = 1)#Data developed from total production quantities FAO
Food.Perc<-read.csv("/Food.Perc.Reg.DomesSupply.All.Crop.csv", row.names = 1)
Food.Perc<-Food.Perc[order(match(rownames(Food.Perc),rownames(Pro.Tot.Val))),]
Food.Perc<-Food.Perc[,order(match(colnames(Food.Perc),colnames(Pro.Tot.Val)))]

Pro.Tot.Val<-Food.Perc*Pro.Tot.Val#Proportion of total production estimated to be available as food. Food data caculated from proportion of each crop available as food from total domestic supply


############################################################################################
###Step 6.2.1: Producer Price Data##########################################################
############################################################################################

PP<-read.csv("Producer_Prices_USD.csv")
PP<-PP[,c(4,8,9,12)]
colnames(PP)<-c("Country", "Crop", "Year", "PP")

CountReg<-read.csv("Countries and Regions.csv")

PP<-merge(PP, CountReg, by="Country")

Cr.Uni<-unique(PP$Crop)
Re.Uni<-unique(PP$Region)


PP.Regions<-matrix(NA,26, 45)
rownames(PP.Regions)<- c("broadbean", "maize", "rye","wheat",
                         "potato","bean","sorghum","oats",
                         "barley","rice",  "lentil",   "chickpea", 
                         "sweetpotato", "banana", "millet","lupin",
                         "cassava","cowpea","pigeonpea", "yam",
                         "taro","plantain",   "buckwheat",  "quinoa",
                         "bambara", "fonio")     

colnames(PP.Regions)<- c("Southeastern.Europe",           "Northern.Africa",               "Caribbean",                     "Southern.South.America",       
                         "Caucasus",                      "Australia",                     "Middle.Europe",                 "Indian.Subcontinent",          
                         "Eastern.Europe",                "Central.America",               "West.Tropical.Africa",          "Western.South.America",        
                         "Southern.Africa",               "Brazil",                        "Malesia",                       "West.Central.Tropical.Africa",
                         "Indo.China",                    "Western.Canada",                "Eastern.Canada",                "Subarctic.America",            
                         "China",                         "Western.Indian.Ocean",          "Southwestern.Pacific",          "Northern.Europe",
                         "Northeast.Tropical.Africa",     "Southwestern.Europe",           "Northern.South.America",        "Western.Asia",
                         "Eastern.Asia",                  "Middle.Asia",                   "East.Tropical.Africa",          "South.Tropical.Africa",     
                         "Mexico",                        "Mongolia",                      "Papuasia",                      "New.Zealand", 
                         "Arabian.Peninsula",             "Siberia",                       "Russian.Far.East",              "Northeastern.U.S.A.",                          
                         "Northwestern.U.S.A.",           "Southeastern.U.S.A.",           "North.Central.U.S.A.",          "South.Central.U.S.A.",                       
                         "Southwestern.U.S.A.") 

for (i in 1:length(Re.Uni))
{
  Reg<-PP[PP$Reg==Re.Uni[i],]
  
  for (j in 1:length(Cr.Uni))
  {
    CropReg<-Reg[Reg$Crop==Cr.Uni[j],]
    CropAve<-mean(CropReg$PP, na.rm=T)
    PP.Regions[j, i]<-CropAve
  } 
}



PP.Regions<-PP.Regions[order(match(rownames(PP.Regions),rownames(Pro.Tot.Val))),]
PP.Regions<-PP.Regions[,order(match(colnames(PP.Regions),colnames(Pro.Tot.Val)))]

PP.Regions<-PP.Regions[,-c(44,45)]###Removing southwestern pacific, and subarctic America

PP.Tot.Val<-as.matrix(Pro.Tot.Val*PP.Regions)

############################################################################################
###Step 6.2.2: NRF Data######################################################################
############################################################################################

#NRF9.3 Data USDA
NRF<-read.csv("NRF9.3_2018.csv", stringsAsFactors=FALSE,row.names = 1)
NRF<-subset(NRF, rownames(NRF) %in% rownames(Pro.Tot.Val))
NRF<-NRF[order(match(rownames(NRF),rownames(Pro.Tot.Val))),]
NRF93<-(NRF[,14]*20000)#Converting NRF values per 50g to values per tonne

NRF93<-(NRF93*Pro.Tot.Val)


############################################################################################
###Step 6.2.3: Calorie Data USDA###############################################################
############################################################################################
kcal<-read.csv("/Users/rhysmanners/Dropbox/Hard Drive/JERM/Results/kcal_protein_per_ha.csv", row.names=1)
kcal<-subset(kcal, rownames(kcal) %in% rownames(Pro.Tot.Val))
kcal<-kcal[order(match(rownames(kcal),rownames(Pro.Tot.Val))),]
kcalTot<-kcal[, 4]
kcalTot<-(kcalTot*Pro.Tot.Val)


############################################################################################
###Step 6.3: Aligning Data##################################################################
############################################################################################
#Citation Data
Citation <- subset(ScTable , rownames(ScTable) %in% rownames(Pro.Tot.Val))
Citation <- Citation[,order(match(colnames(Citation), colnames(Pro.Tot.Val)))]
Citation <- Citation[,-c(44:50)]
Citation<-Citation[order(match(rownames(Citation),rownames(Pro.Tot.Val))),]
Citation <- as.vector(as.matrix(Citation))
Citation[is.na(Citation)] <- 0

#Climate Winner/ Loser Data
Climate <- subset(NModels2016, rownames(NModels2016) %in% rownames(Pro.Tot.Val))
Climate <- Climate[,order(match(colnames(Climate), colnames(Pro.Tot.Val)))]
Climate <- Climate[,-c(44:50)]
Climate <- Climate[order(match(rownames(Climate),rownames(Pro.Tot.Val))),]
Climate <- as.vector(as.numeric(as.matrix(Climate)))


#NRF9.3 Data
Nutrient <- as.vector(as.numeric(as.matrix(NRF93)))

#Total Production Value
ProdVal<-as.vector(as.numeric(as.matrix(PP.Tot.Val)))


#Calorific Data
Calorific <- as.vector(as.matrix(kcalTot))

CropsUnique1<-rownames(Pro.Tot.Val)

Crop <- rep(CropsUnique1, times=dim(kcalTot)[2])
Region <- rep(colnames(Pro.Tot.Val), each=dim(kcalTot)[1])
CropGroup <- as.character(CropData$GROUP[match(Crop, CropData$Crop)])
CropGroup[CropGroup == "vegetables" | CropGroup == "fruit"] <- "FruitsVeg"
CropGroup[CropGroup == "cereals" | CropGroup == "roots and tubers" | Crop == "plantain" | Crop == "banana"] <- "CRTB"


############################################################################################
###Step 7: Statistical Analysis#############################################################
############################################################################################

CRTB <- data[data$CropGroup == "CRTB",]
pulses <- data[data$CropGroup == "pulses",]
CRTBLoser <- tapply(CRTB$Citation[CRTB$Climate == -1], CRTB$Region[CRTB$Climate == -1], function(x) sum(x, na.rm=T))
CRTBWinner <- tapply(CRTB$Citation[CRTB$Climate == 1], CRTB$Region[CRTB$Climate == 1], function(x) sum(x, na.rm=T))
CRTBLoserAcr <- tapply(CRTB$Presence[CRTB$Climate == -1], CRTB$Region[CRTB$Climate == -1], function(x) sum(x, na.rm=T))
CRTBWinnerAcr <- tapply(CRTB$Presence[CRTB$Climate == 1], CRTB$Region[CRTB$Climate == 1], function(x) sum(x, na.rm=T))
t.test(log1p(CRTBLoser / CRTBLoserAcr), log1p(CRTBWinner / CRTBWinnerAcr), paired=TRUE)
plot(as.factor(CRTB$Climate[CRTB$Presence>0.01]), CRTB$Citation[CRTB$Presence>0.01]/CRTB$Presence[CRTB$Presence>0.01])

#data <- data.0.frame(Citation, Climate, Nutrient, Crop, Region, CropGroup)
data.0 <- data.frame(Citation, Climate, Nutrient, Calorific, ProdVal, Crop, Region, CropGroup)
data.0 <- data.0[data.0$Nutrient>0 & (data.0$CropGroup == "CRTB" | data.0$CropGroup == "pulses"),]
data.0<-as.data.frame(data.0[!grepl('NA', rownames(data.0)), ])


m1 <- lm(log1p(Citation) ~ log1p(Nutrient)+log1p(Calorific)+log1p(ProdVal) + Climate + Region + CropGroup , data=data.0)
#m1 <- lm(log1p(Citation) ~ log1p(Nutrient)+ Region + log1p(Calorific) + log1p(ProdVal) + CropGroup , data.0=data.0)
pred1 <- predict(m1, newdata.0=data.0)


m2 <- lm(log1p(Citation) ~ log1p(Nutrient)+log1p(Calorific)+log1p(ProdVal)+ Climate + CropGroup , data=data.0)
#m2 <- lm(log1p(Citation) ~ log1p(Nutrient)+ log1p(Calorific) + log1p(ProdVal) + CropGroup , data.0=data.0)
pred2 <- predict(m2, newdata.0=data.0)

summary (m1)
summary (m2)

step(m1)
step(m2)

##########zero-inflated Poisson GLM
data1 <- data.0
data1$ClimatePositive <- data1$Climate >= 0
#data1$ClimateNegative <- data1$Climate == -1
#data1$ClimateNeutral <- data1$Climate == 0

# Model with zero-inflation -- predicted best by calorific output
model.zip_explanation = zeroinfl(Citation ~ scale(log1p(Nutrient)) + CropGroup + scale(log1p(ProdVal)) + scale(log1p(Calorific)) + ClimatePositive + Region | log1p(Calorific), data=data1)
summary(model.zip_explanation)

# Model without zero-inflation
model.glm_explanation <- glm(Citation ~ log1p(Nutrient) + CropGroup + log1p(ProdVal) + log1p(Calorific) + ClimatePositive + ClimateNegative + Region, data=data1, family = poisson)
step(model.zip_explanation) #all variables retained

# Test if the zero-inflated model is better
vuong(model.zip_explanation, model.glm_explanation)

# Gaps based on nutrient output only
model.zip_gaps_regional = zeroinfl(Citation ~ log1p(Nutrient) + Region | log1p(Calorific), data=data1)
summary(model.zip_gaps_regional)

model.zip_gaps_global = zeroinfl(Citation ~ log1p(Nutrient) + CropGroup | log1p(Calorific), data=data1) # Region does not influence zero inflation
summary(model.zip_gaps_global)

pred2<-predict(model.zip_gaps_global, newdata= data1, type = "count")
pred1<-predict(model.zip_gaps_regional, newdata= data1, type = "count")

plot(data1$Nutrient, data1$Citation+1, pch=20, cex=1, cex.axis=1.25, cex.lab=1.25, axes=FALSE, col=c("black","gray")[as.numeric(data1$CropGroup)], xlab="Nutrient output (total NRF9.3 of each crop-region)", ylab="Research intensity (number of publications)", log="xy")
x <- exp(seq(0, max(log1p(data1$Nutrient)), length.out=nrow(data1)))
nd1 <- data1
nd1$Nutrient <- x
nd1$CropGroup <- "CRTB"
nd2 <- nd1
nd2$CropGroup <- "pulses"
y1 <- predict(model.zip_gaps_global, newdata=nd1)
y2 <- predict(model.zip_gaps_global, newdata=nd2)
lines(cbind(x, y1), lwd=3)
lines(cbind(x, y2), col="gray", lwd=3)
axis(1, at=10^c(7,8,9,10,11,12,13,14), labels=expression(10^7, 10^8, 10^9,
                                                         10^10, 10^11, 10^12, 10^13, 10^14))
axis(2, at=c(1, 11, 101, 1001, 10001), labels=expression(0, 10, 100, 1000, 10000))


# glm1 <- glm.nb(Citation ~ Nutrient + Calorific + ProdVal+ Climate + Region + CropGroup , data=data, family=poisson)
# #m1 <- lm(log1p(Citation) ~ log1p(Nutrient)+ Region + log1p(Calorific) + log1p(ProdVal) + CropGroup , data=data)
# pred1 <- predict(m1, newdata=data)
# 
# summary(glm1)
# summary(glm1)$deviance
# summary(glm1)$df.residual
# 
# 
# 1 - pchisq(summary(glm1)$deviance,
#            summary(glm1)$df.residual)
# 
# cbind(data, 
#       Mean = predict(glm1, newdata = data, type = "response"), 
#       SE = predict(glm1, newdata = data, type = "response", se.fit = T)$se.fit
# )

#######################################################################################

plot(log1p(data$Nutrient), log1p(data$Citation), pch=20, cex=1, cex.axis=1.25, cex.lab=1.25,col=c("black","gray")[as.numeric(data$CropGroup)], xlab="Nutrient Output (log)", ylab="Number of Publications (log)")
abline(coefficients(m2)[1:2], col="black", lwd=2)
abline(coefficients(m2)[1:2]+ c(coefficients(m2)[1], 0), col="gray", lwd=2)

data.reg <- cbind(pred1,data1)#data1 changed to data.reg to coincide with regional predictions
#data.reg <- na.omit(data.reg[(data.reg$pred1) > (data.reg$Citation + 1.5) ,])#No climate change
#data.reg <- na.omit(data.reg[(data.reg$pred1) > (data.reg$Citation + 1.5) & data.reg$Climate == -1,])#CC Loser
data.reg <- na.omit(data.reg[(data.reg$pred1) > (data.reg$Citation + 1.5) & data.reg$Climate == 1,])#CC Winner
#data.reg <- na.omit(data.reg[(data.reg$pred1) < (data.reg$Citation - 1.5) & data.reg$Climate == -1,])#exceed
data.reg <- cbind(data.reg, (data.reg$pred1) - 1 - data.reg$Citation)
data.reg <- data.reg[order(data.reg[,13]),]

data.glob <- cbind(pred2,data1)#data2 changed to data.glob to coincide with global predictions
#data.glob <- na.omit(data.glob[(data.glob$pred) > (data.glob$Citation + 1.5),])#No climate change
#data.glob <- na.omit(data.glob[(data.glob$pred) > (data.glob$Citation + 1.5) & data.glob$Climate == -1,])#CC Loser
data.glob <- na.omit(data.glob[(data.glob$pred) > (data.glob$Citation + 1.5) & data.glob$Climate == 1,])#CC Winner
#data.glob <- na.omit(data.glob[(data.glob$pred) < (data.glob$Citation - 1.5) & data.glob$Climate == -1,])#exceed
data.glob <- cbind(data.glob, (data.glob$pred2) - 1 - data.glob$Citation)
data.glob <- data.glob[order(data.glob[,13]),]

CRTB2 <- data.glob[data.glob$CropGroup == "CRTB",]
CRTB1 <- data.reg[data.reg$CropGroup == "CRTB",]

CRTB1 <- CRTB1[order(CRTB1[,4], decreasing=TRUE),]
CRTB1 <- CRTB1[order(as.numeric(CRTB1$Nutrient), decreasing=TRUE),]
TopCRTB1CitationShortage <- tapply(CRTB1[,13], as.character(CRTB1$Crop), sum)
TopCRTB1Regions <- tapply(as.character(CRTB1$Region), as.character(CRTB1$Crop), function(x) paste(gsub("([.])", " ", x), collapse=", "))
TopCRTB1Regions <- gsub("U S A ", "U.S.A.", TopCRTB1Regions) 
TopCRTB1 <- data.frame(Indicator = TopCRTB1CitationShortage, TopCRTB1Regions)
TopCRTB1 <- TopCRTB1[order(rownames(TopCRTB1), decreasing=FALSE),]

CRTB2 <- CRTB2[order(CRTB2[,4], decreasing=TRUE),]
CRTB2 <- CRTB2[order(as.numeric(CRTB2$Nutrient), decreasing=TRUE),]
TopCRTB2CitationShortage <- tapply(CRTB2[,13], as.character(CRTB2$Crop), sum)
TopCRTB2Acreage <- tapply(CRTB2$Nutrient, as.character(CRTB2$Crop), sum)
TopCRTB2Regions <- tapply(as.character(CRTB2$Region), as.character(CRTB2$Crop), function(x) paste(gsub("([.])", " ", x), collapse=", "))
TopCRTB2Regions <- gsub("U S A ", "U.S.A.", TopCRTB2Regions) 
TopCRTB2 <- data.frame(Indicator = TopCRTB2CitationShortage, TopCRTB2Regions)
TopCRTB2 <- TopCRTB2[order(rownames(TopCRTB2), decreasing=FALSE),]

pulses2 <- data.glob[data.glob$CropGroup == "pulses",]
pulses1 <- data.reg[data.reg$CropGroup == "pulses",]

pulses1 <- pulses1[order(pulses1[,4], decreasing=TRUE),]
pulses1 <- pulses1[order(as.numeric(pulses1$Nutrient), decreasing=TRUE),]
Toppulses1CitationShortage <- tapply(pulses1[,13], as.character(pulses1$Crop), sum)
Toppulses1Acreage <- tapply(pulses1$Nutrient, as.character(pulses1$Crop), sum)
Toppulses1Regions <- tapply(as.character(pulses1$Region), as.character(pulses1$Crop), function(x) paste(gsub("([.])", " ", x), collapse=", "))
Toppulses1Regions <- gsub("U S A ", "U.S.A.", Toppulses1Regions) 
Toppulses1 <- data.frame(Indicator = Toppulses1CitationShortage, Toppulses1Regions)
Toppulses1 <- Toppulses1[order(rownames(Toppulses1), decreasing=FALSE),]

pulses2 <- pulses2[order(pulses2[,4], decreasing=TRUE),]
pulses2 <- pulses2[order(as.numeric(pulses2$Nutrient), decreasing=TRUE),]
Toppulses2CitationShortage <- tapply(pulses2[,13], as.character(pulses2$Crop), sum)
Toppulses2Acreage <- tapply(pulses2$Nutrient, as.character(pulses2$Crop), sum)
Toppulses2Regions <- tapply(as.character(pulses2$Region), as.character(pulses2$Crop), function(x) paste(gsub("([.])", " ", x), collapse=", "))
Toppulses2Regions <- gsub("U S A ", "U.S.A.", Toppulses2Regions) 
Toppulses2 <- data.frame(Indicator = Toppulses2CitationShortage, Toppulses2Regions)
Toppulses2 <- Toppulses2[order(rownames(Toppulses2), decreasing=FALSE),]


write.csv(TopCRTB1, "/GLM_FAOProd_Quant_CRTB1_AllLog_NRF93_CC_Lose_Exceedence.csv")
write.csv(TopCRTB2, "/GLM_FAOProd_Quant_CRTB2_AllLog_NRF93_CC_Lose_Exceedence.csv")
write.csv(Toppulses1, "/GLM_FAOProd_Quant_Pulses1_AllLog_NRF93_CC_Lose_Exceedence.csv")
write.csv(Toppulses2, "/GLM_FAOProd_Quant_Pulses2_AllLog_NRF93_CC_Lose_Exceedence.csv")

############################################################################################
###Step 8: Plotting#########################################################################
############################################################################################


#######################################
######Figure 2#########################
#######################################


pdf(file="/Climate_Change_Map_7x4.pdf", width=8, height=9)#width=1650,height=1250, res=150)
#tiff(file="Test.tiff", width=1650,height=1250, res=150)
par(mfrow=c(7,4), oma=c(0,0,0,0), mar=c(0,0,0,0))

Crops.Uni<-unique(rownames(NModels2016))

for (i in 1:length(Crops.Uni))
{
  Crop<-NModels2016[i,]
  Crop[is.na(Crop)]<-100
  Crop<-as.numeric(Crop)
  regions@data$Crop<-Crop
  Presence<-Ocurrence2016[i,]
  Presence[is.na(Presence)]<-0
  Presence<-as.numeric(Presence)
  regions@data$presence<-Presence
  
  regions@data$COLOUR <- "#FFFFFF"
  regions@data$COLOUR[which(regions$Crop==-1 & regions$presence>=1000)]<-"#ffffbf"
  regions@data$COLOUR[which(regions$Crop==0 & regions$presence>=1000)]<-"#E0E0E0"
  regions@data$COLOUR[which(regions$Crop==1 & regions$presence>=1000)]<-"#998ec3"
  
  plot(regions1, col=regions@data$COLOUR)
  title(sapply(Crops.Uni[i], simpleCap), adj=0.08, line=-1.2, cex.main=1.2)
}  

add_legend("bottom", legend=c("Increase", "Neutral", "Decrease","No Data"), pch=15, 
           col=c("#998ec3", "#E0E0E0", "#ffffbf", "#FFFFFF"),
           horiz=TRUE, bty='n', cex=1.6)
dev.off()


#######################################
######Figure 3#########################
#######################################

png(file="/Figure2_CC_NRF.png", width=1200, height=700, res=80)
plot(data1$Nutrient, data1$Citation+1, pch=20, cex=1, cex.axis=1.25, cex.lab=1.25, axes=FALSE, col=c("black","gray")[as.numeric(data1$CropGroup)], xlab="Nutrient output (total NRF9.3 of each crop-region)", ylab="Research intensity (number of publications)", log="xy")
x <- exp(seq(0, max(log1p(data1$Nutrient)), length.out=nrow(data1)))
nd1 <- data1
nd1$Nutrient <- x
nd1$CropGroup <- "CRTB"
nd2 <- nd1
nd2$CropGroup <- "pulses"
y1 <- predict(model.zip_gaps_global, newdata=nd1)
y2 <- predict(model.zip_gaps_global, newdata=nd2)
lines(cbind(x, y1), lwd=3)
lines(cbind(x, y2), col="gray", lwd=3)
axis(1, at=10^c(7,8,9,10,11,12,13,14), labels=expression(10^7, 10^8, 10^9,
                                                         10^10, 10^11, 10^12, 10^13, 10^14))
axis(2, at=c(1, 11, 101, 1001, 10001), labels=expression(0, 10, 100, 1000, 10000))

dev.off()

#######################################
######Figure 4#########################
#######################################

summary(model.zip_gaps_global)


library(grid)
library(ggplot2)
library(reshape2)

country_intercept<-as.data.frame(summary(model.zip_gaps_regional)$coef$count[3:44,1])
mean_intercept<-mean(country_intercept[1:42,])
model_intercept<-summary(model.zip_gaps_regional)$coef$count[1,1]
rownames(country_intercept)<-sapply(strsplit(rownames(country_intercept), split='Region', fixed=TRUE), function(x) (x[2]))
country_intercept$Crop<-rep("Regional Deviation", 42)
country_intercept$Region<-rownames(country_intercept)
colnames(country_intercept)<-c("Intercept", "Crop", "Region")
arabian<-c(0, "Regional Deviation", "Arabian.Peninsula")
country_intercept<-rbind(country_intercept,arabian)
country_intercept[,1]<-as.numeric(country_intercept[,1])
country_intercept <- cbind(country_intercept, (country_intercept$Intercept) - mean_intercept)

quantile(country_intercept$`(country_intercept$Intercept) - mean_intercept`)

country_intercept$quantile<-0
country_intercept$quantile[which(country_intercept$`(country_intercept$Intercept) - mean_intercept`<=-1.9007636)]<- 2
country_intercept$quantile[which(country_intercept$`(country_intercept$Intercept) - mean_intercept`>-1.9007636 & country_intercept$`(country_intercept$Intercept) - mean_intercept`<0.1450593)]<- 1
country_intercept$quantile[which(country_intercept$`(country_intercept$Intercept) - mean_intercept`>=-0.1450593 & country_intercept$`(country_intercept$Intercept) - mean_intercept`<=0.1450593)]<- 0
country_intercept$quantile[which(country_intercept$`(country_intercept$Intercept) - mean_intercept`>0.1450593 & country_intercept$`(country_intercept$Intercept) - mean_intercept`<=11.407450)]<- -1
country_intercept$quantile[which(country_intercept$`(country_intercept$Intercept) - mean_intercept`>1.9007636)]<- -2


data.plot <- cbind(pred1,data1)
data.plot<-data.plot[!is.na(data.plot$pred1),]
data.plot <- cbind(data.plot, resids_gap_regional)
data.plot$exp<-(data.plot$pred1) - 1 - (data.plot$Citation)

quantile(data.plot$exp)

data.plot$quantile<-0
data.plot$quantile[which(data.plot$exp<=-16.77597)]<- -2
data.plot$quantile[which(data.plot$exp>-16.77597 & data.plot$exp<1.59473)]<- -1
data.plot$quantile[which(data.plot$exp>=-1.59473 & data.plot$exp<=1.59473)]<- 0
data.plot$quantile[which(data.plot$exp>1.59473 & data.plot$exp<=16.77597)]<- 1
data.plot$quantile[which(data.plot$exp>16.77597)]<- 2


col_div_8<-c("#f6e8c3","#dfc27d","#bf812d","#8c510a","#c7eae5","#80cdc1","#35978f","#01665e")
col_div_5<-c("#8c510a", "#d8b365","#f1eef6" , "#5ab4ac", "#01665e", "white")

#ccr<-read.csv("/Users/rhysmanners/Desktop/NModels2017Baseline.csv")
ccr<-read.csv("NModels2017_106crops.csv")
arts<-read.csv("/ScTable.csv")
gcn<-read.csv("/Group Crop Names.csv", header =F)
crn<-read.csv("/Continents.csv")

crn<-as.vector(colnames(crn[,1:44]))
gcn<-as.character(gcn$V1)
arts<-subset(arts, X %in% data.plot$Crop)
arts<-arts[,-c(15, 19, 23, 30, 42, 44, 45)]
colnames(arts)<-colnames(ccr[-c(15,19,23,30,42,44,45)])
ccr<-subset(ccr, X %in% data.plot$Crop)
ccr<-ccr[,-c(15, 19, 23, 30, 42, 44, 45)]

arts<-arts[,crn]#Reordering by continents
ccr<-ccr[,crn]#Reordering by continents
arts$X<-factor(arts$X, levels=c("banana", "barley","buckwheat","cassava","maize","millet","oats","potato","plantain","quinoa","rice", "rye","sorghum","sweetpotato","taro",
                                "wheat","yam","bean","broadbean","chickpea", "cowpea","lentil","lupin","pigeonpea"))
ccr$X<-factor(ccr$X, levels=c("banana", "barley","buckwheat","cassava","maize","millet","oats","potato","plantain","quinoa","rice", "rye","sorghum","sweetpotato","taro",
                              "wheat","yam","bean","broadbean","chickpea", "cowpea","lentil","lupin","pigeonpea"))

arts<-arts[order(arts$X),]#Reordering by crop group
ccr<-ccr[order(ccr$X),]
ccr<-melt(ccr, id="X")
arts<-melt(arts, id="X")
ccr$cc[ccr$value == -1] <- sprintf('\u2193')
ccr$cc[ccr$value== 0] <- sprintf('\u002D')
ccr$cc[ccr$value == 1] <- sprintf('\u2191')
arts$cc<-ccr$cc
colnames(arts)<-c("Crop", "Region", "Publications", "CCSymbol")

CroArtRes<-merge(arts, data.plot, all=T)
CroArtRes<-merge(CroArtRes,country_intercept, all=T)
CroArtRes$quantile[is.na(CroArtRes$quantile)]<- 5
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="banana"] <- "Banana"
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="barley"] <- "Barley"
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="buckwheat"] <- "Buckwheat"
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="cassava"] <- "Cassava"
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="maize"] <- "Maize"
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="millet"] <- "Millet"
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="oats"] <- "Oats"
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="potato"] <- "Potato"
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="plantain"] <- "Plantain"
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="quinoa"] <- "Quinoa"
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="rice"] <- "Rice"
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="rye"] <- "Rye"
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="sorghum"] <- "Sorghum"
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="sweetpotato"] <- "Sweetpotato"
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="taro"] <- "Taro"
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="yam"] <- "Yam"
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="wheat"] <- "Wheat"
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="bean"] <- "Common Bean"
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="broadbean"] <- "Broadbean"
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="chickpea"] <- "Chickpea"
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="cowpea"] <- "Cowpea"
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="lentil"] <- "Lentil"
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="lupin"] <- "Lupin"
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="pigeonpea"] <- "Pigeonpea"
levels(CroArtRes$Crop)[levels(CroArtRes$Crop)=="Regional Deviation"] <- "Regional Deviation"

CroArtRes$Crop<-factor(CroArtRes$Crop, levels=c("Banana", "Barley","Buckwheat","Cassava","Maize","Millet","Oats","Potato","Plantain","Quinoa","Rice", "Rye","Sorghum","Sweetpotato","Taro",
                                                "Wheat","Yam","Broadbean","Chickpea", "Common Bean","Cowpea","Lentil","Lupin","Pigeonpea", "Regional Deviation"))

afr <- textGrob("Africa", gp=gpar(fontsize=10, fontface="bold"))
ame <- textGrob("Americas", gp=gpar(fontsize=10, fontface="bold"))
asi <- textGrob("Asia", gp=gpar(fontsize=10, fontface="bold"))
eur <- textGrob("Europe", gp=gpar(fontsize=10, fontface="bold"))
oce <- textGrob("Oceania", gp=gpar(fontsize=10, fontface="bold"))
sta <- textGrob("Starch\nCrops", gp=gpar(fontsize=9, fontface="bold"))
pul <- textGrob("Pulses", gp=gpar(fontsize=9, fontface="bold"))

png(file="Heat_Map_Quantiles_Regional_Deviation.png", width=1800,height=1150, res=150)
p<-ggplot(CroArtRes, aes(CroArtRes$Region, CroArtRes$Crop)) +
  theme_classic()+
  theme_grey(base_size = 9) +
  theme(legend.position = "right",
        legend.text=element_text(size=8),
        legend.title=element_text(size=9.5),
        legend.title.align=-0.2,
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = -330, hjust = +1),
        plot.margin = unit(c(2,1,1,1), "lines"))+
  geom_tile(aes(fill = as.factor(CroArtRes$quantile)), colour = "white") +
  geom_text(aes(label = as.character(CroArtRes$CCSymbol)), size = 4)+
  geom_vline(xintercept = c(8.5, 23.5, 36.5,41.5),
             color = "black", size=0.75)+
  geom_hline(yintercept = c(1.5,8.5),
             color = "black", size=0.75)+
  scale_colour_manual(breaks = c("-2", "-1", "0", "1", "2", "5"),
                      labels = c("", "Underresearched", "",
                                 "Overresearched", "", ""),
                      values = col_div_5) +
  scale_fill_manual(breaks = c("-2", "-1", "0", "1", "2", "N5"),
                    labels = c("", "Overresearched", "",
                               "Underresearched", "", ""),
                    values = col_div_5, name="Research Intensity")+
  scale_x_discrete("", limits = levels(CroArtRes$Region))+
  scale_y_discrete("", limits = rev(levels(CroArtRes$Crop))) +
  annotation_custom(afr,xmin=4.5,xmax=4.5,ymin=25.9,ymax=25.9) +
  annotation_custom(ame,xmin=16,xmax=16,ymin=25.9,ymax=25.9)+
  annotation_custom(asi,xmin=29.5,xmax=29.5,ymin=25.9,ymax=25.9)+
  annotation_custom(eur,xmin=38.5,xmax=38.5,ymin=25.9,ymax=25.9)+
  annotation_custom(oce,xmin=42.5,xmax=42.5,ymin=25.9,ymax=25.9)+
  annotation_custom(sta,xmin=-2,xmax=-2,ymin=18,ymax=18)+
  annotation_custom(pul,xmin=-2,xmax=-2,ymin=3.5,ymax=3.5)

gt <- ggplot_gtable(ggplot_build(p))

gt$layout$clip[gt$layout$name == "panel"] <- "off"

grid.draw(gt)

dev.off()


#######################################
######Figure 5#########################
#######################################

library(maptools)
library(rgdal)
library(rgeos)

regions <- readShapeSpatial("/Regions_gaps")
regions1 <-gSimplify(regions,tol=0.01, topologyPreserve=TRUE)

Crop_Data<-read.csv("Regional_Gaps_Percentages_Domestic_Supply.csv", row.names = 1)
regions$Starch<-as.numeric(Crop_Data[1,])
regions$Pulse<-as.numeric(Crop_Data[2,])

regions@data$COLOUR <- "#FFFFFF"
#regions@data$COLOUR[which(regions$Starch==NA)]<-"white"
regions@data$COLOUR[which(regions$Starch==0)]<-"#E0E0E0"
regions@data$COLOUR[which(regions$Starch>=0.1 & regions$Starch<=25.1)]<-"#CCE5FF"
regions@data$COLOUR[which(regions$Starch>25.1 & regions$Starch<=50.1)]<-"#66B2FF"
regions@data$COLOUR[which(regions$Starch>50.1 & regions$Starch<=75.1)]<-"#0066CC"
regions@data$COLOUR[which(regions$Starch>=75.1 & regions$Starch<=100)]<-"#003366"

regions@data$COLOUR1 <- "#FFFFFF"
#regions@data$COLOUR1[which(regions$Pulse==NA)]<-"white"
regions@data$COLOUR1[which(regions$Pulse==0)]<-"#E0E0E0"
regions@data$COLOUR1[which(regions$Pulse>=0.1 & regions$Pulse<=25.1)]<-"#E5FFCC"
regions@data$COLOUR1[which(regions$Pulse>25.1 & regions$Pulse<=50.1)]<-"#B2FF66"
regions@data$COLOUR1[which(regions$Pulse>50.1 & regions$Pulse<=75.1)]<-"#66CC00"
regions@data$COLOUR1[which(regions$Pulse>75.1 & regions$Pulse<=100)]<-"#006600"


pdf(file="/Crop_Gap_Sta_Pul_CCWin_NRF_DomeSupp_Poison.pdf", width=11, height=10)#width=1650,height=1250, res=150)
par(mfrow=c(2,1), oma=c(2,1,0,0), mar=c(1,1,1,0))

plot(regions1, col=regions@data$COLOUR)
title("a) Starch Crops", adj=0.08, line=-1.2, cex.main=1.2)
legend(-160, -5,
       legend= c("NA", "0", "0-25", "25-50", "50-75", "75-100"), 
       fill=c("white", "#E0E0E0", "#CCE5FF","#66B2FF","#0066CC", "#003366"), 
       title=c("Crops with Research Gaps (%)"),
       bty="n", # turn off the legend border
       cex=.8,
       y.intersp=0.7) # decrease the font / legend size

plot(regions1, col=regions@data$COLOUR1)
title("b) Pulses", adj=0.10, line=-1.2, cex.main=1.2)
legend(-160, -5, 
       legend= c("NA", "0", "0-25", "25-50", "50-75", "75-100"), 
       fill=c("white", "#E0E0E0", "#E5FFCC","#B2FF66","#66CC00", "#006600"), 
       title=c("Crops with Research Gaps (%)"),
       bty="n", # turn off the legend border
       cex=.8,
       y.intersp=0.7) # decrease the font / legend size

dev.off()
