#Wind Favorability analysis for birds migrating over the Continental US

#This is the second of several R scripts that download and process
#North American Regional Reanalysis data to examine wind conditions
#over the Contenental US during spring and fall bird migration season.

#This script opens NARR data from netcdf files stored on an eternal 
#hard drive and generates & saves data frames for just the relevant data,
#which includes surface level wind data for spring (March 1 to June 15) 
#and fall (August 1 to November 15).

# Some library installs may be necessary
# install.packages(c("ggplot2", "lubridate", "plyr", "reshape2", "scales",
#                    "signal", "stringr", "zoo"))
# install.packages(c("operators", "topmodel", "DEoptim", "XML"))
# install.packages("~/Downloads/EcoHydRology_0.4.12.1.tar.gz", repos = NULL, type = "source")
# install_github("CentreForHydrology/Reanalysis")
# remotes::install_version("rgeos", version = "0.6-4")

#set working directory
setwd("~/Documents/Rcode/NARR")

library(RNetCDF)
library(lubridate)
library(maps)
library(sp)
library(gaiah)
library(dplyr)
library(suncalc)
library(DescTools)

#Define all unique years and months used in the study
years = as.character(1995:2024)
months = c("03","04","05","06","08","09","10","11")
downloads <- expand.grid(years, months)
downloads <- paste0(downloads$Var1, downloads$Var2)
downloads <- downloads[order(downloads)]

#establish list of file names for from downloads
fnv <- paste0("/Volumes/MOON_EXT1/NARR_downloads/vwnd.",downloads,".nc")
fnh <- paste0("/Volumes/MOON_EXT1/NARR_downloads/hwnd.",downloads,".nc")

#Make a data frame template that get's applied to every year/month combination
if(!file.exists("NARR_vwnd_long.csv")) {      #Don't overwrite if file exists
  nc <- RNetCDF::open.nc(fnv[1], write=FALSE)
  nctimes <- RNetCDF::var.get.nc(nc, variable='time', unpack=TRUE)
  secs <- nctimes * 3600.0
  datetimes.utc <- as.POSIXct(secs, origin='1800-01-01', tz='UTC')
  #actual data are in 4D array (grid_cell_x, grid_cell_y, height_layer, dateTime)
  #data <- RNetCDF::var.get.nc(nc, variable='vwnd', unpack=TRUE)
  lons <- RNetCDF::var.get.nc(nc, variable='lon', unpack=TRUE)
  lats <- RNetCDF::var.get.nc(nc, variable='lat', unpack=TRUE)
  lvls <- RNetCDF::var.get.nc(nc, variable='level', unpack=TRUE)
  RNetCDF::close.nc(nc)
  
  #make index for Narr grid
  gridID_df <- data.frame(lon = as.vector(lons), lat = as.vector(lats),
                          id = 1:length(lons), 
                          col_num = base::rep(1:nrow(lons), times=ncol(lons)),
                          row_num = base::rep(1:ncol(lons), each=nrow(lons)))
  
  # #get land surface height data on NARR GRID from a different netcdf File 
  # Source: https://psl.noaa.gov/thredds/catalog/Datasets/NARR/time_invariant/catalog.html?dataset=Datasets/NARR/time_invariant/hgt.sfc.nc
  nc2 <- RNetCDF::open.nc("hgt.sfc.nc", write=FALSE)
  hgt <- RNetCDF::var.get.nc(nc2, variable='hgt', unpack=TRUE)
  hgt <- as.vector(hgt)
  hgt[which(is.na(hgt))] = 0 #set missing values to zero
  RNetCDF::close.nc(nc2)
  
  #get NARR pressure table from csv file
  pTab <- read.csv("PressureTable.csv")
  
  gridID_df$pLev = as.vector(hgt)
  
  getLvl <- function(ht, pTab, levels) {
    indxLvl <- which.min(abs(ht-pTab$meters))
    pressLvl <- pTab$millibars[indxLvl]
    return(which.min(abs(pressLvl-levels)))
  }
  
  gridID_df$pLev <- unlist(lapply(X=gridID_df$pLev, FUN=getLvl, pTab=pTab, levels=lvls))
  gridID_df$pLev <- gridID_df$pLev ##+ 1 # uncomment to go up one pressure level
  
  #my_palette <- topo.colors(15) # in package grDevices
  #plot(gridID_df$lon, gridID_df$lat, col=my_palette[gridID_df$pLev])
  
  # Filter wrld_simpl to North America region
  wrld_simpl <- gaiah::get_wrld_simpl()
  north_america <- wrld_simpl[wrld_simpl$REGION == 19, ]
  # Convert lat and lon to SpatialPointsDataFrame
  df = gridID_df
  coordinates(df) <- ~lon+lat
  proj4string(df) <- proj4string(north_america)
  # flag points that are on land in North America
  on_land <- !is.na(over(df, as(north_america, "SpatialPolygons")))
  gridID_df$land = on_land
  rm(df)
  
  #assign grid points to flyways
  #Flyway designation:
  #WEST (West of the 103rd meridian), 
  #central (between the 103rd and 90th meridian)
  #eastern (east of the 90th meridian; 
  #see La Sorte, Fink, Hochachka, Farnsworth, et al., 2014
  gridID_df$west = gridID_df$lon < -103
  gridID_df$east = gridID_df$lon > -90
  gridID_df$central = (gridID_df$west==F) & (gridID_df$east==F)
  
  #Crop data to land masses and to NEXRAD coverage area
  NARR1 <- gridID_df[which(gridID_df$land),] #land only (US, Canada, Mexico)
  NARR1 <- NARR1[which(NARR1$lat < 49.36),]  #Narrow to north extent of nexrad range
  NARR1 <- NARR1[which(NARR1$lat > 24.54),]   #Narrow to south extent of nexrad range
  NARR1 <- NARR1[which(NARR1$lon < -66.98702),]  #Narrow extent of nexrad range
  NARR1 <- NARR1[which(NARR1$lon > -124.71),]   #Narrow extent of nexrad range
  
  rm(gridID_df) #recover some RAM
  
  #Check plot
  plot(NARR1$lon, NARR1$lat, col="white")
  points(NARR1$lon[which(NARR1$west)], NARR1$lat[which(NARR1$west)], col = "red")
  points(NARR1$lon[which(NARR1$central)], NARR1$lat[which(NARR1$central)], col = "green")
  points(NARR1$lon[which(NARR1$east)], NARR1$lat[which(NARR1$east)], col = "blue")
  
  #Write for future use
  write.csv(NARR1, file = "NARR_vwnd_long.csv", quote = F, row.names = F)
  saveRDS(datetimes.utc, file="NARR_datetimes_utc")
  
  rm(nc) #clear large files
  rm(lons)
  rm(lats)
  rm(NARR1)
} #end if statement for making template.

# Function to get sunrise time for a given latitude, longitude, and date
get_sunset <- function(lat, lon, date) {
  sun <- getSunlightTimes(date = date, lat = lat, lon = lon, keep = "sunset")
  return(sun$sunset)
}

#Make list of all days to process
if(exists("daylist")) {rm(daylist)}
yr = 1995 #any year will work
for(mo in months) {
  dinm <- days_in_month(as.numeric(mo))
  if(mo == "06" || mo == "11") {dinm = 15} #go through June 15 and Nov 15
  dd <- sprintf("%002d", 1:dinm) 
  if(!exists("daylist")) {
    daylist = paste0(mo, dd)
  } else {
    daylist = c(daylist, paste0(mo, dd))
  }
}

NARRV <- read.csv("NARR_vwnd_long.csv") #get the template
NARRV[daylist] <- NA #make columns for all days
  
#Loop through each year and month combination to construct data frames for vertical wind
for(yr in 1995:2024) {
    print(yr)
    fileNARRYr = paste0("NARRWind/NARRwindV", yr, ".csv")
    if(file.exists(fileNARRYr)) {
      print(paste0("Narr file ", fileNARRYr, " exists already"))
      next
    }
    NARRV1 <- NARRV
    fdl1 = fdl2 ="XXX"
    for(dx in daylist) {
      print(dx)
      mo = substring(dx, 1, 2)
      da = substring(dx, 3, 4)
      if(fdl1 != paste0(yr, mo)) {
        print("getting new netcdf file")
        fdl1 <- paste0(yr, mo)
        fileV <- paste0("/Volumes/MOON_EXT1/NARR_downloads/vwnd.", yr, mo,".nc")
        nc <- RNetCDF::open.nc(fileV, write=FALSE)
        dataVW <- RNetCDF::var.get.nc(nc, variable='vwnd', unpack=TRUE)
        nctimes <- RNetCDF::var.get.nc(nc, variable='time', unpack=TRUE)
        secs <- nctimes * 3600
        datetimes.utc <- as.POSIXct(secs, origin='1800-01-01', tz='UTC') #get all date times
        RNetCDF::close.nc(nc)
      }
      sundat <- data.frame(date = as.Date(paste0(yr, "-", mo, "-", da)), lat = NARRV$lat, lon =NARRV$lon)
      sun <- getSunlightTimes(data=sundat, keep = "sunset")
      windv = rep(0, nrow(NARRV1))
      for(ln in 1:nrow(NARRV1)) {
        dateIndex <- which.min(abs(sun$sunset[ln] - datetimes.utc))
        windv[ln] <- dataVW[NARRV1$col_num[ln], NARRV1$row_num[ln], NARRV1$pLev[ln], dateIndex]
      }
      NARRV1[,which(colnames(NARRV1)==dx)] = windv
    }
    write.csv(NARRV1, file = fileNARRYr, row.names = F, quote = F)
    rm(dataVW)
    rm(nc)
  }
  
  #Do loop again for Horizontal 
  for(yr in 1995:2024) {
    print(yr)
    fileNARRYr = paste0("NARRWind/NARRwindH", yr, ".csv")
    if(file.exists(fileNARRYr)) {
      print(paste0("Narr file ", fileNARRYr, " exists already"))
      next
    }
    NARRV1 <- NARRV
    fdl1 = fdl2 ="XXX"
    for(dx in daylist) {
      print(dx)
      mo = substring(dx, 1, 2)
      da = substring(dx, 3, 4)
      if(fdl1 != paste0(yr, mo)) {
        print("getting new netcdf file")
        fdl1 <- paste0(yr, mo)
        fileV <- paste0("/Volumes/MOON_EXT1/NARR_downloads/hwnd.", yr, mo,".nc")
        nc <- RNetCDF::open.nc(fileV, write=FALSE)
        dataVW <- RNetCDF::var.get.nc(nc, variable='uwnd', unpack=TRUE)
        nctimes <- RNetCDF::var.get.nc(nc, variable='time', unpack=TRUE)
        secs <- nctimes * 3600
        datetimes.utc <- as.POSIXct(secs, origin='1800-01-01', tz='UTC') #get all date times
        RNetCDF::close.nc(nc)
      }
      sundat <- data.frame(date = as.Date(paste0(yr, "-", mo, "-", da)), lat = NARRV$lat, lon =NARRV$lon)
      sun <- getSunlightTimes(data=sundat, keep = "sunset")
      windv = rep(0, nrow(NARRV1))
      for(ln in 1:nrow(NARRV1)) {
        dateIndex <- which.min(abs(sun$sunset[ln] - datetimes.utc))
        windv[ln] <- dataVW[NARRV1$col_num[ln], NARRV1$row_num[ln], NARRV1$pLev[ln], dateIndex]
      }
      NARRV1[,which(colnames(NARRV1)==dx)] = windv
    }
    write.csv(NARRV1, file = fileNARRYr, row.names = F, quote = F)
    rm(dataVW)
    rm(nc)
  }

