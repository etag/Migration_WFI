#Wind Favorability analysis for birds migrating over the Continental US

#This is the third of several R scripts that download and process
#North American Regional Reanalysis data to examine wind conditions
#over the Continental US during spring and fall bird migration season.

#This script calculates a wind favorability index for a set of NARR grid
#cells. The index is based on the angular difference between NARR wind
#values and the preferred direction of movement. 

#set working directory
setwd("~/Documents/Rcode/NARR")

library(raster) 
library(lubridate)
library(sp)
library(dplyr)
library(DescTools)
library(numDeriv)
library(pbapply)
library(ggplot2)
library(RColorBrewer)

#Functions 

#get heading from h and v vectors
get_heading <- function(y, x, plot = FALSE) { #Note y before X 
  he_rad <- atan2(y,x)
  he_deg = he_rad*180/pi - 90
  he_deg2 = ifelse(he_deg > 0, (90 - he_deg + 270), abs(he_deg))
  if(plot) {
    ext <- max(abs(c(x,y)))+1
    par(pty="s")
    plot(x=c(-ext,ext), y=c(-ext, ext), col = "white", ylab="", xlab="", asp=1) #dummy plot to define boundaries
    arrows(x = 0, y=0, x1=x, y1=0, length=0.05, col="gray")
    arrows(x = 0, y=0, x1=0, y1=y, length=0.05, col = "blue")
    arrows(x = 0, y=0, x1=x, y1=y, length=0.05, col = "brown")
    rndAng <- round(he_deg2, 1)
    text(x=-ext/1.5, y=ext/1.5, paste0("Wind ang ", rndAng))
  }
  return(he_deg2)
}

#Calculate speed with pythagorian theorem. 
get_speed <- function(y,x) {
  return (sqrt(y^2 + x^2))
}

#Calculate absolute difference between two angles in degrees
calc_delta_angle<- function(a1, a2) {
  ax <- min(c(a1, a2))
  ay <- max(c(a1, a2))
  delta <- ifelse((ay>270 & ax < 90), (360-ay)+ax, ay-ax)
  if(delta > 180) {delta=360-delta}
  return(delta)
}

#Generate wind favorability index value
calc_wind_effect4 <- function(ang, wind, fly, plt=F) { 
  #ang = 93; wind=4.1; fly=10; plt=T 
  ang = ifelse(ang > 180, 180 - ang, ang) #angle should never exceed 180
  angR <- ang/180 * pi #convert degrees to radians
  wy=wind*cos(angR) #Calculate coordinates of wind vector terminus
  wx=wind*sin(angR) #Calculate coordinates of wind vector terminus
  if(ang == 0) {fx = 0; fy = wy + fly}
  if(ang == 180) {fx = 0; fy = fly+wy}
  if(ang != 180 & ang != 0){
    #function for optimization
    fm <- function(angY) { #optimize based on advance minus drift
      (wy+fly*sin(angY)) - abs((wx-fly*cos(angY)))
    }
    grad_h <- function(angY){ #get derivative
      return(grad(fm,angY))
    }
    res <- uniroot(grad_h, c(pi/2,0), tol=1e-10) #find root (optimal value)
    angY = res$root 
    fx = (wx-fly*cos(angY)) #use root to find terminus of flight vector
    fy = (wy+fly*sin(angY)) #use root to find terminus of flight vector
  }
  if(plt) {
    ext <- max(abs(c(fly, wx, wy, fy, fx)))+1
    par(pty="s")
    plot(x=c(-ext,ext), y=c(-ext, ext), col = "white", ylab="", xlab="", asp=1) #dummy plot to define boundaries
    arrows(x = 0, y=0, x1=0, y1=fly, length=0.05, col="gray")
    arrows(x = 0, y=0, x1=wx, y1=wy, length=0.05, col = "blue")
    arrows(x = wx, y=wy, x1=fx, y1=fy, length=0.05, col = "green")
    arrows(x = 0, y=0, x1=fx, y1=fy, length=0.05, col = "brown")
    
    rndAng <- round(90 - angY/pi *180, 1)
    text(x=-ext/1.5, y=ext/1.5, paste0("Flt Ang ", rndAng))
    mvAng <- round(90 - (atan(fy/fx)/pi * 180),1)
    text(x=-ext/1.5, y=0, paste0("Move Ang ", mvAng))
    mvDist= round(sqrt(fx^2+fy^2),1)
    text(x=-ext/1.5, y=-ext/1.5, paste0("Move dist ", mvDist))
         
  }
  indx1 = fy - abs(fx)
  indx2 = indx1/fly
  return(indx2)
}

#test function
#calc_wind_effect4(ang=160, wind=6, fly=10, plt=T) 

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

#Calculation of wind favorability index is slow.
#To speed things up, generate a lookup table for all angles and windspeeds

lookupFile = "wind_fav_lookup_fly10.csv"
if(!file.exists(lookupFile)) {
  lookup = transform(expand.grid(0:180, seq(0, 35, by=0.1)), v=0)
  colnames(lookup) <- c("d", "s", "v")
  lookup$v = pbmapply(FUN=calc_wind_effect4, ang=lookup$d, wind=lookup$s, fly=10, plt=F)
  write.csv(lookup, file=lookupFile, row.names = F, quote = F)
}
lookup <- read.csv(file=lookupFile)

#Read in raster layers for preferred direction of movement
spring_pdm <- raster("spring_pdm.tif") #read in tif of preferred movement directions for spring
fall_pdm <- raster("fall_pdm.tif") #read in tif

#Index column names to distinguish metadata columns as well as spring and fall date columns
NARRV1 = read.csv("NARRWind/NARRwindV1995.csv")
colNameIndx <- colnames(NARRV1)
dateCols = grepl("X", colNameIndx)
colNameIndx[!dateCols] = 99999
colNameIndx[dateCols] = substr(x = colNameIndx[dateCols], 2, 5)
colNameIndx = as.numeric(colNameIndx)
sprInd = which(colNameIndx < 700) #generate an index for all spring columns
fallInd = which(colNameIndx > 700 & colNameIndx != 99999) #generate an index for all spring columns


#extract values from preferred direction of movement raster
pdm_spring <- raster::extract(x = spring_pdm, y = NARRV1[,1:2]) #extract PDM direction for each NARR Cell for Spring
pdm_fall <- raster::extract(x = fall_pdm, y = NARRV1[,1:2])

#Make data frame for all pdms
NARRpdm <- NARRV1
NARRpdm[,sprInd] <- pdm_spring
NARRpdm[,fallInd] <- pdm_fall
#NARRpdm[1:20,115:119] #sanity check

years = 1995:2024

for(yr in years) { #loop through all years.
  print(paste0("processing year ", yr))
  save_index_name <- paste0("NARR_wind_fav/NARRWindFav", yr, ".csv")
  if(file.exists(save_index_name)) {
    print(paste0("year ", yr, " already processed"))
  }
  if(!file.exists(save_index_name)) { 
    
    fnameV <- paste0("NARRWind/NARRwindV", yr, ".csv")
    fnameH <- paste0("NARRWind/NARRwindH", yr, ".csv")
  
    #Get horizontal and vertical wind components
    NARRV1 = read.csv(fnameV) #read in data for all relevant months in a year
    NARRH1 = read.csv(fnameH) #read in data for all relevant months in a year
  
    #Quick scatter plot to view vertical wind values
    # myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
    # ggplot(NARRV1, aes(x = lon, y = lat, color = X0301)) +
    # geom_point(size = .6) +
    # scale_colour_gradientn(colours = myPalette(100), limits=c(min(NARRV1$X0301), max(NARRV1$X0301)))
  
    #Convert matrices to vectors and calculate headings
    hVals = unlist(NARRH1[1:nrow(NARRH1), 11:ncol(NARRH1)]) #get a vector of all horizontal wind values
    vVals = unlist(NARRV1[1:nrow(NARRV1), 11:ncol(NARRV1)]) #get a vector of all horizontal wind values
    wind_dir <- get_heading(vVals, hVals) #make a vector of wind directions
    NARRdir <- NARRV1; 
    NARRdir[,11:ncol(NARRV1)] <- wind_dir #put wind direction back into a data frame.
  
    # #Plot wind directions
    # ggplot(NARRdir, aes(x = lon, y = lat, color = X0820)) +
    # geom_point(size = .6) +
    # scale_colour_gradientn(colours = myPalette(100), limits=c(min(NARRdir$X0301), max(NARRdir$X0301)))

    #calculate wind speeds
    speeds <- get_speed(vVals, hVals)   #make a vector of wind speeds
    speeds <- round(speeds, digits = 1)
    NARRSpeeds <- NARRV1
    NARRSpeeds[,11:ncol(NARRV1)] <- speeds
  
    # ggplot(NARRSpeeds, aes(x = lon, y = lat, color = X0820)) +
    #   geom_point(size = .6) +
    #   scale_colour_gradientn(colours = myPalette(100), limits=c(min(NARRSpeeds$X0301), max(NARRSpeeds$X0301)))

    NARRdelta <- NARRdir
    #calculate angle differences between wind and pdm for spring and fall - loop thorugh columns
    for(col in c(sprInd, fallInd)) {
      a1 = pmin(NARRdir[,col], NARRpdm[,col])
      a2 = pmax(NARRdir[,col], NARRpdm[,col])
      delta <- a2-a1 #get the differences 
      delta <- ifelse(delta > 180, 360-delta, delta) #account for wrap around 
      delta <- round(delta)
      NARRdelta[,col] = delta
    }
  
    #sanity check
    # NARRdir[1570, c("X0820")]
    # NARRpdm[1570, c("X0820")]
    # NARRdelta[1570, c("X0820")]
    
    #Calculate favorabilty index - loop through columns (slow step because of merge)
    Wind_fav <- NARRdelta #template data frame
    for(col in c(sprInd, fallInd)) {
      fav_df <- data.frame(n = 1:nrow(NARRdelta), d = NARRdelta[,col], s = NARRSpeeds[,col])
      fav_df2 <- merge(x=fav_df, y=lookup, by=c('d','s'), all.x=T)
      fav_df2 <- fav_df2[order(fav_df2$n),]
      Wind_fav[,col] = fav_df2$v
    }
  
    #sanity check
    # NARRdir[1570, c("X0820")]
    # NARRpdm[1570, c("X0820")]
    # NARRdelta[1570, c("X0820")]
    # NARRSpeeds[1570, c("X0820")]
    # Wind_fav[1570, c("X0820")]
    # calc_wind_effect4(ang=NARRdelta[1570, c("X0820")], wind = NARRSpeeds[1570, c("X0820")], fly = 10, plt = T)
    
    write.csv(Wind_fav, file=save_index_name, quote = F, row.names = F)
  }
}

#Sensitivity analysis for Wind Favorability index
#How does bird flight speed affect index?
#Generate plots for comparisons of different flight speeds

#build data frame for sensitivity analysis
awf <- data_frame(ang = rep(0:180, 5), 
                  wind = rep(c(4, 7, 10, 13, 16), each=181),
                  fly = 4, fav = 0) #initial flght speed 

#add more rows for different flight speeds
awf2 = awf; awf2$fly = 8
awf3 = awf; awf3$fly = 12
awf4 = awf; awf4$fly = 16
awf <- rbind(awf, awf2, awf3, awf4)
rm(awf2, awf3, awf4)

#Calculate all wind favorability values
awf$fav = mapply(FUN=calc_wind_effect4, awf$ang, awf$wind, awf$fly) 
awf$fly <- paste0(awf$fly, " m/s") #do this for plot labels

#make a facet plot for easy comparison. 
ggplot(awf) + 
  geom_line(aes(x = ang, y = fav, color = as.factor(wind))) +
  facet_grid(factor(fly, levels=c("4 m/s",  "8 m/s",  "12 m/s", "16 m/s"))~., scales = "free_y") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement  = "outside",
        # strip.text.y = element_blank(),
        panel.spacing    = unit(0, "lines"),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position=c(.6,.94),
        legend.background=element_blank(),
        legend.direction="horizontal") +
  guides(color = guide_legend(title="Wind Speed (m/s)", title.position = "top", title.hjust =0.5)) +
  xlab("Delta angle (degrees)") + ylab("Wind Favorabilty Index")
  