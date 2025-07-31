#Wind Favorability analysis for birds migrating over the Continental US

#This is the third of several R scripts that download and process
#North American Regional Reanalysis data to examine wind conditions
#over the Continental US during spring and fall bird migration season.

#This script calculates a wind favorability index for a set of NARR grid
#cells. The index is based on the angular difference between NARR wind
#values and the preferred direction of movement. 

#set working directory
setwd("~/Documents/Rcode/NARR/Migration_WFI/")

library(raster) 
library(lubridate)
library(sf)
library(dplyr)
library(DescTools)
library(numDeriv)
library(ggplot2)
library(RColorBrewer)
library(gstat)

####______Functions______#### 
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
#Slow method - uses numerical approximation to find the root of a funciton.
calc_WFI_slow <- function(ang, wind, fly, plt=F) {  
  #ang = angular offset of wind from preferred flight direction
  #wind velocity (should be same units as fly)
  #fly = assumed flight velocity (groundspeed) 
  
  ang = ifelse(ang > 180, 180 - ang, ang) #angle should never exceed 180
  angR <- abs(ang/180 * pi) #convert degrees to radians - only need positive angles.
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
    par(pty="s", mar=c(2,1,1,0))
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

#Generate wind favorability index value
#Fast method - no numerical approximation. 
calc_WFI_fast <- function(ang, wind, fly, plt=F) { 
  #ang = angular offset of wind from preferred flight direction
  #wind velocity (should be same units as fly)
  #fly = assumed flight velocity (groundspeed) 
  ang = ifelse(ang > 180, 180 - ang, ang) #angle should never exceed 180
  angR <- abs(ang/180 * pi) #convert degrees to radians and get absolute value (no need for negative angles except maybe for graphing)
  wx <- wind*sin(angR) #Calculate x coordinate of wind vector terminus
  wy <- wind*cos(angR) #Calculate y coordinate of wind vector terminus
  if(ang == 0) {fx = 0; fy = wy + fly; angY = 0}
  if(ang == 180) {fx = 0; fy = fly+wy; angY = 0}
  if(ang != 180 & ang != 0){
    angY <- suppressWarnings(expr = acos(wx/fly)) #Calculaate proposed angle back to preferred direction. May be unsolveable (NaN)
    if(is.nan(angY) | angY < pi/4) {angY = pi/4} #If proposed flight angle is unsolved or less than 45, then set angle to 45.
    fx = wx-fly*cos(angY) #use flight angle to find flight vector - subtract from wind x coordinate to find movement terminus
    fy = wy+fly*sin(angY) #use flight angle to find terminus of movement vector
  }
  if(plt) {
    ext <- max(abs(c(fly, wx, wy, fy, fx)))+1
    par(pty="s", mar=c(2,1,1,0))
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

#test/compare WFI functions
# ang1 = 30; wind1=12; fly1=20
# calc_WFI_fast(ang=ang1, wind=wind1, fly=fly1, plt=T) 
# calc_WFI_slow(ang=ang1, wind=wind1, fly=fly1, plt=T) 

####____Generate interpolation raster of wind favorability from radar stations____####
#Interpolate PDM from radar point data
spring_pdm_radar <- read.csv("spring_drift_07_09_25.csv") #read in file with spring PDM for each radar station
fall_pdm_radar <- read.csv("fall_drift_07_09_25.csv") #read in file with fall PDM for each radar station
radars <- read.csv("Radar_Coordinates.csv") #file with all radar coordinates
spring_pdm_radar <- merge(x = spring_pdm_radar, y = radars, by.x="radar_id" ) #merge radar locations with spring PDM data
spring_pdm_radar$radians <- spring_pdm_radar$PDM/180 *pi #convert PDM to radians
spring_pdm_radar$sin <- sin(spring_pdm_radar$radians)    #Get vertical vector
spring_pdm_radar$cos <- cos(spring_pdm_radar$radian)     #Get horizontal vector
NARR1 <- read.csv("NARRWind/NARRwindV1995.csv")  #get some NARR grid coordinates from a wind data file
NARR1 <- NARR1[,1:3] #reduce to lon, lat and id
Narr_sf <- st_as_sf(NARR1, coords = c("lon", "lat"), crs = 4326) #create sf object for NARR grid
spr1 <- st_as_sf(spring_pdm_radar, coords = c("Lon", "Lat"), crs = 4326) #create sf object for spring PDM data
spring_pdm_sin <- idw(formula = sin~1, locations = spr1, newdata = Narr_sf ) #interpolate spring PDM vertical vector
spring_pdm_cos <- idw(formula = cos~1, locations = spr1, newdata = Narr_sf ) #interpolate spring PDM horizontal vector
spring_pdm <- atan2(spring_pdm_sin$var1.pred, spring_pdm_cos$var1.pred) #vectors back to angle
spring_pdm_deg <- ifelse(spring_pdm<0, 360+(spring_pdm/pi*180), spring_pdm/pi*180) #Convert angle to degrees
#Same process for fall PDM
fall_pdm_radar <- merge(x = fall_pdm_radar, y = radars, by.x="radar_id" )
fall_pdm_radar$radians <- fall_pdm_radar$PDM/180 *pi
fall_pdm_radar$sin <- sin(fall_pdm_radar$radians)
fall_pdm_radar$cos <- cos(fall_pdm_radar$radian)
fall1 <- st_as_sf(fall_pdm_radar, coords = c("Lon", "Lat"), crs = 4326) #create sf object
fall_pdm_sin <- idw(formula = sin~1, locations = fall1, newdata = Narr_sf ) #interpolate spring PDM
fall_pdm_cos <- idw(formula = cos~1, locations = fall1, newdata = Narr_sf ) #interpolate spring PDM
fall_pdm <- atan2(fall_pdm_sin$var1.pred, fall_pdm_cos$var1.pred)
fall_pdm_deg <- ifelse(fall_pdm<0, 360+(fall_pdm/pi*180), fall_pdm/pi*180)

NARRpdm = cbind(NARR1, spring_pdm_deg, fall_pdm_deg) #Combine spring and fall PDMs grids into one data frame
write.csv(NARRpdm, file="PDM_spring_fall.csv", quote = FALSE, row.names = FALSE) #Save it


####____Generate Lookup Table ____####
#Calculation of wind favorability index is slow.
#To speed things up, generate a lookup table for all angles and windspeeds

years = 1995:2024

for(yr in years) { #loop through all years.
  if (!dir.exists("NARR_wind_fav")) {dir.create("NARR_wind_fav") }
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
    sprInd <- which(substr(colnames(NARRV1), 1, 3) %in% c("X03", "X04", "X05", "X06"))
    fallInd <- which(substr(colnames(NARRV1), 1, 3) %in% c("X08", "X09", "X10", "X11"))
    for(col in c(sprInd, fallInd)) {
      if (col <= max(sprInd)) {
        a1 = pmin(NARRdir[,col], NARRpdm$spring_pdm)
        a2 = pmax(NARRdir[,col], NARRpdm$spring_pdm) 
      } else {
        a1 = pmin(NARRdir[,col], NARRpdm$fall_pdm)
        a2 = pmax(NARRdir[,col], NARRpdm$fall_pdm) 
      }
      delta <- a2-a1 #get the differences 
      delta <- ifelse(delta > 180, 360-delta, delta) #account for wrap around 
      delta <- round(delta)
      NARRdelta[,col] = delta
    }
    
    #Calculate favorabilty index - loop through columns (slow step)
    Wind_fav <- NARRdelta #template data frame
    for(col in c(sprInd, fallInd)) {
      print(paste0("year ", yr, " column ", colnames(NARRV1)[col]))
      fav_df <- data.frame(n = 1:nrow(NARRdelta), d = NARRdelta[,col], s = NARRSpeeds[,col])
      fav_df$v = mapply(FUN=calc_WFI_fast, ang=fav_df$d, wind=fav_df$s, fly=10, plt=F)
      fav_df <- fav_df[order(fav_df$n),]
      Wind_fav[,col] = fav_df$v
    }
    
    #Change flyway data from three binaries to 1 factor 
    Wind_fav$west[Wind_fav$west] <- "W"
    Wind_fav$west[Wind_fav$east] <- "E"
    Wind_fav$west[Wind_fav$central] <- "C"
    colnames(Wind_fav)[which(colnames(Wind_fav) == "west")] = "flyway"
    Wind_fav$flyway <- as.factor(Wind_fav$flyway)
    Wind_fav <- Wind_fav[,-which(colnames(Wind_fav) %in% c("east", "central"))]
    write.csv(Wind_fav, file=save_index_name, quote = F, row.names = F)
  }
}

 #Sensitivity analysis for Wind Favorability index
#How does bird flight speed affect index?
#Generate plots for comparisons of different flight speeds

#build data frame for sensitivity analysis
awf <- data_frame(ang = rep(0:180, 5), 
                  wind = rep(c(4, 7, 10, 13, 16), each=181),
                  fly = 5, fav = 0) #initial flght speed 

#add more rows for different flight speeds
awf2 = awf; awf2$fly = 10
awf3 = awf; awf3$fly = 15
awf <- rbind(awf, awf2, awf3)
rm(awf2, awf3)

#Calculate all wind favorability values
awf$fav = mapply(FUN=calc_WFI_fast, awf$ang, awf$wind, awf$fly) 
awf$fly <- paste0(awf$fly, " m/s") #do this for plot labels

#make a facet plot for easy comparison. 
p1 <- ggplot(awf) + 
  geom_line(aes(x = ang, y = fav, color = as.factor(wind))) +
  facet_grid(factor(fly, levels=c("5 m/s",  "10 m/s",  "15 m/s"))~., scales = "fixed") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement  = "outside",
        # strip.text.y = element_blank(),
        panel.spacing    = unit(0, "lines"),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position=c(.65,.93),
        legend.key.width = unit(0.1, 'in'),
        legend.background=element_blank(),
        legend.direction="horizontal") +
  guides(color = guide_legend(title="Wind Speed (m/s)", title.position = "top", title.hjust =0.5, title.vjust =-1.3)) +
  xlab("Delta angle (degrees)") + ylab("Wind Favorabilty Index")

pdf(file="WFI_anglePlot.pdf", height = 4, width=4)
  p1
dev.off()

png(file="WFI_anglePlot.png", height = 4, width=4, units="in", res = 300)
p1
dev.off()


