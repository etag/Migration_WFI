#Wind Favorability analysis for birds migrating over the Continental US

#This is one of several R scripts that examines wind favorability for 
#migrating birds using the North American Regional Reanalysis data for
#the Continental US during spring and fall bird migration season.

#This script generates maps of wind favorability with ovelays of bird
#preferred direction of movment as well as wind direction and speed.

# library(shiny)
# library(shinydashboard)
# library(shinyWidgets)
# library(shinyBS)
library(ggplot2)
library(grid)
library(RColorBrewer)
library(sf)
library(dplyr)

#get dataframe of means to generate a plotting template
load("Wind_favorability_US_means.rds")
fav_means <- Wind_fav_means; rm(Wind_fav_means) #rename
load("Wind_favorability_US_stdev.rds")
fav_sd <- Wind_fav_sd; rm(Wind_fav_sd) #rename
load("Wind_Favorability_Shiny/Wind_H_means.rds")
load("Wind_Favorability_Shiny/Wind_V_means.rds")
grd <- read.csv("Wind_Favorability_Shiny/arrow_grid.csv")
load("Wind_Favorability_Shiny/pol_boundaries.rds") #loads to variable "msf"
msf <- st_as_sf(msf)
pdm <- read.csv("PDM_spring_fall.csv")

if(!file.exists("Wind_Favorability_Shiny/tile_id_template.csv")) {
  #Make a regular grid template to force data into grid for plotting
  #define bounds for spatial data
  min_x <- round(min(pdm$lon),2)
  max_x <-  round(max(pdm$lon),2)
  min_y <-  round(min(pdm$lat),2)
  max_y <-  round(max(pdm$lat),2)
  gd <- expand.grid(x = seq(from = min_x, to = max_x, by=0.3), y = seq(from = min_y, to = max_y, by=0.3)) #Make a grid spaced at 0.3 degrees
  t1 <- data.frame(id = pdm$id, lon = pdm$lon, lat = pdm$lat) #get NARR coordinates
  gd$id = 0 #Add in id column
  gd$dist = 999 #Add in distance column
  for(i in 1:nrow(gd)) {
    row <- which.min((gd$x[i]-t1$lon)^2 + (gd$y[i]-t1$lat)^2)
    gd$id[i] = t1$id[row]
    gd$dist[i] = (gd$x[i]-t1$lon[row])^2 + (gd$y[i]-t1$lat[row])^2
  }
  #remove data from grid points that exceed a minimum distance
  gd$id[which(gd$dist>0.15)] <- NA
  #save the template
  gd_save <- gd[,c("x", "y", "id")]
  write.csv(gd_save, "Wind_Favorability_Shiny/tile_id_template.csv", quote = F, row.names = F)
}
tiles <- read.csv("Wind_Favorability_Shiny/tile_id_template.csv")

#calculate some arrow coordinates
pdm <- pdm[which(pdm$id %in% grd$id),]
pdm$springX <- cos(((90-pdm$spring_pdm_deg)/180)*pi)
pdm$springY <- sin(((90-pdm$spring_pdm_deg)/180)*pi)
pdm$fallX <- cos(((90-pdm$fall_pdm_deg)/180)*pi)
pdm$fallY <- sin(((90-pdm$fall_pdm_deg)/180)*pi)

#Plot extent limits
xmin = -124
xmax = -66
ymin = 24
ymax = 50

meansRange <- c(min(fav_means[,12:ncol(fav_means)]), max(fav_means[,12:ncol(fav_means)]))
sdRange <- c(min(fav_sd[,12:ncol(fav_means)]), max(fav_sd[,12:ncol(fav_means)]))

#Size parameters for plot
xyLabSize = 16
tickLabSize = 15
legendSize = 1.5
arrow_scale = 3

#designate info columns and data columns
header_col <- which(!grepl("X", colnames(fav_means)))
val_col <- which(grepl("X", colnames(fav_means)))

#list of all relevant dates
val_col <- which(grepl("X", colnames(fav_means)))
#get max and min for setting data ranges
min_means <- min(fav_means[,val_col])
max_means <- max(fav_means[,val_col])
min_sd <- min(fav_sd[,val_col])
max_sd <- max(fav_sd[,val_col])

# 
# all_dates <- substr(colnames(fav_means)[val_col], start=2, stop = 5)
# dates1 <- seq.Date(as.Date("2021-03-01"), as.Date("2021-06-15"), by = "day")
# dates1 <- c(dates1, seq.Date(as.Date("2021-08-01"), as.Date("2021-11-15"), by = "day"))
# dates1 <- format(dates1, format= "%b %d")
# names(all_dates) <- dates1

#_________Make_plots____________


make_WFI_plot <- function(dateStart = "0501", dateEnd = "0615", plotWhat = "Means", legend=TRUE, xLable = TRUE,
                          SpringPDMArrows = T, FallPDMArrows = F, WindArrows = T) {

  #dateStart = "0501"; dateEnd = "0615"; plotWhat = "Means"; SpringPDMArrows = T; FallPDMArrows = F; WindArrows = T
  
  col1 <- which(colnames(fav_means) == paste0("X", dateStart))
  col2 <- which(colnames(fav_means) == paste0("X", dateEnd))
  colStart = min(col1,col2)
  colEnd = max(col1, col2)
  
  if(plotWhat == "Means") {
    if(colStart==colEnd) {
      row_means <- fav_means[,colStart]
    } else {
      row_means <- rowMeans(fav_means[,colStart:colEnd])
    }
    min_val=min_means; max_val=max_means
  } else {
    if(colStart==colEnd) {
      row_means <- fav_sd[,colStart]
    } else {
      row_means <- rowMeans(fav_sd[,colStart:colEnd])
    }
    min_val=min_sd; max_val=max_sd
  }
  
  dfp <- data.frame(id = fav_means$id, mean = row_means)  
  fav2 <- left_join(tiles, dfp, by = "id")
  
  #print(head(fav2))
  
  if(plotWhat == "Means") {
    scaleLim = meansRange
    colour_breaks <- c(scaleLim[1], 1 , 1.0001, scaleLim[2])
    colours <- c("aquamarine", "orchid", "chartreuse3", "orange")
    leg_title = "Mean\nWind\nFavorability"  
  } else {
    scaleLim = sdRange
    colour_breaks <- scaleLim
    colours <- c("gray80", "navy")
    leg_title <- "Wind\nFavorability\nStandard\nDeviation" 
  }
  
  #Add dummy data to two cells to make range the same for all plots (these points are covered over later)
  fav2$mean[110] = scaleLim[1]; fav2$mean[111] = scaleLim[2]

  m1 = ggplot() +
    geom_raster(data = fav2 , aes(x = x, y = y, fill = mean))+ 
    scale_fill_gradientn(leg_title, na.value = "transparent", 
                         limits = scaleLim,
                         colors = colours[c(1, seq_along(colours), length(colours))],
                         values  = c(0, scales::rescale(colour_breaks, from = range(fav2$mean, na.rm = T)), 1),
                         name = leg_title) +
    geom_sf(data=msf, fill = NA) +
    theme_bw() +
    annotate("rect", xmin = -95, xmax = -90, ymin = 22, ymax = 28, alpha = 1,fill = "white") +
    coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax))
  
    if(xLable) {
      m1 <- m1 + labs(x="Longitude", y="Latitude")
    } else {
      m1 = m1 + labs(x="", y="Latitude")
    }
  
    if(legend) {
      m1 = m1 +
        theme(legend.key.size = unit(legendSize, 'cm'),
          legend.text=element_text(size=tickLabSize),
          axis.text=element_text(size=tickLabSize),
          axis.title = element_text(size=xyLabSize),
          legend.title = element_text(size=xyLabSize),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.margin= unit(c(1, 1, -1, 1), "lines"))
    } else {
      m1 = m1 +
        theme(
              axis.text=element_text(size=tickLabSize),
              axis.title = element_text(size=xyLabSize),
              legend.position = "none",
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.margin= unit(c(1, 1, -1, 1), "lines"))
    }

  #add wind arrows that correspond so selected dates
  if(WindArrows) {
    wind_rows <- which(Wind_H$id %in% grd$id)
    row_means_H <- rowMeans(Wind_H[wind_rows,colStart:colEnd])
    row_means_V <- rowMeans(Wind_V[wind_rows,colStart:colEnd])
    Wind_plot <- cbind(Wind_H[wind_rows, header_col], "u"=row_means_H, "v"= row_means_V)
    Wind_plot <- left_join(Wind_plot, grd[,1:3], by="id")
    Wind_plot$alon = Wind_plot$u/arrow_scale + Wind_plot$x
    Wind_plot$alat = Wind_plot$v/arrow_scale + Wind_plot$y
    m1 = m1 + 
      geom_segment(data = Wind_plot, aes(x = x, y = y, xend = alon, yend = alat),
                   arrow = arrow(length = unit(0.1, "cm")), colour="black")
  }

  #setup for preferred direction of movement arrows  
  pdm2 <- pdm
  pdm2$springU = pdm2$lon + pdm2$springX/(arrow_scale/4)
  pdm2$springV = pdm2$lat + pdm2$springY/(arrow_scale/4)
  pdm2$fallU = pdm2$lon + pdm2$fallX/(arrow_scale/4)
  pdm2$fallV = pdm2$lat + pdm2$fallY/(arrow_scale/4)    

  
  #add spring preferred direciton of movement
  if(SpringPDMArrows) {
      m1 = m1 + 
      geom_segment(data = pdm2, aes(x = lon, y = lat, xend = springU, yend=springV),
                   arrow = arrow(length = unit(0.1, "cm")), colour="dodgerblue")  
  }

  #add fall preferred direction of movement
  if(FallPDMArrows) {
    m1 = m1 +
      geom_segment(data = pdm2, aes(x = lon, y = lat, xend = fallU, yend=fallV),
                   arrow = arrow(length = unit(0.1, "cm")), colour="dodgerblue")
  }
  return(m1)  
}

m1 = make_WFI_plot(dateStart = "0501", dateEnd = "0615", plotWhat = "Means", legend=TRUE, xLable = FALSE,
                   SpringPDMArrows = T, FallPDMArrows = F, WindArrows = T)
m2 = make_WFI_plot(dateStart = "0901", dateEnd = "1015", plotWhat = "Means", legend=FALSE, xLable = TRUE,
                   SpringPDMArrows = F, FallPDMArrows = T, WindArrows = T)

pdf(file = "Wind_Fav_maps.pdf", width = 8, height = 8)
grid.draw(rbind(ggplotGrob(m1), ggplotGrob(m2), size="first"))
dev.off()

png(file = "Wind_Fav_maps.png", width = 8, height = 8, units="in", res = 300)
grid.draw(rbind(ggplotGrob(m1), ggplotGrob(m2), size="first"))
dev.off()



