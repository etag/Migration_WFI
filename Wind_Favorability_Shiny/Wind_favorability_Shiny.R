#Show maps of wind favorabilty data for bird migration

library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyBS)
library(ggplot2)
library(RColorBrewer)
library(sf)
library(dplyr)

#get dataframe of means to generate a plotting template
load("Wind_favorability_US_means.rds")  #Favorability for each grid cell averaged across all years
load("Wind_favorability_US_stdev.rds")  #standard deviation of means across years
load("Wind_H_means.rds")                #Wind horizontal means
load("Wind_V_means.rds")                #Wind vertical means
grd <- read.csv("arrow_grid.csv")       #grid for arrows
load("pol_boundaries.rds")              #polygons (cell boundaries) loads to variable "msf"
msf <- st_as_sf(msf)                    #reformat
tiles <- read.csv("tile_id_template.csv") 
pdm <- read.csv("PDM_spring_and_fall.csv") #preferred direction of movement data

#calculate some arrow coordinates
pdm <- pdm[which(pdm$id %in% grd$id),]
pdm$springX <- cos(((90-pdm$pdmSpring)/180)*pi)
pdm$springY <- sin(((90-pdm$pdmSpring)/180)*pi)
pdm$fallX <- cos(((90-pdm$pdmFall)/180)*pi)
pdm$fallY <- sin(((90-pdm$pdmFall)/180)*pi)

#Make template for main plot
plot_it <- Wind_fav_means[,c("lon", "lat", "id")]

#Plot extent limits
xmin = -124
xmax = -66
ymin = 24
ymax = 50

meansRange <- c(min(Wind_fav_means[,12:ncol(Wind_fav_means)]), max(Wind_fav_means[,12:ncol(Wind_fav_means)]))
sdRange <- c(min(Wind_fav_sd[,12:ncol(Wind_fav_means)]), max(Wind_fav_sd[,12:ncol(Wind_fav_means)]))

#Size parameters for plot
xyLabSize = 16
tickLabSize = 15
legendSize = 1.5
arrow_scale = 3

#designate info columns and data columns
header_col <- which(!grepl("X", colnames(Wind_fav_means)))
val_col <- which(grepl("X", colnames(Wind_fav_means)))

#list of all relevant dates
val_col <- which(grepl("X", colnames(Wind_fav_means)))
min_means <- min(Wind_fav_means[,val_col])
max_means <- max(Wind_fav_means[,val_col])
min_sd <- min(Wind_fav_sd[,val_col])
max_sd <- max(Wind_fav_sd[,val_col])
all_dates <- substr(colnames(Wind_fav_means)[val_col], start=2, stop = 5)
dates1 <- seq.Date(as.Date("2021-03-01"), as.Date("2021-06-15"), by = "day")
dates1 <- c(dates1, seq.Date(as.Date("2021-08-01"), as.Date("2021-11-15"), by = "day"))
dates1 <- format(dates1, format= "%b %d")
names(all_dates) <- dates1

#color palette for plots
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

ui <- fluidPage(

  sidebarLayout(
    
    sidebarPanel( 
      h4("Wind Favorabilty for Bird Migration"),
      helpText("Explore wind conditions for migratory 
               birds by selecting desired date ranges. 
               Data are compiled from the North American 
               Regional Reanalysis for the years 1995 
               to 2024. Values >1 indicated a positive 
               effect of wind, and values <1 mean that 
               wind is a detriment. Under \"Show Arrows\" 
               you can opt to display wind direction and 
               velocity as well as preferred direction of 
               movement (PDM) for birds migrating in the 
               spring and fall respectively"),
      h5("MAP 1"),
      fluidRow(
        awesomeRadio(
          inputId = "PlotType1",
          label = "Choose summary statistic", 
          choices = c("Means", "Standard Deviation"),
          selected = "Means",
          inline = TRUE,
          status = "danger"
        )
      ),
      fluidRow(
        prettyCheckboxGroup(
          inputId = "arrows",
          label = "Show arrows", 
          choices = c("Wind", "Spring PDM", "Fall PDM"),
          selected = NULL,
          inline = TRUE,
          status = "danger"
        )
      ),
      fluidRow(
        column(6, offset=0, pickerInput(inputId = "startDate1", label = "Start Date", choices = all_dates, selected=all_dates[62])),
        column(6, offset=0, pickerInput(inputId = "endDate1", label = "End Date", choices = all_dates, selected=all_dates[81]))
      ),

      helpText("-------------------------"),
      h5("MAP 2"),
      fluidRow(
        awesomeRadio(
          inputId = "PlotType2",
          label = "Choose summary statistic", 
          choices = c("Means", "Standard Deviation"),
          selected = "Means",
          inline = TRUE,
          status = "danger"
        )
      ),
      fluidRow(
        prettyCheckboxGroup(
          inputId = "arrows2",
          label = "Show arrows", 
          choices = c("Wind", "Spring PDM", "Fall PDM"),
          selected = NULL,
          inline = TRUE,
          status = "danger"
        )
      ),
      fluidRow(
        column(6, offset=0, pickerInput(inputId = "startDate2", label = "Start Date", choices = all_dates, selected=all_dates[139])),
        column(6, offset=0, pickerInput(inputId = "endDate2", label = "End Date", choices = all_dates, selected=all_dates[158]))
      )
    ),
    
    mainPanel(
      plotOutput("Map1"),
      plotOutput("Map2"),
    )
  )
)

server <- function(input, output, session) {
  
  rVals <- reactiveValues(  
    r1  = 3,
    r2 = 7
  )
  
  #Plot map 1
  output$Map1 <- renderPlot({
    
    dateStart <- input$startDate1
    dateEnd <- input$endDate1
    
    #print(dateStart)
    #print(dateEnd)
    # dateStart <- "0301"
    # dateEnd <- "0422"
    
    col1 <- which(colnames(Wind_fav_means) == paste0("X", dateStart))
    col2 <- which(colnames(Wind_fav_means) == paste0("X", dateEnd))
    colStart = min(col1,col2)
    colEnd = max(col1, col2)
    
    #print(colStart)
    #print(colEnd)
    
    if(input$PlotType1 == "Means") {
      if(colStart==colEnd) {
        row_means <- Wind_fav_means[,colStart]
      } else {
        row_means <- rowMeans(Wind_fav_means[,colStart:colEnd])
      }
      min_val=min_means; max_val=max_means
    } else {
      if(colStart==colEnd) {
        row_means <- Wind_fav_sd[,colStart]
      } else {
        row_means <- rowMeans(Wind_fav_sd[,colStart:colEnd])
      }
      min_val=min_sd; max_val=max_sd
    }
    
    dfp <- data.frame(id = Wind_fav_means$id, mean = row_means)  
    fav2 <- left_join(tiles, dfp, by = "id")
    #print(head(fav2))
    
    if(input$PlotType1 == "Means") {
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
      theme_bw()+
      theme(legend.key.size = unit(legendSize, 'cm'),
            legend.text=element_text(size=tickLabSize),
            axis.text=element_text(size=tickLabSize),
            axis.title = element_text(size=xyLabSize),
            legend.title = element_text(size=xyLabSize),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
      labs(x="Longitude", y="Latitude") +
      annotate("rect", xmin = -95, xmax = -90, ymin = 22, ymax = 28, alpha = 1,fill = "white") +
      coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax))

    if("Wind" %in% input$arrows) {
      wind_rows <- which(Wind_H$id %in% grd$id)
      row_means_H <- rowMeans(Wind_H[wind_rows,colStart:colEnd])
      row_means_V <- rowMeans(Wind_V[wind_rows,colStart:colEnd])
      Wind_plot <- cbind(Wind_H[wind_rows, header_col], "u"=row_means_H, "v"= row_means_V)
      Wind_plot$alon = Wind_plot$u/arrow_scale + Wind_plot$lon
      Wind_plot$alat = Wind_plot$v/arrow_scale + Wind_plot$lat
      m1 = m1 + 
        geom_segment(data = Wind_plot, aes(x = lon, y = lat, xend = alon, yend = alat),
                     arrow = arrow(length = unit(0.1, "cm")), colour="black")
    }
      
    if("Spring PDM" %in% input$arrows | "Fall PDM" %in% input$arrows) {
      pdm2 <- pdm
      pdm2$springU = pdm2$lon + pdm2$springX/(arrow_scale/4)
      pdm2$springV = pdm2$lat + pdm2$springY/(arrow_scale/4)
      pdm2$fallU = pdm2$lon + pdm2$fallX/(arrow_scale/4)
      pdm2$fallV = pdm2$lat + pdm2$fallY/(arrow_scale/4)    
      }
    
    if("Spring PDM" %in% input$arrows) {
      m1 = m1 + 
        geom_segment(data = pdm2, aes(x = lon, y = lat, xend = springU, yend=springV),
                     arrow = arrow(length = unit(0.1, "cm")), colour="chocolate")  
    }

    if("Fall PDM" %in% input$arrows) {
      m1 = m1 + 
        geom_segment(data = pdm2, aes(x = lon, y = lat, xend = fallU, yend=fallV),
                     arrow = arrow(length = unit(0.1, "cm")), colour="dodgerblue") 
    }
    suppressWarnings(print(m1))
  })
  
  #Plot map 2
  output$Map2 <- renderPlot({
    
    dateStart <- input$startDate2
    dateEnd <- input$endDate2
    
    #print(dateStart)
    #print(dateEnd)
    # dateStart <- "0901"
    # dateEnd <- "0922"
    
    col1 <- which(colnames(Wind_fav_means) == paste0("X", dateStart))
    col2 <- which(colnames(Wind_fav_means) == paste0("X", dateEnd))
    colStart = min(col1,col2)
    colEnd = max(col1, col2)
    
    #print(colStart)
    #print(colEnd)
    
    if(input$PlotType2 == "Means") {
      if(colStart==colEnd) {
        row_means <- Wind_fav_means[,colStart]
      } else {
        row_means <- rowMeans(Wind_fav_means[,colStart:colEnd])
      }
      min_val=min_means; max_val=max_means
    } else {
      if(colStart==colEnd) {
        row_means <- Wind_fav_sd[,colStart]
      } else {
        row_means <- rowMeans(Wind_fav_sd[,colStart:colEnd])
      }
      min_val=min_sd; max_val=max_sd
    }
    
    dfp <- data.frame(id = Wind_fav_means$id, mean = row_means)  
    fav2 <- left_join(tiles, dfp, by = "id")
    #print(head(fav2))
    
    if(input$PlotType2 == "Means") {
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
      theme_bw()+
      theme(legend.key.size = unit(legendSize, 'cm'),
            legend.text=element_text(size=tickLabSize),
            axis.text=element_text(size=tickLabSize),
            axis.title = element_text(size=xyLabSize),
            legend.title = element_text(size=xyLabSize),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
      labs(x="Longitude", y="Latitude") +
      annotate("rect", xmin = -95, xmax = -90, ymin = 22, ymax = 28, alpha = 1,fill = "white") +
      coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax))
    
    if("Wind" %in% input$arrows2) {
      wind_rows <- which(Wind_H$id %in% grd$id)
      row_means_H <- rowMeans(Wind_H[wind_rows,colStart:colEnd])
      row_means_V <- rowMeans(Wind_V[wind_rows,colStart:colEnd])
      Wind_plot <- cbind(Wind_H[wind_rows, header_col], "u"=row_means_H, "v"= row_means_V)
      Wind_plot$alon = Wind_plot$u/arrow_scale + Wind_plot$lon
      Wind_plot$alat = Wind_plot$v/arrow_scale + Wind_plot$lat
      m1 = m1 + 
        geom_segment(data = Wind_plot, aes(x = lon, y = lat, xend = alon, yend = alat),
                     arrow = arrow(length = unit(0.1, "cm")), colour="black")
    }
    
    if("Spring PDM" %in% input$arrows2 | "Fall PDM" %in% input$arrows2) {
      pdm2 <- pdm
      pdm2$springU = pdm2$lon + pdm2$springX/(arrow_scale/4)
      pdm2$springV = pdm2$lat + pdm2$springY/(arrow_scale/4)
      pdm2$fallU = pdm2$lon + pdm2$fallX/(arrow_scale/4)
      pdm2$fallV = pdm2$lat + pdm2$fallY/(arrow_scale/4)    
    }
    
    if("Spring PDM" %in% input$arrows2) {
      m1 = m1 + 
        geom_segment(data = pdm2, aes(x = lon, y = lat, xend = springU, yend=springV),
                     arrow = arrow(length = unit(0.1, "cm")), colour="chocolate")  
    }
    
    if("Fall PDM" %in% input$arrows2) {
      m1 = m1 + 
        geom_segment(data = pdm2, aes(x = lon, y = lat, xend = fallU, yend=fallV),
                     arrow = arrow(length = unit(0.1, "cm")), colour="dodgerblue") 
    }
    suppressWarnings(print(m1))
  })
}
shinyApp(ui, server)





# to generate arrow grid
# 
# density = 2
# grd <- expand.grid(
#   x = seq(from = min(joined$lon), to = max(joined$lon), by = density),
#   y = seq(from = min(joined$lat), to = max(joined$lat), by = density)
# )
# grd <- as.data.frame(grd)
# grd$id = 0
# grd$dist = 99
# 
# for(i in 1:nrow(grd)) {
#   dists <- (abs(grd$x[i]-joined$lon) + abs(grd$y[i]-joined$lat))
#   minDist <- min(dists)
#   nnIndx <- which.min(dists)
#   #print(minDist)
#   grd$id[i] = joined$id[nnIndx]
#   grd$dist[i] = minDist
# }
# grd = grd[which(grd$dist<.35),]
# grd = grd[order(grd$dist),]
# which(duplicated(grd$id))
# write.csv(grd, file ="arrow_grid.csv", row.names = F, quote = F)

