#Wind Favorability analysis for birds migrating over the Continental US

#This is one of several R scripts that examines wind favorability for 
#migrating birds using the North American Regional Reanalysis data for
#the Continental US during spring and fall bird migration season.

#This script calculates summary statistics and generates plots for
#publication.

library(reshape2) 

load("Wind_favorability_US_means.rds")
m1 = Wind_fav_means_con
rm(Wind_fav_means_con)

#Remove unneeded columns
m1 = m1[,-which(names(m1) %in% c("year", "lon", "id", "col_num", "row_num", "pLev", "land", "east", "central"))]
m1L <- reshape2::melt(m1, id.vars = c("lat", "flyway"), variable.name = "day")
colnames(m1L)[1] = "seas"
m1L$day <- as.numeric(substr(m1L$day, start = 2, stop = 5))
m1L$seas[m1L$day < 715] = "spr"
m1L$seas[m1L$day > 715] = "fal"
m1L$seas <- as.factor(m1L$seas)
SeasTest1 = aov(value ~ seas + flyway + seas*flyway, data = m1L)
summary(SeasTest1)
TukeyHSD(SeasTest1, conf.level=.95)
SprLat <- m1L %>%
  group_by(seas) %>%
  summarize_at(.vars = "value", .funs=c(mean, sd))

meanFlySeas <- m1L %>%
  group_by(flyway, seas) %>%
  summarize_at(.vars = "value", .funs=c(mean, sd))


#Start over for latitude test
load("Wind_favorability_US_means.rds")
m1 = Wind_fav_means_con
rm(Wind_fav_means_con)
#Assign data to southern and northern latitudinal bins for quick stats test
m1 = m1 %>% mutate(lat = cut(lat, breaks = c(floor(min(m1$lat)),  40, ceiling(max(m1$lat)) )))
#convert to long and test for latitudinal difference in spring
m1L <- reshape2::melt(m1, id.vars = c("lat", "flyway"), variable.name = "day")
m1L$day <- as.numeric(substr(m1L$day, start = 2, stop = 5))
m1LSpr <- m1L[m1L$day < 715,]
latTest1 = aov(value ~ lat, data = m1LSpr)
summary(latTest1)
SprLat <- m1LSpr %>%
  group_by(lat) %>%
  summarize_at(.vars = "value", .funs=c(mean, sd))
rm(m1LSpr, latTest1, m1L)

#Start over for plotting favorability across months, season, flyway, and latitudes
load("Wind_favorability_US_means.rds")
m1 = Wind_fav_means_con
rm(Wind_fav_means_con)

#Assign data to 1 degree latitudinal bins by rounding
m1$lat = as.factor(round(m1$lat, 0))
#Remove unneeded columns
m1 = m1[,-which(names(m1) %in% c("year", "lon", "id", "col_num", "row_num", "pLev", "land", "east", "central"))]

#convert to long
m1L <- reshape2::melt(m1, id.vars = c("lat", "flyway"), variable.name = "day")
m1L$day <- as.numeric(substr(m1L$day, start = 2, 3))
colnames(m1L)[which(colnames(m1L)=="day")] = "Month"

m2L <- m1L %>% group_by(lat, flyway, Month) %>%
  summarize_at(.vars = "value", .funs = mean )

#get overall means
SAmean <- mean(m2L$value[m2L$Month %in% c("Mar", "Apr", "May", "Jun")]) #mean for entire spring data set
FAmean <- mean(m2L$value[m2L$Month %in% c("Aug", "Sep", "Oct", "Nov")]) #mean for entire fall data set
#determine max and min values for plot range
maxVal <- round(max(m2L$value), 2) + 0.01
minVal <- round(min(m2L$value), 2) - 0.01

m3L <- m1L
m3L$Month = ifelse(m3L$Month < 7, 6.5, 99) #change month numbers to 98 for spring and 99 for fall
m3L <- m3L %>% group_by(lat, flyway, Month) %>%
  summarize_at(.vars = "value", .funs = mean )

m2L = rbind(m2L, m3L)
m2L$Month <- as.factor(m2L$Month)

m4L <- m1L
m4L$flyway = "A"
m4L <- m4L %>% group_by(lat, flyway, Month) %>%
  summarize_at(.vars = "value", .funs = mean )

m5L <- m1L
m5L$flyway = "A"
m5L$Month = ifelse(m5L$Month < 7, 6.5, 99) #change month numbers to 98 for spring and 99 for fall
m5L <- m5L %>% group_by(lat, flyway, Month) %>%
  summarize_at(.vars = "value", .funs = mean )

m4L <- rbind(m4L, m5L)
m4L$Month <- as.factor(m4L$Month)

SE <- m2L[m2L$flyway=="E" & m2L$Month %in% c("3","4","5","6","6.5"),]
SC <- m2L[m2L$flyway=="C" & m2L$Month %in% c("3","4","5","6","6.5"),]
SW <- m2L[m2L$flyway=="W" & m2L$Month %in% c("3","4","5","6","6.5"),]
FE <- m2L[m2L$flyway=="E" & m2L$Month %in% c("8","9","10","11","99"),]
FC <- m2L[m2L$flyway=="C" & m2L$Month %in% c("8","9","10","11","99"),]
FW <- m2L[m2L$flyway=="W" & m2L$Month %in% c("8","9","10","11","99"),]
SA <- m4L[m4L$Month %in% c("3","4","5","6","6.5"),]
FA <- m4L[m4L$Month %in% c("8","9","10","11","99"),]

sc <- hcl.colors(6, palette = "BluGrn") #spring colors
fc <- hcl.colors(6, palette = "PurpOr") #fall colors

SPRE <- ggplot(SE) + 
  lims(y = c(minVal, maxVal)) +
  theme_bw() +
  geom_line(aes(x=lat, y=value, group=Month, colour=Month, linetype = Month)) +
  scale_linetype_manual(values = c(2,3,4,5,1),
                        labels = c("Mar", "Apr","May","Jun", "Spring")) +
  scale_color_manual(values = c(sc[1:4], "black"),
                     labels = c("Mar", "Apr","May","Jun", "Spring")) +
  scale_x_discrete(breaks = seq(25,50,by=2)) +
  theme(legend.position=c(.5,.1),
      legend.margin = margin(0, 0, 0, 0), # turned off for alignment
      legend.direction="horizontal",
      legend.title=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      plot.title=element_text(hjust=0.1, vjust=-7, face='bold')) +
  ggtitle("Spring: Eastern Flyway") +
  ylab("Wind Favorabilty")

SPRW <- ggplot(SW) + 
  lims(y = c(minVal, maxVal)) +
  theme_bw() +
  geom_line(aes(x=lat, y=value, group=Month, colour=Month, linetype = Month)) +
  scale_linetype_manual(values = c(2,3,4,5,1),
                        labels = c("Mar", "Apr","May","Jun", "Spring")) +
  scale_color_manual(values = c(sc[1:4], "black"),
                     labels = c("Mar", "Apr","May","Jun", "Spring")) +
  scale_x_discrete(breaks = seq(25,50,by=2)) +
  theme(legend.position="none",
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        plot.title=element_text(hjust=0.1, vjust=-7, face='bold')) +
  ggtitle("Spring: Western Flyway") +
  ylab("Wind Favorabilty")

SPRC <- ggplot(SC) + 
  lims(y = c(minVal, maxVal)) +
  theme_bw() +
  geom_line(aes(x=lat, y=value, group=Month, colour=Month, linetype = Month)) +
  scale_linetype_manual(values = c(2,3,4,5,1),
                        labels = c("Mar", "Apr","May","Jun", "Spring")) +
  scale_color_manual(values = c(sc[1:4], "black"),
                     labels = c("Mar", "Apr","May","Jun", "Spring")) +
  scale_x_discrete(breaks = seq(25,50,by=2)) +
  theme(legend.position="none",
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        plot.title=element_text(hjust=0.1, vjust=-7, face='bold')) +
  ggtitle("Spring: Central Flyway") +
  ylab("Wind Favorabilty")

SPRA <- ggplot(SA) + 
  lims(y = c(minVal, maxVal)) +
  theme_bw() +
  geom_line(aes(x=lat, y=value, group=Month, colour=Month, linetype = Month)) +
  scale_linetype_manual(values = c(2,3,4,5,1),
                        labels = c("Mar", "Apr","May","Jun", "Spring")) +
  scale_color_manual(values = c(sc[1:4], "black"),
                     labels = c("Mar", "Apr","May","Jun", "Spring")) +
  scale_x_discrete(breaks = seq(25,50,by=2)) +
  theme(legend.position="none",
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        plot.title=element_text(hjust=0.1, vjust=-7, face='bold')) +
  ggtitle("Spring: Continental United States") +
  ylab("Wind Favorabilty")

FALE <- ggplot(FE) + 
  geom_line(aes(x = lat, y = value, group=Month, color = Month, linetype = Month)) + 
  scale_linetype_manual(values = c(2,3,4,5,1),
                        labels = c("Aug", "Sept", "Oct", "Nov", "Fall")) +
  scale_color_manual(values = c(fc[1:4], "black"),
                     labels = c("Aug", "Sept", "Oct", "Nov", "Fall")) +
  lims(y = c(minVal, maxVal)) +
  theme_bw() +
  scale_x_discrete(breaks = seq(25,50,by=2)) +
  theme(legend.position=c(.5,.83),
        legend.margin = margin(0, 0, 0, 0), # turned off for alignment
        legend.direction="horizontal",
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        plot.title=element_text(hjust=0.1, vjust=-7, face='bold')) +
  ggtitle("Fall: Eastern Flyway") +
  ylab("")

# 
# <- ggplot(FE) + 
#   lims(y = c(minVal, maxVal)) +
#   theme_bw() +
#   geom_line(aes(x = lat, y = value, group=Month, color = Month, linetype = Month)) + 
#   scale_linetype_manual(values = c(2,3,4,5,1),
#                         labels = c("Aug", "Sept", "Oct", "Nov", "Fall")) +
#   scale_color_manual(values = c(fc[1:4], "black"),
#                      labels = c("Aug", "Sept", "Oct", "Nov", "Fall")) +
#   scale_x_discrete(breaks = seq(25,50,by=2)) +
#   theme(legend.position=c(.5,.83),
#         legend.margin = margin(0, 0, 0, 0), # turned off for alignment
#         legend.direction="horizontal",
#         legend.title=element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title.y = element_blank(),
#         plot.title=element_text(hjust=0.1, vjust=-7, face='bold')) +
#   ggtitle("Fall: Eastern Flyway") +
#   ylab("")

FALW <- ggplot(FW) + 
  geom_line(aes(x = lat, y = value, group=Month, color = Month, linetype = Month)) + 
  scale_linetype_manual(values = c(2,3,4,5,1),
                        labels = c("Aug", "Sept", "Oct", "Nov", "Fall")) +
  scale_color_manual(values = c(fc[1:4], "black"),
                     labels = c("Aug", "Sept", "Oct", "Nov", "Fall")) +
  lims(y = c(minVal, maxVal)) +
  theme_bw() +
  scale_x_discrete(breaks = seq(25,50,by=2)) +
  theme(legend.position="none",
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        plot.title=element_text(hjust=0.1, vjust=-7, face='bold')) +
  ggtitle("Fall: Western Flyway") +
  ylab("")

FALC <- ggplot(FC) + 
  geom_line(aes(x = lat, y = value, group=Month, color = Month, linetype = Month)) + 
  scale_linetype_manual(values = c(2,3,4,5,1),
                        labels = c("Aug", "Sept", "Oct", "Nov", "Fall")) +
  scale_color_manual(values = c(fc[1:4], "black"),
                     labels = c("Aug", "Sept", "Oct", "Nov", "Fall")) +
  lims(y = c(minVal, maxVal)) +
  theme_bw() +
  scale_x_discrete(breaks = seq(25,50,by=2)) +
  theme(legend.position="none",
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        plot.title=element_text(hjust=0.1, vjust=-7, face='bold')) +
  ggtitle("Fall: Central Flyway") +
  ylab("")

FALA <- ggplot(FA) + 
  geom_line(aes(x = lat, y = value, group=Month, color = Month, linetype = Month)) + 
  scale_linetype_manual(values = c(2,3,4,5,1),
                        labels = c("Aug", "Sept", "Oct", "Nov", "Fall")) +
  scale_color_manual(values = c(fc[1:4], "black"),
                     labels = c("Aug", "Sept", "Oct", "Nov", "Fall")) +
  lims(y = c(minVal, maxVal)) +
  theme_bw() +
  scale_x_discrete(breaks = seq(25,50,by=2)) +
  theme(legend.position="none",
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        plot.title=element_text(hjust=0.1, vjust=-7, face='bold')) +
  ggtitle("Fall: Contenental United States") +
  xlab("Degrees Latitude") +
  ylab("")

library(gridExtra)
pdf(file = "Means_season_flyway_month_lattitude.pdf", width = 8, height = 10)
  grid.arrange(SPRE, FALE, SPRC, FALC, SPRW, FALW, SPRA, FALA , nrow = 4)
dev.off()

png(file = "Means_season_flyway_month_lattitude.png", width = 8, height = 10)
grid.arrange(SPRE, FALE, SPRC, FALC, SPRW, FALW, SPRA, FALA , nrow = 4)
dev.off()




