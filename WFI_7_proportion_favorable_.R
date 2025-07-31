#Wind Favorability analysis for birds migrating over the Continental US

#This script makes calculates the proportion of favorable days
#in a season across years.

library(ggplot2)
library(RColorBrewer)
library(sf)
library(dplyr)
library(reshape2)
library(plyr)

#get dataframe of means to generate a plotting template
load("Wind_favorability_US_Alldata.rds")
fav <- Wind_fav_all; rm(Wind_fav_all)
fav <- fav[,-which(names(fav) %in% c("lon", "lat", "id", "col_num", "row_num", "pLev", "land"))]

#convert to long
favL <- reshape2::melt(fav, id.vars = c("year", "flyway"), variable.name = "day")

#add season as factor
favL$season <- ifelse(as.character(favL$day) < "X0701", "spr", "fal")
favL$season <- as.factor(favL$season)

#Convert favorability index to 1 (favorable) or 0 (not favorable) 
favL$value <- ifelse(favL$value > 1, 1, 0)

#Group by year and season
favYrSea <- favL %>%
  group_by(year, season) %>%
  summarize_at(.vars = "value", .funs=mean)
t1 <- aov(value ~ season, data=favYrSea)
summary(t1)
favYrSeaMeans <- favYrSea %>%
  group_by(season) %>%
  summarize_at(.vars = "value", .funs=c(mean, sd))
favYrSeaMeans

#Group by year, season, and flyway
favYrSeaFly <- favL %>%
  group_by(year, season, flyway) %>%
  summarize_at(.vars = "value", .funs=mean)
t1 <- aov(value ~ season + flyway +season*flyway, data=favYrSeaFly)
summary(t1)
favYrSeaFlyMeans <- favYrSeaFly %>%
  group_by(season, flyway) %>%
  summarize_at(.vars = "value", .funs=c(mean, sd))
favYrSeaFlyMeans
TukeyHSD(t1, conf.level=.95)

#Make wide data format for table with columns for each flyway
favW <- reshape2::dcast(favYrSeaFly, year ~ season + flyway)

#Add in columns for entire study area
favW$fal_A <- favYrSea$value[which(favYrSea$season=="fal")]
favW$spr_A <- favYrSea$value[which(favYrSea$season=="spr")]

#Redo column order
favW <- favW[,c(1,6,3,5,2,7,4,9,8)]

write.csv(favW, file="fav_proportion_table.csv", quote=F, row.names=F)
favW <- read.csv(file="fav_proportion_table.csv")

#Go back to long format for plotting
favPlot <- reshape2::melt(favW, id.vars = c("year"), variable.name = "seas_fly")
t1 <- aov(value ~ seas_fly, data=favPlot)
TukeyHSD(t1, conf.level=.95)

#Dissect season and flyway
favPlot$seas_fly <- as.character(favPlot$seas_fly)
favPlot$flyway <- substr(x = favPlot$seas_fly, 5, 5)
favPlot$season <-  substr(x = favPlot$seas_fly, 1, 3)

#Rename some factors
favPlot$flyway <- mapvalues(favPlot$flyway, from=c("E", "C", "W", "A"), to = c("East", "Central", "West", "All"))
favPlot$season <- mapvalues(favPlot$season, from=c("spr", "fal"), to = c("Spring", "Fall"))

#convert to percentages
favPlot$value <- favPlot$value*100

plotFav <- ggplot(favPlot) + 
  geom_line(aes(x = year, y = value, color = as.factor(season))) +
  facet_grid(factor(flyway, levels=c("E"="East",  "C"="Central",  "W"="West", "A"="All"))~.) +
  scale_color_manual(values = c("Spring" = "black", "Fall" = "gray")) +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.spacing    = unit(0, "lines"),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position=c(.7,.19),
        legend.background=element_blank(),
        legend.direction="horizontal") +
  guides(color = guide_legend(title="Season: ")) + 
  xlab("Year") + ylab("Percentage of Favorable Days")

library(gridExtra)
pdf(file = "Percent_favorable_means_by_season_and_year.pdf", width = 5, height = 5)
plotFav
dev.off()

png(file = "Percent_favorable_means_by_season_and_year.png", width = 5, height = 5, units = "in", res = 300)
plotFav
dev.off()

