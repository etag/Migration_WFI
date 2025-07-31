#Make a regular grid template to force data 
#into grid for plotting

#get a data file 
fav_means <- read.csv("Wind_favorability_US_means.csv")
fav_means <- fav_means[,1:12]

#define bounds for spatial data
min_x <- round(min(fav_means$lon),2)
max_x <-  round(max(fav_means$lon),2)
min_y <-  round(min(fav_means$lat),2)
max_y <-  round(max(fav_means$lat),2)

#Make a grid spaced at 0.3 degrees
gd <- expand.grid(
  x = seq(from = min_x, to = max_x, by=0.3),
  y = seq(from = min_y, to = max_y, by=0.3)
)

#add in some random data for test plots
gd$rand = runif(n=nrow(gd), 1, 10)

#data frame from real data
t1 <- data.frame(id = fav_means$id, lon = fav_means$lon, lat =fav_means$lat)

#Add in id and distance columns
gd$id = 0
gd$dist = 999
for(i in 1:nrow(gd)) {
  dist <- 
    row <- which.min((gd$x[i]-t1$lon)^2 + (gd$y[i]-t1$lat)^2)
  gd$id[i] = t1$id[row]
  gd$dist[i] = (gd$x[i]-t1$lon[row])^2 + (gd$y[i]-t1$lat[row])^2
}

#remove data from grid points that exceed a minimum distance
gd$id[which(gd$dist>0.15)] <- NA
gd$rand[which(gd$dist>0.15)] <- NA

#test plot
ggplot() +
  geom_tile(data = gd , aes(x = x, y = y, color = rand, fill = rand))

#how many duplicates (this counts NAs I think)
sum(duplicated(gd$id))
nrow(gd)

#save the template
gd_save <- gd[,-c(3, 5)]
write.csv(gd_save, "tile_id_template.csv", quote = F, row.names = F)
