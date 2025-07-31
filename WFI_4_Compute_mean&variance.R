#Wind Favorability analysis for birds migrating over the Continental US

#This is the fourth of several R scripts that download and process
#North American Regional Reanalysis data to examine wind conditions
#over the Continental US during spring and fall bird migration season.

#This script uses previously calculated favorability indices 
#to genrate means and standard deviations across grid cells, seasons, and 
#regions (flyways). It also computes means and variances for wind V and H
#components across years. 

#generate long data frame with all fav data
Wind_fav_files <- list.files("NARR_wind_fav/", full.names = T)

#loop through each year to combine files into one
for(i in 1:length(Wind_fav_files)){
  print(paste0("processing ", Wind_fav_files[i]))
  x = Wind_fav_files[i]
  year <- as.numeric(substr(x = x, start = nchar(x)-7, stop = nchar(x)-4))
  Wind_fav2 <- read.csv(x)
  Wind_fav2 <- cbind(year, Wind_fav2)
  if(i==1) { 
    Wind_fav_all <- Wind_fav2
  } else {
    Wind_fav_all <- rbind(Wind_fav_all, Wind_fav2)
  }
}
rm(Wind_fav2)


#Calculate Means across all years
val_col <- which(grepl("X", colnames(Wind_fav_all)))
header_col <- which(!grepl("X", colnames(Wind_fav_all)))
m1 <- Wind_fav_all %>%
  group_by(id) %>%
  summarize_at(.vars = colnames(Wind_fav_all)[val_col], .funs=mean)
header <- Wind_fav_all[which(Wind_fav_all$year==1995), header_col]
Wind_fav_means <- cbind(header, m1[,2:ncol(m1)])
save(Wind_fav_means, file="Wind_favorability_US_means.rds")


#Calculate Standard Deviations across years
sd1 <- Wind_fav_all %>%
  group_by(id) %>%
  summarize_at(.vars = colnames(Wind_fav_all)[val_col], .funs=sd)
Wind_fav_sd <- cbind(Wind_fav_means[,header_col], sd1[,2:ncol(sd1)])
save(Wind_fav_sd, file="Wind_favorability_US_stdev.rds")

#Round all numneric columns to 2 decimal places for smaler file size.
numCol <- as.numeric(which(sapply(Wind_fav_all, class) == "numeric"))
Wind_fav_all[,numCol] <- round(Wind_fav_all[numCol], 2)
save(Wind_fav_all, file="Wind_favorability_US_Alldata.rds") #Save in condensed format
rm(Wind_fav_all)
rm(Wind_fav_sd)
rm(Wind_fav_means)
rm(m1, header)

#Do means for wind H and V components
#generate long data frame with all fav data
Wind_H_files <- list.files(path = "NARRwind", pattern = "NARRwindH", full.names = TRUE)
nFiles <- length(Wind_H_files)
Wind_H1 <- read.csv(Wind_H_files[1])

#n_val_col <- length(val_col)
x = Wind_H_files[1]
year <- as.numeric(substr(x = x, start = nchar(x)-7, stop = nchar(x)-4))
Wind_H1 <- cbind(year, Wind_H1)

for(i in 2:length(Wind_H_files)){
  Wind_H2 <- read.csv(Wind_H_files[i])
  x = Wind_H_files[i]
  year <- as.numeric(substr(x = x, start = nchar(x)-7, stop = nchar(x)-4))
  Wind_H2 <- cbind(year, Wind_H2)
  Wind_H1 <- rbind(Wind_H1, Wind_H2)
}
rm(Wind_H2)


#Calculate Means
val_col <- which(grepl("X", colnames(Wind_H1)))
header_col <- which(!grepl("X", colnames(Wind_H1)))
m1 <- Wind_H1 %>%
  group_by(id) %>%
  summarize_at(.vars = colnames(Wind_H1)[val_col], .funs=mean)
header <- Wind_H1[which(Wind_H1$year==1995), header_col]
Wind_H_means <- cbind(header, m1[,2:ncol(m1)])
save(Wind_H_means, file="Wind_H_means.rds")
rm(m1, header)
rm(Wind_H_means_con)

#_______Do Verticle wind files
Wind_V_files <- list.files(path = "NARRwind", pattern = "NARRwindV", full.names = TRUE)
nFiles <- length(Wind_V_files)
Wind_V1 <- read.csv(Wind_V_files[1])
x = Wind_V_files[1]
year <- as.numeric(substr(x = x, start = nchar(x)-7, stop = nchar(x)-4))
Wind_V1 <- cbind(year, Wind_V1)

for(i in 2:length(Wind_V_files)){
  Wind_V2 <- read.csv(Wind_V_files[i])
  x = Wind_V_files[i]
  year <- as.numeric(substr(x = x, start = nchar(x)-7, stop = nchar(x)-4))
  Wind_V2 <- cbind(year, Wind_V2)
  Wind_V1 <- rbind(Wind_V1, Wind_V2)
}
rm(Wind_V2)

#Calculate Means
val_col <- which(grepl("X", colnames(Wind_V1)))
header_col <- which(!grepl("X", colnames(Wind_V1)))
m1 <- Wind_V1 %>%
  group_by(id) %>%
  summarize_at(.vars = colnames(Wind_V1)[val_col], .funs=mean)
header <- Wind_V1[which(Wind_V1$year==1995), header_col]
Wind_V_means <- cbind(header, m1[,2:ncol(m1)])
save(Wind_V_means, file="Wind_V_means.rds")
rm(m1, header)




