#Wind Favorability analysis for birds migrating over the Continental US

#This is the fourth of several R scripts that download and process
#North American Regional Reanalysis data to examine wind conditions
#over the Continental US during spring and fall bird migration season.

#This script uses previously calculated favorability indices 
#to genrate means and standard deviations across grid cells, seasons, and 
#regions (flyways)

library(ggplot2)
library(RColorBrewer)

#generate long data frame with all fav data
Wind_fav_files <- list.files("NARR_wind_fav/", full.names = T)
nFiles <- length(Wind_fav_files)
Wind_fav1 <- read.csv(Wind_fav_files[1])

#n_val_col <- length(val_col)
x = Wind_fav_files[1] #start with first year (1995)
year <- as.numeric(substr(x = x, start = nchar(x)-7, stop = nchar(x)-4))
Wind_fav1 <- cbind(year, Wind_fav1)

#loop through the rest of the years
for(i in 2:length(Wind_fav_files)){
  print(paste0("processing ", Wind_fav_files[i]))
  Wind_fav2 <- read.csv(Wind_fav_files[i])
  #val_col <- which(grepl("X", colnames(Wind_fav1)))
  #n_val_col <- length(val_col)
  x = Wind_fav_files[i]
  year <- as.numeric(substr(x = x, start = nchar(x)-7, stop = nchar(x)-4))
  Wind_fav2 <- cbind(year, Wind_fav2)
  Wind_fav1 <- rbind(Wind_fav1, Wind_fav2)
}

rm(Wind_fav2)
#write.csv(Wind_fav1, file="Wind_favorability_US_Alldata.csv", quote = F, row.names = F)
#Wind_fav1 <- read.csv(file="Wind_favorability_US_Alldata.csv")
#Save condensed
Wind_fav_all <- Wind_fav1
colnames(Wind_fav_all)[which(colnames(Wind_fav_all)=="col_num")] = "flyway"  #replace longitude colname with flyway
Wind_fav_all$flyway <- "C"  #initialize all flyways as central at first
Wind_fav_all$flyway[which(Wind_fav_all$east)] <- "E"
Wind_fav_all$flyway[which(Wind_fav_all$west)] <- "W"
deleteCols <- which(colnames(Wind_fav_all) %in% c("row_num", "pLev", "land", "west", "east", "central"))
Wind_fav_all <- Wind_fav_all[,-deleteCols]

which(is.numeric(Wind_fav_all))

col(Wind_fav_all)

numCol <- as.numeric(which(sapply(Wind_fav_all, class) == "numeric"))

Wind_fav_all[,numCol] <- round(Wind_fav_all[numCol], 2)
save(Wind_fav_all, file="Wind_favorability_US_Alldata.rds")
rm(Wind_fav_all)


#Calculate Means across all years
val_col <- which(grepl("X", colnames(Wind_fav1)))
header_col <- which(!grepl("X", colnames(Wind_fav1)))
m1 <- Wind_fav1 %>%
  group_by(id) %>%
  summarize_at(.vars = colnames(Wind_fav1)[val_col], .funs=mean)
header <- Wind_fav1[which(Wind_fav1$year==1995), header_col]
Wind_Fav_means <- cbind(header, m1[,2:ncol(m1)])
#write.csv(Wind_Fav_means, file="Wind_favorability_US_means.csv", quote = F, row.names = F)
#Wind_Fav_means <- read.csv(file="Wind_favorability_US_means.csv")
#Save condensed
Wind_fav_means_con <- Wind_Fav_means
colnames(Wind_fav_means_con)[which(colnames(Wind_fav_means_con)=="col_num")] = "flyway"  #replace longitude colname with flyway
Wind_fav_means_con$flyway <- "C"  #initialize all flyways as central at first
Wind_fav_means_con$flyway[which(Wind_fav_means_con$east)] <- "E"
Wind_fav_means_con$flyway[which(Wind_fav_means_con$west)] <- "W"
deleteCols <- which(colnames(Wind_fav_means_con) %in% c("row_num", "pLev", "land", "west", "east", "central"))
Wind_fav_means_con <- Wind_fav_means_con[,-deleteCols]
save(Wind_fav_means_con, file="Wind_favorability_US_means.rds")

load(file="Wind_favorability_US_means.rds")
tail(Wind_fav_means_con)

#Calculate Standard Deviations across years
sd1 <- Wind_fav1 %>%
  group_by(id) %>%
  summarize_at(.vars = colnames(Wind_fav1)[val_col], .funs=sd)
Wind_Fav_sd <- cbind(Wind_Fav_means[,header_col], sd1[,2:ncol(sd1)])
#write.csv(Wind_Fav_sd, file="Wind_favorability_US_stdev.csv", quote = F, row.names = F)
#Wind_Fav_sd <- read.csv(file="Wind_favorability_US_means.csv")
#Save condensed
Wind_fav_sd_con <- Wind_Fav_sd
colnames(Wind_fav_sd_con)[which(colnames(Wind_fav_sd_con)=="col_num")] = "flyway"  #replace longitude colname with flyway
Wind_fav_sd_con$flyway <- "C"  #initialize all flyways as central at first
Wind_fav_sd_con$flyway[which(Wind_fav_sd_con$east)] <- "E"
Wind_fav_sd_con$flyway[which(Wind_fav_sd_con$west)] <- "W"
deleteCols <- which(colnames(Wind_fav_sd_con) %in% c("row_num", "pLev", "land", "west", "east", "central"))
Wind_fav_sd_con <- Wind_fav_sd_con[,-deleteCols]
save(Wind_fav_sd_con, file="Wind_favorability_US_stdev.rds")

rm(Wind_fav_sd_con)
rm(Wind_fav_means_con)
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
  #val_col <- which(grepl("X", colnames(Wind_fav1)))
  #n_val_col <- length(val_col)
  x = Wind_H_files[i]
  year <- as.numeric(substr(x = x, start = nchar(x)-7, stop = nchar(x)-4))
  Wind_H2 <- cbind(year, Wind_H2)
  Wind_H1 <- rbind(Wind_H1, Wind_H2)
}

rm(Wind_H2)
save(Wind_H1, file="Wind_H_Alldata.rds") #is this needed?


#Calculate Means
val_col <- which(grepl("X", colnames(Wind_H1)))
header_col <- which(!grepl("X", colnames(Wind_H1)))
m1 <- Wind_H1 %>%
  group_by(id) %>%
  summarize_at(.vars = colnames(Wind_H1)[val_col], .funs=mean)
header <- Wind_H1[which(Wind_H1$year==1995), header_col]
Wind_H_means <- cbind(header, m1[,2:ncol(m1)])
#write.csv(Wind_H_means, file="Wind_H_means.csv", quote = F, row.names = F)
#Make condensed file
#Wind_H_means <- read.csv(file="Wind_H_means.csv")
Wind_H_means_con <- Wind_H_means
colnames(Wind_H_means_con)[which(colnames(Wind_H_means_con)=="col_num")] = "flyway"  #replace longitude colname with flyway
Wind_H_means_con$flyway <- "C"  #initialize all flyways as central at first
Wind_H_means_con$flyway[which(Wind_H_means_con$east)] <- "E"
Wind_H_means_con$flyway[which(Wind_H_means_con$west)] <- "W"
deleteCols <- which(colnames(Wind_H_means_con) %in% c("row_num", "pLev", "land", "west", "east", "central"))
Wind_H_means_con <- Wind_H_means_con[,-deleteCols]
save(Wind_H_means_con, file="Wind_H_means.rds")

rm(m1, header)
rm(Wind_H_means_con)

#_______Do Verticle wind files
Wind_V_files <- list.files(path = "NARRwind", pattern = "NARRwindV", full.names = TRUE)
nFiles <- length(Wind_V_files)
Wind_V1 <- read.csv(Wind_V_files[1])

#n_val_col <- length(val_col)
x = Wind_V_files[1]
year <- as.numeric(substr(x = x, start = nchar(x)-7, stop = nchar(x)-4))
Wind_V1 <- cbind(year, Wind_V1)

for(i in 2:length(Wind_V_files)){
  Wind_V2 <- read.csv(Wind_V_files[i])
  #val_col <- which(grepl("X", colnames(Wind_fav1)))
  #n_val_col <- length(val_col)
  x = Wind_V_files[i]
  year <- as.numeric(substr(x = x, start = nchar(x)-7, stop = nchar(x)-4))
  Wind_V2 <- cbind(year, Wind_V2)
  Wind_V1 <- rbind(Wind_V1, Wind_V2)
}

rm(Wind_V2)
#write.csv(Wind_V1, file="Wind_V_Alldata.csv", quote = F, row.names = F)

#Calculate Means
val_col <- which(grepl("X", colnames(Wind_V1)))
header_col <- which(!grepl("X", colnames(Wind_V1)))
m1 <- Wind_V1 %>%
  group_by(id) %>%
  summarize_at(.vars = colnames(Wind_V1)[val_col], .funs=mean)
header <- Wind_V1[which(Wind_V1$year==1995), header_col]
Wind_V_means <- cbind(header, m1[,2:ncol(m1)])
#write.csv(Wind_V_means, file="Wind_V_means.csv", quote = F, row.names = F)
#Make condensed file
#Wind_V_means <- read.csv(file="Wind_V_means.csv")
Wind_V_means_con <- Wind_V_means
colnames(Wind_V_means_con)[which(colnames(Wind_V_means_con)=="col_num")] = "flyway"  #replace longitude colname with flyway
Wind_V_means_con$flyway <- "C"  #initialize all flyways as central at first
Wind_V_means_con$flyway[which(Wind_V_means_con$east)] <- "E"
Wind_V_means_con$flyway[which(Wind_V_means_con$west)] <- "W"
deleteCols <- which(colnames(Wind_V_means_con) %in% c("row_num", "pLev", "land", "west", "east", "central"))
Wind_V_means_con <- Wind_V_means_con[,-deleteCols]
save(Wind_V_means_con, file="Wind_V_means.rds")

rm(m1, header)




