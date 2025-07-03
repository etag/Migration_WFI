#Wind Favorability analysis for birds migrating over the Continental US

#This script calculates intervals between days with favorable winds by looping through the 
#wind favorability data for each year. The data compilation takes a couple of hours.

library(dplyr)
library(reshape2)
library(plyr)
library(ggplot2)
library(gridExtra)

#get dataframe with sequences of bad days for each year and combine them
for(yr in 1995:2024){
  #yr = 1995
  read_file_name <- paste0("NARR_wind_fav/NARRWindFav", yr, ".csv")
  fav <- read.csv(file=read_file_name)
  colnames(fav)[which(colnames(fav)=="col_num")] = "flyway"
  colnames(fav)[which(colnames(fav)=="row_num")] = "seas"
  fav$flyway <- "C"  #initialize all flyways as central at first
  fav$flyway[which(fav$east)] <- "E"  #designate eastern flyway
  fav$flyway[which(fav$west)] <- "W"  #designate western flyway
  fav <- fav[,-which(colnames(fav) %in% c("pLev", "land", "west", "east", "central"))] #remove unneeded columns
  
  #convert to long format
  favL <- reshape2::melt(fav, id.vars = c("lon", "lat", "id", "flyway", "seas"), variable.name = "day")
  favL <- favL[order(favL$id, favL$day),] #sort data by cell ID then by day
  days <- as.numeric(substr(x = unique(favL$day), start = 2, stop = 5))
  sprDays <- which(days < 700)
  days[sprDays] <- 1:length(sprDays) + 59 #get day of year starting March 1
  falDays <- which(days > 700)
  days[falDays] <- 1:length(falDays) + 212 #get day of year starting Aug 1
  favL$day <- mapvalues(x=favL$day, from = unique(favL$day), to = days)
  favL$day <- as.numeric(as.character(favL$day)) #have to walk through character vector first??
  favL$seas <- "spr"
  favL$seas[favL$day > 200] <- "fal"
  favL$seas <- as.factor(favL$seas)
  
  favBad <- favL[which(favL$value < 1),]
  favBad$seqLen <- 0
  favBad$seqID <- 0
  
  numRows <- 2*length(unique(favBad$id))
  favSeq <- favBad[1:numRows, which(colnames(favBad) %in% c("id", "flyway", "seas", "seqLen"))]
  favSeq$seas <- "spr" 
  seqCnt = 1
  allID <- unique(favBad$id)
  allID <- allID[order(allID)]
  
  for(i in allID) {
    #i = allID[1]
    if(i %% 500 == 0) {print(i)}
    tmp <- favBad[favBad$id==i,]
    length(which(tmp$seas=="fal"))
    length(which(tmp$seas=="spr"))
    tmp$seqID <- cumsum(c(1, abs(tmp$day[-length(tmp$day)] - tmp$day[-1]) > 1))
    tmp$seqLen <- mapvalues(x = tmp$seqID, from = unique(tmp$seqID), to = table(tmp$seqID))
    tmp <- tmp[!duplicated(tmp$seqID),]
    favSeq$id[c(seqCnt, seqCnt+1)] = i
    favSeq$flyway[c(seqCnt, seqCnt+1)] <- tmp$flyway[1]
    favSeq$seqLen[seqCnt] <- mean(tmp$seqLen[tmp$seas=="spr"])
    favSeq$seqLen[seqCnt+1] <- mean(tmp$seqLen[tmp$seas=="fal"])
    if(is.nan(favSeq$id[seqCnt])) {favSeq$seqLen[seqCnt]=0}     #remove NaN in cases where ther are no unfavorable days 
    if(is.nan(favSeq$id[seqCnt+1])) {favSeq$seqLen[seqCnt+1]=0} #remove NaN in cases where ther are no unfavorable days 
    favSeq$seas[seqCnt+1] <- "fal"
    seqCnt = seqCnt + 2
  }
  
  if(yr == 1995) {
    allSeq <- favSeq
  } else {
    allSeq <- rbind(allSeq, favSeq)
  }
}

allSeq$seqLen[which(is.na(allSeq$seqLen))] = 0
write.csv(allSeq, file = "unfavorable_intervals.csv", row.names = FALSE, quote = FALSE)
allSeq <- read.csv(file = "unfavorable_intervals.csv")

allSeq$flyway <- as.factor(allSeq$flyway)
allSeq$seas <- as.factor(allSeq$seas)
allSeq$seqLen <- as.numeric(allSeq$seqLen)
allSeq$fs <- paste0(allSeq$flyway, allSeq$seas)
allSeq$seas <- mapvalues(x=allSeq$seas, from=c("fal", "spr"), to=c("Fall", "Spring"))
int_means <- allSeq %>% dplyr::group_by(flyway, seas) %>% 
  dplyr::summarise(mean=mean(seqLen), sd = sd(seqLen), median = median(seqLen))

cSeq <- allSeq[allSeq$flyway=="C",]
cMeans <- int_means[int_means$flyway=="C",]
wSeq <- allSeq[allSeq$flyway=="W",]
wMeans <- int_means[int_means$flyway=="W",]
eSeq <- allSeq[allSeq$flyway=="E",]
eMeans <- int_means[int_means$flyway=="E",]


uC <- ggplot(cSeq,aes(x=seqLen, fill=seas, legend=TRUE)) + 
  theme_classic()+
  xlim(0, 10) +
  ylim(0, 0.75) +
  geom_density(alpha=0.25) + 
  geom_vline(data=cMeans, aes(xintercept=median, colour=seas), linetype="dashed", size=0.5) +
  labs(fill = "Season", col = "Season", ) +
  scale_fill_manual(values = c("gray30", "gray90")) +
  scale_color_manual(values = c("gray10", "gray70")) +
  theme(legend.position = "none") +
  xlab("    ") +
  ylab("Density") +
  annotate("text", x = 7, y = 0.75, label="Central Flyway")

uW <- ggplot(wSeq,aes(x=seqLen, fill=seas, legend=TRUE)) + 
  theme_classic()+
  xlim(0, 10) +
  ylim(0, 0.75) +
  geom_density(alpha=0.25) + 
  geom_vline(data=wMeans, aes(xintercept=median, colour=seas), linetype="dashed", size=0.5) +
  scale_fill_manual(values = c("gray30", "gray90")) +
  scale_color_manual(values = c("gray10", "gray70")) +
  theme(legend.position = "none") +
  xlab("Unfavorable Interval (days)") +
  ylab("Density") +
  annotate("text", x = 7, y = 0.75, label="Western Flyway")
  

uE <- ggplot(eSeq,aes(x=seqLen, fill=seas, legend=TRUE)) + 
  theme_classic()+
  xlim(0, 10) +
  ylim(0, 0.75) +
  geom_density(alpha=0.25) + 
  geom_vline(data=eMeans, aes(xintercept=median, colour=seas), linetype="dashed", size=0.5) +
  labs(fill = "Season", col = "Season", ) +
  scale_fill_manual(values = c("gray30", "gray90")) +
  scale_color_manual(values = c("gray10", "gray70")) +
  theme(legend.position = c(0.8, 0.5))+
  xlab("    ") +
  ylab("Density") +
  annotate("text", x = 7, y = 0.75, label="Eastern Flyway")
  
pdf(file = "Favorable_wind_intervals.pdf", width = 3, height = 6)
grid.arrange(uE, uC, uW, nrow = 3)
dev.off()

png(file = "Favorable_wind_intervals.png", width = 3, height = 6, units = "in", pointsize = 300)
grid.arrange(uE, uC, uW, nrow = 3)
dev.off()

