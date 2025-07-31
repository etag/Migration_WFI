#Preferred Directiono of Movement analysis for birds migrating over the Continental US

#This R script generated PDM values for migrating birds in the fall at each NEXRAD radar 
#station based on filtered radar data


library(data.table)
library(dplyr)
files=list.files("/Volumes/nighthawk/merged_data", pattern = "forecast")

#files = files[c(1:48)] #for testing 2 sites

rad_data <- as.list(rep(NA, 1200))
for(f in 1:length(files)){
  print(f)
  rad_data_year=fread(paste("/Volumes/nighthawk/merged_data/",files[f], sep=""))
  rad_data_year = rad_data_year[,c(1:3,6:11,13:14,21:22,25)]
  rad_data[[f]] <- rad_data_year
}

rad_data=rad_data[1:length(files)]
rad_data=rbindlist(rad_data)

rad_data$season = ifelse(rad_data$ordinal.date <= (166), "spring",
                         ifelse(rad_data$ordinal.date >=213, "fall", "summer"))
rad_data=rad_data%>%filter(season !="summer")

rad_data$year = substr(rad_data$time_stamp, 1,4)

s = "fall"
fall=rad_data%>%filter(season== s)

#this will give the airspeed
fall$headingu <- fall$u-fall$uwnd
fall$headingv <- fall$v-fall$vwnd
fall$airspeed<- sqrt(fall$headingu^2+fall$headingv^2)

#this will give us the track direction
fall$track = (180/pi)*(atan2(fall$u,fall$v))

#this will give us the heading direction
fall$heading = (180/pi)*(atan2(fall$headingu,fall$headingv))

hist(fall$track)
hist(fall$heading)

#flip direction to 0 to 360
fall$track=ifelse(fall$track <0, (fall$track+360), fall$track)
fall$heading=ifelse(fall$heading <0, (fall$heading+360), fall$heading)

#changed for fall
fall=fall%>%filter(track<= 270 & track>= (90))
fall=fall%>%filter(airspeed >=5 & airspeed<=(30))
fall=fall%>%filter(rmse <=5)
fall=fall%>%filter(percent_rain <=25)
fall=fall%>%filter(houraftersunset>0)
fall=fall%>%filter(houraftersunset<8)
fall= na.omit(fall)

remaining_height_bins = fall %>% group_by(file) %>%
  summarise(n_height_bins = n())%>%
  filter(n_height_bins>=10)
fall=right_join(fall, remaining_height_bins, by ="file")

#format date
fall$date=as.Date(substr(fall$file, 5, 12), format="%Y%m%d")

fall_cond=fall %>% group_by(file) %>%
  summarise(radar_id= unique(radar_id),
            u = weighted.mean(u, linear_eta, na.rm = T),
            v = weighted.mean(v, linear_eta, na.rm = T),
            uwnd = weighted.mean(uwnd, linear_eta, na.rm = T),
            vwnd = weighted.mean(vwnd, linear_eta, na.rm = T),
            sum_linear_eta = sum(linear_eta),
            time_stamp= unique(time_stamp),
            houraftersunset= mean(houraftersunset),
            lat= mean(lat),
            lon= mean(lon),
            year= unique(year),
            date= unique(date))

fall_samples_per_night=fall_cond %>% group_by(date) %>%
  summarise(number_of_samples_per_night= n()) %>%
  filter(number_of_samples_per_night>10)

fall_cond=right_join(fall_cond, fall_samples_per_night, by ="date")

fall_cond$sampling_period=if_else(as.numeric(substr(fall_cond$file, 14,15))>=0 & as.numeric(substr(fall_cond$file, 14,15))<20,
                                    fall_cond$date-1, fall_cond$date)

fall_cond$airspeed_u <-  fall_cond$u - fall_cond$uwnd
fall_cond$airspeed_v <-  fall_cond$v- (fall_cond$vwnd)
fall_cond$airspeed<- sqrt(fall_cond$airspeed_u^2+fall_cond$airspeed_v^2)
hist(fall_cond$airspeed)

fall_cond$groundspeed<- sqrt(fall_cond$u^2+fall_cond$v^2)
hist(fall_cond$groundspeed)

#track
fall_cond$track= (180/pi)*(atan2(fall_cond$u, fall_cond$v))
hist(fall_cond$track)
#flip direction to 0 to 360
fall_cond$track=if_else(fall_cond$track < 0, (fall_cond$track+360), (fall_cond$track))
hist(fall_cond$track)

#wind direction
fall_cond$wind_dir= (180/pi)*(atan2(fall_cond$uwnd, fall_cond$vwnd))
#flip direction to 0 to 360
fall_cond$wind_dir=if_else(fall_cond$wind_dir < 0, (fall_cond$wind_dir+360), (fall_cond$wind_dir))

#heading
fall_cond$heading = (180/pi)*(atan2(fall_cond$airspeed_u, fall_cond$airspeed_v))
fall_cond$heading=if_else(fall_cond$heading < 0, (fall_cond$heading+360), (fall_cond$heading))

fall_cond=fall_cond%>%filter(airspeed >=5 & airspeed<=(30))

#for spring this is used to flip, but for fall it is commented out
#spring_cond$heading=ifelse(spring_cond$heading <= 180, (spring_cond$heading+180), (spring_cond$heading-180))
#spring_cond$track=ifelse(spring_cond$track <= 180, (spring_cond$track+180), (spring_cond$track-180))
fall_cond$alpha= fall_cond$track- fall_cond$heading
hist(fall_cond$alpha)

fall_cond=subset(fall_cond, alpha>(-120) & alpha<120)
fall_cond$alpha.scaled = scale(fall_cond$alpha,scale=T,center=T)

fall_cond$cbrt_linear_eta=(fall_cond$sum_linear_eta+0.001)^(1/3)
hist(fall_cond$cbrt_linear_eta)

nightly_sampling=fall_cond %>% group_by(radar_id, sampling_period) %>%
  summarise(nsamples_night= n())%>%
  filter(nsamples_night>=10)
fall_cond=right_join(fall_cond, nightly_sampling, by =c("radar_id", "sampling_period"))
fall_cond$ordinal = yday(fall_cond$sampling_period)

fwrite(fall_cond, "/Users/kylehorton/Desktop/drift/fall_drift_model_data.csv", row.names = F)
library(lme4)
mod_fall=lmer(track ~ 0 + radar_id + alpha:radar_id +
                  (alpha.scaled|ordinal) +  (alpha.scaled|radar_id:year) + (alpha.scaled|radar_id:year:ordinal),
                data=fall_cond,REML=T,
                control=lmerControl(optCtrl=list(maxfun=70000)),weight=cbrt_linear_eta)


Vcov <- vcov(mod_fall, useScale = FALSE)
betas <- fixef(mod_fall)
se <- sqrt(diag(Vcov))
zval <- betas / se
pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
values=cbind(betas, se, zval, pval)

station=as.data.frame(substr(row.names(values), 9,12))
colnames(station)= "radar_id"
values=cbind(station,values)
values=values[1:nrow(values),]
mod_values=values

PDM=values[1:143,1:3]
#PDM=values[1:2,1:3] #for testing 2 sites
#PDM$betas=ifelse (PDM$betas <= 180, (PDM$betas+180), (PDM$betas-180)) #this flips the PDMS for spring
colnames(PDM)[2]="PDM"
slopealpha=values[144:286,c(1:3)]
#slopealpha=values[3:4,c(1:3)] #for testing 2 sites
colnames(slopealpha)[2]="alpha"
hist(slopealpha$alpha)
fall_drift=cbind(PDM, slopealpha)

write.csv(fall_drift, "/Users/kylehorton/Desktop/drift/fall_drift_07_09_25.csv", row.names = F)