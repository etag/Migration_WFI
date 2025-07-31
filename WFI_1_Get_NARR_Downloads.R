#Wind Favorability analysis for birds migrating over the Continental US

#This is the first of several R scripts that download and process
#North American Regional Reanalysis data to examine wind conditions
#over the Continental US during spring and fall bird migration season.

#This initial script downloads necessary NARR files to an 
#external hard drive

#set working directory
setwd("~/Documents/Rcode/NARR/Migration_WFI/")

#no libraries needed - all base R

#Get all unique year and month combinations for 1995 to 2023
#for getting monthly NARR files. Date range for spring is 
#March to June 15. For fall - August to November 15.  
years = as.character(1995:2024)
months = c("03","04","05","06","08","09","10","11")
downloads <- expand.grid(years, months)
downloads <- paste0(downloads$Var1, downloads$Var2)
downloads <- downloads[order(downloads)]

#Establish lists of download URLs for vertical and horizontal wind components
vURLs <- paste0("https://downloads.psl.noaa.gov/Datasets/NARR/pressure/vwnd.", downloads, ".nc")
uURLs <- paste0("https://downloads.psl.noaa.gov/Datasets/NARR/pressure/uwnd.", downloads, ".nc")

options(timeout=9000) #prevent download from terminating because of timing out

#establish list of file destination names for file downloads to external hard drive
fnv <- paste0("/Volumes/MOON_EXT1/NARR_downloads/vwnd.",downloads,".nc")
fnh <- paste0("/Volumes/MOON_EXT1/NARR_downloads/hwnd.",downloads,".nc")

#for loop to get each download
for(dl in 1:length(downloads)) {
#for(dl in 1:48) {
  fn1 <- fnv[dl]
    if(file.exists(fn1)) { #skip download if file is already present
    print(paste0(fn1, " already downloaded"))
  } else {
    print(paste0("downloading ", vURLs[dl]))
    download.file(url=vURLs[dl], destfile = fn1)
  }
  fn2 <- fnh[dl]
  if(file.exists(fn2)) {
    print(paste0(fn2, " already downloaded"))
  } else { 
    print(paste0("downloading ", uURLs[dl]))
    download.file(url=uURLs[dl], destfile = fn2)
  }
}
