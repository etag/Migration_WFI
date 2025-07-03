# Migration_WFI
Analysis of wind favorability for migratory birds over the Continental United States
 
This repository contains R scripts and associated files used to generate analyses in "Assessing How Wind Affects Bird Migration Across the Continental United States Using a Novel Wind-Favorability Metric."

R scripts were run using R version 4.4.1.

To replicate our analyses, we recommend working throught the R scripts described below in the order presented:

WFI_1_Get_NARR_Downloads.R
This script downloads all NARR data necessary to calculate a wind favorabiliyt indices. The downloads take considerable time and disk space. Note that we useed and external hard drive for the data storage required by this script.

WFI_2_Generate_NARR_dataframe.R
This script opens NARR data from netcdf files stored on an eternal hard drive and generates and saves a data frames for each year just the relevant data, which includes surface level wind data.

WFI_3_Generate_wind_fav_index.R
This script calculates a wind favorability index (WFI) for NARR grid cells. The index is based on the angular difference between NARR wind values and the preferred direction of movement. Hence, the script also uses preferred direction of movement spatial rasters for spring and fall. The script also contains a sensitivy analysis to see how the WFI is afffected by differences in bird flight speed.

WFI_4_Compute_mean&variance.R
#This script uses previously calculated favorability indices to genrate means and standard deviations across grid cells, seasons, and regions (flyways)

WFI_5_Wind_favorability
