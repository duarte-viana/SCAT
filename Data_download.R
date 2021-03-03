###########################################
# Supporting Information
# "Increasing climatic decoupling of bird abundances and distributions"
# Duarte S. Viana, Jonathan M. Chase
###########################################

# Get BBS and climate data
# CREDITS: code mostly taken from Harris et al. 2018 (https://doi.org/10.7717/peerj.4278)
# Some code was slightly modified

# Load functions
source("Data_management.R")
source("get_prism_data.R")

# Load libraries
library(dplyr)
library(tidyr)
library(stringr)
library(DBI)
library(rdataretriever)
library(prism)
library(raster)
library(sp)
# Issue when running package "retriever" - fix with (https://github.com/weecology/retriever/issues/1390)
# py_install("retriever")


#############################################################################################

############################
# Download species data
# NOT RUN AGAIN
rdataretriever::datasets()
rdataretriever::install_sqlite(dataset='breed-bird-survey',file='bbs_sqlite.db', debug=FALSE, use_cache=TRUE)
############################

############################
# PRISM data
# NOT RUN AGAIN
options(prism.path = "/Users/Viana/Documents/Papers/iDiv/Paper_sDiv/Data/BBS/prismdata")
years_to_use=1981:2018

download_prism()
check_if_prism_files_present(ls_prism_data(),1982:2018)

prism_stacked <- prism_stack(ls_prism_data())
plot(prism_stacked,40)
#Aggregate to 40km. Original cells average ~4.5 across N. America,
#so aggregation of 9x9 creates mostly 40km cells. 
prism_stacked <- raster::aggregate(prism_stacked, fact=9)
writeRaster(prism_stacked,file="prism_aggregated.grd")
############################

#############################################################################################

# Load species data
bbs <- get_bbs_data()

# Get locations of surveys
locations <- dplyr::select(bbs, site_id, long, lat) %>% distinct()
coordinates(locations) <- c("long", "lat")
# set projection of the PRISM data
raster::crs(locations) <- '+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0'

# Load PRISM (past climate) data
prism_stacked <- brick("prism_aggregated.grd")

# Extract climate data for BBS locations
extracted <- raster::extract(prism_stacked, locations)
prism_bbs_data <- data.frame(site_id = locations$site_id, coordinates(locations), extracted)
prism_bbs_data <- prism_bbs_data %>%
  gather(date, value, 4:ncol(prism_bbs_data)) %>%
  tidyr::extract(date, c("clim_var", "year", "month"),
                 "PRISM_([:alpha:]*)_stable_[:alnum:]*_([:digit:]{4})([:digit:]{2})_")
#Format the data a little and load into the sqlite database.
prism_bbs_data$year <- as.numeric(prism_bbs_data$year)
prism_bbs_data$month <- as.numeric(prism_bbs_data$month)

# Calculate BIOCLIM variables
bioclim <- process_bioclim_data(prism_bbs_data)

# Match BBS and climate data
bbs.env <- inner_join(bbs, bioclim, by = c('site_id', 'year')) %>% data.frame()
write.table(bbs.env, file="BBS_BIOCLIM.txt", sep="\t", row.names=FALSE)

















