###########################################
# Supporting Information
# "Increasing climatic decoupling of bird abundances and distributions"
# Duarte S. Viana, Jonathan M. Chase
###########################################

# Data preparation

library(tidyr)
library(sp)
library(rgdal)
library(raster)
library(rgeos)
library(geosphere)
library("mapdata")

# Load BBS data
#bbs.clim <- read.table("BBS_BIOCLIM.txt",sep="\t",header=T)
#save(bbs.clim,file="BBS_BIOCLIM.Rdata")
load("BBS_BIOCLIM.Rdata")
bbs.wide <- spread(bbs.clim, key=species_id, value=abundance)

# Site filtering
ss <- list()
for(i in 1989:2018) ss[[i-1988]] <- unique(bbs.clim$site_id[bbs.clim$year==i])
common <- Reduce(intersect, ss)
clim <- bbs.wide[bbs.wide$site_id %in% common, c("site_id","bio3","bio10","bio15","bio16")]
data.in <- unique(na.exclude(clim)$site_id)
sites.in <- sites.sp[sites$site_id %in% data.in,]

# Map of included sites
quartz(height=3,width=5)
par(mar=c(0,0,0,0))
map('usa')
points(sites.in, col='blue', pch=16, cex=0.7)
length(sites.in) # 179 sites


#######################################################################
# Data selection for subsequent analyses

# Explore different time extents
time.extent <- c(30,25,20) # years

# Common sites across years for each time extent
idata <- list()
ss <- list()
for(i in 1989:2018) ss[[i-1988]] <- unique(bbs.clim$site_id[bbs.clim$year==i])
common <- Reduce(intersect, ss)
clim <- bbs.wide[bbs.wide$site_id %in% common, c("site_id","bio3","bio10","bio15","bio16")]
idata[[1]] <- unique(na.exclude(clim)$site_id)
ss <- list()
for(i in 1994:2018) ss[[i-1993]] <- unique(bbs.clim$site_id[bbs.clim$year==i])
common <- Reduce(intersect, ss)
clim <- bbs.wide[bbs.wide$site_id %in% common, c("site_id","bio3","bio10","bio15","bio16")]
idata[[2]] <- unique(na.exclude(clim)$site_id)
ss <- list()
for(i in 1999:2018) ss[[i-1998]] <- unique(bbs.clim$site_id[bbs.clim$year==i])
common <- Reduce(intersect, ss)
clim <- bbs.wide[bbs.wide$site_id %in% common, c("site_id","bio3","bio10","bio15","bio16")]
idata[[3]] <- unique(na.exclude(clim)$site_id)

# All the points within inner.radius of a center point will be in the test set.
# Everything more than outer.radius away will be in the training set.
# radii are in meters
# Adapted from Harris 2015
l.centers <- list()
jk <- 0
for(j in 1:5){
  centers <-  regularCoordinates(j+14)
  for(k in 1:10){
    jk <- jk+1
    l.centers[[jk]] <- jitter(centers, 70)
  }
}
radius <-  1.8e5

itrain <- list()
itest <- list()
id.sp <- 0

for(s in 1:3){
  if(s==1) years <- 1989:2018
  if(s==2) years <- 1994:2018
  if(s==3) years <- 1999:2018
  sites.in <- idata[[s]]
  bbs.s <- bbs.wide[bbs.wide$site_id %in% sites.in,]
  lonlat <- bbs.s[!duplicated(bbs.s$site_id),c("long","lat")]
  # split training and test data (for 50-fold CV)
  in.train <- list()
  in.test <- list()
  for(j in 1:50){
    centers <-  l.centers[[j]]
    dists <-  pointDistance(centers, lonlat, longlat = TRUE)
    in.train[[j]] <- sites.in[apply(dists, 2,function(x) min(x) > radius)]
    in.test[[j]] <- sites.in[apply(dists, 2, function(x) min(x) <= radius)]
  }
  itrain[[s]] <- in.train
  itest[[s]] <- in.test
  
  for(y in years){
    # All sites surveyed in year y
    bbs.y <- bbs.s[bbs.s$year==y,]
    bbs.y <- bbs.y[!duplicated(bbs.y$site_id),]
    
    # species information
    for(i in 1:424){
      id.sp <- id.sp+1
      bbs.sp.i <- cbind(bbs.y[,c("long","lat","site_id")], abund=bbs.y[,23+i])
      bbs.sp.i$abund[is.na(bbs.sp.i$abund)] <- 0
      Npres.train.j <- c()
      Npres.test.j <- c()
      for(j in 1:5){
        Npres.train.j[j] <- nrow(bbs.sp.i[bbs.sp.i$site_id %in% in.train[[j]] & bbs.sp.i$abund>0,])
        Npres.test.j[j] <- nrow(bbs.sp.i[bbs.sp.i$site_id %in% in.test[[j]] & bbs.sp.i$abund>0,])
      }
      Npres.train <- min(Npres.train.j)
      Npres.test <- min(Npres.test.j)
      
      bbs.sp.i <- bbs.sp.i[bbs.sp.i$abund>0,]
      Npres <- nrow(bbs.sp.i)
      if(Npres<3) data.id[id.sp,] <- c(s, y, names(bbs.y)[23+i], Npres, 
                                       Npres.train, Npres.test, NA,NA,NA,NA,NA,NA)
      if(Npres>=3){
        # calculate occupied area ("range size")
        sites.i <-SpatialPointsDataFrame(bbs.sp.i[,c("long","lat")],
                                         data=data.frame(bbs.sp.i[,"site_id"]),
                                         proj4string=CRS("+proj=longlat +datum=WGS84"))
        hull <- gConvexHull(sites.i)
        hullt <- spTransform(hull,
                             CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"))
        iarea <- gArea(hullt)/1e6 # km2
        # calculate the abundance-weighted centroid (e.g. Ash et al. 2017)
        w.centre.x <- sum(bbs.sp.i$long*bbs.sp.i$abund)/sum(bbs.sp.i$abund)
        w.centre.y <- sum(bbs.sp.i$lat*bbs.sp.i$abund)/sum(bbs.sp.i$abund)
        # abundance metrics
        abund.mean <- mean(bbs.sp.i$abund)
        abund.cv <- sd(bbs.sp.i$abund)/mean(bbs.sp.i$abund)
        abund.total <- sum(bbs.sp.i$abund)
        
        data.id[id.sp,] <- c(s, y, names(bbs.y)[23+i], Npres, Npres.train, Npres.test, 
                             abund.mean, abund.cv, abund.total, 
                             iarea, w.centre.x, w.centre.y)
      }
    }
  }
}


for(i in 1:ncol(data.id)) data.id[,i] <- as.numeric(data.id[,i])

# cases for climatic models
cases <- data.frame(time.extent=rep(c(rep(1,30),rep(2,25),rep(3,20))), year=c(1989:2018,1994:2018,1999:2018))
# Replicate dataset to apply other methods (BRT and GAM)
ni <- nrow(cases)
cases <- cases[rep(1:ni,2),]
cases$data.type <- rep(c("abund","pres"),each=ni)
ni <- nrow(cases)
cases <- cases[rep(1:ni,4),]
cases$method <- rep(c("BRT_lr.01","BRT_lr.001","GAM_kest","GAM_k3"),each=ni)

save(idata,itrain,itest,data.id,cases,file="Data_id.Rdata")


################################################################################
################################################################################

# Get species names and AOU codes

library(DBI)
#library(RSQLite)
library(dplyr)
#library(dbdplyr)

setwd("~/Documents/Papers/iDiv/Paper_sDiv/Data/BBS")

bbs_db <- dbConnect(RSQLite::SQLite(), 'bbs_sqlite.db')
sps <- tbl(bbs_db, "breed_bird_survey_species") %>% data.frame()

sci.names <- paste(sps$genus,sps$species,sep="_")
sps$sci <- sci.names

# Add BirdTree names
birdtree.sp <- read.csv("BLIOCPhyloMasterTax.csv")
sps <- left_join(sps, birdtree.sp[,c("TipLabel","English")], by=c("english_common_name" = "English"))

setwd("~/Documents/Papers/iDiv/Paper_sDiv/BBS_analysis")
save(sps,file="species_names.Rdata")



