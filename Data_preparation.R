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

# Figure S4: Map of included sites
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
clim <- bbs.wide[bbs.wide$site_id %in% common, c("site_id","bio1")]
idata[[1]] <- unique(na.exclude(clim)$site_id)
ss <- list()
for(i in 1994:2018) ss[[i-1993]] <- unique(bbs.clim$site_id[bbs.clim$year==i])
common <- Reduce(intersect, ss)
clim <- bbs.wide[bbs.wide$site_id %in% common, c("site_id","bio1")]
idata[[2]] <- unique(na.exclude(clim)$site_id)
ss <- list()
for(i in 1999:2018) ss[[i-1998]] <- unique(bbs.clim$site_id[bbs.clim$year==i])
common <- Reduce(intersect, ss)
clim <- bbs.wide[bbs.wide$site_id %in% common, c("site_id","bio1")]
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

data.id <- data.frame(time.extent=numeric(0), year=numeric(0), aou=numeric(0), Npres=numeric(0), Npres.train=numeric(0),
                      Npres.test=numeric(0), abund.mean=numeric(0), abund.cv=numeric(0), abund.total=numeric(0))
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
                                       Npres.train, Npres.test, NA,NA,NA)
      if(Npres>=3){
        # abundance metrics
        abund.mean <- mean(bbs.sp.i$abund)
        abund.cv <- sd(bbs.sp.i$abund)/mean(bbs.sp.i$abund)
        abund.total <- sum(bbs.sp.i$abund)
        data.id[id.sp,] <- c(s, y, names(bbs.y)[23+i], Npres, Npres.train, Npres.test, 
                             abund.mean, abund.cv, abund.total)
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


################################################################################
################################################################################

# Descriptive patterns and phylogenetic relatedness

setwd("~/Documents/Papers/iDiv/Paper_BBS_R2-time/EVE")
load("results_v3.Rdata")
load("Data_id_v2.Rdata")

library(dplyr)

# Define the the different approaches
time.extent <- 1:3
data.type <- c("abund","pres")
method <- c("BRT_lr.01","BRT_lr.001","GAM_kest","GAM_k3")
r2.type <- 3:4
meth <- expand.grid(time.extent=time.extent,data.type=data.type,method=method,r2.type=r2.type)


# Responsiveness data (species-climate r2)
l.sp <- list()
for(task.id in 1:nrow(meth)){
  r2s <- as.data.frame(res[which(cases$time.extent==meth$time.extent[task.id] & 
                                   cases$data.type==meth$data.type[task.id] & 
                                   cases$method==meth$method[task.id]),,meth$r2.type[task.id]])
  id.sp <- unique(data.id$aou[data.id$time.extent==meth$time.extent[task.id] & 
                                data.id$Npres.train>=5 & data.id$Npres.test>=5])
  id <- data.id[data.id$time.extent==meth$time.extent[task.id] & data.id$aou %in% id.sp,]
  
  # Prepare data
  lres <- id
  lres$aou <- as.character(lres$aou)
  lres$r2 <- numeric(nrow(lres))
  for(i in unique(lres$aou)) lres$r2[lres$aou==i] <- r2s[,names(r2s)==i]
  lres$aou <- as.character(lres$aou)
  lresi <- na.exclude(lres)
  names(lresi)[2] <- "time"
  N <- tapply(lresi$time,lresi$aou,function(x) length(unique(x)))
  dat.in <- names(N[N>=10]) # at least 10 years
  lresi <- lresi[lresi$aou %in% dat.in,]
  l.sp[[task.id]] <- unique(lresi$aou)
}
lapply(l.sp[meth$time.extent==1],length)
allSame <- function(x) length(unique(x)) == 1
allSame(l.sp[meth$time.extent==1]) # species are all the same across species-climate models

spdf <- data.frame(aou=l.sp[[1]])

# community-level temporal trend in climate matching
years <- 1989:2018
res1 <- as.data.frame(res[which(cases$time.extent==1 & cases$data.type=="abund" & cases$method=="BRT_lr.01"),,4])
res1$year <- years
res2 <- as.data.frame(res[which(cases$time.extent==1 & cases$data.type=="abund" & cases$method=="BRT_lr.001"),,4])
res2$year <- years
res3 <- as.data.frame(res[which(cases$time.extent==1 & cases$data.type=="abund" & cases$method=="GAM_kest"),,4])
res3$year <- years
res4 <- as.data.frame(res[which(cases$time.extent==1 & cases$data.type=="abund" & cases$method=="GAM_k3"),,4])
res4$year <- years

res1g <- gather(res1,key="aou",value="r2",-year)
res2g <- gather(res2,key="aou",value="r2",-year)
res3g <- gather(res3,key="aou",value="r2",-year)
res4g <- gather(res4,key="aou",value="r2",-year)
resg <- rbind(res1g,res2g,res3g,res4g)
resg$method <- rep(c("BRT_lr.01","BRT_lr.001","GAM_kest","GAM_k3"),each=nrow(res1g))
setwd("~/Documents/Papers/iDiv/Paper_BBS_R2-time")
save(resg,file="R2_sp_year.Rdata")


# Mean climate matching for each species
r2.m1 <- apply(res1[,names(res1) %in% spdf$aou],2,mean,na.rm=T)
r2.m2 <- apply(res2[,names(res2) %in% spdf$aou],2,mean,na.rm=T)
r2.m3 <- apply(res3[,names(res3) %in% spdf$aou],2,mean,na.rm=T)
r2.m4 <- apply(res4[,names(res4) %in% spdf$aou],2,mean,na.rm=T)
r2.pooled <- apply(rbind(r2.m1,r2.m2,r2.m3,r2.m4),2,mean,na.rm=T)
r2.pooled <- data.frame(aou=names(r2.pooled),baseline.r2=r2.pooled)
spdf <- left_join(spdf,r2.pooled,by="aou")

#--------------------------------------------------------------

# Phylogenetic relatedness

load("species_names.Rdata")
sps$aou <- as.character(sps$aou)
spdf <- left_join(spdf,sps[,c("aou","sci","sporder","family")],by="aou")

sp.names <- spdf$sci
birdtree.sp <- read.csv("BLIOCPhyloMasterTax.csv")
all(sp.names %in% birdtree.sp$TipLabel)

# Correct unmatched species names
mis.sp <- data.frame(sci=sp.names[!(sp.names %in% birdtree.sp$TipLabel)])
mis.sp <- left_join(mis.sp, sps[,c("TipLabel","sci")], by="sci")
mis.sp$TipLabel <- as.character(mis.sp$TipLabel)
mis.sp$TipLabel[mis.sp$sci=="Setophaga_coronata"] <- "Dendroica_coronata"
spdf$phylo <- spdf$sci
for(i in mis.sp$sci) spdf$phylo[spdf$sci==i] <- mis.sp$TipLabel[mis.sp$sci==i]

# From https://github.com/nicholasjclark/BBS.occurrences/blob/master/Clark_etal_analysis/Appendix_S4_PhyloTraitData.Rmd
# Get species names (for the filtered data) and remove underscore
sp.names <- gsub('_', ' ', spdf$phylo)

#Export the vector as a one-column .csv file 
write.csv(sp.names, "sp.names.csv", row.names = FALSE)

# Once exported, copy the species names and paste into the `Select species` form 
# on the [Phylogeny Subsets](http://birdtree.org/subsets/) page at Birdtree.org. 
# We downloaded 100 trees from the *Ericsson All Species* dataset for our analyses. 
# Once processed, save the resulting .nex file and read in the multiphylo object 
# using functions in the `ape` package
setwd("/Users/Viana/Documents/Papers/iDiv/Paper_BBS_R2-time/Phylo_100_BirdTree")
sp.trees <- ape::read.nexus("output.nex")
# Now check to make sure that all of the species names are represented in the tree. 
# Here, we have to include the underscore once again, as this is included in Birdtree.org phylogeny subsets. 
# This call should return `TRUE` if there are no unmatched names
sp.names.underscore <- gsub(' ', '_', sp.names)
all(sp.names.underscore %in% sp.trees[[1]]$tip.label)

# From https://cran.r-project.org/web/packages/brms/vignettes/brms_phylogenetics.html
# The phylo object contains information on the relationship between species. 
# Using this information, we can construct a covariance matrix of species (Hadfield & Nakagawa, 2010).
# Build a consensus tree
library(phytools)
cons.tree <- consensus.edges(sp.trees,method="mean.edge")
cons.tree <- di2multi(cons.tree)
A <- ape::vcv.phylo(cons.tree)

save(spdf,A,file="phylo_vcv.Rdata")

