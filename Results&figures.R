###########################################
# Supporting Information
# "Increasing climatic decoupling of bird abundances and distributions"
# Duarte S. Viana, Jonathan M. Chase
###########################################

# Results and figures

library(dplyr)
library(tidyr)
library(brms)
#options(buildtools.check = function(action) TRUE )
library(tidybayes)
library(stringr)
library(ggplot2)
library(cowplot)
library(egg)
theme_set(theme_tidybayes() + panel_border())
library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(VoCC)

# Load data
load("Data_id_v2.Rdata")
load("phylo_vcv.Rdata")
load("specialization_index.Rdata")
iucn <- read.csv("IUCN_status.csv")
iucn$sci <- paste(iucn$genusName,iucn$speciesName,sep="_")
traits <- read.table("traits.txt",header=T,sep="\t")
names(traits) <- c("species","lifespan","log.mass","body.length","wingspan",names(traits[,6:ncol(traits)]))
acad <- read.csv("ACAD Global 2021.02.05-filtered.csv")
mig <- read.table("specieslist3_1_migbehav_v1_0.txt",header=T,sep="\t")
hwi <- read.table("Traits_Sheard_HWI.txt",header=T,sep="\t")

# Traits
traits <- spdf

# Add habitat specialization index (from Martin and Fahrig 2018)
traits1.1 <- left_join(traits,ssi[,c("sci","SSI")],by=c("phylo"="sci"))
traits1.2 <- left_join(traits,ssi[,c("sci","SSI")],by="sci")
traits1.1[is.na(traits1.1$SSI),] <- traits1.2[is.na(traits1.1$SSI),]
traits <- as.data.frame(traits1.1)
traits$SSI <- as.numeric(traits$SSI)
a <- which(is.na(traits$SSI))
unique(traits$sci[a])

# from Sheard et al. 2020; https://doi.org/10.1038/s41467-020-16313-6 
hwi$phylo <- gsub(' ', '_', hwi$"Tree.name")
traits <- left_join(traits,hwi[,c("phylo","HWI","Body.mass..log.","Territoriality","Habitat")],by="phylo")
a <- which(is.na(traits$HWI))
unique(traits$phylo[a])
traits$main.habitat <- factor(traits$Habitat)
levels(traits$main.habitat) <- c("dense","semi-open","open")

# From Eyres et al. 2017;  https://doi.org/10.1111/jav.01308
traits <- left_join(traits,mig[,c("IOC3_1_Binomial","Migratory_status")],by=c("sci"="IOC3_1_Binomial"))
a <- which(is.na(traits$Migratory_status))
unique(traits$phylo[a])
traits[traits$phylo %in% unique(traits$phylo[a]),"Migratory_status"] <- mig[mig$IOC3_1_Binomial %in% unique(traits$phylo[a]),"Migratory_status"]
a <- which(is.na(traits$Migratory_status))
unique(traits$phylo[a])

# From ACAD; https://pif.birdconservancy.org/avian-conservation-assessment-database/
acad$sci <- gsub(' ', '_', acad$Scientific.Name)
traits <- left_join(traits,acad[,c("sci","PS.g","BD.g","ND.g","TB.c","TN.c","PT.c","CCS.max","Half.Life")],by="sci")
a <- which(is.na(traits$TB.c))
unique(traits$phylo[a])


# Define the the different approaches
time.extent <- 1:3
data.type <- c("abund","pres")
method <- c("BRT_lr.01","BRT_lr.001","GAM_kest","GAM_k3")
r2.type <- 3:4
#clim.set <- c("clim1","clim2")
meth <- expand.grid(time.extent=time.extent,data.type=data.type,method=method,r2.type=r2.type)
# Define the different species responses
resp <- c("abund.mean","abund.cv", "occupancy")

#-------------------------

# Pool posteriors of the time models of the different approaches

# First, print results of the trends of climate matching estimated with BRT_lr.01 model (model 25)
load("bbs-8862844-25.Rdata")
summary(mm[[1]])

wcc <- which(meth$time.extent==1 & meth$data.type=="abund" & meth$r2.type==4)
res.list <- list()
pop.eff <- list()
for(i in 1:length(wcc)){
  load(paste("bbs-8862844-", wcc[i],".Rdata", sep = ""))
  # Manual extraction of coefficients for r2
  m.time <- mm[[1]]
  t.fe <- as.data.frame(fixef(m.time,summary=FALSE))
  t.fe <- data.frame(b=t.fe$time)
  pop.eff[[i]] <- t.fe$b
  ref <- ranef(m.time,summary=FALSE)
  # species
  t.sp.re <- gather(as.data.frame(ref$species[,,"time"]))
  re.sp <- cbind(t.sp.re,b_time=t.fe[rep(1:nrow(t.fe),length(names(ref$species[1,,2]))),])
  names(re.sp) <- c("species","r_sp","b_year")
  re.sp$coef <- re.sp$b_year + re.sp$r_sp
  res.m <- re.sp[,c("species","coef")]
  names(res.m)[2] <- paste(names(res.m)[2], ".r2", sep="")
  
  for(r in 1:length(resp)){
    m.time <- mm[[r+1]]
    # Manual extraction of coefficients
    t.fe <- as.data.frame(fixef(m.time,summary=FALSE))
    t.fe <- data.frame(b=t.fe$time)
    ref <- ranef(m.time,summary=FALSE)
    # species
    t.sp.re <- gather(as.data.frame(ref$species[,,"time"]))
    re.sp <- cbind(t.sp.re,b_time=t.fe[rep(1:nrow(t.fe),length(names(ref$species[1,,2]))),])
    names(re.sp) <- c("species","r_sp","b_year")
    re.sp$coef <- re.sp$b_year + re.sp$r_sp
    res.m <- cbind(res.m, re.sp[,4])
    names(res.m)[2+r] <- paste("coef.", resp[r], sep="")
    res.list[[i]] <- res.m
  }
}

#save(res.list,pop.eff,file="res_list.Rdata")
#load("res_list.Rdata")

# Combine posteriors
re.sp <- do.call("rbind",res.list)

# Take means and 90% confidence intervals
re.sp.sum <- re.sp %>%
  group_by(species) %>%
  summarize(mean.r2 = mean(coef.r2),lower.r2=quantile(coef.r2,0.05),upper.r2=quantile(coef.r2,0.95),
            mean.abund = mean(coef.abund.mean),lower.abund=quantile(coef.abund.mean,0.05),upper.abund=quantile(coef.abund.mean,0.95),
            mean.abund.cv = mean(coef.abund.cv),lower.abund.cv=quantile(coef.abund.cv,0.05),upper.abund.cv=quantile(coef.abund.cv,0.95),
            mean.occupancy = mean(coef.occupancy),lower.occupancy=quantile(coef.occupancy,0.05),upper.occupancy=quantile(coef.occupancy,0.95))

# Add traits
re.sp.sum <- left_join(re.sp.sum,traits,by=c("species"="aou"))

#-------------------------

# Obtain rates of climate change in the USA

load("BBS_BIOCLIM.Rdata")
bbs.wide <- spread(bbs.clim, key=species_id, value=abundance)

library(prism)
source("forecast-bbs-core.R")
source("get_prism_data.R")
options(prism.path = "/Users/Viana/Documents/Papers/iDiv/Paper_sDiv/Data/BBS/prismdata_annual")
clim_vars <- c("ppt", "tmin", "tmean", "tmax")

# years_to_use=1981:2018
# for (clim_var in clim_vars){
#   get_prism_annual(type=clim_var, year = years_to_use, keepZip=F)
# }
prism_stacked <- pd_stack(prism_archive_ls())

# BBS sites included
sites <- bbs.clim[!duplicated(bbs.clim$site_id),c("long","lat","site_id")]
row.names(sites) <- sites$site_id
sites.sp<-SpatialPointsDataFrame(sites[,1:2],data=data.frame(sites$site_id),
                                 proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
sites.in <- sites.sp[sites$site_id %in% idata[[1]],]
sites.in <- spTransform(sites.in,CRS("+proj=utm +zone=15 +datum=WGS84"))

sites.buf <- gBuffer(sites.in,byid=TRUE,width=50000) # 50km buffer
sites.buf<-spTransform(sites.buf,CRS(proj4string(prism_stacked)))
sites.in$ppt <- NA
sites.in$tmin <- NA
sites.in$tmean <- NA
sites.in$tmax <- NA
for(i in 1:length(sites.buf)){
  raster.i <- crop(prism_stacked,sites.buf[i,])
  for(j in clim_vars){
    raster.ij <- raster::subset(raster.i, grep(j, names(raster.i), value = T))
    tr <- tempTrend(raster.ij, th = 30)
    tr1 <- raster::subset(tr,1)
    sites.in[i,j] <- cellStats(tr1, "mean")
  }
}

re.sp.sum$cc.tmean <- NA
re.sp.sum$cc.ppt <- NA
for(i in re.sp.sum$species){
  spi <- bbs.wide[bbs.wide$site_id %in% idata[[1]],c("site_id",i)]
  spi <- na.exclude(spi)
  spi <- spi[!duplicated(spi$site_id),]
  re.sp.sum$cc.tmean[re.sp.sum$species==i] <- mean(sites.in$tmean[sites.in$sites.site_id %in% spi$site_id])
  re.sp.sum$cc.ppt[re.sp.sum$species==i] <- mean(sites.in$ppt[sites.in$sites.site_id %in% spi$site_id])
}

plot(re.sp.sum$cc.tmean,re.sp.sum$mean.r2)
cor.test(re.sp.sum$cc.tmean,re.sp.sum$mean.r2)
plot(re.sp.sum$cc.ppt,re.sp.sum$mean.r2)
cor.test(re.sp.sum$cc.ppt,re.sp.sum$mean.r2)

#-------------------------

# Obtain rates of land use change in the USA

# load NLCD data (1992)
# nlcd92.1 <- raster("NLCD_USGS/nlcde92_1/nlcde1.tif")
# nlcd92.2 <- raster("NLCD_USGS/nlcde92_2/nlcde2.tif")
# nlcd92.3 <- raster("NLCD_USGS/nlcde92_3/nlcde3.tif")
# nlcd92.4 <- raster("NLCD_USGS/nlcde92_4/nlcde4.tif")
# nlcd92.12 <- merge(nlcd92.1, nlcd92.2)
# nlcd92.34 <- merge(nlcd92.3, nlcd92.4)
# nlcd92 <- merge(nlcd92.12, nlcd92.34)
# plot(nlcd92)
# writeRaster(nlcd92, "NLCD_USGS/nlcd92_merged/nlcd92_merged.tif", format="GTiff") 
nlcd92 <- raster("NLCD_USGS/nlcd92_merged/nlcd92_merged.tif")
# load NLCD data (2004, 2016)
nlcd04 <- stack("NLCD_USGS/NLCD_2004_Land_Cover_L48_20190424/NLCD_2004_Land_Cover_L48_20190424.img")
nlcd16 <- stack("NLCD_USGS/NLCD_2016_Land_Cover_L48_20190424/NLCD_2016_Land_Cover_L48_20190424.img")
# GRS80 is roughly WGS84 for resolutions >1-2m
proj4string(nlcd92) <- proj4string(nlcd04)

# BBS sites included
sites.in2 <- spTransform(sites.in,CRS(proj4string(nlcd04)))

lcc.fun <- function(r1, r2) floor(r2/10)-floor(r1/10)
sites.buf <- gBuffer(sites.in2,byid=TRUE,width=50000)
sites.in$lcc92_04 <- NA
sites.in$lcc04_16 <- NA
sites.in$lcc92_16 <- NA
for(i in 1:length(sites.buf)){
  nlcd92.i <- crop(nlcd92,sites.buf[i,])
  nlcd04.i <- crop(nlcd04,sites.buf[i,])
  nlcd16.i <- crop(nlcd16,sites.buf[i,])
  lcc92_04 <- overlay(nlcd92.i, nlcd04.i, fun=lcc.fun)
  lcc04_16 <- overlay(nlcd04.i, nlcd16.i, fun=lcc.fun)
  lcc92_16 <- overlay(nlcd92.i, nlcd16.i, fun=lcc.fun)
  sites.in$lcc92_04[i] <- 1-(freq(lcc92_04,value=0)/cellStats(lcc92_04, function(i, ...) sum(!is.na(i))))
  sites.in$lcc04_16[i] <- 1-(freq(lcc04_16,value=0)/cellStats(lcc04_16, function(i, ...) sum(!is.na(i))))
  sites.in$lcc92_16[i] <- 1-(freq(lcc92_16,value=0)/cellStats(lcc92_16, function(i, ...) sum(!is.na(i))))
}

re.sp.sum$lcc.mean <- NA
for(i in re.sp.sum$species){
  spi <- bbs.wide[bbs.wide$site_id %in% idata[[1]],c("site_id",i)]
  spi <- na.exclude(spi)
  spi <- spi[!duplicated(spi$site_id),]
  lcc.i <- sites.in$lcc92_16[sites.in$sites.site_id %in% spi$site_id]
  re.sp.sum$lcc.mean[re.sp.sum$species==i] <- mean(lcc.i)
}

plot(re.sp.sum$lcc.mean,re.sp.sum$mean.r2)
cor.test(re.sp.sum$lcc.mean,re.sp.sum$mean.r2)

# load HILDA data
hilda <- raster("hildap_vGLOB-1.0_change-layers/HILDAplus_vGLOB-1.0_luc_change-freq_1960-2019_wgs84.tif")
#hilda.forest <- raster("hildap_vGLOB-1.0_change-layers/HILDAplus_vGLOB-1.0_luc_forest_change_1960-2019_wgs84.tif")
#hilda.open <- raster("hildap_vGLOB-1.0_change-layers/HILDAplus_vGLOB-1.0_luc_pasture-rangeland_change_1960-2019_wgs84.tif")
#hilda.cropland <- raster("hildap_vGLOB-1.0_change-layers/HILDAplus_vGLOB-1.0_luc_cropland_change_1960-2019_wgs84.tif")

# BBS sites included
sites.in3 <- sites.sp[sites$site_id %in% idata[[1]],]

#extract values to points
lcc <- raster::extract(hilda, sites.in3, buffer=50000, fun=mean) # 50km buffer
sites.in$lcc <- lcc

lcc.mean.sites <- mean(sites.in$lcc)
re.sp.sum$lcc.hilda <- NA
for(i in re.sp.sum$species){
  spi <- bbs.wide[bbs.wide$site_id %in% idata[[1]],c("site_id",i)]
  spi <- na.exclude(spi)
  spi <- spi[!duplicated(spi$site_id),]
  lcc.i <- sites.in$lcc[sites.in$sites.site_id %in% spi$site_id]
  re.sp.sum$lcc.hilda[re.sp.sum$species==i] <- mean(lcc.i)
}

plot(re.sp.sum$lcc.hilda,re.sp.sum$mean.r2)
cor.test(re.sp.sum$lcc.hilda,re.sp.sum$mean.r2)

save(re.sp.sum,sites.in,file="re.sp.sum.Rdata")

#-------------------------
# Results and figures

load("res_list.Rdata")
load("re.sp.sum.Rdata")
re.sp.sum <- as.data.frame(re.sp.sum)

# Figure 2: R2 temporal trend according to the different methods

load("R2_sp_year.Rdata")

gg.trends <- resg %>% 
  group_by(method,year) %>%
  summarise(r2 = mean(r2, na.rm=T))  %>%
  ggplot(aes(x=year, y=r2, colour=method)) + 
  geom_point()+
  geom_smooth(method=lm) +
  labs(x="Year", y="Climate matching",colour="") +
  theme(legend.direction="horizontal", legend.position="top",
        legend.title = element_text(size = 10),legend.text = element_text(size = 8))


# Population-level effect
wcc <- which(meth$time.extent==1 & meth$data.type=="abund" & meth$r2.type==4)
pop.est <- do.call("c", pop.eff)
mean(pop.est)
quantile(pop.est,c(0.05,0.95))
lapply(pop.eff,length)
pop.r2 <- data.frame(method=rep(meth$method[wcc],each=8000), b=pop.est)

gg.posterior <- pop.r2 %>%
  mutate(grid_mean = b) %>%
  ggplot(aes(y = method, x = grid_mean)) + # , fill = stat(abs(x) < .9)
  #geom_hline(yintercept = levels(factor(re.sp$sporder)), col="grey") +
  stat_halfeyeh(.width = c(0.50, 0.90)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x =expression(paste(Delta,"climate matching")), y = "Species-climate model") +
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))

# Compose figure 
quartz(height=4,width=7)
ggarrange(gg.trends, gg.posterior, nrow=1)

#-------------------------

# Figure 3: habitat specialiasation

cor.test(re.sp.sum$SSI,re.sp.sum$mean.r2)

gg.time.bm <- re.sp.sum %>%
  ggplot(aes(x=Body.mass..log., y=mean.r2)) + 
  geom_point(col="black", size=1) +
  geom_smooth(method="lm", col="black") +
  labs(x = "Body mass (log)", y = expression(paste(Delta,"climate matching"))) 

gg.time.hs <- re.sp.sum %>%
  ggplot(aes(x=SSI, y=mean.r2)) + 
  geom_point(col="black", size=1) +
  geom_smooth(method="lm", col="black") +
  labs(x = "Habitat specialisation", y = expression(paste(Delta,"climate matching")))

# Compose figure 
quartz(height=3,width=7)
ggarrange(gg.time.bm, gg.time.hs, nrow=1)


#-------------------------

# Figure 4: joint changes in R2 and performance

sig.r2 <- re.sp.sum$lower.r2<=0 & re.sp.sum$upper.r2>=0
sig.abund <- re.sp.sum$lower.abund<=0 & re.sp.sum$upper.abund>=0
sig.r2.abund <- sig.r2==FALSE & sig.abund==FALSE
gg.abund <- re.sp.sum %>%
  ggplot(aes(x=mean.r2, y=mean.abund)) + 
  geom_point(col="darkgrey", size=1) +
  geom_errorbar(aes(ymin=lower.abund, ymax=upper.abund), col="darkgrey",
                position=position_dodge(0.05)) +
  geom_errorbar(aes(xmin=lower.r2, xmax=upper.r2), col="darkgrey",
                position=position_dodge(0.05)) +
  geom_point(data=re.sp.sum[!sig.r2,], aes(x=mean.r2, y=mean.abund), col="red", size=1) +
  geom_point(data=re.sp.sum[!sig.abund,], aes(x=mean.r2, y=mean.abund), col="blue", size=1) +
  geom_point(data=re.sp.sum[sig.r2.abund,], aes(x=mean.r2, y=mean.abund), col="purple", size=1) +
  geom_hline(yintercept = 0, linetype = "dashed", size=1, col="blue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", size=1, col="red", alpha = 0.5) +
  #geom_smooth(method="lm", col="black") +
  labs(x = expression(paste(Delta,"climate matching")), y = expression(paste(Delta,"abundance")))

sig.occup <- re.sp.sum$lower.occupancy<=0 & re.sp.sum$upper.occupancy>=0
sig.r2.occup <- sig.r2==FALSE & sig.occup==FALSE
gg.occupancy <- re.sp.sum %>%
  ggplot(aes(x=mean.r2, y=mean.occupancy)) + 
  geom_point(col="darkgrey", size=1) +
  geom_errorbar(aes(ymin=lower.occupancy, ymax=upper.occupancy), col="darkgrey",
                position=position_dodge(0.05)) +
  geom_errorbar(aes(xmin=lower.r2, xmax=upper.r2), col="darkgrey",
                position=position_dodge(0.05)) +
  geom_point(data=re.sp.sum[!sig.r2,], aes(x=mean.r2, y=mean.occupancy), col="red", size=1) +
  geom_point(data=re.sp.sum[!sig.occup,], aes(x=mean.r2, y=mean.occupancy), col="blue", size=1) +
  geom_point(data=re.sp.sum[sig.r2.occup,], aes(x=mean.r2, y=mean.occupancy), col="purple", size=1) +
  geom_hline(yintercept = 0, linetype = "dashed", size=1, col="blue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", size=1, col="red", alpha = 0.5) +
  #geom_smooth(method="lm", col="black") +
  labs(x = expression(paste(Delta,"climate matching")), y = expression(paste(Delta,"occupancy")))


give.n <- function(x){
  return(c(y = 0.315, label = length(x)))
}

re.sp.sum$PT.c <- factor(re.sp.sum$PT.c)
gg.time.poptrend <- re.sp.sum %>%
  ggplot(aes(x=PT.c, y=mean.r2)) + 
  geom_boxplot() +
  ylim(-0.35,0.32) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x ="Population trend", y = expression(paste(Delta,"climate matching"))) +
  coord_flip() +
  stat_summary(fun.data = give.n, geom = "text", size = 3)

re.sp.sum$TB.c <- factor(re.sp.sum$TB.c)
gg.time.threat <- re.sp.sum %>%
  ggplot(aes(x=TB.c, y=mean.r2)) + 
  geom_boxplot() +
  ylim(-0.35,0.32) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x ="Threat status", y = expression(paste(Delta,"climate matching"))) +
  coord_flip() +
  stat_summary(fun.data = give.n, geom = "text", size = 3)

quartz(height=5.5,width=6)
ggarrange(gg.abund, gg.occupancy, gg.time.poptrend, gg.time.threat, nrow=2)


# Correlation between trends in performance and responsiveness to climate
cor.test(re.sp.sum$mean.r2,re.sp.sum$mean.abund)
cor.test(re.sp.sum$mean.r2,re.sp.sum$mean.occupancy)

#-------------------------

# Extract more results

re.sp.sum$sig.r2 <- sig.r2
re.sp.sum$sig.abund <- sig.abund
re.sp.sum$sig.occup <- sig.occup
re.sp.sum$trend <- NA
re.sp.sum$trend[re.sp.sum$mean.r2<0] <- "negative"
re.sp.sum$trend[re.sp.sum$mean.r2>0] <- "positive"
re.sp.sum$trend[re.sp.sum$lower.r2<=0 & re.sp.sum$upper.r2>=0] <- "null"

# Table S1
tableS1 <- as.data.frame(re.sp.sum[,c("phylo","sporder","mean.r2","sig.r2","mean.abund","sig.abund",
                                      "mean.occupancy","sig.occup","TB.c")])
tableS1$mean.r2 <- format(round(tableS1$mean.r2,3),nsmall=3)
tableS1$mean.abund <- format(round(tableS1$mean.abund,3),nsmall=3)
tableS1$mean.occupancy <- format(round(tableS1$mean.occupancy,3),nsmall=3)
tableS1$mean.r2[!tableS1$sig.r2] <- paste(tableS1$mean.r2[!tableS1$sig.r2],"*",sep="")
tableS1$mean.abund[!tableS1$sig.abund] <- paste(tableS1$mean.abund[!tableS1$sig.abund],"*",sep="")
tableS1$mean.occupancy[!tableS1$sig.occup] <- paste(tableS1$mean.occupancy[!tableS1$sig.occup],"*",sep="")
tableS1 <- tableS1[,c("phylo","sporder","mean.r2","mean.abund","mean.occupancy")]
setwd("/Users/Viana/Documents/Papers/iDiv/Paper_BBS_R2-time")
write.table(tableS1,file="Table_S1.txt",row.names = FALSE, sep="\t")
save(re.sp.sum,file="re.sp.sum.Rdata")


# Figure S1: Phylogeny
setwd("/Users/Viana/Documents/Papers/iDiv/Paper_BBS_R2-time/Phylo_100_BirdTree")
sp.trees <- ape::read.nexus("output.nex")
library(phytools)
trait.y<-setNames(re.sp.sum$mean.r2,re.sp.sum$phylo)
# Plot
library(RColorBrewer)
pal <-  colorRampPalette(c("blue","yellow","red"))
cons.tree <- consensus.edges(sp.trees,method="mean.edge")
cons.tree <- multi2di(cons.tree)
obj<-contMap(cons.tree,trait.y,plot=FALSE)
obj<-setMap(obj,colors=pal(20)) # colors=rev(brewer.pal(11,"RdYlBu"))
quartz(height=5.5,width=5.5)
plot(obj,fsize=c(0.5,0.7),outline=FALSE,lwd=c(2,7),leg.txt="",legend=30,
     type="fan",xlim=c(-140,145),ylim=c(-120,120))


# Figure S2: trends in climate matching vs traits

gg.time.hwi <- re.sp.sum %>%
  ggplot(aes(x=HWI, y=mean.r2)) + 
  geom_point(col="black", size=1) +
  labs(x = "Hand-Wing Index", y = expression(paste(Delta,"climate matching"))) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12))

gg.time.mig <- re.sp.sum %>%
  ggplot(aes(x=Migratory_status, y=mean.r2)) + 
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x ="Migratory status", y = expression(paste(Delta,"climate matching"))) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust = 1))

gg.time.habitat <- re.sp.sum %>%
  ggplot(aes(x=main.habitat, y=mean.r2)) + 
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x ="Main habitat", y = expression(paste(Delta,"climate matching"))) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Compose figure 
quartz(height=4,width=9)
ggarrange(gg.time.hwi,gg.time.mig,gg.time.habitat, nrow=1)

# Correlations traits-trends climate matching
cor.test(re.sp.sum$HWI,re.sp.sum$mean.r2)
summary(aov(mean.r2~Migratory_status,data=re.sp.sum))
summary(aov(mean.r2~main.habitat,data=re.sp.sum))


# Figure S3: baseline values of climate matching vs traits
gg.time.bm <- re.sp.sum %>%
  ggplot(aes(x=Body.mass..log., y=baseline.r2)) + 
  geom_point(col="black", size=1) +
  labs(x = "Body mass (log)", y = "Baseline climate matching") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12))

gg.time.ssi <- re.sp.sum %>%
  ggplot(aes(x=SSI, y=baseline.r2)) + 
  geom_point(col="black", size=1) +
  labs(x = "Habitat specialisation", y = "") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12))

gg.time.hwi <- re.sp.sum %>%
  ggplot(aes(x=HWI, y=baseline.r2)) + 
  geom_point(col="black", size=1) +
  labs(x = "Hand-Wing Index", y = "") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12))

gg.time.mig <- re.sp.sum %>%
  ggplot(aes(x=Migratory_status, y=baseline.r2)) + 
  geom_boxplot() +
  labs(x ="Migratory status", y = "Baseline climate matching") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust = 1))

gg.time.habitat <- re.sp.sum %>%
  ggplot(aes(x=main.habitat, y=baseline.r2)) + 
  geom_boxplot() +
  labs(x ="Main habitat", y = "") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Compose figure 
quartz(height=6,width=9)
ggarrange(gg.time.bm,gg.time.ssi,gg.time.hwi,gg.time.mig,gg.time.habitat, nrow=2)

# Correlations traits-baseline climate matching
cor.test(re.sp.sum$SSI,re.sp.sum$baseline.r2)
cor.test(re.sp.sum$Body.mass..log.,re.sp.sum$baseline.r2)
cor.test(re.sp.sum$HWI,re.sp.sum$baseline.r2)
summary(aov(baseline.r2~Migratory_status,data=re.sp.sum))
summary(aov(baseline.r2~main.habitat,data=re.sp.sum))


# Extract population trends from Sauer et al. 2011
library(tabulizer)
tab <- extract_tables("https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/auk/128/1/10.1525_auk.2010.09220/3/sd1_87.pdf?Expires=1643015249&Signature=kp06HUCzQA0YY9dVysOHdIF2ZmowUgMbT5ahfKBza89n2IRC~S5MatM3V24wlRQ8nd5sK9905sXAzhTLEMSrxA3z~0nA7vBGwDHZFhvwI2lJU-JKRyyuyL6~awm05x7HB~eixwDiaMsIBHZaLnfVYPpwAjTPQDYaRRFKyHs99K~DN6uH1O6ho0kPabW~P3lV5KweJ9lVXUK3SzBL3KUBMMonGBsHCyF2GcUrGpXAXjOcdd7r4gS9NJRI5ziklVb-a~T1ZsDybXXVNHTrGlBhzI~s7OMtRzfFjD9fDczBSJl8U6BmrgTYmJQFFikf5fiTTOq3W-toM4F4ysQmRsZ0Mg__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA",
                      output = "data.frame",pages=1:9)

tab2 <- list()
for(i in 1:9){
  if(i==1){
    header <- 2
    col.trend <- 5
  } 
  if(i!=1){
    header <- 3
    col.trend <- 6
  } 
  dati <- tab[[i]][-c(1:header),]
  trend <- gsub( " .*$", "", dati[,col.trend] )
  blank <- which(trend=="")
  for(j in blank) dati[j-1,1] <- paste(dati[j-1,1],dati[j,1],sep=" ")
  tab2[[i]] <- data.frame(species.Sauer=dati[-blank,1],trend.Sauer=trend[-blank])
}
trends.Sauer <- do.call("rbind",tab2)
trends.Sauer$phylo <- sub("(.*\\()(.*)(\\))", "\\2", trends.Sauer$species.Sauer)
trends.Sauer$phylo <- gsub(' ', '_', trends.Sauer$phylo)
substr(trends.Sauer$trend.Sauer[substr(trends.Sauer$trend.Sauer,1,1)=="–"],1,1) <- "-"
trends.Sauer$trend.Sauer <- as.numeric(trends.Sauer$trend.Sauer)

re.sp.sum <- left_join(re.sp.sum,trends.Sauer,by="phylo")

# Figure S5
quartz(height=3,width=7)
par(mfrow=c(1,2),mar=c(4,4,1,1))
plot(re.sp.sum$trend.Sauer,re.sp.sum$mean.abund,xlim=c(-6,4),xlab="Sauer's et al. estimates",ylab="Our estimates")
plot(re.sp.sum$PT.c,re.sp.sum$mean.abund,xlab="ACAD's estimates (categorised)",ylab="Our estimates")


# Simulation rationale
# The random error (in proportion) decreases as the mean increases
res <- c()
for(i in 1:100) res[i] <- mean(abs(i-rpois(1000,i))/i)
# Figure S6
quartz(height=4,width=5)
par(mar=c(4,4,1,1))
plot(1:100, res, xlab="Mean count", ylab="Proportion of error")

