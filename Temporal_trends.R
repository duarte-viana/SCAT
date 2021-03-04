###########################################
# Supporting Information
# "Increasing climatic decoupling of bird abundances and distributions"
# Duarte S. Viana, Jonathan M. Chase
###########################################

# Estimate temporal trends with "brms"

# Passing arguments to R from command lines
args = commandArgs(trailingOnly=TRUE)
output.file <- args[1]

# try to get SGE_TASK_ID from submit script, otherwise fall back to 1
task.id = as.integer(Sys.getenv("SGE_TASK_ID", "1"))

library(raster)
library(dplyr)
library(tidyr)
library(rstan)
library(brms)
# try to get NSLOTS from submit script, otherwise fall back to 1
slots = as.integer(Sys.getenv("NSLOTS", "1"))
rstan_options(auto_write = TRUE)
options(mc.cores = slots)

# Load data
load("Responsivenes_data.Rdata")
load("Data_id.Rdata")

# Define the the different approaches
time.extent <- 1:3
data.type <- c("abund","pres")
method <- c("BRT_lr.01","BRT_lr.001","GAM_kest","GAM_k3")
r2.type <- 3:4
meth <- expand.grid(time.extent=time.extent,data.type=data.type,method=method,r2.type=r2.type)
# Species performance metrics
resp <- c("abund.mean","range.size")

# output list
mm <- list()

# Responsiveness data (species-climate r2)
r2s <- as.data.frame(res[which(cases$time.extent==meth$time.extent[task.id] & 
                                 cases$data.type==meth$data.type[task.id] & 
                                 cases$method==meth$method[task.id]),,meth$r2.type[task.id]])
# species performance
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
dat.in <- names(N[N>=15]) # at least 15 years
lresi <- lresi[lresi$aou %in% dat.in,]
lresi$species <- factor(lresi$aou)
lresi$r2[lresi$r2<=0] <- 0.0001
lresi$time <- scale(lresi$time)
for(i in unique(lresi$species)){
  lresi$abund.mean[lresi$species==i] <- lresi$abund.mean[lresi$species==i]/max(lresi$abund.mean[lresi$species==i])
  lresi$range.size[lresi$species==i] <- lresi$range.size[lresi$species==i]/max(lresi$range.size[lresi$species==i])
}

# Fit multilevel models
m.time <- brm(r2 ~ time + (time|species),
              data=lresi, family=Beta, chains=4, thin=1, iter=3000, warmup=1000, 
              control = list(adapt_delta = .99), cores=4)
mm[[1]] <- m.time
rn <- 1
for(r in resp){
  rn <- rn + 1
  formula.r <- as.formula(paste(r, " ~ time + (time|species)", sep=""))
  m.time <- brm(formula.r, data=lresi, family=gaussian, chains=4, thin=1, iter=3000, warmup=1000, cores=4)
  mm[[rn]] <- m.time
}

save(mm,file=output.file)

