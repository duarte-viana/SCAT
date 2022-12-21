###########################################
# Supporting Information
# "Increasing climatic decoupling of bird abundances and distributions"
# Duarte S. Viana, Jonathan M. Chase
###########################################

# Fit BRT and GAM models

# Passing arguments to R from command lines
args = commandArgs(trailingOnly=TRUE)
output.file <- args[1]

# try to get SGE_TASK_ID from submit script, otherwise fall back to 1
task.id = as.integer(Sys.getenv("SGE_TASK_ID", "1"))

library(tidyr)
library(gbm)
library(mgcv)
source("R_Functions.R")
# load data info
load("Data_id_v2.Rdata")
# load BBS data
load("BBS_BIOCLIM.Rdata")

# Prepare data

# subset data
wsites <- idata[[cases$time.extent[task.id]]]
dat <- bbs.clim[bbs.clim$year==cases$year[task.id] & bbs.clim$site_id %in% wsites,]
dat.wide <- spread(dat, key=species_id, value=abundance)
in.train <- itrain[[cases$time.extent[task.id]]]
in.test <- itest[[cases$time.extent[task.id]]]
sp.data <- data.id[data.id$time.extent==cases$time.extent[task.id] & data.id$year==cases$year[task.id],]
in.sp <- sp.data$aou[sp.data$Npres.train>=5 & sp.data$Npres.test>=5]
site.id <- dat.wide[dat.wide$site_id %in% wsites,"site_id"]

# Set data matrices
mat.sp0 <- dat.wide[dat.wide$site_id %in% wsites,24:ncol(dat.wide)]
mat.sp0[is.na(mat.sp0)] <- 0
mat.sp <- mat.sp0[,names(mat.sp0) %in% in.sp] # only include species with at least 5 presences in train data and 5 in test data
if(cases$data.type[task.id]=="pres") mat.sp[mat.sp>0] <- 1
# Climatic variables
mat.env <- dat.wide[dat.wide$site_id %in% wsites,c("bio2","bio3","bio5","bio8","bio9","bio15","bio16","bio18")]


# Prepare output object
sp.ids <- unique(bbs.clim$species_id)
res.vp <- array(NA, c(1, length(sp.ids), 4),
                dimnames=list(case=task.id,sp=paste(sp.ids),
                              stat=c("D2","R2","D2.CV","R2.CV")))


# Modelling
try({
  if(cases$method[task.id]=="BRT_lr.01"){
    # BRT
    if(cases$data.type[task.id]=="abund"){
      r2.env <- R2.BRT(Y=mat.sp, X=mat.env, train.data=in.train, test.data=in.test, row.id=site.id,
                                      distr="poisson", inter.depth=2, lr=0.01)
    } 
    if(cases$data.type[task.id]=="pres"){
      r2.env <- R2.BRT(Y=mat.sp, X=mat.env, train.data=in.train, test.data=in.test, row.id=site.id,
                                      distr="bernoulli", inter.depth=2, lr=0.01)
    } 
  }
  
  if(cases$method[task.id]=="BRT_lr.001"){
    # BRT
    if(cases$data.type[task.id]=="abund"){
      r2.env <- R2.BRT(Y=mat.sp, X=mat.env, train.data=in.train, test.data=in.test, row.id=site.id,
                                      distr="poisson", inter.depth=2, lr=0.001)
    } 
    if(cases$data.type[task.id]=="pres"){
      r2.env <- R2.BRT(Y=mat.sp, X=mat.env, train.data=in.train, test.data=in.test, row.id=site.id,
                                      distr="bernoulli", inter.depth=2, lr=0.001)
    } 
  }
  
  if(cases$method[task.id]=="GAM_kest"){
    # GAM with optimised splines 
    if(cases$data.type[task.id]=="abund"){
      r2.env <- R2.GAM(Y=mat.sp, X=mat.env, train.data=in.train, test.data=in.test, row.id=site.id,
                                      family="poisson", env.eff = "splines")
    }
    if(cases$data.type[task.id]=="pres"){
      r2.env <- R2.GAM(Y=mat.sp, X=mat.env, train.data=in.train, test.data=in.test, row.id=site.id,
                                      family="binomial", env.eff = "splines")
    }
  }
  
  if(cases$method[task.id]=="GAM_k3"){
    # GAM with fixed number of knots
    if(cases$data.type[task.id]=="abund"){
      r2.env <- R2.GAM(Y=mat.sp, X=mat.env, train.data=in.train, test.data=in.test, row.id=site.id,
                                      family="poisson", env.eff = "splines", k=3)
    }
    if(cases$data.type[task.id]=="pres"){
      r2.env <- R2.GAM(Y=mat.sp, X=mat.env, train.data=in.train, test.data=in.test, row.id=site.id,
                                      family="binomial", env.eff = "splines", k=3)
    }
  }
  
  # Store results
  res.vp[1,names(mat.sp),1] <- r2.env[[1]]
  res.vp[1,names(mat.sp),2] <- r2.env[[2]]
  res.vp[1,names(mat.sp),3] <- r2.env[[3]]
  res.vp[1,names(mat.sp),4] <- r2.env[[4]]
})

# Save output to a .Rdata file
save(res.vp,file=output.file)

