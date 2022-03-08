###########################################
# Supporting Information
# "Increasing climatic decoupling of bird abundances and distributions"
# Duarte S. Viana, Jonathan M. Chase
###########################################

# Estimate temporal trends with "brms"

# Passing arguments to R from command lines
args = commandArgs(trailingOnly=TRUE)
output.file <- args[1]

# try to get SLURM_ARRAY_TASK_ID  from submit script, otherwise fall back to 1
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

library(dplyr)
library(rstan)
library(brms)
#options(buildtools.check = function(action) TRUE )
# try to get SLURM_CPUS_PER_TASK from submit script, otherwise fall back to 1
cpus_per_task = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
rstan_options(auto_write = TRUE)
options(mc.cores = cpus_per_task)

# Load data
load("results_v3.Rdata")
load("Data_id_v2.Rdata")
load("phylo_vcv.Rdata")

# Define the the different approaches
time.extent <- 1:3
data.type <- c("abund","pres")
method <- c("BRT_lr.01","BRT_lr.001","GAM_kest","GAM_k3")
r2.type <- 3:4
meth <- expand.grid(time.extent=time.extent,data.type=data.type,method=method,r2.type=r2.type)
# Species performance metrics
resp <- c("abund.mean","abund.cv", "occupancy")

# output list
mm <- list()

# Climate matching (species-climate r2)
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
lresi <- left_join(lresi,spdf[,c("aou","phylo")],by="aou")
lresi$species <- factor(lresi$aou)
lresi$phylo <- factor(lresi$phylo)
lresi$r2[lresi$r2<=0] <- 0.0001
lresi$time <- scale(lresi$time)
lresi$occupancy <- lresi$Npres/unlist(lapply(idata,length)[meth$time.extent[task.id]])
for(i in unique(lresi$species)){
  lresi$abund.mean[lresi$species==i] <- lresi$abund.mean[lresi$species==i]/max(lresi$abund.mean[lresi$species==i])
}


# Fit multilevel models
m.time <- brm(r2 ~ time + (1|gr(phylo, cov = A)) + (time|species),
              data=lresi, data2 = list(A = A), family=Beta, chains=4, thin=1, iter=3000, warmup=1000, 
              control = list(adapt_delta = .99), cores=4)
mm[[1]] <- m.time

m.abund.mean <- brm(abund.mean ~ time + (time|species),
                    data=lresi, family=gaussian, chains=4, thin=1, iter=3000, warmup=1000, 
                    control = list(adapt_delta = .99), cores=4)
mm[[2]] <- m.abund.mean

m.abund.cv <- brm(abund.cv ~ time + (time|species),
                  data=lresi, family=gaussian, chains=4, thin=1, iter=3000, warmup=1000, 
                  control = list(adapt_delta = .99), cores=4)
mm[[3]] <- m.abund.cv

m.occup <- brm(occupancy ~ time + (time|species),
               data=lresi, family=Beta, chains=4, thin=1, iter=3000, warmup=1000, 
               control = list(adapt_delta = .99), cores=4)
mm[[4]] <- m.occup

save(mm,file=output.file)

