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

# Load data
load("Data_id.Rdata")
setwd("~/Documents/Papers/iDiv/Paper_BBS_R2-time")
load("species_names.Rdata")
sps$aou <- as.character(sps$aou)
iucn <- read.csv("IUCN_status.csv")
iucn$sci <- paste(iucn$genusName,iucn$speciesName,sep="_")
traits <- read.table("traits.txt",header=T,sep="\t")
names(traits) <- c("species","lifespan","log.mass","body.length","wingspan",names(traits[,6:ncol(traits)]))

# Identify analytical variations
time.extent <- 1:3
data.type <- c("abund","pres")
method <- c("BRT_lr.01","BRT_lr.001","GAM_kest","GAM_k3")
r2.type <- 3:4
meth <- expand.grid(time.extent=time.extent,data.type=data.type,method=method,r2.type=r2.type)
# Define the different species performance metrics
resp <- c("abund.mean","range.size","range.x","range.y")

#-------------------------

# Pool posteriors of the temporal trends in responsiveness (r2) 
# estimated with the different species-climate models

wcc <- which(meth$time.extent==1 & meth$data.type=="abund" & meth$r2.type==4)
res.list <- list()
pop.eff <- list()
for(i in 1:length(wcc)){
  load(paste("bbs-7474408-", wcc[i],".Rdata", sep = "")) # load brms model results (temporal trends)
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

# Combine posteriors
re.sp <- do.call("rbind",res.list)

# Take means and 90% confidence intervals
re.sp.sum <- re.sp %>%
  group_by(species) %>%
  summarize(mean.r2 = mean(coef.r2),lower.r2=quantile(coef.r2,0.05),upper.r2=quantile(coef.r2,0.95),
            mean.abund = mean(coef.abund.mean),lower.abund=quantile(coef.abund.mean,0.05),upper.abund=quantile(coef.abund.mean,0.95),
            mean.range.size = mean(coef.range.size),lower.range.size=quantile(coef.range.size,0.05),upper.range.size=quantile(coef.range.size,0.95),
            mean.range.x = mean(coef.range.x),lower.range.x=quantile(coef.range.x,0.05),upper.range.x=quantile(coef.range.x,0.95),
            mean.range.y = mean(coef.range.y),lower.range.y=quantile(coef.range.y,0.05),upper.range.y=quantile(coef.range.y,0.95))

# Add species names
sp.all <- unique(re.sp.sum$species)
sp.names <- data.frame(aou=sp.all)
sp.names <- left_join(sp.names,sps[,c("aou","sci","sporder")],by="aou")
birdtree.sp <- read.csv("BLIOCPhyloMasterTax.csv")
all(sp.names$sci %in% birdtree.sp$TipLabel)

# Correct unmatched species names
mis.sp <- sp.names[!(sp.names$sci %in% birdtree.sp$TipLabel),]
mis.sp <- left_join(mis.sp, sps[,c("TipLabel","sci")], by="sci")
mis.sp$TipLabel <- as.character(mis.sp$TipLabel)
mis.sp$TipLabel[mis.sp$sci=="Circus_hudsonius"] <- "Circus_hudsonius"
#mis.sp$TipLabel[mis.sp$sci=="Ammospiza_leconteii"] <- "Ammodramus_leconteii"
mis.sp$TipLabel[mis.sp$sci=="Setophaga_coronata"] <- "Dendroica_coronata"
#mis.sp$TipLabel[mis.sp$sci=="Setophaga_nigrescens"] <- "Dendroica_nigrescens"
sp.names$phylo <- sp.names$sci
for(i in mis.sp$sci) sp.names$phylo[sp.names$sci==i] <- mis.sp$TipLabel[mis.sp$sci==i]
# Add species names
re.sp.sum <- left_join(re.sp.sum,sp.names,by=c("species"="aou"))

# Add threat status
iucn$sci <- paste(iucn$genusName,iucn$speciesName,sep="_")
re.sp.sum <- left_join(re.sp.sum,iucn[,c("sci","redlistCategory")],by="sci")
re.sp.sum$redlistCategory <- as.character(re.sp.sum$redlistCategory)
re.sp.sum$redlistCategory[re.sp.sum$phylo=="Picoides_villosus"] <- "Least Concern"
re.sp.sum$redlistCategory[re.sp.sum$phylo=="Dryocopus_pileatus"] <- "Least Concern"

# Add species traits
re.sp.sum1 <- left_join(re.sp.sum,traits,by=c("phylo"="species"))
re.sp.sum2 <- left_join(re.sp.sum,traits,by=c("sci"="species"))
re.sp.sum1[is.na(re.sp.sum1$log.mass),] <- re.sp.sum2[is.na(re.sp.sum1$log.mass),]
re.sp.sum <- as.data.frame(re.sp.sum1)

#-------------------------

# Figure 2: R2 temporal trend according to the different methods

# Population-level effect
wcc <- which(meth$time.extent==1 & meth$data.type=="abund" & meth$r2.type==4)
pop.est <- do.call("c", pop.eff)
mean(pop.est)
quantile(pop.est,c(0.05,0.95))
lapply(pop.eff,length)
pop.r2 <- data.frame(method=rep(meth$method[wcc],each=8000), b=pop.est)

quartz(height=4,width=4)
pop.r2 %>%
  mutate(grid_mean = b) %>%
  ggplot(aes(y = method, x = grid_mean)) + # , fill = stat(abs(x) < .9)
  #geom_hline(yintercept = levels(factor(re.sp$sporder)), col="grey") +
  stat_halfeyeh(.width = c(0.50, 0.90)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x =expression(paste(Delta,"responsiveness to climate")), y = "Species-climate model") +
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))

#-------------------------

# Figure 3: Group R2 temporal trends by threat status and taxonomic order

# Add taxonomic order and threat status
re.sp <- left_join(re.sp,re.sp.sum[,c("species","sporder","redlistCategory")],by="species")

# Figures posterior for order and threat status

gg.time.order <- re.sp %>%
  mutate(grid_mean = coef.r2) %>%
  ggplot(aes(y = sporder, x = grid_mean)) + # , fill = stat(abs(x) < .9)
  #geom_hline(yintercept = levels(factor(re.sp$sporder)), col="grey") +
  stat_halfeyeh(.width = c(0.50, 0.90)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = expression(paste(Delta,"responsiveness to climate")), y = "Order") 

gg.time.threat <- na.exclude(re.sp) %>%
  mutate(grid_mean = coef.r2) %>%
  ggplot(aes(y = redlistCategory, x = grid_mean)) + # , fill = stat(abs(x) < .9)
  #geom_hline(yintercept = levels(factor(re.sp$redlistCategory)), col="grey") +
  stat_halfeyeh(.width = c(0.50, 0.90)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = expression(paste(Delta,"responsiveness to climate")), y = "Threat status") 

gg.time.order.mean <- re.sp.sum %>%
  ggplot(aes(x=sporder, y=mean.r2)) + 
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x ="Order", y = expression(paste(Delta,"responsiveness to climate"))) +
  coord_flip()

gg.time.threat.mean <- re.sp.sum %>%
  ggplot(aes(x=redlistCategory, y=mean.r2)) + 
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x ="Threat status", y = expression(paste(Delta,"responsiveness to climate"))) +
  coord_flip()

# Compose figure 
quartz(height=7,width=7)
ggarrange(gg.time.order, gg.time.threat, gg.time.order.mean, gg.time.threat.mean, nrow=2)

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
  geom_point(data=re.sp.sum[!sig.r2,], aes(x=mean.r2,y=mean.abund), col="red", size=1) +
  geom_point(data=re.sp.sum[!sig.abund,], aes(x=mean.r2,y=mean.abund), col="blue", size=1) +
  geom_point(data=re.sp.sum[sig.r2.abund,], aes(x=mean.r2,y=mean.abund), col="purple", size=1) +
  geom_vline(xintercept = 0, linetype = "dashed", size=1, col="red", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", size=1, col="blue", alpha = 0.5) +
  geom_smooth(method="lm", col="black") +
  labs(x = expression(paste(Delta,"responsiveness to climate")), y = expression(paste(Delta,"abundance")))

sig.rgsi <- re.sp.sum$lower.range.size<=0 & re.sp.sum$upper.range.size>=0
sig.r2.rgsi <- sig.r2==FALSE & sig.rgsi==FALSE
gg.range.size <- re.sp.sum %>%
  ggplot(aes(x=mean.r2, y=mean.range.size)) + 
  geom_point(col="darkgrey", size=1) +
  geom_errorbar(aes(ymin=lower.range.size, ymax=upper.range.size), col="darkgrey",
                position=position_dodge(0.05)) +
  geom_errorbar(aes(xmin=lower.r2, xmax=upper.r2), col="darkgrey",
                position=position_dodge(0.05)) +
  geom_point(data=re.sp.sum[!sig.r2,], aes(x=mean.r2,y=mean.range.size), col="red", size=1) +
  geom_point(data=re.sp.sum[!sig.rgsi,], aes(x=mean.r2,y=mean.range.size), col="blue", size=1) +
  geom_point(data=re.sp.sum[sig.r2.rgsi,], aes(x=mean.r2,y=mean.range.size), col="purple", size=1) +
  geom_vline(xintercept = 0, linetype = "dashed", size=1, col="red", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", size=1, col="blue", alpha = 0.5) +
  geom_smooth(method="lm", col="black") +
  labs(x = expression(paste(Delta,"responsiveness to climate")), y = expression(paste(Delta,"range size")))

sig.rgy <- re.sp.sum$lower.range.y<=0 & re.sp.sum$upper.range.y>=0
sig.r2.rgy <- sig.r2==FALSE & sig.rgy==FALSE
gg.range.y <- re.sp.sum %>%
  ggplot(aes(x=mean.r2, y=mean.range.y)) + 
  geom_point(col="darkgrey", size=1) +
  geom_errorbar(aes(ymin=lower.range.y, ymax=upper.range.y), col="darkgrey",
                position=position_dodge(0.05)) +
  geom_errorbar(aes(xmin=lower.r2, xmax=upper.r2), col="darkgrey",
                position=position_dodge(0.05)) +
  geom_point(data=re.sp.sum[!sig.r2,], aes(x=mean.r2,y=mean.range.y), col="red", size=1) +
  geom_point(data=re.sp.sum[!sig.rgy,], aes(x=mean.r2,y=mean.range.y), col="blue", size=1) +
  geom_point(data=re.sp.sum[sig.r2.rgy,], aes(x=mean.r2,y=mean.range.y), col="purple", size=1) +
  geom_vline(xintercept = 0, linetype = "dashed", size=1, col="red", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", size=1, col="blue", alpha = 0.5) +
  geom_smooth(method="lm", col="black") +
  labs(x = expression(paste(Delta,"responsiveness to climate")), y = expression(paste(Delta,"range latitude")))


quartz(height=3,width=9)
ggarrange(gg.abund, gg.range.size, gg.range.y, nrow=1)


# Correlation between trends in performance and responsiveness to climate
cor.test(re.sp.sum$mean.r2,re.sp.sum$mean.abund)
cor.test(re.sp.sum$mean.r2,re.sp.sum$mean.range.size)
cor.test(re.sp.sum$mean.r2,re.sp.sum$mean.range.y)



#-------------------------

# Table S1
re.sp.sum$sig.r2 <- sig.r2
re.sp.sum$sig.abund <- sig.abund
re.sp.sum$sig.rgsi <- sig.rgsi
re.sp.sum$sig.rgy <- sig.rgy

tableS1 <- as.data.frame(re.sp.sum[,c("phylo","sporder","mean.r2","sig.r2","mean.abund","sig.abund",
                                          "mean.range.size","sig.rgsi","mean.range.y","sig.rgy",
                                          "redlistCategory","main.habitat")])
tableS1$mean.r2 <- format(round(tableS1$mean.r2,3),nsmall=3)
tableS1$mean.abund <- format(round(tableS1$mean.abund,3),nsmall=3)
tableS1$mean.range.size <- format(round(tableS1$mean.range.size,3),nsmall=3)
tableS1$mean.range.y <- format(round(tableS1$mean.range.y,3),nsmall=3)
tableS1$mean.r2[!tableS1$sig.r2] <- paste(tableS1$mean.r2[!tableS1$sig.r2],"*",sep="")
tableS1$mean.abund[!tableS1$sig.abund] <- paste(tableS1$mean.abund[!tableS1$sig.abund],"*",sep="")
tableS1$mean.range.size[!tableS1$sig.rgsi] <- paste(tableS1$mean.range.size[!tableS1$sig.rgsi],"*",sep="")
tableS1$mean.range.y[!tableS1$sig.rgy] <- paste(tableS1$mean.range.y[!tableS1$sig.rgy],"*",sep="")

tableS1 <- tableS1[,c("phylo","sporder","mean.r2","mean.abund","mean.range.size","mean.range.y",
                      "redlistCategory","main.habitat")]
setwd("/Users/Viana/Documents/Papers/iDiv/Paper_BBS_R2-time")
write.table(tableS1,file="Table_S1.txt",row.names = FALSE, sep="\t")


# Figure S1: responsiveness r2 as a function of habitat, migratory status, body mass, and life span
quartz(height=7,width=7)
par(mfrow=c(2,2))
par(mar=c(8,3,1,1))
boxplot(re.sp.sum$mean.r2~as.character(re.sp.sum$main.habitat),las=2,xlab="")
abline(h=0,lty=2)
boxplot(re.sp.sum$mean.r2~as.character(re.sp.sum$migration),las=2,xlab="")
abline(h=0,lty=2)
par(mar=c(5,3,0.1,1))
plot(re.sp.sum$log.mass,re.sp.sum$mean.r2,xlab="Body mass (log)",las=1,pch=16)
plot(re.sp.sum$lifespan,re.sp.sum$mean.r2,xlab="Life span",las=1,pch=16)









