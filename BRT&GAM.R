###########################################
# Supporting Information
# "Increasing climatic decoupling of bird abundances and distributions"
# Duarte S. Viana, Jonathan M. Chase
###########################################

# R functions

# ------------------------------------------------------------------------------
# BOOSTED REGRESSION TREES from package 'gbm'
# ------------------------------------------------------------------------------

R2.BRT <- function(Y, X, train.data=NULL, test.data=NULL, row.id=NULL, distr="poisson", CV = FALSE, inter.depth = 2, lr=0.01)
{
  require(gbm)
  
  # CV folds
  if(CV){
    cv.folds = 5  # 5-fold cross-validation on the training set
  } else cv.folds = 0
  
  # fit the boosted trees
  if(is.null(train.data)){
    d2s <- c()
    r2s <- c()
    for(i in 1:ncol(Y)){
      Y.i <- Y[,i]
      brt <- gbm(Y.i ~ ., data=X,
                 distribution = distr, 
                 shrinkage = lr, # a.k.a. learning rate, or weight of each tree
                 interaction.depth = inter.depth,
                 n.trees = 1000,
                 bag.fraction = 0.75,
                 cv.folds = cv.folds,
                 n.cores = 1) 
      if(CV){
        best.iter <- gbm.perf(brt, method = "cv", plot.it = FALSE)  
        # predicted values, using the best number of trees
      } else{
        best.iter <- 1000
      } 
      # calculate R2s
      # R2 with training data and test (CV) data
      preds <- predict(brt, newdata = X, n.trees = best.iter, type="response")
      if(distr=="poisson"){
        m.null <- mean(Y.i)
        d2s[i] <- D2.Poisson(Y.i, preds, m.null)
        r2s[i] <- cor(Y.i, preds, method="spearman", use='pairwise')^2
      }
      if(distr=="bernoulli"){
        d2s[i] <- D2.binom(Y.i, preds)
        r2s[i] <- cor(Y.i, preds, method="spearman", use='pairwise')^2
      }
    }
    res <- list(d2s,r2s)
  }
  
  if(!is.null(train.data)){
    d2s <- matrix(ncol=ncol(Y), nrow=length(train.data))
    r2s <- matrix(ncol=ncol(Y), nrow=length(train.data))
    d2s.cv <- matrix(ncol=ncol(Y), nrow=length(train.data))
    r2s.cv <- matrix(ncol=ncol(Y), nrow=length(train.data))
    for(j in 1:length(train.data)){
      #print(paste("j=",j))
      Y.train <- Y[site.id %in% train.data[[j]],]
      Y.test <- Y[site.id %in% test.data[[j]],]
      X.train <- X[site.id %in% train.data[[j]],]
      X.test <- X[site.id %in% test.data[[j]],]
      for(i in 1:ncol(Y)){
        #print(paste("i=",i))
        Y.i <- Y.train[,i]
        brt <- gbm(Y.i ~ ., data=X.train,
                   distribution = distr, 
                   shrinkage = lr, # a.k.a. learning rate, or weight of each tree
                   interaction.depth = inter.depth,
                   n.trees = 1000,
                   bag.fraction = 0.75,
                   cv.folds = cv.folds,
                   n.cores = 1) 
        if(CV){
          best.iter <- gbm.perf(brt, method = "cv", plot.it = FALSE)  
          # predicted values, using the best number of trees
        } else{
          best.iter <- 1000
        } 
        # calculate R2s
        # R2 with training data and test (CV) data
        preds.test <- predict(brt, newdata = X.test, n.trees = best.iter, type="response")
        preds <- predict(brt, newdata = X.train, n.trees = best.iter, type="response")
        if(distr=="poisson"){
          m.null <- mean(Y.i)
          d2s.ji <- D2.Poisson(Y.i, preds, m.null)
          r2s.ji <- cor(Y.i, preds, method="spearman", use='pairwise')^2
          d2s.cv.ji <- D2.Poisson(Y.test[,i], preds.test, m.null)
          r2s.cv.ji <- cor(Y.test[,i], preds.test, method="spearman", use='pairwise')^2
          if(is.infinite(d2s.ji) || length(d2s.ji)==0) d2s.ji <- NA
          if(is.infinite(r2s.ji) || length(r2s.ji)==0) r2s.ji <- NA
          if(is.infinite(d2s.cv.ji) || length(d2s.cv.ji)==0) d2s.cv.ji <- NA
          if(is.infinite(r2s.cv.ji) || length(r2s.cv.ji)==0) r2s.cv.ji <- NA
          d2s[j,i] <- d2s.ji
          r2s[j,i] <- r2s.ji
          d2s.cv[j,i] <- d2s.cv.ji
          r2s.cv[j,i] <- r2s.cv.ji
        }
        if(distr=="bernoulli"){
          d2s.ji <- D2.binom(Y.i, preds)
          r2s.ji <- cor(Y.i, preds, method="spearman", use='pairwise')^2
          d2s.cv.ji <- D2.binom(Y.test[,i], preds.test)
          r2s.cv.ji <- R2.Tjur(Y.test[,i], preds.test)
          if(is.infinite(d2s.ji) || length(d2s.ji)==0) d2s.ji <- NA
          if(is.infinite(r2s.ji) || length(r2s.ji)==0) r2s.ji <- NA
          if(is.infinite(d2s.cv.ji) || length(d2s.cv.ji)==0) d2s.cv.ji <- NA
          if(is.infinite(r2s.cv.ji) || length(r2s.cv.ji)==0) r2s.cv.ji <- NA
          d2s[j,i] <- d2s.ji
          r2s[j,i] <- r2s.ji
          d2s.cv[j,i] <- d2s.cv.ji
          r2s.cv[j,i] <- r2s.cv.ji
        }
      }
    }
    # Set negative R2 values to 0
    d2s[d2s<0] <- 0
    r2s[r2s<0] <- 0
    d2s.cv[d2s.cv<0] <- 0
    r2s.cv[r2s.cv<0] <- 0
    # Mean across CV folds
    d2s <- apply(d2s,2,mean,na.rm=TRUE)
    r2s <- apply(r2s,2,mean,na.rm=TRUE)
    d2s.cv <- apply(d2s.cv,2,mean,na.rm=TRUE)
    r2s.cv <- apply(r2s.cv,2,mean,na.rm=TRUE)
    res <- list(d2s,r2s,d2s.cv,r2s.cv)
  }
  return(res)
}


# ------------------------------------------------------------------------------
# GENERALIZED ADDITIVE MODELS (GAM) FROM PACKAGE 'mgcv'
# ------------------------------------------------------------------------------

R2.GAM <- function(Y, X, train.data=NULL, test.data=NULL, row.id=NULL, family = "poisson", env.eff = "splines", k.env=-1)
{
  require(mgcv) 
  
  Y <- data.frame(Y)
  X0 <- data.frame(X)
  xy.indcs <- names(X0) %in% c("long","lat")
  
  # if only x or only y coordinate is provided
  if(sum(xy.indcs) == 1) return("Error, two coordinates must be provided.")
  
  k.spa <- -1
  operator <- " + "
  
  # COMPOSING THE GAM MODEL FORMULA ------------
  
  # if only environmental data are provided
  if(sum(xy.indcs) == 0)
  {
    if(env.eff=="linear"){ 
      env.part <- paste(names(X0), collapse = operator)
    } 
    if(env.eff=="quadratic"){ 
      env.part <- paste("poly(", names(X0), ", 2)", collapse = operator)
    } 
    if(env.eff=="splines"){ 
      env.part <- paste("s(", names(X0), ", k = ", k.env, ")", collapse = operator)
    } 
    
    Formula  <- as.formula(paste("Y.i ~", env.part)) 
  }
  
  # if only xy data are provided
  if(sum(xy.indcs) == ncol(X0))
  {
    Formula <- as.formula(paste("Y.i ~ s(long, lat, k = ", k.spa, ")", sep=""))  
  }
  
  # if xy and environmental data are provided
  if(sum(xy.indcs) == 2 && ncol(X0) > 2)
  {
    env.indcs <- xy.indcs == FALSE
    env.names <-  names(X0)[env.indcs]
    
    if(env.eff=="linear"){ 
      env.part <- paste(env.names, collapse = operator)
    } 
    if(env.eff=="quadratic"){
      env.part <- paste("poly(", env.names, ", 2)", collapse = operator)
    }
    if(env.eff=="splines"){
      env.part <- paste("s(", env.names, ", k = ", k.env, ")", collapse = operator)
    }
    
    Formula <- as.formula(paste(paste("Y.i ~ s(long, lat, k=", k.spa, ")", sep=""), env.part, sep= " + ")) 
  }
  
  # FITTING THE GAM MODELS -----------------
  
  if(is.null(train.data)){
    # the model fitting loop -- doing GAM for each species
    d2s <- c()
    r2s <- c()
    #pvals <- c()
    for(i in 1:ncol(Y)){
      #print(i)
      Y.i <- Y[,i]
      XY <- data.frame(Y.i = Y.i, X0)
      # Fit the GAM model
      m.i <- gam(Formula, data = XY, family=family, method="REML", select=TRUE)
      
      # Return p-value?
      # if(family=="poisson") m0.i <- gam(Y.i ~ 1, family=family, method="REML")
      # if(family=="nb"){
      #   size.m.i <- m.i$family$getTheta(TRUE)
      #   m0.i <- gam(Y.i ~ 1, family=negbin(size.m.i), method="REML")
      # }
      # if(family=="binomial") m0.i <- gam(Y.i ~ 1, family=family, method="REML")
      # pvals[i] <- anova(m0.i, m.i, test="Chisq")$`Pr(>Chi)`[2]
      
      # calculate R2s
      preds <- predict(m.i, type="response")
      if(family=="poisson"){
        m.null <- mean(Y.i)
        d2s[i] <- D2.Poisson(Y.i, preds, m.null)
        r2s[i] <- cor(Y.i, preds, method="spearman", use='pairwise')^2
      }
      if(family=="binomial"){
        d2s[i] <- D2.binom(Y.i, preds)
        r2s[i] <- R2.Tjur(Y.i, preds)
      }
    }
    res <- list(d2s,r2s)
    #res <- list(d2s=d2s, r2s=r2s, pvals=pvals)
  }
  
  if(!is.null(train.data)){
    # the model fitting loop -- doing GAM for each species
    d2s <- matrix(ncol=ncol(Y), nrow=length(train.data))
    r2s <- matrix(ncol=ncol(Y), nrow=length(train.data))
    d2s.cv <- matrix(ncol=ncol(Y), nrow=length(train.data))
    r2s.cv <- matrix(ncol=ncol(Y), nrow=length(train.data))
    for(j in 1:length(train.data)){
      Y.train <- Y[site.id %in% train.data[[j]],]
      Y.test <- Y[site.id %in% test.data[[j]],]
      X.train <- X0[site.id %in% train.data[[j]],]
      X.test <- X0[site.id %in% test.data[[j]],]
      for(i in 1:ncol(Y)){
        #print(i)
        Y.i <- Y.train[,i]
        XY <- data.frame(Y.i = Y.i, X.train)
        # Fit the GAM model
        m.i <- gam(Formula, data = XY, family=family, method="REML", select=TRUE)
        
        # calculate R2s
        # R2 with training data and test (CV) data
        preds.test <- predict(m.i, newdata = X.test, type="response")
        preds <- predict(m.i, type="response")
        if(family=="poisson"){
          m.null <- mean(Y.i)
          d2s.ji <- D2.Poisson(Y.i, preds, m.null)
          r2s.ji <- cor(Y.i, preds, method="spearman", use='pairwise')^2
          d2s.cv.ji <- D2.Poisson(Y.test[,i], preds.test, m.null)
          r2s.cv.ji <- cor(Y.test[,i], preds.test, method="spearman", use='pairwise')^2
          if(is.infinite(d2s.ji) || length(d2s.ji)==0) d2s.ji <- NA
          if(is.infinite(r2s.ji) || length(r2s.ji)==0) r2s.ji <- NA
          if(is.infinite(d2s.cv.ji) || length(d2s.cv.ji)==0) d2s.cv.ji <- NA
          if(is.infinite(r2s.cv.ji) || length(r2s.cv.ji)==0) r2s.cv.ji <- NA
          d2s[j,i] <- d2s.ji
          r2s[j,i] <- r2s.ji
          d2s.cv[j,i] <- d2s.cv.ji
          r2s.cv[j,i] <- r2s.cv.ji
        }
        if(family=="binomial"){
          d2s.ji <- D2.binom(Y.i, preds)
          r2s.ji <- cor(Y.i, preds, method="spearman", use='pairwise')^2
          d2s.cv.ji <- D2.binom(Y.test[,i], preds.test)
          r2s.cv.ji <- R2.Tjur(Y.test[,i], preds.test)
          if(is.infinite(d2s.ji) || length(d2s.ji)==0) d2s.ji <- NA
          if(is.infinite(r2s.ji) || length(r2s.ji)==0) r2s.ji <- NA
          if(is.infinite(d2s.cv.ji) || length(d2s.cv.ji)==0) d2s.cv.ji <- NA
          if(is.infinite(r2s.cv.ji) || length(r2s.cv.ji)==0) r2s.cv.ji <- NA
          d2s[j,i] <- d2s.ji
          r2s[j,i] <- r2s.ji
          d2s.cv[j,i] <- d2s.cv.ji
          r2s.cv[j,i] <- r2s.cv.ji
        }
      }
    }
    # Set negative R2 values to 0
    d2s[d2s<0] <- 0
    r2s[r2s<0] <- 0
    d2s.cv[d2s.cv<0] <- 0
    r2s.cv[r2s.cv<0] <- 0
    # Mean across CV folds
    d2s <- apply(d2s,2,mean,na.rm=TRUE)
    r2s <- apply(r2s,2,mean,na.rm=TRUE)
    d2s.cv <- apply(d2s.cv,2,mean,na.rm=TRUE)
    r2s.cv <- apply(r2s.cv,2,mean,na.rm=TRUE)
    res <- list(d2s,r2s,d2s.cv,r2s.cv)
  }
  return(res)
}



# ------------------------------------------------------------------------------
# Supporting function: Poisson deviance
# Arguments:
# Y - the observations
# mu - the predicted means (of Poisson distribution)
Poisson.deviance <- function(Y, mu){
  2*sum(ifelse(Y == 0, 
               -(Y-mu), # when the observation is 0
               Y*log(Y/mu) - (Y-mu))) # else
}


# ------------------------------------------------------------------------------
# Explained Poisson deviance calculated for each species, then averaged
D2.Poisson <- function(Y.obs, Y.pred, mean.obs){
  pred.D <- Poisson.deviance(Y.obs, Y.pred)
  null.D  <- Poisson.deviance(Y.obs, mean.obs)
  D2 <- 1 - pred.D/null.D
  return(D2)
}

# ------------------------------------------------------------------------------
# Explained Bernoulli deviance calculated for each species, then averaged

D2.binom <- function(Y.obs, Y.pred) {
  null.D = -2*sum(dbinom(x=Y.obs, size = 1, prob = sum(Y.obs/length(Y.obs)), log=TRUE))
  pred.D = -2*sum(dbinom(x=Y.obs, size = 1, prob = Y.pred, log=TRUE))
  D2 <- 1 - pred.D/null.D
  return(D2)
}

# ------------------------------------------------------------------------------
# TJUR PSEUDO R2
# This implementation is taken from  the PseudoR2 function from package DescTools

R2.Tjur <- function(Y.obs, Y.pred){
  unname(diff(tapply(Y.pred, Y.obs, mean, na.rm = TRUE)))
}





