#-------------------------------------------------------------------------------------------------
# Title: Functions to start a simulation and analyses
# Project: Generation time in metacommunities
# Date: 21-12-2021
#-------------------------------------------------------------------------------------------------

# coordinates and environmental variables in space and time
make.land <- function(M = 50, steps = 2200, initial = 200, scale = 500, spa = 10, plot = TRUE){
  require(som.nn)
  require(vegan)
  require(RandomFields)
  # M = number of sites
  # steps = number of steps for simulation
  
  # space (xy)
  coord <- replicate(2, runif(M, 0, 100))
  
  # distances
  distance <- as.matrix(dist.torus(coord))
  
  # environmental variable
  model <- RMexp(var = 0.5, scale = scale)
  sim <- RFsimulate(model = model, x = coord[,1]*spa, y = coord[,2]*spa, T = 0:steps)
  env <- decostand(sim$variable1, "range")
  ecum <- ecdf(env)
  env <- ecum(env)
  
  # landscape
  land <- data.frame("x" = rep(coord[,1], steps+1), "y" = rep(coord[,2], steps+1))
  land$time <- rep(0:steps, each = M)
  land$site <- rep(1:M, steps+1)
  land$env <- env
  land$env[land$time < initial] <- rep(land$env[land$time == initial], initial)
  
  # plot
  if(plot == TRUE){
    par(mar = c(4,4,2,2))
    plot(land$env[land$site == unique(land$site)[1]], type = "l", ylim = c(0,1), xlab = "Time",
       ylab = "Environment")
    for(i in 2:M){
      lines(land$env[land$site == unique(land$site)[i]], col = i)
    }
  }
  
  return(list("distance" = distance, "land" = land))
}

#------------------------------------------------------------------------------------------------
# make traits
make.trait <- function(S = 2, stage = 2, opt = 0.5, breadth = 10, a = 0.4, w = 0,
                       pj = 0.9, m = 0.1, f = 10, pa = 0.1, intra = 1, inter = 0.5) {
  
  # colnames for traits table
  names <- c("species", "stage", "opt", "breadth", "a", "w", "pj", "m", "f", "pa") 
  
  if(S == 1 & stage ==1){
    trait <- as.data.frame(matrix(0, 1, 10))
    colnames(trait) <- names
    trait$species <- 1
    trait$stage <- 1
    trait$opt <- opt
    trait$breadth <- breadth
    trait$a <- a
    trait$w <- w
    trait$f <- f
    
    compet <- matrix(intra)
  }
  
  if(S >= 2 & stage == 1){
    trait <- as.data.frame(matrix(0, S, 10))
    colnames(trait) <- names
    trait$species <- 1:S
    trait$stage <- 1
    trait$opt <- opt
    trait$breadth <- breadth
    trait$a <- a
    trait$w <- w
    trait$f <- f
    
    compet <- matrix(inter, S, S)
    diag(compet) <- intra
  }
  
  if(S == 1 & stage >= 2){
    trait <- as.data.frame(matrix(0, stage, 10))
    colnames(trait) <- names
    trait$species <- 1
    trait$stage <- 1:stage
    trait$opt <- opt
    trait$breadth <- breadth
    trait$a <- a
    trait$w <- w
    trait$pj[trait$stage == 1] <- pj
    trait$m[trait$stage == 1] <- m
    trait$f[trait$stage == 2] <- f
    trait$pa[trait$stage == 2] <- pa
    
    compet <- matrix(intra)
  }
  
  if(S >= 2 & stage >= 2){
    trait <- as.data.frame(matrix(0, S*stage, 10))
    colnames(trait) <- names
    trait$species <- rep(1:S, each = stage)
    trait$stage <- rep(1:stage, S)
    trait$opt <- rep(opt, each = stage)
    trait$breadth <- breadth
    trait$a <- a
    trait$w <- w
    trait$pj[trait$stage == 1] <- pj
    trait$m[trait$stage == 1] <- m
    trait$f[trait$stage == 2] <- f
    trait$pa[trait$stage == 2] <- pa
    
    compet <- matrix(inter, S, S)
    diag(compet) <- intra
  }
  
  ## list of objects
  return(list("trait" = trait, "compet" = compet))
}

#----------------------------------------------------------------------------------------------------
# dispersal kernel

make.w <- function(dist, w){
  listw <- as.list(1:length(w))
  for(i in 1:length(listw)){
    listw[[i]] <- exp(-w[i] * dist^2)
  }
  listw <- do.call(rbind, listw)
  w <- cbind(rep(1:length(w), each = nrow(dist)), listw)
  return(w)
}

#----------------------------------------------------------------------------------------------------
# join species records

join.sp <- function(sp, names, stage){
  join <- as.data.frame(matrix(0, nrow(sp), ncol(sp)/length(unique(stage))))
  for(j in 1:ncol(join)){
    join[,j] <- rowSums(sp[,names == j])
  }
  return(join)
}

#----------------------------------------------------------------------------------------------------
# GAM multi
gam.multi <- function(sp, pred, family = "poisson", select = FALSE, method  = "REML",...){
  sp <- sp[,colSums(sp) !=0]
  fm.pred <- paste("s(", colnames(pred), ")", sep = "", collapse = " + ") 
  r.sq.vec <- rep(0, ncol(sp))
  for(j in 1:ncol(sp)){
    fm <- as.formula(paste("sp[,",j,"] ~ ", fm.pred, sep = "")) #Formula for each species model
    tryCatch({
      r.sq.vec[j] <- summary(mgcv::gam(fm, data = cbind(sp[,j], pred),
                                       family = family, select = select, method = method))$r.sq  
    }, error=function(e){})
    }
  r.sq <- sum(r.sq.vec * apply(sp, 2, var)) / sum(apply(sp, 2, var)) #multivariate r2
  return(r.sq)
}

#------------------------------------------------------------------------------------------------
# plots

plot.lines <- function(x, y, factor, color,...){
  # median, q1 and q3
  tab.median <- as.data.frame(matrix(0, nrow = length(unique(x)) * length(unique(factor)), ncol = 5))
  colnames(tab.median) <- c("factor", "x", "q1", "q2", "q3")
  tab.median$factor <- rep(unique(factor), each = length(unique(x)))
  tab.median$x<- log10(rep(unique(x), length(unique(factor))))
  tab.median$q1 <- tapply(log10(y), interaction(x, factor), quantile, na.rm = TRUE, probs = 0.25)
  tab.median$q2 <- tapply(log10(y), interaction(x, factor), quantile, na.rm = TRUE, probs = 0.50)
  tab.median$q3 <- tapply(log10(y), interaction(x, factor), quantile, na.rm = TRUE, probs = 0.75)
  # plot
  plot(tab.median$x, tab.median$q2, type = "n", axes = FALSE, ...)
  for(i in 1:length(unique(factor))){
    polygon(c(unique(tab.median$x), rev(unique(tab.median$x))),
            c(tab.median$q1[tab.median$factor == unique(tab.median$factor)[i]],
              rev(tab.median$q3[tab.median$factor == unique(tab.median$factor)[i]])),
            border = NA, col = "gray90")
  }
  for(i in 1:length(unique(factor))){
    lines(unique(tab.median$x), c(tab.median$q2[tab.median$factor == unique(tab.median$factor)[i]]),
          col = color[i], ...)
  }
}

#----------------------------------------------------------------------------------------------------------
#plot smooth lines
plot.slines <- function(x, y, factor, color, lg = TRUE, ...){
  if(lg == TRUE){
    x <- log10(x)
  }
  # loess
  lo <- as.list(1:length(unique(factor)))
  for(i in 1:length(lo)){
    lo[[i]] <- loess(y[factor == unique(factor)[i]] ~ x[factor == unique(factor)[i]])
  }
  
  # plot
  plot(x, y, type = "n", axes = FALSE, ...)
  for(i in 1:length(unique(factor))){
    new.x <- seq(min(x[factor == unique(factor)[i]]), max(x[factor == unique(factor)[i]]), 0.1)
    pr <- predict(lo[[i]], newdata = new.x, se = TRUE)
    polygon(c(new.x, rev(new.x)), c(pr$fit + pr$se.fit, rev(pr$fit - pr$se.fit)),
            border = NA, col = "gray90")
  }
  for(i in 1:length(lo)){
    lines(new.x, predict(lo[[i]], newdata = new.x), col = color[i], ...)
  }
}
