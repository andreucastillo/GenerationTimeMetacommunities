#-------------------------------------------------------------------------------------------------
# Title: Results simulations
# Project: Generation times in metacommunities
# Date: 17-12-2021
#-------------------------------------------------------------------------------------------------
rm(list = ls())
source("functionSim.R")

# load results
path <- "results"
set <- list.files(path, pattern = "*ni05", full.names = TRUE)

res <- list(1:length(set))
for(i in 1:length(set)){
  res[[i]] <- read.table(set[i], header = TRUE, sep = ",")
}

# plot settings
library(viridis)
color <- c(viridis(10, direction = -1), "black") 
palette(color)

for(i in 1:2){
  res[[i]] <- res[[i]][res[[i]]$m %in% c(0.1, 0.3, 0.5, 0.7, 0.9, 1),]
}

#---------------------------------------------------------------------------------------------------
# stability (480,460)
par(mfrow = c(3,2), oma = c(2,2,0.5,0.5), mar = c(2,2,0.5,0.5))

for(i in 1:2){
  plot.slines(x = res[[i]]$a, y = res[[i]]$a.stab, factor = paste(res[[i]]$pa, res[[i]]$m),
              color = color, ylim = c(0, 2), xlab = "dispersal", ylab = "a.cv", lwd = 1.5)
  axis(1, at = seq(-5, 0), labels = 10^seq(-5, 0))
  axis(2)
  box()
}

for(i in 1:2){
  plot.slines(x = res[[i]]$a, y = res[[i]]$g.stab, factor = paste(res[[i]]$pa, res[[i]]$m),
              color = color, ylim = c(0, 0.2), xlab = "dispersal", ylab = "g.cv", lwd = 1.5)
  axis(1, at = seq(-5, 0), labels = 10^seq(-5, 0))
  axis(2)
  box()
}

for(i in 1:2){
  plot.slines(x = res[[i]]$a, y = res[[i]]$N, factor = paste(res[[i]]$pa, res[[i]]$m),
              color = color, ylim = c(0,320), xlab = "dispersal", ylab = "N", lwd = 1.5)
  axis(1, at = seq(-5, 0), labels = 10^seq(-5, 0))
  axis(2)
  box()
}

#---------------------------------------------------------------------------------------------------
# diversity (480, 770)
par(mfrow = c(5,2), oma = c(2,2,0.5,0.5), mar = c(2,2,0.5,0.5))
for(i in 1:2){
  plot.slines(x = res[[i]]$a, y = res[[i]]$alpha, factor = paste(res[[i]]$pa, res[[i]]$m),
              color = color, ylim = c(0, 15), xlab = "dispersal", ylab = "alpha", lwd = 1.5)
  axis(1, at = seq(-5, 0), labels = 10^seq(-5, 0))
  axis(2)
  box()
}

for(i in 1:2){
  plot.slines(x = res[[i]]$a, y = res[[i]]$g.space, factor = paste(res[[i]]$pa, res[[i]]$m),
              color = color, ylim = c(0, 25), xlab = "dispersal", ylab = "g.space", lwd = 1.5)
  axis(1, at = seq(-5, 0), labels = 10^seq(-5, 0))
  axis(2)
  box()
}

for(i in 1:2){
  plot.slines(x = res[[i]]$a, y = res[[i]]$b.space, factor = paste(res[[i]]$pa, res[[i]]$m),
              color = color, ylim = c(0, 25), xlab = "dispersal", ylab = "b.space", lwd = 1.5)
  axis(1, at = seq(-5, 0), labels = 10^seq(-5, 0))
  axis(2)
  box()
}

for(i in 1:2){
  plot.slines(x = res[[i]]$a, y = res[[i]]$g.time, factor = paste(res[[i]]$pa, res[[i]]$m),
              color = color, ylim = c(0, 25), xlab = "dispersal", ylab = "g.time", lwd = 1.5)
  axis(1, at = seq(-5, 0), labels = 10^seq(-5, 0))
  axis(2)
  box()
}
for(i in 1:2){
  plot.slines(x = res[[i]]$a, y = res[[i]]$b.time, factor = paste(res[[i]]$pa, res[[i]]$m),
              color = color, ylim = c(0, 25), xlab = "dispersal", ylab = "b.time", lwd = 1.5)
  axis(1, at = seq(-5, 0), labels = 10^seq(-5, 0))
  axis(2)
  box()
}


#----------------------------------------------------------------------------------------------------
# environment (480, 550)
par(mfrow = c(2, 1), oma = c(2,2,0.5,0.5), mar = c(2,2,0.5,0.5))
for(i in 1:2){
  plot.slines(x = res[[i]]$a, y = res[[i]]$R2, factor = paste(res[[i]]$pa, res[[i]]$m),
              color = color, ylim = c(0, 1), xlab = "dispersal", ylab = "R2", lwd = 1.5)
  axis(1, at = seq(-5, 0), labels = 10^seq(-5, 0))
  axis(2)
  box()
}
