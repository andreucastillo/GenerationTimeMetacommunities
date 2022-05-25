#-------------------------------------------------------------------------------------------------
# Title: Analyses of simulations
# Project: Generation time in metacommunities
# Date: 02-12-2021
#-------------------------------------------------------------------------------------------------
rm(list = ls())
source("functionSim.R")
library(mgcv)
library(vegan)

path <- "metacom_ceq_ni10_pa09"
metacom <- paste("metacom/", path, sep = "")
landscape <- "metacom/landscape"

list.land <- list.files(path = landscape, pattern = "*land", recursive = TRUE, full.names = TRUE)
list.sp <- list.files(path = metacom, pattern = "*sp", recursive = TRUE, full.names = TRUE)
list.trait <- list.files(path = metacom, pattern = "*trait", recursive = TRUE, full.names = TRUE)

land <- as.list(1:length(list.land))
for(i in 1:length(land)){
  land[[i]] <- read.table(list.land[i], header = TRUE, sep = ",")
}  

sp <- as.list(1:length(list.sp))
trait <- as.list(1:length(list.trait))

for(i in 1:length(sp)){
  sp[[i]] <- read.table(list.sp[i], header = TRUE, sep = ",")
  trait[[i]] <- read.table(list.trait[i], header = TRUE, sep = ",")
  if(max(trait[[i]]$stage) != 1){
    sp[[i]] <- join.sp(sp[[i]], trait[[i]]$species, trait[[i]]$stage)  
  }
}
rm(list.sp, list.land, list.trait, landscape)

colfunc <- colorRampPalette(c("darkblue", "purple"))
color <- rev(colfunc(10))
palette(color)
par(mar = c(4,4,2,2))

#--------------------------------------------------------------------------------------------------
# Results table
res <- as.data.frame(matrix(0, nrow = length(sp), ncol = 2))
colnames(res) <- c("m", "a")
for(i in 1:length(sp)){
  res$m[i] <- trait[[i]]$m[1]
  res$a[i] <- max(trait[[i]]$a)
}

res$m[res$m == 0] <- 1
res$rep <- 1:15

#---------------------------------------------------------------------------------------------------
# biodiversity
for(i in 1:nrow(res)){
  res$alpha[i] <- mean(specnumber(sp[[i]]))
  res$g.space[i] <- mean(specnumber(sp[[i]], land[[1]]$time))
  res$g.time[i] <- mean(specnumber(sp[[i]], land[[1]]$site))
}
res$b.space <- res$g.space - res$alpha
res$b.time <- res$g.time - res$alpha

plot.slines(x = res$a, y = res$alpha, factor = res$m,
           color = color, xlab = "Dispersal", ylab = "alpha", lwd = 2)
legend("topleft", legend = unique(res$m), col = color, lty = 1, lwd = 2, cex = 0.75)
plot.slines(x = res$a, y = res$g.space, factor = res$m,
           color = color, ylim = c(0,25), xlab = "Dispersal", ylab = "g.space", lwd = 2)
plot.slines(x = res$a, y = res$g.time, factor = res$m,
           color = color, ylim = c(0,25), xlab = "Dispersal", ylab = "g.time", lwd = 2)
plot.slines(x = res$a, y = res$b.space, factor = res$m,
           color = color, xlab = "Dispersal", ylab = "b.space", lwd = 2)
plot.slines(x = res$a, y = res$b.time, factor = res$m,
           color = color, xlab = "Dispersal", ylab = "b.time", lwd = 2)

#----------------------------------------------------------------------------------------------------
# Stability
for(i in 1:nrow(res)){
  n.space <- tapply(rowSums(sp[[i]]), land[[1]]$time, sum)
  res$N[i] <- mean(rowSums(sp[[i]]))
  res$a.stab[i] <- mean(tapply(rowSums(sp[[i]]), land[[1]]$site, sd) /
                          tapply(rowSums(sp[[i]]), land[[1]]$site, mean), na.rm = TRUE)
  res$g.stab[i] <- sd(n.space) / mean(n.space)
}

plot.slines(x = res$a, y = res$N, factor = res$m,
           color = color, xlab = "Dispersal", ylab = "N", lwd = 2)
plot.slines(x = res$a, y = res$a.stab, factor = res$m,
           color = color, xlab = "Dispersal", ylab = "a.cv", lwd = 2, ylim = c(0, 4))
plot.slines(x = res$a, y = res$g.stab, factor = res$m,
           color = color, xlab = "Dispersal", ylab = "g.cv", lwd = 2, ylim = c(0, 1))

#---------------------------------------------------------------------------------------------------
# environmental models
res$R2 <- 0
for(i in 1:nrow(res)){
  print(i)
  res$R2[i] <- gam.multi(sp[[i]], data.frame("env" = land[[res$rep[i]]]$env))  
}

plot.slines(x = res$a, y = res$R2, factor = res$m,
           color = color, xlab = "Dispersal", ylab = "R2", lwd = 2, ylim = c(0, 1))

# save results

write.csv(res, paste("results/", path, ".csv", sep = ""), row.names = FALSE)

