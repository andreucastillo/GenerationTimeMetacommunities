#-------------------------------------------------------------------------------------------------
# Title: Start simulation matrices (prepare multiple simulations)
# Project: Generation times in metacommunities
# Date: 29-12-2021
#-------------------------------------------------------------------------------------------------
rm(list = ls())
source("functionSim.R")

# parameters
competition <- "cst" # ceq or cst
breadth <- 10 # niche breadth: 10 or 0.5
dispersal <- "both" #adult, juvenile, both

m <- rev(seq(0.1,1,0.1)) # juveniles growth (maturity)
a <- c(0.00001, 0.0001, 0.001, 0.01, 0.1, 0.5) # dispersal
pa <- 0.9 # adult survival
pj <- 1 # juvenile survival
rep <- 15

mainpath <- "metacom_full"

S <- 25 # number of species
intra <- 1 * 0.05
R <- 5
w <- 0

#---------------------------------------------------------------------------------------------------
# constant parameters and table

# simulation table
res <- as.data.frame(matrix(0, nrow = length(m) * length(a) * rep, ncol = 1))
colnames(res) <- c("competition")
res$competition <- competition
res$breadth <- breadth
res$dispersal <- dispersal
res$a <- rep(rep(a, each = rep), length(m))
res$w <- w

res$stage <- rep(ifelse(m == 1, 1, 2), each = length(a) * rep)
res$m <- rep(m, each = length(a) * rep)
res$pa <- pa
res$pj <- pj
res$f <- R * ((1-res$pa) * (1 - res$pj * (1-res$m))) / (res$pj * res$m)

res$intra <- intra
res$S <- S
res$rep <- 1:rep

#---------------------------------------------------------------------------------------------------------
# simulation input

path <- paste(mainpath, "/metacom_", competition, "_ni", sub(".", "", as.character(round(breadth,1)), fixed = TRUE),
              "_pa", sub(".", "", as.character(format(round(pa,1),nsmall = 1)), fixed = TRUE), sep = "")
dir.create(path)
write.csv(res, paste(path, "/summary.csv", sep = ""), row.names = FALSE)

for(i in 1:nrow(res)){
  # create directory
  subpath <- paste(path, "/metacom_m",
                   sub(".", "", as.character(format(round(res$m[i],1), nsmall = 1)), fixed = TRUE),
                   "_a", sub(".", "", as.character(format(round(abs(log10(res$a[i]))), nsmall = 1),1), fixed = TRUE),
                   "/", sprintf("r%02d", res$rep[i]), sep = "") 
  dir.create(subpath, recursive = TRUE)
  
  # Other parameters
  opt <- runif(res$S[i], 0, 1) # optimal niche
  if(res$competition[i] == "ceq"){
    inter <- res$intra[i]  
  }
  if(res$competition[i] == "cst"){
    inter <- res$intra[i] * runif(res$S[i] * res$S[i], 0, 0.5)  
  }
  if(res$stage[i] == 1){
    a.mod <- res$a[i]
  }
  if(res$stage[i] == 2){
    if(res$dispersal[i] == "both"){
      a.mod <- res$a[i]
    }
    if(res$dispersal[i] == "juvenile"){
      a.mod <- c(res$a[i], 0)
    }
    if(res$dispersal[i] == "adult"){
      a.mod <- c(0, res$a[i])
    }
  }
  trait <- make.trait(S = res$S[i], stage = res$stage[i], opt = opt, breadth = res$breadth[i], a = a.mod,
                      w = res$w[i], pj = res$pj[i], m = res$m[i], f = res$f[i], pa = res$pa[i],
                      intra = res$intra[i], inter = inter)  
  write.csv(trait$trait, paste(subpath, "/trait.csv", sep = ""), row.names = FALSE)
  write.table(trait$compet, paste(subpath, "/compet.csv", sep = ""), sep = ",",
              row.names = FALSE, col.names = FALSE)
}
