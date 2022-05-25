#-------------------------------------------------------------------------------------------------
# Title: Start landscapes for simulations
# Project: Generation times in metacommunities
# Date: 29-12-2021
#-------------------------------------------------------------------------------------------------
rm(list = ls())
source("functionSim.R")

path <- "metacom/landscape1"
scale <- 500
spa <- scale/10/5
rep <- 5 # number of landscapes
M <- 50
steps <- 2200

unlink(path, recursive = TRUE)
dir.create(path)

for(i in 1:rep){
  print(i)
  subpath <- paste(path, "/", sprintf("r%02d",i), sep = "")
  dir.create(subpath)
  
  land <- make.land(M = M, steps = steps, scale = scale, spa = spa, plot = TRUE)
  write.csv(land$land, paste(subpath, "/land.csv", sep = ""), row.names = FALSE)
  write.csv(land$distance, paste(subpath, "/dist.csv", sep = ""), row.names = FALSE)
}
