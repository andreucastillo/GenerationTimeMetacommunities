#-------------------------------------------------------------------------------------------------
# Title: Cut matrices from the simulations
# Project: Life-history traits in metacommunities
# Date: 02-12-2021
#-------------------------------------------------------------------------------------------------
rm(list = ls())

full <- "metacom_full/metacom_ceq_ni10_pa09" # path with the simulations
cut <- "metacom/metacom_ceq_ni10_pa09" # path for arranged simulations
landscape <- "metacom_full/landscape"
landscape_cut <- "metacom/landscape"

# duration
duration <- 1000
lag <- 20

#-----------------------------------------------------------------------------------------------------
# cut landscapes
list.land <- list.files(path = landscape, pattern = "*land", recursive = TRUE, full.names = TRUE)
land <- read.table(list.land[1], header = TRUE, sep = ",")
# cut time in land matrix
time.sel <- seq(max(land$time)-duration, max(land$time)-1, lag)
land <- land[land$time %in% time.sel,]

#for(i in 1:length(list.land)){
#  land <- read.table(list.land[i], header = TRUE, sep = ",")
#  land <- land[land$time %in% time.sel,] 
#  write.csv(land, file = paste(landscape_cut, "/",
#                               paste(unlist(strsplit(list.land[i], split = "/"))[2:3], collapse = "_"),
#                               sep = ""), row.names = FALSE)
#}

#----------------------------------------------------------------------------------------------------
# list files
list.sp <- list.files(path = full, pattern = "*sp", recursive = TRUE, full.names = TRUE)
list.trait <- list.sp
for(i in 1:length(list.trait)){
  list.trait[i] <- paste(paste(unlist(strsplit(list.sp[i], split ="/"))[-5], collapse = "/"),
                         "trait.csv", sep = "/")
}

# cut
sp <- read.table(list.sp[1], header = FALSE, sep = ",")
trait <- read.table(list.trait[1], header = TRUE, sep = ",")

for(i in 1:length(list.sp)){
  sp <- read.table(list.sp[i], header = FALSE, sep = ",")
  sp <- sp[rownames(sp) %in% rownames(land),]
  trait <- read.table(list.trait[i], header = TRUE, sep = ",")
  # save matrices
  write.csv(sp, file = paste(cut, "/",
                             paste(unlist(strsplit(list.sp[i], split = "/"))[3:5], collapse = "_"),
                             sep = ""), row.names = FALSE)
  write.csv(trait, file = paste(cut, "/",
                               paste(unlist(strsplit(list.trait[i], split = "/"))[3:5], collapse = "_"),
                               sep = ""), row.names = FALSE)
}
