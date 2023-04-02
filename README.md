# GenerationTimeMetacommunities
Simulating metacommunities with stage-structured models

Unpublished study

Authors: Andreu Castillo Escrivà and Francesc Mesquita-Joanes

Contact: acastilloescriva@gmail.com

Metacommunities, as sets of communities linked by dispersal, are structured by a complex combination of ecological processes, such as environmental responses, biotic interactions and dispersal. However, some time-related processes have been usually disregarded in these systems, even when they are considered a fundamental component of change in the population dynamics, such as generation time (i.e., the average time between the birth of an individual and the birth of its offspring). With this code, we analyzed the role of generation time in metacommunities of competing species, by means of simulations with stage-structured models. Stage-structured models might contribute to a better understanding of the spatio-temporal variation of biodiversity, assessing temporal processes in a metacommunity context.

## Files

Settings and landscapes
- startLand.R (create landscapes with one environmental variable for the simulations; random fields)
- startSim.R (set metacommunity paramenters and prepare all the settings for the simulations)

Run simulations
- sim.jl (run simulations in julia)

Analyze simulation output
- cutSim.R (reduce the size of the simulation output for analyzing data)
- anaSim.R (analyze metacommunities; estimate community stability, diversity and environmentl fitting)
- resSim.R (plots)

Useful functions
- functionSim.R (additional R functions used in the other codes)
