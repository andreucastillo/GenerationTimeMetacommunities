#------------------------------------------------------------
# Title: Simulation routine (for multiple simulations)
# Project: Generation times in metacommunities
# Date: 22-12-2021
#------------------------------------------------------------

using Pkg
Pkg.activate(".")

using CSV
using DataFrames
using DelimitedFiles
using Distributions
using Printf
using StatsBase
using DataStructures
using LinearAlgebra

# simulation
pathmetacom = filter(isdir, readdir("metacom_full/metacom_ceq_ni10_pa09/"; join = true))
pathland = "metacom_full/landscape/"
rep = 1:15

# individuals at the beginning settings
rest = [10:10:100;]
P = 0.5

for scenario in 48:length(pathmetacom)
println()
println(pathmetacom[scenario])

for r in rep
print(r,"/")
# load input
land = CSV.read(pathland*@sprintf("r%02.d", r)*"/land.csv", DataFrame)
dist = readdlm(pathland*@sprintf("r%02.d", r)*"/dist.csv", ',', Float64)
trait = CSV.read(pathmetacom[scenario]*@sprintf("/r%02.d", r)*"/trait.csv", DataFrame)
compet = readdlm(pathmetacom[scenario]*@sprintf("/r%02.d", r)*"/compet.csv", ',', Float64)

# start sp matrix
S = length(unique(trait[:,1]))
M = length(unique(land[:,4]))
sp = rand(Poisson(P), M, size(trait)[1])
sptotal = zeros(Int64, size(sp))
sptotal = sptotal + sp
stage = length(unique(trait[:,2]))

# weight distance
wdist = [Array(dist[1:end .!=1,:]) for _ in 1:size(trait)[1]]
for j in 1:length(wdist)
wdist[j] = exp.(-trait[j,6] * wdist[j] .^ 2)
end

# simulation routine
for t in 0:maximum(land[:,3])
# model
# rescue
if(sum(t .== rest) == 1)
sp = sp + rand(Poisson(P), M, size(trait)[1])
end

# dispersal
a = repeat(transpose(trait[:,5]), M, 1)
emigrant = [rand(Binomial(sp[i],a[i])) for i in reshape(1:length(a), size(a))]
immigrantPool = sum(emigrant, dims = 1)
immigrant = zeros(Int64, size(emigrant))
for j in 1:size(emigrant)[2]
if(immigrantPool[j] > 0)
immigrantexp = wdist[j] * emigrant[:,j]
immigrantexp = immigrantexp ./ sum(immigrantexp)
immigrant[:,j] = collect(values(SortedDict(countmap([1:M;
wsample(1:M, immigrantexp, immigrantPool[j])])))) .-1
end
end
sp = sp + immigrant - emigrant

# invasion
sp = sp + rand(Poisson(0.0001), M, size(trait)[1])

# niche and competition
E = land[land[:,3] .== t,5]
niche = zeros(Float64, size(sp))
for n in 1:M
niche[n,:] = exp.(-((E[n] .- trait[:,3]) ./ (2.0*trait[:,4])) .^2.0)
end

if stage == 1
dens = 1 ./ (1 .+ sp * transpose(compet))
for n in 1:M
nc = niche[n,:] .* trait[:,9] .* sp[n,:] .* dens[n,:]
sp[n,:] = [rand(Poisson(lambda)) for lambda in nc]
end
end

if stage >=2
sp_join = zeros(Int64, M, S)
for s in 1:S
sp_join[:,s] = sum(sp[:,trait[:,1] .== s], dims = 2)
end

dens = 1 ./ (1 .+ sp_join * transpose(compet))

for n in 1:M
pj = trait[trait[:,2] .== 1,7]
m = trait[trait[:,2] .== 1,8]
f = trait[trait[:,2] .== 2,9] .* niche[n, trait[:,2] .== 2] .* dens[n,:]
pa = trait[trait[:,2] .== 2,10]
pm = [zeros(Float64, stage, stage) for i in 1:S, j in 1:S]

for s in 1:S
pm[s,s] = [[pj[s] * (1-m[s]), pj[s] * m[s]] [f[s], pa[s]]]
end
pm = reduce(vcat, [reduce(hcat, pm[i, :]) for i in 1:S])

nc = sp[[n],:] * transpose(pm)
sp[n,:] = [rand(Poisson(lambda)) for lambda in nc]
end
end

# join to total species matrix
sptotal = [sptotal; sp]

if t == (maximum(land[:,3])-1)
open(pathmetacom[scenario] * @sprintf("/r%02.d", r) * "/sp.csv", "w") do io
writedlm(io, sptotal, ',')
end
end

end
end

end
