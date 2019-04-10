using Pkg
Pkg.activate()

using Revise
using LinearAlgebra
using StaticArrays
using Makie

include("generatesamples.jl")
using .samples

# testing sphere parameter bitmap

using AbstractPlotting

o = SVector(0.0,5,π);
R = 13;
tsP, tsN = samplesphere(o, R, (4,4));
#sc =  Scene()
#showgeometry(sc, tsP, tsN)

parSp = Array{SArray{Tuple{2},Float64,1,2},1}()

for i in eachindex(tsP)
        vec = tsP[i]-o
        fí = atan(vec[2],vec[1])
        theta = atan(vec[3], vec[1])
        push!(parSp, SVector(fí, theta))
end

#sc2 = scatter(parSp)

hbox(scatter(tsP, markersize=2), scatter(parSp))
