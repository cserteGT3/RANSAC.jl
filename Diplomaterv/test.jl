using Pkg
Pkg.activate()

using Revise
using LinearAlgebra

using Makie

includet("generatesamples.jl")

using .samples

n1 = [1,0,0]
n2 = [0,1,0.5]
nullp = zeros(3)

testpS, testpN = sampleplane(nullp, n1, n2, (15,15), (1,1));

scatter(testpS[1,:],testpS[2,:],testpS[3,:])
# that is lower
meshscatter(testpS[1,:],testpS[2,:],testpS[3,:])

# test this:
#  meshscatter(rand(10), rand(10), rand(10), color=rand(10))
