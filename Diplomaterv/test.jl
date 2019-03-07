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

# Plane
testpS, testpN = sampleplane(nullp, n1, n2, (1,1), (15,15));
plaN = normalsforplot(testpS, testpN);
scatter(testpS[1,:],testpS[2,:],testpS[3,:])
linesegments!(plaN, color = :blue)

# Cylinder
cp, cn = samplecylinder([0,0,1], nullp, 5, 10, (100, 150));
cylN = normalsforplot(cp, cn, 0.2);
scatter(cp[1,:],cp[2,:],cp[3,:])
linesegments!(cylN,  color = :blue)

# that is slower
# meshscatter(testpS[1,:],testpS[2,:],testpS[3,:])
