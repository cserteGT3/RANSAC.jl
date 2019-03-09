using Pkg
Pkg.activate()

using Revise
using LinearAlgebra
using Makie

includet("generatesamples.jl")

using .samples

includet("skecher.jl")

## Easy examples

n1 = [1,0,0]
n2 = [0,1,0.5]
nullp = zeros(3)

# Plane
pV, pN = sampleplane(nullp, n1, n2, (1,1), (4,4));

noiN = noisifynormals(pN, 45);
noiV = noisifyvertices(pV, true, 0.2)

plaN = normalsforplot(pV, pN);
plaNoi = normalsforplot(noiV, noiN);

scatter(pV[1,:], pV[2,:], pV[3,:])
scatter(noiV[1,:], noiV[2,:], noiV[3,:])

linesegments!(plaNoi, color = :blue)
linesegments!(plaN, color = :blue)

# Cylinder
cP, cN = samplecylinder([0,0,1], nullp, 5, 10, (100, 150));

noiNc = noisifynormals(cN, 45);
noiVc = noisifyvertices(cP, true, 0.2);

plaNc = normalsforplot(cP, cN);
plaNcoi = normalsforplot(noiVc, noiNc);

scatter(cP[1,:], cP[2,:], cP[3,:])
scatter(noiVc[1,:], noiVc[2,:], noiVc[3,:])

linesegments!(plaNc,  color = :blue)
linesegments!(plaNcoi,  color = :blue)

## a random sample
vs, ns, normsP = makemeanexample(true)

scatter(vs[1,:], vs[2,:], vs[3,:])
linesegments!(normsP, color = :blue)
