using Pkg
Pkg.activate()

using Revise
using LinearAlgebra
using StaticArrays
using Makie

includet("generatesamples.jl")

using .samples

## Easy examples

n1 = SVector(1,0,0)
n2 = SVector(0,1,0.5)
nullp = SVector(0,0,0)

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

## Regiontree tests

includet("octree.jl")

using .Octree

using StaticArrays
using RegionTrees

r = Cell(SVector(0.0,0,0),SVector(1.0,1,1))
split!(r)
cb1 = child_boundary(r,(1,1,1))
vcb1 = vertices(cb1)

for i in 1:2
   for j in 1:2
       for k in 1:2
           @show (i,j,k)
           @show child_boundary(r,(i,j,k))
       end
   end
end

v2 = collect(1:30);
v3 = rand(Bool,10);
v33 = rand(Bool,10);

v4 = SVector(0,0);
v5 = SVector(-1,1);

using GeometryTypes: Point3f0
r = Cell(SVector(0.0,0,0),SVector(1.0,1,1))
d = OctreeNode(rand(Point3f0, 10), collect(1:10), 1.0)
octref = OctreeRefinery(1)
r2 = Cell(SVector(0.0,0,0),SVector(1.0,1,1), d)

iswithinrectangle(r2.boundary, rand(Point3f0))
map(x -> iswithinrectangle(r2.boundary, x), rand(Point3f0,10))
