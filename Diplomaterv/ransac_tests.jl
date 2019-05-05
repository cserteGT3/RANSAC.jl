## RANSAC
using Pkg
Pkg.activate()

# every include
using LinearAlgebra
using StaticArrays
using RegionTrees
using Random
using Logging
using Revise
using Colors
using Makie

includet("generatesamples.jl")
includet("octree.jl")
includet("utilities.jl")
includet("fitting.jl")
includet("parameterspacebitmap.jl")

using .samples
using .Octree
using .Utilities
using .Fitting
using .ParameterspaceBitmap

include("ransac.jl")

## Test1
# inputs
vs, ns, norms4Plot, shape_s = examplepc2();
# (normalize surface normals if needed)
pcr = PointCloud(vs, ns, 8);
αα = deg2rad(10);
ϵϵ = 0.5;
# number of minimal subsets drawed in one iteration
tt = 30;
# probability that we found shapes
ptt = 0.99
# minimum shape size
ττ = 900
# maximum number of iteration
itermax = 10
# size of the minimal set
draws = 3
include("ransac.jl")
cand, extr = ransac(pcr, αα, ϵϵ, tt, ptt, ττ, itermax, draws, 500, true)
leftover = getrest(pcr);

sc = showshapes(pcr, extr)
sco = scatter(vs)
m = vbox(sco, sc)
# Makie.save("plot.png", m)

scatter(pcr.vertices[leftover])
linesegments!(norms4Plot[leftover], color = :blue)


## Test2
# noised
vs2, ns2, norms4Plot2, shape_s2 = examplepc2(true, all = true, mrotdeg = 10, vertscale = 0.1);
pcr2 = PointCloud(vs2, ns2, 8);
αα = deg2rad(20);
ϵϵ = 0.5;
# number of minimal subsets drawed in one iteration
tt = 30;
# probability that we found shapes
ptt = 0.9
# minimum shape size
ττ = 500
# maximum number of iteration
itermax = 13000

cand, extr = ransac(pcr2, αα, ϵϵ, tt, ptt, ττ, itermax, 3, 500, true)

sc = showshapes(pcr2, extr)
sco = scatter(vs2)
slo = scatter(pcr.vertices[getrest(pcr2)])
m = hbox(slo, vbox(sco, sc))
Makie.save("plot.png", m)
