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
includet("ConfidenceIntervals.jl")

using .samples
using .Octree
using .Utilities
using .ConfidenceIntervals: estimatescore, E
using .Fitting
using .ParameterspaceBitmap

include("ransac.jl")

## Test1
# inputs
vs, ns, norms4Plot, shape_s = examplepc3(true, all = true, vertscale = 0.1);
# (normalize surface normals if needed)
pcr = PointCloud(vs, ns, 8);
αα = deg2rad(20);
ϵϵ = 0.3;
# number of minimal subsets drawed in one iteration
tt = 15;
# probability that we found shapes
ptt = 0.9
# minimum shape size
ττ = 900
# maximum number of iteration
itermax = 20000
# size of the minimal set
draws = 3
include("ransac.jl")
cand, extr = ransac(pcr, αα, ϵϵ, tt, ptt, ττ, itermax, draws, 500, true)
leftover = getrest(pcr);
showtype(extr)
sc = showshapes(pcr, extr)

sco = scatter(vs)
scon = showgeometry(vs, ns)
m = vbox(scon, sc)
# Makie.save("plot.png", m)

scatter(pcr.vertices[leftover])
linesegments!(norms4Plot[leftover], color = :blue)

showgeometry(vs, ns)
## Test2
# noised
vs2, ns2, norms4Plot2, shape_s2 = examplepc2(true, all = true, mrotdeg = 10, vertscale = 0.1);
pcr2 = PointCloud(vs2, ns2, 8);
αα = deg2rad(20);
ϵϵ = 0.5;
# number of minimal subsets drawed in one iteration
tt = 15;
# probability that we found shapes
ptt = 0.9
# minimum shape size
ττ = 900
# maximum number of iteration
itermax = 10000

cand, extr = ransac(pcr2, αα, ϵϵ, tt, ptt, ττ, itermax, 3, 500, true)

sc2 = showshapes(pcr2, extr)
vbox(sc, sc2)
sco = scatter(vs2)
slo = scatter(pcr.vertices[getrest(pcr2)])
m = hbox(slo, vbox(sco, sc))
Makie.save("plot.png", m)


estimatescore(1100, 16000, 1100)
