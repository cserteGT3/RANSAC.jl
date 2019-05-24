## RANSAC
using Pkg
Pkg.activate()

# every include
using LinearAlgebra
using StaticArrays
using Random
using Revise
using Colors
using Makie

includet("RANSAC.jl")
using .RANSAC

## Test1
# inputs
vs, ns, norms4Plot, shape_s = examplepc3();
sc = showgeometry(vs, ns)
sc.center = false
#Makie.save("pure3-normals.png", sc)
# (normalize surface normals if needed)
pcr = PointCloud(vs, ns, 32);
# plane
p_ae = (ϵ = 0.3, α=deg2rad(5));
# cylidner
cy_ae = (ϵ = 0.3, α=deg2rad(5));
# sphere
sp_ae = (ϵ = 0.3, α=deg2rad(5));
one_ae = AlphSilon(sp_ae, p_ae, cy_ae);
# number of minimal subsets drawed in one iteration
tt = 15;
# probability that we found shapes
ptt = 0.9
# minimum shape size
ττ = 900
# maximum number of iteration
itermax = 20
# size of the minimal set
draws = 3
cand, extr = ransac(pcr, one_ae, tt, ptt, ττ, itermax, draws, 500, true)
#leftover = getrest(pcr);
showtype(extr)
sc = showshapes(pcr, extr)

showbytype(pcr, extr)

m2 = vbox(sc, sc2)

sco = scatter(vs)
scon = showgeometry(vs, ns)
m = vbox(scon, sc)
# Makie.save("plot.png", m)

scatter(pcr.vertices[leftover])
linesegments!(norms4Plot[leftover], color = :blue)

showgeometry(vs, ns)
## Test2
# noised
vs2, ns2, norms4Plot2, shape_s2 = examplepc3(true, all = true, mrotdeg = 5, vertscale = 0.3);
pcr2 = PointCloud(vs2, ns2, 32);
# plane
p_ae = (ϵ = 0.8, α=deg2rad(10));
# cylidner
cy_ae = (ϵ = 0.4, α=deg2rad(5));
# sphere
sp_ae = (ϵ = 0.2, α=deg2rad(5));
two_ae = AlphSilon(sp_ae, p_ae, cy_ae);
# number of minimal subsets drawed in one iteration
tt = 15;
# probability that we found shapes
ptt = 0.9
# minimum shape size
ττ = 900
# maximum number of iteration
itermax = 40000

cand, extr = ransac(pcr2, two_ae, tt, ptt, ττ, itermax, 3, 500, true)
showtype(extr)
showbytype(pcr2, extr)
sc2 = showshapes(pcr2, extr)
vbox(sc, sc2)
sco = scatter(vs2)
slo = scatter(pcr.vertices[getrest(pcr2)])
m = hbox(slo, vbox(sco, sc))
Makie.save("plot.png", m)


estimatescore(1100, 16000, 1100)
