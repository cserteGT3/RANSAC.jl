using Pkg
Pkg.activate()

# every include
using LinearAlgebra
using StaticArrays
using Random
using Revise
using Colors
using Makie

includet("../RANSAC.jl")
using .RANSAC

## Test1
# inputs
vs, ns, norms4Plot, shape_s = examplepc3(true, all = true, mrotdeg = 0.1, vertscale = 0.3);
sc = showgeometry(vs, ns)
# sc.center = false
# Makie.save("pure3-normals.png", sc)

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
ptt = 0.8
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

## TEST 2

# inputs
vs, ns, norms4Plot, shape_s = examplepc3(true, all = true, mrotdeg = .7, vertscale = 1.1);

sc = showgeometry(vs, ns)
#sc.center = false
#Makie.save("sc.png", sc)

pcr = PointCloud(vs, ns, 32);
# plane
p_ae = (ϵ = 3, α=deg2rad(10));
# cylidner
cy_ae = (ϵ = 3, α=deg2rad(10));
# sphere
sp_ae = (ϵ = 1, α=deg2rad(5));
one_ae = AlphSilon(sp_ae, p_ae, cy_ae);
# number of minimal subsets drawed in one iteration
tt = 15;
# probability that we found shapes
ptt = 0.8
# minimum shape size
ττ = 900
# maximum number of iteration
itermax = 40000
# size of the minimal set
draws = 3
cand, extr = ransac(pcr, one_ae, tt, ptt, ττ, itermax, draws, 500, true)
#leftover = getrest(pcr);
showtype(extr)
sc = showshapes(pcr, extr)
sc = showbytype(pcr, extr)

sc.center = false
Makie.save("bytye.png", sc)
