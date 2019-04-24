# Outline of the algorithm

using Pkg
Pkg.activate()

# every include
using Revise
using LinearAlgebra
using StaticArrays
using RegionTrees
using Random
using Logging
using Makie
#using GeometryTypes: Point3f0

include("generatesamples.jl")
includet("octree.jl")
includet("utilities.jl")
includet("fitting.jl")
includet("parameterspacebitmap.jl")

using .samples
using .Octree
using .Utilities
using .Fitting
using .ParameterspaceBitmap

# inputs
vs, ns, norms4Plot = makemeanexample(true);
# (normalize surface normals if needed)
pc = PointCloud(vs, ns, 8);
α = deg2rad(20);
ϵ = 0.5;

# build an octree
minV, maxV = findAABB(vs);
octree = Cell(SVector{3}(minV), SVector{3}(maxV), OctreeNode(pc, collect(1:length(vs)), 1));
r = OctreeRefinery(8);
adaptivesampling!(octree, r);
# initialize levelweight vector to 1/d
# TODO: check if levelweight is not empty
fill!(pc.levelweight, 1/length(pc.levelweight))

random_points = randperm(pc.size);
candidates = ShapeCandidate[]
# iterate begin
for i in 1:10
    # select random points
    r1 = popfirst!(random_points)
    # search for r1 in octree
    current_leaf = findleaf(octree, pc.vertices[r1])
    # get all the parents
    cs = getcellandparents(current_leaf)
    # revese the order, cause it's easier to map with levelweight
    reverse!(cs)
    # chosse the level with the highest score
    curr_level = argmax(pc.levelweight[1:length(cs)])
    # choose 3 more from cs[curr_level].data.incellpoints
    sd = shuffle(cs[curr_level].data.incellpoints)[1:3]
    push!(sd, r1)
    # sd: 4 indexes of the actually selected points

    #TODO: do this for every shape
    # fit sphere (and plane) to them
    fp = isplane(pc.vertices[sd], pc.normals[sd], α)
    # if not plane, then continue
    # accept/decline?
    isshape(fp) || continue

    # check compatibility on the first subset
    cps, projs = compatiblesPlane(fp, pc.vertices[pc.subsets[1]], pc.normals[pc.subsets[1]], ϵ, α)
    σs1 = length(findall(cps))
    # estimate the the score for the whole point cloud
    σ_est = estimatescore(length(pc.subsets[1]), pc.size, σs1)
    push!(candidates, ShapeCandidate(fp, σ_est))

    # select the largest score
    best = findhighestscore(candidates)
    # if intervals intersect, then eval more
    best.overlap && @warn "Overlap in confidences."
    # check connected components
    bestcandidate = candidates[best.index]
    # if already bitmapped, then skip bitmapping
    # TODO: not continue but skipping bitmapping
    bestcandidate.bitmapped && continue
    # bitmap it! but on what?

    # check candidate condition: pull from array and finalize

    # update levelweight vector
    # check exit condition
end # iterate end
