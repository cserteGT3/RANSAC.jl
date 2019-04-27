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
using ImageView
using Makie

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

#valamit

function ransac(pc, α, ϵ, t, usegloβ, connekey, pt, τ, itmax)
    # build an octree
    minV, maxV = findAABB(vs);
    octree = Cell(SVector{3}(minV), SVector{3}(maxV), OctreeNode(pc, collect(1:length(vs)), 1));
    r = OctreeRefinery(8);
    adaptivesampling!(octree, r);
    # initialize levelweight vector to 1/d
    # TODO: check if levelweight is not empty
    fill!(pc.levelweight, 1/length(pc.levelweight))
    fill!(pc.levelscore, zero(eltype(pc.levelscore)))

    random_points = randperm(pc.size);
    candidates = ShapeCandidate[]
    extracted = ShapeCandidate[]
    # smallest distance in the pointcloud
    lsd = smallestdistance(pc.vertices)
    @info "Iteration begins."
    # iterate begin
    for k in 1:itermax
        # generate t candidate
        for i in 1:t
            #TODO: that is unsafe, but probably won't interate over the whole pc
            # select a random point
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

            #TODO: this should be something more general
            # fit plane to the selected points
            fp = isplane(pc.vertices[sd], pc.normals[sd], α)
            isshape(fp) && push!(candidates, ShapeCandidate(fp, curr_level))
            # fit sphere to the selected points
            sp = issphere(pc.vertices[sd], pc.normals[sd], ϵ, α)
            isshape(sp) && push!(candidates, ShapeCandidate(sp, curr_level))
        end # for t

        # evaluate the compatible points, currently used as score
        # TODO: do something with octree levels and scores

        for c in candidates
            c.scored && continue
            #TODO: save the bitmmaped parameters for debug
            if isa(c.shape, FittedPlane)
                # plane
                cp, pp = compatiblesPlane(c.shape, pc.vertices, pc.normals, ϵ, α)
                beta = usegloβ ? lsd : smallestdistance(pp)
                bm, idm, _ = bitmapparameters(pp, cp, beta)
                c.inpoints = largestconncomp(bm, idm, connekey)
                c.scored = true
                pc.levelscore[c.octree_lev] = pc.levelscore[c.octree_lev] + length(c.inpoints)
            elseif isa(c.shape, FittedSphere)
                # sphere
                cpl, uo, sp = compatiblesSphere(c.shape, pc.vertices, pc.normals, ϵ, α)
                betau = usegloβ ? lsd : smallestdistance(sp[uo.under])
                betao = usegloβ ? lsd : smallestdistance(sp[uo.over])
                verti = 1:pc.size
                ubm, uid, _ = bitmapparameters(sp[uo.under], cpl[uo.under], betau, verti[uo.under])
                obm, oid, _ = bitmapparameters(sp[uo.over], cpl[uo.over], betao, verti[uo.over])
                upoints = largestconncomp(ubm, uid, connekey)
                opoints = largestconncomp(obm, oid, connekey)
                c.inpoints = length(upoints) >= length(opoints) ? upoints : opoints
                c.scored = true
                pc.levelscore[c.octree_lev] = pc.levelscore[c.octree_lev] + length(c.inpoints)
            else
                # currently nothing else is implemented
                @warn "What the: $c"
            end # if
        end # for c
        if length(candidates) > 0
            # search for the largest score == length(inpoints) (for now)
            best = largestshape(candidates)
            bestsize = length(candidates[best.index].inpoints)
            # if the probability is large enough, extract the shape
            @info "best size: $bestsize"
            @info "best prob: $(prob(bestsize, length(candidates), pc.size, k=4))"
            if prob(bestsize, length(candidates), pc.size, k=4) > pt
                @info "exxxxxxxtraction!!!"
                # invalidate points
                for a in candidates[best.index].inpoints
                    pc.isenabled[a] = false
                end

                # mark candidates that have invalid points
                toremove = Int[]
                for i in eachindex(candidates)
                    i == best.index && continue
                    for a in candidate[i].inpoints
                        if ! pc.isenabled[a]
                            push!(toremove, i)
                            break
                        end
                    end
                end
                # extract the shape
                # TODO: refit
                push!(extracted, candidates[best.index])
                # remove candidates that have invalid points
                # also remove the extracted shape
                push!(toremove, best.index)
                deleteat!(candidates, toremove)
            end # if extract shape
        end # if length(candidates)
        # update octree levels
        updatelevelweight(pc)

        # check exit condition
        prob(τ, length(candidates), pc.size, k=4) > pt && break
        if mod(k,itermax/10) == 0
            @info "Iteration: $k"
        end
    end # iterate end
    return candidates, extracted
end # ransac function

# inputs
vs, ns, norms4Plot = makemeanexample();
# (normalize surface normals if needed)
pcr = PointCloud(vs, ns, 8);
αα = deg2rad(20);
ϵϵ = 0.5;
# number of minimal subsets drawed in one iteration
tt = 6;
# use "global" or "local β"
usegloββ = true
# connetivitiy for connected components
connekeyy = :eight
# probability that we found shapes
ptt = 0.99
# minimum shape size
ττ = 10
# maximum number of iteration
itermax = 200
cand, extr = ransac(pcr, αα, ϵϵ, tt, usegloββ, connekeyy, ptt, ττ, itermax)

function showcandlength(ck)
    for c in ck
        println("candidate length: $(length(c.inpoints))")
    end
end

s = Scene()
