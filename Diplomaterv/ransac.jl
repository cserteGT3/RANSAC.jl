# Outline of the algorithm

#valamit
println("szeva")

function ransac(pc, α, ϵ, t, usegloβ, connekey, pt, τ, itmax, drawN)
    # build an octree
    @assert drawN > 2
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
            if length(random_points)<100
                random_points = randperm(pc.size)
                @warn "Recomputing randperm."
            end
            r1 = popfirst!(random_points)

            while ! pc.isenabled[r1]
                r1 = popfirst!(random_points)
            end
            # search for r1 in octree
            current_leaf = findleaf(octree, pc.vertices[r1])
            # get all the parents
            cs = getcellandparents(current_leaf)
            # revese the order, cause it's easier to map with levelweight
            reverse!(cs)
            # chosse the level with the highest score
            curr_level = argmax(pc.levelweight[1:length(cs)])
            # choose 3 more from cs[curr_level].data.incellpoints
            sdf = shuffle(cs[curr_level].data.incellpoints)
            sd = [r1]
            while length(sd) < 3 && length(sdf) > 0
                n = popfirst!(sdf)
                # don't use the same point twice
                n == r1 && continue
                pc.isenabled[n] && push!(sd, n)
            end
            # sd: 3 indexes of the actually selected points

            # if there's no 3 points, continue to the next draw
            length(sd) < 3 && continue

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
                cp, pp = compatiblesPlane(c.shape, pc.vertices[pc.isenabled], pc.normals[pc.isenabled], ϵ, α)
                c.inpoints = ((1:pc.size)[pc.isenabled])[cp]
                c.scored = true
                pc.levelscore[c.octree_lev] = pc.levelscore[c.octree_lev] + length(c.inpoints)
            elseif isa(c.shape, FittedSphere)
                # sphere
                cpl, uo, sp = compatiblesSphere(c.shape, pc.vertices[pc.isenabled], pc.normals[pc.isenabled], ϵ, α)
                # verti: összes pont indexe, ami enabled és kompatibilis
                verti = ( (1:pc.size)[pc.isenabled] )
                underEn = uo.under .& cpl
                overEn = uo.over .& cpl
                c.inpoints = count(underEn) >= count(overEn) ? verti[underEn] : verti[overEn]
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
            bestsize = best.size
            if mod(k,itermax/5) == 0
                @info "best size: $bestsize"
            end
            # if the probability is large enough, extract the shape
            #@show prob(bestsize, length(candidates), pc.size, k=drawN)

            if prob(bestsize, length(candidates), pc.size, k=drawN) > pt
                @info "Extraction!"
                # invalidate points
                for a in candidates[best.index].inpoints
                    pc.isenabled[a] = false
                end
                # extract the shape and delete from candidates
                push!(extracted, deepcopy(candidates[best.index]))
                deleteat!(candidates, best.index)
                # mark candidates that have invalid points
                toremove = Int[]
                for i in eachindex(candidates)
                    for a in candidates[i].inpoints
                        if ! pc.isenabled[a]
                            push!(toremove, i)
                            break
                        end
                    end
                end
                # extract the shape
                # TODO: refit

                # remove candidates that have invalid points
                deleteat!(candidates, toremove)
            end # if extract shape
        end # if length(candidates)
        # update octree levels
        updatelevelweight(pc)

        # check exit condition
        if prob(τ, length(candidates), pc.size, k=drawN) > pt
            @info "Break out from iteration at: $k"
            break
        end
        if mod(k,itermax/10) == 0
            @info "Iteration: $k"
        end
    end # iterate end
    @warn "Iteration finished."
    return candidates, extracted
end # ransac function

function showcandlength(ck)
    for c in ck
        println("candidate length: $(length(c.inpoints))")
    end
end

function showshapes(s, pointcloud, candidateA)
    colA = [:blue, :black, :darkred, :green, :brown, :yellow, :orange, :lighsalmon1, :goldenrod4, :olivedrab2, :indigo]
    @assert length(candidateA) <= length(colA) "Not enough color in colorarray. Fix it manually. :/"
    for i in 1:length(candidateA)
        ind = candidateA[i].inpoints
        scatter!(s, pointcloud.vertices[ind], color = colA[i])
    end
    s
end

function showshapes(pointcloud, candidateA)
    sc = Scene()
    showshapes(sc, pointcloud, candidateA)
end
