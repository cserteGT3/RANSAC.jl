"""
    ransac(pc, params, setenabled; reset_rand = false)

Run the RANSAC algorithm on a pointcloud with the given parameters.

Returns the candidates (as they are at the end of the iteration),
the extracted primitives
and the time it took to run the algorithm (in seconds).

# Arguments
- `pc::PointCloud`: the point cloud.
- `params::RANSACParameters`: parameters.
- `setenabled::Bool`: if `true`: set every point to enabled.
- `reset_rand::Bool=false`: if `true`, resets the random seed with `Random.seed!(1234)`
"""
function ransac(pc, params, setenabled; reset_rand = false)
    if setenabled
        pc.isenabled = trues(pc.size)
    end
    ransac(pc, params, reset_rand=reset_rand)
end

"""
    ransac(pc, params; reset_rand = false)

Run the RANSAC algorithm on a pointcloud with the given parameters.

Returns the candidates (as they are at the end of the iteration),
the extracted primitives
and the time it took to run the algorithm (in seconds).

# Arguments
- `pc::PointCloud`: the point cloud.
- `params::RANSACParameters`: parameters.
- `reset_rand::Bool=false`: if `true`, resets the random seed with `Random.seed!(1234)`
"""
function ransac(pc, params; reset_rand = false)
    reset_rand && Random.seed!(1234)

    @unpack drawN, minsubsetN, prob_det, τ = params
    @unpack itermax, shape_types = params
    @unpack extract_s, terminate_s = params
    start_time = time_ns()

    # build an octree
    subsetN = length(pc.subsets)
    @logmsg IterLow1 "Building octree."
    minV, maxV = findAABB(pc.vertices)
    octree = Cell(SVector{3}(minV), SVector{3}(maxV), OctreeNode(pc, collect(1:pc.size), 1))
    r = OctreeRefinery(8)
    adaptivesampling!(octree, r)
    @logmsg IterLow1 "Octree finished."
    # initialize levelweight vector to 1/d
    # TODO: check if levelweight is not empty
    fill!(pc.levelweight, 1/length(pc.levelweight))
    fill!(pc.levelscore, zero(eltype(pc.levelscore)))

    random_points = randperm(pc.size)
    
    # FittedShape[] is abstract array, which is bad
    candidates = FittedShape[]
    scoredshapes = ShapeCandidate[]
    extracted = ShapeCandidate[]
    
    # save octree level of a FittedShape
    shape_octree_level = Int[]
    
    # smallest distance in the pointcloud
    #lsd = smallestdistance(pc.vertices)
    
    # allocate for the random selected points
    sd = Vector{Int}(undef, drawN)
    @logmsg IterInf "Iteration begins."

    # for logging
    notifit = itermax > 10 ? div(itermax,10) : 1

    # track the number of candidates for the probabilities
    # 1. :lengthC - number of candidates in a given iteration
    # 2. :allcand - number of candidates that have been ever scored
    # 3. :nofminset - number of minimal sets that have been drawn so far
    countcandidates = [0, 0, 0]

    # iterate begin
    for k in 1:itermax
        if count(pc.isenabled) < τ
            @logmsg IterInf "Break at $k it., because left only: $(count(pc.isenabled))"
            break
        end
        # generate minsubsetN candidate
        for i in 1:minsubsetN
            #TODO: that is unsafe, but probably won't interate over the whole pc
            # select a random point
            if length(random_points)<10
                random_points = randperm(pc.size)
                @logmsg IterLow1 "Recomputing randperm."
            end
            r1 = popfirst!(random_points)

            #TODO: helyettesíteni valami okosabbal,
            # pl mindig az első enabled - ha már nagyon sok ki van véve,
            # akkor az gyorsabb lesz
            while ! pc.isenabled[r1]
                r1 = rand(1:pc.size)
            end
            # search for r1 in octree
            current_leaf = findleaf(octree, pc.vertices[r1])
            # get all the parents
            cs = getcellandparents(current_leaf)
            # revese the order, cause it's easier to map with levelweight
            reverse!(cs)
            # chosse the level with the highest score
            # if multiple maximum, the first=largest cell will be used
            curr_level = argmax(pc.levelweight[1:length(cs)])
            #an indexer array for random indexing
            cell_ind = cs[curr_level].data.incellpoints
            # frome the above, those that are enabled
            bool_cell_ind = @view pc.isenabled[cell_ind]
            enabled_inds = cell_ind[bool_cell_ind]
            # if there's less enabled vertice than needed, skip the rest
            size(enabled_inds, 1) < drawN && continue
            sd[1] = r1

            # made up heuristic
            if size(enabled_inds, 1) < 20*drawN
                # if there's few randoms, then just choose the first ones
                # random index
                sel = 0
                # sd index
                seli = 2
                while true
                    sel += 1
                    # don't choose the first one
                    enabled_inds[sel] == sd[1] && continue
                    sd[seli] = enabled_inds[sel]
                    seli += 1
                    seli == drawN+1 && break
                end
                route = 1
            else
                # if there's enough random, randomly select
                for idk in 2:drawN
                    nexti = rand(1:size(enabled_inds, 1))
                    if sd[1] == enabled_inds[nexti]
                        # try oncemore
                        nexti = rand(1:size(enabled_inds, 1))
                    end
                    #TODO: check if the same
                    #TODO: do something with it
                    sd[idk] = enabled_inds[nexti]
                end
                route = 2
            end

            if !allisdifferent(sd)
                @logmsg IterLow1 "Selected indexes have same element: $sd; route $route was taken."
                continue
            end

            # sd: indexes of the actually selected points
            #TODO: this should be something more general
            # fit plane to the selected points
            f_v = @view pc.vertices[sd]
            f_n = @view pc.normals[sd]

            forcefitshapes!(pc, f_v, f_n, params, candidates, shape_octree_level, curr_level)
        end # for t

        # evaluate the compatible points, currently used as score
        # TODO: do something with octree levels and scores
        which_ = 1
        
        # update number of candidates
        countcandidates[2] += size(candidates,1)
        scorecandidates!(pc, scoredshapes, candidates, which_, params, shape_octree_level)
        # by now every candidate is scored into scoredshapes


        # considered so far k*minsubsetN pcs. minimal sets
        countcandidates[3] = k*minsubsetN
        # length of the candidate array

        countcandidates[1] = size(scoredshapes, 1)
        if ! (size(scoredshapes, 1) < 1)
            # search for the largest score == length(inpoints) (for now)
            best_shape = largestshape(scoredshapes)
            best = findhighestscore(scoredshapes)
            if best.index != best_shape.index
                @logmsg IterLow1 "best: $(best.index), best_shape: $(best_shape.index)"
            end
            bestshape = scoredshapes[best.index]
            # TODO: refine if best.overlap
            scr = E(bestshape.score)
            best_length = size(bestshape.inpoints, 1)

            s = chooseS(countcandidates, extract_s)
            ppp = prob(best_length*subsetN, s, pc.size, drawN)

            #info printing: there's a best
            if k%notifit == 0
                @logmsg IterInf "$k. it, best: $best_length db, score: $scr, prob: $ppp, scored shapes: $(length(scoredshapes)) pcs."
            end

            # if the probability is large enough, extract the shape
            if ppp > prob_det

                # TODO: proper refit, not only getting the points that fit to that shape
                # what do you mean by refit?
                # refit on the whole pointcloud
                if bestshape.shape isa FittedPlane
                    bs = refitplane(bestshape, pc, params)
                elseif bestshape.shape isa FittedSphere
                    bs = refitsphere(bestshape, pc, params)
                elseif bestshape.shape isa FittedCylinder
                    bs = refitcylinder(bestshape, pc, params)
                elseif bestshape.shape isa FittedCone
                    bs = refitcone(bestshape, pc, params)
                end

                scs = size(bs.inpoints,1)
                @logmsg IterInf "Extracting best: $(strt(bs.candidate.shape)) score: $scr, refit length: $scs"
                # invalidate points
                for a in bs.inpoints
                    pc.isenabled[a] = false
                end
                # extract the shape and delete from scoredshapes
                push!(extracted, bs)
                deleteat!(scoredshapes, best.index)
                # mark scoredshapes that have invalid points
                toremove = Int[]
                for i in eachindex(scoredshapes)
                    for a in scoredshapes[i].inpoints
                        if ! pc.isenabled[a]
                            push!(toremove, i)
                            break
                        end
                    end
                end

                # remove scoredshapes that have invalid points
                deleteat!(scoredshapes, toremove)
            end # if extract shape
        else
            #info printing: not best
            if k%notifit == 0
                @logmsg IterInf "$k. it, no candidates."
            end
        end
        # update octree levels
        updatelevelweight(pc)

        # check exit condition
        s = chooseS(countcandidates, terminate_s)
        #if prob(τ/subsetN, s, pc.size, drawN) > prob_det
        if prob(τ, s, pc.size, drawN) > prob_det
            @logmsg IterInf "Break, at this point all shapes should be extracted: $k. iteration."
            break
        end
    end # iterate end
    fint = trunc((time_ns() - start_time)/1_000_000_000, digits=2)
    @logmsg IterInf "Iteration finished in $fint seconds with $(length(extracted)) extracted and $(length(scoredshapes)) scored shapes."
    return scoredshapes, extracted, fint
end # ransac function

"""
    rerunleftover!(pc, nofs, params, sofarextr; reset_rand=true)

For example: `rerunleftover(pcr, 4, p, extr, reset_rand=true)`.
It modifies the list of extracted candidates (`sofarextr`) and `pointcloud.isenabled`.
"""
function rerunleftover!(pc, nofs, params, sofarextr; reset_rand=true)
    lftvr_i = collect(1:pc.size)[pc.isenabled]
    @logmsg IterInf "Rerunning ransac."
    newpc = PointCloud(pc.vertices[lftvr_i], pc.normals[lftvr_i], nofs)
    _, extr2, rtime2 = ransac(newpc, params, true; reset_rand=reset_rand);
    @warn "Additional $rtime2 sec. must be added!!!"
    for i in eachindex(extr2)
        scored = extr2[i]
        old_indexes = lftvr_i[scored.inpoints]
        scored.inpoints = old_indexes
    end

    # set isenabled to false
    for j in eachindex(extr2)
        for i in extr2[j].inpoints
            if pc.isenabled[i]
                pc.isenabled[i] = false
            else
                @warn "$j at $i made mistake."
            end
        end
    end
    append!(sofarextr, extr2)
    return sofarextr, rtime2, extr2
end
