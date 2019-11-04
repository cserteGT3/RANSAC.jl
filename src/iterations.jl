function ransac(pc, params, setenabled; reset_rand = false)
    if setenabled
        pc.isenabled = trues(pc.size)
    end
    ransac(pc, params, reset_rand=reset_rand)
end

function ransac(pc, params; reset_rand = false)
    reset_rand && Random.seed!(1234)

    @unpack drawN, minsubsetN, prob_det, τ = params
    @unpack itermax, leftovers, shape_types = params

    # build an octree
    subsetN = length(pc.subsets)
    @debug "Building octree."
    minV, maxV = findAABB(pc.vertices)
    octree = Cell(SVector{3}(minV), SVector{3}(maxV), OctreeNode(pc, collect(1:pc.size), 1))
    r = OctreeRefinery(8)
    adaptivesampling!(octree, r)
    @debug "Octree finished."
    # initialize levelweight vector to 1/d
    # TODO: check if levelweight is not empty
    fill!(pc.levelweight, 1/length(pc.levelweight))
    fill!(pc.levelscore, zero(eltype(pc.levelscore)))

    random_points = randperm(pc.size)
    candidates = ShapeCandidate[]
    scoredshapes = ScoredShape[]
    extracted = ScoredShape[]
    # smallest distance in the pointcloud
    #lsd = smallestdistance(pc.vertices)
    # allocate for the random selected points
    sd = Vector{Int}(undef, drawN)
    @debug "Iteration begins."

    # for logging
    notifit = itermax > 10 ? div(itermax,10) : 1

    # iterate begin
    for k in 1:itermax
        if count(pc.isenabled) < τ
            @debug "Break at $k, because left only: $(count(pc.isenabled))"
            break
        end
        # generate minsubsetN candidate
        for i in 1:minsubsetN
            #TODO: that is unsafe, but probably won't interate over the whole pc
            # select a random point
            if length(random_points)<10
                random_points = randperm(pc.size)
                @debug "Recomputing randperm."
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
                @debug "Selected indexes have same element: $sd; route $route was taken."
                continue
            end

            # sd: indexes of the actually selected points
            #TODO: this should be something more general
            # fit plane to the selected points
            f_v = @view pc.vertices[sd]
            f_n = @view pc.normals[sd]

            forcefitshapes(f_v, f_n, params, candidates, curr_level)
        end # for t

        # evaluate the compatible points, currently used as score
        # TODO: do something with octree levels and scores
        which_ = 1

        for c in candidates
            #TODO: save the bitmmaped parameters for debug
            push!(scoredshapes, scorecandidate(pc, c, which_, params))
        end # for c
        # by now every candidate is scored into scoredshapes
        empty!(candidates)

        # considered so far k*minsubsetN db. minimal sets
        sofar = k*minsubsetN

        size(scoredshapes, 1) < 1 && @goto endofscoring
        # search for the largest score == length(inpoints) (for now)
        # best = largestshape(scoredshapes)
        best = findhighestscore(scoredshapes)
        bestshape = scoredshapes[best.index]
        # TODO: refine if best.overlap
        scr = E(bestshape.score)
        best_length = length(bestshape.inpoints)

        #ppp = prob(best_length*subsetN, length(scoredshapes), pc.size, drawN)
        ppp = prob(best_length*subsetN, sofar, pc.size, drawN)

        #info printing
        if k%notifit == 0
            @debug "$k. it, best: $best_length db, score: $scr, prob: $ppp, scored shapes: $(length(scoredshapes)) db."
        end


        #TODO: length will be only 1/numberofsubsets
        # if the probability is large enough, extract the shape
        if ppp > prob_det

            # TODO: proper refit, not only getting the points that fit to that shape
            # what do you mean by refit?
            # refit on the whole pointcloud
            if bestshape.candidate.shape isa FittedPlane
                sp = refitplane(bestshape, pc, params)
                sps = size(sp.inpoints,1)
                @debug "Extracting best: $(strt(bestshape.candidate.shape)) score: $scr, refit length: $sps"
            elseif bestshape.candidate.shape isa FittedSphere
                ss = refitsphere(bestshape, pc, params)
                sss = size(ss.inpoints,1)
                @debug "Extracting best: $(strt(bestshape.candidate.shape)) score: $scr, refit length: $sss"
            elseif bestshape.candidate.shape isa FittedCylinder
                sc = refitcylinder(bestshape, pc, params)
                scs = size(sc.inpoints,1)
                @debug "Extracting best: $(strt(bestshape.candidate.shape)) score: $scr, refit length: $scs"
            elseif bestshape.candidate.shape isa FittedCone
                sc = refitcone(bestshape, pc, params)
                scs = size(sc.inpoints,1)
                @debug "Extracting best: $(strt(bestshape.candidate.shape)) score: $scr, refit length: $scs"
            else
                @error "Whatt? panic with $(typeof(bestshape.candidate.shape))"
            end

            # invalidate points
            for a in bestshape.inpoints
                pc.isenabled[a] = false
            end
            # extract the shape and delete from scoredshapes
            push!(extracted, bestshape)
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

        @label endofscoring
        # update octree levels
        updatelevelweight(pc)

        # check exit condition
        # TODO: τ-t is le kéne osztani a subsestek számával
        #if prob(τ/subsetN, length(scoredshapes), pc.size, drawN) > prob_det
        if prob(τ/subsetN, sofar, pc.size, drawN) > prob_det
            @debug "Break, at this point all shapes should be extracted: $k. iteráció."
            break
        end
        #if mod(k,itermax/10) == 0
        #    @debug "Iteration: $k"
        #end
    end # iterate end
    @debug "Iteration finished with $(length(extracted)) extracted and $(length(scoredshapes)) scored shapes."
    return scoredshapes, extracted
end # ransac function
