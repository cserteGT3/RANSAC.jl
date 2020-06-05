"""
    ransac(pc, params, setenabled; reset_rand = false)

Run the RANSAC algorithm on a pointcloud with the given parameters.

Return the extracted primitives and the time it took to run the algorithm (in seconds).

# Arguments
- `pc::RANSACCloud`: the point cloud.
- `params::NamedTuple`: parameters.
- `setenabled::Bool`: if `true`: set every point to enabled.
- `reset_rand::Bool=false`: if `true`, resets the random seed with `Random.seed!(1234)`
"""
function ransac(pc, params, setenabled; reset_rand = false)
    if setenabled
        for i in eachindex(pc.isenabled)
            pc.isenabled[i] = true
        end
    end
    ransac(pc, params, reset_rand=reset_rand)
end

"""
    ransac(pc, params; reset_rand = false)

Run the RANSAC algorithm on a pointcloud with the given parameters.

Return the extracted primitives and the time it took to run the algorithm (in seconds).

# Arguments
- `pc::RANSACCloud`: the point cloud.
- `params::NamedTuple`: parameters.
- `reset_rand::Bool=false`: if `true`, resets the random seed with `Random.seed!(1234)`
"""
function ransac(pc, params; reset_rand = false)
    reset_rand && Random.seed!(1234)

    @extract params : iter_params=iteration
    @extract iter_params : drawN minsubsetN prob_det τ
    @extract iter_params : itermax shape_types
    @extract iter_params : extract_s terminate_s

    # @unpack drawN, minsubsetN, prob_det, τ = params
    # @unpack itermax, shape_types = params
    # @unpack extract_s, terminate_s = params
    start_time = time_ns()

    # FittedShape[] is abstract array, which is bad
    candidates = FittedShape[]
    scoredshapes = IterationCandidates()
    extracted = ExtractedShape[]
    
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
            @logmsg IterInf "Break at $k iteration, because left less points, then τ"
            break
        end
        # generate minsubsetN candidate
        for i in 1:minsubsetN
            sampling_res = samplepointcloud4!(sd, pc, params)
            sampling_res[1] || continue

            # sd: indexes of the actually selected points
            #TODO: this should be something more general
            # fit plane to the selected points
            f_v = @view pc.vertices[sd]
            f_n = @view pc.normals[sd]

            forcefitshapes!(f_v, f_n, params, candidates, shape_octree_level, sampling_res[2], pc)
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
        # length of the candidate struct

        countcandidates[1] = length(scoredshapes)
        if ! (length(scoredshapes) < 1)
            best = findhighestscore(scoredshapes)
            bestshape = scoredshapes.shapes[best.index]
            # TODO: refine if best.overlap
            scr = E(scoredshapes.scores[best.index])
            #best_length = size(bestshape.inpoints, 1)

            s = chooseS(countcandidates, extract_s)
            ppp = prob(scr, s, pc.size, drawN)

            #info printing: there's a best
            if k%notifit == 0
                @logmsg IterInf "$k. it, best: $scr score, prob: $ppp, scored shapes: $(length(scoredshapes)) pcs."
            end

            # if the probability is large enough, extract the shape
            if ppp > prob_det

                # TODO: proper refit, not only getting the points that fit to that shape
                # what do you mean by refit?
                # refit on the whole pointcloud
                extr_shape = refit(bestshape, pc, params)
                
                @logmsg IterInf "Extracting best: $(strt(bestshape)) score: $scr, size: $(size(extr_shape.inpoints,1)) ps."
                # invalidate indexes
                invalidate_indexes!(pc, extr_shape.inpoints)
                # extract the shape and delete from scoredshapes
                push!(extracted, extr_shape)
                deleteat!(scoredshapes, best.index)
                # remove invalid shapes
                removeinvalidshapes!(pc, scoredshapes)
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
    return extracted, fint
end # ransac function
