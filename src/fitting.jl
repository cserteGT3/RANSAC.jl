# primitive fitting API

"""
An abstract type that supertypes all the fitted shapes.
"""
abstract type FittedShape end

"""
Return the default parameters as a `NamedTuple` for the passed `FittedShape`.

# Implementation
Example definition for `MyShape<:FittedShape`:
`defaultshapeparameters(::Type{MyShape})=(myshape=(ϵ=1, α=deg2rad(5),),)`.
"""
function defaultshapeparameters end

"""
Fit a primitive shape to point-normal pairs
(passed as an array of points and an array of normals).
The type of the primitive and a parameter named tuple is also passed.
The RANSACCloud is also passed if needed for more advanced methods.
Return `nothing`, if can't fit the shape to the points.

# Implementation
Signature: `fit(::Type{MyShape}, p, n, pc, params)`.

It should return `nothing` if the fit wasn't succesfull, `MyShape` otherwise.
If the fit results multiple shapes, it can return an array of `MyShape`s as well.
"""
function fit end

"""
Compute the score of a candidate and return its score and the corresponding points.

# Implementation
Signature: `scorecandidate(pc, candidate::MyShape, subsetID, params)`.

`subsetID` is the index of the subset that is used to estimate the score.
You must count the compatible points (that are enabled) in the given subset.
You can get the points by: `pc.vertices[pc.subsets[subsetID]]` and normals similarly.
Use the [`estimatescore`](@ref) function to estimate a score,
then return a  `(score, inpoints)`,
where `inpoints` is the indexes of the points, that are counted.
"""
function scorecandidate end

"""
Refit a primitive to thw whole point cloud, say search for all the compatible points
in the point cloud.

# Implementation
Signature: `refit(s::T, pc, params) where {T<:MyShape}`.

Search all the compatible (and enabled) points in the whole point cloud.
After finding the indexes, return an [`ExtractedShape`](@ref).
"""
function refit end
        
"""
Return a string that tells the "human-readable" type
(plane, spehere, etc.) of a `FittedShape`.

# Implementation
`strt(s::MyShape) = "myshape"`
"""
function strt end

# functions to work with primitives

"""
    ExtractedShape{S<:FittedShape}

Store an extraced primitive (`FittedShape`)
and the points that belong to the shape as `Vector{Int}`.
"""
struct ExtractedShape{S<:FittedShape}
    shape::S
    inpoints::Vector{Int}
end

Base.show(io::IO, x::ExtractedShape{A}) where {A} =
    print(io, "Cand: (", x.shape, "), $(length(x.inpoints)) ps")

"""
    struct IterationCandidates

Store shapes, their score and their points in a struct of arrays.
"""
struct IterationCandidates
    shapes::Vector{FittedShape}
    scores::Vector{ConfidenceInterval}
    inpoints::Vector{Vector{Int}}
end

function IterationCandidates()
    return IterationCandidates(FittedShape[], ConfidenceInterval[], Vector{Vector{Int}}())
end

Base.length(x::IterationCandidates) = size(x.shapes, 1)

Base.show(io::IO, x::IterationCandidates) =
    print(io, "IterationCandidates: $(size(x.shapes,1)) candidates")

"""
    recordscore!(ic::IterationCandidates, shape, score, inpoints)

Record the score and inpoints of a shape to an `IterationCandidates`.
"""
function recordscore!(ic::IterationCandidates, shape, score, inpoints)
    push!(ic.shapes, shape)
    push!(ic.scores, score)
    push!(ic.inpoints, inpoints)
    return ic
end

"""
    deleteat!(ic::IterationCandidates, args...)

Forward `deleteat!` to field arrays.
"""
function deleteat!(ic::IterationCandidates, arg)
    deleteat!(ic.shapes, arg)
    deleteat!(ic.scores, arg)
    deleteat!(ic.inpoints, arg)
    return ic
end

"""
    findhighestscore(A::IterationCandidates)

Find the largest expected value in an `IterationCandidates`.

Indicate if there's an overlap.
"""
function findhighestscore(A::IterationCandidates)
    (length(A) > 0) || return (index = 0, overlap = false)
    ind = 1
    scores = A.scores
    highest = E(scores[1])
    for i in eachindex(scores)
        esc = E(scores[i])
        if esc > highest
            highest = esc
            ind = i
        end
    end

    for i in eachindex(scores)
        i == ind && continue
        isoverlap(scores[i], scores[ind]) && return (index = ind, overlap = true)
    end
    return (index = ind, overlap = false)
end

"""
    forcefitshapes!(points, normals, parameters, candidates, level_array, octree_lev, pc)

Call `fit` based on the given parameters and add valid candidates to
"""
function forcefitshapes!(points, normals, parameters, candidates, level_array, octree_lev, pc)
    @extract parameters.iteration : shape_types
    for s in shape_types
        fitted = fit(s, points, normals, pc, parameters)
        fitted === nothing && continue
        push2candidatesandlevels!(candidates, fitted, level_array, octree_lev)
    end
    return nothing
end

"""
    scorecandidates!(pc, iterationcandidates, candidates, subsetID, params, octree_levels)

Call [`scorecandidate`](@ref) for all candidates.
Reset candidate, and octree level lists (`candidates` and `octree_levels` respectively).
"""
function scorecandidates!(pc, iterationcandidates, candidates, subsetID, params, octree_levels)
    for i in eachindex(candidates)
        sc, ip = scorecandidate(pc, candidates[i], subsetID, params)
        pc.levelscore[octree_levels[i]] += E(sc)
        recordscore!(iterationcandidates, candidates[i], sc, ip)
    end
    empty!(candidates)
    empty!(octree_levels)
    return nothing
end

"""
    invalidate_indexes!(pc, indexlist)

Invalidate indexes of the to-be-extracted shape's points, based on the `indexlist`.
"""
function invalidate_indexes!(pc, indexlist)
    for a in indexlist
        pc.isenabled[a] = false
    end
    return nothing
end

"""
    removeinvalidshapes!(pc, candidates::IterationCandidates)

Remove shapes from `candidates`, which are invalid, say have points that are not enabled.
"""
function removeinvalidshapes!(pc, candidates::IterationCandidates)
    toremove = Int[]
    for i in eachindex(candidates.inpoints)
        for a in candidates.inpoints[i]
            if ! pc.isenabled[a]
                push!(toremove, i)
                break
            end
        end
    end
    deleteat!(candidates, toremove)
    return nothing
end


function samplepointcloud!(indexv, pointcloud, octree, params)
    @extract params : iter_params=iteration
    @extract iter_params : drawN minsubsetN prob_det τ

    r1 = rand(1:pointcloud.size)

    #TODO: helyettesíteni valami okosabbal,
    # pl mindig az első enabled - ha már nagyon sok ki van véve,
    # akkor az gyorsabb lesz
    while ! pointcloud.isenabled[r1]
        r1 = rand(1:pointcloud.size)
    end
    # search for r1 in octree
    current_leaf = findleaf(octree, pointcloud.vertices[r1])
    # get the depth
    max_depth = current_leaf.data.depth
    # get the best level index == best depth
    curr_level_i = argmax(pointcloud.levelweight[1:max_depth])
    # get the best cell, based on the index
    curr_cell = getnthcell(current_leaf, curr_level_i)
    # an indexer array for random indexing
    cell_ind = curr_cell.data.incellpoints
    # frome the above, those that are enabled
    bool_cell_ind = @view pointcloud.isenabled[cell_ind]
    enabled_inds = cell_ind[bool_cell_ind]
    # if there's less enabled vertice than needed, skip the rest
    size(enabled_inds, 1) < drawN && return (false, 0)
    indexv[1] = r1

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
            enabled_inds[sel] == indexv[1] && continue
            indexv[seli] = enabled_inds[sel]
            seli += 1
            seli == drawN+1 && return (false, 1)
        end
        route = 1
    else
        # if there's enough random, randomly select
        for idk in 2:drawN
            nexti = rand(1:size(enabled_inds, 1))
            if indexv[1] == enabled_inds[nexti]
                # try oncemore
                nexti = rand(1:size(enabled_inds, 1))
            end
            #TODO: check if the same
            #TODO: do something with it
            indexv[idk] = enabled_inds[nexti]
        end
        route = 2
    end

    if ! allisdifferent(indexv)
        @logmsg IterLow1 "Selected indexes have same element: $indexv; route $route was taken."
        return (false, 1)
    end
    return (true, curr_level_i)
end


function samplepointcloud2!(indexv, pointcloud, octree, params)
    @extract params : iter_params=iteration
    @extract iter_params : drawN minsubsetN prob_det τ
    @assert drawN > 1 "At least 2 point is expected in this sampling function!"

    r1 = rand(1:pointcloud.size)

    #TODO: helyettesíteni valami okosabbal,
    # pl mindig az első enabled - ha már nagyon sok ki van véve,
    # akkor az gyorsabb lesz
    while ! pointcloud.isenabled[r1]
        r1 = rand(1:pointcloud.size)
    end
    # search for r1 in octree
    current_leaf = findleaf(octree, pointcloud.vertices[r1])
    # get the depth
    max_depth = current_leaf.data.depth
    # get the best level index == best depth
    curr_level_i = argmax(pointcloud.levelweight[1:max_depth])
    # get the best cell, based on the index
    curr_cell = getnthcell(current_leaf, curr_level_i)
    # an indexer array for random indexing
    cell_ind = curr_cell.data.incellpoints
    # frome the above, those that are enabled
    bool_cell_ind = @view pointcloud.isenabled[cell_ind]
    count(bool_cell_ind) < drawN && return (false, 0)

    enabled_inds = cell_ind[bool_cell_ind]
    # if there's less enabled vertice than needed, skip the rest
    indexv[1] = r1

    # if there's enough random, randomly select
    for idk in 2:drawN
        nexti = rand(1:size(enabled_inds, 1))
        if indexv[1] == enabled_inds[nexti]
            # try oncemore
            nexti = rand(1:size(enabled_inds, 1))
        end
        #TODO: check if the same
        #TODO: do something with it
        indexv[idk] = enabled_inds[nexti]
    end

    if ! allisdifferent(indexv)
        @logmsg IterLow1 "Selected indexes have same element: $indexv."
        return (false, 1)
    end
    return (true, curr_level_i)
end

function samplepointcloud3!(indexv, pointcloud, octree, params)
    @extract params : iter_params=iteration
    @extract iter_params : drawN minsubsetN prob_det τ
    @assert drawN > 1 "At least 2 point is expected in this sampling function!"

    r1 = rand(1:pointcloud.size)

    #TODO: helyettesíteni valami okosabbal,
    # pl mindig az első enabled - ha már nagyon sok ki van véve,
    # akkor az gyorsabb lesz
    while ! pointcloud.isenabled[r1]
        r1 = rand(1:pointcloud.size)
    end
    # search for r1 in octree
    current_leaf = findleaf(octree, pointcloud.vertices[r1])
    # get the depth
    max_depth = current_leaf.data.depth
    # get the best level index == best depth
    curr_level_i = argmax(pointcloud.levelweight[1:max_depth])
    # get the best cell, based on the index
    curr_cell = getnthcell(current_leaf, curr_level_i)
    # an indexer array for random indexing
    cell_ind = curr_cell.data.incellpoints
    # frome the above, those that are enabled
    bool_cell_ind = @view pointcloud.isenabled[cell_ind]
    count(bool_cell_ind) < drawN && return (false, 0)

    enabled_inds = cell_ind[bool_cell_ind]
    # if there's less enabled vertice than needed, skip the rest
    indexv[1] = r1

    # if there's enough random, randomly select
    indexv[2:end] = rand(enabled_inds, drawN-1)

    if ! allisdifferent(indexv)
        @logmsg IterLow1 "Selected indexes have same element: $indexv."
        return (false, 1)
    end
    return (true, curr_level_i)
end

function samplepointcloud4!(indexv, pointcloud, params)
    @extract params : iter_params=iteration
    @extract iter_params : drawN minsubsetN prob_det τ
    @assert drawN > 1 "At least 2 point is expected in this sampling function!"

    r1 = rand(1:pointcloud.size)

    #TODO: helyettesíteni valami okosabbal,
    # pl mindig az első enabled - ha már nagyon sok ki van véve,
    # akkor az gyorsabb lesz
    while ! pointcloud.isenabled[r1]
        r1 = rand(1:pointcloud.size)
    end
    # search for r1 in octree
    current_leaf = findleaf(pointcloud.octree, pointcloud.vertices[r1])
    # get the depth
    max_depth = current_leaf.data.depth
    # get the best level index == best depth
    curr_level_i = argmax(pointcloud.levelweight[1:max_depth])
    # get the best cell, based on the index
    curr_cell = getnthcell(current_leaf, curr_level_i)
    # an indexer array for random indexing
    cell_ind = curr_cell.data.incellpoints
    # frome the above, those that are enabled
    enabled_inds = @view cell_ind[pointcloud.isenabled[cell_ind]]
    size(enabled_inds, 1) < drawN && return (false, 0)

    # if there's less enabled vertice than needed, skip the rest
    indexv[1] = r1

    # if there's enough random, randomly select
    for idk in 2:drawN
        nexti = rand(1:size(enabled_inds, 1))
        if indexv[1] == enabled_inds[nexti]
            # try oncemore
            nexti = rand(1:size(enabled_inds, 1))
        end
        #TODO: check if the same
        #TODO: do something with it
        indexv[idk] = enabled_inds[nexti]
    end

    if ! allisdifferent(indexv)
        @logmsg IterLow1 "Selected indexes have same element: $indexv."
        return (false, 1)
    end
    return (true, curr_level_i)
end