"""
An abstract type that supertypes all the fitted shapes.
"""
abstract type FittedShape end

"""
Return a string that tells the type (plane, spehere, etc.) of a `FittedShape`.
"""
function strt end

"""
Return the default parameters as a `NamedTuple` for the given of a `FittedSpahe`.
Example definition for `MyShape<:FittedShape`: `defaultshapeparameters(::Type{MyShape})=(myshape=(ϵ=1),)`
"""
function defaultshapeparameters end

"""
    struct ShapeCandidate{S<:FittedShape}

Store a primitive (`ShapeCandidate`) with its score(`ConfidenceInterval`)
and the points that belong to the shape as `Vector{Int}`.
"""
struct ShapeCandidate{S<:FittedShape}
    shape::S
    score::ConfidenceInterval
    inpoints::Vector{Int}
end

Base.show(io::IO, x::ShapeCandidate{A}) where {A} =
    print(io, "Cand: (", x.shape, "), $(length(x.inpoints)) ps")

getscore(shapecandidate::ShapeCandidate) = shapecandidate.score

"""
    findhighestscore(A)

Find the largest expected value in an array of `ShapeCandidate`s.

Indicate if there's an overlap.
"""
function findhighestscore(A)
    (length(A) > 0) || return (index = 0, overlap = false)
    ind = 1
    highest = E(getscore(A[1]))
    for i in eachindex(A)
        esc = E(getscore(A[i]))
        if esc > highest
            highest = esc
            ind = i
        end
    end

    for i in eachindex(A)
        i == ind && continue
        isoverlap(getscore(A[i]), getscore(A[ind])) && return (index = ind, overlap = true)
    end
    return (index = ind, overlap = false)
end

"""
    largestshape(A)

Find the largest shape in an array of `ShapeCandidate`s.
"""
function largestshape(A)
    length(A) < 1 && return (index = 0, size = 0)
    bestscore = length(A[1].inpoints)
    bestind = 1
    for i in eachindex(A)
        length(A[i].inpoints) <= bestscore && continue
        bestscore = length(A[i].inpoints)
        bestind = i
    end
    return (index = bestind, size = bestscore)
end

function forcefitshapes!(pc, points, normals, parameters, candidates, level_array, octree_lev)
    @extract parameters.iteration : shape_types
    for s in shape_types
        fitted = fit(s, points, normals, parameters)
        fitted === nothing && continue
        push!(candidates, fitted)
        push!(level_array, octree_lev)
    end
    return nothing
end

"""
    scorecandidates!(pc, scored_cands, candidates, subsetID, params, octree_levels)

Call `scorecandidate` for all candidates.
Reset candidate, and octree level lists (`candidates` and `octree_levels` respectively).
"""
function scorecandidates!(pc, scored_cands, candidates, subsetID, params, octree_levels)
    for i in eachindex(candidates)
        sc = scorecandidate(pc, candidates[i], subsetID, params)
        pc.levelscore[octree_levels[i]] += E(sc.score)
        push!(scored_cands, sc)
    end
    empty!(candidates)
    empty!(octree_levels)
    return nothing
end

"""
    invalidate_indexes!(pc, shape)

Invalidate indexes of the to-be-extracted shape's points.
"""
function invalidate_indexes!(pc, shape)
    for a in shape.inpoints
        pc.isenabled[a] = false
    end
    return nothing
end

"""
    removeinvalidshapes!(pc, shapelist)

Remove shapes from `shapelist`, which are invalid, say have points that are not enabled.
"""
function removeinvalidshapes!(pc, shapelist)
    toremove = Int[]
    for i in eachindex(shapelist)
        for a in shapelist[i].inpoints
            if ! pc.isenabled[a]
                push!(toremove, i)
                break
            end
        end
    end
    deleteat!(shapelist, toremove)
    return nothing
end
