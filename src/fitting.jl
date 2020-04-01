"""
An abstract type that supertypes all the fitted shapes.
"""
abstract type FittedShape end

"""
    struct ShapeCandidate{S<:FittedShape}

Store a primitive candidate and the octree level, that it is from.
"""
struct ShapeCandidate{S<:FittedShape}
    shape::S
    octree_lev::Int
end

Base.show(io::IO, x::ShapeCandidate{S}) where {S} =
    print(io, "Cand: (", x.shape, ")")

Base.show(io::IO, ::MIME"text/plain", x::ShapeCandidate{S}) where {S} =
    print(io, "ShapeCandidate{$S}\n", x, ", octree: $(x.octree_lev)")

"""
    mutable struct ScoredShape{A<:AbstractArray}

Store a primitive (`ShapeCandidate`) with its score(`ConfidenceInterval`)
and the points that belong to the shape.
"""
mutable struct ScoredShape{A<:AbstractArray}
    candidate::ShapeCandidate
    score::ConfidenceInterval
    inpoints::A
end

Base.show(io::IO, x::ScoredShape{A}) where {A} =
    print(io, "Scored: (", x.candidate, "), $(length(x.inpoints)) ps")

Base.show(io::IO, ::MIME"text/plain", x::ScoredShape{A}) where {A} =
    print(io, "ScoredShape{$A}\n", x, ", score: ($(x.score))")

getscore(scoredshape::ScoredShape) = scoredshape.score

"""
    findhighestscore(A)

Find the largest expected value in an array of `ScoredShape`s.

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

Find the largest shape in an array of `ScoredShape`s.
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

"""
This is "dummy" vector for type stability.
"""
const NaNVec = SVector(0.0,0.0,0.0)

function forcefitshapes!(pc, points, normals, parameters, candidates, octree_lev)
    @unpack shape_types = parameters
    for s in shape_types
        if s === :plane
            fitted = fitplane(points, normals, parameters)
        elseif s === :sphere
            fitted = fitsphere(points, normals, parameters)
        elseif s === :cylinder
            fitted = fitcylinder(points, normals, parameters)
        elseif s === :cone
            fitted = fitcone(points, normals, parameters)
        else
            error("$s is not recognized as valid shape type.")
        end
        fitted === nothing || push!(candidates, ShapeCandidate(fitted, octree_lev))
    end
    return candidates
end
