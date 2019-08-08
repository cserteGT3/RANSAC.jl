"""
An abstract type that wraps the fitted shapes.
"""
abstract type FittedShape end

struct ShapeCandidate{S<:FittedShape}
    shape::S
    octree_lev::Int
end

mutable struct ScoredShape{A<:AbstractArray}
    candidate::ShapeCandidate
    score::ConfidenceInterval
    inpoints::A
end

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
