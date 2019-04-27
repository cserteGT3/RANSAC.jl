module Fitting

include("utilities.jl")

using StaticArrays: SVector, MVector
using LinearAlgebra: cross, ×, dot, normalize, normalize!, norm
using ZChop: zchop, zchop!

using .Utilities

export FittedShape, isshape
export FittedPlane, isplane
export FittedSphere, issphere
export ShapeCandidate, findhighestscore
export largestshape

"""
An abstract type that wraps the fitted shapes.
"""
abstract type FittedShape end

struct FittedPlane{A<:AbstractArray} <: FittedShape
    isplane::Bool
    point::A
    normal::A
end

function isshape(shape::FittedPlane)
    return shape.isplane
end

mutable struct ShapeCandidate{S<:FittedShape, A<:AbstractArray}
    shape::S
    score::Union{ConfidenceInterval, Nothing}
    inpoints::A
    scored::Bool
    octree_lev::Int
end

ShapeCandidate(shape, score, octlev) = ShapeCandidate(shape, score, [], false, octlev)
ShapeCandidate(shape, octlev) = ShapeCandidate(shape, nothing, [], false, octlev)

"""
    findhighestscore(A)

Find the largest expected value in an array of `ShapeCandidate`s.

Indicate if there's an overlap.
"""
function findhighestscore(A)
    (length(A) > 0) || return (index = 0, overlap = false)
    ind = 1
    highest = E(A[1].score)
    for i in eachindex(A)
        esc = E(A[i].score)
        if esc > highest
            highest = esc
            ind = i
        end
    end
    overlaps = [isoverlap(A[ind].score, a.score) for a in A]
    overlap[ind] = false
    # no overlap:
    overlaps == falses(length(A)) && return (index = ind, overlap = false)
    # overlap:
    return (index = ind, overlap = true)
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

"""
This is "dummy" vector for type stability.
"""
const NaNVec = SVector(0.0,0.0,0.0)

"""
    isplane(p, n, alpharad, collin_threshold = 0.2)

Fit a plane to 3 points. Their and additional point's normals are used to validate the fit.

A collinearity check is used to not filter out points on one line.
# Arguments:
- `alpharad::Real`: maximum difference between the normals (in radians).
- `collin_threshold::Real=0.2`: threshold for the collinearity check (lower is stricter).
"""
function isplane(p, n, alpharad, collin_threshold = 0.2)
    lp = length(p)
    @assert length(p) > 2 "At least 3 point is needed."
    @assert lp == length(n) "Size must be the same."
    crossv = MVector{3}(normalize(cross(p[2]-p[1], p[3]-p[1])))
    zchop!(crossv)
    # how many point's normal must be checked
    if norm(crossv) < collin_threshold
        # The points are collinear
        # solution: use 1 more point
        if length(p) < 4
            # return with false, cause no more points can be used
            return FittedPlane(false, NaNVec, NaNVec)
        end
        crossv = MVector{3}(normalize(cross(p[2]-p[4], p[3]-p[1])))
        if norm(crossv) < collin_threshold
            # return false, cause these are definitely on one line
            return FittedPlane(false, NaNVec, NaNVec)
        end
    end
    # here we have the normal of the theoretical plane
    norm_ok = falses(lp)
    invnorm_ok = falses(lp)

    thr = cos(alpharad)
    for i in 1:lp
        dotp = dot(crossv, zchop!(MVector{3}(normalize(n[i]))))
        norm_ok[i] = dotp > thr
        invnorm_ok[i] = dotp < -thr
    end
    norm_ok == trues(lp) && return FittedPlane(true, p[1], SVector{3}(crossv))
    invnorm_ok == trues(lp) && return FittedPlane(true, p[1], SVector{3}(-1*crossv))
    return FittedPlane(false, NaNVec, NaNVec)
end

struct FittedSphere{A<:AbstractArray, R<:Real} <: FittedShape
    issphere::Bool
    center::A
    radius::R
    outwards::Bool
end

function isshape(shape::FittedSphere)
    return shape.issphere
end

function fitsphere(v, n)
    n1n = normalize(n[1])
    n2n = normalize(n[2])

    b = false
    if abs(dot(n1n, n2n)) > 0.98
        # parallel normals
        # fit it, if bad, will fall out later
        centerp = (v[1]+v[2])/2
        return FittedSphere(true, centerp, norm(centerp-v[1]), b)
    else
        # not parallel normals
        # first check if they intersect
        g = v[2]-v[1]
        h = cross(n2n, g)
        k = cross(n2n, n1n)
        nk = norm(k)
        nh = norm(h)

        if nk < 0.02 || nh < 0.02
            # if nk == 0 || nh == 0, but it's numeric
            # no intersection
            n2 = n2n × (n1n × n2n)
            n1 = n1n × (n2n × n1n)
            c1 = v[1] + dot( (v[2]-v[1]), n2 )/dot( n[1], n2 ) * n[1]
            c2 = v[2] + dot( (v[1]-v[2]), n1 )/dot( n[2], n1 ) * n[2]
            centerp = (c1+c2)/2
            r = (norm(v[1]-centerp) + norm(v[1]-centerp))/2
            return FittedSphere(true, centerp, r, b)
        else
            # intersection
            if dot(h, k) > 0
                # point the same direction -> +
                M = v[1] + nh/nk * n1n
                return FittedSphere(true, SVector{3}(zchop!(MVector{3}(M))), norm(M-v[1]), b)
            else
                # point in different direction -> -
                M = v[1] - nh/nk * n1n
                return FittedSphere(true, SVector{3}(zchop!(MVector{3}(M))), norm(M-v[1]), b)
            end
        end
    end
end

function setsphereOuterity(sp, b)
    FittedSphere(sp.issphere, sp.center, sp.radius, b)
end

"""
    issphere(p, n, epsilon, alpharad)

Fit a sphere to 2 points. Additional points and their normals are used to validate the fit.

# Arguments:
- `epsilon::Real`: maximum distance-difference between the fitted and measured spheres.
- `alpharad::Real`: maximum difference between the normals (in radians).
"""
function issphere(p, n, epsilon, alpharad)
    pl = length(p)
    @assert pl == length(n) "Size must be the same."
    @assert pl > 2 "Size must be at least 3."
    # "forcefit" a sphere
    sp = fitsphere(p, n)
    # check if real sphere
    sp.issphere || return sp
    thr = cos(alpharad)
    vert_ok = falses(pl)
    norm_ok = falses(pl)
    invnorm_ok = falses(pl)
    for i in 1:pl
        # vertice check
        vert_ok[i] = abs(norm(p[i]-sp.center)-sp.radius) < epsilon
        # normal check
        dotp = dot( normalize(p[i]-sp.center), normalize(n[i]) )
        norm_ok[i] = dotp > thr
        invnorm_ok[i] = dotp < -thr
    end
    vert_ok == trues(pl) || return FittedSphere(false, NaNVec, 0, false)
    norm_ok == trues(pl) && return setsphereOuterity(sp, true)
    invnorm_ok == trues(pl) && return setsphereOuterity(sp, false)
    return FittedSphere(false, NaNVec, 0, false)
end

end #module
