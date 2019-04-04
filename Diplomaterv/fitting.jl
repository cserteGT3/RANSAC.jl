module Fitting

using StaticArrays: SVector, MVector
using LinearAlgebra: cross, ×, dot, normalize, normalize!, norm
using ZChop: zchop, zchop!

export FittedPlane, isplane
export FittedSphere, issphere

struct FittedPlane{A<:AbstractArray}
    isplane::Bool
    point::A
    normal::A
end

"""
This is "dummy" vector for type stability.
"""
const NaNVec = SVector(0.0,0.0,0.0)

"""
    isplane(p, n, alpharad, collin_threshold = 0.2)

Fit a plane to 3 or 4 points. Their normals are used to validate the fit.

A collinearity check is used to not filter out points on one line.
# Arguments:
- `alpharad::Real`: maximum difference between the normals (in radians).
- `collin_threshold::Real=0.2`: threshold for the collinearity check (lower is stricter).
"""
function isplane(p, n, alpharad, collin_threshold = 0.2)
    @assert length(p) > 2 "At least 3 point is needed."
    crossv = MVector{3}(normalize(cross(p[2]-p[1], p[3]-p[1])))
    zchop!(crossv)
    # how many point's normal must be checked
    tc = 3
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
        # we use 4 points so 4 normals must be checked
        tc = 4
    end
    # here we have the normal of the theoretical plane

    norm_ok = falses(tc)
    invnorm_ok = falses(tc)

    thr = cos(alpharad)
    for i in 1:tc
        dotp = dot(crossv, zchop!(MVector{3}(normalize(n[i]))))
        norm_ok[i] = dotp > thr
        invnorm_ok[i] = dotp < -thr
    end
    norm_ok == trues(tc) && return FittedPlane(true, p[1], SVector{3}(crossv))
    invnorm_ok == trues(tc) && return FittedPlane(true, p[1], SVector{3}(-1*crossv))
    return FittedPlane(false, NaNVec, NaNVec)
end

struct FittedSphere{A<:AbstractArray, R<:Real}
    issphere::Bool
    center::A
    radius::R
end

function fitsphere(v, n)
    n1n = normalize(n[1])
    n2n = normalize(n[2])

    if abs(dot(n1n, n2n)) > 0.98
        # parallel normals
        # fit it, if bad, will fall out later
        centerp = (v[1]+v[2])/2
        return FittedSphere(true, centerp, norm(centerp-v[1]))
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
            return FittedSphere(true, centerp, r)
        else
            # intersection
            if dot(h, k) > 0
                # point the same direction -> +
                M = v[1] + nh/nk * n1n
                return FittedSphere(true, SVector{3}(zchop!(MVector{3}(M))), norm(M-v[1]))
            else
                # point in different direction -> -
                M = v[1] - nh/nk * n1n
                return FittedSphere(true, SVector{3}(zchop!(MVector{3}(M))), norm(M-v[1]))
            end
        end
    end
end

function issphere(p, n, epsilon, alpharad)
    @assert length(p) == length(n) "Size must be the same."
    @assert length(p) > 2 "Size must be at least 3."
    sp = fitsphere(p, n)
    return sp
end

end #module
