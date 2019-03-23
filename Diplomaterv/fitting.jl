module Fitting

using StaticArrays: SVector, MVector
using LinearAlgebra: cross, dot, normalize, normalize!, norm
using ZChop: zchop, zchop!

export FittedPlane, isplane

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

end #module
