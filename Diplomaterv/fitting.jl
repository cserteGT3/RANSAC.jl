module Fitting

using StaticArrays: SVector
using LinearAlgebra: cross, dot, normalize, normalize!, norm
using ZChop: zchop

export FittedPlane, isplane

struct FittedPlane{A<:AbstractArray}
    isplane::Bool
    point::A
    normal::A
end

function isplane(p, n, alpharad, collin_threshold = 0.2)
    @assert length(p) > 2 "At least 3 point is needed."
    crossv = cross(p[2]-p[1], p[3]-p[1])
    # how many point's normal must be checked
    tc = 3
    if norm(crossv) < collin_threshold
        # The points are collinear
        # solution: use 1 more point
        if length(p) < 4
            # return with false, cause no more points can be used
            return FittedPlane(false, [], [])
        end
        crossv = cross(p[2]-p[4], p[3]-p[1])
        if norm(crossv) < collin_threshold
            # return false, cause these are definitely on one line
            return FittedPlane(false, [], [])
        end
        # we use 4 points so 4 normals must be checked
        tc = 4
    end
    # here we have the normal of the theoretical plane
    calc_norm = zchop(normalize!(crossv))

    norm_ok = falses(tc)
    invnorm_ok = falses(tc)

    thr = cos(alpharad)
    for i in 1:tc
        dotp = dot(calc_norm, p[i])
        norm_ok[i] = dotp > thr
        invnorm_ok[i] = dotp < -thr
    end
    norm_ok == trues(4) && return FittedPlane(true, p1[1], calc_norm)
    invnorm_ok == trues(4) && return FittedPlane(true, -1 .*p1[1], calc_norm)
    return FittedPlane(false, [], [])
end

end #module
