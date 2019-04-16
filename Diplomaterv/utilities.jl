module RodriguesRotations

using LinearAlgebra
using ZChop: zchop

export rodriguesdeg, rodriguesrad

"""
    rodrigues(nv, Θ)

Create a rotation matrix from a normalized axis and an angle (in radian).
Near-zero elements will be chopped to zero.
"""
function rodrigues(nv, Θ)
    et = eltype(nv)
    R = zeros(et,3,3)
    R = nv*nv' + cos(Θ).*(Matrix{et}(I, 3,3) - nv*nv') + sin(Θ).*crossprodtensor(nv)
    return R
end

"""
    crossprodtensor(v)

Create the cross product tensor of 3 dimensional vector.
"""
function crossprodtensor(v)
    [0 -v[3] v[2];
    v[3] 0 -v[1];
    -v[2] v[1] 0]
end

"""
    rodriguesrad(nv, ϑ)

Create a rotation matrix from an axis (`nv`) and an angle (`ϑ`) in radians.
"""
function rodriguesrad(nv, ϑ)
    nvn = normalize(nv)
    return rodrigues(nvn, ϑ)
end

"""
    rodriguesdeg(nv, ϑ)

Create a rotation matrix from an axis (`nv`) and an angle (`ϑ`) in degrees.
"""
function rodriguesdeg(nv, ϑ)
    nvn = normalize(nv)
    return rodrigues(nvn, deg2rad(ϑ))
end

end

"""
Collection of useful functions.
"""
module Utilities

using LinearAlgebra: normalize, normalize!, cross, dot
using StaticArrays: SVector

export arbitrary_orthogonal
export isparallel
export findAABB
export chopzaxis

"""
    arbitrary_orthogonal(vec)

Create an arbitrary orthogonal vector to `vec`.
"""
function arbitrary_orthogonal(vec)
    @assert size(vec,1) == 3 "Implemented only for 3 dimensional."
    v = normalize(vec)
    b0 = (v[1] <  v[2]) && (v[1] <  v[3]);
    b1 = (v[2] <= v[1]) && (v[2] <  v[3]);
    b2 = (v[3] <= v[1]) && (v[3] <= v[2]);
    rv = normalize!([convert(Int, b0), convert(Int, b1), convert(Int, b2)])
    return convert(typeof(vec), cross(v, rv))
end


"""
    isparallel(v1, v2, alpharad)

Check if `v1` and `v2` are parallel in an `alpharad` angle and point towards the same direction.

The vectors considered to be normalized.
"""
function isparallel(v1, v2, alpharad)
    dot(v1, v2) > cos(alpharad)
end


"""
    findAABB(points)

Find the axis aligned minimum bounding box of the given pointcloud.
"""
function findAABB(points)
    minV = convert(Array, points[1])
    maxV = convert(Array, points[1])
    for i in 1:length(points)
        a = points[i]
        for j in 1:length(minV)
            minV[j] = minV[j] > a[j] ? a[j] : minV[j]
            maxV[j] = maxV[j] < a[j] ? a[j] : maxV[j]
        end
    end
    return convert(eltype(points), minV), convert(eltype(points), maxV)
end


"""
    chopzaxis(points)

Chop the 3rd element of every point.
"""
function chopzaxis(points)
    @assert length(points[1]) == 3 "Implemented only for 3 long vectors."
    [ SVector(a[1], a[2]) for a in points]
end

end # module
