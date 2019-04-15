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

export arbitrary_orthogonal
export isparallel

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

Check if `v1` and `v2` are parallel in an `alpharad` angle.

The vectors considered to be normalized.
"""
function isparallel(v1, v2, alpharad)
    dot(v1, v2) > cos(alpharad)
end

end # module
