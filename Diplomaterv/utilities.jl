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

using LinearAlgebra: normalize, normalize!, cross, dot, norm
using StaticArrays: SVector

export arbitrary_orthogonal
export isparallel
export findAABB
export chopzaxis
export unitdisk2square
export ConfidenceInterval, notsoconfident
export estimatescore
export smallestdistance
export prob

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
    minV = convert(Array, deepcopy(points[1]))
    maxV = convert(Array, deepcopy(points[1]))
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


"""
    unitdisk2square(p)

Transform point on unit disk to unit sqaure.

Maps a unit circle to a unit sqaure [0,1]^2.
Source: shirley1997 - A low distortion map between disk and square
"""
function unitdisk2square(p)
    r = norm(p)
    phi = atan(p[2], p[1])
    if phi < -π/4
        phi += 2*π # in range [-pi/4,7pi/4]
    end
    if phi < π/4 # region 1
        a = r
        b = phi*a/(π/4)
    elseif phi < 3*π/4 # region 2
        b = r
        a = -(phi-π/2)*b/(π/4)
    elseif phi < 5*π/4 # region 3
        a = -r
        b = (phi - π)*a/(π/4)
    else # region 4
        b = -r
        a = -(phi-3*π/2)*b/(π/4)
    end
    # TODO:
    # decision: squre: [-1,1]^2 or [0,1]^2
    convert(typeof(p), [(a+1)/2, (b+1)/2])
end


struct ConfidenceInterval{F<:Real}
    min::F
    max::F
    ConfidenceInterval(x, y) = x > y ? error("out of order") : new{typeof(x)}(x, y)
end

"""
    notsoconfident(x, y)

Return with a new `ConfidenceInterval` parameterized correctly with x and y.
"""
function notsoconfident(x, y)
    ConfidenceInterval(min(x, y), max(x, y))
end

"""
    isoverlap(i1::A, i2::B) where {A<:ConfidenceInterval, B<:ConfidenceInterval}

Test if the two `ConfidenceInterval`s are overlapping.
"""
function isoverlap(i1::A, i2::B) where {A<:ConfidenceInterval, B<:ConfidenceInterval}
    i1.min == i2.min && return true
    if i1.min < i2.min
        return i2.min <= i1.max
    else
        return isoverlap(i2, i1)
    end
end

"""
    E(x::ConfidenceInterval)

Expected value of a [ConfidenceInterval](@ref).
"""
E(x::ConfidenceInterval) = (x.min+x.max)/2

"""
    hypergeomdev(N, x, n)

Standard deviation of the hypergeometric distribution.

Just copied from the article (Efficient RANSAC).
"""
function hypergeomdev(N, x, n)
    sq = sqrt((x*n*(N-x)*(N-n))/(N-1))
    (min=(x*n+sq)/N, max=(x*n-sq)/N)
end

"""
    estimatescore(S1length, Plength, σS1)

Give an estimate score for the whole pointcloud.
"""
function estimatescore(S1length, Plength, σS1)
    gd = hypergeomdev(-2-S1length, -2-Plength, -1-σS1)
    return notsoconfident(-1-gd.min, -1-gd.max)
end

"""
    smallestdistance(points)

Find the smallest distance between the points.
"""
function smallestdistance(points)
    @assert length(points) > 1 "At least two point is needed for that."
    ld = norm(points[2]-points[1])
    for i in eachindex(points)
        for j in eachindex(points)
            if i!=j
                d = norm(points[i]-points[j])
                ld = d < ld ? d : ld
            end
        end
    end
    ld
end

"""
    prob(n, s, N; k = 4)

The probability of successful detection of a shape sized `n`,
from a point cloud size of `N`, with `k` size of minimal sets and `s` draws.

# Arguments:
- `n`: size of the shape.
- `s`: number of candidates that have been drawn.
- `N`: size of the point cloud.
- `k`: size of the minimal set required to define a shape candidate.
"""
prob(n, s, N; k = 4) = 1-(1-(n/N)^k)^s

end # module
