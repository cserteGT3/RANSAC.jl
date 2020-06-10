"""
    rodrigues(nv, Θ)

Create a rotation matrix from a normalized axis and an angle (in radian).
"""
function rodrigues(nv, Θ)
    et = eltype(nv)
    #R = nv*nv' + cos(Θ).*(Matrix{et}(I, 3,3) - nv*nv') + sin(Θ).*crossprodtensor(nv)
    R = nv*nv' + cos(Θ).*(Matrix(I, 3,3) - nv*nv')
    pluscrossprod!(R, sin(Θ), nv)
    return R
end

"""
    rodrigues(nv, Θ)

Create a rotation matrix from a normalized axis and an angle (in radian).
"""
function rodrigues(nv::SA, Θ) where {SA<:StaticArray}
    #R = nv*nv' + cos(Θ).*(Matrix{et}(I, 3,3) - nv*nv') + sin(Θ).*crossprodtensor(nv)
    R = MMatrix{3,3}(nv*nv' + cos(Θ).*(Matrix(I, 3,3) - nv*nv'))
    pluscrossprod!(R, sin(Θ), nv)
    return SMatrix{3,3}(R)
end

"""
    pluscrossprod!(A, value, v)

Add a crossproduct tensor to `A`.
Equals to: `A+value .* crossprodtensor(v)`
"""
function pluscrossprod!(A, value, v)
    A[1,2] -= value*v[3]
    A[1,3] += value*v[2]
    
    A[2,1] += value*v[3]
    A[2,3] -= value*v[1]

    A[3,1] -= value*v[2]
    A[3,2] += value*v[1]

    return A
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

## Collection of useful functions.


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
    arbitrary_orthogonal2(vec)

Create an arbitrary orthogonal vector to `vec`.
"""
function arbitrary_orthogonal2(vec)
    @assert size(vec,1) == 3 "Implemented only for 3 dimensional."
    v = normalize(vec)
    rv = normalize(rand(3))
    return convert(typeof(vec), cross(v, rv))
end


"""
    isparallel(v1, v2, alpharad)

Check if `v1` and `v2` are parallel in an `alpharad` angle
    and point towards the same direction.

The vectors considered to be normalized.
"""
function isparallel(v1, v2, alpharad)
    dot(v1, v2) > cos(alpharad)
end


"""
    findAABB(points)

Find the axis aligned minimum bounding box of the given points.
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
    [SVector(a[1], a[2]) for a in points]
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
    minmaxdistance(points)

Find the smallest, largest and average distance between the points.
"""
function minmaxdistance(points)
    @assert length(points) > 1 "At least two point is needed for that."
    ld = norm(points[2]-points[1])
    maxid = norm(points[2]-points[1])
    sumd = 0.
    sumi = 0
    for i in eachindex(points)
        for j in eachindex(points)
            if i!=j
                d = norm(points[i]-points[j])
                sumd += d
                sumi += 1
                ld = d < ld ? d : ld
                maxid = d > maxid ? d : maxid
            end
        end
    end
    avg_d = sumd/sumi
    return (mind=ld, maxd=maxid, avgd=avg_d)
end

"""
    not0minmaxdistance(points)

Find the smallest and largest distance between the points.
Zero distance is not taken into account.
"""
function not0minmaxdistance(points)
    @assert length(points) > 1 "At least two point is needed for that."
    ld = norm(points[2]-points[1])
    maxid = norm(points[2]-points[1])
    for i in eachindex(points)
        for j in eachindex(points)
            if i!=j
                d = norm(points[i]-points[j])
                d == 0 && continue
                ld = d < ld ? d : ld
                maxid = d > maxid ? d : maxid
            end
        end
    end
    return (mind=ld, maxd=maxid)
end

"""
    prob(n, s, N, k)

The probability of successful detection of a shape sized `n`,
from a point cloud size of `N`, with `k` size of minimal sets and `s` draws.

# Arguments:
- `n`: size of the shape.
- `s`: number of candidates that have been drawn.
- `N`: size of the point cloud.
- `k`: size of the minimal set required to define a shape candidate.
"""
prob(n, s, N, k) = 1-(1-(n/N)^k)^s


"""
    havesameelement(A, B)

Return true if `A` and `B` have at least one common element.
"""
function havesameelement(A, B)
    for a in A
        for b in B
            a == b && return true
        end
    end
    false
end

"""
    allisdifferent(a::AbstractArray{T, 1}) where {T}

Check if all the elements of an array are different.
Supports only one dimensional arrays.
"""
function allisdifferent(a::AbstractArray{T, 1}) where {T}
    #TODO: make it more abstract with:
    #https://docs.julialang.org/en/v1/base/iterators/
    size(a, 1) < 2 && return true
    for i in 2:size(a, 1)
        for j in 1:i-1
            a[i] == a[j] && return false
        end
    end
    return true
end

function chooseS(A, k)
    symtoint = Dict(:lengthC=>1, :allcand=>2, :nofminset=>3)
    return A[symtoint[k]]
end

"""
    defaultiterationparameters(shape_types)

Construct a named tuple with the default iteration parameters.
`shape_types` is an array of `FittedShape`s,
that controls which primitives you want to fit to the point cloud.

# Examples

```julia-repl
julia> RANSAC.defaultiterationparameters([FittedPlane])
(iteration = (drawN = 3, minsubsetN = 15, prob_det = 0.9,
shape_types = UnionAll[FittedPlane], τ = 900, itermax = 1000,
extract_s = :nofminset, terminate_s = :nofminset),)

julia> RANSAC.defaultiterationparameters([FittedPlane, FittedSphere, FittedCone])
(iteration = (drawN = 3, minsubsetN = 15, prob_det = 0.9,
shape_types = UnionAll[FittedPlane, FittedSphere, FittedCone], τ = 900, itermax = 1000,
extract_s = :nofminset, terminate_s = :nofminset),)
```

# Implementation
- `drawN`: number of points to be sampled (length of a minimal subset).
- `minsubsetN`: number of minimal sets sampled in one iteration.
- `prob_det`: probability of detection.
- `τ`: minimal shape size.
- `itermax`: maximum number of iteration.
- `shape_types`: shapes that are fitted to the point cloud (array of types).
- `extract_s`, `terminate_s`: they are for easier testing, do not delete or modify it.
"""
function defaultiterationparameters(shape_types)
    # `drawN`: number of points to be sampled (length of a minimal subset)
    # `minsubsetN`: number of minimal sets sampled in one iteration
    # `prob_det`: probability of detection
    # `τ`: minimal shape size
    # `itermax`: maximum number of iteration
    # `shape_types`: shapes that are fitted to the point cloud (array of types)
    # `extract_s`, `terminate_s`:
        # track the number of candidates for the probabilities
        # which to use to terminate the algorithm
    # 1. :lengthC - number of candidates in a given iteration
    # 2. :allcand - number of candidates that have been ever scored
    # 3. :nofminset - number of minimal sets that have been drawn so far
    ip = (drawN=3, minsubsetN=15, prob_det=0.9, shape_types=shape_types, τ=900, itermax=1000, extract_s=:nofminset, terminate_s=:nofminset)
    return (iteration=ip,)
end

"""
    defaultcommonparameters()

Construct a `NamedTuple` with the default common parameters.

# Examples
```julia-repl
julia> defaultcommonparameters()
(common = (collin_threshold = 0.2, parallelthrdeg = 1.0),)
```

# Implementation
This section describes the role of the common parameters.
- `collin_threshold`: 3 points can be nearly collinear,
    in some cases they must be filtered. See the code of:
    ` fit(::Type{FittedPlane}, p, n, params)`.
- `parallelthrdeg`: threshold for two vectos being parallel, in degrees.
    If `abs(dot(a,b))>cosd(parallelthrdeg)`, `a` and `b` are considered to be parallel.
"""
function defaultcommonparameters()
    #`collin_threshold`: threshold of points being collinear
    # `parallelthrdeg`: threshold of two vectors being parallel (in degrees)
    cp = (collin_threshold = 0.2, parallelthrdeg = 1.)
    return (common=cp,)
end

"""
    defaultparameters(shape_types::Vector{T}) where {T}

Construct a `NamedTuple` with the given shape types and the default parameters.

# Examples
```julia-repl
julia> defaultparameters([FittedSphere, FittedPlane])
(iteration = (drawN = 3, minsubsetN = 15, prob_det = 0.9,
shape_types = UnionAll[FittedSphere, FittedPlane], τ = 900, itermax = 1000,
extract_s = :nofminset, terminate_s = :nofminset),
common = (collin_threshold = 0.2, parallelthrdeg = 1.0),
sphere = (ϵ = 0.3, α = 0.08726646259971647, sphere_par = 0.02),
plane = (ϵ = 0.3, α = 0.08726646259971647))
```
"""
function defaultparameters(shape_types::Vector{T}) where {T}
    def_iter_pars = defaultiterationparameters(shape_types)
    def_iter_pars = merge(def_iter_pars, defaultcommonparameters())
    def_pars = defaultshapeparameters.(shape_types)
    for a in def_pars
        def_iter_pars = merge(def_iter_pars, a)
    end
    def_iter_pars
end

"""
    ransacparameters(p::T=DEFAULT_PARAMETERS; kwargs...) where {T<:NamedTuple}

Construct a `NamedTuple` based on a previous one,
defaulting to `DEFAULT_PARAMETERS` and override it with the kwargs.
Check the docs and examples for more.

# Examples
```julia-repl
julia> p1 = ransacparameters()
(iteration = (drawN = 3, minsubsetN = 15, prob_det = 0.9),
plane = (ϵ = 0.3, α = 0.08726646259971647),
cone = (ϵ = 0.3, α = 0.08726646259971647, minconeopang = 0.03490658503988659),
cylinder = (ϵ = 0.3, α = 0.08726646259971647),
sphere = (ϵ = 0.3, α = 0.08726646259971647, sphere_par = 0.1))

julia> p2 = ransacparameters(p1; sphere=(ϵ=0.9, α=deg2rad(1),), plane=(ϵ=1.0,))
(iteration = (drawN = 3, minsubsetN = 15, prob_det = 0.9),
plane = (ϵ = 1.0, α = 0.08726646259971647),
cone = (ϵ = 0.3, α = 0.08726646259971647, minconeopang = 0.03490658503988659),
cylinder = (ϵ = 0.3, α = 0.08726646259971647),
sphere = (ϵ = 0.9, α = 0.017453292519943295, sphere_par = 0.1))
```
"""
function ransacparameters(p::T=DEFAULT_PARAMETERS; kwargs...) where {T<:NamedTuple}
    newp = p
    for a in keys(kwargs)
        oldpar = get(p, a, kwargs[a])
        newpar = merge(oldpar, kwargs[a])
        newp = merge(newp, (a=>newpar,))
    end
    return newp
end

"""
    ransacparameters(p::Array{T}; kwargs...) where {T<:UnionAll}

Construct a `NamedTuple` for a given types of shapes 
using [`defaultparameters`](@ref) and override it with the kwargs.
Check the docs and examples for more.

# Examples
```julia-repl
julia> p1 = ransacparameters([FittedSphere, FittedCylinder])
(iteration = (drawN = 3, minsubsetN = 15, prob_det = 0.9,
shape_types = UnionAll[FittedSphere, FittedCylinder], τ = 900, itermax = 1000,
extract_s = :nofminset, terminate_s = :nofminset),
common = (collin_threshold = 0.2, parallelthrdeg = 1.0),
sphere = (ϵ = 0.3, α = 0.08726646259971647, sphere_par = 0.02),
cylinder = (ϵ = 0.3, α = 0.08726646259971647))

julia> p2 = ransacparameters([FittedSphere, FittedCylinder], sphere=(ϵ=0.01,), cylinder=(α=0.02,))
(iteration = (drawN = 3, minsubsetN = 15, prob_det = 0.9,
shape_types = UnionAll[FittedSphere, FittedCylinder], τ = 900, itermax = 1000,
extract_s = :nofminset, terminate_s = :nofminset),
common = (collin_threshold = 0.2, parallelthrdeg = 1.0),
sphere = (ϵ = 0.01, α = 0.08726646259971647, sphere_par = 0.02),
cylinder = (ϵ = 0.3, α = 0.02))
```
"""
function ransacparameters(p::Array{T}; kwargs...) where {T<:UnionAll}
    dpars = defaultparameters(p)
    return ransacparameters(dpars; kwargs...)
end

"""
    push2candidatesandlevels!(candidates, candidate, levels, current_level)

Push `candidate` to `candidates` and `current_level` to `levels`.
Handle the case, when `candidate` is an array of candidates.
Used in [`forcefitshapes!`](@ref).
"""
function push2candidatesandlevels!(candidates, candidate, levels, current_level)
    push!(candidates, candidate)
    push!(levels, current_level)
end

function push2candidatesandlevels!(candidates, candidate::AbstractArray, levels, current_level)
    append!(candidates, candidate)
    append!(levels, fill(current_level, size(candidate)))
end

"""
    setfloattype(nt, T)

Convert every `Real` but not `Integer` value to type `T` in a `NamedTuple` recursively.
"""
function setfloattype(nt, T)
    Keys = []
    Values = []
    for (k,v) in pairs(nt)
        push!(Keys, k)
        v = nt[k]
        if v isa NamedTuple
            push!(Values, setfloattype(v,T))
        elseif (v isa Number) && ! (v isa Integer)
            # v is a Real or a float
            push!(Values, convert(T, v))
        else
            push!(Values, v)
        end
    end
    return (; zip(Keys,Values)...)
end