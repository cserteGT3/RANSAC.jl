## extend RegionTrees

"""
    getcellandparents(cell::Cell)

Collect the cell and it's parents into an array.
"""
function getcellandparents(cell::Cell)
    p = allparents(cell)
    aoc = [i for i in p]
    pushfirst!(aoc, cell)
    aoc
end

"""
    getnthcell(c::Cell, n)

Get the cell at the n-th level (at the depth `n`) from the parents of a cell.
If `n` is larger than the depth of `c`, return nothing
(it means, that the asked depth is below the current level).
Also return nothing, if the asked level is below the minimum level, currently set to 1.
"""
function getnthcell(c::Cell, n)
    #n < 1 && throw(BoundsError(c, n))
    n < 1 && return nothing
    maxd = c.data.depth
    n == maxd && return c
    for cell in allparents(c)
        cell.data.depth == n && return cell
    end
    return nothing
end

"""
    mutable struct PointCloud{A<:AbstractArray, B<:AbstractArray, C<:AbstractArray}

A struct to wrap a point cloud. Stores the vertices, the normals,
    the subsets (as an array of arrays of vertice indexes),
    an array to indicate if a point is part of an already extracted primitive,
    the size of the point cloud (number of vertices),
    the weight of each octree level (populated during the construction of the octree,
    by copying the given element),
    an array to store the sum of the score of already
    extracted primitives for each octree level.
"""
mutable struct PointCloud{A<:AbstractArray, B<:AbstractArray, C<:AbstractArray}
    vertices::A
    normals::A
    subsets::B
    isenabled::BitArray{1}
    size::Int
    levelweight::C
    levelscore::C
end

Base.show(io::IO, x::PointCloud{A,B,C}) where {A,B,C} =
    print(io, "PointCloud of size $(x.size) & $(length(x.subsets)) subsets")

Base.show(io::IO, ::MIME"text/plain", x::PointCloud{A,B,C}) where {A,B,C} =
    print(io, "PointCloud{$A,$B,$C}\n", x)

"""
    PointCloud(vertices, normals, subsets, isenabled, size, levelweight, levelscore)

Constructor that converts vertices and normals to array of `SVector{3,Float64}`.

# Arguments
- `vertices::AbstractArray`: an array of vertices.
- `normals::AbstractArray`: an array of surface normals.
- `subsets::AbstractArray`: an array consisting of index-arrays,
    that describe which point is in which subset.
- `isenabled::BitArray{1}`: an array indicating if the given point is enabled
    (`true`) or already has been extracted (`false`).
- `size::Int`: number of points in the cloud.
- `levelweight::Vector{Float}`: the weight of every octree level.
- `levelscore::Vector{Float}`: for every octree level,
    store the sum of the scores of primitives that have been extracted.
"""
function PointCloud(vertices, normals, subsets, isenabled, size, levelweight, levelscore)
    vc = [SVector{3,Float64}(v) for v in vertices]
    nc = [SVector{3,Float64}(v) for v in normals]
    PointCloud(vc, nc, subsets, isenabled, size, levelweight, levelscore)
end

"""
    PointCloud(vertices, normals, subsets, isenabled)

Construct a `PointCloud` with the given fields.
Not defined fields default to: `size=size(vertices,1)`,
`levelweight` and `levelscore` defaults to `[0]` with appropriate type.
"""
function PointCloud(vertices, normals, subsets, isenabled)
    za = zeros(eltype(eltype(vertices)), 1)
    return PointCloud(vertices, normals, subsets, isenabled, size(vertices,1),za,za)
end

"""
    PointCloud(vertices, normals, subsets)

Construct a `PointCloud` with the given fields.
Not defined fields default to: `size=size(vertices,1)`,
`levelweight` and `levelscore` defaults to `[0]` with appropriate type,
and all vertices are enabled.
"""
function PointCloud(vertices, normals, subsets)
    return PointCloud(vertices, normals, subsets, trues(size(vertices,1)))
end

"""
    PointCloud(vertices, normals, numofsubsets::Int)

Construct a `PointCloud` with `numofsubsets` number of random subsets.
Every point is enabled, and other fields default to: `size=size(vertices,1)`,
`levelweight` and `levelscore` defaults to `[0]` with appropriate type,
and all vertices are enabled.
"""
function PointCloud(vertices, normals, numofsubsets::Int)
    @assert numofsubsets > 0 "At least 1 subset please!"
    # subset length
    function makesubset(l, n)
        ssl = fld(l, n)
        alls = randperm(l)
        subsets = [ alls[(i-1)*ssl+1:i*ssl] for i in 1:n-1]
        push!(subsets, alls[(n-1)*ssl+1:end])
        subsets
    end
    subs = makesubset(length(vertices), numofsubsets)
    return PointCloud(vertices, normals, subs)
end

"""
    struct OctreeNode{B<:AbstractArray}

Store the data attached to one cell of the octree.
Stores the reference of the point cloud,
the indices of the points that are inside of the cell,
and its own depth (depth of the root is 1).
"""
struct OctreeNode{B<:AbstractArray, I<:Integer}
    pointcloud::PointCloud
    incellpoints::B
    depth::I
end

Base.show(io::IO, x::OctreeNode{A}) where {A} =
    print(io, "OctreeNode: $(length(x.incellpoints)) ps, $(x.depth) d")

Base.show(io::IO, ::MIME"text/plain", x::OctreeNode{A}) where {A} =
    print(io, """OctreeNode{$A}\n with a $(x.pointcloud.size) pointcloud, at depth $(x.depth)""")

struct OctreeRefinery <: AbstractRefinery
    count::Int64
end

function needs_refinement(r::OctreeRefinery, cell)
    length(cell.data.incellpoints) > r.count
end

function refine_data(r::OctreeRefinery, cell::Cell, indices)
    boundary = child_boundary(cell, indices)
    # a visszatérési érték azoknak a pontoknak az indexe, amelyik a boundary-ban van
    # de hogyan szerzem meg a pontokat az index alapján?-> pointcloud
    points = @view cell.data.pointcloud.vertices[cell.data.incellpoints]
    bolarr = map(x -> iswithinrectangle(boundary, x), points)
    # new cell inherits the weight and ++ of the depth
    d = cell.data.depth + 1
    lw = cell.data.pointcloud.levelweight
    if d > length(lw)
        push!(lw, convert(eltype(lw), d))
        push!(cell.data.pointcloud.levelscore, convert(eltype(lw), d))
    end
    OctreeNode(cell.data.pointcloud, cell.data.incellpoints[bolarr], d)
end

"""
    iswithinrectangle(rect::HyperRectangle, p)

Decide if p (3 dimensional point) is within rect.

Works only with those prisms which axes are aligned with the coordinate frame.
"Bottom/left" gives `false` and "top/right" gives `true`.
"""
function iswithinrectangle(rect::HyperRectangle, p)
    vs = vertices(rect)
    vmin = vs[1,1,1]
    vmax = vs[2,2,2]
    for i in 1:3
        vmin[i] < p[i] || return false
        vmax[i] >= p[i] || return false
    end
    return true
end

function updatelevelweight(pc, x = 0.9)
    P = pc.levelweight
    σ = pc.levelscore
    w = sum( σ[i]/P[i] for i in eachindex(P) )
    for i in eachindex(P)
        pc.levelweight[i] = x*σ[i]/(w*P[i]) + (1-x)/length(P)
    end
end

"""
    octreedepth(pc::PointCloud, n::Int=8)

Construct an octree with `OctreeRefinery(n)` and compute its depth.
(It's only an octree if `n==8`, but nevermind.)
"""
function octreedepth(pc::PointCloud, n::Int=8)
    minV, maxV = findAABB(pc.vertices)
    octree=Cell(SVector{3}(minV), SVector{3}(maxV), OctreeNode(pc, collect(1:pc.size), 1))
    r = OctreeRefinery(n)
    adaptivesampling!(octree, r)
    return octreedepth(octree)
end

"""
    octreedepth(cell::Cell)

Iterate through the leafes and search the maximum depth (`cell.data.depth`).
"""
function octreedepth(cell::Cell)
    maxdepth = cell.data.depth
    alll = allleaves(cell)
    for a in alll
        if a.data.depth > maxdepth
            maxdepth = a.data.depth
        end
    end
    return maxdepth
end
