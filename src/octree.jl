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
"""
function PointCloud(vertices, normals, subsets, isenabled, size, levelweight, levelscore)
    vc = [SVector{3,Float64}(v) for v in vertices]
    nc = [SVector{3,Float64}(v) for v in normals]
    PointCloud(vc, nc, subsets, isenabled, size, levelweight, levelscore)
end

"""
    PointCloud(vertices, normals, subsets, isenabled)

Construct a `PointCloud` with filling it's `size` and `levelweight` fields.
"""
function PointCloud(vertices, normals, subsets, isenabled)
    return PointCloud(vertices, normals, subsets, isenabled, length(vertices), [1.0], [1.0])
end

"""
    PointCloud(vertices, normals, subsets)

Construct a `PointCloud`, filling it's other fields.

Every point is enabled, and the weight vector defaults to `[1.0]`.
"""
function PointCloud(vertices, normals, subsets)
    return PointCloud(vertices, normals, subsets, trues(length(vertices)), length(vertices), [1.0], [1.0])
end

"""
    PointCloud(vertices, normals, numofsubsets::Int)

Construct a `PointCloud` with `numofsubsets` random subsets.

Every point is enabled, and the weight vector defaults to `[1.0]`.
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
    return PointCloud(vertices, normals, subs, trues(length(vertices)), length(vertices), [1.0], [1.0])
end

struct OctreeNode{B<:AbstractArray}
    pointcloud::PointCloud
    incellpoints::B
    depth::Int64
end

Base.show(io::IO, x::OctreeNode{A}) where {A} =
    print(io, "OctreeNode: $(length(x.incellpoints)) ps, $(x.depth) d")

Base.show(io::IO, ::MIME"text/plain", x::OctreeNode{A}) where {A} =
    print(io, """OctreeNode{$A}\n$(x.pointcloud.size) pointcloud, $(x.depth) deep, $(length(x.incellpoints)) ps""")

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
