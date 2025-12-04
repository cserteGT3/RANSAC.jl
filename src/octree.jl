## extend RegionTrees

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
    while true
        c = parent(c)
        c === nothing && return nothing
        c.data.depth == n && return c
    end
    return nothing
end

"""
    struct RANSACCloud{A<:AbstractArray, B<:AbstractArray, C<:AbstractArray}

A struct to wrap a point cloud. Stores the vertices, the normals,
    the octree of the vertices,
    the subsets (as an array of arrays of vertice indexes),
    an array to indicate if a point is part of an already extracted primitive,
    the size of the point cloud (number of vertices),
    the weight of each octree level (populated during the construction of the octree,
    by copying the given element),
    an array to store the sum of the score of already
    extracted primitives for each octree level.
"""
struct RANSACCloud{A<:AbstractArray, B<:AbstractArray, C<:AbstractArray, T}
    vertices::A
    normals::B
    octree::Cell{T}
    subsets::Vector{Vector{Int}}
    isenabled::BitArray{1}
    size::Int
    levelweight::C
    levelscore::C
    #= 
    # Arguments
    - `vertices::AbstractArray`: an array of vertices.
    - `normals::AbstractArray`: an array of surface normals.
    - `subsets::Vector{Vector{Int}}`: an array consisting of index-arrays,
        that describe which point is in which subset.
    - `isenabled::BitArray{1}`: an array indicating if the given point is enabled
        (`true`) or already has been extracted (`false`).
    - `size::Int`: number of points in the cloud.
    - `levelweight::Vector{Float}`: the weight of every octree level.
    - `levelscore::Vector{Float}`: for every octree level,
        store the sum of the scores of primitives that have been extracted.
    =#
end

Base.show(io::IO, x::RANSACCloud{A,B,C}) where {A,B,C} =
    print(io, "RANSACCloud of size $(x.size)")

Base.show(io::IO, ::MIME"text/plain", x::RANSACCloud{A,B,C}) where {A,B,C} =
    print(io, "RANSACCloud of size $(x.size) & $(length(x.subsets)) subsets")

"""
    nomodRANSACCloud(vertices, normals, subsets)

Construct a `RANSACCloud` without touching the vertices, normals and subsets.
Other fields are computed.

# Arguments
- `vertices`: an array of vertices.
- `normals`: an array of surface normals.
- `subsets::Vector{Vector{Int}}`: a list of indexes for each subset.
"""
function RANSACCloud(vertices, normals, subsets)
    octree = buildoctree(vertices)
    octree_d = octreedepth(octree)
    s = size(vertices, 1)
    levelscore = zeros(Float64, (octree_d,))
    levelweight = fill!(similar(levelscore), 1/octree_d)
    RANSACCloud(vertices, normals, octree, subsets, trues(s), s, levelscore, levelweight)
end

"""
    RANSACCloud(vertices, normals, numofsubsets::Int)

Construct a `RANSACCloud` with `numofsubsets` number of random subsets.
Vertices and normals are converted to array of `SVector`s.

# Arguments
- `vertices`: an array of vertices.
- `normals`: an array of surface normals.
- `numofsubsets::Int`: number of subsets.
"""
function RANSACCloud(vertices, normals, numofsubsets::Int)
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
    return RANSACCloud(vertices, normals, subs)
end

"""
    struct OctreeNode{A}

Store the data attached to one cell of the octree.
Stores the indices of the points that are inside of the cell,
and its own depth (depth of the root is 1).
"""
struct OctreeNode{A}
    geometries::Vector{A}
    nodeindices::Vector{Int}
    depth::Int
end

Base.show(io::IO, x::OctreeNode{A}) where {A} =
    print(io, "OctreeNode: $(length(x.incellpoints)) ps, $(x.depth) d")

Base.show(io::IO, ::MIME"text/plain", x::OctreeNode{A}) where {A} =
    print(io, """OctreeNode{$A}\n with $(length(x.incellpoints)) points, at depth $(x.depth)""")

struct OctreeRefinery <: AbstractRefinery
    count::Int
end

function needs_refinement(r::OctreeRefinery, cell)
    length(cell.data.nodeindices) > r.count
end

hr2b(hr::RegionTrees.HyperRectangle) = boundingbox([sv2p(hr.origin), sv2p(hr.origin + hr.widths)])

hr2b(c::RegionTrees.Cell) = hr2b(c.boundary)

function refine_data(r::OctreeRefinery, cell::Cell, indices)
    # new boundary
    boundary = child_boundary(cell, indices)
    newbox = hr2b(boundary)
    # points that are in the current cell
    geoms = @view cell.data.geometries[cell.data.nodeindices]
    # points that are in the new cell
    newinds = Int[]
    for (g, i) in zip(geoms, cell.data.nodeindices)
        intersects(g, newbox) && push!(newinds, i)
    end
    # new cell inherits ++ of the depth
    d = cell.data.depth + 1
    OctreeNode(cell.data.geometries, newinds, d)
end

function updatelevelweight(pc, x = 9//10)
    P = pc.levelweight
    σ = pc.levelscore
    w = sum( σ[i]/P[i] for i in eachindex(P) )
    for i in eachindex(P)
        pc.levelweight[i] = x*σ[i]/(w*P[i]) + (1-x)/length(P)
    end
end

"""
    octreedepth(pc::RANSACCloud)

Compute the depth of the octree stored in the cloud.
"""
function octreedepth(pc::RANSACCloud)
    return octreedepth(pc.octree)
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

"""
    buildoctree(vertices)

Build octree of the vertices.
"""
function buildoctree(vertices)
    b = boundingbox(vertices)
    allinds = collect(1:size(vertices, 1))
    r = OctreeRefinery(8)
    width = b.max - b.min
    root = Cell(p2sv(b.min), v2sv(width), OctreeNode(vertices, allinds, 1))
    adaptivesampling!(root, r)
    return root
end
