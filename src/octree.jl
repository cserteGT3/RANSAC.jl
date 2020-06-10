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
struct RANSACCloud{A<:AbstractArray, B<:AbstractArray, C<:AbstractArray}
    vertices::A
    normals::B
    octree::Cell
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
function nomodRANSACCloud(vertices, normals, subsets)
    octree = buildoctree(vertices)
    octree_d = octreedepth(octree)
    s = size(vertices, 1)
    levelscore = zeros(eltype(eltype(vertices)), (octree_d,))
    levelweight = fill!(similar(levelscore), 1/octree_d)
    RANSACCloud(vertices, normals, octree, subsets, trues(s), s, levelscore, levelweight)
end

"""
    RANSACCloud(vertices, normals, subsets; force_eltype::Union{Nothing,DataType}=nothing)

Construct a `RANSACCloud` with the given `subsets`.
Vertices and normals are converted to array of `SVector`s.

# Arguments
- `vertices`: an array of vertices.
- `normals`: an array of surface normals.
- `subsets::Vector{Vector{Int}}`: a list of indexes for each subset.
- `force_eltype::Union{Nothing,DataType}=nothing`:
    an element type can be forced for the normals and vertices
    (in practice `Float64` or `Float32`). If not specified,
    the element type of the passed `vertices` and `normals` will be used.
"""
function RANSACCloud(vertices, normals, subsets; force_eltype::Union{Nothing,DataType}=nothing)
    @assert size(vertices) == size(normals) "Every point must have a normal."
    v_type = force_eltype === nothing ? eltype(eltype(vertices)) : force_eltype
    n_type = force_eltype === nothing ? eltype(eltype(normals)) : force_eltype
    vc = [SVector{3,v_type}(v) for v in vertices]
    nc = [SVector{3,n_type}(v) for v in normals]
    nomodRANSACCloud(vc, nc, subsets)
end

"""
    RANSACCloud(vertices, normals, numofsubsets::Int)

Construct a `RANSACCloud` with `numofsubsets` number of random subsets.
Vertices and normals are converted to array of `SVector`s.

# Arguments
- `vertices`: an array of vertices.
- `normals`: an array of surface normals.
- `numofsubsets::Int`: number of subsets.
- `force_eltype::Union{Nothing,DataType}=nothing`:
    an element type can be forced for the normals and vertices
    (in practice `Float64` or `Float32`). If not specified,
    the element type of the passed `vertices` and `normals` will be used.
"""
function RANSACCloud(vertices, normals, numofsubsets::Int; force_eltype::Union{Nothing,DataType}=nothing)
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
    return RANSACCloud(vertices, normals, subs; force_eltype=force_eltype)
end

"""
    struct OctreeNode{B<:AbstractArray}

Store the data attached to one cell of the octree.
Stores the indices of the points that are inside of the cell,
and its own depth (depth of the root is 1).
"""
struct OctreeNode{B<:AbstractArray, I<:Integer}
    incellpoints::B
    depth::I
end

Base.show(io::IO, x::OctreeNode{A}) where {A} =
    print(io, "OctreeNode: $(length(x.incellpoints)) ps, $(x.depth) d")

Base.show(io::IO, ::MIME"text/plain", x::OctreeNode{A}) where {A} =
    print(io, """OctreeNode{$A}\n with $(length(x.incellpoints)) points, at depth $(x.depth)""")

struct OctreeRefinery{A<:AbstractArray} <: AbstractRefinery
    count::Int
    vertices::A
end

function needs_refinement(r::OctreeRefinery, cell)
    length(cell.data.incellpoints) > r.count
end

function refine_data(r::OctreeRefinery, cell::Cell, indices)
    # new boundary
    boundary = child_boundary(cell, indices)
    # points that are in the current cell
    points = @view r.vertices[cell.data.incellpoints]
    # points that are in the new cell
    bolarr = map(x -> iswithinrectangle(boundary, x), points)
    # new cell inherits ++ of the depth
    d = cell.data.depth + 1
    OctreeNode(cell.data.incellpoints[bolarr], d)
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
    l = size(vertices, 1)
    minV, maxV = findAABB(vertices)
    root=Cell(SVector{3}(minV), SVector{3}(maxV), OctreeNode(collect(1:l), 1))
    r = OctreeRefinery(8, vertices)
    adaptivesampling!(root, r)
    return root
end
