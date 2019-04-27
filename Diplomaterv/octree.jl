module Octree

import RegionTrees: AbstractRefinery, needs_refinement, refine_data

using RegionTrees: child_boundary, Cell
using RegionTrees: HyperRectangle, vertices
using Random: randperm

export OctreeNode
export iswithinrectangle, OctreeRefinery
export PointCloud
export getcellandparents
export allparents
export updatelevelweight

## extend RegionTrees
"""
    cellandparents(cell::Cell)

Returns all parents of a cell.
"""
function allparents(cell::Cell)
    Channel() do c
        queue = [cell]
        while !isempty(queue)
            current = pop!(queue)
            p = parent(current)
            if ! (p === nothing)
                put!(c, p)
                push!(queue, p)
            end

        end
    end
end

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

end #module

module ExamplePointCloud
# construct an example pointcloud

using StaticArrays: SVector
export buildexamplecloud

const maxV = [5.5, 10, 7.8];

"""
    buildexamplecloud(maxp = [5.5, 10, 7.8])

Create an example point cloud for octree tests and visualization.

One cornerpoint of the prism is the origin and the other can be given as `maxp`.
Return tuple's fields: `corners`, `pc`, `minP`, `maxP`.
"""
function buildexamplecloud(maxp = [5.5, 10, 7.8])
    @assert length(maxp) == 3
    minp = [0.0, 0, 0];
    mmV = [minp maxp];

    # corner points
    corners=[SVector(mmV[1,i], mmV[2,j], mmV[3,k]) for i in 1:2 for j in 1:2 for k in 1:2]

    full8p=[SVector(0.75*maxp[1]+(i/12)*maxp[1], 0.75*maxp[2]+(j/12)*maxp[2], 0.75*maxp[3]+(k/12)*maxp[3]) for i in 0:2 for j in 0:2 for k in 0:2]
    quadrant=[SVector((i/3)*maxp[1], (j/3)*maxp[2], (k/3)*maxp[3]) for i in 0:1 for j in 0:1 for k in 0:1]
    online=[SVector((3/4)*maxp[1], (j/2)*maxp[2], (k/2)*maxp[3]) for j in 0:0.5:1 for k in 0:0.5:1]

    pc = vcat(full8p, quadrant, online)
    return (corners = corners, pc = pc, minP = minp, maxP = maxp)
end

end #module
