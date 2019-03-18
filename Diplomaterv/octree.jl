module Octree

import RegionTrees: AbstractRefinery, needs_refinement, refine_data

using RegionTrees: child_boundary, Cell
using RegionTrees: HyperRectangle, vertices

export OctreeNode
export iswithinrectangle, OctreeRefinery

mutable struct OctreeNode{A<:AbstractArray, B<:AbstractArray}
    cloudarray::A
    incellpoints::B
    weight::Float64
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
    # de hogyan szerzem meg a pontokat az index alapján?-> cloudarray
    points = cell.data.cloudarray[cell.data.incellpoints]
    bolarr = map(x -> iswithinrectangle(boundary, x), points)
    # new cell inherits the weight and ++ of the depth
    d = cell.data.depth + 1
    OctreeNode(cell.data.cloudarray, cell.data.incellpoints[bolarr], cell.data.weight, d)
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
