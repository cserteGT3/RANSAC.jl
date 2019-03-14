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

end
