module Octree

import RegionTrees: AbstractRefinery, needs_refinement, refine_data

using RegionTrees: child_boundary, Cell
using RegionTrees: HyperRectangle, vertices

export OctreeNode
export iswithinrectangle

mutable struct OctreeNode{A<:AbstractArray}
    indexes::A
    weight::Float64
end

struct OctreeRefinery <: AbstractRefinery
    count::Int64
end

function needs_refinement(r::OctreeRefinery, cell)
    cell.data.indexes > r.count
end

function refine_data(r::OctreeRefinery, cell::Cell, indices)
    boundary = child_boundary(cell, indices)
    # a visszatérési érték azoknak a pontoknak az indexe, amelyik a boundary-ban van
    indexes = cell.data.indexes
    # de hogyan szerzem meg a pontokat az index alapján?
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
