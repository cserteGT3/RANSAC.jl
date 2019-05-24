using Pkg
Pkg.activate()

# every include
using LinearAlgebra
using StaticArrays
using Random
using Revise
using Colors
using Makie
using RegionTrees: adaptivesampling!, Cell, allleaves
import GeometryTypes

includet("../RANSAC.jl")
using .RANSAC

## Construct an example pointcloud

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

function testvisoctree()
    ## visualize
    pcd = buildexamplecloud();
    corners = pcd.corners
    cloud = pcd.pc

    s = Scene()
    scatter!(s, corners, color=:blue)
    scatter!(s, cloud, color= :green)

    # build octree

    pc = PointCloud(cloud, cloud, 1)
    root = Cell(SVector{3}(pcd.minP), SVector{3}(pcd.maxP), OctreeNode(pc, collect(1:length(cloud)), 1))

    r = OctreeRefinery(8)

    adaptivesampling!(root, r)

    tp = 0.1
    @show "let's get started"
    i = 1
    for leaf in allleaves(root)
        ku = GeometryTypes.HyperRectangle(Vec3(leaf.boundary.origin...), Vec3(leaf.boundary.widths...))
        if leaf.data.depth == 2
            mesh!(s, ku, color = (:red, tp), transparency = true)
            text!(s, "$i,2", position = Vec3(leaf.boundary.origin...), textsize = 0.5 )
        elseif leaf.data.depth == 3
            mesh!(s, ku, color = (:green, tp), transparency = true)
            text!(s, "$i,3", position = Vec3(leaf.boundary.origin...), textsize = 0.5 )
        else
            println("please no")
        end
        @show leaf.data.depth
        i +=1
    end
    s, pc
end

scene, pc = testvisoctree()
