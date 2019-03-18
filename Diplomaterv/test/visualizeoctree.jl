## Octree test

# construct a pointcloud

using StaticArrays: SVector

minV = [0.0, 0, 0];
maxV = [5.5, 10, 7.8];
mmV = [minV maxV];

# corner points
corners=[SVector(mmV[1,i], mmV[2,j], mmV[3,k]) for i in 1:2 for j in 1:2 for k in 1:2]

full8p=[SVector(0.75*maxV[1]+(i/12)*maxV[1], 0.75*maxV[2]+(j/12)*maxV[2], 0.75*maxV[3]+(k/12)*maxV[3]) for i in 0:2 for j in 0:2 for k in 0:2 ];
quadrant=[SVector((i/3)*maxV[1], (j/3)*maxV[2], (k/3)*maxV[3]) for i in 0:1 for j in 0:1 for k in 0:1 ];
online = [SVector((3/4)*maxV[1], (j/2)*maxV[2], (k/2)*maxV[3]) for j in 0:0.5:1 for k in 0:0.5:1 ];

pc = vcat(full8p, quadrant, online)


## visualize

using Makie
s = Scene()
scatter!(corners, color=:blue)
scatter!(full8p, color=:red)
scatter!(quadrant, color= :green)
scatter!(online, color=:brown)

# build octree
using RegionTrees: adaptivesampling!, Cell, allleaves

include("../octree.jl")
using .Octree

root = Cell(SVector{3}(minV), SVector{3}(maxV), OctreeNode(pc, collect(1:length(pc)), 0.0, 1))

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
    global i +=1
end
s
using GeometryTypes
kubi = GeometryTypes.HyperRectangle(Vec3(root.boundary.origin...), Vec3(root.boundary.widths...))

s22  = mesh(kubi, color = (:green, 0.2), transparency = true)
text!(s22, "1", position = Vec3(root.boundary.origin...), textsize = 0.5)
