## Octree test
using StaticArrays: Svector

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

scatter(corners, color=:blue)
scatter!(full8p, color=:red)
scatter!(quadrant, color= :green)
scatter!(online, color=:brown)
