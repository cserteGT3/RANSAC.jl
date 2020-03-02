# RANSAC.jl

This package implements the efficient RANSAC algorithm for point clouds.
Paper can be found [here](https://cg.cs.uni-bonn.de/en/publications/paper-details/schnabel-2007-efficient/).

>  R. Schnabel, R. Wahl, R. Klein
>	"Efficient RANSAC for Point-Cloud Shape Detection",
>	in Computer Graphics Forum, Vol. 26, No. 2, pages 214-226,
>	Blackwell Publishing, June 2007

## Efficient RANSAC

The efficient RANSAC algorithm is used to detect primitive shapes in point clouds.
Primitive shapes can be: plane, sphere, cylinder, cone and torus, though torus is not implemented yet.

The input of the algorithm is a point cloud with associated surface normals.
The output is a set of primitive shapes with corresponding sets of points, and the rest of the points that do not belong to any primitives.

## Getting started

1. Install the package:
```julia
] add https://github.com/cserteGT3/RANSAC.jl
```
2. Get an example point cloud. Use your own, or download a public dataset (for example the one used in `Schnabel2007`, read more about [here](@ref Example))
3. You can use [FileIO] and [MeshIO] to load models.
```julia
using FileIO
m = load("testm.obj")
pc = PointCloud(m, 2)
rparams = RANSACParameters{Float64}()
```

## The Efficient RANSAC algorithm

A [full page](@ref Efficient RANSAC) is dedicated to describe the algorithm and to help to understand the parameters.
