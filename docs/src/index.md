# RANSAC.jl

This package implements the efficient RANSAC algorithm for point clouds.
Paper can be found here: [`Schnabel2007`](https://cg.cs.uni-bonn.de/en/publications/paper-details/schnabel-2007-efficient/).

>  R. Schnabel, R. Wahl, R. Klein
>	"Efficient RANSAC for Point-Cloud Shape Detection",
>	in Computer Graphics Forum, Vol. 26, No. 2, pages 214-226,
>	Blackwell Publishing, June 2007

## The Efficient RANSAC paradigm

The efficient RANSAC algorithm is used to detect primitive shapes in point clouds.
Primitive shapes can be: plane, sphere, cylinder, cone and torus, though torus is not implemented yet.

The input of the algorithm is a point cloud with associated surface normals.
The output is a set of primitive shapes with corresponding sets of points, and the rest of the points that do not belong to any primitives.

## Quick tour

1. Install the package:
```julia
] add https://github.com/cserteGT3/RANSAC.jl
```
2. Get an example point cloud. Use your own, or download a public dataset (for example the one used in `Schnabel2007`, read more about on the [Example](@ref) page)
3. You can use [FileIO](https://github.com/JuliaIO/FileIO.jl) and [MeshIO](https://github.com/JuliaIO/MeshIO.jl) to load models.
```julia
using FileIO
m = load("testm.obj")
using RANSAC
pc = PointCloud(m, 2)
rparams = RANSACParameters{Float64}()
```

## The Efficient RANSAC algorithm

A full page ([Efficient RANSAC](@ref)) is dedicated to describe the algorithm and to help to understand the parameters.
If something is not clear, please open an issue or check the original paper.

## Differences from the reference implementation

* no bitmap
* separate parameters for each shape
* no tori
