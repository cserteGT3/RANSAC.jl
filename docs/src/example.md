# Example

This page guides you through the use of RANSAC.jl.
A public dataset is used, that was published by Schnabel et al. and is accessible [here](https://cg.cs.uni-bonn.de/en/publications/paper-details/schnabel-2009-completion/).
> R. Schnabel, P. Degener, R. Klein
> "Completion and Reconstruction with Primitive Shapes",
> in Computer Graphics Forum (Proc. of Eurographics), Vol. 28, No. 2, pages 503-512


Read the algorithm description here to get a full understanding of the algorithm and its parameters.


## Loading the data

As MeshIO and other softwares had troubles opening the `_input.obj` files, I used [MeshLab](http://www.meshlab.net/) to open them and export to non-binary encoded PLY.

```julia
julia> using FileIO

julia> m = load("fandisk_input.obj")
HomogenousMesh(
    faces: 17513xGeometryTypes.Face{3,GeometryTypes.OffsetInteger{-1,UInt32}},     vertices: 8935xGeometryTypes.Point{3,Float32},     normals: 8935xGeometryTypes.Normal{3,Float32}, )
```

## Constructing a `PointCloud`

Currently the package only handles vectors of `Float64`.

```julia
julia> using RANSAC

julia> pc = PointCloud(m.vertices, m.normals, 8)
PointCloud{Array{StaticArrays.SArray{Tuple{3},Float64,1,3},1},Array{Array{Int64,1},1},Array{Float64,1}}
PointCloud of size 8935 & 8 subsets
```

## Set the parameters

Check the [Short description](@ref) and [Parameters](@ref) pages to understand the parameters.

```julia
julia> p = RANSACParameters{Float64}();

julia> p = setalphas(p, deg2rad(10));

julia> p = setepsilons(p, 0.05);

julia> p = RANSACParameters(p, τ=50, itermax=100_000);
```

## Run!

The `ransac()` function runs the iteration.

```julia
julia> _, extr, _ = ransac(pc, p, true, reset_rand=true);
```

## See the results

There's the [RANSACVisualizer](https://github.com/cserteGT3/RANSACVisualizer.jl) package that provides a few utility functions to check the results.
