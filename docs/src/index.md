# RANSAC.jl

This package implements the efficient RANSAC algorithm for point clouds.
Paper can be found here: [`Schnabel2007`](https://cg.cs.uni-bonn.de/en/publications/paper-details/schnabel-2007-efficient/).

>  R. Schnabel, R. Wahl, R. Klein
>	"Efficient RANSAC for Point-Cloud Shape Detection",
>	in Computer Graphics Forum, Vol. 26, No. 2, pages 214-226,
>	Blackwell Publishing, June 2007

A full page ([Efficient RANSAC](@ref)) is dedicated to describe the algorithm and to help to understand the parameters.
If something is not clear, please open an issue or check the original paper.

## Quick tour

The efficient RANSAC algorithm is used to segment and fit primitive shapes (sphere, plane, cylinder, torus, cone) to point clouds.
Up to now mostly C++ and Python implementations have been published, this is the first one in Julia (as far as I know).

### Main features

* easy-to-use primitive recognition
* extensible: it's easy to add new primitive shapes
* fast (to be honest, it's not yet comparable with the reference and PCL implementations)

### Differences from the reference implementation

* no bitmap
* separate parameters for each shape
* no tori

### Install the package

```julia
] add https://github.com/cserteGT3/RANSAC.jl
```

### Load a point cloud

Get an example point cloud. Use your own, or download a public dataset (for example the one used in `Schnabel2009`, read more about on the [Example](@ref) page)
You can use [MeshIO.jl](https://github.com/JuliaIO/MeshIO.jl) to load models.

```julia
using FileIO
m = load("testm.obj")
```

### Construct a `RANSACCloud`

```julia
using RANSAC
pc = RANSACCloud(m.position, m.normals, 2)
rparams = ransacparameters()
```

### Run the iteration

```julia
extr, _ = ransac(pc, rparams, true);
```

See the [Example](@ref) page for a detailed tour.

## Short algorithm description

The input of the algorithm is a point cloud of size ``N`` with points and associated surface normals.
The output is a set of primitive shapes with corresponding sets of points, and the rest of the points that do not belong to any primitives.
Primitive shapes can be: plane, sphere, cylinder, cone and torus (though torus is not implemented yet).

In every iteration, new shape candidates are created by fitting primitives to randomly sampled minimal sets.
Every shape primitive is generated for every minimal set and the valid ones are continuously collected in set ``C``.
Then every candidate is scored and the one with the highest score is considered the best.
A candidate is only extracted if the probability (``p_t``) that no better candidates are in ``C`` is high enough.
If the best candidate is chosen to be extracted, its points are removed from the point cloud, and every other candidate that has a removed point is also deleted from ``C``.
The iteration continues until the probability that every at least ``\tau`` sized shape is found is larger than a parameter threshold.
The score of a candidate is defined by the number of compatible points.
A point is compatible if it is in the ``\epsilon`` band of the shape, and its normal does not deviate from the surface normal more than an ``\alpha`` angle.
Also, only those points are considered that count towards the largest connecting component in the parameter space bitmap of the shape.