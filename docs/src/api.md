# Public API

Note, that not all exported functions are considered as part of the public API.
The private API is not mature yet, expect it to change.

## Representing a point cloud

```@docs
PointCloud
PointCloud(vertices, normals, subsets)
PointCloud(vertices, normals, subsets, isenabled)
PointCloud(vertices, normals, subsets, isenabled, size, levelweight, levelscore)
```

## Parameters

For parameters nested named tuples are used, because it's easy to construct them, change their values or extend them.
Earlier Parameters.jl was used, but I could not solve the extension part, then came the named tuples.
You can construct it by hand, use the exported `ransacparameters()` function or load from a YAML file.

### Structure of the parameters

The easiest way to construct the desired named tuple is to use the `ransacparameters()` function.

```@docs
ransacparameters
```

As the docstring shows, you can construct a new one based on an old one, and give keyword arguments that will overwrite the old values.
Note, that the key-values that are in the keyword arguments will be overwritten, not the named tuple itself (so the values not listed in the keyword argument will not change).

As you can see in the above examples, the parameter must have two fields:`iteration`, `common` and the primitive types that you want to fit (`sphere`, `plane`, etc.).
Note, that `p.iteration.shape_types` field controls which primitives are fitted.
Another important thing regarding the `ransacparameter()` function, that the keyword named tuple must have a trailing comma, so this is good:

```julia
p2 = ransacparameters([FittedSphere, FittedCylinder], sphere=(ϵ=0.01,))
```

, but the following is NOT:

```julia
p2 = ransacparameters([FittedSphere, FittedCylinder], sphere=(ϵ=0.01))
```

The `ransacparameters()` function uses the not exported `default...parameters()` function, whose docstrings describes which parameters control what:

```@docs
RANSAC.defaultcommonparameters
RANSAC.defaultiterationparameters
```

Check [`RANSAC.defaultshapeparameters`](@ref) as well.

The `defaultparameters` function joins these together:

```@docs
RANSAC.defaultparameters
```

### Construct by hand

As parameters are plain named tuples, one can easily construct their own.
The above functions make heavy use of the `merge` function.
Check the code, if you wish.

### Parse from YAML file

With the help of [YAML.jl](https://github.com/BioJulia/YAML.jl) one can easily read the parameters from a YAML file.
As shown below, you can specify which parameters you want to change (the others are going to be the default ones).

An example file (`config.yml`):

```yml
plane:
  - ϵ: 0.1
  - α: 0.01

sphere:
  - ϵ: 0.2
  - α: 0.05
  # parameter in sphere fitting
  - sphere_par: 0.01

cone:
  - ϵ: 1.
  - α: 3.14
  # filter those cones, whose opening angle is less than `minconeopang` radians
  - minconeopang: 1.

iteration:
  # number of points to be sampled (length of a minimal subset)
  - drawN: 9
  # number of minimal sets sampled in one iteration
  - minsubsetN: 2
  # probability of detection
  - prob_det: 0.999
  # minimal shape size
  - τ: 10000
  # maximum number of iteration
  - itermax: 100000
  # shapes that are fitted to the point cloud
  - shape_types:
    - plane
    - sphere
    - cone

common:
  # threshold of two vectors being parallel (in degrees)
  - parallelthrdeg: 0.5
  # threshold of points being collinear
  - collin_threshold: 0.3
```

Then you can use the `readconfig` function to read the file:

```julia
julia> readconfig("config.yml")
(iteration = (drawN = 9, minsubsetN = 2, prob_det = 0.999,
shape_types = UnionAll[FittedPlane, FittedSphere, FittedCone], τ = 10000,
itermax = 100000, extract_s = :nofminset, terminate_s = :nofminset),
common = (collin_threshold = 0.3, parallelthrdeg = 0.5),
plane = (ϵ = 0.1, α = 0.01), cone = (ϵ = 1.0, α = 3.14, minconeopang = 1.0),
cylinder = (ϵ = 0.3, α = 0.08726646259971647),
sphere = (ϵ = 0.2, α = 0.05, sphere_par = 0.01))
```

```@docs
readconfig
RANSAC.DEFAULT_PARAMETERS
RANSAC.DEFAULT_SHAPE_DICT
```

## Primitives

```@docs
FittedShape
FittedPlane
FittedSphere
FittedCylinder
FittedCone
```

## Iteration

The `ransac()` function does the iteration.

```@docs
ransac
```

## Exporting the results

With the help of [JSON.jl](https://github.com/JuliaIO/JSON.jl), the resulted shapes can be easily saved to JSON files.
For this purpose, the `exportJSON()` function can be used.
Note that `io` must be specified, the "default" fallback to `stdout` is not implemented.

```@docs
exportJSON
```

A few examples:

```julia
julia> using RANSAC, StaticArrays

julia> s1 = FittedPlane(SVector(0,0,1.), SVector(12.5, 7, 24))
FittedPlane{SArray{Tuple{3},Float64,1,3}}
normal: [12.5, 7.0, 24.0], point: [0.0, 0.0, 1.0]

julia> exportJSON(stdout, s1)
{"point":[0.0,0.0,1.0],"normal":[12.5,7.0,24.0],"type":"plane"}

julia> exportJSON(stdout, s1, 2)
{
  "point": [
    0.0,
    0.0,
    1.0
  ],
  "normal": [
    12.5,
    7.0,
    24.0
  ],
  "type": "plane"
}
```

It is advised to export shapes in an array for easier processing (though I'm not a JSON expert):

```julia
julia> exportJSON(stdout, [s1])
{"primitives":[{"point":[0.0,0.0,1.0],"normal":[12.5,7.0,24.0],"type":"plane"}]}

julia> exportJSON(stdout, [s1], 2)
{
  "primitives": [
    {
      "point": [
        0.0,
        0.0,
        1.0
      ],
      "normal": [
        12.5,
        7.0,
        24.0
      ],
      "type": "plane"
    }
  ]
}
```

Works of course for different primitives:

```julia
julia> s2 = FittedSphere(SVector(1.2, 3., 5.), 1.5, true)
FittedSphere{SArray{Tuple{3},Float64,1,3}, Float64}
center: [1.2, 3.0, 5.0], R: 1.5, outwards

julia> exportJSON(stdout, [s1, s2], 2)
{
  "primitives": [
    {
      "point": [
        0.0,
        0.0,
        1.0
      ],
      "normal": [
        12.5,
        7.0,
        24.0
      ],
      "type": "plane"
    },
    {
      "outwards": true,
      "radius": 1.5,
      "center": [
        1.2,
        3.0,
        5.0
      ],
      "type": "sphere"
    }
  ]
}
```

As can be seen above, in these cases an array of `"primitives"` is printed.
Under the hood, the `toDict()` function does the job of converting the primitives to `Dict`s.

```@docs
RANSAC.toDict(s::FittedShape)
RANSAC.toDict(::Vector{T}) where {T<:Union{FittedShape,ExtractedShape}}
```
