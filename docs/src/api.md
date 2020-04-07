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

With the help of [Parameters.jl](https://github.com/mauro3/Parameters.jl) it's easy to parameterize the algorithm.
The `RANSACParameters` type collects all the parameters, though its fields are subject to change, the current fields and default values are listed below.
Note that fields that are float types must use float literals (so `1.0` is ok, but `1` is not).

### Constructing from code

```julia
@with_kw struct RANSACParameters{R<:Real} @deftype R
    ϵ_plane = 0.3
    α_plane = deg2rad(5)

    ϵ_sphere = 0.3
    α_sphere = deg2rad(5)

    ϵ_cylinder = 0.3
    α_cylinder = deg2rad(5)

    ϵ_cone = 0.3
    α_cone = deg2rad(5)
    # filter those cones, whose opening angle is less than `minconeopang` radians
    minconeopang = deg2rad(2)

    ϵ_torus = 0.3
    α_torus = deg2rad(5)

    # number of points to be sampled (length of a minimal subset)
    drawN::Int = 3; @assert drawN>2
    # number of minimal sets sampled in one iteration
    minsubsetN::Int = 15; @assert minsubsetN>0
    # probability of detection
    prob_det = 0.9
    # minimal shape size
    τ::Int = 900
    # maximum number of iteration
    itermax::Int = 1000

    # threshold of two vectors being parallel (in degrees)
    parallelthrdeg = 1.
    # threshold of points being collinear
    collin_threshold = 0.2
    # parameter in sphere fitting
    sphere_par = 0.02

    # shapes that are fitted to the point cloud
    shape_types::Array{Symbol,1} = [:sphere, :plane, :cylinder, :cone]
end
```

You can use the following functions to set every ``\alpha`` or ``\epsilon`` parameters to a certain value:

```@docs
setalphas
setepsilons
```

### Parse from YAML file

With the help of [YAML.jl](https://github.com/BioJulia/YAML.jl) one can easily read the parameters from a YAML file.
As shown below, you can specify which parameters you want to change (the others are going to be the default ones).
Note, that in the YAML file, you must specify fields with float types as float-literals.

An example file (`config.yml`):
```yml
# ϵ_plane is float, so it can't be `1`
ϵ_plane: 1.0
α_plane: 0.0872

# maximum number of iteration
itermax: 100

# threshold of two vectors being parallel (in degrees)
parallelthrdeg: 1.2
# threshold of points being collinear
collin_threshold: 0.22
# parameter in sphere fitting
sphere_par: 0.025

# shapes that are fitted to the point cloud
# must be an array
shape_types:
  - plane
  - cone
```
Then you can use the `readconfig` function to read the file:
```julia
julia> readconfig("config.yml")
RANSACParameters{Float64}
  ϵ_plane: Float64 1.0
  α_plane: Float64 0.0872
  ϵ_sphere: Float64 0.3
  α_sphere: Float64 0.08726646259971647
  ϵ_cylinder: Float64 0.3
  α_cylinder: Float64 0.08726646259971647
  ϵ_cone: Float64 0.3
  α_cone: Float64 0.08726646259971647
  minconeopang: Float64 0.03490658503988659
  ϵ_torus: Float64 0.3
  α_torus: Float64 0.08726646259971647
  drawN: Int64 3
  minsubsetN: Int64 15
  prob_det: Float64 0.9
  τ: Int64 900
  itermax: Int64 100
  parallelthrdeg: Float64 1.2
  collin_threshold: Float64 0.22
  β: Float64 1.0
  sphere_par: Float64 0.025
  transl_conn: Symbol eight
  shape_types: Array{Symbol}((2,))
  extract_s: Symbol nofminset
  terminate_s: Symbol nofminset
```
You can also specify the type of the parameter struct (can be `Float32` as well).
```julia
julia> readconfig("config.yml", RANSACParameters{Float32})
RANSACParameters{Float32}
  ϵ_plane: Float32 1.0f0
  α_plane: Float32 0.0872f0
  ϵ_sphere: Float32 0.3f0
  α_sphere: Float32 0.08726646f0
  ϵ_cylinder: Float32 0.3f0
  α_cylinder: Float32 0.08726646f0
  ϵ_cone: Float32 0.3f0
  α_cone: Float32 0.08726646f0
  minconeopang: Float32 0.034906585f0
  ϵ_torus: Float32 0.3f0
  α_torus: Float32 0.08726646f0
  drawN: Int64 3
  minsubsetN: Int64 15
  prob_det: Float32 0.9f0
  τ: Int64 900
  itermax: Int64 100
  parallelthrdeg: Float32 1.2f0
  collin_threshold: Float32 0.22f0
  β: Float32 1.0f0
  sphere_par: Float32 0.025f0
  transl_conn: Symbol eight
  shape_types: Array{Symbol}((2,))
  extract_s: Symbol nofminset
  terminate_s: Symbol nofminset
```

```@docs
readconfig
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
For this purpose, the `exportJSON()` function can be used. Note that `io` must be specified, the "default" fallback to `stdout` is not implemented.
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
RANSAC.toDict(::Vector{T}) where {T<:Union{FittedShape,ShapeCandidate,ScoredShape}}
```
