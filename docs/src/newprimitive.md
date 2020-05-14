# Introducing a new primitive shape

!!! warning "API under development"
    The primitive fitting API is still under development, so please be careful.
    If you have any issues, or questions, feel free to open issues/PR-s.

## Fitting API

If you want to extend the package with a new primitive, this page will guide you what to define/extend.
You can also check the "built-in" primitives (´src/shapes/...´), which follow the same API.

The primitive must be subtype of `FittedShape` (`struct MyShape<:FittedShape`), and you must define the following functions:

- `defaultshapeparameters()`
- `fit()`
- `scorecandidate()`
- `refit!()`
- `strt()`

Docs of the functions:

```@docs
RANSAC.defaultshapeparameters
RANSAC.fit
RANSAC.scorecandidate
RANSAC.estimatescore
RANSAC.refit
RANSAC.strt
```

## Reading parameters from YAML file

If you want to read the parameters from a YAML file, you need to pass a dictionary to the [`readconfig`](@ref) function, that contains all the types, that you wish to use.
If you have `MyShape<:FittedShape` and the following YAML file

```yml
plane:
  - ϵ: 0.1
  - α: 0.01

myshape:
  - ϵ: 0.1
  - α: 0.01
  - anotherpar: 42

iteration:
  - prob_det: 0.999
  - τ: 10000
  - itermax: 100000
  - shape_types:
    - plane
    - myshape
```

you should use the following script to read the parameters:

```julia
julia> myd = Dict("plane"=>FittedPlane, "myshape"=>MyShape)
Dict{String,Type} with 2 entries:
  "plane"   => FittedPlane
  "myshape" => MyShape

julia> p = readconfig("config.yml", shapedict=myd)
(iteration = (drawN = 3, minsubsetN = 15, prob_det = 0.999,
shape_types = Type[FittedPlane, MyShape], τ = 10000, itermax = 100000,
extract_s = :nofminset, terminate_s = :nofminset),
common = (collin_threshold = 0.2, parallelthrdeg = 1.0),
plane = (ϵ = 0.1, α = 0.01),
cone = (ϵ = 0.3, α = 0.08726646259971647, minconeopang = 0.03490658503988659),
cylinder = (ϵ = 0.3, α = 0.08726646259971647),
sphere = (ϵ = 0.3, α = 0.08726646259971647, sphere_par = 0.02),
myshape = (ϵ = 0.1, α = 0.01, anotherpar = 42))
```

The parameters for `MyShape` can be accessed with `p.myshape`.
In any functions you need, the `@extract` macro can help from the [ExtractMacro.jl](https://github.com/carlobaldassi/ExtractMacro.jl) package.

An example:

```julia
function fit(::Type{MyShape}, p, n, params)
    @extract params: params_myshape=myshape
    @extract params_myshape : α, anotherpar
    # you can use α, anotherpar as you wish
    ...
end
```
