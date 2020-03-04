# Public API

Note, that not all exported functions are considered as part of the public API.

## Representing a point cloud



## Parameters

With the help of [Parameters.jl](https://github.com/mauro3/Parameters.jl) it's easy to parameterize the algorithm.
The `RANSACParameters` type collects all the parameters, though its fields are subject to change, the current fields and default values are listed below.

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
    # if the number of enabled points fall under `leftover`,
    # the iteration terminates
    leftovers::Int = 1

    # threshold of two vectors being parallel (in degrees)
    parallelthrdeg = 1
    # threshold of points being collinear
    collin_threshold = 0.2
    # parameter in sphere fitting
    sphere_par = 0.02

    # shapes that are fitted to the point cloud
    shape_types::Array{Symbol,1} = [:sphere, :plane, :cylinder, :cone]
end
```
