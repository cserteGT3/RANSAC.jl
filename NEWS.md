# NEWS

In this file I keep track of the changes from release to release.

## v0.6.0

#### Bugfixes

- handling negative sqrt error - closes #29

#### Changes

- progress is shown with the ProgressMeter.jl
- documentation is updated to use the new RANSACVisualizer package with Makie recipes
- the Readme have been extended with a related projects section

## v0.5.0

#### Changes

- **BREAKING**: pass `RANSACCloud` too when fitting shapes - [docs](https://csertegt3.github.io/RANSAC.jl/stable/newprimitive/)

### v0.5.1

#### Changes

- enable returning multiple shapes when fitting - [docs](https://csertegt3.github.io/RANSAC.jl/stable/newprimitive/)

### v0.5.2

#### Changes

- I thought that there's a license file in the repo. I was wrong so I fixed that.

## v0.4.0

#### Changes

- complete rework of the fitting and `FittedShape` API - [docs](https://csertegt3.github.io/RANSAC.jl/stable/newprimitive/)
- rework parameters: change Parameters.jl to named tuples (and `@extract` macro) - [docs](https://csertegt3.github.io/RANSAC.jl/stable/api/#Parameters-1)
- replace `PointCloud` with `RANSACCloud` and merge octree into it - [docs](https://csertegt3.github.io/RANSAC.jl/stable/api/#Representing-a-point-cloud-1)
- speed optimizations
- removed the sampling functions, that generate shapes for testing - see [RANSACBenchmark.jl](https://github.com/cserteGT3/RANSACBenchmark.jl)

## v0.3.0

#### New feature

- `exportJSON` function to export the reconstructed shapes to JSON.

#### Changes

- removed the `is_shape_` field of `Fitted_Shape_`s, which is a breaking change, but should not cause issues, unless you constructed such types manually. Then you should remove the first argument from the constructor calls.
- with the above change `isshape()` has been removed and `fit_shape_` functions now return `nothing` if it can't fit the given primitive.

### v0.3.1

#### New feature

- `readconfig` function that reads parameters from a YAML file

#### Changes

- cleared the parameter struct, unused parameters were removed

#### Bugfixes

- resolved #5

## v0.2.0

I removed some undocumented features (reconstructing translational surfaces) that had unregistered dependencies.
This way I don't have to check in the Manifest.

## v0.1.0

First release. The code is in the state as I finished my master's thesis.
Documentation is still WIP.

### v0.1.1

Finishing touches to the documentation.
