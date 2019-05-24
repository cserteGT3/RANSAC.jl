# Roadmap

## General TODOs

- [ ] consistent use of `zchop`
- [ ] consistent use of thresholds
- [ ] define assertion of normed surface normals
* - [ ] decision
* - [ ] consistent implementation: delete unnecessary `normalize()`
- [x] add `cellandparents()` to `RegionTrees.jl`
* - [ ] writed test for it

## Generate examples

### Primitives
- [x] plane
- [x] sphere
- [x] cylinder
- [ ] cone
- [ ] torus

### Noisify
- [x] add gaussian noise to vertices
- [x] add gaussian noise to normals
- [ ] add outliers

### Visualize
- [x] fire up Makie
- [x] plot normals too

## Octree
- [ ] working Octree
- [x] `OctreeNode`
- [ ] Octree tests
- [x] `iswithinrectangle()` tests
- [ ] generalize `iswithinrectangle()`
- [ ] visualize Octree

### Tasks
- [x] use `Array` of `StaticArrays` (instead of 2dim array)
- [ ] check if `Point3f0` and 64 bit `Svector` is fine
- [ ] does changing the original cloudarray changes things?
- [ ] `struct PointCloud`: vertices, normals, enabling bits

## Fitting
- [ ] `isa...` for every shape
* - [x] plane
* - [x] sphere
* - [ ] cylinder
* - [ ] cone
* - [ ] torus
- [ ] add tests for every `isa...`
* - [ ] plane
* - [ ] sphere
* - [ ] cylinder
* - [ ] cone
* - [ ] torus
* - [ ] generalized interface for fitting

## Main loop

- [ ] add tests for `ConfidenceInterval`
- [ ] use of `Distributions.Hypergeometric`(?)
