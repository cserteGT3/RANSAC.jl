# Roadmap

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
* - [ ] sphere
* - [ ] cylinder
* - [ ] cone
* - [ ] torus
- [ ] add tests for every `isa...`
* - [ ] plane
* - [ ] sphere
* - [ ] cylinder
* - [ ] cone
* - [ ] torus
