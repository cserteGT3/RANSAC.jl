# Roadmap

## Generate examples

### Primitives
- [x] plane
- [ ] sphere
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
- [ ] `OctreeNode`
- [ ] Octree tests
- [ ] `iswithinrectangle()` test

### Tasks
- [ ] use `Array` of `StaticArrays` (instead of 2dim array)
- [ ] check if `Point3f0` and 64 bit `Svector` is fine
- [ ] does changing the original cloudarray changes things?
- [ ] `struc PointCloud`: vertices, normals, enabling bits
