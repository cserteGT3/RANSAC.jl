# RANSAC.jl

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg)-->
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://cserteGT3.github.io/RANSAC.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://cserteGT3.github.io/RANSAC.jl/dev)
[![Build Status](https://travis-ci.com/cserteGT3/RANSAC.jl.svg?branch=master)](https://travis-ci.com/cserteGT3/RANSAC.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/0wwq0nr9jhj2shq3/branch/master?svg=true)](https://ci.appveyor.com/project/cserteGT3/ransac-jl/branch/master)
[![codecov.io](http://codecov.io/github/cserteGT3/RANSAC.jl/coverage.svg?branch=master)](http://codecov.io/github/cserteGT3/RANSAC.jl?branch=master)


This package implements the efficient RANSAC algorithm for point clouds.
Paper can be found [here](https://cg.cs.uni-bonn.de/en/publications/paper-details/schnabel-2007-efficient/).

> 	R. Schnabel, R. Wahl, R. Klein
	"Efficient RANSAC for Point-Cloud Shape Detection",
	in Computer Graphics Forum, Vol. 26, No. 2, pages 214-226,
	Blackwell Publishing, June 2007

## Efficient RANSAC

The efficient RANSAC algorithm is used to segment and fit primitive shapes (sphere, plane, cylinder, torus, cone) to point clouds.
Up to now mostly C++ and Python implementations have been published, this is the first one in Julia (as far as I know).

### Main features

* easy-to-use primitive recognition
* extensible: it's easy to add new primitive shapes
* fast (to be honest, it's not yet comparable with the reference and CGAL implementations)

### Differences from the reference implementation

* no bitmap
* separate parameters for each shape
* no tori

Install the package by:

```julia
] add https://github.com/cserteGT3/RANSAC.jl
```

The input of the algorithm is a point cloud with associated surface normals.
The output is a set of primitive shapes with corresponding sets of points, and the rest of the points that do not belong to any primitives.

Here's an example with a point cloud and the detected primitives according to their type.
![RANSAC example](img/ransac_example.png)
