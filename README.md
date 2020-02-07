# RANSAC.jl

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://cserteGT3.github.io/PropertyFiles.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://cserteGT3.github.io/PropertyFiles.jl/dev)-->
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

The efficient RANSAC algorithm is used to detect primitive shapes in point clouds.
Primitive shapes can be: plane, sphere, cylinder, cone and torus, though torus is not implemented yet.

The input of the algorithm is a point cloud with associated surface normals.
The output is a set of primitive shapes with corresponding sets of points, and the rest of the points that do not belong to any primitives.

Here's an example with a point cloud and the detected primitives according to their type.
![RANSAC example](img/ransac_example.png)
