# Example

This page guides you through the use of RANSAC.jl.
A public dataset is used, that was published by Schnabel et al. and is accessible [here](https://cg.cs.uni-bonn.de/en/publications/paper-details/schnabel-2009-completion/).
> R. Schnabel, P. Degener, R. Klein
> "Completion and Reconstruction with Primitive Shapes",
> in Computer Graphics Forum (Proc. of Eurographics), Vol. 28, No. 2, pages 503-512


Read the algorithm description here to get a full understanding of the algorithm and its parameters.


## Loading the data

As MeshIO and other softwares had troubles opening the `_input.obj` files, I used [MeshLab](http://www.meshlab.net/) to open them and export to non-binary encoded PLY.

## Constructing a `PointCloud`

Currently the package only handles vectors of `Float64`.

## Parameters

## Run!

## See the results

There's the [RANSACVisualizer](https://github.com/cserteGT3/RANSACVisualizer.jl) package that provides a few utility functions to check the results.
