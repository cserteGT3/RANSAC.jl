# Example

This page guides you through the use of RANSAC.jl.
A public dataset is used, that was published by Schnabel et al. in ... and is accessible here.
Read the algorithm description here to get a full understanding of the algorithm and its parameters.

## Loading the data

As MeshIO and other softwares had troubles opening the `_input.obj` files, I used [MeshLab]() to open them and export to non-binary encoded PLY.

## Constructing a `PointCloud`

Currently the package only handles vectors of `Float64`.

## Parameters

## Run!

## See the results

There's the [RANSACVisualizer]() package that provides a few utility functions to check the results.
