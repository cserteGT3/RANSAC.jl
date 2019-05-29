module RANSAC

using LinearAlgebra
using Random: randperm, shuffle
using Logging
using StaticArrays: SVector, MVector
using RegionTrees
using ZChop: zchop, zchop!
using AbstractPlotting: Point3f0
using Makie: scatter, linesegments!, Scene
using Images: label_components, component_lengths

import RegionTrees: AbstractRefinery, needs_refinement, refine_data

export  rodriguesdeg,
        rodriguesrad

export  arbitrary_orthogonal,
        isparallel,
        findAABB,
        chopzaxis,
        unitdisk2square,
        prob,
        smallestdistance,
        havesameelement,
        AlphSilon

export  ConfidenceInterval,
        notsoconfident,
        estimatescore,
        smallestdistance,
        E

export  OctreeNode,
        iswithinrectangle,
        OctreeRefinery,
        PointCloud,
        getcellandparents,
        allparents,
        updatelevelweight

export  FittedShape,
        isshape,
        FittedPlane,
        isplane,
        FittedSphere,
        issphere,
        FittedCylinder,
        iscylinder,
        ShapeCandidate,
        findhighestscore,
        ScoredShape,
        largestshape

export  sampleplane,
        samplecylinder,
        samplesphere,
        normalsforplot,
        noisifyvertices,
        noisifynormals,
        makemeanexample,
        examplepc2,
        examplepc3,
        examplepc4,
        examplepc5,
        showgeometry

export  project2plane,
        compatiblesPlane,
        bitmapparameters,
        compatiblesSphere,
        compatiblesCylinder,
        largestconncomp,
        refitsphere,
        refitplane,
        refitcylinder

export  ransac,
        showcandlength,
        showshapes,
        getrest,
        showtype,
        showbytype

include("utilities.jl")
include("ConfidenceIntervals.jl")
include("octree.jl")
include("fitting.jl")
include("sampling.jl")
include("parameterspacebitmap.jl")
include("iterations.jl")

end #module