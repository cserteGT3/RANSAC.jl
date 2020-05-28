module RANSAC

using LinearAlgebra
using Random
using Logging
using Logging: default_logcolor
using StaticArrays: SVector, SMatrix, MMatrix, StaticArray
using RegionTrees
#using Images: label_components, component_lengths, component_subscripts
using ExtractMacro
import JSON
import YAML

import RegionTrees: AbstractRefinery, needs_refinement, refine_data
import Base: show, length, deleteat!

# debuuugggggg
#using Infiltrator
#using Makie: scatter, vbox

export  rodriguesdeg,
        rodriguesrad

export  arbitrary_orthogonal,
        arbitrary_orthogonal2,
        isparallel,
        findAABB,
        unitdisk2square,
        smallestdistance,
        prob,
        havesameelement

export  ConfidenceInterval,
        notsoconfident,
        isoverlap,
        E

export  OctreeNode,
        OctreeRefinery,
        RANSACCloud

export  FittedShape,
        FittedPlane,
        FittedSphere,
        FittedCylinder,
        FittedCone,
        ExtractedShape,
        ransacparameters

export  project2plane,
        refit

export  ransac

export  nosource_debuglogger,
        nosource_infologger,
        sourced_debuglogger,
        sourced_infologger,
        nosource_IterInflog,
        nosource_IterLow1log,
        nosource_IterLow2log,
        nosource_Compute1log,
        nosource_Compute2log,
        nosource_Errorlog

export  exportJSON,
        readconfig

include("utilities.jl")
include("confidenceintervals.jl")
include("octree.jl")
include("fitting.jl")
include("shapes/plane.jl")
include("shapes/sphere.jl")
include("shapes/cylinder.jl")
include("shapes/cone.jl")
#include("parameterspacebitmap.jl")
include("iterations.jl")
include("logging.jl")
include("orientedbox_.jl")
include("json-yaml.jl")

"""
`const DEFAULT_PARAMETERS = 
    defaultparameters([FittedPlane, FittedCone, FittedCylinder, FittedSphere])`
"""
const DEFAULT_PARAMETERS = defaultparameters([FittedPlane, FittedCone, FittedCylinder, FittedSphere])

"""
`const DEFAULT_SHAPE_DICT = Dict("plane"=>FittedPlane,
    "cone"=>FittedCone, "cylinder"=>FittedCylinder, "sphere"=>FittedSphere)`
"""
const DEFAULT_SHAPE_DICT = Dict("plane"=>FittedPlane, "cone"=>FittedCone, "cylinder"=>FittedCylinder, "sphere"=>FittedSphere)

end #module
