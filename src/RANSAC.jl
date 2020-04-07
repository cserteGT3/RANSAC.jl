module RANSAC

using LinearAlgebra
using Random
using Logging
using Logging: default_logcolor
using StaticArrays: SVector, SMatrix
using RegionTrees
#using Images: label_components, component_lengths, component_subscripts
using Parameters
import JSON
import YAML

import RegionTrees: AbstractRefinery, needs_refinement, refine_data
import Base: show

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
        havesameelement,
        #p2table,
        RANSACParameters,
        setepsilons,
        setalphas

export  ConfidenceInterval,
        notsoconfident,
        isoverlap,
        E

export  OctreeNode,
        iswithinrectangle,
        OctreeRefinery,
        PointCloud,
        getcellandparents,
        octreedepth,
        updatelevelweight

export  FittedShape,
        FittedPlane,
        FittedSphere,
        FittedCylinder,
        FittedCone,
        ShapeCandidate,
        findhighestscore,
        ScoredShape,
        largestshape,
        forcefitshapes

export  sampleplane,
        sampleplanefromcorner,
        samplecylinder,
        samplesphere,
        samplecone,
        normalsforplot,
        noisifyvertices,
        noisifynormals,
        makemeanexample,
        examplepc2,
        examplepc3,
        examplepc4,
        examplepc5,
        examplepc6

export  project2plane,
        refitsphere,
        refitplane,
        refitcylinder,
		refitcone

export  ransac,
        rerunleftover!

export  nosource_debuglogger,
        nosource_infologger,
        sourced_debuglogger,
        sourced_infologger,
        nosource_IterInflog,
        nosource_IterLow1log,
        nosource_IterLow2log,
        nosource_Compute1log,
        nosource_Compute2log,
        nosource_Errorlog,
        nosource_metafmt

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
include("sampling.jl")
#include("parameterspacebitmap.jl")
include("iterations.jl")
include("logging.jl")
include("orientedbox_.jl")
include("json-yaml.jl")

end #module
