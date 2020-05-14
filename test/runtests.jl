using RANSAC
using Test
using LinearAlgebra
using StaticArrays: SVector
import JSON

include("testhelpers.jl")

@testset "octree tests" begin
    include("octree.jl")
end

@testset "dummy sphere test" begin
    include("dummyspheretest.jl")
end

@testset "utility tests" begin
    include("utilitytests.jl")
end

@testset "json export" begin
    include("json.jl")
end

@testset "yaml import" begin
    include("yaml.jl")
end

@testset "confidence interval" begin
    include("confidenceintervals.jl")
end

#=
@testset "parameterspace bitmap" begin
    include("parameterspacebitmap.jl")
end
=#
