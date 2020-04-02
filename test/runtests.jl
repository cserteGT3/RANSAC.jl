using RANSAC
using Test
using LinearAlgebra
using StaticArrays: SVector
import JSON

@testset "iswithinrectangle test" begin
    include("rectangletest.jl")
end

@testset "dummy sphere test" begin
    include("dummyspheretest.jl")
end

@testset "utility tests" begin
    include("utilitytests.jl")
end

@testset "cone" begin
    include("cone.jl")
end

@testset "json export" begin
    include("json.jl")
end

#=
@testset "parameterspace bitmap" begin
    include("parameterspacebitmap.jl")
end
=#
