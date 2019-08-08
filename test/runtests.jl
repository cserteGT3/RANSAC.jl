using Test
using StaticArrays: SVector

#include("../src/RANSAC.jl")

@testset "iswithinrectangle test" begin
    include("rectangletest.jl")
end

@testset "dummy sphere test" begin
    include("dummyspheretest.jl")
end

@testset "utility tests" begin
    include("utilitytests.jl")
end
