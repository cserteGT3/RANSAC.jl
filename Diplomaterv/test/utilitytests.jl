# Test functions from the Utilities module

include("../utilities.jl")

using .Utilities: findAABB
using StaticArrays: SVector
using Test

@testset "rand 3D points" begin
    rs3 = [SVector{3}(rand(3)) for _ in 1:30]
    iV3 = SVector(-1.0,-1,-1)
    aV3 = SVector(2.0,2.0,2.0)
    push!(rs3, iV3)
    push!(rs3, aV3)

    minV3, maxV3 = findAABB(rs3)
    @test minV3 == iV3
    @test maxV3 == aV3
end

@testset "rand 2D points" begin
    rs2 = [SVector{2}(rand(2)) for _ in 1:30]
    iV2 = SVector(-1.0,-1)
    aV2 = SVector(2.0,2.0)
    push!(rs2, iV2)
    push!(rs2, aV2)

    minV2, maxV2 = findAABB(rs2)
    @test minV2 == iV2
    @test maxV2 == aV2
end
