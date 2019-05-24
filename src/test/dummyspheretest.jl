# Dummy test of sphere fitting

include("../fitting.jl")

using .Fitting
using StaticArrays: SVector
using Test

const EPSI = 0.1
const ALFI = deg2rad(10)

@testset "true sphere 1" begin
    tv1 = [SVector(0,-1,0.0), SVector(0,0,-1.0), SVector(1,0,0.0), SVector(0,1,0.0)]
    tn1 = [SVector(0,-1,0.0), SVector(0,0,-1.0), SVector(1,0,0.0), SVector(0,1,0.0)]
    fs = issphere(tv1, tn1, EPSI, ALFI)
    fp = isplane(tv1, tn1, π/2)
    @test fs.issphere == true
    @test fp.isplane == false
end

@testset "true sphere 2" begin
    tv1 = [SVector(0,-0.99,0.0), SVector(0,0,-1.0), SVector(1.01,0,0.0), SVector(0,1,0.0)]
    tn1 = [SVector(0,-1,0.0), SVector(0,0,-1.0), SVector(1,0,0.0), SVector(0,1,0.0)]
    fs1 = issphere(tv1, tn1, EPSI, ALFI) # true
    fs2 = issphere(tv1, tn1, 0.01, ALFI) # false, cause ϵ
    fp = isplane(tv1, tn1, π/2)
    @test fs1.issphere == true
    @test fs2.issphere == false
    @test fp.isplane == false
end

@testset "false sphere 1" begin
    tv1 = [SVector(0,1,0.0), SVector(0,0,-1.0), SVector(1,0,0.0), SVector(0,1,0.0)]
    tn1 = [SVector(0,-1,0.0), SVector(0,0,-1.0), SVector(1,0,0.0), SVector(0,1,0.0)]
    fs1 = issphere(tv1, tn1, EPSI, ALFI)
    # should be false even with large thresholds
    fs2 = issphere(tv1, tn1, 10, π/2)
    fp = isplane(tv1, tn1, π/2)
    @test fs1.issphere == false
    @test fs2.issphere == false
    @test fp.isplane == false
end
